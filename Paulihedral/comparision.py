import sys, time
from benchmark.offline import *
from benchmark.molecule import *
from benchmark.ising import *
from parallel_bl import *
from qiskit import QuantumCircuit, transpile
from qiskit_nature.mappers.second_quantization import JordanWignerMapper, ParityMapper, BravyiKitaevMapper
from qiskit_nature.problems.second_quantization.lattice import (
    BoundaryCondition,
    FermiHubbardModel,
    HyperCubicLattice,
    Lattice,
    LatticeDrawStyle,
    LineLattice,
    SquareLattice,
    TriangularLattice,
)
import synthesis_SC
import synthesis_FT
from tools import *
from arch import *
from pytket.pauli import Pauli, QubitPauliString
from pytket.circuit import Qubit, Node
from pytket import Circuit, OpType
from pytket.circuit import PauliExpBox
from pytket.passes import DecomposeBoxes, PauliSimp, SynthesiseTket, FullPeepholeOptimise, RoutingPass, DecomposeSwapsToCXs
from pytket.qasm import circuit_from_qasm, circuit_from_qasm_str, circuit_from_qasm_io, circuit_to_qasm_str, circuit_to_qasm
import time, sys, os
from pytket.routing import place_with_map
from pytket.routing import Placement, LinePlacement, GraphPlacement, NoiseAwarePlacement
from config import test_scale
import re
import numpy as np
import random


old_cwd = os.getcwd()
set_cwd()
ctime = time.time

local2Regex = re.compile(r'[XYZxyz]\d+[XYZxyz]\d+')
qubitRegex = re.compile(r'\d+')
gateRegex = re.compile(r'[XYZ]')


def gene_FermiHubbard_oplist(rows, cols):
    qubit_converter = QubitConverter(mapper=BravyiKitaevMapper())
    square_lattice = SquareLattice(rows=rows, cols=cols, boundary_condition=BoundaryCondition.PERIODIC)
    t = -1.0  # the interaction parameter
    v = 0.0  # the onsite potential
    u = 5.0  # the interaction parameter U
    fhm = FermiHubbardModel.uniform_parameters(
        lattice=square_lattice,
        uniform_interaction=t,
        uniform_onsite_potential=v,
        onsite_interaction=u,
    )
    ham = qubit_converter.convert(fhm.second_q_ops(display_format="sparse"))
    oplist = []
    for i in ham.primitive:
        oplist.append([pauliString(str(i.paulis[0]), coeff=1.0)])
    return oplist

def gene_random_oplist(qubit, order=3, seed=10):
    s = min(np.math.factorial(qubit), int(qubit**order), int(4**qubit))
    random.seed(seed)
    oplist = []
    ps = ['I', 'X', 'Y', 'Z']
    while len(oplist) < s:
        pss = ""
        for i in range(qubit):
            pss += random.choice(ps)
        if pss not in oplist:
            oplist.append(pss)
    res = []
    for i in oplist:
        res.append([pauliString(i, coeff=1.0)])
    return res

def gene_cond_random_oplist(qubit, num, seed=10):
    s = min(np.math.factorial(qubit), num, int(4**qubit))
    random.seed(seed)
    oplist = []
    ps = "XYZ"
    while len(oplist) < s:
        t = random.randint(1, qubit//2)
        t0 = random.randint(0, t)
        pss = ""
        for i in range(qubit-2*t):
            pss += random.choice(ps)
        pss = 2*t0*'I' + pss + 2*(t-t0)*'I'
        if pss not in oplist:
            oplist.append(pss)
    res = []
    for i in oplist:
        res.append([pauliString(i, coeff=1.0)])
    return res

# PH no device method
def PH(parr):
    print('PH passes, Our schedule, Our synthesis', flush=True)
    nq = len(parr[0][0])
    length = nq//2 # `length' is a hyperparameter, and can be adjusted for best performance
    t0 = ctime()
    a1 = depth_oriented_scheduling(parr, length=length, maxiter=30)
    qc = synthesis_FT.block_opt_FT(a1)
    print('PH, Time costed:', ctime()-t0, flush=True)
    qc1 = transpile(qc, basis_gates=['u3', 'cx'], optimization_level=0)
    t0 = ctime()
    qc2 = transpile(qc, basis_gates=['u3', 'cx'], optimization_level=3)
    print_qc(qc2)
    print('Qiskit L3, Time costed:', ctime()-t0, flush=True)
    ansatz = circuit_from_qasm_str(qc1.qasm())
    t0 = ctime()
    FullPeepholeOptimise().apply(ansatz)
    SynthesiseTket().apply(ansatz)
    print(f"CNOT: {ansatz.n_gates_of_type(OpType.CX)}, Single: {ansatz.n_gates-ansatz.n_gates_of_type(OpType.CX)}, Total: {ansatz.n_gates}, Depth: {ansatz.depth()}")
    print("TKet O2:", ctime()-t0)
# TKet no device method
def TK(parr):
    print('TK passes', flush=True)
    n = len(parr[0][0])
    q = [Qubit(i) for i in range(n)]
    oplist = {}
    def to_pauli_list(ps):
        r = []
        for i in ps:
            if i == 'I':
                r.append(Pauli.I)
            elif i == 'X':
                r.append(Pauli.X)
            elif i == 'Y':
                r.append(Pauli.Y)
            elif i == 'Z':
                r.append(Pauli.Z)
        return r
    for i in parr:
        for j in i:
            op = QubitPauliString(q, to_pauli_list(j.ps))
            oplist[op] = 1/3.14
    def add_excitation(circ, term_dict, param=1.0):
        for term, coeff in term_dict.items():
            qubits, paulis = zip(*term.to_dict().items())
            pbox = PauliExpBox(paulis, coeff * param)
            circ.add_pauliexpbox(pbox, qubits)
    ansatz = Circuit(n)
    t0 = ctime()
    add_excitation(ansatz, oplist)
    PauliSimp().apply(ansatz)
    print(f"TK Pauli Simp: {ctime()-t0}")
    t0 = ctime()
    FullPeepholeOptimise().apply(ansatz)
    print(f"CNOT: {ansatz.n_gates_of_type(OpType.CX)}, Single: {ansatz.n_gates-ansatz.n_gates_of_type(OpType.CX)}, Total: {ansatz.n_gates}, Depth: {ansatz.depth()}")
    print("TK O2:", ctime()-t0)
    ansatz = Circuit(n)
    add_excitation(ansatz, oplist)
    PauliSimp().apply(ansatz)
    SynthesiseTket().apply(ansatz)
    qstr = circuit_to_qasm_str(ansatz)
    qc = QuantumCircuit.from_qasm_str(qstr)
    t0 = ctime()
    qc = transpile(qc, basis_gates=['cx', 'u3'], optimization_level=3)
    print_qc(qc)
    print("Qiskit L3:", ctime()-t0)

# generate Hamiltonian from string
def gene_oplist_strlist(sl):
    oplist = []
    if len(qubitRegex.findall(str(sl))) > 0:
        pslist = local2Regex.findall(str(sl))
        qubitlist = qubitRegex.findall(str(sl))
        max_qubit = max(list(map(int, qubitlist)))
        for i in pslist:
            gates = gateRegex.findall(str(i))
            qubitindex = qubitRegex.findall(str(i))
            index = list(map(int, qubitindex))
            ps = index[0]*'I' + gates[0] + (index[1]-index[0]-1)*'I' + gates[1] + (max_qubit-index[1])*'I'
            oplist.append([pauliString(ps, coeff=1.0)])
    else:
        for i in sl:
            oplist.append([pauliString(i, coeff=1.0)])
    return oplist




#parr = [[pauliString('XXXII', 1.0),pauliString('XXIII', 1.0)],[pauliString('XXIXI', 1.0),pauliString('XXIIX', 1.0)]]
#parr1 = gene_molecule_oplist(H2)
#heisen_parr = [gene_dot_1d(29, interaction='Z')+gene_dot_1d(29, interaction='X')+gene_dot_1d(29, interaction='Y'), gene_dot_2d(4,5, interaction='Z')+gene_dot_2d(4,5, interaction='Y')+gene_dot_2d(4,5, interaction='X'), gene_dot_3d(1,2,4, interaction='Z')+gene_dot_3d(1,2,4, interaction='Y')+gene_dot_3d(1,2,4, interaction='X')]
#print(heisen_parr[0])
with open('benchmark/circuits/randomcircuit_q05_06.txt') as file:
    localdata = file.read()
local2gate = local2Regex.findall(localdata)
parr2 = gene_oplist_strlist(local2gate)
#print(parr2)
#parr = load_oplist('MgH', benchmark='uccsd')
#parrraw = load_oplist('LiH', benchmark='molecule')
#parr = gene_FermiHubbard_oplist(3,3)
#parr3 = gene_random_oplist(10)
#print(parrraw[0][0].coeff)
#PH(parr2)
#TK(parrraw)