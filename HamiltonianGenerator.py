from qiskit_nature.drivers import UnitsType, Molecule
from qiskit_nature.drivers.second_quantization import (
    ElectronicStructureDriverType,
    ElectronicStructureMoleculeDriver,
)
from qiskit_nature.problems.second_quantization import ElectronicStructureProblem
from qiskit_nature.problems.second_quantization.lattice import FermiHubbardModel
from qiskit_nature.converters.second_quantization import QubitConverter
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
import os
import pickle
import re
import numpy as np
import random



gateRegex = re.compile(r'[IXYZ]')
pauliRegex = re.compile(r'[XYZ]\d+')
qubitRegex = re.compile(r'\d+')
package_directory = os.path.dirname(os.path.abspath(__file__))

H2 = [['H', [0., 0, 0]],['H', [0, 0, -1.5]]]
He2 = [['He', [0., 0, 0]],['He', [0, 0, 2.6740000]]]
HF = [['F', [0., 0, 0]],['H', [0, 0, 0.9153800]]]
HCl = [['Cl', [0., 0., 0.]],['H', [0, 0, 1.2744500]]]
ClF = [['F', [0., 0., 0.]],['Cl', [0., 0., 1.6303300]]]
LiH = [['Li', [0,0,0]],['H', [0,0,0.76]]]
N2 = [['N', [0., 0, 0]],['N', [0, 0, -1.5]]]
NH2 = [['N', [0.,0.,0.]],['H', [0.,0.,1.0240500]],['H', [0.9971600,0.,-0.2331400]]]  # NH2 charge = 0, multiplicity=2
NH3 = [['N', [0.,0.,0.1128900]],['H', [0., 0.9380200,-0.2634100]],
    ['H', [0.8123500,-0.4690100,-0.2634100]],
    ['H', [-0.8123500,-0.4690100,-0.2634100]]]
H2O = [['O', [0.,0.,0.]],['H', [0.95,-0.55,0.]],['H', [-0.95,-0.55,0.]]]
H2S = [['S', [0., 0, 0]],['H', [0, 0, -1.5]],['H', [0, 0, 1.5]]]
CH2 = [['C', [0.,0.,0.]],['H', [0.,0.,1.1077100]],['H', [1.0837800,0.,-0.2289900]]]
CH4 = [['C', [0.,1.0,1.0]],['H', [0.051054399,0.051054399,0.051054399]],
    ['H', [1.948945601,1.948945601,0.051054399]],
    ['H', [0.051054399,1.948945601,1.948945601]],
    ['H', [1.948945601,0.051054399,1.948945601]]]
MgO = [['Mg', [0., 0, 0]],['O', [0, 0, -1.5]]]
CO2 = [['O', [1.4, 0., 0.]],['C', [0., 0., 0.]],['O', [-1.4, 0., 0.]]]
NaCl = [['Na', [0., -1.5, -1.5]],['Cl', [1.5, -1.5, -1.5]]]
KOH = [['K', [3.4030, 0.2500, 0.0]],['O', [2.5369, -0.2500, -1.5]],['H', [2.0,0.06,0.0]]]
FeO = [['Fe', [0., 0, 0]],['O', [0, 0, -1.5]]]


def load_oplist(name, benchmark='uccsd'):
    if benchmark == 'molecule':
        fth = os.path.join(package_directory, 'data', name + '.pickle')
        with open(fth, 'rb') as f:
            entry = pickle.load(f)
        oplist = []
        paulis = []
        for i in entry:
            paulis.extend(i)
        for i in paulis:
            oplist.append(i.ps)
        return oplist
    if benchmark == 'uccsd':
        fth = os.path.join(package_directory, 'data', name + '_UCCSD.pickle')
        with open(fth, 'rb') as f:
            entry = pickle.load(f)
        oplist = []
        paulis = []
        for i in entry:
            paulis.extend(i)
        for i in paulis:
            oplist.append(i.ps)
        return oplist

def get_qubit_op(geo, molecule_charge, molecule_multiplicity):
    molecule = Molecule(geometry=geo, charge=molecule_charge, multiplicity=molecule_multiplicity)
    driver = ElectronicStructureMoleculeDriver(
    molecule, basis="sto3g", driver_type=ElectronicStructureDriverType.PYSCF)
    es_problem = ElectronicStructureProblem(driver)
    second_q_op = es_problem.second_q_ops()
    # JordanWigner mapping
    qubit_converter = QubitConverter(mapper=BravyiKitaevMapper())
    qubit_op = qubit_converter.convert(second_q_op[0])
    # Parity mapping
    '''qubit_converter = QubitConverter(mapper=ParityMapper(), two_qubit_reduction=True)
    qubit_op = qubit_converter.convert(second_q_op[0], num_particles=es_problem.num_particles)'''

    return qubit_op

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
        oplist.append(str(i.paulis[0]))
    return oplist


def gene_molecule_oplist(atom_geo, molecule_charge=0, molecule_multiplicity=1):
    qubit_op = get_qubit_op(atom_geo, molecule_charge, molecule_multiplicity)
    oplist = []
    for i in qubit_op.primitive:
        oplist.append(str(i.paulis[0]))
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
        res.append(i)
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
        res.append(i)
    return res

def gene_dot_1d(w, interaction='Z'):
    nq = w + 1
    oplist = []
    for i in range(nq - 1):
        ps = i*'I' + interaction + interaction + (nq-2-i)*'I'
        oplist.append(ps)
    return oplist

# w >= 2, h >= 2
def gene_dot_2d(w, h, offset=0, numq=0, interaction='Z'):
    if numq == 0:
        nq = (w+1)*(h+1)
    else:
        nq = numq
    oplist = []
    for i in range(w):
        for j in range(h):
            k = (w+1)*j + i + offset
            ps = ['I']*nq
            ps[k] = interaction
            ps[k+1] = interaction
            oplist.append("".join(ps))
            ps = ['I']*nq
            ps[k] = interaction
            ps[k+1+w] = interaction
            oplist.append("".join(ps))
    for j in range(h):
        k = (w+1)*j + w + offset
        ps = ['I']*nq
        ps[k] = interaction
        ps[k+1+w] = interaction
        oplist.append("".join(ps))
    for i in range(w):
        k = h*(w+1) + i + offset
        ps = ['I']*nq
        ps[k] = interaction
        ps[k+1] = interaction
        oplist.append("".join(ps))
    return oplist

# w, h, l >= 2
# w*(h+1)*(l+1) + (w+1)*h*(l+1) + (w+1)*(h+1)*l
def gene_dot_3d(w, h, l, interaction='Z'):
    oplist = []
    nq = (w+1)*(h+1)*(l+1)
    for i in range(l+1):
        oplist += gene_dot_2d(w, h, numq=nq, offset=i*(w+1)*(h+1), interaction=interaction)
    for i in range(l):
        for j in range((w+1)*(h+1)):
            k = i*(w+1)*(h+1) + j
            ps = ['I']*nq
            ps[k] = interaction
            ps[k+(w+1)*(h+1)] = interaction
            oplist.append("".join(ps))
    return oplist

