from qiskit import QuantumCircuit, QuantumRegister, transpile
from qiskit.circuit import Parameter
import numpy as np
import sys, time, os
import re
from pytket.pauli import Pauli, QubitPauliString
from pytket.circuit import Qubit, Node
from pytket import Circuit, OpType
from pytket.circuit import PauliExpBox
from pytket.passes import DecomposeBoxes, PauliSimp, SynthesiseTket, FullPeepholeOptimise, RoutingPass, DecomposeSwapsToCXs
from pytket.qasm import circuit_from_qasm, circuit_from_qasm_str, circuit_from_qasm_io, circuit_to_qasm_str, circuit_to_qasm
from Paulihedral.parallel_bl import *
import Paulihedral.synthesis_FT as synthesis_FT
from Paulihedral.tools import *
from Paulihedral.arch import *
from Paulihedral.benchmark.mypauli import *
from TemplateMatching import *


ctime = time.time
#----------------------------------------------------------------------------------------
#                Parameter of 2-local Hamiltonian Evolution Gate (γ)
#----------------------------------------------------------------------------------------
gamma = Parameter("$\\gamma$")
gamma = 1
#----------------------------------------------------------------------------------------
#                           Define useful Regex compilers
#----------------------------------------------------------------------------------------
# Regex to find local gates and qubit index
local2Regex = re.compile(r'[XYZ]\d+[XYZ]\d+')
qubitRegex = re.compile(r'\d+')
gateRegex = re.compile(r'[XYZ]')
sgateRegex = re.compile(r'[XYZ]\d+')
reduceRegex = re.compile(r'\d+[lr]*')
pauliRegex = re.compile(r'[XYZ]\d+')
allgateRegex = re.compile(r'[IXYZ]')
#----------------------------------------------------------------------------------------
#          The parameter is set to 1 by default since it must be fixed and has no
#                          influence on Optimization Algorithm.
#----------------------------------------------------------------------------------------


# 2local Circuit Types
def GateXX(nqubits,i,j,first='none',second='none',cut='no'):
    if cut not in ['no','yes'] :
        print('Unexpected Arguments: '+cut)
        sys.exit(0)
    if first not in ['l','r','lr','none'] or second not in ['l','r','lr','none'] :
        if first not in ['l','r','lr','none'] and second not in ['l','r','lr','none'] :
            print('Unexpected Arguments: '+first+' , '+second)
        elif first not in ['l','r','lr','none'] and second in ['l','r','lr','none']:
            print('Unexpected Arguments: '+first)
        elif second not in ['l','r','lr','none'] and first in ['l','r','lr','none']:
            print('Unexpected Arguments: '+second)       
        sys.exit(0)

    qc_xx = QuantumCircuit(nqubits)
    if (first == 'l' or first == 'lr') and (second == 'r' or second == 'none') :
        qc_xx.h(nqubits[j])
    elif (first == 'r' or first == 'none') and (second == 'l' or second == 'lr') :
        qc_xx.h(nqubits[i])
    elif (first == 'r' or first == 'none') and (second == 'r' or second == 'none') :
        qc_xx.h(nqubits[i])
        qc_xx.h(nqubits[j])

    if cut == 'no' :
        qc_xx.cx(nqubits[i],nqubits[j])
    qc_xx.rz(gamma,nqubits[j])
    qc_xx.cx(nqubits[i],nqubits[j])

    if (first == 'r' or first == 'lr') and (second == 'l' or second == 'none') :
        qc_xx.h(nqubits[j])
    elif (first == 'l' or first == 'none') and (second == 'r' or second == 'lr') :
        qc_xx.h(nqubits[i])
    elif (first == 'l' or first == 'none') and (second == 'l' or second == 'none') :
        qc_xx.h(nqubits[i])
        qc_xx.h(nqubits[j])       

    return qc_xx

def GateXY(nqubits,i,j,first='none',second='none',cut='no'):
    if cut not in ['no','yes'] :
        print('Unexpected Arguments: '+cut)
        sys.exit(0)
    if first not in ['l','r','lr','none'] or second not in ['l','r','lr','none'] :
        if first not in ['l','r','lr','none'] and second not in ['l','r','lr','none'] :
            print('Unexpected Arguments: '+first+' , '+second)
        elif first not in ['l','r','lr','none'] and second in ['l','r','lr','none']:
            print('Unexpected Arguments: '+first)
        elif second not in ['l','r','lr','none'] and first in ['l','r','lr','none']:
            print('Unexpected Arguments: '+second)       
        sys.exit(0)

    qc_xy = QuantumCircuit(nqubits)
    if (first == 'l' or first == 'lr') and (second == 'r' or second == 'none') :
        qc_xy.rx(-np.pi/2,nqubits[j])
    elif (first == 'r' or first == 'none') and (second == 'l' or second == 'lr') :
        qc_xy.h(nqubits[i])
    elif (first == 'r' or first == 'none') and (second == 'r' or second == 'none') :
        qc_xy.h(nqubits[i])
        qc_xy.rx(-np.pi/2,nqubits[j])
    
    if cut == 'no' :
        qc_xy.cx(nqubits[i],nqubits[j])
    qc_xy.rz(gamma,nqubits[j])
    qc_xy.cx(nqubits[i],nqubits[j])

    if (first == 'r' or first == 'lr') and (second == 'l' or second == 'none') :
        qc_xy.rx(np.pi/2,nqubits[j])
    elif (first == 'l' or first == 'none') and (second == 'r' or second == 'lr') :
        qc_xy.h(nqubits[i])
    elif (first == 'l' or first == 'none') and (second == 'l' or second == 'none') :
        qc_xy.h(nqubits[i])
        qc_xy.rx(np.pi/2,nqubits[j])

    return qc_xy

def GateXZ(nqubits,i,j,first='none',second='none',cut='no'):
    if cut not in ['no','yes'] :
        print('Unexpected Arguments: '+cut)
        sys.exit(0)
    if first not in ['l','r','lr','none'] or second not in ['l','r','lr','none'] :
        if first not in ['l','r','lr','none'] and second not in ['l','r','lr','none'] :
            print('Unexpected Arguments: '+first+' , '+second)
        elif first not in ['l','r','lr','none'] and second in ['l','r','lr','none']:
            print('Unexpected Arguments: '+first)
        elif second not in ['l','r','lr','none'] and first in ['l','r','lr','none']:
            print('Unexpected Arguments: '+second)       
        sys.exit(0)

    qc_xz = QuantumCircuit(nqubits)
    if (first == 'r' or first == 'none') and (second == 'l' or second == 'lr') :
        qc_xz.h(nqubits[i])
    elif (first == 'r' or first == 'none') and (second == 'r' or second == 'none') :
        qc_xz.h(nqubits[i])

    if cut == 'no' :
        qc_xz.cx(nqubits[i],nqubits[j])
    qc_xz.rz(gamma,nqubits[j])
    qc_xz.cx(nqubits[i],nqubits[j])

    if (first == 'l' or first == 'none') and (second == 'r' or second == 'lr') :
        qc_xz.h(nqubits[i])
    elif (first == 'l' or first == 'none') and (second == 'l' or second == 'none') :
        qc_xz.h(nqubits[i])

    return qc_xz

def GateYX(nqubits,i,j,first='none',second='none',cut='no'):
    if cut not in ['no','yes'] :
        print('Unexpected Arguments: '+cut)
        sys.exit(0)
    if first not in ['l','r','lr','none'] or second not in ['l','r','lr','none'] :
        if first not in ['l','r','lr','none'] and second not in ['l','r','lr','none'] :
            print('Unexpected Arguments: '+first+' , '+second)
        elif first not in ['l','r','lr','none'] and second in ['l','r','lr','none']:
            print('Unexpected Arguments: '+first)
        elif second not in ['l','r','lr','none'] and first in ['l','r','lr','none']:
            print('Unexpected Arguments: '+second)       
        sys.exit(0)

    qc_yx = QuantumCircuit(nqubits)
    if (first == 'l' or first == 'lr') and (second == 'r' or second == 'none') :
        qc_yx.h(nqubits[j])
    elif (first == 'r' or first == 'none') and (second == 'l' or second == 'lr') :
        qc_yx.rx(-np.pi/2,nqubits[i])
    elif (first == 'r' or first == 'none') and (second == 'r' or second == 'none') :
        qc_yx.rx(-np.pi/2,nqubits[i])
        qc_yx.h(nqubits[j])

    if cut == 'no' :
        qc_yx.cx(nqubits[i],nqubits[j])
    qc_yx.rz(gamma,nqubits[j])
    qc_yx.cx(nqubits[i],nqubits[j])

    if (first == 'r' or first == 'lr') and (second == 'l' or second == 'none') :
        qc_yx.h(nqubits[j])
    elif (first == 'l' or first == 'none') and (second == 'r' or second == 'lr') :
        qc_yx.rx(np.pi/2,nqubits[i])
    elif (first == 'l' or first == 'none') and (second == 'l' or second == 'none') :
        qc_yx.rx(np.pi/2,nqubits[i])
        qc_yx.h(nqubits[j])

    return qc_yx

def GateYY(nqubits,i,j,first='none',second='none',cut='no'):
    if cut not in ['no','yes'] :
        print('Unexpected Arguments: '+cut)
        sys.exit(0)
    if first not in ['l','r','lr','none'] or second not in ['l','r','lr','none'] :
        if first not in ['l','r','lr','none'] and second not in ['l','r','lr','none'] :
            print('Unexpected Arguments: '+first+' , '+second)
        elif first not in ['l','r','lr','none'] and second in ['l','r','lr','none']:
            print('Unexpected Arguments: '+first)
        elif second not in ['l','r','lr','none'] and first in ['l','r','lr','none']:
            print('Unexpected Arguments: '+second)       
        sys.exit(0)

    qc_yy = QuantumCircuit(nqubits)
    if (first == 'l' or first == 'lr') and (second == 'r' or second == 'none') :
        qc_yy.rx(-np.pi/2,nqubits[j])
    elif (first == 'r' or first == 'none') and (second == 'l' or second == 'lr') :
        qc_yy.rx(-np.pi/2,nqubits[i])
    elif (first == 'r' or first == 'none') and (second == 'r' or second == 'none') :
        qc_yy.rx(-np.pi/2,nqubits[i])
        qc_yy.rx(-np.pi/2,nqubits[j])

    if cut == 'no' :
        qc_yy.cx(nqubits[i],nqubits[j])
    qc_yy.rz(gamma,nqubits[j])
    qc_yy.cx(nqubits[i],nqubits[j])

    if (first == 'r' or first == 'lr') and (second == 'l' or second == 'none') :
        qc_yy.rx(np.pi/2,nqubits[j])
    elif (first == 'l' or first == 'none') and (second == 'r' or second == 'lr') :
        qc_yy.rx(np.pi/2,nqubits[i])
    elif (first == 'l' or first == 'none') and (second == 'l' or second == 'none') :
        qc_yy.rx(np.pi/2,nqubits[i])
        qc_yy.rx(np.pi/2,nqubits[j])

    return qc_yy

def GateYZ(nqubits,i,j,first='none',second='none',cut='no'):
    if cut not in ['no','yes'] :
        print('Unexpected Arguments: '+cut)
        sys.exit(0)
    if first not in ['l','r','lr','none'] or second not in ['l','r','lr','none'] :
        if first not in ['l','r','lr','none'] and second not in ['l','r','lr','none'] :
            print('Unexpected Arguments: '+first+' , '+second)
        elif first not in ['l','r','lr','none'] and second in ['l','r','lr','none']:
            print('Unexpected Arguments: '+first)
        elif second not in ['l','r','lr','none'] and first in ['l','r','lr','none']:
            print('Unexpected Arguments: '+second)       
        sys.exit(0)

    qc_yz = QuantumCircuit(nqubits)
    if (first == 'r' or first == 'none') and (second == 'l' or second == 'lr') :
        qc_yz.rx(-np.pi/2,nqubits[i])
    elif (first == 'r' or first == 'none') and (second == 'r' or second == 'none') :
        qc_yz.rx(-np.pi/2,nqubits[i])

    if cut == 'no' :
        qc_yz.cx(nqubits[i],nqubits[j])
    qc_yz.rz(gamma,nqubits[j])
    qc_yz.cx(nqubits[i],nqubits[j])

    if (first == 'l' or first == 'none') and (second == 'r' or second == 'lr') :
        qc_yz.rx(np.pi/2,nqubits[i])
    elif (first == 'l' or first == 'none') and (second == 'l' or second == 'none') :
        qc_yz.rx(np.pi/2,nqubits[i])

    return qc_yz

def GateZX(nqubits,i,j,first='none',second='none',cut='no'):
    if cut not in ['no','yes'] :
        print('Unexpected Arguments: '+cut)
        sys.exit(0)
    if first not in ['l','r','lr','none'] or second not in ['l','r','lr','none'] :
        if first not in ['l','r','lr','none'] and second not in ['l','r','lr','none'] :
            print('Unexpected Arguments: '+first+' , '+second)
        elif first not in ['l','r','lr','none'] and second in ['l','r','lr','none']:
            print('Unexpected Arguments: '+first)
        elif second not in ['l','r','lr','none'] and first in ['l','r','lr','none']:
            print('Unexpected Arguments: '+second)       
        sys.exit(0)

    qc_zx = QuantumCircuit(nqubits)
    if (first == 'l' or first == 'lr') and (second == 'r' or second == 'none') :
        qc_zx.h(nqubits[j])
    elif (first == 'r' or first == 'none') and (second == 'r' or second == 'none') :
        qc_zx.h(nqubits[j])

    if cut == 'no' :
        qc_zx.cx(nqubits[i],nqubits[j])
    qc_zx.rz(gamma,nqubits[j])
    qc_zx.cx(nqubits[i],nqubits[j])

    if (first == 'r' or first == 'lr') and (second == 'l' or second == 'none') :
        qc_zx.h(nqubits[j])
    elif (first == 'l' or first == 'none') and (second == 'l' or second == 'none') :
        qc_zx.h(nqubits[j])

    return qc_zx

def GateZY(nqubits,i,j,first='none',second='none',cut='no'):
    if cut not in ['no','yes'] :
        print('Unexpected Arguments: '+cut)
        sys.exit(0)
    if first not in ['l','r','lr','none'] or second not in ['l','r','lr','none'] :
        if first not in ['l','r','lr','none'] and second not in ['l','r','lr','none'] :
            print('Unexpected Arguments: '+first+' , '+second)
        elif first not in ['l','r','lr','none'] and second in ['l','r','lr','none']:
            print('Unexpected Arguments: '+first)
        elif second not in ['l','r','lr','none'] and first in ['l','r','lr','none']:
            print('Unexpected Arguments: '+second)       
        sys.exit(0)

    qc_zy = QuantumCircuit(nqubits)
    if (first == 'l' or first == 'lr') and (second == 'r' or second == 'none') :
        qc_zy.rx(-np.pi/2,nqubits[j])
    elif (first == 'r' or first == 'none') and (second == 'r' or second == 'none') :
        qc_zy.rx(-np.pi/2,nqubits[j])

    if cut == 'no' :
        qc_zy.cx(nqubits[i],nqubits[j])
    qc_zy.rz(gamma,nqubits[j])
    qc_zy.cx(nqubits[i],nqubits[j])

    if (first == 'r' or first == 'lr') and (second == 'l' or second == 'none') :
        qc_zy.rx(np.pi/2,nqubits[j])
    elif (first == 'l' or first == 'none') and (second == 'l' or second == 'none') :
        qc_zy.rx(np.pi/2,nqubits[j])

    return qc_zy

def GateZZ(nqubits,i,j,first='none',second='none',cut='no'):
    if cut not in ['no','yes'] :
        print('Unexpected Arguments: '+cut)
        sys.exit(0)
    if first not in ['l','r','lr','none'] or second not in ['l','r','lr','none'] :
        if first not in ['l','r','lr','none'] and second not in ['l','r','lr','none'] :
            print('Unexpected Arguments: '+first+' , '+second)
        elif first not in ['l','r','lr','none'] and second in ['l','r','lr','none']:
            print('Unexpected Arguments: '+first)
        elif second not in ['l','r','lr','none'] and first in ['l','r','lr','none']:
            print('Unexpected Arguments: '+second)       
        sys.exit(0)

    qc_zz = QuantumCircuit(nqubits)

    if cut == 'no' :
        qc_zz.cx(nqubits[i],nqubits[j])
    qc_zz.rz(gamma,nqubits[j])
    qc_zz.cx(nqubits[i],nqubits[j])

    return qc_zz


# construct QuantumCircuit
def ConstructCircuit(inputcircuit):
    # find number of qubits
    qubitindex = qubitRegex.findall(str(inputcircuit))
    qubits = int(max(list(map(int, qubitindex)))) + 1
    nqubits = QuantumRegister(qubits)
    qc = QuantumCircuit(nqubits)
    for gates in inputcircuit:
        gatecontent = gateRegex.findall(str(gates))
        gateindex = qubitRegex.findall(str(gates))
        idx1 = int(gateindex[0])
        idx2 = int(gateindex[1])

        # add gates to QuantumCircuit
        # Note here we apply a unitary to circuit and then put its inverse, ideally we should get the outcome the same as input gate
        # However, in the presence of noise, we cannot obtain the same input state with probability 1, so when we execute UU^{dagger},
        # the probability of error is multiplied by 2.
        if gatecontent[0] == 'X' and gatecontent[1] == 'X':
            qc.append(GateXX(nqubits,idx1,idx2),[nqubits[i] for i in range(0, qubits)])
        
        elif gatecontent[0] == 'X' and gatecontent[1] == 'Y':
            qc.append(GateXY(nqubits,idx1,idx2),[nqubits[i] for i in range(0, qubits)])
   
        elif gatecontent[0] == 'X' and gatecontent[1] == 'Z':
            qc.append(GateXZ(nqubits,idx1,idx2),[nqubits[i] for i in range(0, qubits)])
     
        elif gatecontent[0] == 'Y' and gatecontent[1] == 'X':
            qc.append(GateYX(nqubits,idx1,idx2),[nqubits[i] for i in range(0, qubits)])
    
        elif gatecontent[0] == 'Y' and gatecontent[1] == 'Y':
            qc.append(GateYY(nqubits,idx1,idx2),[nqubits[i] for i in range(0, qubits)])

        elif gatecontent[0] == 'Y' and gatecontent[1] == 'Z':
            qc.append(GateYZ(nqubits,idx1,idx2),[nqubits[i] for i in range(0, qubits)])
 
        elif gatecontent[0] == 'Z' and gatecontent[1] == 'X':
            qc.append(GateZX(nqubits,idx1,idx2),[nqubits[i] for i in range(0, qubits)])
     
        elif gatecontent[0] == 'Z' and gatecontent[1] == 'Y':
            qc.append(GateZY(nqubits,idx1,idx2),[nqubits[i] for i in range(0, qubits)])
    
        elif gatecontent[0] == 'Z' and gatecontent[1] == 'Z':
            qc.append(GateZZ(nqubits,idx1,idx2),[nqubits[i] for i in range(0, qubits)])
 
        else:
            print('Unexpected Gates')
            sys.exit(0)

    return qc


# l r Identification
def Identify(content):
    if 'l' in content[:] and 'r' in content[:]:
        re = 'lr'
    elif 'l' in content[:] and 'r' not in content[:]:
        re = 'l'
    elif 'l' not in content[:] and 'r' in content[:]:
        re = 'r'
    else:
        re = 'none'
    return re


# Reconstruct optimized circuit with input as ordered circuit
def Optimize(inputcircuit,method='none'):
    # find number of qubits
    qubitindex = qubitRegex.findall(str(inputcircuit))
    qubits = int(max(list(map(int, qubitindex)))) + 1
    nqubits = QuantumRegister(qubits)
    if method not in ['none','DFS','EC'] :
        print('Unexpected Arguments: '+method)
        sys.exit(0)

    gatedict = {}
    for i in range(qubits):
        gatedict[str(i)] = []

    # record gates in every qubit
    for gates in inputcircuit:
        gatecontent = sgateRegex.findall(str(gates))
        gateindex = qubitRegex.findall(str(gates))
        gate1 = gatecontent[0]
        gate2 = gatecontent[1]
        idx1 = int(gateindex[0])
        idx2 = int(gateindex[1])
        for i in range(qubits):
            if idx1 == i:
                gatedict[str(i)].append(gate1)
            if idx2 == i:
                gatedict[str(i)].append(gate2)

    # find adjacent gates that can be reduced and mark them
    for i in range(qubits):
        for (idx1,idx2) in zip( range(len(gatedict[str(i)])-1) , range(1,len(gatedict[str(i)])) ):
            if  gatedict[str(i)][idx1][0] == gatedict[str(i)][idx2][0]:
                gatedict[str(i)][idx1] = gatedict[str(i)][idx1] + 'r'
                gatedict[str(i)][idx2] = gatedict[str(i)][idx2] + 'l'

    ReducedCircuit = []
    countindex = []
    for i in range(qubits):
        countindex.append(0)

    # renew the previous gates with a additional mark
    for gates in inputcircuit:
        gateindex = qubitRegex.findall(str(gates))
        idx1 = int(gateindex[0])
        idx2 = int(gateindex[1])
        for i in range(qubits):
            if idx1 == i:
                gate1 = gatedict[str(i)][countindex[i]]
                countindex[i] = countindex[i] + 1
            if idx2 == i:
                gate2 = gatedict[str(i)][countindex[i]]
                countindex[i] = countindex[i] + 1
        ReducedCircuit.append(gate1+gate2)

    qc = QuantumCircuit(nqubits)
    visited = []
    optimize = True
    for gates in ReducedCircuit:
        gatecontent = gateRegex.findall(str(gates))
        gateindex = qubitRegex.findall(str(gates))
        reducecontent = reduceRegex.findall(str(gates))
        idx1 = int(gateindex[0])
        idx2 = int(gateindex[1])
        if method == 'none' :
            optimize = False
        if method == 'EC' and (idx1 in visited or idx2 in visited) :
            optimize = False
        if method == 'DFS' :
            if idx2 in visited :
                optimize = False
            else:
                optimize = True
        # add gates to QuantumCircuit
        # Note here we apply a unitary to circuit and then put its inverse, ideally we should get the outcome the same as input gate
        # However, in the presence of noise, we cannot obtain the same input state with probability 1, so when we execute UU^{dagger},
        # the probability of error is multiplied by 2.
        if gatecontent[0] == 'X' and gatecontent[1] == 'X':
            r1 = Identify(reducecontent[0])
            r2 = Identify(reducecontent[1])
            if optimize:
                opt = 'yes'
            else:
                opt = 'no'
            qc.append(GateXX(nqubits,idx1,idx2,r1,r2,opt),[nqubits[i] for i in range(0, qubits)])

        elif gatecontent[0] == 'X' and gatecontent[1] == 'Y':
            r1 = Identify(reducecontent[0])
            r2 = Identify(reducecontent[1])
            if optimize:
                opt = 'yes'
            else:
                opt = 'no'
            qc.append(GateXY(nqubits,idx1,idx2,r1,r2,opt),[nqubits[i] for i in range(0, qubits)])
      
        elif gatecontent[0] == 'X' and gatecontent[1] == 'Z':
            r1 = Identify(reducecontent[0])
            r2 = Identify(reducecontent[1])
            if optimize:
                opt = 'yes'
            else:
                opt = 'no'
            qc.append(GateXZ(nqubits,idx1,idx2,r1,r2,opt),[nqubits[i] for i in range(0, qubits)])
      
        elif gatecontent[0] == 'Y' and gatecontent[1] == 'X':
            r1 = Identify(reducecontent[0])
            r2 = Identify(reducecontent[1])
            if optimize:
                opt = 'yes'
            else:
                opt = 'no'
            qc.append(GateYX(nqubits,idx1,idx2,r1,r2,opt),[nqubits[i] for i in range(0, qubits)])
     
        elif gatecontent[0] == 'Y' and gatecontent[1] == 'Y':
            r1 = Identify(reducecontent[0])
            r2 = Identify(reducecontent[1])
            if optimize:
                opt = 'yes'
            else:
                opt = 'no'
            qc.append(GateYY(nqubits,idx1,idx2,r1,r2,opt),[nqubits[i] for i in range(0, qubits)])
    
        elif gatecontent[0] == 'Y' and gatecontent[1] == 'Z':
            r1 = Identify(reducecontent[0])
            r2 = Identify(reducecontent[1])
            if optimize:
                opt = 'yes'
            else:
                opt = 'no'
            qc.append(GateYZ(nqubits,idx1,idx2,r1,r2,opt),[nqubits[i] for i in range(0, qubits)])
   
        elif gatecontent[0] == 'Z' and gatecontent[1] == 'X':
            r1 = Identify(reducecontent[0])
            r2 = Identify(reducecontent[1])
            if optimize:
                opt = 'yes'
            else:
                opt = 'no'
            qc.append(GateZX(nqubits,idx1,idx2,r1,r2,opt),[nqubits[i] for i in range(0, qubits)])
       
        elif gatecontent[0] == 'Z' and gatecontent[1] == 'Y':
            r1 = Identify(reducecontent[0])
            r2 = Identify(reducecontent[1])
            if optimize:
                opt = 'yes'
            else:
                opt = 'no'
            qc.append(GateZY(nqubits,idx1,idx2,r1,r2,opt),[nqubits[i] for i in range(0, qubits)])
      
        elif gatecontent[0] == 'Z' and gatecontent[1] == 'Z':
            r1 = Identify(reducecontent[0])
            r2 = Identify(reducecontent[1])
            if optimize:
                opt = 'yes'
            else:
                opt = 'no'
            qc.append(GateZZ(nqubits,idx1,idx2,r1,r2,opt),[nqubits[i] for i in range(0, qubits)])
     
        else:
            print('Unexpected Gates')
            sys.exit(0)

        visited.append(idx1)
        visited.append(idx2)

    return qc

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

# Qiskit no device method
def Qis(circuit):
    parr = gene_oplist_strlist(circuit)
    qiskitcircuit = {}
    nq = len(parr[0][0])
    length = nq//2 # `length' is a hyperparameter, and can be adjusted for best performance
    a1 = depth_oriented_scheduling(parr, length=length, maxiter=30)
    qc = synthesis_FT.qiskit_synthesis(a1)
    qc2 = transpile(qc, basis_gates=['u3', 'cx'], optimization_level=3)
    qiskitcircuit['gateCount'] = sum(qc2.count_ops().values())
    qiskitcircuit['CNOTCount'] = count_CX(list(qc2.count_ops().keys()),list(qc2.count_ops().values()))
    qiskitcircuit['depth'] = qc2.depth()

    return qiskitcircuit

# PH no device method
def PH(circuit):
    #print('PH passes, Our schedule, Our synthesis', flush=True)
    parr = gene_oplist_strlist(circuit)
    phcircuit = {}
    nq = len(parr[0][0])
    length = nq//2 # `length' is a hyperparameter, and can be adjusted for best performance
    t0 = ctime()
    a1 = depth_oriented_scheduling(parr, length=length, maxiter=30)
    qc = synthesis_FT.block_opt_FT(a1)
    #print('PH, Time costed:', ctime()-t0, flush=True)
    qc1 = transpile(qc, basis_gates=['u3', 'cx'], optimization_level=0)
    t0 = ctime()
    qc2 = transpile(qc, basis_gates=['u3', 'cx'], optimization_level=3)
    #print_qc(qc2)
    phcircuit['gateCount'] = sum(qc2.count_ops().values())
    phcircuit['CNOTCount'] = count_CX(list(qc2.count_ops().keys()),list(qc2.count_ops().values()))
    phcircuit['depth'] = qc2.depth()
    #print('Qiskit L3, Time costed:', ctime()-t0, flush=True)
    '''ansatz = circuit_from_qasm_str(qc1.qasm())
    t0 = ctime()
    FullPeepholeOptimise().apply(ansatz)
    SynthesiseTket().apply(ansatz)
    print(f"CNOT: {ansatz.n_gates_of_type(OpType.CX)}, Single: {ansatz.n_gates-ansatz.n_gates_of_type(OpType.CX)}, Total: {ansatz.n_gates}, Depth: {ansatz.depth()}")
    print("TKet O2:", ctime()-t0)'''

    return phcircuit

# TKet no device method
def TK(circuit):
    #print('TK passes', flush=True)
    parr = gene_oplist_strlist(circuit)
    tkcircuit = {}
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
            qubits, paulis = zip(*term.map.items())
            pbox = PauliExpBox(paulis, coeff * param)
            circ.add_pauliexpbox(pbox, qubits)
    '''ansatz = Circuit(n)
    t0 = ctime()
    add_excitation(ansatz, oplist)
    PauliSimp().apply(ansatz)
    print(f"TK Pauli Simp: {ctime()-t0}")
    t0 = ctime()
    FullPeepholeOptimise().apply(ansatz)
    print(f"CNOT: {ansatz.n_gates_of_type(OpType.CX)}, Single: {ansatz.n_gates-ansatz.n_gates_of_type(OpType.CX)}, Total: {ansatz.n_gates}, Depth: {ansatz.depth()}")
    print("TK O2:", ctime()-t0)'''
    ansatz = Circuit(n)
    add_excitation(ansatz, oplist)
    PauliSimp().apply(ansatz)
    SynthesiseTket().apply(ansatz)
    qstr = circuit_to_qasm_str(ansatz)
    qc = QuantumCircuit.from_qasm_str(qstr)
    t0 = ctime()
    qc = transpile(qc, basis_gates=['cx', 'u3'], optimization_level=3)
    #print_qc(qc)
    tkcircuit['gateCount'] = sum(qc.count_ops().values())
    tkcircuit['CNOTCount'] = count_CX(list(qc.count_ops().keys()),list(qc.count_ops().values()))
    tkcircuit['depth'] = qc.depth()
    #print("Qiskit L3:", ctime()-t0)

    return tkcircuit