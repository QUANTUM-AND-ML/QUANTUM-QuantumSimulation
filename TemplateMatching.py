from qiskit import *
from qiskit.circuit import Parameter
import numpy as np

theta = Parameter("$\\theta$")

# To obtain the order of qubits and the number of qubits required for the action of the Two local gate list.
def extract_bit(local2gate_simplification):
    '''
    :param local2gate_simplification:Two local gates，example: ['Z0X1', 'Z0Y1', 'Z1X2', 'Y1X2', 'Y1Z2', 'X1X2', 'Y2X3', 'Y2Z3', 'Z2Z3', 'Z0Z2', 'Z1X3', 'Y1Z3', 'X1X3', 'Z1Z3']
    :return types_of_gates: the role of Two local doors of type
    '''
    local2gate_initial_mark = []  
    a = ''  
    for j in range(0, len(local2gate_simplification), 1):
        for i in local2gate_simplification[j]:  
            if str.isdigit(i):  
                a += i  
            else:
                a += " " 
    num_list = a.split(" ")  
    for i in num_list:  
        try:  
            if int(i) >= 0:
                local2gate_initial_mark.append(int(i))  
            else:
                pass 
        except:
            pass
    #print("local2gate_initial is \n", local2gate_initial_mark)
    count = max(local2gate_initial_mark)
    b = []
    for i in range(0, len(local2gate_initial_mark), 2):
        b.append(local2gate_initial_mark[i:i + 2])
    local2gate_initial_mark = b
    #print(local2gate_initial_mark)
    types_of_gates = []
    for i in range(0,len(local2gate_simplification),1):
        s = ''.join([x for x in local2gate_simplification[i] if x.isalpha()])
        types_of_gates.append(s)

    return types_of_gates, local2gate_initial_mark, count

def neighbor(ansatz,order):
    Control_Bit = [-1,-1]
    Controlled_Bit = [-1,-1]
    Quantum_Gate_Marking = [0,0,0,0]

    # Determine the neighbor numbers before and after Rx,Rz
    if ansatz[order][0] == 'rx' or ansatz[order][0] == 'rz' or ansatz[order][0] == 'h':
        # Find front neighbor number
        if order == 0:
            Control_Bit[0] = -1
            Controlled_Bit[0] = -1
        else:
            for i in range(order-1,-1,-1):
                if ansatz[i][0] == 'rx' or ansatz[i][0] == 'rz' or ansatz[i][0] == 'h':
                    if ansatz[i][1] == ansatz[order][1]:
                        Control_Bit[0] = i
                        Quantum_Gate_Marking[0] = 1
                        break
                elif ansatz[i][0] == 'cx':
                    if ansatz[i][1] == ansatz[order][1]:
                        Control_Bit[0] = i
                        Quantum_Gate_Marking[0] = 2
                        break
                    elif ansatz[i][2] == ansatz[order][1]:
                        Control_Bit[0] = i
                        Quantum_Gate_Marking[0] = 3
                        break
                else:
                    print('Appearance of unknown quantum gates')
                    break
        # Find the neighbor number behind
        if order == len(ansatz):
            Control_Bit[1] = -1
            Controlled_Bit[1] = -1
        else:
            for i in range(order+1,len(ansatz),1):
                if ansatz[i][0] == 'rx' or ansatz[i][0] == 'rz' or ansatz[i][0] == 'h':
                    if ansatz[i][1] == ansatz[order][1]:
                        Control_Bit[1] = i
                        Quantum_Gate_Marking[1] = 1
                        break
                elif ansatz[i][0] == 'cx':
                    if ansatz[i][1] == ansatz[order][1]:
                        Control_Bit[1] = i
                        Quantum_Gate_Marking[1] = 2
                        break
                    elif ansatz[i][2] == ansatz[order][1]:
                        Control_Bit[1] = i
                        Quantum_Gate_Marking[1] = 3
                        break
                else:
                    print('Appearance of unknown quantum gates')
                    break

    # Determine the neighbor numbers before and after Cx
    elif ansatz[order][0] == 'cx':
        # Find the neighbor number of the previous control bit
        if order == 0:
            Control_Bit[0] = -1
            Controlled_Bit[0] = -1
        else:
            for i in range(order-1,-1,-1):
                if ansatz[i][0] == 'rx' or ansatz[i][0] == 'rz' or ansatz[i][0] == 'h':
                    if ansatz[i][1] == ansatz[order][1]:
                        Control_Bit[0] = i
                        Quantum_Gate_Marking[0] = 1
                        break
                if ansatz[i][0] == 'cx':
                    if ansatz[i][1] == ansatz[order][1]:
                        Control_Bit[0] = i
                        Quantum_Gate_Marking[0] = 2
                        break
                    elif ansatz[i][2] == ansatz[order][1]:
                        Control_Bit[0] = i
                        Quantum_Gate_Marking[0] = 3
                        break
        # Find the neighbor number of the control bit that follows
        if order == len(ansatz):
            Control_Bit[1] = -1
            Controlled_Bit[1] = -1
        else:
            for i in range(order+1,len(ansatz),1):
                if ansatz[i][0] == 'rx' or ansatz[i][0] == 'rz' or ansatz[i][0] == 'h':
                    if ansatz[i][1] == ansatz[order][1]:
                        Control_Bit[1] = i
                        Quantum_Gate_Marking[1] = 1
                        break
                if ansatz[i][0] == 'cx':
                    if ansatz[i][1] == ansatz[order][1]:
                        Control_Bit[1] = i
                        Quantum_Gate_Marking[1] = 2
                        break
                    elif ansatz[i][2] == ansatz[order][1]:
                        Control_Bit[1] = i
                        Quantum_Gate_Marking[1] = 3
                        break
        # Find the neighbor number of the previous controlled bit
        if order == 0:
            Control_Bit[0] = -1
            Controlled_Bit[0] = -1
        else:
            for i in range(order-1,-1,-1):
                if ansatz[i][0] == 'rx' or ansatz[i][0] == 'rz' or ansatz[i][0] == 'h':
                    if ansatz[i][1] == ansatz[order][2]:
                        Controlled_Bit[0] = i
                        Quantum_Gate_Marking[2] = 1
                        break
                if ansatz[i][0] == 'cx':
                    if ansatz[i][1] == ansatz[order][2]:
                        Controlled_Bit[0] = i
                        Quantum_Gate_Marking[2] = 2
                        break
                    elif ansatz[i][2] == ansatz[order][2]:
                        Controlled_Bit[0] = i
                        Quantum_Gate_Marking[2] = 3
                        break
        # Find the neighbor number of the controlled bit that follows
        if order == len(ansatz):
            Control_Bit[1] = -1
            Controlled_Bit[1] = -1
        else:
            for i in range(order+1,len(ansatz),1):
                if ansatz[i][0] == 'rx' or ansatz[i][0] == 'rz' or ansatz[i][0] == 'h':
                    if ansatz[i][1] == ansatz[order][2]:
                        Controlled_Bit[1] = i
                        Quantum_Gate_Marking[3] = 1
                        break
                if ansatz[i][0] == 'cx':
                    if ansatz[i][1] == ansatz[order][2]:
                        Controlled_Bit[1] = i
                        Quantum_Gate_Marking[3] = 2
                        break
                    elif ansatz[i][2] == ansatz[order][2]:
                        Controlled_Bit[1] = i
                        Quantum_Gate_Marking[3] = 3
                        break

    else:
        print('The search is for unknown quantum gates')
    return Control_Bit,Controlled_Bit,Quantum_Gate_Marking

def Search_parameters(ansatz, order):
    new_order = 0
    for i in range(0,order,1):
        if ansatz[i][0] == 'rx' or ansatz[i][0] == 'rz':
            new_order = new_order + 1
    return new_order

# Delete two duplicate CNOT doors
def rule_CX(ansatz):
    new_ansatz = list(ansatz)
    Count = len(new_ansatz)
    i = 0
    while Count > 0:
        if new_ansatz[i][0] == 'cx':
            Control_Bit, Controlled_Bit, Quantum_Gate_Marking = neighbor(new_ansatz, i)
            if Control_Bit[1] == Controlled_Bit[1] and Control_Bit[1] != -1:
                new_ansatz.pop(Control_Bit[1])
                new_ansatz.pop(i)
                i = i - 1
                Count = Count - 1
        i = i + 1
        Count = Count - 1
    return new_ansatz

# Delete two duplicate H-gates
def rule_H(ansatz):
    new_ansatz = list(ansatz)
    Count = len(new_ansatz)
    i = 0
    while Count > 0:
        if new_ansatz[i][0] == 'h':
            Control_Bit, Controlled_Bit, Quantum_Gate_Marking = neighbor(new_ansatz, i)
            if Control_Bit[1] != -1:
                if new_ansatz[Control_Bit[1]][0] == 'h':
                    new_ansatz.pop(Control_Bit[1])
                    new_ansatz.pop(i)
                    i = i - 1
                    Count = Count - 1
        i = i + 1
        Count = Count - 1
    return new_ansatz

# Delete two RX gates with opposite phases
def rule_RX(ansatz,parameters):
    new_ansatz = list(ansatz)
    new_parameters = list(parameters)
    Count = len(new_ansatz)
    i = 0
    while Count > 0:
        if new_ansatz[i][0] == 'rx':
            Control_Bit, Controlled_Bit, Quantum_Gate_Marking = neighbor(new_ansatz, i)
            if Control_Bit[1] != -1:
                if new_ansatz[Control_Bit[1]][0] == 'rx' and new_parameters[Search_parameters(new_ansatz,i)] + new_parameters[Search_parameters(new_ansatz,Control_Bit[1])] == 0:
                    new_parameters.pop(Search_parameters(new_ansatz, Control_Bit[1]))
                    new_parameters.pop(Search_parameters(new_ansatz,i))
                    new_ansatz.pop(Control_Bit[1])
                    new_ansatz.pop(i)
                    i = i - 1
                    Count = Count - 1
        i = i + 1
        Count = Count - 1
    return new_ansatz,new_parameters

# To implement the statistics of the number of CNOT doors, enter two lists, corresponding to the type and number of doors
def count_CX(types,counts):
    cnotcount = 0
    for i in range(0,len(types),1):
        if types[i] == 'cx':
            cnotcount = cnotcount + counts[i]
    return cnotcount

def XY(circuit,Affected_bit):
    gates = []
    phase = []
    circuit.h(Affected_bit[0])
    gates.append(['h', Affected_bit[0]])

    circuit.rx(np.pi / 2, Affected_bit[1])
    gates.append(['rx', Affected_bit[1]])
    phase.append(np.pi / 2)

    circuit.cx(Affected_bit[0], Affected_bit[1])
    gates.append(['cx', Affected_bit[0], Affected_bit[1]])

    circuit.rz(theta, Affected_bit[1])
    gates.append(['rz', Affected_bit[1]])
    phase.append(theta)

    circuit.cx(Affected_bit[0], Affected_bit[1])
    gates.append(['cx', Affected_bit[0], Affected_bit[1]])

    circuit.h(Affected_bit[0])
    gates.append(['h', Affected_bit[0]])

    circuit.rx(-np.pi / 2, Affected_bit[1])
    gates.append(['rx', Affected_bit[1]])
    phase.append(-np.pi / 2)

    return gates,phase

def XY1(circuit,Affected_bit):
    gates = []
    phase = []
    circuit.h(Affected_bit[0])
    gates.append(['h', Affected_bit[0]])

    circuit.cx(Affected_bit[0], Affected_bit[1])
    gates.append(['cx', Affected_bit[0],Affected_bit[1]])

    circuit.rx(np.pi / 2, Affected_bit[1])
    gates.append(['rx', Affected_bit[1]])
    phase.append(np.pi / 2)

    circuit.rz(theta, Affected_bit[1])
    gates.append(['rz', Affected_bit[1]])
    phase.append(theta)

    circuit.rx(-np.pi / 2, Affected_bit[1])
    gates.append(['rx', Affected_bit[1]])
    phase.append(-np.pi / 2)

    circuit.cx(Affected_bit[0], Affected_bit[1])
    gates.append(['cx', Affected_bit[0],Affected_bit[1]])

    circuit.h(Affected_bit[0])
    gates.append(['h', Affected_bit[0]])
    return gates,phase

def XY2(circuit,Affected_bit):
    gates = []
    phase = []
    circuit.h(Affected_bit[1])
    gates.append(['h', Affected_bit[1]])

    circuit.cx(Affected_bit[1], Affected_bit[0])
    gates.append(['cx', Affected_bit[1], Affected_bit[0]])

    circuit.rz(np.pi / 2, Affected_bit[1])
    gates.append(['rz', Affected_bit[1]])
    phase.append(np.pi / 2)

    circuit.rx(theta, Affected_bit[1])
    gates.append(['rx', Affected_bit[1]])
    phase.append(theta)

    circuit.rz(-np.pi / 2, Affected_bit[1])
    gates.append(['rz', Affected_bit[1]])
    phase.append(-np.pi / 2)

    circuit.cx(Affected_bit[1], Affected_bit[0])
    gates.append(['cx', Affected_bit[1], Affected_bit[0]])

    circuit.h(Affected_bit[1])
    gates.append(['h', Affected_bit[1]])
    return gates,phase

def XZ(circuit,Affected_bit):
    gates = []
    phase = []
    circuit.h(Affected_bit[0])
    gates.append(['h', Affected_bit[0]])

    circuit.cx(Affected_bit[0], Affected_bit[1])
    gates.append(['cx', Affected_bit[0], Affected_bit[1]])

    circuit.rz(theta, Affected_bit[1])
    gates.append(['rz', Affected_bit[1]])
    phase.append(theta)

    circuit.cx(Affected_bit[0], Affected_bit[1])
    gates.append(['cx', Affected_bit[0], Affected_bit[1]])

    circuit.h(Affected_bit[0])
    gates.append(['h', Affected_bit[0]])
    return gates,phase

def XZ1(circuit,Affected_bit):
    gates = []
    phase = []
    circuit.h(Affected_bit[1])
    gates.append(['h', Affected_bit[1]])

    circuit.cx(Affected_bit[1], Affected_bit[0])
    gates.append(['cx', Affected_bit[1], Affected_bit[0]])

    circuit.rx(theta, Affected_bit[1])
    gates.append(['rx', Affected_bit[1]])
    phase.append(theta)

    circuit.cx(Affected_bit[1], Affected_bit[0])
    gates.append(['cx', Affected_bit[1], Affected_bit[0]])

    circuit.h(Affected_bit[1])
    gates.append(['h', Affected_bit[1]])
    return gates,phase

def XZ2(circuit,Affected_bit):
    gates = []
    phase = []
    circuit.h(Affected_bit[1])
    gates.append(['h', Affected_bit[1]])

    circuit.cx(Affected_bit[0], Affected_bit[1])
    gates.append(['cx', Affected_bit[0], Affected_bit[1]])

    circuit.rx(theta, Affected_bit[0])
    gates.append(['rx', Affected_bit[0]])
    phase.append(theta)

    circuit.cx(Affected_bit[0], Affected_bit[1])
    gates.append(['cx', Affected_bit[0], Affected_bit[1]])

    circuit.h(Affected_bit[1])
    gates.append(['h', Affected_bit[1]])
    return gates,phase

def ZX(circuit,Affected_bit):
    gates = []
    phase = []
    circuit.h(Affected_bit[1])
    gates.append(['h', Affected_bit[1]])

    circuit.cx(Affected_bit[0], Affected_bit[1])
    gates.append(['cx', Affected_bit[0], Affected_bit[1]])

    circuit.rz(theta, Affected_bit[1])
    gates.append(['rz', Affected_bit[1]])
    phase.append(theta)

    circuit.cx(Affected_bit[0], Affected_bit[1])
    gates.append(['cx', Affected_bit[0], Affected_bit[1]])

    circuit.h(Affected_bit[1])
    gates.append(['h', Affected_bit[1]])
    return gates,phase

def ZX1(circuit,Affected_bit):
    gates = []
    phase = []
    circuit.h(Affected_bit[0])
    gates.append(['h', Affected_bit[0]])

    circuit.cx(Affected_bit[0], Affected_bit[1])
    gates.append(['cx', Affected_bit[0], Affected_bit[1]])

    circuit.rx(theta, Affected_bit[0])
    gates.append(['rx', Affected_bit[0]])
    phase.append(theta)

    circuit.cx(Affected_bit[0], Affected_bit[1])
    gates.append(['cx', Affected_bit[0], Affected_bit[1]])

    circuit.h(Affected_bit[0])
    gates.append(['h', Affected_bit[0]])
    return gates,phase

def ZX2(circuit,Affected_bit):
    gates = []
    phase = []
    circuit.h(Affected_bit[0])
    gates.append(['h', Affected_bit[0]])

    circuit.cx(Affected_bit[1], Affected_bit[0])
    gates.append(['cx', Affected_bit[1], Affected_bit[0]])

    circuit.rx(theta, Affected_bit[1])
    gates.append(['rx', Affected_bit[1]])
    phase.append(theta)

    circuit.cx(Affected_bit[1], Affected_bit[0])
    gates.append(['cx', Affected_bit[1], Affected_bit[0]])

    circuit.h(Affected_bit[0])
    gates.append(['h', Affected_bit[0]])
    return gates,phase

def YX(circuit,Affected_bit):
    gates = []
    phase = []
    circuit.rx(np.pi / 2, Affected_bit[0])
    gates.append(['rx', Affected_bit[0]])
    phase.append(np.pi / 2)

    circuit.h(Affected_bit[1])
    gates.append(['h', Affected_bit[1]])

    circuit.cx(Affected_bit[0], Affected_bit[1])
    gates.append(['cx', Affected_bit[0], Affected_bit[1]])

    circuit.rz(theta, Affected_bit[1])
    gates.append(['rz', Affected_bit[1]])
    phase.append(theta)

    circuit.cx(Affected_bit[0], Affected_bit[1])
    gates.append(['cx', Affected_bit[0], Affected_bit[1]])

    circuit.rx(-np.pi / 2, Affected_bit[0])
    gates.append(['rx', Affected_bit[0]])
    phase.append(-np.pi / 2)

    circuit.h(Affected_bit[1])
    gates.append(['h', Affected_bit[1]])
    return gates,phase

def YX1(circuit,Affected_bit):
    gates = []
    phase = []
    circuit.h(Affected_bit[1])
    gates.append(['h', Affected_bit[1]])

    circuit.cx(Affected_bit[1], Affected_bit[0])
    gates.append(['cx', Affected_bit[1], Affected_bit[0]])

    circuit.rx(np.pi / 2, Affected_bit[0])
    gates.append(['rx', Affected_bit[0]])
    phase.append(np.pi / 2)

    circuit.rz(theta, Affected_bit[0])
    gates.append(['rz', Affected_bit[0]])

    circuit.rx(-np.pi / 2, Affected_bit[0])
    gates.append(['rx', Affected_bit[0]])
    phase.append(-np.pi / 2)

    circuit.cx(Affected_bit[1], Affected_bit[0])
    gates.append(['cx', Affected_bit[1], Affected_bit[0]])

    circuit.h(Affected_bit[1])
    gates.append(['h', Affected_bit[1]])
    return gates,phase

def YX2(circuit,Affected_bit):
    gates = []
    phase = []
    circuit.h(Affected_bit[0])
    gates.append(['h', Affected_bit[0]])

    circuit.cx(Affected_bit[0], Affected_bit[1])
    gates.append(['cx', Affected_bit[0], Affected_bit[1]])

    circuit.rz(np.pi / 2, Affected_bit[0])
    gates.append(['rz', Affected_bit[0]])
    phase.append(np.pi / 2)

    circuit.rx(theta, Affected_bit[0])
    gates.append(['rx', Affected_bit[0]])
    phase.append(theta)

    circuit.rz(-np.pi / 2, Affected_bit[0])
    gates.append(['rz', Affected_bit[0]])
    phase.append(-np.pi / 2)

    circuit.cx(Affected_bit[0], Affected_bit[1])
    gates.append(['cx', Affected_bit[0], Affected_bit[1]])

    circuit.h(Affected_bit[0])
    gates.append(['h', Affected_bit[0]])
    return gates,phase

def YY(circuit,Affected_bit):
    gates = []
    phase = []
    circuit.rx(np.pi / 2, Affected_bit[0])
    gates.append(['rx', Affected_bit[0]])
    phase.append(np.pi / 2)

    circuit.rx(np.pi / 2, Affected_bit[1])
    gates.append(['rx', Affected_bit[1]])
    phase.append(np.pi / 2)

    circuit.cx(Affected_bit[0], Affected_bit[1])
    gates.append(['cx', Affected_bit[0], Affected_bit[1]])

    circuit.rz(theta, Affected_bit[1])
    gates.append(['rz', Affected_bit[1]])
    phase.append(theta)

    circuit.cx(Affected_bit[0], Affected_bit[1])
    gates.append(['cx', Affected_bit[0], Affected_bit[1]])

    circuit.rx(-np.pi / 2, Affected_bit[0])
    gates.append(['rx', Affected_bit[0]])
    phase.append(-np.pi / 2)

    circuit.rx(-np.pi / 2, Affected_bit[1])
    gates.append(['rx', Affected_bit[1]])
    phase.append(-np.pi / 2)
    return gates,phase

def YY1(circuit,Affected_bit):
    gates = []
    phase = []
    circuit.rx(np.pi / 2, Affected_bit[0])
    gates.append(['rx', Affected_bit[0]])
    phase.append(np.pi / 2)

    circuit.cx(Affected_bit[0], Affected_bit[1])
    gates.append(['cx', Affected_bit[0], Affected_bit[1]])

    circuit.rx(np.pi / 2, Affected_bit[1])
    gates.append(['rx', Affected_bit[1]])
    phase.append(np.pi / 2)

    circuit.rz(theta, Affected_bit[1])
    gates.append(['rz', Affected_bit[1]])
    phase.append(theta)

    circuit.rx(-np.pi / 2, Affected_bit[1])
    gates.append(['rx', Affected_bit[1]])
    phase.append(-np.pi / 2)

    circuit.cx(Affected_bit[0], Affected_bit[1])
    gates.append(['cx', Affected_bit[0], Affected_bit[1]])

    circuit.rx(-np.pi / 2, Affected_bit[0])
    gates.append(['rx', Affected_bit[0]])
    phase.append(-np.pi / 2)
    return gates,phase

def YY2(circuit,Affected_bit):
    gates = []
    phase = []
    circuit.rx(np.pi / 2, Affected_bit[1])
    gates.append(['rx', Affected_bit[1]])
    phase.append(np.pi / 2)

    circuit.cx(Affected_bit[1], Affected_bit[0])
    gates.append(['cx', Affected_bit[1], Affected_bit[0]])

    circuit.rx(np.pi / 2, Affected_bit[0])
    gates.append(['rx', Affected_bit[0]])
    phase.append(np.pi / 2)

    circuit.rz(theta, Affected_bit[0])
    gates.append(['rz', Affected_bit[0]])
    phase.append(theta)

    circuit.rx(-np.pi / 2, Affected_bit[0])
    gates.append(['rx', Affected_bit[0]])
    phase.append(-np.pi / 2)

    circuit.cx(Affected_bit[1], Affected_bit[0])
    gates.append(['cx', Affected_bit[1], Affected_bit[0]])

    circuit.rx(-np.pi / 2, Affected_bit[1])
    gates.append(['rx', Affected_bit[1]])
    phase.append(-np.pi / 2)
    return gates,phase

def YZ(circuit,Affected_bit):
    gates = []
    phase = []
    circuit.rx(np.pi / 2, Affected_bit[0])
    gates.append(['rx', Affected_bit[0]])
    phase.append(np.pi / 2)

    circuit.cx(Affected_bit[0], Affected_bit[1])
    gates.append(['cx', Affected_bit[0], Affected_bit[1]])

    circuit.rz(theta, Affected_bit[1])
    gates.append(['rz', Affected_bit[1]])
    phase.append(theta)

    circuit.cx(Affected_bit[0], Affected_bit[1])
    gates.append(['cx', Affected_bit[0], Affected_bit[1]])

    circuit.rx(-np.pi / 2, Affected_bit[0])
    gates.append(['rx', Affected_bit[0]])
    phase.append(-np.pi / 2)
    return gates,phase

def YZ1(circuit,Affected_bit):
    gates = []
    phase = []
    circuit.cx(Affected_bit[1], Affected_bit[0])
    gates.append(['cx', Affected_bit[1], Affected_bit[0]])

    circuit.rx(np.pi / 2, Affected_bit[0])
    gates.append(['rx', Affected_bit[0]])
    phase.append(np.pi / 2)

    circuit.rz(theta, Affected_bit[0])
    gates.append(['rz', Affected_bit[0]])
    phase.append(theta)

    circuit.rx(-np.pi / 2, Affected_bit[0])
    gates.append(['rx', Affected_bit[0]])
    phase.append(-np.pi / 2)

    circuit.cx(Affected_bit[1], Affected_bit[0])
    gates.append(['cx', Affected_bit[1], Affected_bit[0]])
    return gates,phase

def ZY(circuit,Affected_bit):
    gates = []
    phase = []
    circuit.rx(np.pi / 2, Affected_bit[1])
    gates.append(['rx', Affected_bit[1]])
    phase.append(np.pi / 2)

    circuit.cx(Affected_bit[0], Affected_bit[1])
    gates.append(['cx', Affected_bit[0], Affected_bit[1]])

    circuit.rz(theta, Affected_bit[1])
    gates.append(['rz', Affected_bit[1]])
    phase.append(theta)

    circuit.cx(Affected_bit[0], Affected_bit[1])
    gates.append(['cx', Affected_bit[0], Affected_bit[1]])

    circuit.rx(-np.pi / 2, Affected_bit[1])
    gates.append(['rx', Affected_bit[1]])
    phase.append(-np.pi / 2)
    return gates,phase

def ZY1(circuit,Affected_bit):
    gates = []
    phase = []
    circuit.cx(Affected_bit[0], Affected_bit[1])
    gates.append(['cx', Affected_bit[0], Affected_bit[1]])

    circuit.rx(np.pi / 2, Affected_bit[1])
    gates.append(['rx', Affected_bit[1]])
    phase.append(np.pi / 2)

    circuit.rz(theta, Affected_bit[1])
    gates.append(['rz', Affected_bit[1]])
    phase.append(theta)

    circuit.rx(-np.pi / 2, Affected_bit[1])
    gates.append(['rx', Affected_bit[1]])
    phase.append(-np.pi / 2)

    circuit.cx(Affected_bit[0], Affected_bit[1])
    gates.append(['cx', Affected_bit[0], Affected_bit[1]])
    return gates,phase

def XX(circuit,Affected_bit):
    gates = []
    phase = []
    circuit.h(Affected_bit[0])
    gates.append(['h', Affected_bit[0]])

    circuit.h(Affected_bit[1])
    gates.append(['h', Affected_bit[1]])

    circuit.cx(Affected_bit[0], Affected_bit[1])
    gates.append(['cx', Affected_bit[0], Affected_bit[1]])

    circuit.rz(theta, Affected_bit[1])
    gates.append(['rz', Affected_bit[1]])
    phase.append(theta)

    circuit.cx(Affected_bit[0], Affected_bit[1])
    gates.append(['cx', Affected_bit[0], Affected_bit[1]])

    circuit.h(Affected_bit[0])
    gates.append(['h', Affected_bit[0]])

    circuit.h(Affected_bit[1])
    gates.append(['h', Affected_bit[1]])
    return gates,phase

def XX1(circuit,Affected_bit):
    gates = []
    phase = []
    circuit.cx(Affected_bit[0], Affected_bit[1])
    gates.append(['cx', Affected_bit[0], Affected_bit[1]])

    circuit.rx(theta, Affected_bit[0])
    gates.append(['rx', Affected_bit[0]])
    phase.append(theta)

    circuit.cx(Affected_bit[0], Affected_bit[1])
    gates.append(['cx', Affected_bit[0], Affected_bit[1]])
    return gates,phase

def XX2(circuit,Affected_bit):
    gates = []
    phase = []
    circuit.cx(Affected_bit[1], Affected_bit[0])
    gates.append(['cx', Affected_bit[1], Affected_bit[0]])

    circuit.rx(theta, Affected_bit[1])
    gates.append(['rx', Affected_bit[1]])
    phase.append(theta)

    circuit.cx(Affected_bit[1], Affected_bit[0])
    gates.append(['cx', Affected_bit[1], Affected_bit[0]])
    return gates, phase

def ZZ(circuit,Affected_bit):
    gates = []
    phase = []
    circuit.cx(Affected_bit[0], Affected_bit[1])
    gates.append(['cx', Affected_bit[0], Affected_bit[1]])

    circuit.rz(theta, Affected_bit[1])
    gates.append(['rz', Affected_bit[1]])
    phase.append(theta)

    circuit.cx(Affected_bit[0], Affected_bit[1])
    gates.append(['cx', Affected_bit[0], Affected_bit[1]])
    return gates,phase

def ZZ1(circuit,Affected_bit):
    gates = []
    phase = []
    circuit.cx(Affected_bit[1], Affected_bit[0])
    gates.append(['cx', Affected_bit[1], Affected_bit[0]])

    circuit.rz(theta, Affected_bit[0])
    gates.append(['rz', Affected_bit[0]])
    phase.append(theta)

    circuit.cx(Affected_bit[1], Affected_bit[0])
    gates.append(['cx', Affected_bit[1], Affected_bit[0]])
    return gates,phase

def simplification(local2gate_simplification):
    theta = Parameter("$\\theta$")
    gates = []
    phase = []
    types_of_gates, local2gate_simplification_bit, count = extract_bit(local2gate_simplification)
    circuit = QuantumCircuit(count + 1)
    circuit_simplification = QuantumCircuit(count + 1)
    circuit_simplification_full = QuantumCircuit(count + 1)

    # Original circuit
    for i in range(0,len(types_of_gates),1):
        if types_of_gates[i] == 'XY':
            a1,b1 = XY(circuit,local2gate_simplification_bit[i])
        elif types_of_gates[i] == 'XZ':
            a2,b2 = XZ(circuit,local2gate_simplification_bit[i])
        elif types_of_gates[i] == 'ZX':
            a3,b3 = ZX(circuit,local2gate_simplification_bit[i])
        elif types_of_gates[i] == 'YX':
            a4,b4 = YX(circuit,local2gate_simplification_bit[i])
        elif types_of_gates[i] == 'YY':
            a5,b5 = YY(circuit,local2gate_simplification_bit[i])
        elif types_of_gates[i] == 'YZ':
            a6,b6 = YZ(circuit,local2gate_simplification_bit[i])
        elif types_of_gates[i] == 'ZY':
            a7,b7 = ZY(circuit,local2gate_simplification_bit[i])
        elif types_of_gates[i] == 'XX':
            a8,b8 = XX(circuit,local2gate_simplification_bit[i])
        elif types_of_gates[i] == 'ZZ':
            a9,b9 = ZZ(circuit,local2gate_simplification_bit[i])
        else:
            pass
        circuit.barrier()

    # Transformed circuit
    j = 0
    yzxx = []
    while j < len(types_of_gates):
        if types_of_gates[j] == 'XY':
            a10,b10 = XY1(circuit_simplification,local2gate_simplification_bit[j])
            gates += a10
            phase += b10
        elif types_of_gates[j] == 'XZ':
            a11,b11 = XZ(circuit_simplification,local2gate_simplification_bit[j])
            gates += a11
            phase += b11
        elif types_of_gates[j] == 'ZX':
            a12,b12 = ZX1(circuit_simplification,local2gate_simplification_bit[j])
            gates += a12
            phase += b12
        elif types_of_gates[j] == 'YX':
            a13,b13 = YX2(circuit_simplification,local2gate_simplification_bit[j])
            gates += a13
            phase += b13
        elif types_of_gates[j] == 'YY':
            a14,b14 = YY1(circuit_simplification,local2gate_simplification_bit[j])
            gates += a14
            phase += b14
        else:
            if j >= 1 and j < len(types_of_gates) - 1:
                if types_of_gates[j] == 'YZ':
                    if local2gate_simplification_bit[j] == local2gate_simplification_bit[j - 1] and local2gate_simplification_bit[j] == local2gate_simplification_bit[j + 1] and types_of_gates[j] == 'YZ' and types_of_gates[j - 1] != 'YY' and types_of_gates[j +1] != 'ZY' and (types_of_gates[j + 1] == 'XX' or types_of_gates[j + 1] == 'ZZ'):
                        a15,b15 = YZ1(circuit_simplification, local2gate_simplification_bit[j])
                        gates += a15
                        phase += b15
                    else:
                        a16,b16 = YZ(circuit_simplification, local2gate_simplification_bit[j])
                        gates += a16
                        phase += b16
                elif types_of_gates[j] == 'ZY':
                    a17,b17 = ZY1(circuit_simplification, local2gate_simplification_bit[j])
                    gates += a17
                    phase += b17
                elif types_of_gates[j] == 'XX':
                    if types_of_gates[j - 1] == 'YZ' and local2gate_simplification_bit[j] == local2gate_simplification_bit[j - 1]:
                        a18,b18 = XX2(circuit_simplification, local2gate_simplification_bit[j])
                        gates += a18
                        phase += b18
                        yzxx = local2gate_simplification_bit[j]
                        #print(yzxx)
                    else:
                        a19,b19 = XX1(circuit_simplification,local2gate_simplification_bit[j])
                        gates += a19
                        phase += b19
                elif types_of_gates[j] == 'ZZ':
                    if types_of_gates[j - 1] == 'YZ' and local2gate_simplification_bit[j] == local2gate_simplification_bit[j - 1]:
                        a20,b20 = ZZ1(circuit_simplification, local2gate_simplification_bit[j])
                        gates += a20
                        phase += b20
                    elif yzxx == local2gate_simplification_bit[j]:
                        a21,b21 = ZZ1(circuit_simplification, local2gate_simplification_bit[j])
                        gates += a21
                        phase += b21
                    else:
                        a22,b22 = ZZ(circuit_simplification, local2gate_simplification_bit[j])
                        gates += a22
                        phase += b22
                else:
                    pass
            else:
                if types_of_gates[j] == 'YZ':
                    a23,b23 = YZ(circuit_simplification, local2gate_simplification_bit[j])
                    gates += a23
                    phase += b23
                elif types_of_gates[j] == 'ZY':
                    a24,b24 = ZY1(circuit_simplification, local2gate_simplification_bit[j])
                    gates += a24
                    phase += b24
                elif types_of_gates[j] == 'XX':
                    if types_of_gates[j - 1] == 'YZ' and local2gate_simplification_bit[j] == \
                            local2gate_simplification_bit[j - 1]:
                        a25,b25 = XX2(circuit_simplification, local2gate_simplification_bit[j])
                        gates += a25
                        phase += b25
                    else:
                        a26,b26 = XX1(circuit_simplification, local2gate_simplification_bit[j])
                        gates += a26
                        phase += b26
                elif types_of_gates[j] == 'ZZ':
                    if types_of_gates[j - 1] == 'YZ' and local2gate_simplification_bit[j] == local2gate_simplification_bit[j - 1]:
                        a27,b27 = ZZ1(circuit_simplification, local2gate_simplification_bit[j])
                        gates += a27
                        phase += b27
                    elif yzxx == local2gate_simplification_bit[j]:
                        a28,b28 = ZZ1(circuit_simplification, local2gate_simplification_bit[j])
                        gates += a28
                        phase += b28
                    else:
                        a29,b29 = ZZ(circuit_simplification, local2gate_simplification_bit[j])
                        gates += a29
                        phase += b29
                else:
                    pass
        #circuit_simplification.barrier()
        j = j + 1
    #circuit.draw(output='mpl')
    #plt.show

    #circuit_simplification.draw(output='mpl')
    #plt.show()

    # Input for testing
    count = 0
    for i in range(len(gates)):
        if gates[i][0] == 'rx' or gates[i][0] == 'rz':
            count = count + 1
    if count != len(phase):
        print('Input structure and phase mismatch')

    # Application merges identical door operations and generates lines
    for i in range(3):
        gates = rule_H(gates)
        gates,phase = rule_RX(gates,phase)
        gates = rule_CX(gates)
    for i in range(0,len(gates),1):
        if gates[i][0] == 'h':
            circuit_simplification_full.h(gates[i][1])
        elif gates[i][0] == 'rx':
            circuit_simplification_full.rx(phase[Search_parameters(gates,i)],gates[i][1])
        elif gates[i][0] == 'rz':
            circuit_simplification_full.rz(phase[Search_parameters(gates,i)],gates[i][1])
        elif gates[i][0] == 'cx':
            circuit_simplification_full.cx(gates[i][1],gates[i][2])
        else:
            pass
    
    return circuit_simplification , circuit_simplification_full
