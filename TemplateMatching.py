from qiskit import *
from qiskit.circuit import Parameter
import numpy as np

theta = Parameter("$\\theta$")

# To obtain the order of qubits and the number of qubits required for the action of the Two local gate list.
def extract_bit(local2gate_simplification):
    '''
    :param local2gate_simplification:输入的列表形式Two local门，example: ['Z0X1', 'Z0Y1', 'Z1X2', 'Y1X2', 'Y1Z2', 'X1X2', 'Y2X3', 'Y2Z3', 'Z2Z3', 'Z0Z2', 'Z1X3', 'Y1Z3', 'X1X3', 'Z1Z3']
    :return types_of_gates:提取作用的Two local门的类型；local2gate_initial_mark:提取的作用比特列表，按照作用对的形式排列；count：作用比特的最大值，确定需要量子比特的数量
    '''
    local2gate_initial_mark = []  # 新建空列表，用以存储提取的数值
    a = ''  # 将空值赋值给a
    for j in range(0, len(local2gate_simplification), 1):
        for i in local2gate_simplification[j]:  # 将字符串进行遍历
            if str.isdigit(i):  # 判断i是否为数字，如果“是”返回True，“不是”返回False
                a += i  # 如果i是数字格式，将i以字符串格式加到a上
            else:
                a += " "  # 如果i不是数字格式，将“ ”（空格）加到a上
    # 数字与数字之间存在许多空格，所以需要对字符串a按''进行分割。
    num_list = a.split(" ")  # 按''进行分割，此时a由字符串格式编程列表
    #print("num_list is \n", num_list)
    for i in num_list:  # 对列表a，进行遍历
        try:  # try 结构体，防止出错后直接退出程序
            if int(i) >= 0:
                local2gate_initial_mark.append(int(i))  # 如果列表a的元素为数字，则赋值给num_list_new
            else:
                pass  # 如果不是数字，pass
        except:
            pass
    #print("local2gate_initial is \n", local2gate_initial_mark)
    count = max(local2gate_initial_mark)
    # 打印出的结果
    #print("len(local2gate_initial):", len(local2gate_initial_mark))  # 作为验证，可以数一下列表元素个数
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
    '''
    :param ansatz:输入的量子线路结构，列表格式; Example: [('rx', 0), ('rz', 1), ('cx', 0, 1)];
    :param order: 所查找量子门所处量子线路结构的列表标记; Example:第一个量子门，所处列表标记为0
    :return:Control_Bit:初始为[-1,-1],单比特门只考虑这个，列表的第一位是前面相邻门的序号，两比特门是控制位对应前面序号，列表第二位是后面相邻门序号，两比特类似,
            Controlled_Bit:初始为[-1，-1]，两比特门受控位对应相邻门的序号，列表第一位是前面相邻列表的序号，列表第二位是后面相邻门序号,
            Quantum_Gate_Marking:初始为[0，0，0，0]，共四位的列表，前两个代表控制比特前后邻居门的类型，1表示单比特门，2表示两比特门的控制位，3表示两比特门的受控位，后两个表示受控比特的前后邻居门的类型;
    '''
    Control_Bit = [-1,-1]
    Controlled_Bit = [-1,-1]
    Quantum_Gate_Marking = [0,0,0,0]

    #判断Rx,Rz前后的邻居编号
    if ansatz[order][0] == 'rx' or ansatz[order][0] == 'rz' or ansatz[order][0] == 'h':
        #寻找前面邻居编号
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
                    print('出现未知量子门')
                    break
        #寻找后面邻居编号
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
                    print('出现未知量子门')
                    break

    #判断Cx前后的邻居编号
    elif ansatz[order][0] == 'cx':
        #寻找前面控制比特的邻居编号
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
        #寻找后面控制比特的邻居编号
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
        #寻找前面受控制比特的邻居编号
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
        #寻找后面受控制比特的邻居编号
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
        print('查找的是未知的量子门')
    return Control_Bit,Controlled_Bit,Quantum_Gate_Marking

def Search_parameters(ansatz, order):
    new_order = 0
    for i in range(0,order,1):
        if ansatz[i][0] == 'rx' or ansatz[i][0] == 'rz':
            new_order = new_order + 1
    return new_order

#删除两个重复的CNOT门
def rule_CX(ansatz):
    '''
    :param ansatz: 输入的量子线路结构，列表格式; Example: [('rx', 0), ('rz', 1), ('cx', 0, 1)];
    :return: 使用规则一后的量子线路结构;
    '''
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
    #print('使用规则CX次数:',len(ansatz) - len(new_ansatz))
    return new_ansatz

#删除两个重复的H门
def rule_H(ansatz):
    '''
    :param ansatz: 输入的量子线路结构，列表格式; Example: [('rx', 0), ('rz', 1), ('cx', 0, 1)];
    :return: 使用规则一后的量子线路结构;
    '''
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
    #print('使用规则H次数:',len(ansatz) - len(new_ansatz))
    return new_ansatz

#删除两个相位相反的RX门
def rule_RX(ansatz,parameters):
    '''
    :param ansatz: 输入的量子线路结构，列表格式; Example: [('rx', 0), ('rz', 1), ('cx', 0, 1)];
    :return: 使用规则一后的量子线路结构;
    '''
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
    #print('使用规则H次数:',len(ansatz) - len(new_ansatz))
    return new_ansatz,new_parameters

#实现CNOT门数目的统计，输入两个列表，对应门的类型和数目
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

    #原始线路图
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

    #变换后线路图
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

    #测试输入
    count = 0
    for i in range(len(gates)):
        if gates[i][0] == 'rx' or gates[i][0] == 'rz':
            count = count + 1
    if count != len(phase):
        print('输入结构和相位不匹配')

    #应用合并相同的门操作并生成线路
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
    #circuit_simplification_full.draw(output='mpl')
    #plt.show()

    return circuit_simplification , circuit_simplification_full
