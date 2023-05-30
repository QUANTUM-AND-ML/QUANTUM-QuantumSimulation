import matplotlib.pyplot as plt
from CircuitConstruct import *
from CircuitOptimize import *
from HamiltonianGenerator import *
import os
import numpy as np 
from tqdm import tqdm
from TemplateMatching import *


#------------------------------------------------------------------------------------------------
#                     The Function used for benchmarking a specific circuit
#                               (return a dic stores results)
#------------------------------------------------------------------------------------------------
def Benchmarking(circuit):
    #----------------------------------------------------------------------------------------  
    #                    outputCircuit after our Optimization Algorithm
    #     (Firstly using subcircuit Decomposition method to decompose k-local Pailis into 
    #              2-locals, then conduct the greedy based circuit optimization)
    #---------------------------------------------------------------------------------------- 
    local2gate = HamiltonianDecomposition(circuit)
    # store benchmarking results
    results = {}
    results['gateCount'] = {}
    results['depth'] = {}
    results['CNOTCount'] = {}
    #----------------------------------------------------------------------------------------
    # Original Circuit
    #OriginalCircuit = ConstructCircuit(circuit)
    #OriginalCircuit = OriginalCircuit.decompose()
    # OriginalCircuit.draw('mpl')
    # get Ordered Circuit
    DFSCircuit = DFSOrdering(local2gate)
    ECCircuit = ECOrdering(local2gate)
    #---------( The following two ConstructCircuit functions are used for inspection before final output, thet're not necessary )---------
    # optimized circuits of two different strategies 
    # DFSOptimizedCircuit = ConstructCircuit(DFSCircuit)
    # DFSOptimizedCircuit.decompose().draw('mpl')
    # ECOptimizedCircuit = ConstructCircuit(ECCircuit)
    # ECOptimizedCircuit.decompose().draw('mpl')
    # optimized circuit
    #---------------------------------------------------
    Outputcircuit_DFS = Optimize(DFSCircuit,'DFS')
    Outputcircuit_DFS = Outputcircuit_DFS.decompose()
    #Outputcircuit_DFS.draw('mpl')
    Outputcircuit_EC = Optimize(ECCircuit,'EC')
    Outputcircuit_EC = Outputcircuit_EC.decompose()
    #Outputcircuit_EC.draw('mpl')
    #plt.show()
    #---------------------------------------------------------------------------------------- 
    #                   Circuit after optimization of Paulihedral Compiler
    #----------------------------------------------------------------------------------------
    paulihedral_result = PH(circuit)
    #---------------------------------------------------------------------------------------- 
    #                     Circuit after optimization of Qiskit Optimizer
    #---------------------------------------------------------------------------------------- 
    qiskit_result = Qis(circuit)
    #---------------------------------------------------------------------------------------- 
    #              Circuit after optimization of t|ket> quantum simulation programs
    #---------------------------------------------------------------------------------------- 
    tk_result = TK(circuit)
    #---------------------------------------------------------------------------------------- 
    #                    Compute Gate count and Circuit depth
    #---------------------------------------------------------------------------------------- 
    results['gateCount']['DFS'] = sum(Outputcircuit_DFS.count_ops().values())
    results['depth']['DFS'] = Outputcircuit_DFS.depth()
    results['CNOTCount']['DFS'] = count_CX(list(Outputcircuit_DFS.count_ops().keys()),list(Outputcircuit_DFS.count_ops().values()))

    results['gateCount']['EC'] = sum(Outputcircuit_EC.count_ops().values())
    results['depth']['EC'] = Outputcircuit_EC.depth()
    results['CNOTCount']['EC'] = count_CX(list(Outputcircuit_EC.count_ops().keys()),list(Outputcircuit_EC.count_ops().values()))

    results['gateCount']['paulihedral'] = paulihedral_result['gateCount']
    results['depth']['paulihedral'] = paulihedral_result['depth']
    results['CNOTCount']['paulihedral'] = paulihedral_result['CNOTCount']

    results['gateCount']['qiskit'] = qiskit_result['gateCount']
    results['depth']['qiskit'] = qiskit_result['depth']
    results['CNOTCount']['qiskit'] = qiskit_result['CNOTCount']

    results['gateCount']['tket'] = tk_result['gateCount']
    results['depth']['tket'] = tk_result['depth']
    results['CNOTCount']['tket'] = tk_result['CNOTCount']
    

    for optimizer in ['DFS','EC','paulihedral','qiskit','tket']:
        print('   ' + optimizer)
        print('gate count:\t' + str(results['gateCount'][optimizer]))
        print('CNOT count:\t',results['CNOTCount'][optimizer])
        print('depth:\t\t' + str(results['depth'][optimizer]) + '\n')

    return results,local2gate



#-------------------------------------------------------------------------------------------------------
#     ⭐   Benchmarking of our algorithm and other state-of-the-art Compilers and Optimizers    ⭐
#-------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    #----------------------------------------------------------------------------------------
    #              ~ Calculate benchmarking results from corresponding data ~
    #----------------------------------------------------------------------------------------
    #--------------------------------Models Benchmarking-------------------------------------
    
    #----------------------------------------------------------------------------------------
    #                ~ Plot benchmarking results from benchmark models ~
    #----------------------------------------------------------------------------------------
    
    #----------------------------------------------------------------------------------------
    #                          ~ Benchmark a single circuit ~
    #     (output Order of the circuit before and after optimization, output benchmarking
    #              results, draw the circuit by plotting qiskit QuantumCircuit)
    #----------------------------------------------------------------------------------------
    #parr = gene_random_oplist(12)
    Hamiltonian = CH2 # Can be modified to (HF, LiH, H2O, NH2, CH2, NH3, CH4, UCCSD-8, UCCSD-12, UCCSD-16, etc.)
    parr = gene_molecule_oplist(Hamiltonian)
    #heisen_parr = [gene_dot_1d(29, interaction='Z')+gene_dot_1d(29, interaction='X')+gene_dot_1d(29, interaction='Y'), gene_dot_2d(4,5, interaction='Z')+gene_dot_2d(4,5, interaction='Y')+gene_dot_2d(4,5, interaction='X'), gene_dot_3d(1,2,4, interaction='Z')+gene_dot_3d(1,2,4, interaction='Y')+gene_dot_3d(1,2,4, interaction='X')]
    #print(heisen_parr[0])
    #moles = ['LiH', 'BeH2', 'CH4', 'MgH', 'LiCl', 'CO2']
    #parr = gene_FermiHubbard_oplist(5,5)
    #parr = load_oplist('CH4', benchmark='uccsd')

    
    print('~ Comparing to state-of-the-art Quantum Compilers ~')
    # output benchmarking
    result,local2gate = Benchmarking(parr)
    #print(local2gate)
    local2gate_initial = list(set(local2gate))
    #print(local2gate_initial)
    # 字符串赋值
    local2gate_initial_mark = []  
    a = ''  
    for j in range(0,len(local2gate_initial),1):
        for i in local2gate_initial[j]: 
            if str.isdigit(i): 
                a += i  
            else:
                a += " "  
    num_list = a.split(" ")  
    #print("num_list is \n", num_list)
    for i in num_list:  
        try:  
            if int(i) >= 0:
                local2gate_initial_mark.append(int(i)) 
            else:
                pass  
        except:
            pass
    b = []
    for i in range(0, len(local2gate_initial_mark), 2):
        b.append(local2gate_initial_mark[i:i + 2])
    local2gate_initial_mark = b
    #print(local2gate_initial_mark)
    values = [0 for x in range(len(local2gate_initial))]

    #Prioritize the sorting of target bits
    for i in range(0,len(local2gate_initial),1):
        values[i] = abs(local2gate_initial_mark[i][1] - local2gate_initial_mark[i][0])*10000000 +local2gate_initial_mark[i][0]*100
        #for j in range(0,abs(local2gate_initial_mark[i][1] - local2gate_initial_mark[i][0]))
        #if abs(local2gate_initial_mark[i][1] - local2gate_initial_mark[i][0]) == 1 and local2gate_initial_mark[i][0] % 2 != 0:
            #values[i] = values[i] + 50000
        values[i] = values[i] + (local2gate_initial_mark[i][0] % (abs(local2gate_initial_mark[i][1] - local2gate_initial_mark[i][0]) * 2))*10000

        s = ''.join([x for x in local2gate_initial[i] if x.isalpha()])
        if s == 'XY':
            values[i] = values[i] + 1
        elif s == 'XZ':
            values[i] = values[i] + 2
        elif s == 'ZX':
            values[i] = values[i] + 3
        elif s == 'YX':
            values[i] = values[i] + 4
        elif s == 'YY':
            values[i] = values[i] + 5
        elif s == 'YZ':
            values[i] = values[i] + 6
        elif s == 'ZY':
            values[i] = values[i] + 7
        elif s == 'XX':
            values[i] = values[i] + 8
        elif s == 'ZZ':
            values[i] = values[i] + 9
        else:
            pass

    dictionary = dict(zip(local2gate_initial, values))
    local2gate_sort = dict(sorted(dictionary.items(), key=lambda x: x[1], reverse=False)).keys() # Sort the dictionary in ascending order by the size of the values and output the sorted keys
    local2gate_simplification = list(local2gate_sort) # Get the list of Two local doors after sorting
    #print(local2gate_simplification)
    circuit2, circuit1 = simplification(local2gate_simplification)
    print('   ' + 'TM')
    print('gate count:\t' + str(sum(circuit1.count_ops().values())))
    print('CNOT count:\t' + str(count_CX(list(circuit1.count_ops().keys()),list(circuit1.count_ops().values()))))
    print('depth:\t\t' + str(circuit1.depth()) + '\n')
