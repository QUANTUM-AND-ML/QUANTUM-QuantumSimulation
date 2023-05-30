from collections import Counter
from CircuitConstruct import *


# find index of first gate
def takeIndex(elem):
    gatecontent = gateRegex.findall(str(elem))
    gateindex = qubitRegex.findall(str(elem))
    temp = gateindex[0] + gatecontent[0] + gateindex[1] + gatecontent[1]
    return temp

def qubitIndex(elem):
    gatecontent = gateRegex.findall(str(elem))
    gateindex = qubitRegex.findall(str(elem))
    temp = gateindex[0] + gatecontent[0]
    return temp

# Parallelize all possible operations in DFSOrdering
def DFSParallelization(Circuit,layerindex):
    oddterm = []
    eventerm = []
    #print(layerindex)
    #print(Circuit)
    if len(layerindex) !=1 :   
        for i in range(layerindex[-2],layerindex[-1],2):
            oddterm.append(Circuit[i])
        for i in range(layerindex[-2]+1,layerindex[-1],2):
            eventerm.append(Circuit[i])
        Circuit[layerindex[-2]:layerindex[-1]] = oddterm + eventerm
    #print(oddterm + eventerm)
    #print(Circuit)

# Parallelize all possible operations in ECOrdering
def ECParallelization(Circuit,layerindex):
    oddterm = []
    eventerm = []
    if len(layerindex) ==1 :
        for i in range(0,layerindex[-1],2):
            oddterm.append(Circuit[i])
        for i in range(1,layerindex[-1],2):
            eventerm.append(Circuit[i])
        Circuit[0:layerindex[-1]] = oddterm + eventerm
    else:    
        for i in range(layerindex[-2],layerindex[-1],2):
            oddterm.append(Circuit[i])
        for i in range(layerindex[-2]+1,layerindex[-1],2):
            eventerm.append(Circuit[i])
        Circuit[layerindex[-2]:layerindex[-1]] = oddterm + eventerm
    


#----------------------------------------------------------------------------------------
#          Rearrange circuit by Greedy Algorithm to reducing single qubit gates
#----------------------------------------------------------------------------------------
                      
#---------------------------------------------------------------------------------------- 
#                      Edge Coloring ordering algorithm Strategy
#         (Edge Coloring for all layers to make it as much parallel as possible,
#                concatenation-gate ordering in one layer takes priority)
def ECOrdering(inputcircuit):
    gate1pool = []
    gate2pool = []
    lastlayergate1pool = []
    lastlayergate2pool = []
    llgate1index = []
    llgate2index = []
    OrderedCircuit = []
    SortedOrderedCircuit = []
    lonelygate = []
    rest = []
    last = ''
    top = True
    # layerindex stores the indexes in OrderedCircuit from which new layers start, every element in layerindex list labels
    # the index ( in OrderedCircuit ) of the first gate in one layer 
    layerindex = []
    # sort all local2gates and group them into gate1pool and gate2pool
    for gates in inputcircuit:
            gatecontent = sgateRegex.findall(str(gates))
            gate1pool.append(gatecontent[0])
            gate2pool.append(gatecontent[1])

    totalgate = gate1pool + gate2pool
    lonelygatecount = dict(Counter(totalgate))

    # find lonely gates that cannot be reduced
    for gates in inputcircuit:
            gatecontent = sgateRegex.findall(str(gates))
            lonelygate1 = gatecontent[0]
            lonelygate2 = gatecontent[1]
            if lonelygatecount[lonelygate1] == 1 and lonelygatecount[lonelygate2] == 1 :
                lonelygate.append(gates)
            
    Sortedlocal2gate_R = sorted(list(set(inputcircuit) - set(lonelygate)),key=takeIndex)

    # find the superior order that maximize the amount of reduced single qubit gates
    while SortedOrderedCircuit != Sortedlocal2gate_R:
        for gates in Sortedlocal2gate_R:
            gatecontent = sgateRegex.findall(str(gates))
            gate1 = gatecontent[0]
            gate2 = gatecontent[1]

            # rest gates not been optimized
            rest = list(set(Sortedlocal2gate_R)-set(OrderedCircuit))
            rest.sort(key=takeIndex)
            restgate1pool = []
            restgate2pool = []
            lastlayercombo = []
            currentgate2opt = []
            
            
            for restgate in rest:
                restcontent = sgateRegex.findall(str(restgate))
                restgate1pool.append(restcontent[0])
                restgate2pool.append(restcontent[1])
                if restcontent[0] == gate1:
                    currentgate2opt.append(restcontent[1])

            # find qubit index of gate2 in last layer
            for lastlayergate1, lastlayergate2 in zip(lastlayergate1pool,lastlayergate2pool) :
                llgate1content = qubitRegex.findall(str(lastlayergate1))
                llgate2content = qubitRegex.findall(str(lastlayergate2))
                llgate1index.append(llgate1content[0])
                llgate2index.append(llgate2content[0])
            
            llgate1index = list(set(llgate1index))
            llgate2index = list(set(llgate2index))
            llgate1index.sort()
            llgate2index.sort()
            lastlayergate1pool.sort(key=qubitIndex)
            lastlayergate2pool.sort(key=qubitIndex)
            lastlayergatespool = list( set(lastlayergate1pool) | set(lastlayergate2pool) )
            lastlayergatespool.sort(key=qubitIndex)

            # find all possible gate combinations that can be optimized
            for lastlayergate1 in lastlayergatespool:
                for lastlayergate2 in lastlayergatespool:
                    lastlayercombo.append(lastlayergate1+lastlayergate2)

            if not OrderedCircuit:
                last = restgate1pool[0]

            if gates in OrderedCircuit:
                continue    

            judge_rlc = bool(set(rest) & set(lastlayercombo))
            judge_c2l2 = bool(set(currentgate2opt) & set(lastlayergate2pool))
            judge_c2mr2 = bool(set(currentgate2opt) - set(restgate2pool))
            judge_c2ar1 = bool(set(currentgate2opt) & set(restgate1pool))
            if top:
                # choose gate that can reduce gate of last layer twice
                if gate1 in lastlayergatespool and gate2 in lastlayergatespool :
                    OrderedCircuit.append(gates)
                    top = False
                    rest = list(set(Sortedlocal2gate_R)-set(OrderedCircuit))
                    rest.sort(key=takeIndex)
                    last = gate2
                    # if it can't be reduced anymore, switch to next layer
                    if gate2 not in restgate1pool:   
                        restgate1pool = []       
                        for restgate in rest:
                            restcontent = sgateRegex.findall(str(restgate))
                            restgate1pool.append(restcontent[0])
                        if rest:
                            last = restgate1pool[0]
                        top = True
                        layerindex.append(len(OrderedCircuit))
                        ECParallelization(OrderedCircuit,layerindex)
                # choose gate that can reduce gate of last layer once
                elif not judge_rlc and gate1 == last and gate2 in lastlayergate2pool :
                    OrderedCircuit.append(gates)
                    top = False
                    if gate1[1:] not in llgate1index:
                        lastlayergate1pool.append(gate1)
                    else:
                        lastlayergate1pool[llgate1index.index(gate1[1:])] = gate1
                    rest = list(set(Sortedlocal2gate_R)-set(OrderedCircuit))
                    rest.sort(key=takeIndex)
                    last = gate2
                    # if it can't be reduced anymore, switch to next layer
                    if gate2 not in restgate1pool:   
                        restgate1pool = []       
                        for restgate in rest:
                            restcontent = sgateRegex.findall(str(restgate))
                            restgate1pool.append(restcontent[0])
                        if rest:
                            last = restgate1pool[0]
                        top = True
                        layerindex.append(len(OrderedCircuit))
                        ECParallelization(OrderedCircuit,layerindex)
                # on the top of the layer, if current gate2 cannot reduce the gate2 at same position in the last layer, order the gates whose gate2 can 
                # be reduced by following gate2 on the tail , choose the one that cannot as priority
                elif not judge_rlc and not judge_c2l2 and gate1 == last and gate2 not in restgate2pool :
                    OrderedCircuit.append(gates)
                    top = False
                    if gate1[1:] not in llgate1index:
                        lastlayergate1pool.append(gate1)
                    else:
                        lastlayergate1pool[llgate1index.index(gate1[1:])] = gate1
                    if gate2[1:] not in llgate2index:
                        lastlayergate2pool.append(gate2)
                    else:
                        lastlayergate2pool[llgate2index.index(gate2[1:])] = gate2
                    rest = list(set(Sortedlocal2gate_R)-set(OrderedCircuit))
                    rest.sort(key=takeIndex)
                    last = gate2
                    # if it can't be reduced anymore, switch to next layer
                    if gate2 not in restgate1pool:
                        restgate1pool = []      
                        for restgate in rest:
                            restcontent = sgateRegex.findall(str(restgate))
                            restgate1pool.append(restcontent[0])
                        if rest:
                            last = restgate1pool[0]
                        top = True
                        layerindex.append(len(OrderedCircuit))
                        ECParallelization(OrderedCircuit,layerindex)
                # rest cases
                elif not judge_rlc and not judge_c2mr2 and gate1 == last :
                    OrderedCircuit.append(gates)
                    top = False
                    if gate1[1:] not in llgate1index:
                        lastlayergate1pool.append(gate1)
                    else:
                        lastlayergate1pool[llgate1index.index(gate1[1:])] = gate1
                    if gate2[1:] not in llgate2index:
                        lastlayergate2pool.append(gate2)
                    else:
                        lastlayergate2pool[llgate2index.index(gate2[1:])] = gate2
                    rest = list(set(Sortedlocal2gate_R)-set(OrderedCircuit))
                    rest.sort(key=takeIndex)
                    last = gate2
                    # if it can't be reduced anymore, switch to next layer
                    if gate2 not in restgate1pool:
                        restgate1pool = []          
                        for restgate in rest:
                            restcontent = sgateRegex.findall(str(restgate))
                            restgate1pool.append(restcontent[0])
                        if rest:
                            last = restgate1pool[0]
                        top = True  
                        layerindex.append(len(OrderedCircuit))
                        ECParallelization(OrderedCircuit,layerindex)
            else:
                # choose gate that can reduce gate of last layer primarily
                if gate1 == last and gate2 in lastlayergate2pool :
                    OrderedCircuit.append(gates)
                    if gate1[1:] not in llgate1index:
                        lastlayergate1pool.append(gate1)
                    else:
                        lastlayergate1pool[llgate1index.index(gate1[1:])] = gate1
                    rest = list(set(Sortedlocal2gate_R)-set(OrderedCircuit))
                    rest.sort(key=takeIndex)
                    last = gate2
                    # if it can't be reduced anymore, switch to next layer
                    if gate2 not in restgate1pool:
                        restgate1pool = []          
                        for restgate in rest:
                            restcontent = sgateRegex.findall(str(restgate))
                            restgate1pool.append(restcontent[0])
                        if rest:
                            last = restgate1pool[0]
                        top = True 
                        layerindex.append(len(OrderedCircuit))
                        ECParallelization(OrderedCircuit,layerindex)
                # if current gate2 cannot reduce the gate2 at same position in the last layer, choose the one that can reduce next gate1
                elif not judge_c2l2 and gate1 == last and gate2 in restgate1pool :
                    OrderedCircuit.append(gates)
                    if gate1[1:] not in llgate1index:
                        lastlayergate1pool.append(gate1)
                    else:
                        lastlayergate1pool[llgate1index.index(gate1[1:])] = gate1
                    if gate2[1:] not in llgate2index:
                        lastlayergate2pool.append(gate2)
                    else:
                        lastlayergate2pool[llgate2index.index(gate2[1:])] = gate2
                    last = gate2
                # the bottom of the current layer
                elif not judge_c2ar1 and gate1 == last :
                    OrderedCircuit.append(gates)
                    if gate1[1:] not in llgate1index:
                        lastlayergate1pool.append(gate1)
                    else:
                        lastlayergate1pool[llgate1index.index(gate1[1:])] = gate1
                    if gate2[1:] not in llgate2index:
                        lastlayergate2pool.append(gate2)
                    else:
                        lastlayergate2pool[llgate2index.index(gate2[1:])] = gate2
                    rest = list(set(Sortedlocal2gate_R)-set(OrderedCircuit))
                    rest.sort(key=takeIndex)
                    restgate1pool = []          
                    for restgate in rest:
                        restcontent = sgateRegex.findall(str(restgate))
                        restgate1pool.append(restcontent[0])
                    if rest:
                            last = restgate1pool[0]
                    top = True
                    layerindex.append(len(OrderedCircuit))
                    ECParallelization(OrderedCircuit,layerindex)
        SortedOrderedCircuit = sorted(OrderedCircuit,key=takeIndex)

    
    # Parallelize the lonely gates if possible, then merge the total ordered circuit. 
    OrderedCircuit.extend(lonelygate)
    for l in range(2,len(OrderedCircuit)-1):
        gatel = OrderedCircuit[l]
        gatelindex = qubitRegex.findall(str(OrderedCircuit[l]))
        nextone = False
        #for i in range(len(OrderedCircuit)-1):
        for i in range(l):
            gatei = OrderedCircuit[i]
            gateiplus = OrderedCircuit[i+1]
            gateiindex = qubitRegex.findall(str(gatei))
            gateiplusindex = qubitRegex.findall(str(gateiplus))
            indexlonelyi = bool(set(gateiindex) & set(gatelindex))
            indexiiplus = bool(set(gateiindex) & set(gateiplusindex))
            indexliplus = bool(set(gatelindex) & set(gateiplusindex))
            # if a gate in lonely gates don't occupy same qubits with OrderedCircuit[i], we insert the gate into the position i+1 for paralleling,
            # then push backward the rest gates in OrderedCircuit.
            if nextone:
                break
            else:    
                if i > 0 :
                    gateiminus = OrderedCircuit[i-1]
                    gateiminusindex = qubitRegex.findall(str(gateiminus))
                    indexiiminus = bool(set(gateiindex) & set(gateiminusindex))
                    indexliminus = bool(set(gatelindex) & set(gateiminusindex))
                    if not indexlonelyi and indexiiplus and indexiiminus:
                        for index in range(l,i+1,-1):
                            OrderedCircuit[index] = OrderedCircuit[index-1]
                        OrderedCircuit[i+1] = gatel
                        nextone = True
                    elif not indexlonelyi and not indexliplus and not indexliminus:
                        for index in range(l,i+1,-1):
                            OrderedCircuit[index] = OrderedCircuit[index-1]
                        OrderedCircuit[i+1] = gatel
                        nextone = True
                else:
                    if not indexlonelyi and indexiiplus:
                        for index in range(l,i+1,-1):
                            OrderedCircuit[index] = OrderedCircuit[index-1]
                        OrderedCircuit[i+1] = gatel
                        nextone = True
                    elif not indexlonelyi and not indexliplus:
                        for index in range(l,i+1,-1):
                            OrderedCircuit[index] = OrderedCircuit[index-1]
                        OrderedCircuit[i+1] = gatel
                        nextone = True

    return OrderedCircuit
#----------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------
#                         DFS ordering algorithm Strategy
#                (DFS for the first layer and Edge Coloring for the rest,
#                concatenation-gate ordering in one layer takes priority)
def DFSOrdering(inputcircuit):
    # find number of qubits
    qubitindex = qubitRegex.findall(str(inputcircuit))
    qubits = int(max(list(map(int, qubitindex)))) + 1
    gate1pool = []
    gate2pool = []
    lastlayergate1pool = []
    lastlayergate2pool = []
    llgate1index = []
    llgate2index = []
    OrderedCircuit = []
    SortedOrderedCircuit = []
    lonelygate = []
    rest = []
    last = ''
    top = True
    # layerindex stores the indexes in OrderedCircuit from which new layers start, every element in layerindex list labels
    # the index ( in OrderedCircuit ) of the first gate in one layer 
    layerindex = []
    # sort all local2gates and group them into gate1pool and gate2pool
    for gates in inputcircuit:
            gatecontent = sgateRegex.findall(str(gates))
            gate1pool.append(gatecontent[0])
            gate2pool.append(gatecontent[1])

    totalgate = gate1pool + gate2pool
    lonelygatecount = dict(Counter(totalgate))

    # find lonely gates that cannot be reduced
    for gates in inputcircuit:
            gatecontent = sgateRegex.findall(str(gates))
            lonelygate1 = gatecontent[0]
            lonelygate2 = gatecontent[1]
            if lonelygatecount[lonelygate1] == 1 and lonelygatecount[lonelygate2] == 1 :
                lonelygate.append(gates)
            
    Sortedlocal2gate_R = sorted(list(set(inputcircuit) - set(lonelygate)),key=takeIndex)

    # find the superior order that maximize the amount of reduced single qubit gates
    while SortedOrderedCircuit != Sortedlocal2gate_R:
        for gates in Sortedlocal2gate_R:
            gatecontent = sgateRegex.findall(str(gates))
            gate1 = gatecontent[0]
            gate2 = gatecontent[1]

            # rest gates not been optimized
            rest = list(set(Sortedlocal2gate_R)-set(OrderedCircuit))
            rest.sort(key=takeIndex)
            restgate1pool = []
            restgate2pool = []
            lastlayercombo = []
            currentgate2opt = []
            
            
            for restgate in rest:
                restcontent = sgateRegex.findall(str(restgate))
                restgate1pool.append(restcontent[0])
                restgate2pool.append(restcontent[1])
                if restcontent[0] == gate1:
                    currentgate2opt.append(restcontent[1])

            # find qubit index of gate2 in last layer
            for lastlayergate1, lastlayergate2 in zip(lastlayergate1pool,lastlayergate2pool) :
                llgate1content = qubitRegex.findall(str(lastlayergate1))
                llgate2content = qubitRegex.findall(str(lastlayergate2))
                llgate1index.append(llgate1content[0])
                llgate2index.append(llgate2content[0])
            
            llgate1index = list(set(llgate1index))
            llgate2index = list(set(llgate2index))
            llgate1index.sort()
            llgate2index.sort()
            lastlayergate1pool.sort(key=qubitIndex)
            lastlayergate2pool.sort(key=qubitIndex)
            lastlayergatespool = list( set(lastlayergate1pool) | set(lastlayergate2pool) )
            lastlayergatespool.sort(key=qubitIndex)

            # find all possible gate combinations that can be optimized
            for lastlayergate1 in lastlayergatespool:
                for lastlayergate2 in lastlayergatespool:
                    lastlayercombo.append(lastlayergate1+lastlayergate2)

            if not OrderedCircuit:
                last = restgate1pool[0]

            if gates in OrderedCircuit:
                continue    

            judge_rlc = bool(set(rest) & set(lastlayercombo))
            judge_c2l2 = bool(set(currentgate2opt) & set(lastlayergate2pool))
            judge_c2mr2 = bool(set(currentgate2opt) - set(restgate2pool))
            judge_c2ar1 = bool(set(currentgate2opt) & set(restgate1pool))
            if top:
                # choose gate that can reduce gate of last layer twice
                if gate1 in lastlayergatespool and gate2 in lastlayergatespool :
                    OrderedCircuit.append(gates)
                    top = False
                    rest = list(set(Sortedlocal2gate_R)-set(OrderedCircuit))
                    rest.sort(key=takeIndex)
                    last = gate2
                    # if it can't be reduced anymore, switch to next layer
                    if gate2 not in restgate1pool:   
                        restgate1pool = []       
                        for restgate in rest:
                            restcontent = sgateRegex.findall(str(restgate))
                            restgate1pool.append(restcontent[0])
                        if rest:
                            last = restgate1pool[0]
                        top = True
                        layerindex.append(len(OrderedCircuit))
                        DFSParallelization(OrderedCircuit,layerindex)
                # choose gate that can reduce gate of last layer once
                elif not judge_rlc and gate1 == last and gate2 in lastlayergate2pool :
                    OrderedCircuit.append(gates)
                    top = False
                    if gate1[1:] not in llgate1index:
                        lastlayergate1pool.append(gate1)
                    else:
                        lastlayergate1pool[llgate1index.index(gate1[1:])] = gate1
                    rest = list(set(Sortedlocal2gate_R)-set(OrderedCircuit))
                    rest.sort(key=takeIndex)
                    last = gate2
                    # if it can't be reduced anymore, switch to next layer
                    if gate2 not in restgate1pool:   
                        restgate1pool = []       
                        for restgate in rest:
                            restcontent = sgateRegex.findall(str(restgate))
                            restgate1pool.append(restcontent[0])
                        if rest:
                            last = restgate1pool[0]
                        top = True
                        layerindex.append(len(OrderedCircuit))
                        DFSParallelization(OrderedCircuit,layerindex)
                # on the top of the layer, if current gate2 cannot reduce the gate2 at same position in the last layer, order the gates whose gate2 can 
                # be reduced by following gate2 on the tail , choose the one that cannot as priority
                elif not judge_rlc and not judge_c2l2 and gate1 == last and gate2 not in restgate2pool :
                    OrderedCircuit.append(gates)
                    top = False
                    if gate1[1:] not in llgate1index:
                        lastlayergate1pool.append(gate1)
                    else:
                        lastlayergate1pool[llgate1index.index(gate1[1:])] = gate1
                    if gate2[1:] not in llgate2index:
                        lastlayergate2pool.append(gate2)
                    else:
                        lastlayergate2pool[llgate2index.index(gate2[1:])] = gate2
                    rest = list(set(Sortedlocal2gate_R)-set(OrderedCircuit))
                    rest.sort(key=takeIndex)
                    last = gate2
                    # if it can't be reduced anymore, switch to next layer
                    if gate2 not in restgate1pool:
                        restgate1pool = []      
                        for restgate in rest:
                            restcontent = sgateRegex.findall(str(restgate))
                            restgate1pool.append(restcontent[0])
                        if rest:
                            last = restgate1pool[0]
                        top = True
                        layerindex.append(len(OrderedCircuit))
                        DFSParallelization(OrderedCircuit,layerindex)
                # rest cases
                elif not judge_rlc and not judge_c2mr2 and gate1 == last :
                    OrderedCircuit.append(gates)
                    top = False
                    if gate1[1:] not in llgate1index:
                        lastlayergate1pool.append(gate1)
                    else:
                        lastlayergate1pool[llgate1index.index(gate1[1:])] = gate1
                    if gate2[1:] not in llgate2index:
                        lastlayergate2pool.append(gate2)
                    else:
                        lastlayergate2pool[llgate2index.index(gate2[1:])] = gate2
                    rest = list(set(Sortedlocal2gate_R)-set(OrderedCircuit))
                    rest.sort(key=takeIndex)
                    last = gate2
                    # if it can't be reduced anymore, switch to next layer
                    if gate2 not in restgate1pool:
                        restgate1pool = []          
                        for restgate in rest:
                            restcontent = sgateRegex.findall(str(restgate))
                            restgate1pool.append(restcontent[0])
                        if rest:
                            last = restgate1pool[0]
                        top = True  
                        layerindex.append(len(OrderedCircuit))
                        DFSParallelization(OrderedCircuit,layerindex)
            else:
                # choose gate that can reduce gate of last layer primarily
                if gate1 == last and gate2 in lastlayergate2pool :
                    OrderedCircuit.append(gates)
                    if gate1[1:] not in llgate1index:
                        lastlayergate1pool.append(gate1)
                    else:
                        lastlayergate1pool[llgate1index.index(gate1[1:])] = gate1
                    rest = list(set(Sortedlocal2gate_R)-set(OrderedCircuit))
                    rest.sort(key=takeIndex)
                    last = gate2
                    # if it can't be reduced anymore, switch to next layer
                    if gate2 not in restgate1pool:
                        restgate1pool = []          
                        for restgate in rest:
                            restcontent = sgateRegex.findall(str(restgate))
                            restgate1pool.append(restcontent[0])
                        if rest:
                            last = restgate1pool[0]
                        top = True 
                        layerindex.append(len(OrderedCircuit))
                        DFSParallelization(OrderedCircuit,layerindex)
                # if current gate2 cannot reduce the gate2 at same position in the last layer, choose the one that can reduce next gate1
                elif not judge_c2l2 and gate1 == last and gate2 in restgate1pool :
                    OrderedCircuit.append(gates)
                    if gate1[1:] not in llgate1index:
                        lastlayergate1pool.append(gate1)
                    else:
                        lastlayergate1pool[llgate1index.index(gate1[1:])] = gate1
                    if gate2[1:] not in llgate2index:
                        lastlayergate2pool.append(gate2)
                    else:
                        lastlayergate2pool[llgate2index.index(gate2[1:])] = gate2
                    last = gate2
                # the bottom of the current layer
                elif not judge_c2ar1 and gate1 == last :
                    OrderedCircuit.append(gates)
                    if gate1[1:] not in llgate1index:
                        lastlayergate1pool.append(gate1)
                    else:
                        lastlayergate1pool[llgate1index.index(gate1[1:])] = gate1
                    if gate2[1:] not in llgate2index:
                        lastlayergate2pool.append(gate2)
                    else:
                        lastlayergate2pool[llgate2index.index(gate2[1:])] = gate2
                    rest = list(set(Sortedlocal2gate_R)-set(OrderedCircuit))
                    rest.sort(key=takeIndex)
                    restgate1pool = []          
                    for restgate in rest:
                        restcontent = sgateRegex.findall(str(restgate))
                        restgate1pool.append(restcontent[0])
                    if rest:
                            last = restgate1pool[0]
                    top = True
                    layerindex.append(len(OrderedCircuit))
                    DFSParallelization(OrderedCircuit,layerindex)
        SortedOrderedCircuit = sorted(OrderedCircuit,key=takeIndex)  
    
    # Parallelize the lonely gates if possible, then merge the total ordered circuit. 
    OrderedCircuit.extend(lonelygate)
    DFSlayer = []
    for l in range(2,len(OrderedCircuit)-1):
        gatel = OrderedCircuit[l]
        gatelindex = qubitRegex.findall(str(OrderedCircuit[l]))
        if int(gatelindex[1]) == (qubits-1) :
            DFSlayer.append(l)
        nextone = False
        for i in range(l):
            gatei = OrderedCircuit[i]
            gateiplus = OrderedCircuit[i+1]
            gateiindex = qubitRegex.findall(str(gatei))
            gateiplusindex = qubitRegex.findall(str(gateiplus))
            indexlonelyi = bool(set(gateiindex) & set(gatelindex))
            indexiiplus = bool(set(gateiindex) & set(gateiplusindex))
            indexliplus = bool(set(gatelindex) & set(gateiplusindex))
            # if a gate in lonely gates don't occupy same qubits with OrderedCircuit[i], we insert the gate into the position i+1 for paralleling,
            # then push backward the rest gates in OrderedCircuit.
            if DFSlayer:
                if i > DFSlayer[0]:
                    gateiminus = OrderedCircuit[i-1]
                    gateiminusindex = qubitRegex.findall(str(gateiminus))
                    indexiiminus = bool(set(gateiindex) & set(gateiminusindex))
                    indexliminus = bool(set(gatelindex) & set(gateiminusindex))
                    if not indexlonelyi and indexiiplus and indexiiminus:
                        for index in range(l,i+1,-1):
                            OrderedCircuit[index] = OrderedCircuit[index-1]
                        OrderedCircuit[i+1] = gatel
                        nextone = True
                    elif not indexlonelyi and not indexliplus and not indexliminus:
                        for index in range(l,i+1,-1):
                            OrderedCircuit[index] = OrderedCircuit[index-1]
                        OrderedCircuit[i+1] = gatel
                        nextone = True
            if nextone:
                break

    return OrderedCircuit
#---------------------------------------------------------------------------------------- 

def concatenate(list):
    concatenation = ''
    for terms in list:
        concatenation += terms
    return concatenation

def Depth4block(sub1,sub2):
    depth4block = []
    depth4block.extend(sub1)
    depth4block.extend(sub2)
    depth4block.extend(sub1)
    depth4block.extend(sub2)
    return depth4block

def Depth5block(sub1,sub2):
    depth5block = []
    depth5block.extend(sub2)
    depth5block.extend(sub1)
    depth5block.extend(sub2)
    depth5block.extend(sub1)
    depth5block.extend(sub2)
    return depth5block


#----------------------------------------------------------------------------------------
#                         Subcircuit Decomposition algorithm 
#      (SubcircuitDecomposition function can decompose k-local pauli strings
#      in the form of Pauli+Index(e.g. X0Y1X2Z3) into 2-locals (e.g. X0Z1,X1X2).
#       Originla paper describing techniques used in algorithm can be found at:
#              https://www.nature.com/articles/s41467-021-25196-0 )
def SubcircuitDecomposition(localn_Hamiltonian):
    if len(pauliRegex.findall(str(localn_Hamiltonian))) < 3:
        #print('2-local circuit does not need to be decomposed')
        return [localn_Hamiltonian]
    depth5decomposition = []
    temp1 = ''
    temp2 = ''
    firstlocal = []
    secondlocal = []
    paulik = str(localn_Hamiltonian)
    layers = len(pauliRegex.findall(str(localn_Hamiltonian))) - 2
    for i in range(layers):
        localnblock = pauliRegex.findall(paulik)
        secondpauli = gateRegex.findall(str(localnblock[1]))
        secondqubit = qubitRegex.findall(str(localnblock[1]))
        if secondpauli[0] == 'X':
            temp1 = 'Y'
            temp2 = 'Z'
        elif secondpauli[0] == 'Y':
            temp1 = 'Z'
            temp2 = 'X'
        else:
            temp1 = 'X'
            temp2 = 'Y'    
        firstlocal.append(localnblock[0] + temp1 + secondqubit[0])
        restlocal = concatenate(localnblock[2:])
        paulik = temp2 + secondqubit[0] + restlocal
        secondlocal.append(temp2 + secondqubit[0] + restlocal)
    secondblock = Depth4block([firstlocal[-1]],[secondlocal[-1]])
    if len(firstlocal) > 1:
        for localterm in firstlocal[-2::-1]:
            depth5decomposition = Depth5block([localterm],secondblock)
            secondblock = depth5decomposition
    else:
        depth5decomposition = secondblock

    return depth5decomposition

#    Decomposing a k-local Hamiltonian into 2-local Hamiltonians using Subcircuit Decomposition algorithm
#         Input : a list of k-local Paulistrings such as XYZIIXXIIYY or X0Y1Z2X5X6Y9Y10
#        output : a list of 2-local Paulistrings such as Y1Y2, X2Z3, ...
def HamiltonianDecomposition(klocal_Hamiltonian):
    decomposedklocal = []
    decomposedoutput = []
    if len(qubitRegex.findall(str(klocal_Hamiltonian))) > 0:
        for term in klocal_Hamiltonian:
            decomposedoutput.extend(SubcircuitDecomposition(term))
    else:
        for term in klocal_Hamiltonian:
            klocalblockraw = allgateRegex.findall(str(term))
            for i in range(len(klocalblockraw)):
                klocalblockraw[i] = klocalblockraw[i] + str(i)
            klocalblock = pauliRegex.findall(str(klocalblockraw))
            klocal = ''
            if len(klocalblock) > 1:
                for paulis in klocalblock:
                    klocal += paulis
                decomposedklocal.extend([klocal])
        for i in range(len(decomposedklocal)):
            decomposedoutput.extend(SubcircuitDecomposition(decomposedklocal[i]))

    return decomposedoutput
#----------------------------------------------------------------------------------------