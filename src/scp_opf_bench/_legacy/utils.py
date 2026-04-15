# OPF 구현
#%%
#import networks
import gurobipy as gp
from gurobipy import GRB
import pandapower as pp
import pandapower.networks as pn
import numpy as np
import os
import pickle
GlobalSbase = 1e6

def InitializeConstants(net, return_shunt = False, Sbase = 1e6, update = False):
    if os.path.exists('networkdata/initialize/' + net.name + '.pkl') and (not update):
        with open('networkdata/initialize/' + net.name + '.pkl', "rb") as f:
                Vbase, Ibase, Sbase, Zbase, Buses, SubstationBuses, NonsubstationBuses, BidirectionalLines, UnidirectionalLines, UnidirectionalTrafo, NonTrafoUnidirectionalLines, TrafoTappos, TrafoTapStepPercent, LineTupleToIndex, TrafoTupleToIndex, R, X, AdjacencyList, DirectedLinesForEachTree, UnidirectionalLineTupleToIndex = pickle.load(f)
                return Vbase, Ibase, Sbase, Zbase, Buses, SubstationBuses, NonsubstationBuses, BidirectionalLines, UnidirectionalLines, UnidirectionalTrafo, NonTrafoUnidirectionalLines, TrafoTappos, TrafoTapStepPercent, LineTupleToIndex, TrafoTupleToIndex, R, X, AdjacencyList, DirectedLinesForEachTree, UnidirectionalLineTupleToIndex

    Vbase = 4.16*1000
    if len(net.ext_grid) > 0:
        Vbase = net.bus.loc[net.ext_grid.bus.iloc[0], 'vn_kv']*1000 # Vbase is the first ext_grid's vn_kv 
    else:
        bus = net.gen.loc[net.gen.slack].loc[0, 'bus']
        Vbase = net.gen.loc[net.gen.slack].loc[0, 'vm_pu']*net.bus.loc[bus, 'vn_kv']*1000
    Sbase = Sbase # trafo sn_mva # 이거 바꾸면 RefreshPQ도 바꿔줘야함
    Ibase = Sbase/Vbase 
    Zbase = Vbase/Ibase
    Buses = list(net.bus.index)
    SubstationBuses = None
    if len(net.ext_grid) > 0:
        SubstationBuses = list(net.ext_grid.bus)
    else:
        SubstationBuses = list(net.gen.loc[net.gen.slack].bus)

    NonsubstationBuses = list(net.bus.index)
    for bus in SubstationBuses:
        NonsubstationBuses.remove(bus)

    OpenedLines = []
    for switch in net.switch.loc[(net.switch.et == 'l') & (~net.switch.closed)].index:
        OpenedLines.append(net.switch.loc[switch, 'element'])

    BidirectionalLines = []
    LineTupleToIndex = dict()
    UnidirectionalLineTupleToIndex = dict()
    ###################################################################################
    ########## R, X, Lines (including trafo)/양방향 다 고려 #############################
    R, X, G, C, GB, CB  = dict(), dict(), dict(), dict(), dict(), dict()
    for bus in net.bus.index:
        GB[bus], CB[bus] = 0, 0 

    for line in net.line.index:
        if line in OpenedLines: # Check if it is opened
            continue
        from_bus, to_bus = net.line.loc[line, 'from_bus'], net.line.loc[line, 'to_bus']
        r, x = net.line.loc[line, 'r_ohm_per_km']*net.line.loc[line, 'length_km']/Zbase, net.line.loc[line, 'x_ohm_per_km']*net.line.loc[line, 'length_km']/Zbase
        g, c = net.line.loc[line, 'g_us_per_km']*net.line.loc[line, 'length_km']*1e-6*Zbase/2, net.line.loc[line, 'c_nf_per_km']*2*3.141592*60*1e-9*net.line.loc[line, 'length_km']*net.line.loc[line, 'parallel']/2*Zbase
        BidirectionalLines.append((from_bus, to_bus)), BidirectionalLines.append((to_bus, from_bus))
        R[(from_bus, to_bus)], R[(to_bus, from_bus)] = r, r
        X[(from_bus, to_bus)], X[(to_bus, from_bus)] = x, x
        G[(from_bus, to_bus)], G[(to_bus, from_bus)] = g, g
        C[(from_bus, to_bus)], C[(to_bus, from_bus)] = c, c
        GB[from_bus] += g
        GB[to_bus] += g
        CB[from_bus] += c
        CB[to_bus] += c
        LineTupleToIndex[(from_bus, to_bus)] = line
        UnidirectionalLineTupleToIndex[(from_bus, to_bus)] = line



    for switch in net.switch.loc[(net.switch.et=='b') & (net.switch.closed)].index:
        from_bus, to_bus = net.switch.loc[switch, 'bus'], net.switch.loc[switch, 'element']
        BidirectionalLines.append((from_bus, to_bus)), BidirectionalLines.append((to_bus, from_bus))
        R[(from_bus, to_bus)], R[(to_bus, from_bus)] = 0, 0
        X[(from_bus, to_bus)], X[(to_bus, from_bus)] = 0, 0

    TrafoTapStepPercent = dict()
    for tf in net.trafo.index:
        from_bus, to_bus = net.trafo.loc[tf, 'hv_bus'], net.trafo.loc[tf, 'lv_bus']
        Vrated = net.trafo.loc[tf, 'vn_hv_kv']*1000
        Irated = net.trafo.loc[tf, 'sn_mva']*1000000/Vrated

        Zreftrafo = (Vrated/1000)**2*net.sn_mva/net.trafo.loc[tf, 'sn_mva']

        #z = net.trafo.loc[tf, 'vk_percent']/100*net.sn_mva/net.trafo.loc[tf, 'sn_mva']*Zreftrafo/Zbase
        #r = net.trafo.loc[tf, 'vkr_percent']/100*net.sn_mva/net.trafo.loc[tf, 'sn_mva']*Zreftrafo/Zbase
        z = net.trafo.loc[tf, 'vk_percent']/100*Vrated/Irated/Zbase
        r = net.trafo.loc[tf, 'vkr_percent']/100*Vrated/Irated/Zbase
        x = (z**2 - r**2)**(1/2)
        R[(from_bus, to_bus)], R[(to_bus, from_bus)] = r, r
        X[(from_bus, to_bus)], X[(to_bus, from_bus)] = x, x

        BidirectionalLines.append((from_bus, to_bus)), BidirectionalLines.append((to_bus, from_bus))
        TrafoTapStepPercent[(from_bus, to_bus)] = net.trafo.loc[tf, 'tap_step_percent']
        TrafoTapStepPercent[(to_bus, from_bus)] = TrafoTapStepPercent[(from_bus, to_bus)]
    
    #################################################
    ########## Adj Mat ###########################
    AdjacencyList = {}
    for bus in Buses:
        AdjacencyList[bus] = []
    for node1, node2 in BidirectionalLines:
        if node1 not in AdjacencyList:
            AdjacencyList[node1] = []
        if node2 not in AdjacencyList:
            AdjacencyList[node2] = [] 
        AdjacencyList[node2].append(node1)  # If the graph is undirected
    #################################################

    ########## Directed Graph (Not used) #########################
    DirectedLinesForEachTree = dict()
    Check = set()
    def DepthSearch(bus, substation):
        for adjbus in AdjacencyList[bus]:
            if not adjbus in Check:
                DirectedLinesForEachTree[substation].append((bus,adjbus)), Check.add(adjbus)
                DepthSearch(adjbus, substation)

    for substation in SubstationBuses:
        DirectedLinesForEachTree[substation] = []
        Check.add(substation)
        DepthSearch(substation, substation)
    ##############################################################

    ########### Trafo Line ################################
    UnidirectionalTrafo = []
    TrafoTupleToIndex = dict()
    TrafoTappos = dict()
    for trafo in net.trafo.index:
        from_bus, to_bus = net.trafo.loc[trafo, 'hv_bus'], net.trafo.loc[trafo, 'lv_bus']
        min_pos, max_pos = net.trafo.loc[trafo, 'tap_min'], net.trafo.loc[trafo, 'tap_max']
        UnidirectionalTrafo.append((from_bus, to_bus))
        TrafoTupleToIndex[(from_bus, to_bus)] = trafo
        UnidirectionalLineTupleToIndex[(from_bus, to_bus)] = trafo
        TrafoTappos[(from_bus, to_bus)] = net.trafo.loc[trafo, 'tap_pos']
    ######################################################
    
    ########### only one side of two direction is included #############
    NonTrafoUnidirectionalLines = []

    for line in BidirectionalLines:
        line1, line2 = line, tuple(reversed(line))
        if not (line1 in NonTrafoUnidirectionalLines + UnidirectionalTrafo) and not (line2 in NonTrafoUnidirectionalLines + UnidirectionalTrafo):
            NonTrafoUnidirectionalLines.append(line)
    ####################################################################

    UnidirectionalLines = UnidirectionalTrafo + NonTrafoUnidirectionalLines

    with open('networkdata/initialize/' + net.name + '.pkl', 'wb') as file:
        pickle.dump([Vbase, Ibase, Sbase, Zbase, Buses, SubstationBuses, NonsubstationBuses, BidirectionalLines, UnidirectionalLines, UnidirectionalTrafo, NonTrafoUnidirectionalLines, TrafoTappos, TrafoTapStepPercent, LineTupleToIndex, TrafoTupleToIndex, R, X, AdjacencyList, DirectedLinesForEachTree, UnidirectionalLineTupleToIndex], file)
    if not return_shunt:
        return Vbase, Ibase, Sbase, Zbase, Buses, SubstationBuses, NonsubstationBuses, BidirectionalLines, UnidirectionalLines, UnidirectionalTrafo, NonTrafoUnidirectionalLines, TrafoTappos, TrafoTapStepPercent, LineTupleToIndex, TrafoTupleToIndex, R, X, AdjacencyList, DirectedLinesForEachTree, UnidirectionalLineTupleToIndex
    if return_shunt:
        return Vbase, Ibase, Sbase, Zbase, Buses, SubstationBuses, NonsubstationBuses, BidirectionalLines, UnidirectionalLines, UnidirectionalTrafo, NonTrafoUnidirectionalLines, TrafoTappos, TrafoTapStepPercent, LineTupleToIndex, TrafoTupleToIndex, R, X, G, C, GB, CB, AdjacencyList, DirectedLinesForEachTree, UnidirectionalLineTupleToIndex

#%%
def RefreshPQ(net, Sbase = 1e6):
    P, Q = dict(), dict() # Initialization
    Sbase = Sbase
    Buses = list(net.bus.index)
    for bus in Buses:
        P[bus], Q[bus] = 0, 0 

    for load in net.load.index:
        bus = net.load.loc[load, 'bus']
        p, q = net.load.loc[load, 'p_mw']*net.load.loc[load, 'scaling']*1000000/Sbase, net.load.loc[load, 'q_mvar']*net.load.loc[load, 'scaling']*1000000/Sbase
        P[bus], Q[bus] = P[bus] - p, Q[bus] - q

    for sgen in net.sgen.index:
        bus = net.sgen.loc[sgen, 'bus']
        p = net.sgen.loc[sgen, 'p_mw']*net.sgen.loc[sgen, 'scaling']*1000000/Sbase
        P[bus]= P[bus] + p
    return P, Q

def makeYbusmatrix(net, Sbase):
    Vbase, Ibase, Sbase, Zbase, Buses, SubstationBuses, NonsubstationBuses, BidirectionalLines, UnidirectionalLines, UnidirectionalTrafo, NonTrafoUnidirectionalLines, TrafoTappos, TrafoTapStepPercent, LineTupleToIndex, TrafoTupleToIndex, R, X, AdjacencyList, DirectedLinesForEachTree, UnidirectionalLineTupleToIndex = InitializeConstants(net, Sbase = Sbase)

    #pp.runpp(net)
    #Ybus = net["_ppc"]["internal"]["Ybus"]InitializeConstants
    N = len(net.bus.index)
    Y = np.zeros(shape = (N,N), dtype =complex)
    
    for key, value in zip(R.keys(), R.values()): ## only for SDP/ Switch가 zero impedance 가질 경우 안풀림 
        if value < 1e-10:
            R[key] = 1e-3
            X[key] = 1e-3

    for i, busi in enumerate(net.bus.index):
        for j, busj in enumerate(net.bus.index):
            if (busi, busj) in R.keys():
                z = complex(R[busi,busj], X[busi, busj])
                y = 1/z
                Y[i,j] = -y
            
            #Y[i,j] = Ybus[busi, busj]
    for i in range(Y.shape[0]):
        for j in range(Y.shape[1]):
            if i != j:
                Y[i,i] += -Y[i][j] 
    return Y