#%%
import gurobipy as gp
from gurobipy import GRB
import pandapower as pp
import pandapower.networks as pn
from copy import deepcopy
import pandas as pd
import numpy as np
from math import sqrt
import os
from pandapower.networks.power_system_test_cases import case1888rte
#%%
# Mesh 네트워크 컨버전시 주의사항
# (1) gen -> sgen 변경시 pp.runpp 돌리면 gen마다 결정되는 Q가 있는데, 이걸 제일 나중에 설정해줘야함 (trafo 등 다른 요소 변경 후에)
# (2) trafo -> line 변경시 임피던스 좀 줄이는게 좋음

# 일반적인 프로세스
# 1 - length_km 줄여가면서 trafo -> line 변경 시도 / runpp 수렴 안하면 length_km 더 줄임
# 2 - 1이 끝난 후 gen -> sgen 변경 시도 / runpp 수렴 안하면 length_km 더 줄여서 1번부터
def ieee1888(mode = 'ppf', sgen = True):
    net = case1888rte()
    net.name = 'ieee1888'

    net.line.c_nf_per_km = 0.0
    net.ext_grid.vm_pu = 1
    net.ext_grid.va_degree = 0
    #net.line.length_km = net.line.length_km*0.1
    net.shunt.in_service = False
#    net.load.q_mvar = abs(net.load.q_mvar)

    for sgen in net.sgen.index:
        net.sgen.drop(sgen, inplace = True)

    net.trafo.tap_pos = 0
    net.line.length_km = net.line.length_km*0.5
    
    check = dict()
    for idx in net.line.index:
        bus1, bus2 = net.line.loc[idx, 'from_bus'], net.line.loc[idx, 'to_bus']
        if ((bus1, bus2) in check) or ((bus2, bus1) in check):
            net.line.drop(idx, inplace= True)
        else:
            check[(bus1, bus2)] = 1
    if mode == 'opf':
        return net
    elif mode == 'vvc':
        pass
    elif mode == 'qdispatch':
        change_trafo_to_line(net)
    else:
        change_trafo_to_line(net)
        pass
    if not sgen:
        for sgen in net.sgen.index:
            net.sgen.drop(sgen, inplace = True)
    check = dict()
    for idx in net.line.index:
        bus1, bus2 = net.line.loc[idx, 'from_bus'], net.line.loc[idx, 'to_bus']
        if ((bus1, bus2) in check) or ((bus2, bus1) in check):
            net.line.drop(idx, inplace= True)
        else:
            check[(bus1, bus2)] = 1
    change_gen_to_sgen(net)
    return net

#%%
def change_gen_to_sgen(net):
    pp.runpp(net)
    gens = list(net.gen.index)
    for gen in gens:
        bus, p, q = net.gen.loc[gen, 'bus'], net.res_gen.loc[gen, 'p_mw'], net.res_gen.loc[gen, 'q_mvar']
        
        sn_mva = (p**2 + q**2)**(1/2)
        sgen = pp.create_sgen(net, bus, p, q, sn_mva)
        net.gen.drop(gen, inplace = True)
        #try: < 이걸 하나하나 체크하면 중간에 뭐 q 세팅이 바뀌는지 안돔
        #    pp.runpp(net)
        #except:
        #    print("Nonconvergence in adding sgen ", sgen)
    return net

def change_trafo_to_line(net):
    trafos = list(net.trafo.index)
    for trafo in trafos:
        vk_percent, vkr_percent, sn_mva = net.trafo.loc[trafo, 'vk_percent'], net.trafo.loc[trafo, 'vkr_percent'], net.trafo.loc[trafo, 'sn_mva']
        z, r = vk_percent/100.*net.sn_mva/sn_mva, vkr_percent/100.*net.sn_mva/sn_mva
        x = (z**2 - r**2)**(1/2)
        bus1, bus2 = net.trafo.loc[trafo, 'hv_bus'], net.trafo.loc[trafo, 'lv_bus']
        net.trafo.drop(trafo, inplace = True)
        line = pp.create_line_from_parameters(net, bus1, bus2, 1, r, x, 0, 9999)
        
        # Check if the replacement is valid:
#        try:
#            pp.runpp(net)
#        except:
#            print("Nonconvergence in adding line ", line)
    ext_grid_bus = net.ext_grid.bus[0]
    vn_kv = net.bus.loc[ext_grid_bus, 'vn_kv']
    net.bus.vn_kv = vn_kv
    return net

#%%

# %%
#net.load.p_mw = net.load.p_mw*1.5
# %%
#pp.runpp(net)

# %%
#net.load.p_mw = net.load.p_mw*1.01
#net.load.q_mvar = net.load.q_mvar*1.01
#pp.runpp(net)

# %%
#net.load.p_mw = net.load.p_mw/1.01
#net.load.q_mvar = net.load.q_mvar/1.01
#pp.runpp(net)

# %%
def test():
    net = ieee1888()
    pp.runpp(net, max_iteration = 10)
    net.res_bus
    net.name = 'ieee1888'

    net.line.c_nf_per_km = 0.0
    net.ext_grid.vm_pu = 1
    net.ext_grid.va_degree = 0
    #net.line.length_km = net.line.length_km*0.1
    net.shunt.in_service = False
    #    net.load.q_mvar = abs(net.load.q_mvar)

    for sgen in net.sgen.index:
        net.sgen.drop(sgen, inplace = True)

    check = dict()
    for idx in net.line.index:
        bus1, bus2 = net.line.loc[idx, 'from_bus'], net.line.loc[idx, 'to_bus']
        if ((bus1, bus2) in check) or ((bus2, bus1) in check):
            net.line.drop(idx, inplace= True)
        else:
            check[(bus1, bus2)] = 1
    net.trafo.tap_pos = 0
    net.line.length_km = net.line.length_km*0.5
    change_trafo_to_line(net)
    pp.runpp(net)
    net.res_gen
    net.res_bus
    change_gen_to_sgen(net)
    pp.runpp(net, max_iteration = 10)
    # %%
