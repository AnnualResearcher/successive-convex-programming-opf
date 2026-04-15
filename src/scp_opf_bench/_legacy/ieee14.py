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
from pandapower.networks.power_system_test_cases import case14
#%%
def ieee14(mode = 'qdispatch', sgen = True):
    net = case14()
    net.name = 'ieee14'
    net.line.c_nf_per_km = 0.0
    net.ext_grid.vm_pu = 1
    net.ext_grid.va_degree = 0
    net.shunt.in_service = False

    if mode == 'opf':
        return net
    elif mode == 'vvc':
        change_gen_to_sgen(net)
    elif mode == 'qdispatch':
        change_gen_to_sgen(net)
        change_trafo_to_line(net)
    else:
        change_gen_to_sgen(net)
        change_trafo_to_line(net)

    if not sgen:
        for sgen in net.sgen.index:
            net.sgen.drop(sgen, inplace = True)
    return net

# %%
def change_gen_to_sgen(net):
    gens = list(net.gen.index)
    for gen in gens:
        bus, p, q = net.gen.loc[gen, 'bus'], net.gen.loc[gen, 'p_mw'], net.gen.loc[gen, 'max_q_mvar']
        
        sn_mva = (p**2 + q**2)**(1/2)
        sgen = pp.create_sgen(net, bus, p, 0, sn_mva)
        net.gen.drop(gen, inplace = True)
        try:
            pp.runpp(net)
        except:
            print("Nonconvergence in adding sgen ", sgen)
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
        try:
            pp.runpp(net)
        except:
            print("Nonconvergence in adding line ", line)
    ext_grid_bus = net.ext_grid.bus[0]
    vn_kv = net.bus.loc[ext_grid_bus, 'vn_kv']
    net.bus.vn_kv = vn_kv
    return net




# %%
