import gurobipy as gp
from gurobipy import GRB
import pandapower as pp
import pandapower.networks as pn
from copy import deepcopy
import pandas as pd
import numpy as np
from math import sqrt
import os

def three_bus_example():
    net = pp.create_empty_network()
    buses = [0,1,2]
    net.name = '3_bus_net'
    for bus in buses:
        pp.create_bus(net, vn_kv = 4.16)
    pp.create_ext_grid(net, 0)
    pp.create_line_from_parameters(net, 0,1, 1, 1., 54.,0,1000)
    pp.create_line_from_parameters(net, 1,2, 1, 1., 54.,0,1000)
    pp.create_line_from_parameters(net, 0,2, 1, 1., 54.,0,1000)

    pp.create_load(net, 1, p_mw = 10, q_mvar = 20)
    pp.create_load(net, 2, p_mw = 30, q_mvar = 10)
    pp.create_sgen(net, 2, p_mw = 0, q_mvar = 10, sn_mva = 50)
    net.line.length_km = 0.001
    return net


def ieee9_radial2(ppf = False, opf = False, vvc = False, use_sgen = False, shunt = False):
    if int(ppf) + int(opf) + int(vvc) >= 2 or int(ppf) + int(opf) + int(vvc) == 0:
        raise('Only one among pp, opf, vvc needs to be set True')
    
    net = pn.case9() # Example network
    net.name = 'ieee9_radial'
    pp.create_switch(net, 0, 0, et = 'l', closed = True)
    pp.create_switch(net, 3, 1, et = 'l', closed = True)
    pp.create_switch(net, 0, 7, et = 'b', closed= False)
    net.line.drop(4, inplace = True)
    net.ext_grid.loc[0, 'vm_pu'] = 1.05
    pp.create_load(net,1, p_mw = 100, q_mvar = 35)
    pp.create_load(net,2, p_mw = 100, q_mvar = 35)
    net.line.length_km = 0.1
    net.line.c_nf_per_km = 0
    net.line.r_ohm_per_km = net.line.r_ohm_per_km +1
    net.line.loc[3,'from_bus'] = 5
    net.line.loc[3,'to_bus'] = 2
    net.line.loc[8,'from_bus'] = 3
    net.line.loc[8,'to_bus'] = 8
    net.line.loc[7,'from_bus'] = 8
    net.line.loc[7,'to_bus'] = 7 
    net.line.loc[5,'from_bus'] = 7
    net.line.loc[5,'to_bus'] = 6 
    if ppf:
        for i in net.gen.index:
            net.gen.drop(i, inplace = True)
        for i in net.poly_cost.index:
            net.poly_cost.drop(i, inplace = True)
        for i in net.shunt.index:
            net.shunt.drop(i, inplace = True)
    if opf:
        net.gen.loc[0, 'slack'] = True
        for i in net.ext_grid.index:
            net.ext_grid.drop(i, inplace = True)        
        for i in net.poly_cost[net.poly_cost.et == 'ext_grid'].index:
            net.poly_cost.drop(i, inplace = True)
        pass
   
    if vvc:
        for i in net.gen.index:
            net.gen.drop(i, inplace = True)
        for i in net.poly_cost.index:
            net.poly_cost.drop(i, inplace = True)   

        net.line.drop(8, inplace = True)     
        trafo_std_type = {
            'sn_mva' : 10,
            'vn_hv_kv' : net.bus.vn_kv[3],
            'vn_lv_kv' : net.bus.vn_kv[8],
            'vk_percent' : 0.008,
            'vkr_percent' : 0.001,
            'pfe_kw' : 0,
            'i0_percent' : 0,
            'shift_degree' : 0,
            'tap_side' : 'lv',
            'tap_min' : -9,
            'tap_max' : 9,
            'tap_neutral' : 0,
            'tap_step_percent' : 1.5,
            'tap_step_degree' : 0,
            'tap_phase_shifter' : False,
            'tap_dependency_table' : False,
            'tap_changer_type' : 'Ratio'
        }
        pp.create_std_type(net, trafo_std_type, name = 'OLTC', element= 'trafo')

        pp.create_transformer(net, hv_bus = 3, lv_bus=8, std_type = 'OLTC', tap_pos = 0)
    if use_sgen:
        pp.create_sgen(net,6,p_mw = 150, q_mvar = 50)
        pp.create_sgen(net,5,p_mw = 120, q_mvar = 100) 
        net.sgen.sn_mva = (net.sgen.p_mw**2 + net.sgen.q_mvar**2)**(1/2)*1.35
    if shunt:
        net.line.c_nf_per_km = 2100.
        net.line.g_us_per_km = 20.
    return net

def ieee9_mesh2(ppf = False, opf = False, vvc = False, use_sgen = False, shunt = False):
    if int(ppf) + int(opf) + int(vvc) >= 2 or int(ppf) + int(opf) + int(vvc) == 0:
        raise('Only one among pp, opf, vvc needs to be set True')
    
    net = pn.case9() # Example network
    pp.create_switch(net, 0, 0, et = 'l', closed = True)
    pp.create_switch(net, 3, 1, et = 'l', closed = True)
    pp.create_switch(net, 0, 7, et = 'b', closed= False)
    #net.line.drop(4, inplace = True)
    net.ext_grid.loc[0, 'vm_pu'] = 1.05
    pp.create_load(net,1, p_mw = 100, q_mvar = 35)
    pp.create_load(net,2, p_mw = 100, q_mvar = 35)
    net.line.length_km = 0.1
    net.line.c_nf_per_km = 0
    net.line.r_ohm_per_km = net.line.r_ohm_per_km +1
    net.line.loc[3,'from_bus'] = 5
    net.line.loc[3,'to_bus'] = 2
    net.line.loc[8,'from_bus'] = 3
    net.line.loc[8,'to_bus'] = 8
    net.line.loc[7,'from_bus'] = 8
    net.line.loc[7,'to_bus'] = 7 
    net.line.loc[5,'from_bus'] = 7
    net.line.loc[5,'to_bus'] = 6 
    if ppf:
        for i in net.gen.index:
            net.gen.drop(i, inplace = True)
        for i in net.poly_cost.index:
            net.poly_cost.drop(i, inplace = True)
        for i in net.shunt.index:
            net.shunt.drop(i, inplace = True)
    if opf:
        net.gen.loc[0, 'slack'] = True
        for i in net.ext_grid.index:
            net.ext_grid.drop(i, inplace = True)        
        for i in net.poly_cost[net.poly_cost.et == 'ext_grid'].index:
            net.poly_cost.drop(i, inplace = True)
        pass
   
    if vvc:
        for i in net.gen.index:
            net.gen.drop(i, inplace = True)
        for i in net.poly_cost.index:
            net.poly_cost.drop(i, inplace = True)   

        #net.line.drop(8, inplace = True)     
        trafo_std_type = {
            'sn_mva' : 10,
            'vn_hv_kv' : net.bus.vn_kv[3],
            'vn_lv_kv' : net.bus.vn_kv[8],
            'vk_percent' : 0.008,
            'vkr_percent' : 0.001,
            'pfe_kw' : 0,
            'i0_percent' : 0,
            'shift_degree' : 0,
            'tap_side' : 'lv',
            'tap_min' : -9,
            'tap_max' : 9,
            'tap_neutral' : 0,
            'tap_step_percent' : 1.5,
            'tap_step_degree' : 0,
            'tap_phase_shifter' : False,
            'tap_dependency_table' : False,
            'tap_changer_type' : 'Ratio'
        }
        #pp.create_std_type(net, trafo_std_type, name = 'OLTC', element= 'trafo')

        #pp.create_transformer(net, hv_bus = 3, lv_bus=8, std_type = 'OLTC', tap_pos = 0)
    if use_sgen:
        pp.create_sgen(net,6,p_mw = 150, q_mvar = 50)
        pp.create_sgen(net,5,p_mw = 120, q_mvar = 100) 
        net.sgen.sn_mva = (net.sgen.p_mw**2 + net.sgen.q_mvar**2)**(1/2)*1.35

    return net

def ieee118():
    net = pn.case118()

    #for i in net.gen.index:
    #    net.gen.drop(i, inplace = True)
    #for i in net.poly_cost.index:
    #    net.poly_cost.drop(i, inplace = True)
    
    net.trafo.tap_min = -9
    net.trafo.tap_max = 9
    net.trafo.tap_step_percent = 1.5
    net.trafo.tap_pos = 0
    net.trafo.tap_side = 'lv'
    net.trafo.tap_neutral = 0
    return net

def ieee33(ppf = False, opf = False, vvc = False, use_sgen = False):
    net = pp.networks.case33bw()
    net.name = 'ieee33'

    Dropindex = [32, 33, 34, 35, 36] # 그림의 index보다 하나 작음
    ShuntBus = [17, 32]
    TrafoBus = [[6, 7]]#,[25,26]] #[6, 7],[25,26]

    net.line.drop(Dropindex, inplace= True)
    for bus in ShuntBus:
        pp.create_shunt(net, bus, q_mvar = -0.5)

    trafo_std_type = {
        'sn_mva' : 10,
        'vn_hv_kv' : net.bus.vn_kv[0],
        'vn_lv_kv' : net.bus.vn_kv[0],
        'vk_percent' : 0.8,
        'vkr_percent' : 0.2,
        'pfe_kw' : 0,
        'i0_percent' :0,
        'shift_degree' : 0,
        'tap_side' : 'lv',
        'tap_min' : -9,
        'tap_max' : 9,
        'tap_neutral' : 0,
        'tap_step_percent' : 0,
        'tap_step_degree' : 0,
        'tap_phase_shifter' : False,
        'tap_dependency_table' : False,
        'tap_changer_type' : 'Ratio'

    }

    pp.create_std_type(net, trafo_std_type, name = 'OLTC', element= 'trafo')


    for idx in net.line.index:
        pp.create_switch(net, net.line.from_bus[idx], idx, et = 'l', closed = True, name = 'l' + str(idx))

    # tie lines
    #pp.create_switches(net, [5, 11, 24, 28], [20, 21, 28, 11], et = 'b', closed = False, name = ['b0', 'b1', 'b2', 'b3'])
    pp.create_switches(net, [24,7,11], [28,20,21], et = 'b', closed = False, name = ['b0', 'b1', 'b2'])

    for i in net.switch.index:
        net.switch.drop(i,inplace = True)
    if ppf:
        for i in net.gen.index:
            net.gen.drop(i, inplace = True)
        for i in net.poly_cost.index:
            net.poly_cost.drop(i, inplace = True)
        for i in net.trafo.index:
            net.trafo.drop(i, inplace = True)
        for i in net.shunt.index:
            net.shunt.drop(i, inplace = True)

    if opf: # need to add generator (but wont be used)
        net.gen.loc[0, 'slack'] = True
        for i in net.ext_grid.index:
            net.ext_grid.drop(i, inplace = True)        
        for i in net.poly_cost[net.poly_cost.et == 'ext_grid'].index:
            net.poly_cost.drop(i, inplace = True)
        for i in net.trafo.index:
            net.trafo.drop(i, inplace = True)
        for i in net.shunt.index:
            net.shunt.drop(i, inplace = True)

    if vvc:
        for i in net.gen.index:
            net.gen.drop(i, inplace = True)
        for i in net.poly_cost.index:
            net.poly_cost.drop(i, inplace = True)
        # switch 정의하기 전에 정의해라 무조건
        for trafo in TrafoBus:
            for idx in net.line.index:
                if net.line.from_bus[idx] == trafo[0] and net.line.to_bus[idx] == trafo[1]:
                    net.line.drop(idx, inplace = True)               
                    tfindex = pp.create_transformer(net, trafo[0], trafo[1], 'OLTC', tap_pos = 0)

                elif net.line.to_bus[idx] == trafo[0] and net.line.from_bus[idx] == trafo[1]:
                    net.line.drop(idx, inplace = True)            
                    tfindex = pp.create_transformer(net, trafo[1], trafo[0], 'OLTC', tap_pos = 0)

        #for i in net.shunt.index:
        #    net.shunt.drop(i, inplace = True)

#            pp.create_switch(net, trafo[0], tfindex, et = 't', closed = True)

    if use_sgen:
        pp.create_sgen(net, 24,p_mw = 0.1, q_mvar = 0.2, sn_mva = 1, type='PV')
        pp.create_sgen(net, 32,p_mw = 0.1, q_mvar = 0.2, sn_mva = 1, type='PV')
        pp.create_sgen(net, 17,p_mw = 0.1, q_mvar = 0.2, sn_mva = 0.7, type='PV')
        pp.create_sgen(net, 21,p_mw = 0.1, q_mvar = 0.2, sn_mva = 0.7, type='PV')
    
    
    net.line.max_i_ka = 0.02
    net.load.p_mw = net.load.p_mw*5*0.17
    net.load.q_mvar = net.load.q_mvar*5*0.17

    #multiplier = [1.8*3, 1.1*2.5, 0.7, 0.7, 0.7, 0.4, 0.4, 0.4, 0.3, 0.3, 1, 1, 1.5, 1.4, 1.4, 1.1, 0.95, 0.95, 1, 1, 0.6*2.5, 0.6*4, 0.4*4.8, 0.8*2, 0.8, 0.8, 0.8, 0.8, 1, 1, 1.1, 1]
    #ika = deepcopy(net.line.max_i_ka)
    #for idx1, idx in enumerate(ika.index):
    #    net.line.loc[idx, 'max_i_ka'] = ika[idx]*multiplier[idx1]*3/1.8
    net.sgen.sn_mva = (net.sgen.p_mw**2 + net.sgen.q_mvar**2)**(1/2)*1.35

    return net

def ieee123(ppf = False, opf = False, vvc = False, use_sgen = False):
    if int(ppf) + int(opf) + int(vvc) >= 2 or int(ppf) + int(opf) + int(vvc) == 0:
        raise('Only one among pp, opf, vvc needs to be set True')
    
    previous_dir = os.getcwd()
    os.chdir(os.path.dirname(os.path.realpath(__file__)))
    line_data = pd.read_excel('IEEE123Node/line data.xls')
    line_config = pd.read_excel('IEEE123Node/line config.xls')

    net = pp.create_empty_network(name = 'IEEE123')
    for i in range(len(line_config)):
        pp.create_std_type(net, dict(line_config[line_config.columns[1:]].iloc[i]), name = line_config.name[i], element = 'line')

    bus_set = set(line_data[line_data.columns[0]]).union(set(line_data[line_data.columns[1]]))
    for i in bus_set:
        pp.create_bus(net, vn_kv = 4.16, name = f'{i}', id = i)
        
    ## external grid
    ext1 = pp.create_bus(net, vn_kv = 4.16, name = '150', id = 150)
    pp.create_ext_grid(net, ext1, vm_pu = 1.05, s_sc_max_mva = 3)

    ## singlebus
    sb1 = pp.create_bus(net, vn_kv = 4.16, name = '610', id = 610)
    sb2 = pp.create_bus(net, vn_kv = 4.16, name = '350', id = 350)

    bus_to_index = pd.Series(net.bus.index, index = net.bus.name)         ## 이는 그림 상 나온 bus 이름 (int)을 받으면 그를 pandapower 상 index로 변환
    bus_to_index.index.name = None
    bus_to_index.index = bus_to_index.index.map(int)                

    ## capacitor
    capacitor_data = pd.read_excel('IEEE123Node/capacitor data.xls')
    for i in range(len(capacitor_data)):
        pp.create_shunt(net, bus_to_index.loc[capacitor_data[capacitor_data.columns[0]][i]], capacitor_data[capacitor_data.columns[1]][i]/1000.)

    ## line
    for i in range(len(line_data)): 
        pp.create_line(net, bus_to_index.loc[line_data[line_data.columns[0]][i]], bus_to_index.loc[line_data[line_data.columns[1]][i]], length_km = line_data[line_data.columns[2]][i]*0.0003048, std_type = line_data[line_data.columns[3]][i], name = f'{i+1}')
        pp.create_switch(net, bus_to_index.loc[line_data[line_data.columns[0]][i]], i, et = "l", closed = True) # 모든 line에 switch 연결
    pp.create_line(net, bus_to_index.loc[350], bus_to_index.loc[300], length_km = 0.0003, std_type = line_data[line_data.columns[3]][0], name = f'{i+2}')
    pp.create_line(net, bus_to_index.loc[61], bus_to_index.loc[610], length_km = 0.0003, std_type = line_data[line_data.columns[3]][0], name = f'{i+3}')

    load_data = pd.read_excel('IEEE123Node/spot loads data.xls')
    load_data['P'] = (load_data['P1'] + load_data['P2'] + load_data['P3'])/1000.
    load_data['Q'] = (load_data['Q1'] + load_data['Q2'] + load_data['Q3'])/1000.
    for i in range(len(load_data)):
        pp.create_load(net, bus = bus_to_index[load_data[load_data.columns[0]][i]], p_mw = load_data['P'][i], q_mvar = load_data['Q'][i], sn_mva = sqrt(load_data['P'][i]**2 + load_data['Q'][i]**2), name = f'{i+1}')
 
    switch_data = pd.read_excel('IEEE123Node/switch data.xls')
    for i in range(len(switch_data)):
        from_bus = switch_data[switch_data.columns[0]][i]
        to_bus = switch_data[switch_data.columns[1]][i]
        closed = switch_data[switch_data.columns[2]][i]
        if closed == "closed":
            closed = True
        else:
            closed = False
        pp.create_switch(net, bus_to_index.loc[from_bus], bus_to_index.loc[to_bus], et = "b", closed = closed)
    #pp.create_switch(net, bus_to_index.loc[33], bus_to_index.loc[233], et = "b", closed = False) ## excel에 없음
    #pp.create_switch(net, bus_to_index.loc[13], 12, et = "l", closed = True) ## excel에 없음
    pp.create_switch(net, bus_to_index.loc[39], bus_to_index.loc[66], et = "b", closed = False) ## excel에 없음
    pp.create_switch(net, bus_to_index.loc[17], bus_to_index.loc[96], et = "b", closed = False) ## excel에 없음
    #pp.create_switch(net, bus_to_index.loc[76], 76, et = "l", closed = True) ## excel에 없음

    ## 참고 : 부하의 유효전력 합은 3.49 mw, 무효전력 합은 1.92 mvar

    sgen_data = pd.read_excel('IEEE123Node/sgen data.xls')
  # p_mw를 크게 올리면 Q error 커지긴 함
    if use_sgen:
        for i in range(len(sgen_data)):
            pp.create_sgen(net, bus_to_index[sgen_data[sgen_data.columns[0]][i]], p_mw=sgen_data[sgen_data.columns[1]][i]/3, q_mvar=sgen_data[sgen_data.columns[2]][i]/3, sn_mva=sgen_data[sgen_data.columns[3]][i]/3, scaling=sgen_data[sgen_data.columns[4]][i], name=sgen_data[sgen_data.columns[5]][i], type=sgen_data[sgen_data.columns[6]][i])
        net.sgen.sn_mva = (net.sgen.p_mw**2 + net.sgen.q_mvar**2)**(1/2)*1.35


    print("IEEE network is loaded")
    net.sgen.loc[net.sgen.type == 'WT', 'type'] = 'WP'
    os.chdir(previous_dir)

    net.line.max_i_ka = net.line.max_i_ka*0.7

    net.sgen.p_mw = net.sgen.p_mw
    net.sgen.q_mvar = net.sgen.q_mvar
    net.shunt.in_service =  False
    net.shunt.q_mvar = 0.1
    net.line.c_nf_per_km = 0
    net.line.g_us_per_km = 0
    
    if ppf:
        for i in net.gen.index:
            net.gen.drop(i, inplace = True)
        for i in net.poly_cost.index:
            net.poly_cost.drop(i, inplace = True)
        for i in net.shunt.index:
            net.shunt.drop(i, inplace = True)

    if opf:
        pp.create_gen(net,84, p_mw = 10, vm_pu = 1, max_p_mw= 10, min_p_mw = 0.1, max_q_mvar=10, min_q_mvar=-10)
        pp.create_gen(net,123, p_mw = 10, vm_pu = 1, max_p_mw= 10, min_p_mw = 0.1, max_q_mvar=10, min_q_mvar=-10)
        pp.create_gen(net,41, p_mw = 10, vm_pu = 1, max_p_mw= 10, min_p_mw = 0.1, max_q_mvar=10, min_q_mvar=-10)
        net.gen.loc[0, 'slack'] = True
        for i in net.ext_grid.index:
            net.ext_grid.drop(i, inplace = True)        
        for i in net.poly_cost[net.poly_cost.et == 'ext_grid'].index:
            net.poly_cost.drop(i, inplace = True)
        pp.create_poly_cost(net, 0, et = 'gen', cp0_eur = 400, cp1_eur_per_mw= 0.8, cp2_eur_per_mw2= 0.01)
        pp.create_poly_cost(net, 1, et = 'gen', cp0_eur = 200, cp1_eur_per_mw= 1.2, cp2_eur_per_mw2= 0.075)
        pp.create_poly_cost(net, 2, et = 'gen', cp0_eur = 300, cp1_eur_per_mw= 1, cp2_eur_per_mw2= 0.012)
      
    if vvc:
        for i in net.gen.index:
            net.gen.drop(i, inplace = True)
        for i in net.poly_cost.index:
            net.poly_cost.drop(i, inplace = True)   
        net.switch.drop(122, inplace = True)
        rng2 = np.random.default_rng(seed=42)

        trafo_std_type = {
            'sn_mva' : 10,
            'vn_hv_kv' : net.bus.vn_kv[123],
            'vn_lv_kv' : net.bus.vn_kv[115],
            'vk_percent' : 0.08,
            'vkr_percent' : 0.01, # 바꿈
            'pfe_kw' : 0,
            'i0_percent' : 0,
            'shift_degree' : 0,
            'tap_side' : 'lv',
            'tap_min' : -9,
            'tap_max' : 9,
            'tap_neutral' : 0,
            'tap_step_percent' : 1.5,
            'tap_step_degree' : 0,
            'tap_phase_shifter' : False,
            'oltc' : True,
            'tap_dependency_table' : False,
            'tap_changer_type' : 'Ratio'
        }
        pp.create_std_type(net, trafo_std_type, name = 'OLTC', element= 'trafo')
        pp.create_transformer(net, hv_bus = 123, lv_bus=115, std_type = 'OLTC', tap_pos = 0)

        number_of_trafo = 4
        trafo_line = rng2.choice(list(net.line.index), 10, replace = False)
        trafo_line = trafo_line[:number_of_trafo]
        trafo_line.sort()

        for linee in trafo_line:
            bus1, bus2 = net.line.from_bus[linee], net.line.to_bus[linee]
            net.line.drop(linee, inplace= True)
            pp.create_transformer(net, hv_bus = bus1, lv_bus=bus2, std_type = 'OLTC', tap_pos = 0)
            net.switch.drop(net.switch.index[net.switch.element == linee][0], inplace = True)


    net.bus['max_vm_pu'] = 1.1
    net.bus['min_vm_pu'] = 0
    net.name = 'ieee123'
    return net

def ieee9_mesh():
    net = pn.case9() # Example network
    pp.create_switch(net, 0, 0, et = 'l', closed = True)
    pp.create_switch(net, 3, 1, et = 'l', closed = True)
    pp.create_switch(net, 0, 7, et = 'b', closed= False)
    pp.create_load(net,1, p_mw = 100, q_mvar = 35)
    pp.create_load(net,2, p_mw = 100, q_mvar = 35)
    net.line.length_km = 0.1
    net.line.c_nf_per_km = 0
    net.line.r_ohm_per_km = net.line.r_ohm_per_km +1
    pp.create_sgen(net,6,p_mw = 150, q_mvar = 50)
    pp.create_sgen(net,5,p_mw = 120, q_mvar = 100) 
    net.sgen.sn_mva = (net.sgen.p_mw**2 + net.sgen.q_mvar**2)**(1/2)*1.35
    net.ext_grid.drop(0, inplace = True)
    net.poly_cost.drop(0, inplace = True)
    net.gen.loc[0, 'slack'] = True
    return net

def ieee9_radial():
    net = pn.case9() # Example network
    pp.create_switch(net, 0, 0, et = 'l', closed = True)
    pp.create_switch(net, 3, 1, et = 'l', closed = True)
    pp.create_switch(net, 0, 7, et = 'b', closed= False)
    net.line.drop(4, inplace = True)
    pp.create_load(net,1, p_mw = 100, q_mvar = 35)
    pp.create_load(net,2, p_mw = 100, q_mvar = 35)
    net.line.length_km = 0.1
    net.line.c_nf_per_km = 0
    net.line.r_ohm_per_km = net.line.r_ohm_per_km +1
    pp.create_sgen(net,6,p_mw = 150, q_mvar = 50)
    pp.create_sgen(net,5,p_mw = 120, q_mvar = 100) 
    net.ext_grid.drop(0, inplace = True)
    net.poly_cost.drop(0, inplace = True)

    net.sgen.sn_mva = (net.sgen.p_mw**2 + net.sgen.q_mvar**2)**(1/2)*1.35
    net.gen.loc[0, 'slack'] = True
    return net

def ieee123_trafo(use_sgen = True):
    previous_dir = os.getcwd()
    os.chdir(os.path.dirname(os.path.realpath(__file__)))
    line_data = pd.read_excel('IEEE123Node/line data.xls')
    line_config = pd.read_excel('IEEE123Node/line config.xls')

    net = pp.create_empty_network(name = 'IEEE123')
    for i in range(len(line_config)):
        pp.create_std_type(net, dict(line_config[line_config.columns[1:]].iloc[i]), name = line_config.name[i], element = 'line')

    bus_set = set(line_data[line_data.columns[0]]).union(set(line_data[line_data.columns[1]]))
    for i in bus_set:
        pp.create_bus(net, vn_kv = 4.16, name = f'{i}', id = i)
        
    ## external grid
    ext1 = pp.create_bus(net, vn_kv = 4.16, name = '150', id = 150)
    #ext2 = pp.create_bus(net, vn_kv = 4.16, name = '233', id = 233)
#    ext3 = pp.create_bus(net, vn_kv = 4.16, name = '451', id = 451)
#    ext4 = pp.create_bus(net, vn_kv = 4.16, name = '195', id = 195)

    pp.create_ext_grid(net, ext1, vm_pu = 1.03, s_sc_max_mva = 3)
    #pp.create_ext_grid(net, ext2, vm_pu = 1.03, s_sc_max_mva = 3)
#    pp.create_ext_grid(net, ext3, s_sc_max_mva = 5)
#    pp.create_ext_grid(net, ext4, s_sc_max_mva = 5)

    ## singlebus
    sb1 = pp.create_bus(net, vn_kv = 4.16, name = '610', id = 610)
    sb2 = pp.create_bus(net, vn_kv = 4.16, name = '350', id = 350)

    bus_to_index = pd.Series(net.bus.index, index = net.bus.name)         ## 이는 그림 상 나온 bus 이름 (int)을 받으면 그를 pandapower 상 index로 변환
    bus_to_index.index.name = None
    bus_to_index.index = bus_to_index.index.map(int)                

    ## capacitor
    
    capacitor_data = pd.read_excel('IEEE123Node/capacitor data.xls')
    for i in range(len(capacitor_data)):
        pp.create_shunt(net, bus_to_index.loc[capacitor_data[capacitor_data.columns[0]][i]], capacitor_data[capacitor_data.columns[1]][i]/1000.)

    ## line
    for i in range(len(line_data)): 
        pp.create_line(net, bus_to_index.loc[line_data[line_data.columns[0]][i]], bus_to_index.loc[line_data[line_data.columns[1]][i]], length_km = line_data[line_data.columns[2]][i]*0.0003048, std_type = line_data[line_data.columns[3]][i], name = f'{i+1}')
        pp.create_switch(net, bus_to_index.loc[line_data[line_data.columns[0]][i]], i, et = "l", closed = True) # 모든 line에 switch 연결
    pp.create_line(net, bus_to_index.loc[350], bus_to_index.loc[300], length_km = 0.0003, std_type = line_data[line_data.columns[3]][0], name = f'{i+2}')
    pp.create_line(net, bus_to_index.loc[61], bus_to_index.loc[610], length_km = 0.0003, std_type = line_data[line_data.columns[3]][0], name = f'{i+3}')



    load_data = pd.read_excel('IEEE123Node/spot loads data.xls')
    load_data['P'] = (load_data['P1'] + load_data['P2'] + load_data['P3'])/1000.
    load_data['Q'] = (load_data['Q1'] + load_data['Q2'] + load_data['Q3'])/1000.
    for i in range(len(load_data)):
        pp.create_load(net, bus = bus_to_index[load_data[load_data.columns[0]][i]], p_mw = load_data['P'][i], q_mvar = load_data['Q'][i], sn_mva = sqrt(load_data['P'][i]**2 + load_data['Q'][i]**2), name = f'{i+1}')
 

    switch_data = pd.read_excel('IEEE123Node/switch data.xls')
    for i in range(len(switch_data)):
        from_bus = switch_data[switch_data.columns[0]][i]
        to_bus = switch_data[switch_data.columns[1]][i]
        closed = switch_data[switch_data.columns[2]][i]
        if closed == "closed":
            closed = True
        else:
            closed = False
        pp.create_switch(net, bus_to_index.loc[from_bus], bus_to_index.loc[to_bus], et = "b", closed = closed)
    #pp.create_switch(net, bus_to_index.loc[33], bus_to_index.loc[233], et = "b", closed = False) ## excel에 없음
    #pp.create_switch(net, bus_to_index.loc[13], 12, et = "l", closed = True) ## excel에 없음
    pp.create_switch(net, bus_to_index.loc[39], bus_to_index.loc[66], et = "b", closed = False) ## excel에 없음
    pp.create_switch(net, bus_to_index.loc[17], bus_to_index.loc[96], et = "b", closed = False) ## excel에 없음
    #pp.create_switch(net, bus_to_index.loc[76], 76, et = "l", closed = True) ## excel에 없음

    ## 참고 : 부하의 유효전력 합은 3.49 mw, 무효전력 합은 1.92 mvar

    sgen_data = pd.read_excel('IEEE123Node/sgen data.xls')

    if use_sgen is True:
        for i in range(len(sgen_data)):
            pp.create_sgen(net, bus_to_index[sgen_data[sgen_data.columns[0]][i]], p_mw=sgen_data[sgen_data.columns[1]][i], q_mvar=sgen_data[sgen_data.columns[2]][i], sn_mva=sgen_data[sgen_data.columns[3]][i], scaling=sgen_data[sgen_data.columns[4]][i], name=sgen_data[sgen_data.columns[5]][i], type=sgen_data[sgen_data.columns[6]][i])

    print("IEEE network is loaded")
    net.sgen.type[net.sgen.type == 'WT'] = 'WP'
    os.chdir(previous_dir)

    net.line.max_i_ka = net.line.max_i_ka*0.7

    net.load.p_mw = net.load.p_mw*2
    net.sgen.p_mw = net.sgen.p_mw
    net.sgen.q_mvar = net.sgen.q_mvar
    net.sgen.sn_mva = net.sgen.p_mw*1
    net.shunt.in_service =  False
    net.shunt.q_mvar = 0.001
    net.line.c_nf_per_km = 0
    net.line.g_us_per_km = 0

    net.switch.drop(122, inplace = True)

    trafo_std_type = {
        'sn_mva' : 10,
        'vn_hv_kv' : net.bus.vn_kv[123],
        'vn_lv_kv' : net.bus.vn_kv[115],
        'vk_percent' : 8,
        'vkr_percent' : 0,
        'pfe_kw' : 0,
        'i0_percent' : 0,
        'shift_degree' : 0,
        'tap_side' : 'lv',
        'tap_min' : -9,
        'tap_max' : 9,
        'tap_neutral' : 0,
        'tap_step_percent' : 1.5,
        'tap_step_degree' : 0,
        'tap_phase_shifter' : False,
        'oltc' : True,
        'tap_dependency_table' : False,
        'tap_changer_type' : 'Ratio'
    }
    pp.create_std_type(net, trafo_std_type, name = 'OLTC', element= 'trafo')

    pp.create_transformer(net, hv_bus = 123, lv_bus=115, std_type = 'OLTC', tap_pos = 0)

    return net

def ieee8500(ppf = True, use_sgen = True, vvc = False , seed = 42, number_of_trafo = 20):
    if os.path.exists('ieee8500.p'):
        net = pp.from_pickle('ieee8500.p')
        number_of_sgen = 1213 #1213
        rng = np.random.default_rng(seed=seed)
        #sgen_bus = rng.choice(range(3699, 4875), number_of_sgen)
        sgen_bus = rng.choice(4875, number_of_sgen, replace = False)
        
        sgen_bus.sort()
        load_pw_mean, load_qvar_mean = net.load.p_mw.mean(), net.load.q_mvar.mean()
        net.load.p_mw = net.load.p_mw*1.5
        if use_sgen:
            for bus in sgen_bus:
                p, q = load_pw_mean*(rng.random())*1, load_qvar_mean*(1+ (rng.random()-0.5))*0.2
                #p, q = load_pw_mean*(rng.random())*0.2, load_qvar_mean*(1+ (rng.random()-0.5))*0.2
                pp.create_sgen(net, bus, p_mw = p, q_mvar = q, sn_mva = p*4 + q*1.35)

        if vvc:
            number_of_trafo = number_of_trafo # 115
            trafo_line = rng.choice(list(net.line.index), 115, replace = False)
            trafo_line = trafo_line[:number_of_trafo]
            trafo_line.sort()
            trafo_std_type = {
                'sn_mva' : 10,
                'vn_hv_kv' : net.bus.vn_kv[0],
                'vn_lv_kv' : net.bus.vn_kv[0],
                'vk_percent' : 2, 
                'vkr_percent' : 1,
                'pfe_kw' : 0,
                'i0_percent' : 0,
                'shift_degree' : 0,
                'tap_side' : 'lv',
                'tap_min' : -9,
                'tap_max' : 9,
                'tap_neutral' : 0,
                'tap_step_percent' : 1,
                'tap_step_degree' : 0,
                'tap_phase_shifter' : False,
                'oltc' : True,
                'tap_dependency_table' : False,
                'tap_changer_type' : 'Ratio'
            }
            pp.create_std_type(net, trafo_std_type, name = 'OLTC', element= 'trafo')

            for linee in trafo_line:
                bus1, bus2 = net.line.from_bus[linee], net.line.to_bus[linee]
                net.line.drop(linee, inplace= True)
                pp.create_transformer(net, hv_bus = bus1, lv_bus=bus2, std_type = 'OLTC', tap_pos = 0)

            #number_of_shunt = 738 
            shunt_bus = rng.choice(list(net.bus.index), 12, replace = False)
            for bus in shunt_bus:
                pp.create_shunt(net, bus, q_mvar = 0.01)
        net.name = 'ieee8500'
        return net
    else:
        return create_ieee8500()

def create_ieee8500():
    DSSObject = comtypes.client.CreateObject('OpenDSSEngine.DSS')

    DSSObject.Text.Command = 'redirect "C:/Program Files/OpenDSS/IEEETestCases/8500-Node/Run_8500Node.dss"'

    DSSText = DSSObject.Text
    DSSCircuit = DSSObject.ActiveCircuit
    DSSActiveBus = DSSCircuit.ActiveBus
    DSSSolution = DSSCircuit.Solution
    DSSLines = DSSCircuit.Lines
    DSSLoads = DSSCircuit.Loads
    DSSBuses = DSSCircuit.ActiveBus
    DSSActiveElement = DSSCircuit.ActiveCktElement

    DSSLoads.kW
    DSSLoads.kva
    #dss.ActiveCircuit.ActiveBus.kVBase
    #dss.ActiveCircuit.Vsources.AllNames
    # https://dss-extensions.org/python_apis.html < 여기 제일 아래에 DSS-Python List


    DSSObject.ActiveCircuit.Vsources.AllNames
    DSSObject.ActiveCircuit.Transformers.AllNames
    DSSObject.ActiveCircuit.Loads.AllNames
    DSSObject.ActiveCircuit.Lines.AllNames
    DSSObject.ActiveCircuit.AllBusNames
    DSSActiveBus.AllPDEatBus
    DSSActiveBus.AllPCEatBus

    net = pp.create_empty_network()


    BusNameToIndex = dict()
    BusIndexToName = dict()

    for idx, name in enumerate(DSSObject.ActiveCircuit.AllBusNames):
        BusNameToIndex[name] = idx
        BusIndexToName[idx+1] = name
        pp.create_bus(net, 66.4, index = idx)

    check_if_multiple_line = dict()
    DSSLines = DSSCircuit.Lines

    idx = DSSLines.First
    while idx != 0:

        Bus1idx, Bus2idx = BusNameToIndex[DSSLines.Bus1.split(".")[0]], BusNameToIndex[DSSLines.Bus2.split(".")[0]]
        if (Bus1idx, Bus2idx) in check_if_multiple_line.keys():
            print('multiple lines')
        else:
            r, x, c, i, dist = DSSLines.R1, DSSLines.X1, DSSLines.C1, DSSLines.EmergAmps, DSSLines.length
            pp.create_line_from_parameters(net, Bus1idx, Bus2idx, dist, r, x, 0, i/1000.)
            check_if_multiple_line[(Bus1idx, Bus2idx)], check_if_multiple_line[(Bus2idx, Bus1idx)] = True, True
        idx = DSSLines.Next


    DSSObject.ActiveCircuit.Transformers.AllNames

    idx = DSSObject.ActiveCircuit.Transformers.First
    while idx != 0:
        Bus1idx, Bus2idx = BusNameToIndex[DSSObject.ActiveCircuit.ActiveElement.BusNames[0].split(".")[0]], BusNameToIndex[DSSObject.ActiveCircuit.ActiveElement.BusNames[1].split(".")[0]] # 어떤 ActiveElement건 Busname 반환
        if (Bus1idx, Bus2idx) in check_if_multiple_line.keys():
            print('multiple lines')
        else:
            r, x, c, i, dist = 0.001, 0.01, 0., 600, 0.001
            pp.create_line_from_parameters(net, Bus1idx, Bus2idx, dist, r, x, 0, i/1000.)
            check_if_multiple_line[(Bus1idx, Bus2idx)], check_if_multiple_line[(Bus2idx, Bus1idx)] = True, True
        idx = DSSObject.ActiveCircuit.Transformers.Next

    idx = DSSObject.ActiveCircuit.Vsources.First
    while idx != 0:
        Bus1idx = BusNameToIndex[DSSObject.ActiveCircuit.ActiveElement.BusNames[0].split('.')[0]]
        v_pu = DSSObject.ActiveCircuit.Vsources.pu
        pp.create_ext_grid(net, Bus1idx, v_pu)
        idx = DSSObject.ActiveCircuit.Vsources.Next

    DSSObject.ActiveCircuit.Reactors.Count
    idx = DSSObject.ActiveCircuit.Reactors.First
    while idx != 0:
        Bus1idx, Bus2idx = BusNameToIndex[DSSObject.ActiveCircuit.ActiveElement.BusNames[0].split(".")[0]], BusNameToIndex[DSSObject.ActiveCircuit.ActiveElement.BusNames[1].split(".")[0]] # 어떤 ActiveElement건 Busname 반환
        if (Bus1idx, Bus2idx) in check_if_multiple_line.keys():
            print('multiple lines')
        else:
            r, x, c, i, dist = 0.001, 0.01, 0., 600, 0.001
            pp.create_line_from_parameters(net, Bus1idx, Bus2idx, dist, r, x, 0, i/1000.)
            check_if_multiple_line[(Bus1idx, Bus2idx)], check_if_multiple_line[(Bus2idx, Bus1idx)] = True, True
        idx = DSSObject.ActiveCircuit.Reactors.Next
    idx = DSSObject.ActiveCircuit.Loads.First
    while idx != 0:
        Bus1idx = BusNameToIndex[DSSObject.ActiveCircuit.ActiveElement.BusNames[0].split(".")[0]]
        p_kw, q_kvar = DSSObject.ActiveCircuit.Loads.kW, DSSObject.ActiveCircuit.Loads.kvar
        pp.create_load(net, Bus1idx, p_kw/1000., q_kvar/1000.)
        idx = DSSObject.ActiveCircuit.Loads.Next

    net.line.length_km = 40*net.line.length_km 
    net.name = 'ieee8500'
    pp.to_pickle(net, 'ieee8500.p')
    return net
