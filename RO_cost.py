# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 02:38:36 2020
RO COST MODEL
Following a structure similar to VAGMD_cost.py
@author: adama

Battery cost
"""



import numpy as np
import math
from desalinationcosts import CR_factor


class RO_cost(object):
    def __init__(self,
                 Capacity = [1000], # Desalination plant capacity (m3/day)
                 Prod = 328500, # Annual permeate production (m3)
                 Area = 40.8 , # Membrane area (m2)
                 fuel_usage = 30, # Total fuel usage (%)
                 # Pflux = 1.1375,  # Permeate flow per module (m3/h/module)
                 num_modules = 100, # Number of module (from design model) 
                 membrane_cost=30, # cost per unit area of membrane (USD/m2)
                 pressure_vessel_cost= 1000, # cost per pressure vessel (USD/vessel)-- just guessing 1000 to fill value for now: 2/25/2020 
#                 FFR  = 17, # Feed flow rate per module (m3/h/module)
#                 GOR = 10.475,  # Gained output ratio
                # downtime = 0.1, # Yearly downtime of the plant (ratio)
                 yrs = 20, # Expected plant lifetime
                 int_rate = 0.04 , # Average interest rate
                 coe = 0.05,  # Unit cost of electricity ($/kWh)
                 sam_coe = 0.02,
                 ave_grid_coe = 0.04, # average cost of grid energy used ($/kWh) - calcualted from the simulation model
                 monthly_grid_coe = [], 
                 solar_coe =  0,
                 chem_cost=0.05, # specific chemical cost ($/m3)
                 labor_cost=0.1,  # specific labor cost ($/m3)
                 rep_rate=5,    # membrane replacement rate (yrs)
                 equip_cost_method='general', # Option 1:'general' capex breakdown by cost factor; Option 2: 'specify' equipment costs
                 sec= 2.5, # specific energy consumption (kwh/m3)
                 HP_pump_pressure=60, # high pressure pump discharge pressure (bar)
                 HP_pump_flowrate=41.67,   # high pressure pump flowrate (m3/hr)
                 BP_pump_pressure=10.5, # booster pump discharge pressure (bar)
                 BP_pump_flowrate=62.5, # booster pump flowrate (m3/hr)
                 ERD_flowrate=62.5,     # ERD flowrate set to brine flowrate(m3/hr)
                 ERD_pressure=49  ,      # ERD pressure set to brine pressure entering (bar)  
                 disposal_cost=0.03,    # specific waste disposal cost($/m3)
                 IX_cost = 0.15, # ion-exchange cost ($/m3)
                 insurance = 0.5, # insurance (percentage of CAPEX)
                 #coh = 0.01 , # Unit cost of heat ($/kWh)
                 #sam_coh = 0.02, # Unit cost of heat from SAM ($/kWh)
                 #solar_inlet = 85, # Solar field inlet temperature
                 #solar_outlet = 95, # Solar field outlet temperature
                 #HX_eff = 0.85, # Heat exchanger efficiency
                 #cost_module_re = 0.220 , # Cost of module replacement ($/m3)
                 unit_capex = [1000],
                 unit_capex_main= 1647,  # total EPC cost, USD/(m3/day)
                 unit_capex_passes=668,  # total EPC cost, USD/(m3/day)              
                 downtime = 10, # downtime percentage
                 cost_storage = 150 , # Cost of battery ($/kWh)
                 storage_cap = 0, # Capacity of battery (kWh)
                 monthly_input = True,
                 monthly_elec_curtailed = [],
                 monthly_grid_consumed = [],

                 brackish_depth = 0, # ft
                 well_cost = 650, # $/ft per 3.25 MGD capacity (~15000 m3/day)
                 RO_rr = 50,
                 ):
        self.HP_pump_pressure=HP_pump_pressure
        self.HP_pump_flowrate=HP_pump_flowrate
        self.BP_pump_pressure=BP_pump_pressure
        self.BP_pump_flowrate=BP_pump_flowrate
        self.ERD_flowrate=ERD_flowrate
        self.ERD_pressure=ERD_pressure
        self.ann_prod=Prod * (1-downtime /100)
        self.chem_cost=chem_cost 
        self.labor_cost=labor_cost 
        self.operation_hour = 24 #* (1-downtime) # Average daily operation hour (h/day)
        # self.Pflux = Pflux
        self.Area = Area
        # self.PF_module = self.Pflux * self.Area
        self.num_modules = num_modules
        self.total_area = self.num_modules*self.Area
        self.membrane_cost=membrane_cost
        self.pressure_vessel_cost=pressure_vessel_cost
        self.SEC=sec
        self.capacity=Capacity
        self.equip_cost_method=equip_cost_method
        self.monthly_grid_coe = monthly_grid_coe
        self.monthly_input = monthly_input
        self.monthly_elec_curtailed = monthly_elec_curtailed
        self.monthly_grid_consumed = monthly_grid_consumed

        
        # calculate membrane replacement cost
        rep_yr = rep_rate
        self.replacement_rate = 0
        while rep_yr < yrs:
            self.replacement_rate += 1 / (1+int_rate) ** rep_yr
            rep_yr += rep_rate              
        self.replacement_rate *= CR_factor(yrs,int_rate)
        self.membrane_replacement_cost=self.membrane_cost*self.total_area*self.replacement_rate/self.ann_prod
        
        
        self.disposal_cost=disposal_cost 
#        self.HX_eff = HX_eff
        
        self.unit_capex = unit_capex
        self.unit_capex_main = unit_capex_main
        self.unit_capex_passes = unit_capex_passes
        self.downtime = downtime / 100
#        self.th_module = th_module
#        self.FFR = FFR
        self.IX_cost = IX_cost
        self.insurance = insurance / 100
        self.coe = coe
        if solar_coe is None:
            self.sam_coe = sam_coe
        else:
            self.sam_coe = float(solar_coe)

        self.fuel_usage = fuel_usage / 100
#        self.cost_module_re = cost_module_re
        self.yrs = yrs
        self.int_rate = int_rate
        self.cost_storage = cost_storage
        self.storage_cap = storage_cap
        self.ave_grid_coe = ave_grid_coe

        self.brackish_depth = brackish_depth
        self.well_cost = well_cost
        self.RO_rr = RO_rr 
                
    def lcow(self):
        if self.equip_cost_method=='specify':
            self.total_module_cost = self.total_area*self.membrane_cost + self.NV*self.pressure_vessel_cost
            self.HPpump_cost=53*self.HP_pump_flowrate*self.HP_pump_pressure
#            self.test=pump_cost(HP_pump_pressure,HP_pump_flowrate)
            self.BPpump_cost=53*self.BP_pump_flowrate*self.BP_pump_pressure
    #        self.ERD_cost=
#            self.other_equip_cost=
#            self.equip_cost=self.HPpump_cost + self.BPpump_cost +self.ERD_cost + self.other_equip_cost
            self.CAPEX = self.total_module_cost*CR_factor(self.yrs,self.int_rate) / self.ann_prod
#           (self.cost_sys*1000*self.int_rate*(1+self.int_rate)**self.yrs) / ((1+self.int_rate)**self.yrs-1) / self.Prod
            self.unit_capex=self.total_module_cost/self.capacity  
#        elif self.equip_cost_method=='general':
        else:    
            for i in range(len(self.unit_capex)):
                if i == 0 :
                    if self.unit_capex[i]:
                        self.unit_capex_empirical = [self.unit_capex[i]]
                    else:
                        self.unit_capex_empirical = [3726.1 * self.capacity[i] ** (-0.071)]
                else:
                    if self.unit_capex[i]:
                        self.unit_capex_empirical.append(self.unit_capex[i])
                    else:
                        self.unit_capex_empirical.append(808.39 * self.capacity[i] ** (-0.017))
            
            self.EPC_cost = sum([self.capacity[i] * self.unit_capex_empirical[i] for i in range(len(self.capacity)) ])
            self.drilling_cost = self.well_cost * self.brackish_depth * self.capacity[0] / 3785.4 / 3.26
            
            self.CAPEX =(self.EPC_cost +  self.cost_storage * self.storage_cap + self.drilling_cost)*CR_factor(self.yrs,self.int_rate) / self.ann_prod

        
        # cost of electricity is modified so that it takes the hourly cost from the grid
        # self.cost_elec = self.SEC * (self.fuel_usage * self.coe + (1-self.fuel_usage) * self.sam_coe)
        if self.monthly_input:
            self.monthly_sam_coe = [ self.sam_coe / (1-self.monthly_elec_curtailed[i]) for i in range(len(self.monthly_elec_curtailed))]
            self.monthly_ave_price = [ self.monthly_grid_consumed[i] * self.monthly_grid_coe[i] + (1-self.monthly_grid_consumed[i]) * self.monthly_sam_coe[i] for i in range(len(self.monthly_sam_coe))]
            Month = [31,28,31,30,31,30,31,31,30,31,30,31]
            self.cost_elec = self.SEC * sum([Month[i]/365*self.monthly_ave_price[i] for i in range(len(Month))])

        else:
            
            self.pump_brackish = self.brackish_depth * 0.3048 * 9.8 * 1000 * 2.778 * 1e-7 / self.RO_rr / 0.85# kWh/m3 

            self.cost_elec = (self.SEC + self.pump_brackish) * ( self.fuel_usage * self.ave_grid_coe + (1-self.fuel_usage) * self.sam_coe)
        
        # print(round(self.cost_elec,4), round(self.fuel_usage,4), round(self.ave_grid_coe,4), round(self.sam_coe,4))

        self.insurance_cost = self.EPC_cost * self.insurance / self.ann_prod
        self.OPEX = self.disposal_cost + self.chem_cost +self.labor_cost + self.cost_elec + self.insurance_cost + self.IX_cost  + self.membrane_replacement_cost#maintenance and membrane replacement

        self.LCOW = self.CAPEX + self.OPEX

#        self.test=(self.total_capex*self.int_rate*(1+self.int_rate)**self.yrs) / ((1+self.int_rate)**self.yrs-1) / self.ann_prod
        cost_output = []
        cost_output.append({'Name':'Desal Annualized CAPEX','Value':self.CAPEX,'Unit':'$/m3'})
        cost_output.append({'Name':'Desal OPEX','Value':self.OPEX,'Unit':'$/m3'})
        cost_output.append({'Name':'Levelized cost of water','Value':self.LCOW,'Unit':'$/m3'})
        cost_output.append({'Name':'Levelized cost of electricity (from grid)','Value':self.ave_grid_coe,'Unit':'$/m3'})
        cost_output.append({'Name':'Levelized cost of electricity (from solar field)','Value':self.sam_coe,'Unit':'$/m3'})
        cost_output.append({'Name':'Energy cost','Value':self.cost_elec,'Unit':'$/m3'})    
        

        return cost_output
    
if __name__ == 'main':
    grid_rate = [0.060899999999999996, 0.067, 0.0638, 0.0645, 0.0638, 0.08070000000000001, 0.0802, 0.08070000000000001, 0.0793, 0.0679, 0.0708, 0.0695]
    monthly_grid_consumed = [0.15 for _ in range(12)]
    monthly_elec_curtailed = [0 for _ in range(12)]
    test = RO_cost(monthly_grid_coe = grid_rate,
                   monthly_elec_curtailed= monthly_elec_curtailed,
                   monthly_grid_consumed= monthly_grid_consumed)
    test.lcow()
#    def pump_cost(pumppressure,pumpflowrate):
#        pumpcapex=53*pumppressure*pumpflowrate
#        return pumpcapex
#    def cost_method(self):
# RO_new_cap, prod, grid_perc, SEC, lcoe_fcr, grid_rate = 19930.256162875, \
#                                                         7274543.499450075,\
#                                                         66.97247499821341,\
#                                                         3.776546404259085,\
#                                                         0.06981718655101347, 0.136
# rocost=RO_cost(Capacity = [RO_new_cap], Prod = prod, fuel_usage=grid_perc, sec = SEC,
#                        sam_coe=lcoe_fcr, solar_coe=None, ave_grid_coe = grid_rate).lcow()
# print(rocost)

# test = RO_cost(Capacity = [15000], Prod = 365000*15,unit_capex=[None]).lcow()
# print(test)
