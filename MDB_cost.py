# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 14:04:17 2020

@author: zzrfl
"""

import numpy as np
import math

class MDB_cost(object):
    def __init__(self,
                 crystallizer = None,
                 salt_recovery_revenue = 70, # $/ton NaCl
                 cryst_water_prod = 0, # m3
                 cryst_salt_prod = 0,  # kg
                 Capacity = 1000, # Desalination plant capacity (m3/day)
                 Prod = 328500, # Annual permeate production (m3)
                 downtime = 10,
                 fuel_usage = 0, # %
                 RR = 30, # %
                 Area = 25.92 , # Membrane area (m2)
                 Pflux = 2.832,  # Permeate flux per module (l/h/module)
                 TCO  = 76.288, # Condenser channel outlet temperature (ºC)
                 TEI  = 80, # Condenser channel inlet temperature (C)
                 FFR  = 1100, # Feed flow rate per module (l/h/module)
                 th_module  = 4.560, # Thermal power supplied per module (kW(th)/module)
                 STEC  = 62.180, # Specific thermal energy consumption (kWh(th)/m3) 
                 SEEC = 1.8, # Specific electric energy consumption (kWh(e)/m3) 
#                 GOR = 10.475,  # Gained output ratio
                # downtime = 0.1, # Yearly downtime of the plant (ratio)
                 yrs = 20, # Expected plant lifetime
                 int_rate = 0.04 , # Average interest rate
                 coe = 0.05 , # Unit cost of electricity ($/kWh)
                 coh = 0.01 , # Unit cost of heat ($/kWh)
                 sam_coh = 0.02, # Unit cost of heat from SAM ($/kWh)
                 solar_coh = '',
                 solar_inlet = 85, # Solar field inlet temperature
                 solar_outlet = 95, # Solar field outlet temperature
                 HX_eff = 0.85, # Heat exchanger efficiency
                 cost_module_re = 0.220 , # Cost of module replacement ($/m3)

                 MD_membrane = 0.075*1.2, # Base price of AGMD membrane (k$/m2)
                 HX = 0.35 * 1.2, # Base price of heat exchanger (k$/m2)
                 
                 MD_module = 1.95* 1.2,# Base price of AGMD module assembly (k$/base capacity)
                 MD_module_capacity = 3, # Base capacity of module assembly (modules)   
                                  
                 endplates = 0.85* 1.2,# Base price of HX endplates (k$/base capacity)
                 endplates_capacity =10, # Base capacity of housing rack (m2)
                 heat_cool = 5* 1.2, # Base price of heating/cooling installation (k$/base capacity)
                 heat_cool_capacity =10, # Base capacity of housing rack (m3/h)
                 
                 h_r = 5* 1.2,  # Base price of housing rack (k$/base capacity)
                 h_r_capacity = 3, # Base capacity of housing rack (modules)
                 tank = 5* 1.2,  # Base price of tank (with plumbing) (k$/base capacity)
                 tank_capacity = 3, # Base capacity of tank (modules)
                 pump = 3* 1.2,  # Base price of pump (k$/base capacity)
                 pump_capacity = 5, # Base capacity of pump (m3/h)
                 other = 15* 1.2,  # Base price of controller, cabling and programming (k$/base capacity)
                 other_capacity = 3, # Base capacity of controller, cabling and programming (modules)
                 cost_storage = 26 , # Cost of thermal storage ($/kWh)
                 storage_cap = 13422, # Capacity of thermal storage (kWh)    
                 insurance = 0.5,  # insurance (percentage of CAPEX)
                 maintenance = 1.3, # percentage of CAPEX (%)
                 operational = 0.1 # Specific operational cost ($/m3)    
                 
                 ):
        
        self.operation_hour = 24 #* (1-downtime) # Average daily operation hour (h/day)
        self.Pflux = Pflux
        self.Area = Area
        self.cost_storage = cost_storage
        self.storage_cap = storage_cap
        self.fuel_usage = fuel_usage/100
        self.PF_module = self.Pflux * self.Area
        self.num_modules = math.ceil(Capacity *1000 / self.PF_module / self.operation_hour) # Number of module required
        self.TEI = TEI
        self.TCO = TCO
        self.HX_eff = HX_eff
        self.solar_inlet = solar_inlet
        self.solar_outlet = solar_outlet
        self.th_module = th_module
        self.FFR = FFR
        self.coe = coe
        self.Capacity = Capacity
        self.coh = coh
        self.downtime = downtime / 100
        
        if solar_coh != '':
            self.sam_coh = float(solar_coh)
        else:
            self.sam_coh = sam_coh

        self.cost_module_re = cost_module_re
        self.Prod = Prod
        self.STEC = STEC
        self.yrs = yrs
        self.int_rate = int_rate
        self.SEEC = SEEC
        
        self.MD_membrane = MD_membrane
        self.HX = HX
        self.MD_module = MD_module
        self.MD_module_capacity = MD_module_capacity
        self.endplates = endplates
        self.endplates_capacity =endplates_capacity
        self.heat_cool = heat_cool
        self.heat_cool_capacity =heat_cool_capacity
        
        self.h_r = h_r
        self.h_r_capacity = h_r_capacity
        self.tank = tank
        self.tank_capacity = tank_capacity
        self.pump = pump
        self.pump_capacity = pump_capacity
        self.other = other
        self.other_capacity = other_capacity
        self.RR = RR / 100
        
        self.insurance = insurance / 100
        self.maintenance = maintenance / 100 
        self.operationalcost = operational
        self.crystallier = crystallizer
        self.salt_recovery_revenue = salt_recovery_revenue
        self.cryst_water_prod = cryst_water_prod
        self.cryst_salt_prod = cryst_salt_prod
        
    def lcow(self):

        self.Prod *= (1-self.downtime)
        if self.crystallier:
            self.Prod += self.cryst_water_prod
        self.module_cost = (self.MD_module*self.MD_module_capacity*(self.num_modules/self.MD_module_capacity)**0.8 + self.MD_membrane * self.Area * self.num_modules) 
        self.delta_T2 = self.solar_outlet - self.TEI
        self.delta_T1 = self.solar_inlet - self.TCO
        try:
            self.LMTD = (self.delta_T2 - self.delta_T1) / math.log(self.delta_T2/self.delta_T1)
        except:
            self.LMTD = 10
        self.HX_area = self.num_modules * self.th_module / self.HX_eff / 2.5 / self.LMTD
        self.HX_cost = 2 * ( self.endplates * (self.HX_area/self.endplates_capacity)**0.6 + self.HX * self.HX_area) 
        self.Feed = self.Capacity / 24 / self.RR
        self.other_cap = ((self.h_r * (self.num_modules/self.h_r_capacity)**0.6 
                        + self.tank * (self.num_modules/self.tank_capacity)**0.5
                        + self.pump*(self.Feed/self.pump_capacity)**0.6 
                        + self.other *(self.num_modules/self.other_capacity)**0.3 
                        + 0.25*5 
                        + 2*self.heat_cool*(self.Feed/self.heat_cool_capacity)**0.6) )
        self.cost_sys = (self.module_cost + self.HX_cost + self.other_cap + self.cost_storage * self.storage_cap / 1000)

        self.CAPEX = (self.cost_sys*1000*self.int_rate*(1+self.int_rate)**self.yrs) / ((1+self.int_rate)**self.yrs-1) / (self.Prod+0.1) 

        
        self.cost_elec = self.SEEC * self.coe
        self.maintenancecost = self.cost_sys*1000   * self.maintenance / (self.Prod)
        self.insurancecost = self.cost_sys*1000 * self.insurance / (self.Prod) 
        self.cost_th = self.STEC * (self.fuel_usage * self.coh + (1-self.fuel_usage) * self.sam_coh)
        self.OPEX = self.cost_elec + self.cost_th + self.cost_module_re + self.maintenancecost + self.insurancecost + self.operationalcost

        
        self.LCOW = self.CAPEX + self.OPEX
        
        cost_output = []
        cost_output.append({'Name':'Desal CAPEX','Value':self.CAPEX,'Unit':'$/m3'})
        cost_output.append({'Name':'Desal OPEX','Value':self.OPEX,'Unit':'$/m3'})
        cost_output.append({'Name':'Levelized cost of water','Value':self.LCOW,'Unit':'$/m3'})
        cost_output.append({'Name':'Levelized cost of heat (from fossile fuel)','Value':self.coh,'Unit':'$/kWh'})
        cost_output.append({'Name':'Levelized cost of heat (from solar field)','Value':self.sam_coh,'Unit':'$/kWh'})
        cost_output.append({'Name':'Cost of heat per unit','Value':self.cost_th,'Unit':'$/m3'})   
        cost_output.append({'Name':'Energy cost','Value':self.cost_th + self.cost_elec,'Unit':'$/m3'})   
        if self.crystallier:
            self.salt_recovery_revenue = self.cryst_salt_prod * self.salt_recovery_revenue / 100 / self.Prod
            cost_output.append({'Name':'Salt recovery revenue','Value':self.salt_recovery_revenue,'Unit':'$/m3 fresh water'})        
        return cost_output
