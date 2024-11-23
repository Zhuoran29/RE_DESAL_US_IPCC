# -*- coding: utf-8 -*-
"""
Model integrated from WaterTAP
Created on Dec 01 2022

@author: Zhuoran Zhang
"""

from pyomo.environ import (
    ConcreteModel,
    SolverFactory,
    TerminationCondition,
    assert_optimal_termination,
    value,
    Constraint,
    Var,
    Objective,
    Expression,
)
from pyomo.environ import units as pyunits
from pyomo.util.check_units import (
    assert_units_consistent,
    assert_units_equivalent,
    check_units_equivalent,
)
import pyomo.util.infeasible as infeas
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_activated_constraints,
    number_unfixed_variables_in_activated_equalities,
    number_activated_equalities,
    number_unused_variables,
)

import idaes.core.util.model_statistics as stats
from idaes.core.util.constants import Constants
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core.solvers import get_solver
from idaes.core import UnitModelCostingBlock

from watertap.property_models import cryst_prop_pack as props
from watertap.unit_models.crystallizer import Crystallization
from watertap.costing import WaterTAPCosting, CrystallizerCostType

from io import StringIO
from pyomo.util.infeasible import (
    log_active_constraints,
    log_close_to_bounds,
    log_infeasible_bounds,
    log_infeasible_constraints,
)
from pyomo.common.log import LoggingIntercept
import logging

class Crystallizer(object):
    
    def __init__(self,
         specified_capacity       =  500, # Userspecified capacity (m3/day)
         desal_brine_capacity     = 500,  # Use desal brine flow rate (m3/day)
         specified_salinity  =    70     ,        # Feed salinity (g/L)
         desal_brines_salinity = 70,
         feed_pressure        =  101325  , # Feed pressure (Pa)
         feed_temperature     =  20,    # Feed temperature (oC)
         crystallizer_temperature  = 55 , # Crystallizer temperature (oC)
         crystallizer_yield        = 0.4 , # Crystallizer yield
         salt_revenue = 70, # Salt recovery value ($/ton)
         Fossil_f = 1 # Fossil fuel fraction
         ):
        if specified_capacity:
            self.feed_flow_mass = specified_capacity / 3.6 / 24 # Convert to kg/s
        else:
            self.feed_flow_mass = desal_brine_capacity / 3.6 / 24
        
        if specified_salinity:
            self.feed_mass_frac_NaCl = specified_salinity / 1000 # Convert to mass fraction
        else:
            self.feed_mass_frac_NaCl = desal_brines_salinity / 1000

        self.feed_pressure = feed_pressure
        self.feed_temperature = feed_temperature
        self.crystallizer_temperature = crystallizer_temperature
        self.crystallizer_yield = crystallizer_yield
        self.salt_revenue = salt_revenue
        self.Fossil_f  = Fossil_f

    def design(self):

        m = self.initialization()

        # fully specify system
        feed_flow_mass = self.feed_flow_mass
        feed_mass_frac_NaCl = self.feed_mass_frac_NaCl
        feed_mass_frac_H2O = 1 - self.feed_mass_frac_NaCl
        feed_pressure = self.feed_pressure
        feed_temperature = 273.15 + self.feed_temperature
        eps = 1e-6
        crystallizer_temperature = 273.15 + self.crystallizer_temperature
        crystallizer_yield = self.crystallizer_yield

        # Fully define feed
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            feed_flow_mass * feed_mass_frac_NaCl
        )
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Sol", "NaCl"].fix(eps)
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(eps)
        m.fs.unit.inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.inlet.temperature[0].fix(feed_temperature)

        # Update scaling factors
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e-0*feed_flow_mass, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e-1*feed_flow_mass, index=("Liq", "NaCl")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e-0*feed_flow_mass, index=("Vap", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e-1*feed_flow_mass, index=("Sol", "NaCl")
        )

        # Solve the model
        solver = get_solver()
        results = solver.solve(m, tee=True)
        assert_optimal_termination(results)

        b = m.fs.unit
        solid_mass = b.solids.flow_mass_phase_comp[0, "Sol", "NaCl"].value *  m.fs.costing.utilization_factor.value
        flow_rate = 365 * 86400 * sum(b.properties_in[0].flow_vol_phase[p].value for p in b.properties_in[0].params.phase_list)
        solid_revenue = solid_mass * 365 * 86400 * (self.salt_revenue/100)     
        vapor_mass_flow_rate = b.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"].value #kg/s
        self.vapor_vol_flow_rate = vapor_mass_flow_rate / 1000 * 3600   # m3/hr

        self.P_req = value(b.work_mechanical[0])
        self.salt_mass_flow = value(b.crystallization_yield["NaCl"]* b.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"]) * 3600 # kg/hr
        self.design_output = []    
        self.design_output.append({'Name':'Crystallizer thermal energy requirement','Value':value(b.work_mechanical[0]),'Unit':'kW'})
        self.design_output.append({'Name':'Specific thermal energy consumption','Value':value(b.work_mechanical[0])/3.6/self.feed_flow_mass,'Unit':'kWh/m3'})
        self.design_output.append({'Name':'Solid mass in solids stream','Value':solid_mass,'Unit':'kg/s'})
        self.design_output.append({'Name':'Solid mass in liquid stream','Value':value(b.outlet.flow_mass_phase_comp[0, "Liq", "NaCl"]),'Unit':'kg/s'})
        self.design_output.append({'Name':'Liquid stream solvent flow','Value':value(b.outlet.flow_mass_phase_comp[0, "Liq", "H2O"]),'Unit':'kg/s'})
        self.design_output.append({'Name':'Liquid stream composition','Value':value(b.properties_out[0].mass_frac_phase_comp["Liq", "NaCl"]),'Unit':''})
        self.design_output.append({'Name':'Fresh water production rate','Value':vapor_mass_flow_rate,'Unit':'kg/s'})
        self.design_output.append({'Name':'Saturation pressure','Value':value(b.pressure_operating),'Unit':'Pa'})
        self.design_output.append({'Name':'Crystallizer diameter','Value':value(b.diameter_crystallizer),'Unit':'m'})
        self.design_output.append({'Name':'Minimum active volume','Value':value(b.volume_suspension),'Unit':'m3'})
        self.design_output.append({'Name':'Residence time','Value':value(b.t_res),'Unit':'hr'})
        self.design_output.append({'Name':'Total capital cost','Value':value(m.fs.costing.aggregate_capital_cost), 'Unit':'$'})
        # self.design_output.append({'Name':'Salt recovery revenue','Value':solid_revenue, 'Unit':'$'})


        # self.cost_output = []    
        # flow_rate = 365 * 86400 * sum(m.fs.unit.properties_in[0].flow_vol_phase[p].value for p in m.fs.unit.properties_in[0].params.phase_list) * m.fs.costing.utilization_factor.value
        # OPEX = m.fs.costing.total_operating_cost.value / flow_rate / m.fs.costing.utilization_factor.value
        # CAPEX = m.fs.costing.aggregate_capital_cost.value * m.fs.costing.factor_capital_annualization.value / flow_rate / m.fs.costing.utilization_factor.value
        # self.cost_output.append({'Name':'LCOW','Value':m.fs.costing.LCOW.value,'Unit':'$/m3 brine'})
        # self.cost_output.append({'Name':'OPEX','Value':OPEX,'Unit':'$/m3 brine'})
        # self.cost_output.append({'Name':'CAPEX','Value':m.fs.costing.aggregate_capital_cost.value,'Unit':'$/m3 brine'})

        return self.design_output

    def initialization(self):
        solver = get_solver()

        # setup model
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = props.NaClParameterBlock()
        m.fs.costing = WaterTAPCosting()

        # create uit models
        m.fs.unit = Crystallization(property_package=m.fs.properties)
        m.fs.unit.crystal_growth_rate.fix()
        m.fs.unit.souders_brown_constant.fix()
        m.fs.unit.crystal_median_length.fix()
        
        # fully specify system
        feed_flow_mass = 1
        feed_mass_frac_NaCl = 0.185
        feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
        feed_pressure = 101325
        feed_temperature = 273.15 + 20
        eps = 1e-6
        crystallizer_temperature = 273.15 + 55
        crystallizer_yield = 0.40

        # Fully define feed
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            feed_flow_mass * feed_mass_frac_NaCl
        )
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Sol", "NaCl"].fix(eps)
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(eps)
        m.fs.unit.inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.inlet.temperature[0].fix(feed_temperature)

        # Define operating conditions
        m.fs.unit.temperature_operating.fix(crystallizer_temperature)
        m.fs.unit.crystallization_yield["NaCl"].fix(crystallizer_yield)

        # scaling
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e-0, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e-1, index=("Liq", "NaCl")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e-0, index=("Vap", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e-1, index=("Sol", "NaCl")
        )

        iscale.calculate_scaling_factors(m.fs)

        # initialize units
        m.fs.unit.initialize()

        # solve
        results = solver.solve(m, tee=False)

        # costing
        m.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method_arguments={"cost_type": CrystallizerCostType.mass_basis},
        )
        m.fs.costing.cost_process()
        m.fs.costing.add_annual_water_production(m.fs.unit.properties_in[0].flow_vol)
        m.fs.costing.add_LCOW(m.fs.unit.properties_in[0].flow_vol)
        m.fs.costing.add_specific_energy_consumption(m.fs.unit.properties_in[0].flow_vol)

        # consistent units
        assert_units_consistent(m)

        # re-solve model
        results = solver.solve(m, tee=False)

        return m
        
    def simulation(self, gen, storage, desal_load):
        self.thermal_load = self.P_req + desal_load # kW
        self.max_prod = self.vapor_vol_flow_rate  # m3/h
        self.storage_cap = storage * self.thermal_load # kWh
        
        to_desal = [0 for i in range(len(gen))]
        to_storage =  [0 for i in range(len(gen))]
        storage_load =  [0 for i in range(len(gen))]
        storage_cap_1 =  [0 for i in range(len(gen))]
        storage_cap_2 = [0 for i in range(len(gen))]
        storage_status =  [0 for i in range(len(gen))]
        solar_loss =  [0 for i in range(len(gen))]
        load =  [0 for i in range(len(gen))]
        prod =  [0 for i in range(len(gen))]
        salt =  [0 for i in range(len(gen))]
        fuel =  [0 for i in range(len(gen))]   
        energy_consumption =  [0 for i in range(len(gen))] 
        for i in range(len(gen)):
            to_desal[i] = min(self.thermal_load, gen[i])
            to_storage[i] = abs(gen[i] - to_desal[i])
            storage_load[i] = gen[i] - self.thermal_load
            if i != 0:
                storage_cap_1[i] = storage_status[i-1]
            storage_cap_2[i] = max(storage_load[i] + storage_cap_1[i], 0)
            storage_status[i] = min(storage_cap_2[i] , self.storage_cap)
            solar_loss[i] = abs(storage_status[i] - storage_cap_2[i])
            load[i] = to_desal[i] + max(0, storage_cap_1[i] - storage_cap_2[i])
            if load[i] / self.thermal_load < self.Fossil_f:
                fuel[i] = self.thermal_load - load[i]
            energy_consumption[i] = fuel[i]+load[i]     
            prod[i] = (fuel[i]+load[i] )/ self.thermal_load * self.max_prod
            salt[i] = (fuel[i]+load[i] )/ self.thermal_load * self.salt_mass_flow
        Month = [0,31,59,90,120,151,181,212,243,273,304,334,365]
        Monthly_prod = [ sum( prod[Month[i]*24:(Month[i+1]*24)] ) for i in range(12) ]
               
        simu_output = []
        simu_output.append({'Name':'Water production','Value':prod,'Unit':'m3'})
        simu_output.append({'Name':'Salt production','Value':salt,'Unit':'kg'})
        simu_output.append({'Name':'Storage status','Value':storage_status,'Unit':'kWh'})
        simu_output.append({'Name':'Storage Capacity','Value':self.storage_cap,'Unit':'kWh'})
        simu_output.append({'Name':'Fossil fuel usage','Value':fuel,'Unit':'kWh'})
        simu_output.append({'Name':'Total water production','Value':sum(prod),'Unit':'m3'})
        simu_output.append({'Name':'Total salt production','Value':sum(salt),'Unit':'kg'})
        simu_output.append({'Name':'Monthly water production','Value': Monthly_prod,'Unit':'m3'})
        simu_output.append({'Name':'Total fossil fuel usage','Value':sum(fuel),'Unit':'kWh'})
        simu_output.append({'Name':'Percentage of fossil fuel consumption','Value':sum(fuel)/max(1,sum(energy_consumption))*100,'Unit':'%'})        
        simu_output.append({'Name':'Solar energy curtailment','Value':solar_loss,'Unit':'kWh'})
        simu_output.append({'Name':'Curtailed solar thermal energy','Value':(sum(gen) - sum(load)) / 1000000 ,'Unit':'GWh'})   
        simu_output.append({'Name':'Percentage of curtailed energy','Value':(sum(gen) - sum(load)) / sum(gen) * 100 ,'Unit':'%'})
               
        return simu_output

#%%
if __name__ == 'main':
    case = Crystallizer(         
         specified_capacity       =  5000, #86.4 , # m3/day
         specified_salinity  =  70, #212.6     , # Feed salinity (g/L)
         feed_pressure        =  101325  , # Feed pressure (Pa)
         feed_temperature     =  20,    # Feed temperature (oC)
         crystallizer_temperature  = 55 , # Crystallizer temperature (oC)
         crystallizer_yield        = 0.4 , # Crystallizer yield
         Fossil_f = 1 # Fossil fuel fraction
         )
    output, m = case.design()
    for i in output:
        print(i)

# %%