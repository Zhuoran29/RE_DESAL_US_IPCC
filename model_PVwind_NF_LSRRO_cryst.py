from VAGMD_batch.VAGMD_batch import VAGMD_batch
from MDB_cost import MDB_cost
from StaticCollector_flatplate import StaticCollector_fp
from SC_ETC_cost import ETC_cost
from PVWatts import PVwatts
from wind import wind

import matplotlib.pyplot as plt
import numpy as np
from pyomo.environ import value, units as pyunits
from idaes.core.util.constants import Constants
import sys
import os
import time
import csv
import pandas as pd
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    badly_scaled_var_generator,
)
import idaes.core.util.scaling as iscale
from urllib.request import urlopen
# from lxml import etree
from bs4 import BeautifulSoup
import requests
import json
from Wind_resource.wind_toolkit_api import get_wind_resource
from lookup_openei_rates import lookup_rates
from pyomo.environ import (
    ConcreteModel,
    TerminationCondition,
    SolverStatus,
    Objective,
    Expression,
    maximize,
    value,
    Set,
    Var,
    log,
    units as pyunits,
)
from watertap_reflo.src.watertap_contrib.reflo.analysis.example_flowsheets.multi_effect_NaCl_crystallizer import (
    build_fs_multi_effect_crystallizer,
    add_costings,
    multi_effect_crystallizer_initialization,
    get_model_performance,
)
import watertap_reflo.src.watertap_contrib.reflo.property_models.cryst_prop_pack as props
from watertap.core.solvers import get_solver
from watertap_reflo.src.watertap_contrib.reflo.costing import TreatmentCosting, CrystallizerCostType

from watertap.unit_models.zero_order import NanofiltrationZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.costing.zero_order_costing import ZeroOrderCosting


from fs_NF_LSRRO_crys import fs_NF_LSRRO_crys
from fs_NF_LSRRO_crys_simplified import fs_NF_LSRRO_crys_simple

from Wind_resource import wind25, wind50, wind75, wind100, wind125, wind150, wind175, wind200, wind225, wind250

solver = get_solver()

import warnings
warnings.filterwarnings("ignore")


def PV_wind_cal():
    TX_counties = 'TX_SEEC.csv'
    US_solar_directory = '/Users/zhuoranzhang/Documents/DOE_Chile/DOE_CSP_PROJECT/SAM_flatJSON/solar_resource/'
    TX_wind_directory = '/Users/zhuoranzhang/Documents/Github_repos/Modeling_code/Wind_resource/'
    test_csv = "test2.csv"

    df = pd.read_csv(test_csv, encoding = 'unicode_escape')
    csv_outfile1  = 'TX_PV_wind_100MW.csv' 
    csv_outfile2  = 'TX_PV_wind_150MW.csv' 
    csv_outfile3  = 'TX_PV_wind_200MW.csv' 
    csv_outfile4  = 'TX_PV_wind_250MW.csv' 
    csv_outfile5  = 'TX_PV_wind_300MW.csv' 
    data_file1 = open(csv_outfile1, 'w', newline='') 
    data_file2 = open(csv_outfile2, 'w', newline='') 
    data_file3 = open(csv_outfile3, 'w', newline='') 
    data_file4 = open(csv_outfile4, 'w', newline='') 
    data_file5 = open(csv_outfile5, 'w', newline='') 
    csv_writer1 = csv.writer(data_file1)
    csv_writer2 = csv.writer(data_file2)
    csv_writer3 = csv.writer(data_file3)
    csv_writer4 = csv.writer(data_file4)
    csv_writer5 = csv.writer(data_file5)
    csv_writers = [
                   csv_writer1, 
                   csv_writer2, 
                   csv_writer3, 
                   csv_writer4, 
                   csv_writer5, 
                #    csv_writer6
                   ]
    for i in csv_writers:
        i.writerow([
                    "County", "State_ID", "County_ID", "Latitude", "Longitude",
                    "SEEC (kWh/m3)", "Desal capacity (m3/day)",
                    "LCOE solar ($/kWh)", "LCOE solar adjusted ($/kWh)",
                    "LCOE wind ($/kWh)", "LCOE wind adjusted ($/kWh)",
                    "LCOW ($/m3)", "Solar capacity (MW)", "Wind capacity (MW)",
                    "Solar_perc (%)", "Solar_curt (%)", "Wind_perc (%)", "Wind_curt (%)",
                    "Renewable_curt (%)", "Capacity factor", "Avg water prod (m3/day)",
                    ])
        
    energy_demand = 100 * 1000 # kW
    failed_county = []

    # Solar capacities in different combinations
    MW100 = [25, 50, 75]
    MW150 = [25, 50, 75, 100, 125]
    MW200 = [50, 75, 100, 125, 150]
    MW250 = [50, 75, 100, 125, 150, 175, 200]
    MW300 = [75, 100, 125, 150, 175, 200, 225]
    combinations = [MW100, MW150, MW200, MW250, MW300]
    demands = [100, 150, 200, 250, 300]

    wind_functions = {25: wind25,
                      50: wind50,
                      75: wind75,
                      100: wind100,
                      125: wind125,
                      150: wind150,
                      175: wind175,
                      200: wind200,
                      225: wind225,
    }

    start_time = time.time()
    progress = 0
    for index, row in df.iterrows():
        # try:     
        county = row["County_ID"]
        wind_file = f"wind_{county}.csv"
        weatherfile =  row["filename"]
        SEEC = row["NF SEEC (kWh-e/m3)"] + row["Pumping energy (kWh-e/m3)"] + row["LSRRO SEEC (kWh-e/m3)"] # kWh/m3
        water_prod = energy_demand / SEEC  # m3/hour
        capex_lcow = 0.45 # $/m3
        fixed_opex_lcow = 0.25 # $/m3
        STEC = row["Cryst STEC (kWh-th/m3)"]
        cost_salt = row["Salt recovery ($/m3)"]
        salt_prod = cost_salt / (-0.024) * water_prod # kg/hr
        # print(water_prod)
        
        weatherfile_path = US_solar_directory + weatherfile 
        wind_resource_path = TX_wind_directory + wind_file


        for i in range(len(combinations)):
            demand = demands[i]
            results = {
                "LCOW ($/m3)": float("inf"),
                "Solar capacity (MW)": 0,
                "Wind capacity (MW)": 0,
                "Solar_perc (%)": 0,
                "Solar_curt (%)": 0,
                "Wind_perc (%)": 0,
                "Wind_curt (%)": 0,
                "Renewable_curt (%)":0,
                "Capacity factor": 0,
                "Avg water prod (m3/day)":0,

            }
            for solar in combinations[i]:

                wind = demand - solar
                solar_gen, capacity_factor, kwh_per_kw, lcoe_solar, monthly_energy = PVwatts(
                                                                weatherfile_path,
                                                                solar*1000,
                                                                0)
                wind_gen, lcoe_wind = wind_functions[wind](wind_resource=wind_resource_path)

                if len(solar_gen) != len(wind_gen):
                    if len(solar_gen) < len(wind_gen):
                        wind_gen = wind_gen[0: len(solar_gen)-1]
                    else:
                        solar_gen = solar_gen[0: len(wind_gen)-1]
                
                flag = lcoe_solar < lcoe_wind

                solar_load, wind_load, to_desal, solar_curt, wind_curt, curt, water_gen, salt_gen = simulation(
                                        solar_gen,
                                        wind_gen,
                                        power = 100000, # kW
                                        water_prod = water_prod, # m3/hr,
                                        salt_prod = salt_prod, # kg/hr
                                        flag = flag,
                                        storage = 0,
                                        )
                solar_perc = sum(solar_load) / 100000 / 8760
                wind_perc = sum(wind_load) / 100000 / 8760

                solar_curt_perc = sum(solar_curt)/sum(solar_gen)
                wind_curt_perc = sum(wind_curt)/sum(wind_gen)
                curt_perc = sum(curt) / (sum(solar_gen) + sum(wind_gen))
                cf = sum(water_gen)/water_prod/8760

                lcoe_solar_adjusted = round(lcoe_solar / (1 - solar_curt_perc),3)
                lcoe_wind_adjusted = round(lcoe_wind / (1 - wind_curt_perc),3)

               # Calculate new LCOW
                cost_th = STEC * 0.02 # $/m3
                cost_elec = (sum(solar_load) * lcoe_solar_adjusted + sum(wind_load) * lcoe_wind_adjusted) / sum(water_gen)
                cost_salt = sum(salt_gen) * (-0.024) / sum(water_gen)

                designed_annual_water_prod = water_prod * 8760 * 0.9
                annualized_CAPEX = capex_lcow * designed_annual_water_prod

                LCOW = (annualized_CAPEX/sum(water_gen)
                        + fixed_opex_lcow
                        + cost_th
                        + cost_elec
                        + cost_salt
                        )                
                # print(f"{solar} solar + {wind} wind: {round(solar_perc*100,1)}%, {round(wind_perc*100,1)}%, {round(100*curt_perc,1)}%, {round(100*cf,1)}%, {round(lcoe_solar_adjusted,3)}, {round(lcoe_wind_adjusted,3)}, {round(LCOW,3)};")
   
                if LCOW < results["LCOW ($/m3)"]:
                    results = {
                    "LCOW ($/m3)": LCOW,
                    "Solar capacity (MW)": solar,
                    "Wind capacity (MW)": wind,
                    "Solar_perc (%)": solar_perc,
                    "Solar_curt (%)": solar_curt_perc,
                    "Wind_perc (%)": wind_perc,
                    "Wind_curt (%)": wind_curt_perc,
                    "Renewable_curt (%)": curt_perc,
                    "Capacity factor": cf,
                    "Avg water prod (m3/day)": sum(water_gen) / 365,
                }

                    "County", "State_ID", "County_ID", "Latitude", "Longitude",
                    "SEEC (kWh/m3)", "Desal capacity (m3/day)",
                    "LCOE solar ($/kWh)", "LCOE solar adjusted ($/kWh)",
                    "LCOE wind ($/kWh)", "LCOE wind adjusted ($/kWh)",
                    "LCOW ($/m3)", "Solar capacity (MW)", "Wind capacity (MW)",
                    "Solar_perc (%)", "Solar_curt (%)", "Wind_perc (%)", "Wind_curt (%)",
                    "Renewable_curt (%)", "Capacity factor", "Avg water prod (m3/day)",
            
            csv_writers[i].writerow([
                row["County"], row["State_ID"], row["County_ID"], row["Latitude"], row["Longitude"],
                SEEC, water_prod*24,
                lcoe_solar, lcoe_solar_adjusted,
                lcoe_wind, lcoe_wind_adjusted, 
                results["LCOW ($/m3)"], results["Solar capacity (MW)"], results["Wind capacity (MW)"],
                results["Solar_perc (%)"], results["Solar_curt (%)"], results["Wind_perc (%)"], results["Wind_curt (%)"],
                results["Renewable_curt (%)"], results["Capacity factor"], results["Avg water prod (m3/day)"],
            ])

        progress += 1
        print('Progress: ', round(progress/len(df.index)*100,1), '%', 'Processing time: ', round((time.time()-start_time)/60,1), 'min')
            
        # except:
        #     failed_county.append(county)

def PV_wind_battery_cal():
    TX_counties = 'TX_SEEC.csv'
    US_solar_directory = '/Users/zhuoranzhang/Documents/DOE_Chile/DOE_CSP_PROJECT/SAM_flatJSON/solar_resource/'
    TX_wind_directory = '/Users/zhuoranzhang/Documents/Github_repos/Modeling_code/Wind_resource/'
    test_csv = "test2.csv"
    df = pd.read_csv(test_csv, encoding = 'unicode_escape')
    
    # csv_outfile1  = 'TX_PV_wind_100MW_batt.csv' 
    # data_file1 = open(csv_outfile1, 'w', newline='') 
    # csv_writer1 = csv.writer(data_file1)

    # csv_outfile2  = 'TX_PV_wind_150MW_batt.csv' 
    # data_file2 = open(csv_outfile2, 'w', newline='') 
    # csv_writer2 = csv.writer(data_file2)

    csv_outfile3  = 'TX_PV_wind_200MW_batt.csv' 
    data_file3 = open(csv_outfile3, 'w', newline='') 
    csv_writer3 = csv.writer(data_file3)

    csv_outfile4  = 'TX_PV_wind_250MW_batt.csv' 
    data_file4 = open(csv_outfile4, 'w', newline='')  
    csv_writer4 = csv.writer(data_file4)

    csv_outfile5  = 'TX_PV_wind_300MW_batt.csv' 
    data_file5 = open(csv_outfile5, 'w', newline='') 
    csv_writer5 = csv.writer(data_file5)

    csv_outfile6  = 'TX_PV_wind_350MW_batt.csv' 
    data_file6 = open(csv_outfile6, 'w', newline='')
    csv_writer6 = csv.writer(data_file6)

    csv_writers = [
                #    csv_writer1, 
                #    csv_writer2, 
                   csv_writer3, 
                   csv_writer4, 
                   csv_writer5, 
                   csv_writer6,
                   ]
    for i in csv_writers:
        i.writerow([
                    "County", "State_ID", "County_ID", "Latitude", "Longitude",
                    "SEEC (kWh/m3)", "Desal capacity (m3/day)",
                    "LCOE solar ($/kWh)", "LCOE solar adjusted ($/kWh)",
                    "LCOE wind ($/kWh)", "LCOE wind adjusted ($/kWh)",
                    "LCOW ($/m3)", "Solar capacity (MW)", "Wind capacity (MW)", "Battery (hr)",
                    "Solar_perc (%)", "Solar_curt (%)", "Wind_perc (%)", "Wind_curt (%)",
                    "Renewable_perc (%)", "Renewable_curt (%)", "Capacity factor", "Avg water prod (m3/day)",
                    ])
        
    energy_demand = 100 * 1000 # kW
    failed_county = []

    # Solar capacities in different combinations
    MW100 = [25, 50, 75]
    MW150 = [25, 50, 75, 100, 125]
    MW200 = [50, 75, 100, 125, 150]
    MW250 = [50, 75, 100, 125, 150, 175, 200]
    MW300 = [75, 100, 125, 150, 175, 200, 225]
    MW350 = [100, 125, 150, 175, 200, 225, 250]
    combinations = [
        # MW100, 
        # MW150, 
        MW200,
        MW250, 
        MW300, 
        MW350]
    demands = [
        # 100,
        # 150, 
        200, 
        250, 
        300,
        350]
    battery_power_cost = {0: 0,
                          2: 710,
                          4: 1036,
                          6: 1362,
                          8: 1688,
                          10: 2014,
                          }
    wind_functions = {25: wind25,
                      50: wind50,
                      75: wind75,
                      100: wind100,
                      125: wind125,
                      150: wind150,
                      175: wind175,
                      200: wind200,
                      225: wind225,
                      250: wind250,
    }

    start_time = time.time()
    progress = 0
    for index, row in df.iterrows():
        # try:     
        county = row["County_ID"]
        wind_file = f"wind_{county}.csv"
        weatherfile =  row["filename"]
        SEEC = row["NF SEEC (kWh-e/m3)"] + row["Pumping energy (kWh-e/m3)"] + row["LSRRO SEEC (kWh-e/m3)"] # kWh/m3
        water_prod = energy_demand / SEEC  # m3/hour
        capex_lcow = 0.45 # $/m3
        fixed_opex_lcow = 0.25 # $/m3
        STEC = row["Cryst STEC (kWh-th/m3)"]
        cost_salt = row["Salt recovery ($/m3)"]
        salt_prod = cost_salt / (-0.024) * water_prod # kg/hr
        # print(water_prod)
        
        weatherfile_path = US_solar_directory + weatherfile 
        wind_resource_path = TX_wind_directory + wind_file

        for i in range(len(combinations)):
            demand = demands[i]
            results = {
                "LCOW ($/m3)": float("inf"),
                "Solar capacity (MW)": 0,
                "Wind capacity (MW)": 0,
                "Battery (hr)": 0,
                "Solar_perc (%)": 0,
                "Solar_curt (%)": 0,
                "Wind_perc (%)": 0,
                "Wind_curt (%)": 0,
                "Renewable_perc (%)": 0,
                "Renewable_curt (%)":0,
                "Capacity factor": 0,
                "Avg water prod (m3/day)":0,
            }
            for solar in combinations[i]:
                wind = demand - solar
                solar_gen, capacity_factor, kwh_per_kw, lcoe_solar, monthly_energy = PVwatts(
                                                                weatherfile_path,
                                                                solar*1000,
                                                                0)
                wind_gen, lcoe_wind = wind_functions[wind](wind_resource=wind_resource_path)

                if len(solar_gen) != len(wind_gen):
                    if len(solar_gen) < len(wind_gen):
                        wind_gen = wind_gen[0: len(solar_gen)-1]
                    else:
                        solar_gen = solar_gen[0: len(wind_gen)-1]
                
                flag = lcoe_solar < lcoe_wind

                for storage in [0,2,4,6,8,10]:
                    solar_load, wind_load, to_desal, solar_curt, wind_curt, curt, water_gen, salt_gen = simulation_battery(
                                            solar_gen,
                                            wind_gen,
                                            power = 100000, # kW
                                            water_prod = water_prod, # m3/hr,
                                            salt_prod = salt_prod, # kg/hr
                                            flag = flag,
                                            storage = storage,
                                            )
                    solar_perc = sum(solar_load) / 100000 / 8760
                    wind_perc = sum(wind_load) / 100000 / 8760

                    solar_curt_perc = sum(solar_curt)/sum(solar_gen)
                    wind_curt_perc = sum(wind_curt)/sum(wind_gen)
                    curt_perc = sum(curt) / (sum(solar_gen) + sum(wind_gen))
                    cf = sum(water_gen)/water_prod/8760

                    lcoe_solar_adjusted = round(lcoe_solar / (1 - solar_curt_perc),3)
                    lcoe_wind_adjusted = round(lcoe_wind / (1 - wind_curt_perc),3)

                # Calculate new LCOW
                    cost_th = STEC * 0.02 # $/m3
                    cost_elec = (sum(solar_load) * lcoe_solar_adjusted + sum(wind_load) * lcoe_wind_adjusted) / sum(water_gen)
                    cost_salt = sum(salt_gen) * (-0.024) / sum(water_gen)

                    battery_cost = (157 * storage + battery_power_cost[storage]) * 100000

                    designed_annual_water_prod = water_prod * 8760 * 0.9
                    annualized_CAPEX = capex_lcow * designed_annual_water_prod + battery_cost * 0.1

                    LCOW = (annualized_CAPEX/sum(water_gen)
                            + fixed_opex_lcow
                            + cost_th
                            + cost_elec
                            + cost_salt
                            )                
                    print(f"{solar} solar + {wind} wind + {storage} batt: {round(solar_curt_perc*100,1)}%, {round(wind_curt_perc*100,1)}%, {round(100*curt_perc,1)}%, {round(100*cf,1)}%, {round(lcoe_solar_adjusted,3)}, {round(lcoe_wind_adjusted,3)}, {round(LCOW,3)};")
   
                    if LCOW < results["LCOW ($/m3)"]:
                        results = {
                        "LCOW ($/m3)": LCOW,
                        "Solar capacity (MW)": solar,
                        "Wind capacity (MW)": wind,
                        "Battery (hr)": storage,
                        "Solar_perc (%)": solar_perc,
                        "Solar_curt (%)": solar_curt_perc,
                        "Wind_perc (%)": wind_perc,
                        "Wind_curt (%)": wind_curt_perc,
                        "Renewable_perc (%)": solar_perc + wind_perc,
                        "Renewable_curt (%)": curt_perc,
                        "Capacity factor": cf,
                        "Avg water prod (m3/day)": sum(water_gen) / 365,
                    }
            
            csv_writers[i].writerow([
                row["County"], row["State_ID"], row["County_ID"], row["Latitude"], row["Longitude"],
                SEEC, water_prod*24,
                lcoe_solar, lcoe_solar_adjusted,
                lcoe_wind, lcoe_wind_adjusted, 
                results["LCOW ($/m3)"], results["Solar capacity (MW)"], results["Wind capacity (MW)"], results["Battery (hr)"],
                results["Solar_perc (%)"], results["Solar_curt (%)"], results["Wind_perc (%)"], results["Wind_curt (%)"],
                results["Renewable_perc (%)"], results["Renewable_curt (%)"], results["Capacity factor"], results["Avg water prod (m3/day)"],
            ])

        progress += 1
        print('Progress: ', round(progress/len(df.index)*100,1), '%', 'Processing time: ', round((time.time()-start_time)/60,1), 'min')
            
def PV_wind_battery_grid_cal():
    TX_counties = 'TX_SEEC.csv'
    US_solar_directory = '/Users/zhuoranzhang/Documents/DOE_Chile/DOE_CSP_PROJECT/SAM_flatJSON/solar_resource/'
    TX_wind_directory = '/Users/zhuoranzhang/Documents/Github_repos/Modeling_code/Wind_resource/'
    test_csv = "test2.csv"

    df = pd.read_csv(TX_counties, encoding = 'unicode_escape')
    
    # df = pd.read_csv(test_csv, encoding = 'unicode_escape')
    # csv_outfile1  = 'TX_PV_wind_100MW_batt_grid.csv' 
    # data_file1 = open(csv_outfile1, 'w', newline='') 
    # csv_writer1 = csv.writer(data_file1)

    # csv_outfile2  = 'TX_PV_wind_150MW_batt_grid.csv' 
    # data_file2 = open(csv_outfile2, 'w', newline='') 
    # csv_writer2 = csv.writer(data_file2)

    csv_outfile3  = 'TX_PV_wind_200MW_batt_grid.csv' 
    data_file3 = open(csv_outfile3, 'w', newline='') 
    csv_writer3 = csv.writer(data_file3)

    csv_outfile4  = 'TX_PV_wind_250MW_batt_grid.csv' 
    data_file4 = open(csv_outfile4, 'w', newline='')  
    csv_writer4 = csv.writer(data_file4)

    csv_outfile5  = 'TX_PV_wind_300MW_batt_grid.csv' 
    data_file5 = open(csv_outfile5, 'w', newline='') 
    csv_writer5 = csv.writer(data_file5)

    csv_outfile6  = 'TX_PV_wind_350MW_batt_grid.csv' 
    data_file6 = open(csv_outfile6, 'w', newline='')
    csv_writer6 = csv.writer(data_file6)

    csv_writers = [
                #    csv_writer1, 
                #    csv_writer2, 
                   csv_writer3, 
                   csv_writer4, 
                   csv_writer5, 
                   csv_writer6,
                   ]
    for i in csv_writers:
        i.writerow([
                    "County", "State_ID", "County_ID", "Latitude", "Longitude",
                    "SEEC (kWh/m3)", "Desal capacity (m3/day)",
                    "LCOE solar ($/kWh)", "LCOE solar adjusted ($/kWh)",
                    "LCOE wind ($/kWh)", "LCOE wind adjusted ($/kWh)",
                    "LCOW ($/m3)", "Solar capacity (MW)", "Wind capacity (MW)", "Battery (hr)",
                    "Solar_perc (%)", "Solar_curt (%)", "Wind_perc (%)", "Wind_curt (%)", 
                    "Renewable_perc (%)", "Renewable_curt (%)", "Capacity factor", "Avg water prod (m3/day)",
                    ])
        
    energy_demand = 100 * 1000 # kW
    failed_county = []

    # Solar capacities in different combinations
    MW100 = [25, 50, 75]
    MW150 = [25, 50, 75, 100, 125]
    MW200 = [50, 75, 100, 125, 150]
    MW250 = [50, 75, 100, 125, 150, 175, 200]
    MW300 = [75, 100, 125, 150, 175, 200, 225]
    MW350 = [100, 125, 150, 175, 200, 225, 250]
    combinations = [
                    # MW100, 
                    # MW150, 
                    MW200, 
                    MW250, 
                    MW300, 
                    MW350]
    demands = [
        # 100, 
        # 150, 
        200, 
        250, 
        300, 
        350]
    battery_power_cost = {0: 0,
                          2: 710,
                          4: 1036,
                          6: 1362,
                          8: 1688,
                          10: 2014,
                          }
    wind_functions = {25: wind25,
                      50: wind50,
                      75: wind75,
                      100: wind100,
                      125: wind125,
                      150: wind150,
                      175: wind175,
                      200: wind200,
                      225: wind225,
                      250: wind250,
    }
    grid_prices = {
        4 : 0.0786,
        6 : 0.1709,
        12: 0.0916,
        35: 0.0656,
        48: 0.0713,
    }

    start_time = time.time()
    progress = 0
    for index, row in df.iterrows():
        # try:     
        county = row["County_ID"]
        wind_file = f"wind_{county}.csv"
        weatherfile =  row["filename"]
        SEEC = row["NF SEEC (kWh-e/m3)"] + row["Pumping energy (kWh-e/m3)"] + row["LSRRO SEEC (kWh-e/m3)"] # kWh/m3
        water_prod = energy_demand / SEEC  # m3/hour
        capex_lcow = 0.45 # $/m3
        fixed_opex_lcow = 0.25 # $/m3
        STEC = row["Cryst STEC (kWh-th/m3)"]
        cost_salt = row["Salt recovery ($/m3)"]
        salt_prod = cost_salt / (-0.024) * water_prod # kg/hr
        state = row["State_ID"]
        # print(water_prod)
        
        weatherfile_path = US_solar_directory + weatherfile 
        wind_resource_path = TX_wind_directory + wind_file

        for i in range(len(combinations)):
            demand = demands[i]
            results = {
                "LCOW_no_grid": float("inf"),
                "LCOW ($/m3)": float("inf"),
                "Solar capacity (MW)": 0,
                "Wind capacity (MW)": 0,
                "Battery (hr)": 0,
                "LCOE solar": 0,
                "LCOE wind": 0,
                "LCOE solar adjusted": 0,
                "LCOE wind adjusted": 0,
                "Solar_perc (%)": 0,
                "Solar_curt (%)": 0,
                "Wind_perc (%)": 0,
                "Wind_curt (%)": 0,
                "Renewable_perc (%)": 0,
                "Renewable_curt (%)":0,
                "Capacity factor": 0,
                "Avg water prod (m3/day)":0,
            }
            for solar in combinations[i]:
                wind = demand - solar
                solar_gen, capacity_factor, kwh_per_kw, lcoe_solar, monthly_energy = PVwatts(
                                                                weatherfile_path,
                                                                solar*1000,
                                                                0)
                wind_gen, lcoe_wind = wind_functions[wind](wind_resource=wind_resource_path)

                if len(solar_gen) != len(wind_gen):
                    if len(solar_gen) < len(wind_gen):
                        wind_gen = wind_gen[0: len(solar_gen)-1]
                    else:
                        solar_gen = solar_gen[0: len(wind_gen)-1]
                
                flag = lcoe_solar < lcoe_wind
                # print('')
                for storage in [0,2,4,6,8,10]:
                    solar_load, wind_load, to_desal, solar_curt, wind_curt, curt, water_gen, salt_gen, fuel = simulation_battery_grid(
                                            solar_gen,
                                            wind_gen,
                                            power = 100000, # kW
                                            water_prod = water_prod, # m3/hr,
                                            salt_prod = salt_prod, # kg/hr
                                            flag = flag,
                                            storage = storage,
                                            )
                    solar_perc = sum(solar_load) / 100000 / 8760
                    wind_perc = sum(wind_load) / 100000 / 8760
                    fuel_perc = sum(fuel) / 100000 / 8760

                    solar_curt_perc = sum(solar_curt)/sum(solar_gen)
                    wind_curt_perc = sum(wind_curt)/sum(wind_gen)
                    curt_perc = sum(curt) / (sum(solar_gen) + sum(wind_gen))
                    cf = sum(water_gen)/water_prod/8760

                    lcoe_solar_adjusted = round(lcoe_solar / (1 - solar_curt_perc),3)
                    lcoe_wind_adjusted = round(lcoe_wind / (1 - wind_curt_perc),3)

                # Calculate new LCOW
                    LCOH = 0.02
                    grid_price = grid_prices[state]
                    cost_th = STEC * LCOH # $/m3
                    cost_elec = (sum(solar_load) * lcoe_solar_adjusted + sum(wind_load) * lcoe_wind_adjusted) / sum(water_gen)
                    cost_salt = sum(salt_gen) * (-0.024) / sum(water_gen)

                    battery_cost = (157 * storage + battery_power_cost[storage]) * 100000

                    designed_annual_water_prod = water_prod * 8760 * 0.9
                    annualized_CAPEX = capex_lcow * designed_annual_water_prod + battery_cost * 0.1

                    LCOW_no_grid = (annualized_CAPEX/ (sum(water_gen)*0.9)
                            + fixed_opex_lcow
                            + cost_th
                            + cost_elec
                            + cost_salt
                            )   
                    
                    # Calculate LCOW with grid participation
                    cost_elec = (sum(solar_load) * lcoe_solar_adjusted + sum(wind_load) * lcoe_wind_adjusted + sum(fuel) * grid_price) / water_prod / 8760
                    cost_salt = salt_prod * (-0.024) / water_prod
                    cf = 1
                    LCOW = (annualized_CAPEX/ (water_prod * 8760 * 0.9)
                            + fixed_opex_lcow
                            + cost_th
                            + cost_elec
                            + cost_salt
                            )                
                    # print(f"{solar} solar + {wind} wind + {storage} batt: {round(100*cf,1)}%, {round(100*fuel_perc,1)}%, {round(lcoe_solar_adjusted,3)}, {round(lcoe_wind_adjusted,3)}, {round(LCOW_no_grid,3)}, {round(LCOW,3)};")
                    if LCOW_no_grid < results["LCOW_no_grid"]:
                        results = {
                        "LCOW_no_grid": LCOW_no_grid,
                        "LCOW ($/m3)": LCOW,
                        "Solar capacity (MW)": solar,
                        "Wind capacity (MW)": wind,
                        "Battery (hr)": storage,
                        "LCOE solar": lcoe_solar,
                        "LCOE wind": lcoe_wind,
                        "LCOE solar adjusted": lcoe_solar_adjusted,
                        "LCOE wind adjusted": lcoe_wind_adjusted,
                        "Solar_perc (%)": solar_perc,
                        "Solar_curt (%)": solar_curt_perc,
                        "Wind_perc (%)": wind_perc,
                        "Wind_curt (%)": wind_curt_perc,
                        "Renewable_perc (%)": solar_perc + wind_perc,
                        "Renewable_curt (%)": curt_perc,
                        "Capacity factor": cf,
                        "Avg water prod (m3/day)": sum(water_gen) / 365,
                    }
            
            csv_writers[i].writerow([
                row["County"], row["State_ID"], row["County_ID"], row["Latitude"], row["Longitude"],
                SEEC, water_prod*24,
                results["LCOE solar"], results["LCOE solar adjusted"],
                results["LCOE wind"], results["LCOE wind adjusted"], 
                results["LCOW ($/m3)"], results["Solar capacity (MW)"], results["Wind capacity (MW)"], results["Battery (hr)"],
                results["Solar_perc (%)"], results["Solar_curt (%)"], results["Wind_perc (%)"], results["Wind_curt (%)"],
                results["Renewable_perc (%)"], results["Renewable_curt (%)"], results["Capacity factor"], results["Avg water prod (m3/day)"],
            ])

        progress += 1
        print('Progress: ', progress,"/",len(df), 'Processing time: ', round((time.time()-start_time)/60,1), 'min')
          

def simulation( solar_gen, 
                wind_gen,
                power, # kW
                water_prod,  # m3/hr
                salt_prod,
                flag,
                storage=0):

    gen = [solar_gen[i] + wind_gen[i] for i in range(len(solar_gen))]
    
    
    to_desal = [0 for i in range(len(gen))]
    to_storage =  [0 for i in range(len(gen))]
    storage_load =  [0 for i in range(len(gen))]
    storage_cap_1 =  [0 for i in range(len(gen))]
    storage_cap_2 = [0 for i in range(len(gen))]
    storage_status =  [0 for i in range(len(gen))]
    solar_loss =  [0 for i in range(len(gen))]
    load =  [0 for i in range(len(gen))]
    salt =  [0 for i in range(len(gen))]
    prod =  [0 for i in range(len(gen))]
    fuel =  [0 for i in range(len(gen))]
    energy_consumption =  [0 for i in range(len(gen))]
    actual_load =  [0 for i in range(len(gen))]

    curt = [0 for i in range(len(gen))]
    wind_curt = [0 for i in range(len(gen))]
    solar_curt = [0 for i in range(len(gen))]
    water_gen = [0 for i in range(len(gen))]
    solar_load = [0 for i in range(len(gen))]
    wind_load = [0 for i in range(len(gen))]
    salt_gen = [0 for i in range(len(gen))]
    for i in range(len(gen)):
        
        to_desal[i] = min(power, gen[i])
        curt[i] = gen[i] - to_desal[i]

        if flag:
            solar_curt[i] =min(solar_gen[i], curt[i])
            wind_curt[i] = curt[i] - solar_curt[i]
        else:
            wind_curt[i] = min(wind_gen[i], curt[i])
            solar_curt[i] = curt[i] - wind_curt[i]
        
        solar_load[i] = solar_gen[i] - solar_curt[i]
        wind_load [i] = wind_gen[i] - wind_curt[i]

        water_gen[i] = to_desal[i] / power * water_prod
        salt_gen[i] = to_desal[i] / power * salt_prod
    
    return solar_load, wind_load, to_desal, solar_curt, wind_curt, curt, water_gen, salt_gen

def simulation_battery( solar_gen, 
                wind_gen,
                power, # kW
                water_prod,  # m3/hr
                salt_prod,
                flag,
                storage=0):

    gen = [solar_gen[i] + wind_gen[i] for i in range(len(solar_gen))]

    storage_cap = storage * power # kWh
    
    to_desal = [0 for i in range(len(gen))]
    to_storage =  [0 for i in range(len(gen))]
    storage_load =  [0 for i in range(len(gen))]
    storage_cap_1 =  [0 for i in range(len(gen))]
    storage_cap_2 = [0 for i in range(len(gen))]
    storage_status =  [0 for i in range(len(gen))]
    solar_loss =  [0 for i in range(len(gen))]
    load =  [0 for i in range(len(gen))]
    salt =  [0 for i in range(len(gen))]
    prod =  [0 for i in range(len(gen))]
    fuel =  [0 for i in range(len(gen))]
    energy_consumption =  [0 for i in range(len(gen))]
    actual_load =  [0 for i in range(len(gen))]

    curt = [0 for i in range(len(gen))]
    wind_curt = [0 for i in range(len(gen))]
    solar_curt = [0 for i in range(len(gen))]
    water_gen = [0 for i in range(len(gen))]
    solar_load = [0 for i in range(len(gen))]
    wind_load = [0 for i in range(len(gen))]
    salt_gen = [0 for i in range(len(gen))]
    for i in range(len(gen)):
        to_desal[i] = min(power, gen[i])
        to_storage[i] = abs(gen[i] - to_desal[i])
        storage_load[i] = gen[i] - power
        if i != 0:
            storage_cap_1[i] = storage_status[i-1]
        storage_cap_2[i] = max(storage_load[i] + storage_cap_1[i], 0)
        storage_status[i] = min(storage_cap_2[i] , storage_cap)
        curt[i] = abs(storage_status[i] - storage_cap_2[i])

        load[i] = to_desal[i] + max(0, storage_cap_1[i] - storage_cap_2[i])
        if max(0,load[i] / power) < 0:
            fuel[i] = power - load[i]
        
        if flag:
            solar_curt[i] =min(solar_gen[i], curt[i])
            wind_curt[i] = curt[i] - solar_curt[i]
        else:
            wind_curt[i] = min(wind_gen[i], curt[i])
            solar_curt[i] = curt[i] - wind_curt[i]

        solar_load[i] = solar_gen[i] - solar_curt[i]
        wind_load [i] = wind_gen[i] - wind_curt[i]

        actual_load[i] = max(0,load[i])
        energy_consumption[i] = fuel[i]+load[i]
        water_gen[i] = (fuel[i]+load[i] )/ power * water_prod  
        salt_gen[i] = (fuel[i]+load[i] ) / power * salt_prod
    
    return solar_load, wind_load, to_desal, solar_curt, wind_curt, curt, water_gen, salt_gen

def simulation_battery_grid( solar_gen, 
                wind_gen,
                power, # kW
                water_prod,  # m3/hr
                salt_prod,
                flag,
                storage=0):

    gen = [solar_gen[i] + wind_gen[i] for i in range(len(solar_gen))]

    storage_cap = storage * power # kWh
    
    to_desal = [0 for i in range(len(gen))]
    to_storage =  [0 for i in range(len(gen))]
    storage_load =  [0 for i in range(len(gen))]
    storage_cap_1 =  [0 for i in range(len(gen))]
    storage_cap_2 = [0 for i in range(len(gen))]
    storage_status =  [0 for i in range(len(gen))]
    solar_loss =  [0 for i in range(len(gen))]
    load =  [0 for i in range(len(gen))]
    salt =  [0 for i in range(len(gen))]
    prod =  [0 for i in range(len(gen))]
    fuel =  [0 for i in range(len(gen))]
    energy_consumption =  [0 for i in range(len(gen))]
    actual_load =  [0 for i in range(len(gen))]

    curt = [0 for i in range(len(gen))]
    wind_curt = [0 for i in range(len(gen))]
    solar_curt = [0 for i in range(len(gen))]
    water_gen = [0 for i in range(len(gen))]
    solar_load = [0 for i in range(len(gen))]
    wind_load = [0 for i in range(len(gen))]
    salt_gen = [0 for i in range(len(gen))]
    for i in range(len(gen)):
        to_desal[i] = min(power, gen[i])
        to_storage[i] = abs(gen[i] - to_desal[i])
        storage_load[i] = gen[i] - power
        if i != 0:
            storage_cap_1[i] = storage_status[i-1]
        storage_cap_2[i] = max(storage_load[i] + storage_cap_1[i], 0)
        storage_status[i] = min(storage_cap_2[i] , storage_cap)
        curt[i] = abs(storage_status[i] - storage_cap_2[i])

        load[i] = to_desal[i] + max(0, storage_cap_1[i] - storage_cap_2[i])
        if max(0,load[i] / power) < 1:
            fuel[i] = power - load[i]
        
        if flag:
            solar_curt[i] =min(solar_gen[i], curt[i])
            wind_curt[i] = curt[i] - solar_curt[i]
        else:
            wind_curt[i] = min(wind_gen[i], curt[i])
            solar_curt[i] = curt[i] - wind_curt[i]

        solar_load[i] = solar_gen[i] - solar_curt[i]
        wind_load [i] = wind_gen[i] - wind_curt[i]

        water_gen[i] = (load[i] )/ power * water_prod  
        salt_gen[i] = (load[i] ) / power * salt_prod
    
    return solar_load, wind_load, to_desal, solar_curt, wind_curt, curt, water_gen, salt_gen, fuel


def retrieve_elec_price(lat, lon):
    # url_test = "https://apps.openei.org/IURDB/rate/view/5bbe36f35457a3b467135b51#3__Energy"
    url, eneryg_rates = lookup_rates(lat, lon)
    print(url)
    print(energy_rates)
    if url:
        try:
            page = requests.get(url)

            soup = BeautifulSoup(page.text, 'lxml')

            table1 = soup.find('section', id = "energy_rate_strux_table")
            table2 = list(table1.find_all("div", {"class": "strux_view_cell"}))

            rates = []
            adj = []
            for i in range(len(table2)):
                if i > 7 and i%7 == 4:
                    rates.append(float(table2[i].string))    
                if i > 7 and i%7 == 5:
                    adj.append(float(table2[i].string))
                    # Average elec rates   # From EIA?  # Tiers
            return [sum(rates)/len(rates), 'Yes',       len(rates)]
        except:
            return [0.05, 'No', 1]
    
    else:
        return [0.05, 'No', 1 ]

if __name__ == "__main__":
    # PV_wind_cal()
    # PV_wind_battery_cal()
    PV_wind_battery_grid_cal()
