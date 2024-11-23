from VAGMD_batch.VAGMD_batch import VAGMD_batch
from MDB_cost import MDB_cost
from StaticCollector_flatplate import StaticCollector_fp
from SC_ETC_cost import ETC_cost
from PVWatts import PVwatts

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

solver = get_solver()

import warnings
warnings.filterwarnings("ignore")


def NF_LSRRO_cryst_cal():
    state_elec_price_path = 'state_ind_retail_price_2021.csv'
    df_elec = pd.read_csv(state_elec_price_path, encoding= 'unicode_escape')

    state_elec_price = {}
    for index, row in df_elec.iterrows():
        state_elec_price[row["State_ID"]] = row['Industrial retail price ($/kWh)']

    US_solar_directory = '/Users/zhuoranzhang/Documents/DOE_Chile/DOE_CSP_PROJECT/SAM_flatJSON/solar_resource/'
    files = os.listdir(US_solar_directory)

    test_csv = "test.csv"
    major_ions_csv = 'TDS_county_filtered.csv'
    desal_csv = "county_center_dry_ssp5_TableToExcel.csv"
    failed_csv = "failed_counties.csv"

    df = pd.read_csv(test_csv, encoding= 'unicode_escape')   

    # print(df.keys())

    csv_outfile1 = 'county_dry_ssp5_2050_PV_NF_LSRRO_cryst_sm1.csv'
    # csv_outfile2= 'county_dry_ssp5_2070_IPH_MD_cryst_sm14.csv' 
    csv_outfile3 = 'county_dry_ssp5_2050_PV_NF_LSRRO_cryst_sm18.csv' 
    # csv_outfile4 = 'county_dry_ssp5_2070_IPH_MD_cryst_sm22.csv' 
    # csv_outfile5 = 'county_dry_ssp5_2070_IPH_MD_cryst_sm26.csv'
    csv_outfile6 = 'county_dry_ssp5_2050_PV_NF_LSRRO_cryst_sm3.csv'

    data_file1 = open(csv_outfile1, 'w', newline='') 
    # data_file2 = open(csv_outfile2, 'w', newline='') 
    data_file3 = open(csv_outfile3, 'w', newline='') 
    # data_file4 = open(csv_outfile4, 'w', newline='') 
    # data_file5 = open(csv_outfile5, 'w', newline='') 
    data_file6 = open(csv_outfile6, 'w', newline='') 

    csv_writer1 = csv.writer(data_file1)
    # csv_writer2 = csv.writer(data_file2)
    csv_writer3 = csv.writer(data_file3)
    # csv_writer4 = csv.writer(data_file4)
    # csv_writer5 = csv.writer(data_file5)
    csv_writer6 = csv.writer(data_file6)

    csv_writers = [csv_writer1, 
                #    csv_writer2, 
                   csv_writer3, 
                #    csv_writer4, 
                #    csv_writer5, 
                   csv_writer6]

    for i in csv_writers:
        i.writerow(["County", "State_ID", "County_ID", "Latitude", "Longitude", 
                    "Brackish_TDS (ppm)",  "Brackish_depth (feet)", "Brackish_availability (%)", 
                    "NF_Capacity (m3/day)", "LSRRO_recovery", "LSRRO_Capacity (m3/day)", 
                    "Cryst capacity (m3/day)", "PV_Capacity (kW)", "Solar multiple", "Storage/battery (hour)", 
                    "LCOE from solar field($/kWh)", "LCOE after curtailment adjustment ($/kWh)", 
                    "LCOE from grid", "From EIA", "LCOW ($/m3)", "PV_perc (%)", "Grid_perc (%)", "Curtail_perc (%)",
                    "NF SEEC (kWh-e/m3)", "LSRRO SEEC (kWh-e/m3)", 
                    "Cryst SEEC (kWh-e/m3)", "Cryst STEC (kWh-th/m3)",  "Pumping energy (kWh-e/m3)",
                    "NF CAPEX ($)", "LSRRO CAPEX ($)", "Cryst CAPEX ($)", "Drilling cost ($)", 
                    "Elec cost ($/m3)", "Heat cost ($/m3)", "Salt recovery ($/m3)"
                    ])

    start_time = time.time()
    progress = 0
    skipped_county = []
    failed_county = []
    half_failed = []
    for index, row in df.iterrows():
        try:      
            latitude = row["county_lat"]
            longitude = row["county_lon"]
            county = row["NAME"]
            county_id = row["GEOID"]
            # state = row["State"]
            state_id = row["STATEFP"]
            weatherfile =  row["filename"]
            # RO capacity is to provide additional water demand
            desal_cap = (row["dp_8_dry_Y2050"] - row["dp_8_dry_Y2015"]) * 3785.4

            depth = row["brackish_availability_new_depth_ft"]
            brackish = row["brackish_availability_new_volume_mgd"]

            if not desal_cap > 0:
                # print('county', county, 'skipped')
                progress += 1
                skipped_county.append(county_id)
                continue

            tds = row["MEAN_TDS_mgL"]

            solute_dict = {
                'Na': row["MEAN_Na_mgL"],
                'Mg': row["MEAN_Mg_mgL"],
                'Ca': row["MEAN_Ca_mgL"],
                'K': row["MEAN_K_mgL"],
                'Cl': row["MEAN_Cl_mgL"],
                'CO3': row["MEAN_CO3_mgL"],
                'HCO3': row["MEAN_HCO3_mgl"],
                'SO4': row["MEAN_SO4_mgL"],
            }


            if tds == 0:
                tds = 5001
                depth = 1000
                continue

            # Model simulation

            # ADD DESAL
            m, lsrro, cryst, data_table, overall_performance,capacities = fs_NF_LSRRO_crys(
                                                                            desal_cap,
                                                                            solute_dict)


            SEEC_cryst = value(pyunits.convert(
            (
                cryst.fs.eff_1.magma_circulation_flow_vol
                * cryst.fs.eff_1.dens_mass_slurry
                * Constants.acceleration_gravity
                * cryst.fs.eff_1.costing.costing_package.crystallizer.pump_head_height
                / cryst.fs.eff_1.costing.costing_package.crystallizer.efficiency_pump
            ),
            to_units=pyunits.kW,
            )) / (overall_performance["Total water prod vol (m3/s)"] * 3600) # m3/h condensation
            
            # brine_flow_rate * 24  # kWh-e/m3 feed
            lsrro_rr = capacities[-1]
            STEC_cryst = overall_performance["Overall STEC (kWh/m3 feed)"] * (1 - lsrro_rr) / (0.753 + 0.247 * lsrro_rr)

            SEEC_lsrro = value(lsrro.fs.costing.specific_energy_consumption)

            lsrro_prod = capacities[2] # m3/day
            cryst_prod = overall_performance["Total water prod vol (m3/s)"] * 86400 # m3/day

            total_water_prod = lsrro_prod + cryst_prod
            SEEC_nf = value(m.fs.zo_costing.aggregate_flow_electricity) / value(m.fs.NF.properties_treated[0].flow_vol) / 3600

            # Pump energy
            brackish_avail = brackish * 3785.4 / (desal_cap / lsrro_rr / 0.85) * 100
                                                                            #  NF RR # pump eff
            pump_energy_brackish = depth * 0.3048 * 9.8 * 1000 * 2.778 * 1e-7 / 0.85 / 0.85 # kWh/m3 
            SEEC_avg = (lsrro_prod * SEEC_lsrro + cryst_prod * SEEC_cryst) / total_water_prod + SEEC_nf + pump_energy_brackish

            # Total power
            Elec_req = SEEC_avg * desal_cap / 24 + pump_energy_brackish / 24 * desal_cap # kW
            Thermal_req = overall_performance["Initial thermal energy consumption (kW)"] # kW

            weatherfile_path = US_solar_directory + weatherfile        

            try:
                grid_rate = 0.07 #state_elec_price[state_id] / 100
            except:
                grid_rate = state_elec_price[0] / 100

            from_EIA = "Yes"
            tiers = 1
            # grid_rate, from_EIA, tiers = retrieve_elec_price(latitude, longitude)

            info = (county, state_id, county_id, 
            latitude, longitude, tds, depth, 
            brackish_avail, desal_cap, 
            capacities[3], Thermal_req, grid_rate, from_EIA,
            )
            solar_multiple_optimizer(
                                    m, 
                                    lsrro,
                                    cryst,
                                    data_table,
                                    overall_performance, 
                                    weatherfile_path, 
                                    SEEC_avg,
                                    [SEEC_nf, SEEC_lsrro, SEEC_cryst, pump_energy_brackish],
                                    Thermal_req, 
                                    Elec_req,
                                    lsrro_prod,
                                    total_water_prod,
                                    latitude,  
                                    grid_rate, 
                                    desal_cap, 
                                    depth, 
                                    STEC_cryst,
                                    cryst_prod,
                                    capacities,
                                    csv_writers,
                                    info,
                                    )


            progress += 1
            print('Progress: ', round(progress/len(df.index)*100,1), '%', 'Processing time: ', round((time.time()-start_time)/60,1), 'min', 'failed: ',len(failed_county) )
            

        except BaseException as e:
            # print(county_id, desal_cap)
            failed_county.append(county_id)
            progress += 1
            continue
    
    return skipped_county, failed_county, half_failed

def solar_multiple_optimizer(
                            m, 
                            lsrro,
                            cryst,
                            data_table,
                            overall_performance, 
                            weatherfile_path, 
                            SEEC_avg,
                            SEECs,
                            Thermal_req, 
                            Elec_req,
                            lsrro_prod,
                            total_water_prod,
                            latitude, 
                            grid_rate, 
                            desal_cap, 
                            depth, 
                            STEC_cryst,
                            cryst_prod,
                            capacities,
                            csv_writers,
                            info,
                              ):
    solar_multiples = [1, 1.8, 3]
    storage_hours = [ 0, 4, 8, 12, 16]

    SEEC_nf, SEEC_lsrro, SEEC_cryst, pump_energy_brackish = SEECs
        
    try:
        for i in range(len(solar_multiples)):
            sm = solar_multiples[i]
            results = {"Curtailment percentage (%)": 0, "LCOW": float("inf"), "Solar multiple": 1,  "Storage": 0, 
                "LCOE from solar field": float("inf"), "LCOE adjusted from curtailment":float("inf") , "Grid percentage": 0,
                "NF SEEC (kWh-e/m3)": 0, "LSRRO SEEC (kWh-e/m3)": 0,
                "Cryst SEEC (kWh-e/m3)": 0, "Cryst STEC (kWh-th/m3)": 0, "Pumping energy (kWh-e/m3)": 0,
                "NF CAPEX ($)": 0, "LSRRO CAPEX ($)": 0, "Cryst CAPEX ($)": 0, "Drilling cost ($)": 0,
                "Elec cost ($/m3)": 0 , "Heat cost ($/m3)": 0, "Salt recovery ($/m3)": 0}

            PV_capacity = Elec_req * sm * 1.1
            gen, capacity_factor, kwh_per_kw, lcoe_fcr, monthly_energy = PVwatts(weatherfile_path,
                                                                                PV_capacity,
                                                                                0)
            for storage in storage_hours:
                simu_output = simulation(
                                        gen = gen, 
                                        power = Elec_req, # kWh
                                        prod  = total_water_prod,  # m3/day
                                        overall_performance = overall_performance,
                                        storage= storage)

                curt_perc = simu_output[-1]["Value"]
                grid_perc = simu_output[-3]["Value"]
                water_prod = simu_output[5]["Value"]
                salt_prod = simu_output[6]["Value"]
                grid_use = simu_output[-4]["Value"]

                feed_vol, nf_cap, lsrro_cap, cryst_cap, lsrro_rr = capacities

                #  CAPEX
                NF_capex = value(pyunits.convert(value(pyunits.convert(m.fs.zo_costing.total_capital_cost, to_units=pyunits.USD_2018)) * pyunits.USD_2018, to_units = pyunits.USD_2021))  / ( value(m.fs.NF.properties_treated[0].flow_vol) * 86400) * nf_cap
                LSRRO_capex = value(lsrro.fs.costing.total_capital_cost) / (value(lsrro.fs.product.properties[0].flow_vol) * 86400) * lsrro_cap
                cryst_capex = value(cryst.fs.costing.total_capital_cost)
                drilling_cost = 650 * depth * desal_cap / 3785.4 / 3.26

                # Fixed OPEX
                NF_fix_opex = value(pyunits.convert(value(pyunits.convert(m.fs.zo_costing.total_fixed_operating_cost, to_units=pyunits.USD_2018 / pyunits.year)) * pyunits.USD_2018, to_units = pyunits.USD_2021)) / ( value(m.fs.NF.properties_treated[0].flow_vol) * 86400) * nf_cap
                LSRRO_fix_opex = value(cryst.fs.costing.total_fixed_operating_cost) / (value(lsrro.fs.product.properties[0].flow_vol) * 86400) * lsrro_cap
                cryst_fix_opex = value(cryst.fs.costing.total_fixed_operating_cost)

                # variable OPEX
                lcoh = 0.05 # $/kWh
                NaCl_rec = -0.024 # $/kg

                cost_th = STEC_cryst * lcoh  # $/m3
                cost_elec = SEEC_avg * (grid_perc/100 * grid_rate + (1- grid_perc/100) * lcoe_fcr)
                cost_NaCl = (overall_performance['Total solids collected (kg/s)'] * 86400) / desal_cap * NaCl_rec

                # Water Prod (m3/day)

                total_capex = NF_capex + LSRRO_capex + cryst_capex + drilling_cost
                total_fix_opex = NF_fix_opex + LSRRO_fix_opex + cryst_fix_opex


                LCOW = ((total_capex * 0.1 
                    + total_fix_opex) 
                    / ((lsrro_cap + cryst_cap * 0.753 ) * 365 * 0.9)
                    + cost_th
                    + cost_elec
                    + cost_NaCl
                )


                lcoe_after_curtail = lcoe_fcr / (100 - curt_perc) * 100


                # Update value if a lower LCOW is achieved
                if LCOW < results["LCOW"]:
                    results = {"Curtailment percentage (%)": curt_perc, "LCOW": LCOW, "Solar multiple": sm,  "Storage": storage, 
                                "LCOE from solar field": lcoe_fcr, "LCOE adjusted from curtailment":lcoe_after_curtail , "Grid percentage": grid_perc,
                                "NF SEEC (kWh-e/m3)": SEEC_nf, "LSRRO SEEC (kWh-e/m3)": SEEC_lsrro,
                                "Cryst SEEC (kWh-e/m3)": SEEC_cryst,"Cryst STEC (kWh-th/m3)": STEC_cryst,  "Pumping energy (kWh-e/m3)": pump_energy_brackish,
                                "NF CAPEX ($)": NF_capex, "LSRRO CAPEX ($)": LSRRO_capex, "Cryst CAPEX ($)": cryst_capex, "Drilling cost ($)": drilling_cost,
                                "Elec cost ($/m3)": cost_elec , "Heat cost ($/m3)": cost_th, "Salt recovery ($/m3)": cost_NaCl
                                }

            curt_perc = results["Curtailment percentage (%)"]
            LCOW = results["LCOW"]
            storage = results["Storage"]
            sm = results["Solar multiple"]
            LCOE_solar = results["LCOE from solar field"]
            LCOE_curt = results["LCOE adjusted from curtailment"]
            grid_perc = results["Grid percentage"]
            NF_SEEC = results["NF SEEC (kWh-e/m3)"]
            LSRRO_SEEC = results["LSRRO SEEC (kWh-e/m3)"]
            cryst_STEC = results["Cryst STEC (kWh-th/m3)"]
            cryst_SEEC = results["Cryst SEEC (kWh-e/m3)"]
            pumping_energy = results["Pumping energy (kWh-e/m3)"]
            NF_capex = results["NF CAPEX ($)"]
            LSRRO_capex = results["LSRRO CAPEX ($)"]
            cryst_capex = results["Cryst CAPEX ($)"]
            drilling_cost = results["Drilling cost ($)"]

            (county, state_id, county_id, 
            latitude, longitude, tds, depth, 
            brackish_avail, desal_cap, 
            cryst_cap, Elec_req, grid_rate, from_EIA,
            ) = info

            csv_writers[i].writerow([county, state_id, county_id, latitude, longitude, 
                                    tds, depth, brackish_avail, 
                                    nf_cap, lsrro_rr, desal_cap, 
                                    cryst_cap, Elec_req*sm, sm, storage, 
                                    LCOE_solar, LCOE_curt,
                                    grid_rate, from_EIA, LCOW, 100-grid_perc, grid_perc, curt_perc,
                                    NF_SEEC, LSRRO_SEEC,
                                    cryst_SEEC, cryst_STEC, pumping_energy, 
                                    NF_capex, LSRRO_capex, cryst_capex, drilling_cost,
                                    cost_elec, cost_th, cost_NaCl
                                    ])


    except:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno)

    return results


def simulation( gen, 
                power, # kWh
                prod,  # m3/day
                overall_performance, # cryst performance
                storage=0):

    cryst_load = overall_performance["Initial thermal energy consumption (kW)"]
    salt_prod = overall_performance["Total solids collected (kg/s)"] * 3600 # kg/hr    
    
    power_load = power # kWh
    max_prod = prod / 24 # m3/h
    storage_cap = storage * power_load # kWh
    
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
    for i in range(len(gen)):
        to_desal[i] = min(power_load, gen[i])
        to_storage[i] = abs(gen[i] - to_desal[i])
        storage_load[i] = gen[i] - power_load
        if i != 0:
            storage_cap_1[i] = storage_status[i-1]
        storage_cap_2[i] = max(storage_load[i] + storage_cap_1[i], 0)
        storage_status[i] = min(storage_cap_2[i] , storage_cap)
        solar_loss[i] = abs(storage_status[i] - storage_cap_2[i])
        load[i] = to_desal[i] + max(0, storage_cap_1[i] - storage_cap_2[i])
        if max(0,load[i] / power_load) < 1:
            fuel[i] = power_load - load[i]
        
        actual_load[i] = max(0,load[i])
        energy_consumption[i] = fuel[i]+load[i]
        prod[i] = (fuel[i]+load[i] )/ power_load * max_prod  
        salt[i] = (fuel[i]+load[i] )/ power_load * salt_prod
    
    Month = [0,31,59,90,120,151,181,212,243,273,304,334,365]
    Monthly_prod = [ sum( prod[Month[i]*24:(Month[i+1]*24)] ) for i in range(12) ]
    monthly_grid_consumed = [ sum( fuel[Month[i]*24:(Month[i+1]*24)] ) / sum( energy_consumption[Month[i]*24:(Month[i+1]*24)] )  for i in range(12)]
    monthly_elec_curtailed = [ (sum( gen[Month[i]*24:(Month[i+1]*24)]) - sum(load[Month[i]*24:(Month[i+1]*24)] ))/sum( gen[Month[i]*24:(Month[i+1]*24)] ) for i in range(12)]

    simu_output = []

    simu_output.append({'Name':'Water production','Value':prod,'Unit':'m3'})
    simu_output.append({'Name':'Salt production','Value':prod,'Unit':'kg'})
    simu_output.append({'Name':'Storage status','Value':storage_status,'Unit':'kWh'})
    simu_output.append({'Name':'Storage Capacity','Value':storage_cap,'Unit':'kWh'})
    simu_output.append({'Name':'Fossil fuel usage','Value':fuel,'Unit':'kWh'})
    simu_output.append({'Name':'Total water production','Value':sum(prod),'Unit':'m3'})
    simu_output.append({'Name':'Total salt production','Value':sum(salt),'Unit':'kg'})
    simu_output.append({'Name':'Monthly water production','Value': Monthly_prod,'Unit':'m3'})
    simu_output.append({'Name':'Total fossil fuel usage','Value':sum(fuel),'Unit':'kWh'})
    simu_output.append({'Name':'Percentage of fossil fuel consumption','Value':sum(fuel)/max(1,sum(energy_consumption))*100,'Unit':'%'})          
    simu_output.append({'Name':'Curtailed solar electric energy','Value':max(0, (sum(gen) - sum(load)) / 1000000) ,'Unit':'GWh'})   
    simu_output.append({'Name':'Percentage of curtailed energy','Value':max(0, (sum(gen) - sum(load))) / sum(gen) * 100 ,'Unit':'%'}) 
    
    return simu_output


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

# Retreive monthly retail price for each state (Industrial sector)
def EIA_api():
    url = 'https://api.eia.gov/v2/electricity/retail-sales/data/?api_key=NrYesSRYplmSurhuPtNPqbqC0CwAKtifeBQ7dXX0'
    x = '&data[]=price&frequency=monthly&facets[sectorid][]=IND&start=2020-12-31&end=2021-12-31'
    page = requests.get(url+x)
    data = json.loads(page.text)
    

    state_elec_price_path = 'state_ind_retail_price_2021.csv'
    df_elec = pd.read_csv(state_elec_price_path, encoding= 'unicode_escape')

    idx = ['2021-01', '2021-02','2021-03','2021-04','2021-05','2021-06','2021-07','2021-08',
           '2021-09','2021-10','2021-11','2021-12',]
    monthly_df = pd.DataFrame(index = idx,  columns = df_elec['Name'].tolist())
    for i in data['response']['data']:
        try:
            monthly_df.at[i['period'], i['stateDescription']] = i['price'] / 100
        except:
            print(i['period'], i['stateDescription'])

    return monthly_df

if __name__ == "__main__":
    skipped, failed, half_failed = NF_LSRRO_cryst_cal()

    print('skipped', len(skipped))
    print('failed', len(failed))
    print('half failed', len(half_failed))

