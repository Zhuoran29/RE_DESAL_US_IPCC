from OARO import OARO
from OARO_cost import OARO_cost
from PVWatts import PVwatts

import os
import time
import csv
import pandas as pd
import sys

from urllib.request import urlopen
from lxml import etree
from bs4 import BeautifulSoup
import requests
import json

from lookup_openei_rates import lookup_rates

def PV_OARO_cal():
    state_elec_price_path = 'state_ind_retail_price_2021.csv'
    df_elec = pd.read_csv(state_elec_price_path, encoding= 'unicode_escape') 

    state_elec_price = {}
    for index, row in df_elec.iterrows():
        state_elec_price[row["Name"]] = row['Industrial retail price ($/kWh)']

    US_solar_directory = 'D:/PhD/DOE/Sensitivity/US_solar_data/'
    files = os.listdir(US_solar_directory)

    test_csv = "test.csv"
    desal_csv = "desal_all_sectors.csv"

    df = pd.read_csv(desal_csv, encoding= 'unicode_escape')   

    # print(df.keys())

    csv_outfile = 'D:/PhD/Papers/New/Modeling_code/PV_OARO_results.csv' 
    data_file = open(csv_outfile, 'w', newline='') 
    csv_writer = csv.writer(data_file)
    
    csv_writer.writerow(["Location", "State","Latitude", "Longitude", "Cap_increase", "OARO_Capacity (m3/day)", \
                         "PV_Capacity (kW)", "Solar multiple", "Storage/battery (hour)", "STEC (kWh/m3 feed)", \
                         "Thermal power (MW)", "LCOE from solar field($/kWh)", "LCOE after curtailment adjustment ($/kWh)", \
                         "LCOE from grid", "From EIA", "LCOW ($/m3)", "PV_perc (%)", "Grid_perc (%)", "Curtail_perc (%)", "Sector"\
                         ])
    start_time = time.time()
    progress = 0
    for index, row in df.iterrows():
        try:
            if not row["Capacity__"]:
                continue
            elif "Municipalities" in row["Customer_t"] or "Tourist" in row["Customer_t"]:
                try:
                    cap_increase = row["dp_dry_50"]
                    sector = "Municipality"
                except:
                    continue
            elif "Industry" in row["Customer_t"]:
                try:
                    cap_increase = row["i_dry_50"]
                    sector = "Industry"
                except:
                    continue                
            elif "Power" in row["Customer_t"]:
                try:
                    cap_increase = row["th_dry_50"]
                    sector = "Thermoelectricity"
                except:
                    continue     
            elif "Irrigation" in row["Customer_t"]:
                try:
                    cap_increase = row["ir_dry_50"]
                    sector = "Irrigation"
                except:
                    continue   
            else:
                try:
                    cap_increase = row["all_dry_50"]
                    sector = "Others"
                except:
                    continue                                

            latitude = row["Latitude"]
            longitude = row["Longitude"]
            location = row["Location"]
            state = row["State_Regi"]
            weatherfile =  row["filename"]
            RO_cap = row["Capacity__"]
            # Check water utility price projections

            # if "Seawater" in row["Feedwater"]:
            #     tds = 30
            # elif "Brine" in row["Feedwater"]:
            #     tds = 50
            # elif "Pure" in row["Feedwater"]:
            #     tds = 5
            # elif "Brackish" in row["Feedwater"]:
            #     tds = 20
            # elif "River" in row["Feedwater"]:
            #     tds = 10
            # else:
            #     tds = 30

            if "Brackish" not in row["Feedwater"]:
                continue

            # Model simulation
            OARO_new_cap = RO_cap * (1 + cap_increase)
            OARO_sys = OARO(FeedC_r = 20, # Feed concentration (g/L)
                            Capacity = OARO_new_cap, # System capacity  (m3/day)
                            rr = 1,)

            OARO_design_outputs, OARO_cost_outputs = OARO_sys.design()
            brine_volume = OARO_new_cap * (1 - 0.92) # m3/day
            crystallizer_P_req = brine_volume * 157 / 24 / 1000   # MW-th
            STEC = crystallizer_P_req * 24 * 1000 / OARO_new_cap  # kWh-th / m3 feed

            SEC = OARO_design_outputs[-2]["Value"]
            P_req = OARO_design_outputs[1]["Value"]  #MW

            weatherfile_path = US_solar_directory + weatherfile        

            try:
                grid_rate = state_elec_price[state] / 100
            except:
                grid_rate = state_elec_price["U.S. Total"] / 100

            from_EIA = "Yes"
            tiers = 1
            # grid_rate, from_EIA, tiers = retrieve_elec_price(latitude, longitude

            results = solar_multiple_optimizer(OARO_sys, weatherfile_path, P_req, latitude, SEC, grid_rate, OARO_new_cap, OARO_cost_outputs, False)

            curt_perc = results["Curtailment percentage (%)"]
            LCOW = results["LCOW"]
            storage = results["Storage"]
            sm = results["Solar multiple"]
            LCOE_solar = results["LCOE from solar field"]
            LCOE_curt = results["LCOE adjusted from curtailment"]
            grid_perc = results["Grid percentage"]

            csv_writer.writerow([location, state, latitude, row["Longitude"], cap_increase, OARO_new_cap, P_req*sm, sm, storage, STEC, crystallizer_P_req,  LCOE_solar, LCOE_curt,\
                                grid_rate, from_EIA, LCOW, 100-grid_perc, grid_perc, curt_perc, sector])

            progress += 1
            print('Progress: ', round(progress/len(df.index)*100,1), '%', 'Processing time: ', round((time.time()-start_time)/60,1), 'min')
    

        except BaseException as e:
            print('Fail')
            print(row["Location"] + row["State_Regi"] + str(e))
            continue


def solar_multiple_optimizer(OARO_sys, weatherfile_path, P_req, latitude, SEC, grid_rate, OARO_new_cap, OARO_cost_outputs, monthly):
    solar_multiples = [1,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3]
    storage_hours = [0] #[0, 4, 8, 12]
    results = {"Curtailment percentage (%)": 0, "LCOW": float("inf"), "Solar multiple": 1, "Storage": 0, "LCOE from solar field": float("inf"), "Grid percentage": 0}

    try:
        for sm in solar_multiples:
            for storage in storage_hours:
                PV_capacity = P_req * sm
                gen, capacity_factor, kwh_per_kw, lcoe_fcr, monthly_energy = PVwatts(weatherfile_path, PV_capacity, latitude)
                RO_simulation, monthly_grid_consumed, monthly_elec_curtailed = OARO_sys.simulation(gen = gen, storage = storage)

                curt_perc = RO_simulation[-1]["Value"]
                grid_perc = RO_simulation[-3]["Value"]
                prod = RO_simulation[4]["Value"]
                grid_use = RO_simulation[-4]["Value"]

                lcoe_after_curtail = lcoe_fcr / (100 - curt_perc) * 100
                # print(OARO_new_cap, prod, grid_perc, SEC, lcoe_fcr, grid_rate)
                # print('cost started')
                cost_model = OARO_cost(Capacity = OARO_new_cap, Prod = prod, fuel_usage=grid_perc,
                            oaro_area  = OARO_cost_outputs['oaro_area'], ro_area  = OARO_cost_outputs['ro_area'], 
                            sec = SEC, 
                            monthly_input = monthly,
                            sam_coe=lcoe_fcr, 
                            solar_coe=None, 
                            coe = grid_rate, 
                            # monthly_grid_coe = grid_rate, 
                            storage_cap = OARO_sys.storage_cap,
                            monthly_grid_consumed = monthly_grid_consumed,
                            monthly_elec_curtailed = monthly_elec_curtailed)
                            
                # print('cost initialized')
                cost = cost_model.lcow()
                # print('cost ended')
                LCOW = cost[2]['Value']
                # print('sm = ', sm, 'LCOW = ', round(LCOW, 4), 'elec_cost = ', round(cost_model.cost_elec,4), 'LCOE = ', round(lcoe_after_curtail ,4), 'grip_perc = ', round(grid_perc,1))

                # Update value if a lower LCOW is achieved
                if LCOW < results["LCOW"]:
                    results = {"Curtailment percentage (%)": curt_perc, "LCOW": LCOW, "Solar multiple": sm,  "Storage": storage, "LCOE from solar field": lcoe_fcr, "LCOE adjusted from curtailment":lcoe_after_curtail , "Grid percentage": grid_perc}

                    # if curt_perc < 15:
                    #     if LCOW < results["LCOW"]:
                    #         results = {"Curtailment percentage (%)": curt_perc, "LCOW": LCOW, "Solar multiple": sm, "LCOE from solar field": lcoe_fcr, "Grid percentage": grid_perc}
                    # else:
                    #     if results["LCOW"] > 1e9:
                    #         results = {"Curtailment percentage (%)": curt_perc, "LCOW": LCOW, "Solar multiple": sm, "LCOE from solar field": lcoe_fcr, "Grid percentage": grid_perc}
    except:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno)

    return results

if __name__ == "__main__":
    PV_OARO_cal()
    # PV_RO_monthly_cal()
    # print(retrieve_elec_price(40.7472, -79.9959))