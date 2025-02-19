from RO_Fixed_Load import RO
from RO_cost import RO_cost
from PVWatts import PVwatts

import sys
import os
import time
import csv
import pandas as pd

from urllib.request import urlopen
# from lxml import etree
from bs4 import BeautifulSoup
import requests
import json

from lookup_openei_rates import lookup_rates

def PV_RO_cal():
    state_elec_price_path = 'state_ind_retail_price_2021.csv'
    df_elec = pd.read_csv(state_elec_price_path, encoding= 'unicode_escape')

    state_elec_price = {}
    for index, row in df_elec.iterrows():
        state_elec_price[row["State_ID"]] = row['Industrial retail price ($/kWh)']

    US_solar_directory = '/Users/zhuoranzhang/Documents/DOE_Chile/DOE_CSP_PROJECT/SAM_flatJSON/solar_resource'
    files = os.listdir(US_solar_directory)
    coastal_csv = "coastal_county.csv"

    test_csv = "test.csv"
    desal_csv = "county_center_dry_ssp5_TableToExcel.csv"

    df = pd.read_csv(coastal_csv, encoding= 'unicode_escape')   

    # print(df.keys())

    csv_outfile = 'coastal_PV_RO_2050_new.csv' 
    data_file = open(csv_outfile, 'w', newline='') 
    csv_writer = csv.writer(data_file)
    
    csv_writer.writerow(["County", "State_ID", "County_ID", "Latitude", "Longitude", 
                         "Brackish_TDS (ppm)",  "Brackish_depth (feet)", "Brackish_availability (%)", "RO_recovery",
                         "RO_Capacity (m3/day)", "PV_Capacity (kW)", "Solar multiple", "Storage/battery (hour)", 
                         "LCOE from solar field($/kWh)", "LCOE after curtailment adjustment ($/kWh)", 
                         "LCOE from grid", "From EIA", "LCOW ($/m3)", "PV_perc (%)", "Grid_perc (%)", "Curtail_perc (%)",
                         "RO SEC (kWh/m3)", "Pumping energy (kWh/m3)", "RO CAPEX ($)", "Drilling cost ($)", 
                         ])
    start_time = time.time()
    progress = 0
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
            RO_cap = (row["dp_8_dry_Y2050"] - row["dp_8_dry_Y2015"]) * 3785.4
            
            if not RO_cap > 0:
                # print('county', county, 'skipped')
                progress += 1
                continue
            tds = 35166.4

            if tds == 0:
                tds = 5001

            depth = 0
            brackish = row["brackish_availability_new_volume_mgd"]
            # Model simulation

            # Rocovery rate is selected to reach 60 g/L brine salinity
            # if brackish > 0:
            RO_rr = 1 - tds / 1000 / 70
            # else:
            #     RO_rr = 1 - 5 / 50

            brackish_avail = brackish * 3785.4 / (RO_cap / RO_rr) * 100
            pump_energy_brackish = depth * 0.3048 * 9.8 * 1000 * 2.778 * 1e-7 / RO_rr / 0.85 # kWh/m3 

            RO_sys = RO(nominal_daily_cap_tmp = RO_cap,
                        FeedC_r = tds / 1000,
                        R1 = RO_rr * 100)

            RO_outputs = RO_sys.RODesign()
            SEC = RO_outputs[-1]["Value"]
            P_req = RO_outputs[-2]["Value"] + pump_energy_brackish / 24 * RO_cap

            weatherfile_path = US_solar_directory + weatherfile        

            try:
                grid_rate = state_elec_price[state_id] / 100
            except:
                grid_rate = state_elec_price[0] / 100

            from_EIA = "Yes"
            tiers = 1
            # grid_rate, from_EIA, tiers = retrieve_elec_price(latitude, longitude)

            results = solar_multiple_optimizer(RO_sys, weatherfile_path, P_req, latitude, SEC, grid_rate, RO_cap, False, RO_rr, depth)

            curt_perc = results["Curtailment percentage (%)"]
            LCOW = results["LCOW"]
            storage = results["Storage"]
            sm = results["Solar multiple"]
            LCOE_solar = results["LCOE from solar field"]
            LCOE_curt = results["LCOE adjusted from curtailment"]
            grid_perc = results["Grid percentage"]
            RO_SEC = results["RO SEC (kWh/m3)"]
            pumping_energy = results["Pumping energy (kWh/m3)"]
            RO_capex = results["RO EPC ($)"]
            drilling_cost = results["Drilling cost ($)"]

            csv_writer.writerow([county, state_id, county_id, latitude, longitude, 
                                    tds, depth, brackish_avail, RO_rr,
                                    RO_cap, P_req*sm, sm, storage, 
                                    LCOE_solar, LCOE_curt,
                                    grid_rate, from_EIA, LCOW, 100-grid_perc, grid_perc, curt_perc,
                                    RO_SEC, pumping_energy, RO_capex, drilling_cost,
                                    ])

            progress += 1
            print('Progress: ', round(progress/len(df.index)*100,1), '%', 'Processing time: ', round((time.time()-start_time)/60,1), 'min', 'sm:', sm,'gridprice: ',grid_rate )
            

        except BaseException as e:
            print('Fail')
            print(row["cb_2017_us_county_500k_NAME"] + str(e))
            progress += 1
            continue

def solar_multiple_optimizer(RO_sys, weatherfile_path, P_req, latitude, SEC, grid_rate, RO_new_cap, monthly, RO_rr, depth ):
    solar_multiples = [3]
    storage_hours = [0, 2,4,6,8,10]
    results = {"Curtailment percentage (%)": 0, "LCOW": float("inf"), "Solar multiple": 1,  "Storage": 0, 
                "LCOE from solar field": float("inf"), "LCOE adjusted from curtailment":float("inf") , "Grid percentage": 0,
                "RO SEC (kWh/m3)": 0, "Pumping energy (kWh/m3)": 0,
                "RO EPC ($)": 0, "Drilling cost ($)": 0}
    # try:
    for sm in solar_multiples:
        PV_capacity = P_req * sm
        gen, capacity_factor, kwh_per_kw, lcoe_fcr, monthly_energy = PVwatts(weatherfile_path, PV_capacity, 0)

        for storage in storage_hours:
            RO_simulation, monthly_grid_consumed, monthly_elec_curtailed = RO_sys.simulation(gen = gen, storage = storage)

            curt_perc = RO_simulation[-1]["Value"]
            grid_perc = RO_simulation[-3]["Value"]
            prod = RO_simulation[4]["Value"]
            grid_use = RO_simulation[-4]["Value"]

            lcoe_after_curtail = lcoe_fcr / (100 - curt_perc) * 100 * 0.391
            # print(RO_new_cap, prod, grid_perc, SEC, lcoe_fcr, grid_rate)
            cost_model = RO_cost(Capacity = [RO_new_cap], 
                                Prod = prod, 
                                fuel_usage=grid_perc, 
                                sec = SEC, 
                                monthly_input = monthly,
                                sam_coe=lcoe_after_curtail, 
                                solar_coe=None, 
                                ave_grid_coe = grid_rate, 
                                monthly_grid_coe = grid_rate, 
                                storage_cap = RO_sys.storage_cap,
                                monthly_grid_consumed = monthly_grid_consumed, 
                                monthly_elec_curtailed = monthly_elec_curtailed, 
                                RO_rr = RO_rr, 
                                brackish_depth=depth,
                                cost_storage = 159,)
                        
            cost = cost_model.lcow()

            LCOW = cost[2]['Value']
            # print('sm = ', sm, 'LCOW = ', round(LCOW, 4), 'elec_cost = ', round(cost_model.cost_elec,4), 'LCOE = ', round(lcoe_after_curtail ,4), 'grip_perc = ', round(grid_perc,1))

            # Update value if a lower LCOW is achieved
            if LCOW < results["LCOW"]:
                results = {"Curtailment percentage (%)": curt_perc, "LCOW": LCOW, "Solar multiple": sm,  "Storage": 0, 
                            "LCOE from solar field": lcoe_fcr, "LCOE adjusted from curtailment":lcoe_after_curtail , "Grid percentage": grid_perc,
                            "RO SEC (kWh/m3)": cost_model.SEC, "Pumping energy (kWh/m3)": cost_model.pump_brackish,
                            "RO EPC ($)": cost_model.EPC_cost, "Drilling cost ($)": cost_model.drilling_cost}

                # if curt_perc < 15:
                #     if LCOW < results["LCOW"]:
                #         results = {"Curtailment percentage (%)": curt_perc, "LCOW": LCOW, "Solar multiple": sm, "LCOE from solar field": lcoe_fcr, "Grid percentage": grid_perc}
                # else:
                #     if results["LCOW"] > 1e9:
                #         results = {"Curtailment percentage (%)": curt_perc, "LCOW": LCOW, "Solar multiple": sm, "LCOE from solar field": lcoe_fcr, "Grid percentage": grid_perc}
    # except:
    #     exc_type, exc_obj, exc_tb = sys.exc_info()
    #     fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    #     print(exc_type, fname, exc_tb.tb_lineno)
    return results



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
    PV_RO_cal()
    # PV_RO_monthly_cal()
    # print(retrieve_elec_price(40.7472, -79.9959))