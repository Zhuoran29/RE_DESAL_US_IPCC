from VAGMD_batch.VAGMD_batch import VAGMD_batch
from MDB_cost import MDB_cost
from StaticCollector_flatplate import StaticCollector_fp
from SC_ETC_cost import ETC_cost
from IPH_LF import IPH_LF

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

from watertap_reflo.src.watertap_contrib.reflo.analysis.example_flowsheets.multi_effect_NaCl_crystallizer import (
    build_fs_multi_effect_crystallizer,
    add_costings,
    multi_effect_crystallizer_initialization,
    get_model_performance,
)
import watertap_reflo.src.watertap_contrib.reflo.property_models.cryst_prop_pack as props
from watertap.core.solvers import get_solver
from watertap_reflo.src.watertap_contrib.reflo.costing import TreatmentCosting, CrystallizerCostType

solver = get_solver()

import warnings
warnings.filterwarnings("ignore")


def IPH_MD_cryst_cal():
    state_elec_price_path = 'state_ind_retail_price_2021.csv'
    df_elec = pd.read_csv(state_elec_price_path, encoding= 'unicode_escape')

    state_elec_price = {}
    for index, row in df_elec.iterrows():
        state_elec_price[row["State_ID"]] = row['Industrial retail price ($/kWh)']

    US_solar_directory = '/Users/zhuoranzhang/Documents/DOE_Chile/DOE_CSP_PROJECT/SAM_flatJSON/solar_resource/'
    files = os.listdir(US_solar_directory)

    test_csv = "test.csv"
    desal_csv = "county_center_dry_ssp5_TableToExcel.csv"
    failed_csv = "failed_counties.csv"

    df = pd.read_csv(desal_csv, encoding= 'unicode_escape')   

    # print(df.keys())

    csv_outfile1 = 'county_dry_ssp5_2050_IPH_MD_cryst_sm1.csv'
    # csv_outfile2= 'county_dry_ssp5_2070_IPH_MD_cryst_sm14.csv' 
    csv_outfile3 = 'county_dry_ssp5_2050_IPH_MD_cryst_sm18.csv' 
    # csv_outfile4 = 'county_dry_ssp5_2070_IPH_MD_cryst_sm22.csv' 
    # csv_outfile5 = 'county_dry_ssp5_2070_IPH_MD_cryst_sm26.csv'
    csv_outfile6 = 'county_dry_ssp5_2050_IPH_MD_cryst_sm3.csv'

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
                    "Brackish_TDS (ppm)",  "Brackish_depth (feet)", "Brackish_availability (%)", "RO_recovery",
                    "MD_Capacity (m3/day)", "Cryst capacity (m3/day)", "FPC_Capacity (kW)", "Solar multiple", "Actual SM", "Thermal storage (hour)", 
                    "LCOH from solar field($/kWh)", "LCOH after curtailment adjustment ($/kWh)", 
                    "LCOE from grid", "From EIA", "LCOW ($/m3)", "PV_perc (%)", "Grid_perc (%)", "Curtail_perc (%)",
                    "MD STEC (kWh-th/m3)", "MD SEEC (kWh-e/m3)", 
                    "Cryst STEC (kWh-th/m3)", "Cryst SEEC (kWh-e/m3)", "Pumping energy (kWh-e/m3)",
                    "MD CAPEX ($)", "Cryst CAPEX ($)", "Drilling cost ($)", 
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
            county = row["cb_2017_us_county_500k_NAME"]
            county_id = row["cb_2017_us_county_500k_GEOID"]
            # state = row["State"]
            state_id = row["cb_2017_us_county_500k_STATEFP"]
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
            tds = row["brackish_availability_new_TDS_ppm"]



            if tds == 0:
                tds = 5001
                depth = 1000


            # Model simulation

            # Rocovery rate is selected to reach 180 g/L brine salinity
            # if brackish > 0:
            desal_rr = 1 - tds / 1000 / 180
            # else:
            #     desal_rr = 1 - 5 / 50

            brackish_avail = brackish * 3785.4 / (desal_cap / desal_rr) * 100
            pump_energy_brackish = depth * 0.3048 * 9.8 * 1000 * 2.778 * 1e-7 / desal_rr / 0.85 # kWh/m3 

            MDB = VAGMD_batch(
                            module = 1,  # '0' for AS7C1.5L module and '1' for AS26C2.7L module
                            TEI_r  = 80, # Evaporator channel inlet temperature (ºC)
                            TCI_r  = 25, # Condenser channel inlet temperature (ºC)
                            FFR_r  = 1100, # Feed flow rate (l/h)
                            FeedC_r= tds/1000,  # Feed concentration (g/L)
                            V0     = 50,  # Initial batch volume (L)
                            RR     = desal_rr * 100,  # Recovery rate (%)
                            
                            Capacity = desal_cap, # System Capcity (m3/day)
                            Fossil_f = 1, # Fossil fuel fraction
                            
                            j = 'c', # cooling system: 'c' for closed, 'o' for open
                            Ttank = 25, # Initial temeprature of the saline feed
                            TCoolIn = 15, # Initial tmeeprature of the cooling water
                            dt = 60, # Time step for the simulation (< 480 second)
                            )

            MDB_outputs = MDB.design()

            brine_flow_rate = desal_cap * (1 - desal_rr) # m3/day
            brine_temp = MDB.Ttank[-1]
            brine_density = 1110.50 # kg/m3 @ 180 g/L, 35 C


            # Add crystallizer
            m = build_fs_multi_effect_crystallizer(
                operating_pressure_eff1=0.4455,  # bar
                operating_pressure_eff2=0.2758,  # bar
                operating_pressure_eff3=0.1651,  # bar
                operating_pressure_eff4=0.095,  # bar
                feed_flow_mass= 1,  # kg/s
                feed_mass_frac_NaCl=0.1613,
                feed_pressure=101325,  # Pa
                feed_temperature=273.15 + brine_temp,  # K
                crystallizer_yield=0.8,
                steam_pressure=1.5,  # bar (gauge pressure)
            )
            add_costings(m)

            # Negative value for salt recovery value ($/kg)
            
            m.fs.costing.crystallizer.steam_cost.fix(0)
            m.fs.costing.crystallizer.NaCl_recovery_value.fix(-0.024)

            multi_effect_crystallizer_initialization(m)

            results = solver.solve(m)

            eff_1 = m.fs.eff_1

            eff_1.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].unfix()
            eff_1.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].unfix()
            eff_1.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"].fix()


            @m.Constraint(eff_1.flowsheet().time, doc="total of feed flow")
            def total_feed(b, t):
                effs = [m.fs.eff_1, m.fs.eff_2, m.fs.eff_3, m.fs.eff_4]
                return (
                    sum(i.properties_in[0].flow_vol_phase["Liq"] for i in effs)
                    == brine_flow_rate / 86400 * pyunits.m**3 / pyunits.s
                )
            try:
                results = solver.solve(m)
            except:
            # Update scaling factors
                for eff in [m.fs.eff_1, m.fs.eff_2, m.fs.eff_3, m.fs.eff_4]:

                    iscale.set_scaling_factor(eff.properties_in[0].flow_mass_phase_comp["Liq", "H2O"], 1e-1/1000)
                    iscale.set_scaling_factor(eff.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"], 1e-1/1000)
                    iscale.set_scaling_factor(eff.properties_in[0].flow_mass_phase_comp["Vap", "H2O"], 1e-1/1000)
                    iscale.set_scaling_factor(eff.properties_in[0].flow_mass_phase_comp["Sol", "NaCl"], 1e-1/1000)
                    iscale.set_scaling_factor(eff.energy_flow_superheated_vapor, 1e-4)
                    for ind, c in eff.eq_vapor_energy_constraint.items():
                        sf = iscale.get_scaling_factor(eff.energy_flow_superheated_vapor)
                        iscale.constraint_scaling_transform(c, sf)

                    iscale.set_scaling_factor(eff.pressure_operating, 1e-3)
                    iscale.set_scaling_factor(eff.properties_out[0].pressure, 1e-5)
                    iscale.set_scaling_factor(eff.properties_solids[0].pressure, 1e-5)
                    iscale.set_scaling_factor(eff.properties_vapor[0].pressure, 1e-3)
                    iscale.set_scaling_factor(eff.properties_out[0].flow_vol_phase["Liq"], 1e3)
                    iscale.set_scaling_factor(eff.properties_out[0].flow_vol_phase["Vap"], 1e8)
                    iscale.set_scaling_factor(eff.properties_out[0].pressure_sat, 1e-3)
                    iscale.set_scaling_factor(eff.properties_solids[0].flow_vol_phase["Vap"], 1e8)
                    iscale.set_scaling_factor(eff.properties_solids[0].flow_vol_phase["Sol"], 1e8)
                    iscale.set_scaling_factor(
                        eff.properties_vapor[0].dens_mass_solvent["Vap"], 1e1
                    )
                    iscale.set_scaling_factor(eff.properties_out[0].flow_vol_phase["Sol"], 1e12)
                    iscale.set_scaling_factor(
                        eff.properties_solids[0].flow_vol_phase["Liq"], 1e11
                    )
                    iscale.set_scaling_factor(eff.properties_vapor[0].flow_vol_phase["Liq"], 1e12)
                    iscale.set_scaling_factor(eff.properties_vapor[0].flow_vol_phase["Sol"], 1e12)

                results = solver.solve(m)
                half_failed.append(county_id)


            data_table, overall_performance = get_model_performance(m)

            elec_cryst = value(pyunits.convert(
            (
                m.fs.eff_1.magma_circulation_flow_vol
                * m.fs.eff_1.dens_mass_slurry
                * Constants.acceleration_gravity
                * m.fs.eff_1.costing.costing_package.crystallizer.pump_head_height
                / m.fs.eff_1.costing.costing_package.crystallizer.efficiency_pump
            ),
            to_units=pyunits.kW,
            )) / brine_flow_rate * 24  # kWh-e/m3 feed


            STEC = MDB_outputs[-4]["Value"] + overall_performance["Overall STEC (kWh/m3 feed)"] * (1 - desal_rr)
            SEEC = MDB_outputs[-3]["Value"] + elec_cryst

            Thermal_req = MDB.P_req + overall_performance["Initial thermal energy consumption (kW)"] # kW
            Elec_req = SEEC * desal_cap / 24 + pump_energy_brackish / 24 * desal_cap

            cryst_water_prod = overall_performance["Total water prod vol (m3/s)"] * 86400 # m3/day

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
            brackish_avail, desal_rr, desal_cap, 
            brine_flow_rate, Thermal_req, grid_rate, from_EIA,
            )
            solar_multiple_optimizer(MDB, 
                                    m, 
                                    overall_performance, 
                                    weatherfile_path, 
                                    Thermal_req, 
                                    latitude, 
                                    STEC, 
                                    grid_rate, 
                                    desal_cap, 
                                    False, 
                                    desal_rr,
                                    depth, 
                                    elec_cryst,
                                    cryst_water_prod,
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

def solar_multiple_optimizer(MDB, 
                             cryst, 
                             overall_performance, 
                             weatherfile_path, 
                             P_req, 
                             latitude, 
                             STEC, 
                             grid_rate, 
                             desal_cap, 
                             monthly, 
                             desal_rr, 
                             depth, 
                             cryst_elec,
                             cryst_water_prod,
                             csv_writers,
                             info,
                              ):
    solar_multiples = [1, 1.8, 3]
    storage_hours = [ 0, 4, 8, 12, 16]


        
    try:
        for i in range(len(solar_multiples)):
            sm = solar_multiples[i]
            results = {"Curtailment percentage (%)": 0, "LCOW": float("inf"), "Solar multiple": 1,  "actual sm": 1, "Storage": 0, 
                "LCOH from solar field": float("inf"), "LCOH adjusted from curtailment":float("inf") , "Grid percentage": 0,
                "MD STEC (kWh-th/m3)": 0, "MD SEEC (kWh-e/m3)": 0,
                "Cryst STEC (kWh-th/m3)": 0, "Cryst SEEC (kWh-e/m3)": 0, "Pumping energy (kWh-e/m3)": 0,
                "MD EPC ($)": 0, "Cryst EPC ($)": 0, "Drilling cost ($)": 0}
            IPH_capacity = P_req * sm

            energy, capacity_factor, lcoe_fcr, annual_energy, solarm = IPH_LF(
                                                        weatherfile_path,
                                                        capacity = P_req,
                                                        specified_solar_multiple = sm,
                                                        )
            for storage in storage_hours:

                simu_output = MDB.simulation(gen = energy, storage = storage, overall_performance = overall_performance)

                curt_perc = simu_output[-1]["Value"]
                grid_perc = simu_output[-3]["Value"]
                water_prod = simu_output[5]["Value"]
                salt_prod = simu_output[6]["Value"]
                grid_use = simu_output[-4]["Value"]

                lcoh = lcoe_fcr
                lcoh_after_curtail = lcoh / (100 - curt_perc) * 100


                HX_capex = value(cryst.fs.capex_heat_exchanger + cryst.fs.capex_end_plates)

                cost_model = MDB_cost(Capacity = desal_cap, 
                                        Prod = water_prod, 
                                        fuel_usage=grid_perc, 
                                        Area = MDB.Area,
                                        Pflux = MDB.PFlux_avg,
                                        RR = MDB.R[-1],
                                        TCO = sum(MDB.TCO) / len(MDB.TCO),
                                        TEI = MDB.TEI_r,
                                        FFR = MDB.FFR_r,
                                        th_module = sum(MDB.ThPower) / len(MDB.ThPower),
                                        STEC = MDB.STEC[-1],
                                        SEEC = MDB.SEEC[-1],
                                        sam_coh=lcoh_after_curtail, 
                                        coh = 0.03,
                                        coe =  grid_rate,
                                        cost_storage = 15,
                                        storage_cap = MDB.storage_cap,
                                        brackish_depth=depth,
                                        cryst = cryst,
                                        cryst_elec = cryst_elec,
                                        cryst_th = overall_performance["Overall STEC (kWh/m3 feed)"],
                                        HX_capex = HX_capex,
                                        salt_prod = salt_prod,
                                        salt_cost = -0.024,
                                        cryst_water_prod = cryst_water_prod,
                                        )
                            
                cost = cost_model.lcow()

                LCOW = cost[2]['Value']

                # Update value if a lower LCOW is achieved
                if LCOW < results["LCOW"]:
                    results = {"Curtailment percentage (%)": curt_perc, "LCOW": LCOW, "Solar multiple": sm, "actual sm": solarm, "Storage": storage, 
                                "LCOH from solar field": lcoh, "LCOH adjusted from curtailment":lcoh_after_curtail , "Grid percentage": grid_perc,
                                "MD STEC (kWh-th/m3)": cost_model.STEC, "MD SEEC (kWh-e/m3)": cost_model.SEEC,
                                "Cryst STEC (kWh-th/m3)": overall_performance["Overall STEC (kWh/m3 feed)"], 
                                "Cryst SEEC (kWh-e/m3)": cryst_elec, "Pumping energy (kWh-e/m3)": cost_model.pump_brackish,
                                "MD EPC ($)": cost_model.cost_sys*1000, "Cryst EPC ($)": cost_model.cryst_CAPEX, "Drilling cost ($)": cost_model.drilling_cost
                                }
                    # if curt_perc < 15:
                    #     if LCOW < results["LCOW"]:
                    #         results = {"Curtailment percentage (%)": curt_perc, "LCOW": LCOW, "Solar multiple": sm, "LCOE from solar field": lcoe_fcr, "Grid percentage": grid_perc}
                    # else:
                    #     if results["LCOW"] > 1e9:
                    #         results = {"Curtailment percentage (%)": curt_perc, "LCOW": LCOW, "Solar multiple": sm, "LCOE from solar field": lcoe_fcr, "Grid percentage": grid_perc}

            curt_perc = results["Curtailment percentage (%)"]
            LCOW = results["LCOW"]
            storage = results["Storage"]
            sm = results["Solar multiple"]
            actual_sm = results['actual sm']
            LCOH_solar = results["LCOH from solar field"]
            LCOH_curt = results["LCOH adjusted from curtailment"]
            grid_perc = results["Grid percentage"]
            MD_STEC = results["MD STEC (kWh-th/m3)"]
            MD_SEEC = results["MD SEEC (kWh-e/m3)"]
            cryst_STEC = results["Cryst STEC (kWh-th/m3)"]
            cryst_SEEC = results["Cryst SEEC (kWh-e/m3)"]
            pumping_energy = results["Pumping energy (kWh-e/m3)"]
            MD_capex = results["MD EPC ($)"]
            cryst_capex = results["Cryst EPC ($)"]
            drilling_cost = results["Drilling cost ($)"]

            (county, state_id, county_id, 
            latitude, longitude, tds, depth, 
            brackish_avail, desal_rr, desal_cap, 
            brine_flow_rate, Thermal_req, grid_rate, from_EIA,
            ) = info

            csv_writers[i].writerow([county, state_id, county_id, latitude, longitude, 
                                    tds, depth, brackish_avail, desal_rr,
                                    desal_cap, brine_flow_rate, Thermal_req*sm, sm, actual_sm, storage, 
                                    LCOH_solar, LCOH_curt,
                                    grid_rate, from_EIA, LCOW, 100-grid_perc, grid_perc, curt_perc,
                                    MD_STEC, MD_SEEC,
                                    cryst_STEC, cryst_SEEC, pumping_energy, 
                                    MD_capex, cryst_capex, drilling_cost,
                                    ])


    except:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno)

    # return results



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
    skipped, failed, half_failed = IPH_MD_cryst_cal()

    print('skipped', len(skipped))
    print('failed', len(failed))
    print('half failed', len(half_failed))

