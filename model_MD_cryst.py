from VAGMD_batch.VAGMD_batch import VAGMD_batch
from MDB_cost import MDB_cost
from StaticCollector_flatplate import StaticCollector_fp
from SC_ETC_cost import ETC_cost
from idaes.core import FlowsheetBlock
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
    constraint_scaling_transform,
)
import idaes.core.util.scaling as iscale
from urllib.request import urlopen
# from lxml import etree
from bs4 import BeautifulSoup
import requests
import json

from idaes.core import UnitModelCostingBlock
from lookup_openei_rates import lookup_rates
from pyomo.environ import (
    ConcreteModel,
    TerminationCondition,
    check_optimal_termination,
    SolverStatus,
    Objective,
    Expression,
    Constraint,
    maximize,
    value,
    Set,
    Var,
    log,
    units as pyunits,
)
from watertap.unit_models.zero_order import NanofiltrationZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.costing.zero_order_costing import ZeroOrderCosting

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


def FPC_MD_cryst_cal():
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
    _3sectors_csv = 'TDS_county_3_sectors.csv'
    coastal_csv = "coastal_county.csv"

    df = pd.read_csv(_3sectors_csv, encoding= 'unicode_escape')   

    # print(df.keys())

    # csv_outfile1 = 'county_dry_ssp5_2050_FPC_MD_cryst_sm1.csv'
    # csv_outfile2= 'county_dry_ssp5_2070_IPH_MD_cryst_sm14.csv' 
    csv_outfile3 = 'dp_FPC_MD_cryst_2070_sm18.csv' 
    # csv_outfile4 = 'county_dry_ssp5_2070_IPH_MD_cryst_sm22.csv' 
    # csv_outfile5 = 'county_dry_ssp5_2070_IPH_MD_cryst_sm26.csv'
    csv_outfile6 = 'dp_FPC_MD_cryst_2070_sm3.csv'

    # data_file1 = open(csv_outfile1, 'w', newline='') 
    # data_file2 = open(csv_outfile2, 'w', newline='') 
    data_file3 = open(csv_outfile3, 'w', newline='') 
    # data_file4 = open(csv_outfile4, 'w', newline='') 
    # data_file5 = open(csv_outfile5, 'w', newline='') 
    data_file6 = open(csv_outfile6, 'w', newline='') 

    # csv_writer1 = csv.writer(data_file1)
    # csv_writer2 = csv.writer(data_file2)
    csv_writer3 = csv.writer(data_file3)
    # csv_writer4 = csv.writer(data_file4)
    # csv_writer5 = csv.writer(data_file5)
    csv_writer6 = csv.writer(data_file6)

    csv_writers = [
                    # csv_writer1, 
                #    csv_writer2, 
                   csv_writer3, 
                #    csv_writer4, 
                #    csv_writer5, 
                   csv_writer6]

    for i in csv_writers:
        i.writerow(["County", "State_ID", "County_ID", "Latitude", "Longitude", 
                    "Brackish_TDS (g/L)",  "Brackish_depth (feet)", "Brackish_availability (%)",
                    "MDB_recovery", "MDB_Capacity (m3/day)", 
                    "Cryst capacity (m3/day)", "FPC_Capacity (kW)", "Solar multiple", "Storage (hour)", 
                    "LCOH from solar field($/kWh)", "LCOH after curtailment adjustment ($/kWh)", 
                    "LCOE from grid", "From EIA", "LCOW ($/m3)", "PV_perc (%)", "Grid_perc (%)", "Curtail_perc (%)",
                    "MDB SEEC (kWh-e/m3)", "MDB STEC (kWh-th/m3)", 
                    "Cryst SEEC (kWh-e/m3)", "Cryst STEC (kWh-th/m3)",  "Pumping energy (kWh-e/m3)",
                    "MDB CAPEX ($)", "Cryst CAPEX ($)", "Drilling cost ($)", 
                    "Elec cost ($/m3)", "Heat cost ($/m3)", "Salt recovery ($/m3)",
                    "CAPEX LCOW ($/m3)","Fixed OPEX LCOW ($/m3)", "Energy LCOW ($/m3)"
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
            desal_cap = (row["dp_8_dry_Y2070"] - row["dp_8_dry_Y2015"]) * 3785.4
            # desal_cap = (row["3_sectors_Y2070"] - row["3_sectors_Y2015"]) * 3785.4

            depth = row["brackish_availability_new_depth_ft"]
            brackish = row["brackish_availability_new_volume_mgd"]

            if not desal_cap > 0:
                # print('county', county, 'skipped')
                progress += 1
                skipped_county.append(county_id)
                continue

            tds = row["MEAN_TDS_mgL"]
            # tds = 35000

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
                depth = 2000 

            # Model simulation

            # Rocovery rate is selected to reach 180 g/L brine salinity
            # if brackish > 0:
            # desal_rr = 1 - tds / 1000 / 180
            # else:
            #     desal_rr = 1 - 5 / 50

            NF_permeate_tds = tds / 1000

            # Calculate capacity of each components
            cbrine = 180
            MDB_rr = 1 - NF_permeate_tds / cbrine
            
            MDB_cap = desal_cap / (0.247 + 0.753 / MDB_rr) # m3/day
            cryst_cap = MDB_cap * (1 - MDB_rr) / MDB_rr # m3/day

            nf_cap = MDB_cap / MDB_rr
            feed_vol = nf_cap / 0.85

            # Add MDB component
            MDB = VAGMD_batch(
                            module = 1,  # '0' for AS7C1.5L module and '1' for AS26C2.7L module
                            TEI_r  = 80, # Evaporator channel inlet temperature (ºC)
                            TCI_r  = 25, # Condenser channel inlet temperature (ºC)
                            FFR_r  = 1100, # Feed flow rate (l/h)
                            FeedC_r= NF_permeate_tds,  # Feed concentration (g/L)
                            V0     = 50,  # Initial batch volume (L)
                            RR     = MDB_rr * 100,  # Recovery rate (%)
                            
                            Capacity = MDB_cap, # System Capcity (m3/day)
                            Fossil_f = 1, # Fossil fuel fraction
                            
                            j = 'c', # cooling system: 'c' for closed, 'o' for open
                            Ttank = 25, # Initial temeprature of the saline feed
                            TCoolIn = 15, # Initial tmeeprature of the cooling water
                            dt = 15, # Time step for the simulation (< 480 second)
                            )

            MDB_outputs = MDB.design()

            brine_flow_rate = MDB_cap * (1 - MDB_rr) / MDB_rr# m3/day
            brine_temp = MDB.Ttank[-1]
            brine_density = 1110.50 # kg/m3 @ 180 g/L, 35 C


            # Add crystallizer
            m2 = build_fs_multi_effect_crystallizer(
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
            add_costings(m2)

            # Negative value for salt recovery value ($/kg)
            
            m2.fs.costing.crystallizer.steam_cost.fix(0)
            m2.fs.costing.crystallizer.NaCl_recovery_value.fix(-0.024)

            multi_effect_crystallizer_initialization(m2)

            results = solver.solve(m2)

            eff_1 = m2.fs.eff_1

            eff_1.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].unfix()
            eff_1.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].unfix()
            eff_1.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"].fix()


            @m2.Constraint(eff_1.flowsheet().time, doc="total of feed flow")
            def total_feed(b, t):
                effs = [m2.fs.eff_1, m2.fs.eff_2, m2.fs.eff_3, m2.fs.eff_4]
                return (
                    sum(i.properties_in[0].flow_vol_phase["Liq"] for i in effs)
                    == brine_flow_rate / 86400 * pyunits.m**3 / pyunits.s
                )
                
            try:
                results = solver.solve(m2)
            except:
            # Update scaling factors
                for eff in [m2.fs.eff_1, m2.fs.eff_2, m2.fs.eff_3, m2.fs.eff_4]:

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

                results = solver.solve(m2)
                half_failed.append(county_id)

            data_table, overall_performance = get_model_performance(m2)

            # Collect SEEC and STEC
            SEEC_cryst = value(pyunits.convert(
            (
                m2.fs.eff_1.magma_circulation_flow_vol
                * m2.fs.eff_1.dens_mass_slurry
                * Constants.acceleration_gravity
                * m2.fs.eff_1.costing.costing_package.crystallizer.pump_head_height
                / m2.fs.eff_1.costing.costing_package.crystallizer.efficiency_pump
            ),
            to_units=pyunits.kW,
            )) / (overall_performance["Total water prod vol (m3/s)"] * 3600) # m3/h condensation
            
            STEC_cryst = overall_performance["Overall STEC (kWh/m3 feed)"] * (1 - MDB_rr) / (0.753 + 0.247 * MDB_rr)

            STEC_MDB = MDB_outputs[-4]["Value"]
            SEEC_MDB = MDB_outputs[-3]["Value"] 


            # Pump energy
            brackish_avail = brackish * 3785.4 / (desal_cap / MDB_rr / 0.99) * 100
                                                                            #  NF RR # pump eff
            pump_energy_brackish = depth * 0.3048 * 9.8 * 1000 * 2.778 * 1e-7 / 0.99 / 0.85 # kWh/m3 

            # Get SEEC/STEC wrt total produced water
            cryst_prod = overall_performance["Total water prod vol (m3/s)"] * 86400 # m3/day

            SEEC_avg = (MDB_cap * SEEC_MDB + cryst_prod * SEEC_cryst) / (cryst_prod + MDB_cap) + pump_energy_brackish
            STEC_avg = STEC_MDB + STEC_cryst

            # Total power
            Thermal_req = MDB.P_req + overall_performance["Initial thermal energy consumption (kW)"] # kW
            Elec_req = SEEC_avg * desal_cap / 24  # kW

            cryst_water_prod = overall_performance["Total water prod vol (m3/s)"] * 86400 # m3/day


            weatherfile_path = US_solar_directory + weatherfile        

            try:
                grid_rate = state_elec_price[state_id] / 100
            except:
                grid_rate = state_elec_price[0] / 100

            from_EIA = "Yes"
            tiers = 1
            # grid_rate, from_EIA, tiers = retrieve_elec_price(latitude, longitude)

            info = (county, state_id, county_id, 
            latitude, longitude, tds, depth, 
            brackish_avail, desal_cap, 
            cryst_cap, Thermal_req, grid_rate, from_EIA,
            )
            solar_multiple_optimizer(MDB, 
                                    m2, 
                                    overall_performance, 
                                    weatherfile_path, 
                                    latitude, 
                                    grid_rate, 
                                    desal_cap, 
                                    False, 
                                    MDB_rr,
                                    depth, 
                                    SEEC_avg,
                                    STEC_avg,
                                    [SEEC_MDB, SEEC_cryst, pump_energy_brackish],
                                    [STEC_MDB, STEC_cryst],
                                    MDB_cap,
                                    cryst_cap,
                                    nf_cap,
                                    Thermal_req,
                                    Elec_req,
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
                             latitude, 
                             grid_rate, 
                             desal_cap, 
                             monthly, 
                             MDB_rr, 
                             depth,
                             SEEC_avg,
                             STEC_avg,
                             SEECs,
                             STECs, 
                             MDB_cap,
                             cryst_cap,
                             nf_cap,
                             Thermal_req,
                             Elec_req,
                             csv_writers,
                             info,
                              ):
    solar_multiples = [1.8 , 3]
    storage_hours = [ 0, 4, 8, 12, 16]


    try:
        for i in range(len(solar_multiples)):
            sm = solar_multiples[i]
            results = {"Curtailment percentage (%)": 0, "LCOW": float("inf"), "Solar multiple": 1,  "Storage": 0, 
                "LCOH from solar field": float("inf"), "LCOH adjusted from curtailment":float("inf") , "Grid percentage": 0,
                "MDB SEEC (kWh-e/m3)": 0, "MDB STEC (kWh-th/m3)": 0,
                "Cryst SEEC (kWh-e/m3)": 0, "Cryst STEC (kWh-th/m3)": 0, "Pumping energy (kWh-e/m3)": 0,
                "LSRRO CAPEX ($)": 0, "Cryst CAPEX ($)": 0, "Drilling cost ($)": 0,
                "Elec cost ($/m3)": 0 , "Heat cost ($/m3)": 0, "Salt recovery ($/m3)": 0}
            
            FPC_capacity = Thermal_req * sm
            
            FPC = StaticCollector_fp(desal_thermal_power_req = FPC_capacity/1000,
                                        file_name= weatherfile_path,
                                        tilt_angle = latitude)

            heat_gen, sc_output = FPC.design()

            SEEC_MDB, SEEC_cryst, pump_energy_brackish = SEECs
            STEC_MDB, STEC_cryst = STECs

            for storage in storage_hours:

                simu_output = MDB.simulation(gen = heat_gen, storage = storage, overall_performance = overall_performance)

                curt_perc = simu_output[-1]["Value"]
                grid_perc = simu_output[-3]["Value"]
                water_prod = simu_output[5]["Value"]
                salt_prod = simu_output[6]["Value"]
                grid_use = simu_output[-4]["Value"]

                FPC_cost = ETC_cost(aperture_area= sc_output[3]['Value'],
                        thermal_gen=sum(heat_gen),
                        P_req=MDB.design_output[-6]["Value"],
                        yrs = 25,
                        capacity = FPC_capacity,
                        )

                lcoh = FPC_cost.lcoh()[-1]["Value"]
                lcoh_after_curtail = lcoh / (100 - curt_perc) * 100

                # Calculate MDB cost
                cost_model = MDB_cost(Capacity = MDB_cap, 
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
                                        cost_storage = 0,
                                        storage_cap = MDB.storage_cap,
                                        )
                cost_model.lcow()

                #  CAPEX
                MDB_capex = cost_model.cost_sys*1000
                cryst_capex = value(cryst.fs.costing.total_capital_cost)
                drilling_cost = 650 * depth * desal_cap / 3785.4 / 3.26
                TES_cost = 15 * storage * Thermal_req 

                # Fixed OPEX
                cryst_fix_opex = value(cryst.fs.costing.total_fixed_operating_cost)

                # variable OPEX
                external_lcoh = 0.02 # $/kWh
                NaCl_rec = -0.024 # $/kg

                cost_th = STEC_avg * (grid_perc/100 * external_lcoh + (1- grid_perc/100) * lcoh_after_curtail)  # $/m3
                cost_elec = SEEC_avg * grid_rate
                cost_NaCl = (overall_performance['Total solids collected (kg/s)'] * 86400) / desal_cap * NaCl_rec

                MDB_var_opex = (cost_model.cost_module_re + cost_model.maintenancecost + cost_model.insurancecost + cost_model.operationalcost)

                # Water Prod (m3/day)

                total_capex =  MDB_capex + cryst_capex + drilling_cost + TES_cost
                total_fix_opex =  cryst_fix_opex

                # Calculate LCOW ($/m3)
                LCOW = ((total_capex * 0.1 
                    + total_fix_opex) 
                    / ((MDB_cap + cryst_cap * 0.753 ) * 365 * 0.9)
                    + cost_th
                    + cost_elec
                    + cost_NaCl
                    + MDB_var_opex
                )

                capex_lcow = total_capex * 0.1 / ((MDB_cap + cryst_cap * 0.753 ) * 365 * 0.9)
                fixed_opex_lcow = total_fix_opex / ((MDB_cap + cryst_cap * 0.753 ) * 365 * 0.9)
                opex_lcow = (cost_th
                    + cost_elec
                    + cost_NaCl + MDB_var_opex)

                # Update value if a lower LCOW is achieved
                if LCOW < results["LCOW"]:
                    results = {"Curtailment percentage (%)": curt_perc, "LCOW": LCOW, "Solar multiple": sm,  "Storage": storage, 
                                "LCOH from solar field": lcoh, "LCOH adjusted from curtailment":lcoh_after_curtail , "Grid percentage": grid_perc,
                                "MDB SEEC (kWh-e/m3)": SEEC_MDB, "MDB STEC (kWh-th/m3)": STEC_MDB,
                                "Cryst SEEC (kWh-e/m3)": SEEC_cryst,"Cryst STEC (kWh-th/m3)": STEC_cryst,  "Pumping energy (kWh-e/m3)": pump_energy_brackish,
                                "MDB CAPEX ($)": MDB_capex, "Cryst CAPEX ($)": cryst_capex, "Drilling cost ($)": drilling_cost,
                                "Elec cost ($/m3)": cost_elec , "Heat cost ($/m3)": cost_th, "Salt recovery ($/m3)": cost_NaCl
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
            LCOH_solar = results["LCOH from solar field"]
            LCOH_curt = results["LCOH adjusted from curtailment"]
            grid_perc = results["Grid percentage"]
            MDB_SEEC = results["MDB SEEC (kWh-e/m3)"]
            MDB_STEC = results["MDB STEC (kWh-th/m3)"]
            cryst_STEC = results["Cryst STEC (kWh-th/m3)"]
            cryst_SEEC = results["Cryst SEEC (kWh-e/m3)"]
            pumping_energy = results["Pumping energy (kWh-e/m3)"]
            MDB_capex = results["MDB CAPEX ($)"]
            cryst_capex = results["Cryst CAPEX ($)"]
            drilling_cost = results["Drilling cost ($)"]

            (county, state_id, county_id, 
            latitude, longitude, brackish_tds, depth, 
            brackish_avail, desal_cap, 
            cryst_cap, Thermal_req, grid_rate, from_EIA,
            ) = info

            csv_writers[i].writerow([county, state_id, county_id, latitude, longitude, 
                                    brackish_tds/1000, depth, brackish_avail, 
                                    MDB_rr, desal_cap, 
                                    cryst_cap, Thermal_req*sm, sm, storage, 
                                    LCOH_solar, LCOH_curt,
                                    grid_rate, from_EIA, LCOW, 100-grid_perc, grid_perc, curt_perc,
                                    MDB_SEEC, MDB_STEC,
                                    cryst_SEEC, cryst_STEC, pumping_energy, 
                                    MDB_capex, cryst_capex, drilling_cost,
                                    cost_elec, cost_th, cost_NaCl,
                                    capex_lcow, fixed_opex_lcow, opex_lcow
                                    ])


    except:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno)




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
    skipped, failed, half_failed = FPC_MD_cryst_cal()

    print('skipped', len(skipped))
    print('failed', len(failed))
    print('half failed', len(half_failed))

