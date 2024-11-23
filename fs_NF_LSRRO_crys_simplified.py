from pyomo.environ import (
    Block,
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    value,
    Var,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from watertap.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core import UnitModelCostingBlock
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_unfixed_variables,
)

from watertap.unit_models.zero_order import NanofiltrationZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.costing.zero_order_costing import ZeroOrderCosting

from parameter_sweep import LinearSample, parameter_sweep
from watertap.flowsheets.lsrro import lsrro
from watertap.flowsheets.lsrro.lsrro import ABTradeoff, ACase, BCase
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    constraint_scaling_transform,
    unscaled_variables_generator,
    unscaled_constraints_generator,
    badly_scaled_var_generator,
)
import idaes.logger as idaeslog
from watertap_reflo.src.watertap_contrib.reflo.analysis.example_flowsheets.multi_effect_NaCl_crystallizer import (
    build_fs_multi_effect_crystallizer,
    add_costings,
    multi_effect_crystallizer_initialization,
    get_model_performance,
)
import idaes.core.util.scaling as iscale
from pyomo.environ import (
    units as pyunits,
    Expression,
)
from watertap.unit_models.reverse_osmosis_1D import (
    ReverseOsmosis1D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from parameter_sweep import LinearSample, parameter_sweep
from watertap.flowsheets.lsrro import lsrro
from watertap.flowsheets.lsrro.lsrro import build, set_operating_conditions, optimize_set_up, ABTradeoff, ACase, BCase, run_lsrro_case, initialize, solve

solver = get_solver()

def fs_NF_LSRRO_crys_simple(capacity, # m3/day
                     solute_dict,
                     tds):  # mg/L
    m = ConcreteModel()
    m.db = Database()

    if sum(solute_dict.values()) == 0:
        solute_dict["Na"] = tds / 58.5 * 23
        solute_dict["Cl"] = tds / 58.5 * 35.5

    # Avoid 0 for zero_order_base_properties
    for key,v in solute_dict.items():
        if v == 0:
            solute_dict[key] = 1e-8

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.params = WaterParameterBlock(solute_list=[
                                                    "magnesium",
                                                    "potassium",
                                                    "sodium",
                                                    "sulfate",
                                                    "carbonate",
                                                    "bicarbonate",
                                                    "chloride",
                                                    "calcium",
                                                    ])

    m.fs.NF = NanofiltrationZO(property_package=m.fs.params, database=m.db)
    feed_flow_rate = 10

    m.fs.NF.inlet.flow_mass_comp[0, "H2O"].fix(feed_flow_rate * 1000)
    m.fs.NF.inlet.flow_mass_comp[0, "sulfate"].fix(solute_dict["SO4"] * feed_flow_rate / 1000)
    m.fs.NF.inlet.flow_mass_comp[0, "sodium"].fix(solute_dict["Na"] * feed_flow_rate / 1000)
    m.fs.NF.inlet.flow_mass_comp[0, "potassium"].fix(solute_dict["K"] * feed_flow_rate / 1000)
    m.fs.NF.inlet.flow_mass_comp[0, "magnesium"].fix(solute_dict["Mg"] * feed_flow_rate / 1000)
    m.fs.NF.inlet.flow_mass_comp[0, "carbonate"].fix(solute_dict["CO3"] * feed_flow_rate / 1000)
    m.fs.NF.inlet.flow_mass_comp[0, "bicarbonate"].fix(solute_dict["HCO3"] * feed_flow_rate / 1000)
    m.fs.NF.inlet.flow_mass_comp[0, "chloride"].fix(solute_dict["Cl"] * feed_flow_rate / 1000)
    m.fs.NF.inlet.flow_mass_comp[0, "calcium"].fix(solute_dict["Ca"] * feed_flow_rate / 1000) 

    # Load parameters
    data = m.db.get_unit_operation_parameters("nanofiltration")
    m.fs.NF.load_parameters_from_database()
    m.fs.NF.initialize(outlvl = idaeslog.CRITICAL)

    # electroneutrality
    charge_comp = {"sodium": 1, "calcium": 2, "magnesium": 2, "sulfate": -2, "chloride": -1,
                    "potassium": 1, "carbonate": -2, "bicarbonate": -1,}
    molar_mass = {"sodium": 23, "calcium": 40, "magnesium": 24, "sulfate": 96, "chloride": 35.5,
                    "potassium": 39, "carbonate": 60, "bicarbonate": 61,}

    m.fs.NF.eq_electroneutrality = Constraint(
        expr=0
        == sum(
            charge_comp[j]
            * m.fs.NF.treated.flow_mass_comp[0, j] * 1000 / molar_mass[j]
            for j in charge_comp.keys()
        )
    )   

    constraint_scaling_transform(m.fs.NF.eq_electroneutrality, 1e5)

    # Create TDS item for permeate and brine
    m.fs.NF_treated_tds = Expression(
        expr=(sum(m.fs.NF.treated.flow_mass_comp[0, j] for j in charge_comp) 
        / m.fs.NF.properties_treated[0].flow_vol), doc="TDS of NF permeate"
    )
    m.fs.NF_brine_tds = Expression(
        expr=(sum(m.fs.NF.byproduct.flow_mass_comp[0, j] for j in charge_comp) 
        / m.fs.NF.properties_byproduct[0].flow_vol), doc="TDS of NF brine"
    )

    try:
        m.fs.NF.removal_frac_mass_comp[0, "chloride"].unfix()

        results = solver.solve(m)
        assert check_optimal_termination(results)  

    except:
        try:
            m.fs.NF.removal_frac_mass_comp[0, "chloride"].fix(0.15)
            m.fs.NF.removal_frac_mass_comp[0, "sodium"].unfix()

            results = solver.solve(m)
            assert check_optimal_termination(results) 
        except:
            try:
                m.fs.NF.removal_frac_mass_comp[0, "sodium"].fix(0.1)
                m.fs.NF.removal_frac_mass_comp[0, "calcium"].unfix()
                results = solver.solve(m)
                assert check_optimal_termination(results) 

            except:
                try:
                    m.fs.NF.removal_frac_mass_comp[0, "calcium"].fix(0.79)
                    m.fs.NF.removal_frac_mass_comp[0, "sulfate"].unfix()
                    results = solver.solve(m)
                    assert check_optimal_termination(results) 
                except:
                    solute_dict["Na"] = tds / 58.5 * 23
                    solute_dict["Cl"] = tds / 58.5 * 35.5
                    m.fs.NF.inlet.flow_mass_comp[0, "sulfate"].fix(1e-8)
                    m.fs.NF.inlet.flow_mass_comp[0, "sodium"].fix(solute_dict["Na"] * feed_flow_rate / 1000)
                    m.fs.NF.inlet.flow_mass_comp[0, "potassium"].fix(1e-8)
                    m.fs.NF.inlet.flow_mass_comp[0, "magnesium"].fix(1e-8)
                    m.fs.NF.inlet.flow_mass_comp[0, "carbonate"].fix(1e-8)
                    m.fs.NF.inlet.flow_mass_comp[0, "bicarbonate"].fix(1e-8)
                    m.fs.NF.inlet.flow_mass_comp[0, "chloride"].fix(solute_dict["Cl"] * feed_flow_rate / 1000)
                    m.fs.NF.inlet.flow_mass_comp[0, "calcium"].fix(1e-8)


                    m.fs.NF.removal_frac_mass_comp[0, "sulfate"].fix()
                    m.fs.NF.removal_frac_mass_comp[0, "chloride"].unfix()

                    results = solver.solve(m)
                    assert check_optimal_termination(results)                     

    # Get TDS (kg/m3) and flow rate (m3/s) of NF permeate
    NF_permeate_tds = value(pyunits.convert(m.fs.NF_treated_tds, to_units = pyunits.kg/ pyunits.m**3))
    NF_brine_tds = value(pyunits.convert(m.fs.NF_brine_tds, to_units = pyunits.kg/ pyunits.m**3))
    NF_permeate_vol = value(m.fs.NF.properties_treated[0].flow_vol)


    # Add cost
    m.fs.zo_costing = ZeroOrderCosting()
    m.fs.NF.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.zo_costing)
    
    m.fs.zo_costing.cost_process()
    # m.fs.zo_costing.initialize()

    # m.fs.zo_costing.electricity_cost.fix(0.07)
    @m.Expression()
    def LCOW(b):
        return (
            b.fs.zo_costing.total_capital_cost * b.fs.zo_costing.capital_recovery_factor
            + m.fs.zo_costing.total_operating_cost
        ) / (
            pyunits.convert(
                b.fs.NF.properties_treated[0].flow_vol,
                to_units=pyunits.m**3 / pyunits.year,
            )
            * 0.9 # b.fs.zo_costing.utilization_factor
        )
    
    results = solver.solve(m)
    try:    
        assert check_optimal_termination(results) 
    except:
        NF_permeate_tds = tds / 1000

    # Calculate capacity of each components
    cbrine = 180
    lsrro_rr = 1 - NF_permeate_tds / cbrine
    
    lsrro_cap = capacity / (0.247 + 0.753 / lsrro_rr) # m3/day
    cryst_cap = lsrro_cap * (1 - lsrro_rr) / lsrro_rr # m3/day

    nf_cap = lsrro_cap / lsrro_rr
    feed_vol = nf_cap / 0.85

    capacities = [feed_vol, nf_cap, lsrro_cap, cryst_cap, NF_permeate_tds, lsrro_rr ]

        
    # Add lsrro
    # lsrro, results = run_lsrro_case(
    #     number_of_stages=3,
    #     # water_recovery=0.9861,
    #     Cin=NF_permeate_tds,  # inlet NaCl conc kg/m3,
    #     Qin=1e-3,  # inlet feed flowrate m3/s
    #     Cbrine=180,  # brine conc kg/m3
    #     A_case=ACase.optimize,
    #     B_case=BCase.optimize,
    #     AB_tradeoff=ABTradeoff.equality_constraint,
    #     # A_value=4.2e-12, #membrane water permeability coeff m/s-Pa
    #     has_NaCl_solubility_limit=True,
    #     has_calculated_concentration_polarization=True,
    #     has_calculated_ro_pressure_drop=True,
    #     permeate_quality_limit=500e-6,
    #     AB_gamma_factor=1,
    #     B_max=3.5e-6,
    #     number_of_RO_finite_elements=10,
    #     set_default_bounds_on_module_dimensions=True,
    # )
    lsrro = None

    # Add crystallizers
    cryst = build_fs_multi_effect_crystallizer(
        operating_pressure_eff1=0.4455,  # bar
        operating_pressure_eff2=0.2758,  # bar
        operating_pressure_eff3=0.1651,  # bar
        operating_pressure_eff4=0.095,  # bar
        feed_flow_mass= 1,  # kg/s
        feed_mass_frac_NaCl=0.1613,
        feed_pressure=101325,  # Pa
        feed_temperature=273.15 + 25,  # K
        crystallizer_yield=0.8,
        steam_pressure=1.5,  # bar (gauge pressure)
    )
    add_costings(cryst)

    # Negative value for salt recovery value ($/kg)
    cryst.fs.costing.crystallizer.steam_cost.fix(0)
    cryst.fs.costing.crystallizer.NaCl_recovery_value.fix(-0.024)

    multi_effect_crystallizer_initialization(cryst)

    results = solver.solve(cryst)

    eff_1 = cryst.fs.eff_1

    eff_1.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].unfix()
    eff_1.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].unfix()
    eff_1.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"].fix()


    @cryst.Constraint(eff_1.flowsheet().time, doc="total of feed flow")
    def total_feed(b, t):
        effs = [cryst.fs.eff_1, cryst.fs.eff_2, cryst.fs.eff_3, cryst.fs.eff_4]
        return (
            sum(i.properties_in[0].flow_vol_phase["Liq"] for i in effs)
            == cryst_cap / 86400 * pyunits.m**3 / pyunits.s
        )
    try:
        results = solver.solve(cryst)
    except:
    # Update scaling factors
        for eff in [cryst.fs.eff_1, cryst.fs.eff_2, cryst.fs.eff_3, cryst.fs.eff_4]:

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

        results = solver.solve(cryst)
        # half_failed.append(county_id)


    data_table, overall_performance = get_model_performance(cryst)    


    return m, lsrro, cryst, data_table, overall_performance, NF_brine_tds, capacities


if __name__ == "__main__":
    capacity = 100
    feed_flow_rate = 10
    solute_dict = {
        'Na': 553.9776611,
        'Mg': 8.724057656,
        'Ca': 10.46886919,
        'K': 15.70330378,
        'Cl': 678.0746445,
        'CO3': 0,
        'HCO3': 396.0722176,
        'SO4': 6.979246125,
    }

    m, lsrro, cryst, data_table, overall_performance, NF_permeate_tds, capacities = fs_NF_LSRRO_crys_simple(capacity,
                                                                        solute_dict,
                                                                        3)

    # results = solver.solve(m)
    # # Check for optimal solution
    # assert check_optimal_termination(results)

    df = m.fs.NF._get_stream_table_contents()
    p = m.fs.NF._get_performance_contents()


    charge_comp = {"sodium": 1, "calcium": 2, "magnesium": 2, "sulfate": -2, "chloride": -1,
                    "potassium": 1, "carbonate": -2, "bicarbonate": -1,}
    molar_mass = {"sodium": 23, "calcium": 40, "magnesium": 24, "sulfate": 96, "chloride": 35.5,
                    "potassium": 39, "carbonate": 60, "bicarbonate": 61,}
