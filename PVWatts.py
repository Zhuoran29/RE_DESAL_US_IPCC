from PySSC import PySSC
def PVwatts(weatherfile_path, capacity, tilt):
	ssc = PySSC()
	ssc.module_exec_set_print(0)
	data = ssc.data_create()
	ssc.data_set_string( data, b'solar_resource_file', b'' + weatherfile_path.encode("ascii", "backslashreplace") )
	ssc.data_set_number( data, b'system_capacity', capacity )
	ssc.data_set_number( data, b'module_type', 0 )
	ssc.data_set_number( data, b'dc_ac_ratio', 1.1 )
	ssc.data_set_number( data, b'array_type', 2 )
	ssc.data_set_number( data, b'tilt', tilt )
	ssc.data_set_number( data, b'azimuth', 180 )
	ssc.data_set_number( data, b'gcr', 0.40000000000000002 )
	ssc.data_set_number( data, b'losses', 14.075660705566406 )
	ssc.data_set_number( data, b'en_snowloss', 0 )
	ssc.data_set_number( data, b'inv_eff', 96 )
	ssc.data_set_number( data, b'batt_simple_enable', 0 )
	ssc.data_set_number( data, b'adjust:constant', 0 )
	module = ssc.module_create(b'pvwattsv7')	
	ssc.module_exec_set_print( 0 )
	if ssc.module_exec(module, data) == 0:
		print ('pvwattsv7 simulation error')
		idx = 1
		msg = ssc.module_log(module, 0)
		while (msg != None):
			print ('	: ' + msg.decode("utf - 8"))
			msg = ssc.module_log(module, idx)
			idx = idx + 1
		SystemExit( "Simulation Error" );
	ssc.module_free(module)
	ssc.data_set_number( data, b'enable_interconnection_limit', 0 )
	ssc.data_set_number( data, b'grid_interconnection_limit_kwac', 100000 )
	ssc.data_set_array_from_csv( data, b'grid_curtailment', b'/Users/zhuoranzhang/Documents/Github_repos/Modeling_code/grid_curtailment.csv');
	module = ssc.module_create(b'grid')	
	ssc.module_exec_set_print( 0 )
	if ssc.module_exec(module, data) == 0:
		print ('grid simulation error')
		idx = 1
		msg = ssc.module_log(module, 0)
		while (msg != None):
			print ('	: ' + msg.decode("utf - 8"))
			msg = ssc.module_log(module, idx)
			idx = idx + 1
		SystemExit( "Simulation Error" )
	ssc.module_free(module)
	ssc.data_set_number( data, b'capital_cost', 825*capacity )
	ssc.data_set_number( data, b'fixed_operating_cost', 14*capacity )
	ssc.data_set_number( data, b'variable_operating_cost', 0 )
	ssc.data_set_number( data, b'fixed_charge_rate', 0.1 )
	module = ssc.module_create(b'lcoefcr')	
	ssc.module_exec_set_print( 0 )
	if ssc.module_exec(module, data) == 0:
		print ('lcoefcr simulation error')
		idx = 1
		msg = ssc.module_log(module, 0)
		while (msg != None):
			print ('	: ' + msg.decode("utf - 8"))
			msg = ssc.module_log(module, idx)
			idx = idx + 1
		SystemExit( "Simulation Error" )
	ssc.module_free(module)
	energy = ssc.data_get_array(data, b'gen')
	capacity_factor = ssc.data_get_number(data, b'capacity_factor')
	kwh_per_kw = ssc.data_get_number(data, b'kwh_per_kw')
	lcoe_fcr = ssc.data_get_number(data, b'lcoe_fcr')
	monthly_energy = ssc.data_get_array(data, b'monthly_energy')
	ssc.data_free(data)

	return [energy, capacity_factor, kwh_per_kw, lcoe_fcr, monthly_energy]

if __name__ == "__main__":
	res = (PVwatts(weatherfile_path = "/Users/zhuoranzhang/Documents/DOE_Chile/DOE_CSP_PROJECT/SAM_flatJSON/solar_resource/USA FL Miami Intl Ap (TMY3).csv",
			capacity = 15,
			tilt = 33))
	for i in res:
		print(i)