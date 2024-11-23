from PySSC import PySSC

def IPH_LF(weatherfile_path = '/Users/zhuoranzhang/Documents/DOE_Chile/DOE_CSP_PROJECT/SAM_flatJSON/solar_resource/USA AZ Tucson International Ap (TMY3).csv', 
		   capacity = 200000,
		   specified_solar_multiple = 1,
		   ):
	
	q_pb_des = capacity / 1000
	I_bn_des = 950

	if capacity < 10000:
		nModBoil = int(capacity / 2000) + 1
	else:
		nModBoil = 6

	nLoops, solarm = dependency(q_pb_des,
								I_bn_des,
								nModBoil,
								specified_solar_multiple,
								)
	# print('nloops:', nLoops)
	# print('actual sm', solarm)

	ssc = PySSC()
	ssc.module_exec_set_print(0)
	data = ssc.data_create()
	ssc.data_set_string( data, b'file_name', b'' + weatherfile_path.encode("ascii", "backslashreplace") )
	ssc.data_set_number( data, b'I_bn_des', I_bn_des )
	ssc.data_set_number( data, b'T_cold_ref', 100 )
	ssc.data_set_number( data, b'P_turb_des', 20 )
	ssc.data_set_number( data, b'x_b_des', 0.75 )
	ssc.data_set_number( data, b'q_pb_des', q_pb_des )
	ssc.data_set_number( data, b'use_quality_or_subcooled', 0 )
	ssc.data_set_number( data, b'deltaT_subcooled', 10 )
	ssc.data_set_number( data, b'fP_hdr_c', 0.01 )
	ssc.data_set_number( data, b'fP_sf_boil', 0.074999999999999997 )
	ssc.data_set_number( data, b'fP_hdr_h', 0.025000000000000001 )
	ssc.data_set_number( data, b'nModBoil', nModBoil )
	ssc.data_set_number( data, b'nLoops', nLoops )
	ssc.data_set_number( data, b'eta_pump', 0.84999999999999998 )
	ssc.data_set_number( data, b'theta_stow', 10 )
	ssc.data_set_number( data, b'theta_dep', 10 )
	ssc.data_set_number( data, b'T_fp', 10 )
	ssc.data_set_number( data, b'Pipe_hl_coef', 0.0035000000000000001 )
	ssc.data_set_number( data, b'SCA_drives_elec', 0.20000000000000001 )
	ssc.data_set_number( data, b'ColAz', 0 )
	ssc.data_set_number( data, b'e_startup', 2.7000000000000002 )
	ssc.data_set_number( data, b'T_amb_des_sf', 42 )
	ssc.data_set_number( data, b'V_wind_max', 20 )
	ssc.data_set_number( data, b'csp.lf.sf.water_per_wash', 0.02 )
	ssc.data_set_number( data, b'csp.lf.sf.washes_per_year', 12 )
	A_aperture = [[ 513.60000000000002 ], [ 0 ]]
	ssc.data_set_matrix( data, b'A_aperture', A_aperture )
	L_col = [[ 44.799999999999997 ], [ 0 ]]
	ssc.data_set_matrix( data, b'L_col', L_col )
	OptCharType = [[ 1 ], [ 0 ]]
	ssc.data_set_matrix( data, b'OptCharType', OptCharType )
	IAM_T = [[ 0.98960000000000004,   0.043999999999999997,   -0.072099999999999997,   -0.23269999999999999,   0 ], [ 0,   0,   0,   0,   0 ]]
	ssc.data_set_matrix( data, b'IAM_T', IAM_T )
	IAM_L = [[ 1.0031000000000001,   -0.22589999999999999,   0.53680000000000005,   -1.6434,   0.72219999999999995 ], [ 0,   0,   0,   0,   0 ]]
	ssc.data_set_matrix( data, b'IAM_L', IAM_L )
	TrackingError = [[ 1 ], [ 0 ]]
	ssc.data_set_matrix( data, b'TrackingError', TrackingError )
	GeomEffects = [[ 0.83999999999999997 ], [ 0 ]]
	ssc.data_set_matrix( data, b'GeomEffects', GeomEffects )
	rho_mirror_clean = [[ 0.93500000000000005 ], [ 0 ]]
	ssc.data_set_matrix( data, b'rho_mirror_clean', rho_mirror_clean )
	dirt_mirror = [[ 0.94999999999999996 ], [ 0 ]]
	ssc.data_set_matrix( data, b'dirt_mirror', dirt_mirror )
	error = [[ 1 ], [ 0 ]]
	ssc.data_set_matrix( data, b'error', error )
	HLCharType = [[ 1 ], [ 0 ]]
	ssc.data_set_matrix( data, b'HLCharType', HLCharType )
	HL_dT = [[ 0,   0.67200000000000004,   0.0025560000000000001,   0,   0 ], [ 0,   0,   0,   0,   0 ]]
	ssc.data_set_matrix( data, b'HL_dT', HL_dT )
	HL_W = [[ 1,   0,   0,   0,   0 ], [ 0,   0,   0,   0,   0 ]]
	ssc.data_set_matrix( data, b'HL_W', HL_W )
	D_2 = [[ 0.066000000000000003 ], [ 0 ]]
	ssc.data_set_matrix( data, b'D_2', D_2 )
	D_3 = [[ 0.070000000000000007 ], [ 0 ]]
	ssc.data_set_matrix( data, b'D_3', D_3 )
	D_4 = [[ 0.115 ], [ 0 ]]
	ssc.data_set_matrix( data, b'D_4', D_4 )
	D_5 = [[ 0.12 ], [ 0 ]]
	ssc.data_set_matrix( data, b'D_5', D_5 )
	D_p = [[ 0 ], [ 0 ]]
	ssc.data_set_matrix( data, b'D_p', D_p )
	Rough = [[ 4.5000000000000003e-05 ], [ 0 ]]
	ssc.data_set_matrix( data, b'Rough', Rough )
	Flow_type = [[ 1 ], [ 0 ]]
	ssc.data_set_matrix( data, b'Flow_type', Flow_type )
	AbsorberMaterial = [[ 1 ], [ 0 ]]
	ssc.data_set_matrix( data, b'AbsorberMaterial', AbsorberMaterial )
	HCE_FieldFrac = [[ 0.98499999999999999,   0.01,   0.0050000000000000001,   0 ], [ 0,   0,   0,   0 ]]
	ssc.data_set_matrix( data, b'HCE_FieldFrac', HCE_FieldFrac )
	alpha_abs = [[ 0.95999999999999996,   0.95999999999999996,   0.80000000000000004,   0 ], [ 0,   0,   0,   0 ]]
	ssc.data_set_matrix( data, b'alpha_abs', alpha_abs )
	b_eps_HCE1 = [[ 0 ], [ 0.1384 ]]
	ssc.data_set_matrix( data, b'b_eps_HCE1', b_eps_HCE1 )
	b_eps_HCE2 = [[ 0 ], [ 0.65000000000000002 ]]
	ssc.data_set_matrix( data, b'b_eps_HCE2', b_eps_HCE2 )
	b_eps_HCE3 = [[ 0 ], [ 0.65000000000000002 ]]
	ssc.data_set_matrix( data, b'b_eps_HCE3', b_eps_HCE3 )
	b_eps_HCE4 = [[ 0 ], [ 0.1384 ]]
	ssc.data_set_matrix( data, b'b_eps_HCE4', b_eps_HCE4 )
	sh_eps_HCE1 = [[ 0 ], [ 0 ]]
	ssc.data_set_matrix( data, b'sh_eps_HCE1', sh_eps_HCE1 )
	sh_eps_HCE2 = [[ 0 ], [ 0 ]]
	ssc.data_set_matrix( data, b'sh_eps_HCE2', sh_eps_HCE2 )
	sh_eps_HCE3 = [[ 0 ], [ 0 ]]
	ssc.data_set_matrix( data, b'sh_eps_HCE3', sh_eps_HCE3 )
	sh_eps_HCE4 = [[ 0 ], [ 0 ]]
	ssc.data_set_matrix( data, b'sh_eps_HCE4', sh_eps_HCE4 )
	alpha_env = [[ 0.02,   0.02,   0,   0 ], [ 0,   0,   0,   0 ]]
	ssc.data_set_matrix( data, b'alpha_env', alpha_env )
	EPSILON_4 = [[ 0.85999999999999999,   0.85999999999999999,   1,   0 ], [ 0,   0,   0,   0 ]]
	ssc.data_set_matrix( data, b'EPSILON_4', EPSILON_4 )
	Tau_envelope = [[ 0.96299999999999997,   0.96299999999999997,   1,   0 ], [ 0,   0,   0,   0 ]]
	ssc.data_set_matrix( data, b'Tau_envelope', Tau_envelope )
	GlazingIntactIn = [[ 1,   1,   0,   1 ], [ 0,   0,   0,   0 ]]
	ssc.data_set_matrix( data, b'GlazingIntactIn', GlazingIntactIn )
	AnnulusGas = [[ 1,   1,   1,   1 ], [ 0,   0,   0,   0 ]]
	ssc.data_set_matrix( data, b'AnnulusGas', AnnulusGas )
	P_a = [[ 0.0001,   750,   750,   0 ], [ 0,   0,   0,   0 ]]
	ssc.data_set_matrix( data, b'P_a', P_a )
	Design_loss = [[ 150,   1100,   1500,   0 ], [ 0,   0,   0,   0 ]]
	ssc.data_set_matrix( data, b'Design_loss', Design_loss )
	Shadowing = [[ 0.95999999999999996,   0.95999999999999996,   0.95999999999999996,   0 ], [ 0,   0,   0,   0 ]]
	ssc.data_set_matrix( data, b'Shadowing', Shadowing )
	Dirt_HCE = [[ 0.97999999999999998,   0.97999999999999998,   1,   0 ], [ 0,   0,   0,   0 ]]
	ssc.data_set_matrix( data, b'Dirt_HCE', Dirt_HCE )
	b_OpticalTable = [[ -180,   -160,   -140,   -120,   -100,   -80,   -60,   -40,   -20,   0,   20,   40,   60,   80,   100,   120,   140,   160,   180,   -999.89999999999998 ], [ 0,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1 ], [ 10,   0.97999999999999998,   0.97444500000000001,   0.97197599999999995,   0.97284700000000002,   0.97690999999999995,   0.97690999999999995,   0.97284700000000002,   0.97197599999999995,   0.97444500000000001,   0.97999999999999998,   0.97444500000000001,   0.97197599999999995,   0.97284700000000002,   0.97690999999999995,   0.97690999999999995,   0.97284700000000002,   0.97197599999999995,   0.97444500000000001,   0.97999999999999998 ], [ 20,   0.93000000000000005,   0.92297600000000002,   0.92893000000000003,   0.94600499999999998,   0.95401899999999995,   0.95401899999999995,   0.94600499999999998,   0.92893000000000003,   0.92297600000000002,   0.93000000000000005,   0.92297600000000002,   0.92893000000000003,   0.94600499999999998,   0.95401899999999995,   0.95401899999999995,   0.94600499999999998,   0.92893000000000003,   0.92297600000000002,   0.93000000000000005 ], [ 30,   0.83999999999999997,   0.83861799999999997,   0.87069099999999999,   0.91302099999999997,   0.94091100000000005,   0.94091100000000005,   0.91302099999999997,   0.87069099999999999,   0.83861799999999997,   0.83999999999999997,   0.83861799999999997,   0.87069099999999999,   0.91302099999999997,   0.94091100000000005,   0.94091100000000005,   0.91302099999999997,   0.87069099999999999,   0.83861799999999997,   0.83999999999999997 ], [ 40,   0.71999999999999997,   0.72994700000000001,   0.80368700000000004,   0.86696099999999998,   0.90003900000000003,   0.90003900000000003,   0.86696099999999998,   0.80368700000000004,   0.72994700000000001,   0.71999999999999997,   0.72994700000000001,   0.80368700000000004,   0.86696099999999998,   0.90003900000000003,   0.90003900000000003,   0.86696099999999998,   0.80368700000000004,   0.72994700000000001,   0.71999999999999997 ], [ 50,   0.55000000000000004,   0.59125499999999998,   0.70745400000000003,   0.79350900000000002,   0.83955999999999997,   0.83955999999999997,   0.79350900000000002,   0.70745400000000003,   0.59125499999999998,   0.55000000000000004,   0.59125499999999998,   0.70745400000000003,   0.79350900000000002,   0.83955999999999997,   0.83955999999999997,   0.79350900000000002,   0.70745400000000003,   0.59125499999999998,   0.55000000000000004 ], [ 60,   0.34000000000000002,   0.43217800000000001,   0.59747799999999995,   0.66400599999999999,   0.69351099999999999,   0.69351099999999999,   0.66400599999999999,   0.59747799999999995,   0.43217800000000001,   0.34000000000000002,   0.43217800000000001,   0.59747799999999995,   0.66400599999999999,   0.69351099999999999,   0.69351099999999999,   0.66400599999999999,   0.59747799999999995,   0.43217800000000001,   0.34000000000000002 ], [ 70,   0.13,   0.26525399999999999,   0.42558600000000002,   0.46449600000000002,   0.47710599999999997,   0.47710599999999997,   0.46449600000000002,   0.42558600000000002,   0.26525399999999999,   0.13,   0.26525399999999999,   0.42558600000000002,   0.46449600000000002,   0.47710599999999997,   0.47710599999999997,   0.46449600000000002,   0.42558600000000002,   0.26525399999999999,   0.13 ], [ 80,   0.01,   0.113694,   0.20891000000000001,   0.23325499999999999,   0.23882800000000001,   0.23882800000000001,   0.23325499999999999,   0.20891000000000001,   0.113694,   0.01,   0.113694,   0.20891000000000001,   0.23325499999999999,   0.23882800000000001,   0.23882800000000001,   0.23325499999999999,   0.20891000000000001,   0.113694,   0.01 ], [ 90,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ]]
	ssc.data_set_matrix( data, b'b_OpticalTable', b_OpticalTable )
	sh_OpticalTable = [[ 0,   0 ], [ 0,   0 ]]
	ssc.data_set_matrix( data, b'sh_OpticalTable', sh_OpticalTable )
	ssc.data_set_number( data, b'heat_sink_dP_frac', 0.01 )
	ssc.data_set_number( data, b'adjust:constant', 4 )
	ssc.data_set_number( data, b'electricity_rate', 0.059999999999999998 )
	ssc.data_set_number( data, b'fixed_operating_cost', 5.6 * capacity * solarm  )
	ssc.data_set_number( data, b'capital_cost', 392 * capacity * solarm )
	ssc.data_set_number( data, b'variable_operating_cost', 0.001 )
	ssc.data_set_number( data, b'fixed_charge_rate', 0.1071793755750878 )
	module = ssc.module_create(b'linear_fresnel_dsg_iph')	
	ssc.module_exec_set_print( 0 )
	if ssc.module_exec(module, data) == 0:
		print ('linear_fresnel_dsg_iph simulation error')
		idx = 1
		msg = ssc.module_log(module, 0)
		while (msg != None):
			print ('	: ' + msg.decode("utf - 8"))
			msg = ssc.module_log(module, idx)
			idx = idx + 1
		SystemExit( "Simulation Error" )
	ssc.module_free(module)
	module = ssc.module_create(b'iph_to_lcoefcr')	
	ssc.module_exec_set_print( 0 )
	if ssc.module_exec(module, data) == 0:
		print ('iph_to_lcoefcr simulation error')
		idx = 1
		msg = ssc.module_log(module, 0)
		while (msg != None):
			print ('	: ' + msg.decode("utf - 8"))
			msg = ssc.module_log(module, idx)
			idx = idx + 1
		SystemExit( "Simulation Error" )
	ssc.module_free(module)
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
	lcoe_fcr = ssc.data_get_number(data, b'lcoe_fcr')
	annual_energy = ssc.data_get_number(data, b'annual_energy')
	ssc.data_free(data)
	return [energy, capacity_factor, lcoe_fcr, annual_energy, solarm]


def dependency(
			   q_pb_des,
			   I_bn_des,
			   nModBoil,
			   specified_solar_multiple,

			   ):
	
	loop_aperture = nModBoil * 513.6

	loop_opt_eff = 0.746
	loop_therm_eff = 0.990

	total_loop_conv_eff = loop_opt_eff * loop_therm_eff

	sm1_aperture = q_pb_des / (I_bn_des * total_loop_conv_eff) * 1e6
	sm1_numloops = sm1_aperture / loop_aperture

	nLoops = int(specified_solar_multiple * sm1_aperture / loop_aperture) + 1

	actual_aper = loop_aperture * nLoops
	solarm = actual_aper / sm1_aperture

	return nLoops, solarm



if __name__ ==  '__main__':
	energy, capacity_factor, lcoe_fcr,annual_energy, solarm = IPH_LF(capacity = 200000,
																  specified_solar_multiple=1.8)
	print('LCOH: ', lcoe_fcr)
	print('CF: ',capacity_factor)
	print('actual sm', solarm)
	print('Annual energy: ',annual_energy)