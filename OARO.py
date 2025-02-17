# Need to pip install APMonitor and add dependency in beta version
from APMonitor.apm import *
import os
from pathlib import Path


class OARO(object):
    def __init__(self,
                 FeedC_r = 35, # Feed concentration (g/L)
                 Capacity = 1000, # System capacity  (m3/day)
                 rr = 25,
                 ):
        self.FeedC_r = int(FeedC_r)
        self.Capacity = Capacity
        self.rr = rr
        self.base_path = Path(__file__).resolve().parent.absolute() 
    def design(self):
        s = 'http://byu.apmonitor.com'
        
        # Application name
        a = 'DesalinationModels'
        
        # Clear previous application
        apm(s,a,'clear all')

        #APM Model Files
        ROOT_DIR = self.base_path / "EnhancedRO_APM_Models"
        
        
        #LSRRO Cases: 70 g/l Feed, 75% recovery; 35 g/l Feed, 85% recovery; 20 g/l Feed, 92% recovery
        
        # LSRRO_70g_L_75rec=ROOT_DIR+'/LSRRO_8stage_cin_70_r_75.apm'
        # LSRRO_35g_L_85rec=ROOT_DIR+'/LSRRO_6stage_cin_35_r_85.apm'
        # LSRRO_20g_L_92rec=ROOT_DIR+'/LSRRO_6stage_cin_20_r_92.apm'
        
        # #COMRO 
        # COMRO_70g_L_75rec=ROOT_DIR+'/COMRO_4stage_cin_70_r_75.apm' # relatively long to find solution
        # COMRO_35g_L_85rec=ROOT_DIR+'/COMRO_3stage_cin_35_r_85.apm'
        # COMRO_20g_L_92rec=ROOT_DIR+'/COMRO_4stage_cin_20_r_92.apm'
        
        #OARO 
        # recrate = str(self.rr)
        if self.rr == 1:
            if self.FeedC_r == 70:
                recrate = '75'
            elif self.FeedC_r == 35:
                recrate = '85'
            elif self.FeedC_r == 20:
                recrate = '92'
            elif self.FeedC_r == 125:
                recrate = '50'
        
        elif self.rr == 2:
            if self.FeedC_r == 70:
                recrate = '50'
            elif self.FeedC_r == 35:
                recrate = '75'
            elif self.FeedC_r == 20:
                recrate = '85'
            elif self.FeedC_r == 125:
                recrate = '25'        
            
        OARO_70g_L_50rec=ROOT_DIR/'OARO_2stage_cin_70_r_50.apm'     
        OARO_70g_L_75rec=ROOT_DIR/'OARO_4stage_cin_70_r_75.apm'
        OARO_35g_L_85rec=ROOT_DIR/'OARO_2stage_cin_35_r_85.apm' # LCOW result is a little bit higher than result in MATLAB (1.5 vs 1.4). SEC is also a little bit higher (5.3 vs 5.2)
        OARO_20g_L_92rec=ROOT_DIR/'OARO_2stage_cin_20_r_92.apm' # not converging to same sol as MATLAB; closer when changing otol from default of 1e-6 to 1e-3
        OARO_125g_L_50rec=ROOT_DIR/'OARO_5stage_cin_125_r_50.apm'
        OARO_125g_L_25rec=ROOT_DIR/'OARO_3stage_cin_125_r_25.apm'        
        OARO_20g_L_85rec=ROOT_DIR/'OARO_2stage_cin_20_r_85.apm'
        OARO_35g_L_75rec=ROOT_DIR/'OARO_2stage_cin_35_r_75.apm'
        
        feedsal = str(self.FeedC_r)
        tech = 'OARO'
        model_file= eval(tech + "_" + feedsal + "g_L_"+ recrate + "rec")
        
        # csv_load(s,a,csvfilename) 
        
        apm_load(s,a,model_file)
        
        
        #apm_option(s,a,'nlc.diaglevel',10)
        apm_option(s,a,'apm.imode',3) # chooses operation mode; 3= steady-state optimization
        apm_option(s,a,'apm.max_iter',1000)  # max number of iterations
        apm_option(s,a,'apm.scaling',1) # activate scaling
        apm_option(s,a,'apm.reduce',5) # activates reduction
        if (tech=='OARO'): 
            if feedsal=='70': 
                apm_option(s,a,'apm.otol',1e-4) # default tolerance is 1e-6
            elif feedsal=='20':
                apm_option(s,a,'apm.otol',1e-3) # default tolerance is 1e-6
            elif (feedsal=='125') & (recrate=='50'):
                apm_option(s,a,'apm.otol',1e-7) # default tolerance is 1e-6
        else:
            apm_option(s,a,'apm.otol',1e-6) # default tolerance is 1e-6
        output = apm(s,a,'solve');
        
        SEC = 0
        rtq = 0
        if (apm_tag(s,a,'apm.appstatus')==1):
            # retrieve solution if successful
            sol = apm_sol(s,a)
            #z = sol.x 
        
            rtq=sol['rtq']
            LCOW=sol['lcow']
            # oaro_area = sol['Amksum']
            # ro_area = sol['Amksum']
            # print(tech)
            # if tech=='COMRO':
            #     print("Cin=",sol['cpin'])
            # else:
            #     print("Cin=",sol['cin'])
            # print("Water Recovery Rate=",sol['rt'])
            # print("Solution recovery rate=",rtq)
            # print("LCOW=",LCOW,"$/m3")
            # print("area=", sol['amj'],'m2')
            # print("Specialized membrane area =",sol['amk'])
            # print('pumpcost=', sol['pumptotalcostsum[{}]'.format(int(sol['nstage']))])
            # print('erdcost=', sol['erdcostsum[{}]'.format(int(sol['nstage']))])
            # print('opex=', sol['opex'])
            # print('capex=',sol['capex'] )
            # print('cr_factor=', sol['cr_factor'])
            # print('energycost=',sol['energycost'] )
            # print('salmakeupcost=',sol['salmakeupcost'] )
            # print('Qpk=', sol['qpk'])


            #print(sol['mpj[0]'])
           
           
            if sol['nstage']==2:
                if tech=='OARO':
                    SEC=sol['secsum'] # OARO 2stage
                    Cb=sol['cb']      # OARO 2stage
                    Amj=sol['amjsum']
                    Amk=sol['amksum']
                    # SEC= sol['sec']
                elif tech=='LSRRO':
                    SEC= sol['sec']
                    Cb = sol['cfjout']
                    
           
                print("Brine concentration=",Cb,"g/L")
        
          
            elif sol['nstage']>2:
                SEC= sol['sec']
               #print("SEC=",SEC,"kWh/m3")
        
                if tech=='LSRRO':
                    Cbstring=""""sol['cfjout[{}]']".format(int(sol['ncc']))"""
                    Cb=eval(eval(Cbstring))
                    Amj=eval(eval(""""sol['amjsum[{}]']".format(int(sol['ncc']))""")) 
                    Amk=sol['amksum']           
              
                elif tech=='COMRO':
                    Cbstring=""""sol['cfoutj[{}]']".format(int(sol['nstage']))"""   ########## may want to change variable name for final brine concentration to match in all technologies (right now, cfoutj vs cfjout)
                    Cb=eval(eval(Cbstring))
                    Amj=eval(eval(""""sol['amjsum[{}]']".format(int(sol['nstage']))""")) 
                    Amk=eval(eval(""""sol['amksum[{}]']".format(int(sol['nstage']))""")) 
                    
                elif tech=='OARO':
                    Cb=sol['cfout']
                    Amj=eval(eval(""""sol['amjsum[{}]']".format(int(sol['ncc']))""")) 
                    Amk=sol['amksum']
            elif (tech=='COMRO') & (sol['nstage']==1):
                Amj=sol['amj1']
                Amk=sol['amk1']       
                SEC=sol['sec']
                Cb=sol['cfoutj']               

            specialized_area = Amj
            ro_area = Amk           
            
            print("Brine concentration=",Cb,"g/L")       
            print("Cost of electricity=",sol['ecost'],"kWh/m3")
            print("Number of stages=",sol['nstage'])
            print(model_file)            
        else:
            # % not successful, retrieve infeasiblilities report
            apm_get(s,a,'infeasibilities.txt')
            # disp(output)
        self.ThPower = SEC * (468/24*rtq)
        self.num_modules = (self.Capacity / (468/24*rtq) /24 )
        self.rtq = rtq
        
                    
        costout = {'oaro_area': specialized_area * self.num_modules,
                   'ro_area':   ro_area * self.num_modules,
                   'pumpcost':  sol['pumptotalcostsum[{}]'.format(int(sol['nstage']))] * self.num_modules ,
                   'erdcost': sol['erdcostsum[{}]'.format(int(sol['nstage']))] * self.num_modules,
                   # 'ma':   sol['ma'] * self.num_modules,
                   'sec':  SEC
                  }         
        brine_s = self.FeedC_r / (1 - int(recrate)/100)
        self.design_output = []
        # self.design_output.append({'Name':'Number of modules required','Value':self.num_modules,'Unit':''})
        self.design_output.append({'Name':'Permeate flow rate','Value':(468/24*rtq) * self.num_modules  *24,'Unit':'m3/day'})    
        self.design_output.append({'Name':'Electric energy requirement','Value':self.ThPower * self.num_modules / 1000,'Unit':'MW(e)'})
        self.design_output.append({'Name':'Brine concentration','Value':brine_s,'Unit':'g/L'})
        
        self.design_output.append({'Name':'Specialized membrane area','Value': specialized_area * self.num_modules,'Unit':'m2'})    
        self.design_output.append({'Name':'Conventional RO membrane area','Value':ro_area * self.num_modules,'Unit':'m2'})
        self.design_output.append({'Name':'Specific electricity consumption','Value': SEC,'Unit':'kWh(e)/m3'})
        self.design_output.append({'Name':'Recovery ratio','Value': rtq*100 ,'Unit':'%'})    
        
        return self.design_output, costout
        
    def simulation(self, gen, storage):
        self.thermal_load = self.ThPower * self.num_modules # kWh
        self.max_prod = (468/24*self.rtq) * self.num_modules # m3/h

        self.storage_cap = storage * self.thermal_load # kWh
        self.Fossil_f = 1
        
        to_desal = [0 for i in range(len(gen))]
        to_storage =  [0 for i in range(len(gen))]
        storage_load =  [0 for i in range(len(gen))]
        storage_cap_1 =  [0 for i in range(len(gen))]
        storage_cap_2 = [0 for i in range(len(gen))]
        storage_status =  [0 for i in range(len(gen))]
        solar_loss =  [0 for i in range(len(gen))]
        load =  [0 for i in range(len(gen))]
        prod =  [0 for i in range(len(gen))]
        fuel =  [0 for i in range(len(gen))]
        energy_consumption =  [0 for i in range(len(gen))]
        actual_load =  [0 for i in range(len(gen))]
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
            actual_load[i] = max(0,load[i])
            energy_consumption[i] = fuel[i]+load[i] 
            prod[i] = (fuel[i]+load[i] )/ self.thermal_load * self.max_prod
        
        Month = [0,31,59,90,120,151,181,212,243,273,304,334,365]
        Monthly_prod = [ sum( prod[Month[i]*24:(Month[i+1]*24)] ) for i in range(12) ]
        monthly_grid_consumed = [ sum( fuel[Month[i]*24:(Month[i+1]*24)] ) / sum( energy_consumption[Month[i]*24:(Month[i+1]*24)] )  for i in range(12)]
        monthly_elec_curtailed = [ (sum( gen[Month[i]*24:(Month[i+1]*24)]) - sum(load[Month[i]*24:(Month[i+1]*24)] ))/sum( gen[Month[i]*24:(Month[i+1]*24)] ) for i in range(12)]


        simu_output = []

        simu_output.append({'Name':'Water production','Value':prod,'Unit':'m3'})
        simu_output.append({'Name':'Storage status','Value':storage_status,'Unit':'kWh'})
        simu_output.append({'Name':'Storage Capacity','Value':self.storage_cap,'Unit':'kWh'})
        simu_output.append({'Name':'Grid electricity usage','Value':fuel,'Unit':'kWh'})
        simu_output.append({'Name':'Total water production','Value':sum(prod),'Unit':'m3'})
        simu_output.append({'Name':'Monthly water production','Value': Monthly_prod,'Unit':'m3'})
        simu_output.append({'Name':'Total fossil fuel usage','Value':sum(fuel),'Unit':'kWh'})
        simu_output.append({'Name':'Percentage of fossil fuel consumption','Value':sum(fuel)/max(1,sum(energy_consumption))*100,'Unit':'%'})           
        simu_output.append({'Name':'Curtailed solar electric energy','Value':max(0, (sum(gen) - sum(load)) / 1000000) ,'Unit':'GWh'})   
        simu_output.append({'Name':'Percentage of curtailed energy','Value':max(0, (sum(gen) - sum(load))) / sum(gen) * 100 ,'Unit':'%'})    
        # Add brine volume and concentration (using 100% rejection(make it a variable))
        
        return simu_output, monthly_grid_consumed, monthly_elec_curtailed
                    
        
    
        
