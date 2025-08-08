# %% Import
import time
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from copy import deepcopy

from scipy.integrate import solve_ivp
from scipy.optimize import shgo,dual_annealing,minimize
from scipy.optimize import approx_fprime

import json
from function_simulation import function_simulation

# %% Simulator
def simulate(time_initial,time_final,DTWIN_state,DTWIN_design,DTWIN_config):
    NEW_DTWIN_state=deepcopy(DTWIN_state)
    
    nn=0   
    for i1 in DTWIN_config['Brxtor_list']:
        ts0=np.array([time_initial,time_final])
        Xo0=np.array([])
        for i2 in DTWIN_config['Species_list']:
            Xo0=np.append(Xo0,DTWIN_state[i1]['Current'][i2])
        u0=np.array([DTWIN_design[i1]['Glucose_feed'],nn,DTWIN_config['number_br'],DTWIN_design[i1]['Induction_time'],DTWIN_design[i1]['Inductor_conc'],
                     DTWIN_design[i1]['Dextrine_feed'],DTWIN_design[i1]['Enzyme_feed']])
        # THs=np.array(DTWIN_config['Params'])
        THs=DTWIN_config['Params']
        # print(THs)
        D0=DTWIN_design[i1]['Pulses'].copy()
        
        t,y=function_simulation(ts0,Xo0,u0,THs,D0)
        
        nn2=0
        for i2 in DTWIN_config['Species_list']:
            NEW_DTWIN_state[i1]['All'][i2]['time']=np.append(np.array(NEW_DTWIN_state[i1]['All'][i2]['time']),t[1:].flatten()).tolist()
            NEW_DTWIN_state[i1]['All'][i2]['Value']=np.append(np.array(NEW_DTWIN_state[i1]['All'][i2]['Value']),y[1:,nn2].flatten()).tolist()
            
            NEW_DTWIN_state[i1]['Current'][i2]=float(y[-1,nn2])

            nn2=nn2+1
 
        nn=nn+1
    return NEW_DTWIN_state
# %% Sampler
def sample(time_initial,time_final,DTWIN_state,DTWIN_design,DTWIN_config):
    NEW_DTWIN_state=deepcopy(DTWIN_state)

    for i1 in DTWIN_config['Brxtor_list']:
        for i2 in DTWIN_config['Species_list']:
            if i2=='DOT':
                ts_sample_all=np.array(DTWIN_design[i1]['time_sample'][i2])
            else:
                ts_sample_all=np.array(DTWIN_design[i1]['time_sample'][i2])*(1+np.random.normal(0,1,size=len(np.array(DTWIN_design[i1]['time_sample'][i2])))*DTWIN_config['Noise_time']*0)
           
            ts_sample=ts_sample_all[(ts_sample_all>time_initial) & (ts_sample_all<=time_final)]

            tX_state=DTWIN_state[i1]['All'][i2]['time']
            X_state=DTWIN_state[i1]['All'][i2]['Value']
            
            if i2=='DOT':
                X_interp=np.interp(ts_sample,tX_state,X_state)
                X_interp=X_interp*(1+np.random.normal(0,1,size=len(X_interp))*1/5*DTWIN_config['Noise_concentration'])
            else:
                X_interp=np.interp(ts_sample,tX_state,X_state)
                X_interp=X_interp*(1+np.random.normal(0,1,size=len(X_interp))*DTWIN_config['Noise_concentration'])

            NEW_DTWIN_state[i1]['Sample'][i2]['time']=np.append(np.array(NEW_DTWIN_state[i1]['Sample'][i2]['time']),ts_sample).tolist()
            NEW_DTWIN_state[i1]['Sample'][i2]['Value']=np.append(np.array(NEW_DTWIN_state[i1]['Sample'][i2]['Value']),X_interp).tolist()    
    
    return NEW_DTWIN_state
# %% Write
def write(filename,time_initial,time_final,DTWIN_state,DTWIN_design,DTWIN_config):
    
    with open(filename) as json_file:   
        File_dict = json.load(json_file)

    time_samples_analysis=DTWIN_config['time_samples_analysis']        
    for i1 in DTWIN_config['Brxtor_list']:
      
        for i2 in DTWIN_config['Species_list']:
            if i2=='Xv':
                # tsf=(np.array(list(File_dict[i1]['measurements_aggregated']['OD600']['measurement_time'].values()))).tolist()
                # Xsf=list(File_dict[i1]['measurements_aggregated']['OD600']['OD600'].values())
                
                ts_new=np.array(DTWIN_state[i1]['Sample'][i2]['time'])
                Xsf_new=np.array(DTWIN_state[i1]['Sample'][i2]['Value'])*2.7027
                
                time_samples_iter=np.array(time_samples_analysis)
                time_samples_iter=time_samples_iter[time_samples_iter<=(time_final-1)]
                if len(time_samples_iter)>0:
                    time_analysis=time_samples_iter[-1]
                    Xsf_new=Xsf_new[ts_new<=time_analysis]
                    ts_new=ts_new[ts_new<=time_analysis]
                    tsf=(ts_new*3600).tolist()
                    Xsf=Xsf_new.tolist()
                    
                else:
                    tsf=[]
                
    
            
                File_dict[i1]['measurements_aggregated']['OD600']['measurement_time']={}
                File_dict[i1]['measurements_aggregated']['OD600']['OD600']={}
                for i3 in range(len(tsf)):
                    File_dict[i1]['measurements_aggregated']['OD600']['measurement_time'][str(i3)]=tsf[i3]
                    
                    File_dict[i1]['measurements_aggregated']['OD600']['OD600'][str(i3)]=Xsf[i3]
                    
            elif i2=='DOT':
                tsf=(np.array(list(File_dict[i1]['measurements_aggregated'][i2]['measurement_time'].values()))).tolist()
                Xsf=list(File_dict[i1]['measurements_aggregated'][i2][i2].values())

                ts_new=np.array(DTWIN_state[i1]['Sample'][i2]['time'])
                Xsf_new=np.array(DTWIN_state[i1]['Sample'][i2]['Value'])
                
                Xsf_new=Xsf_new[(ts_new>time_initial) & (ts_new<=time_final)]
                ts_new=ts_new[(ts_new>time_initial) & (ts_new<=time_final)]
                
                tsf=tsf+(ts_new*3600).tolist()
                Xsf=Xsf+Xsf_new.tolist()

                File_dict[i1]['measurements_aggregated'][i2]['measurement_time']={}
                File_dict[i1]['measurements_aggregated'][i2][i2]={}
                for i3 in range(0,len(tsf)):
                    File_dict[i1]['measurements_aggregated'][i2]['measurement_time'][str(i3)]=tsf[i3]
                    File_dict[i1]['measurements_aggregated'][i2][i2][str(i3)]=Xsf[i3]        
            else:
                # tsf=(np.array(list(File_dict[i1]['measurements_aggregated'][i2]['measurement_time'].values()))).tolist()
                # Xsf=list(File_dict[i1]['measurements_aggregated'][i2][i2].values())

                ts_new=np.array(DTWIN_state[i1]['Sample'][i2]['time'])
                Xsf_new=np.array(DTWIN_state[i1]['Sample'][i2]['Value'])
                
                time_samples_iter=np.array(time_samples_analysis)
                time_samples_iter=time_samples_iter[time_samples_iter<=(time_final-1)]
                if len(time_samples_iter)>0:
                    time_analysis=time_samples_iter[-1]
                    Xsf_new=Xsf_new[ts_new<=time_analysis]
                    ts_new=ts_new[ts_new<=time_analysis]
                    tsf=(ts_new*3600).tolist()
                    Xsf=Xsf_new.tolist()
                    
                else:
                    tsf=[]
                # Xsf_new=Xsf_new[(ts_new>time_initial) & (ts_new<=time_final)]
                # ts_new=ts_new[(ts_new>time_initial) & (ts_new<=time_final)]
                
                # tsf=tsf+(ts_new*3600).tolist()
                # Xsf=Xsf+Xsf_new.tolist()
                try:
                    File_dict[i1]['measurements_aggregated'][i2]['measurement_time']={}
                    File_dict[i1]['measurements_aggregated'][i2][i2]={}
                except:
                    File_dict[i1]['measurements_aggregated'][i2]={'measurement_time':{},i2:{}}
                for i3 in range(0,len(tsf)):
                    File_dict[i1]['measurements_aggregated'][i2]['measurement_time'][str(i3)]=tsf[i3]
                    File_dict[i1]['measurements_aggregated'][i2][i2][str(i3)]=Xsf[i3]
                            
        # Feed
        ts_pulse=(np.array(list(File_dict[i1]['measurements_aggregated']['Cumulated_feed_volume_glucose']['measurement_time'].values()))).tolist()
        F_pulse=list(File_dict[i1]['measurements_aggregated']['Cumulated_feed_volume_glucose']['Cumulated_feed_volume_glucose'].values())
        ts_dextrine=(np.array(list(File_dict[i1]['measurements_aggregated']['Cumulated_feed_volume_dextrine']['measurement_time'].values()))).tolist()
        F_dextrine=list(File_dict[i1]['measurements_aggregated']['Cumulated_feed_volume_dextrine']['Cumulated_feed_volume_dextrine'].values())
        ts_enzyme=(np.array(list(File_dict[i1]['measurements_aggregated']['Cumulated_feed_volume_enzyme']['measurement_time'].values()))).tolist()
        F_enzyme=list(File_dict[i1]['measurements_aggregated']['Cumulated_feed_volume_enzyme']['Cumulated_feed_volume_enzyme'].values())

        ts_new_pulse=np.array(DTWIN_design[i1]['Pulses']['time_pulse'])#np.array(DTWIN_state[i1]['Sample'][i2]['time'])
        F_new_pulse=np.array(DTWIN_design[i1]['Pulses']['Feed_pulse'])#np.array(DTWIN_state[i1]['Sample'][i2]['Value'])
        ts_new_dextrine=np.array(DTWIN_design[i1]['Pulses']['time_dextrine'])#np.array(DTWIN_state[i1]['Sample'][i2]['time'])
        F_new_dextrine=np.array(DTWIN_design[i1]['Pulses']['Feed_dextrine'])#np.array(DTWIN_state[i1]['Sample'][i2]['Value'])
        ts_new_enzyme=np.array(DTWIN_design[i1]['Pulses']['time_enzyme'])#np.array(DTWIN_state[i1]['Sample'][i2]['time'])
        F_new_enzyme=np.array(DTWIN_design[i1]['Pulses']['Feed_enzyme'])#np.array(DTWIN_state[i1]['Sample'][i2]['Value'])
        
        F_new_pulse=F_new_pulse[(ts_new_pulse>time_initial) & (ts_new_pulse<=time_final)]
        ts_new_pulse=ts_new_pulse[(ts_new_pulse>time_initial) & (ts_new_pulse<=time_final)].copy()
        F_new_dextrine=F_new_dextrine[(ts_new_dextrine>time_initial) & (ts_new_dextrine<=time_final)]
        ts_new_dextrine=ts_new_dextrine[(ts_new_dextrine>time_initial) & (ts_new_dextrine<=time_final)].copy()
        F_new_enzyme=F_new_enzyme[(ts_new_enzyme>time_initial) & (ts_new_enzyme<=time_final)]
        ts_new_enzyme=ts_new_enzyme[(ts_new_enzyme>time_initial) & (ts_new_enzyme<=time_final)].copy()

        ts_pulse=ts_pulse+(ts_new_pulse*3600).tolist()
        if len(F_pulse)==1:
            last_pulse=F_pulse[0]+0
        elif len(F_pulse)>1:
            last_pulse=F_pulse[-1]+0
        else:
            last_pulse=0
        F_pulse_all=F_pulse+(np.cumsum(F_new_pulse)+last_pulse).tolist()

        ts_dextrine=ts_dextrine+(ts_new_dextrine*3600).tolist()
        if len(F_dextrine)==1:
            last_pulse=F_dextrine[0]+0
        elif len(F_dextrine)>1:
            last_pulse=F_dextrine[-1]+0
        else:
            last_pulse=0
        F_dextrine_all=F_dextrine+(np.cumsum(ts_new_dextrine)+last_pulse).tolist()

        ts_enzyme=ts_enzyme+(ts_new_enzyme*3600).tolist()
        if len(F_enzyme)==1:
            last_pulse=F_enzyme[0]+0
        elif len(F_enzyme)>1:
            last_pulse=F_enzyme[-1]+0
        else:
            last_pulse=0
        F_enzyme_all=F_enzyme+(np.cumsum(F_new_enzyme)+last_pulse).tolist()
            
        File_dict[i1]['measurements_aggregated']['Cumulated_feed_volume_glucose']['measurement_time']={}
        File_dict[i1]['measurements_aggregated']['Cumulated_feed_volume_glucose']['Cumulated_feed_volume_glucose']={}
        File_dict[i1]['measurements_aggregated']['Cumulated_feed_volume_dextrine']['measurement_time']={}
        File_dict[i1]['measurements_aggregated']['Cumulated_feed_volume_dextrine']['Cumulated_feed_volume_dextrine']={}
        File_dict[i1]['measurements_aggregated']['Cumulated_feed_volume_enzyme']['measurement_time']={}
        File_dict[i1]['measurements_aggregated']['Cumulated_feed_volume_enzyme']['Cumulated_feed_volume_enzyme']={}
        for i3 in range(0,len(ts_pulse)):
            File_dict[i1]['measurements_aggregated']['Cumulated_feed_volume_glucose']['measurement_time'][str(i3)]=ts_pulse[i3]+0
            File_dict[i1]['measurements_aggregated']['Cumulated_feed_volume_glucose']['Cumulated_feed_volume_glucose'][str(i3)]=F_pulse_all[i3]+0
        for i3 in range(0,len(ts_dextrine)):
            File_dict[i1]['measurements_aggregated']['Cumulated_feed_volume_dextrine']['measurement_time'][str(i3)]=ts_dextrine[i3]+0
            File_dict[i1]['measurements_aggregated']['Cumulated_feed_volume_dextrine']['Cumulated_feed_volume_dextrine'][str(i3)]=F_dextrine_all[i3]+0
        for i3 in range(0,len(ts_enzyme)):
            File_dict[i1]['measurements_aggregated']['Cumulated_feed_volume_enzyme']['measurement_time'][str(i3)]=ts_enzyme[i3]+0
            File_dict[i1]['measurements_aggregated']['Cumulated_feed_volume_enzyme']['Cumulated_feed_volume_enzyme'][str(i3)]=F_enzyme_all[i3]+0
        
    
    with open(filename, "w") as outfile:
        json.dump(File_dict, outfile) 
        
    return File_dict
# %% Read
def read(filename,DTWIN_design,DTWIN_config):
    NEW_DTWIN_design=deepcopy(DTWIN_design)
  
    with open(filename) as json_file:   
                File_dict = json.load(json_file)

       
    for i1 in DTWIN_config['Brxtor_list']:
        f0_pulse=np.array(list(File_dict[i1]['setpoints']['Feed_glc_cum_setpoints']['Feed_glc_cum_setpoints'].values()))
        tf_pulse=np.array(list(File_dict[i1]['setpoints']['Feed_glc_cum_setpoints']['setpoint_time'].values()))/3600
        f0_dextrine=np.array(list(File_dict[i1]['setpoints']['Feed_dextrine_cum_setpoints']['Feed_dextrine_cum_setpoints'].values()))
        tf_dextrine=np.array(list(File_dict[i1]['setpoints']['Feed_dextrine_cum_setpoints']['setpoint_time'].values()))/3600
        f0_enzyme=np.array(list(File_dict[i1]['setpoints']['Feed_enzyme_cum_setpoints']['Feed_enzyme_cum_setpoints'].values()))
        tf_enzyme=np.array(list(File_dict[i1]['setpoints']['Feed_enzyme_cum_setpoints']['setpoint_time'].values()))/3600
        # f0_pulse=np.array(File_dict[i1]['setpoints']['Feed_glc_cum_setpoints']['Feed_glc_cum_setpoints'])
        # tf_pulse=np.array(File_dict[i1]['setpoints']['Feed_glc_cum_setpoints']['setpoint_time'])/3600
        # f0_dextrine=np.array(File_dict[i1]['setpoints']['Feed_dextrine_cum_setpoints']['Feed_dextrine_cum_setpoints'])
        # tf_dextrine=np.array(File_dict[i1]['setpoints']['Feed_dextrine_cum_setpoints']['setpoint_time'])/3600
        # f0_enzyme=np.array(File_dict[i1]['setpoints']['Feed_enzyme_cum_setpoints']['Feed_enzyme_cum_setpoints'])
        # tf_enzyme=np.array(File_dict[i1]['setpoints']['Feed_enzyme_cum_setpoints']['setpoint_time'])/3600

        # f_pulse=np.hstack([f0_pulse[0],np.diff(f0_pulse)])
        # f_dextrine=np.hstack([f0_dextrine[0],np.diff(f0_dextrine)])
        # f_enzyme=np.hstack([f0_enzyme[0],np.diff(f0_enzyme)])
        f_pulse=np.hstack([f0_pulse[0],np.diff(f0_pulse)])
        f_dextrine=np.hstack([f0_dextrine[0],np.diff(f0_dextrine)])
        f_enzyme=np.hstack([f0_enzyme[0],np.diff(f0_enzyme)])
        # mask_f0=f0!=None
        # f_acc=f0[mask_f0]
        # f=np.hstack([f_acc[0],np.diff(f_acc)])
        
        # tf=tf0[mask_f0]/3600

        NEW_DTWIN_design[i1]['Pulses']['time_pulse']=tf_pulse.tolist()
        NEW_DTWIN_design[i1]['Pulses']['Feed_pulse']=f_pulse.tolist()
        NEW_DTWIN_design[i1]['Pulses']['time_dextrine']=tf_dextrine.tolist()
        NEW_DTWIN_design[i1]['Pulses']['Feed_dextrine']=f_dextrine.tolist()
        NEW_DTWIN_design[i1]['Pulses']['time_enzyme']=tf_enzyme.tolist()
        NEW_DTWIN_design[i1]['Pulses']['Feed_enzyme']=f_enzyme.tolist()
        
        return NEW_DTWIN_design
        
