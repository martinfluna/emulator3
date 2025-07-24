import json
import numpy as np
import matplotlib.pyplot as plt

import time
from  scipy.signal import savgol_filter as smooth_filter 

# %% 
   
with open('db_output.json') as json_file: 
    data = json.load(json_file)
with open('Config_dot.json') as json_file:     
    Config_dot = json.load(json_file)
with open('Feed_dot.json') as json_file:
    Feed_profile = json.load(json_file)


# %% Processing data

list_br=list(data.keys())

t_last=np.zeros(len(list_br))

nn=0
for i in list_br:
    nn=nn+1
    
    tt=np.array(list(data[i]['measurements_aggregated']['DOT']['measurement_time'].values()))
    ddot_raw=np.array(list(data[i]['measurements_aggregated']['DOT']['DOT'].values()))
    
    try:
        ddot=smooth_filter(ddot_raw,window_length=5,polyorder=0) #smooth
    except:
        ddot=ddot_raw.copy()

    
    tt=tt/3600 #to hours
    t_last[nn-1]=float(tt[-1])
    

# %%  check batch
    if (tt[-1]<Config_dot['DOT_node_start']):
        Config_dot['check_batch'][i]=1
    else:
        ddot=ddot[tt>Config_dot['forget_before']] #forget first X hours
        tt=tt[tt>Config_dot['forget_before']] #forget first X hours

    # plt.plot(tt,ddot)  

    if (Config_dot['check_batch'][i]==1) and (tt[-1]>=Config_dot['DOT_node_start']):
        i_t_min=np.argmin(ddot)
        
        criteria_1=tt[i_t_min]<tt[-1]

        criteria_2=(max(ddot)-ddot[i_t_min])*Config_dot['DOT_threshold']<(ddot[-1]-ddot[i_t_min])
        
        if (criteria_1) and (criteria_2) and (tt[-1]<Config_dot['time_batchEnd'][i]):
            time_batch=round(tt[-1]*6)/6
            Config_dot['check_batch'][i]=0
            print(time_batch,'criteria 2 ',criteria_2,i,'shortened')

        elif tt[-1]>=(Config_dot['time_batch'][i]-10/60):
            time_batch=round(Config_dot['time_batch'][i]*6)/6+10/60#+20/60?
            print(time_batch,i,'delayed')

        else:
            time_batch=round(Config_dot['time_batch'][i]*6)/6
            print(time_batch,i,'no change')

        Config_dot['time_batch'][i]=time_batch+0
        Config_dot['time_batchEnd'][i]=time_batch+0
        delta_time_batch=time_batch-Feed_profile[i]['measurement_time'][0]

        Feed_profile[i]['measurement_time']=(np.array(Feed_profile[i]['measurement_time'])+delta_time_batch).tolist()#check
       
        Config_dot['time_induction'][i]=Config_dot['time_induction'][i]+delta_time_batch
    

        
# %% check pulse
    if (tt[-1]>Config_dot['time_batchEnd'][i]):
        print('pulse check')
        tt=np.array(list(data[i]['measurements_aggregated']['DOT']['measurement_time'].values()))/3600
        ddot=ddot_raw[tt>Config_dot['forget_before']]
        tt=tt[tt>Config_dot['forget_before']]
        plt.plot(tt,ddot)
        plt.show()
        if min(ddot)<Config_dot['O2min']:
            time_uu_corrected=np.array(Feed_profile[i]['measurement_time'])
            uu_corrected=np.array(Feed_profile[i]['feed_profile'])
            uu_corrected[time_uu_corrected>tt[-1]]=5
            Feed_profile[i]['feed_profile']=uu_corrected.tolist()
            Feed_profile[i]['setpoint_value']=np.cumsum(uu_corrected).tolist()#check
            # Feed_profile['NEXT']['optimize_feed'][nn-1]=0

            print('feedback on '+i)


        

# %%     
Config_dot['time_last_exec'] =  max(t_last)
plt.show()  
        
with open('Config_dot.json', "w") as outfile:
    json.dump(Config_dot, outfile)
            
with open('Feed_dot.json', "w") as outfile:
    json.dump(Feed_profile, outfile)
        

for i in list_br:
    # convert measurement time to seconds
    Feed_profile[i]['measurement_time'] = [int(m * 3600) for m in Feed_profile[i]['measurement_time']]

    # delete extra column for pulses
    del Feed_profile[i]['feed_profile']

    
with open('Feed.json', "w") as outfile:
    json.dump(Feed_profile, outfile)
    
