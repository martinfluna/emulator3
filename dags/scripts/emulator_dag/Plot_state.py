import numpy as np
import pandas as pd
import json
import time

import matplotlib.pyplot as plt
# %% 


with open('EMULATOR_state.json') as json_file:   
        EMULATOR_state = json.load(json_file)
with open('EMULATOR_design.json') as json_file:   
        EMULATOR_design = json.load(json_file)
with open('EMULATOR_config.json') as json_file:   
        EMULATOR_config = json.load(json_file)

mbr='19419'
species='Dextrine_R' # ['Xv','Glucose','Acetate','DOT','Fluo_RFP','Volume','Dextrine_S','Dextrine_R','Enzyme']

t=EMULATOR_state[mbr]['All'][species]['time']
y=EMULATOR_state[mbr]['All'][species]['Value']        

plt.plot(t,y,'.')