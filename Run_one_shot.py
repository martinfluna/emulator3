import time
import Node_run_emulator
import Node_start_emulator
# import numpy as np
# import json

# %%
# Node_design_file.create()
# Node_createDesign_json.create()

Node_start_emulator.start_emu()
time_a=time.time()
# time_b=time.time()
time.sleep(0.934) 
# %%
time_c=time.time()   
# Node_run_emulator.run_emu()
Node_run_emulator.run_emu()
time_d=time.time()    

print(time_d-time_c)