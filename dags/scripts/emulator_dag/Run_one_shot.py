import time
import Node_run_emulator
import Node_start_emulator

# %%
Node_start_emulator.start_emu()
time_a=time.time()
time.sleep(0.934) 
# %%
time_c=time.time()   
Node_run_emulator.run_emu()
time_d=time.time()    

print(time_d-time_c)