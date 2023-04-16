import ccr
# import ccr.incised as inc
import os
import flopy
from flopy.utils import ZoneBudget
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
# import pyemu
import subprocess

mf6 = False
base_dir = 'test1_newBCs_Sce2_refined'#'mf6_sce2_refined'

if mf6 == True:
    modelname = 'incised'
    exe = 'mf6.exe'

    #To get the zone locs based on the KH
    temp_dir = 'mf6_sce1_refined'
    ws = os.path.join(os.getcwd(), temp_dir)
    sim = flopy.mf6.MFSimulation.load(sim_name=f'{modelname}', version='mf6', exe_name=exe, sim_ws=ws)
    sim.simulation_data.mfpath.set_sim_path(ws)
    m = sim.get_model()
    bas = m.dis.idomain.array
    bas_shape = bas.shape
    npf = m.npf.k.array

   #To get the zone budget
    ws = os.path.join(os.getcwd(), base_dir)
    sim = flopy.mf6.MFSimulation.load(sim_name=f'{modelname}', version='mf6', exe_name=exe, sim_ws=ws)
    sim.simulation_data.mfpath.set_sim_path(ws)
    m = sim.get_model()

    cbb_file = os.path.join(base_dir, modelname + ".cbb")
    zon = np.ones(bas_shape, dtype=int)
    # upgradien zone
    zon[:, 38:77, 1] = 4
    # River zone:
    zon[:, 38:77, -1] = 2
    # Pond zone:
    zon_rows, zon_cols = np.nonzero(bas[0, :, :])
    # zon[0,np.nonzero(bas.ibound[0,:,:])] = 3
    # zon[0,zon_rows,zon_cols] = 3
    # zon[pond_cells] = 3
    zon[npf[:,:,:]<0.3] = 3
    # Collector well:
    # zon[-1,well_cells] = 4
    # zon[-3, well_cells[0], well_cells[1]] = 4
    zonbud = m.output.zonebudget(zon)
    zonbud.change_model_ws(ws)
    zonbud.write_input()
    zbud6_exe = os.path.join(base_dir, 'zbud6.exe')
    zbud_nam = os.path.join(base_dir, 'zonebud.zbnam')

    # subprocess.call([f'{zbud6_exe}', f'{zbud_nam}'])

    zb = zonbud.get_budget()
    df = zonbud.get_dataframes().reset_index()
    flow_DNG = df['ZONE_2'][5] * 28.3168  # SS rate into river (L/day)
    flow_UPG = df['ZONE_4'][5] * 28.3168  # SS rate into river (L/day)
    flow_pond = df['ZONE_3'][5] * 28.3168  # SS effluent rate from pond (L/day)
    # flow_well = 0.5 * (df['ZONE_1'][5] + df['ZONE_4'][2]) * 28.3168  # SS flow through well (L/day)
    print(df.to_string())
    print('Flow into river: ' + str(flow_DNG) + ' L/day.')
    print('Flow from pond: ' + str(flow_pond) + ' L/day.')
    print('Flow through UPG: ' + str(flow_UPG) + ' L/day.')
    print('Percent: ' + str(flow_pond / flow_DNG * 100.0))
    exit()

else:
    c_riv, c_well, flows = ccr.incised_vertical.model_output(base_dir=base_dir)
    print('River concentration: ' + str(c_riv))
    print('Well concentration: ' + str(c_well))
    print(flows)



