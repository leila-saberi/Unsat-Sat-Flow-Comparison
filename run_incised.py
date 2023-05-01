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

mf6 = True
base_dir = 'mf6_sce2_refined'#'test1_newBCs_Sce2_refined'#'mf6_sce2_refined'

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
    zon[:, 27:54, 1] = 4
    # River zone:
    zon[:, 27:54, -1] = 2
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
    zbud6_exe = os.path.join(os.getcwd(), base_dir)
    # zbud_nam = os.path.join(base_dir, 'zonebud.zbnam')

    # args = f'{zbud6_exe}', f'{zbud_nam}'

    # zbud6_exe = r'C:/Users/lsaberi/GitHub/Unsat-Sat-Flow-Comparison/mf6_sce0_refined/zonebudget.bat'

    p = subprocess.Popen(["zonebudget.bat"], shell=True, cwd=zbud6_exe)

    stdout, stderr = p.communicate()
    print(p.returncode)  # is 0 if success

    if p.returncode == 0:
        zb = zonbud.get_budget()
        df = zonbud.get_dataframes().reset_index()
        flow_DNG = df['ZONE_2'][7] * 28.3168  # SS rate into river (L/day)
        flow_UPG = df['ZONE_4'][7] * 28.3168  # SS rate into river (L/day)
        flow_pond = df['ZONE_3'][7] * 28.3168  # SS effluent rate from pond (L/day)
        flow_DNG_UZF = df['ZONE_2'][2] * 28.3168  # UZF SS rate into river (L/day)
        flow_UPG_UZF = df['ZONE_4'][2] * 28.3168  # UZF SS rate into river (L/day)
        flow_pond_UZF = df['ZONE_3'][2] * 28.3168  # UZF SS effluent rate from pond (L/day)
        # flow_well = 0.5 * (df['ZONE_1'][5] + df['ZONE_4'][2]) * 28.3168  # SS flow through well (L/day)
        print(df.to_string())
        print('Flow into river: ' + str(flow_DNG) + ' L/day.')
        print('Flow from pond: ' + str(flow_pond) + ' L/day.')
        print('Flow through UPG: ' + str(flow_UPG) + ' L/day.')
        print('UZF Flow into river: ' + str(flow_DNG_UZF) + ' L/day.')
        print('UZF Flow from pond: ' + str(flow_pond_UZF) + ' L/day.')
        print('UZF Flow through UPG: ' + str(flow_UPG_UZF) + ' L/day.')
        print('Percent: ' + str(flow_pond / flow_DNG * 100.0))
        exit()

else:
    c_riv, c_well, flows = ccr.incised_vertical.model_output(base_dir=base_dir)
    print('River concentration: ' + str(c_riv))
    print('Well concentration: ' + str(c_well))
    print(flows)



