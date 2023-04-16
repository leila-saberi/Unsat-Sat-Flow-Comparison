"""
Script to convert modflow-nwt to modflow6
Authors: Ahmad Askar + Lisbon team
Modified for EPRI by: Michelle Pedrazas

"""

import os
import shutil
import flopy
import sys

# Specifying the directories of the old and new models
cwd = os.getcwd()
org_d = os.path.join(cwd, "test1_newBCs_Sce2_refined") #MODFLOW-NWT
org_trans_d = os.path.join(cwd, "test1_newBCs_Sce2_refined")
new_d = os.path.join(cwd, "mf6_sce2_refined") #MODFLOW6
new_trans_d = os.path.join(cwd, "mf6trans_sce2_refined")
no_et_mass_out = False

# removing the new model directory if exist
if os.path.exists(new_d):
    shutil.rmtree(new_d)
if os.path.exists(new_trans_d):
    shutil.rmtree(new_trans_d)

# reading the nwt model
mf_old = flopy.modflow.Modflow.load("incised.nam", model_ws=org_d,
                                    verbose=False,
                                    check=False,
                                    exe_name="mfnwt.exe")

# reading the mt3d model
mt_old = flopy.mt3d.mt.Mt3dms.load("incised_mt3.nam", model_ws=org_trans_d,
                                   verbose=False)

############################################# MF6 Simulation
# Creating the MF6 Simulation
new_sim = flopy.mf6.MFSimulation(sim_name="incised", version="mf6", exe_name="mf6.exe", sim_ws=new_d)

# TDIS Package for MF6 Simulation
flopy.mf6.ModflowTdis(new_sim,
                      nper=mf_old.dis.nper,
                      perioddata=[mf_old.dis.get_final_totim(), #PERLEN: length of a stress period
                                  mf_old.dis.nstp.array[0],     #NSTP: number of time steps in a stress period
                                  mf_old.dis.tsmult.array[0]],  #TSMULT: multiplier for the length of successive time steps
                      filename=f"{new_sim.name}.tdis",          #filename: File name for this package
                      pname = f"{new_sim.name}")                #pname: Package name for this package

############################################## Groundwater Flow Model
#Build GWF model
flow_model = flopy.mf6.ModflowGwf(new_sim, modelname=f"{new_sim.name}", model_nam_file=f"{new_sim.name}.nam")
# DIS Package for GWF model (Gwfdis)
flopy.mf6.ModflowGwfdis(flow_model, nlay=mf_old.dis.nlay,
                        nrow=mf_old.dis.nrow, ncol=mf_old.dis.ncol, delr=mf_old.dis.delr, delc=mf_old.dis.delc,
                        top=mf_old.dis.top.array, botm=mf_old.dis.botm.array, idomain=mf_old.modelgrid.idomain)

# flopy.mf6.ModflowGwfupw(flow_model)

# IC Package for GWF model (Gwfic)
flopy.mf6.ModflowGwfic(flow_model, strt=mf_old.bas6.strt.array)

# NPF Package for GWF model (Gwfnpf)
flopy.mf6.ModflowGwfnpf(flow_model, save_flows=True, save_specific_discharge=True,
                        save_saturation=True, k=mf_old.upw.hk.array, k33=mf_old.upw.vka.array, icelltype=1)

# OC package for GWF model (Gwfoc)
flopy.mf6.ModflowGwfoc(flow_model, budget_filerecord=f"{new_sim.name}.cbb",head_filerecord=f"{new_sim.name}.hds",
                       saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
                       printrecord=[("BUDGET", "ALL")])

#CHD Package
Old_chd=mf_old.chd.stress_period_data.get_dataframe()
chd_period = {}
for per in range(mf_old.dis.nper):
    chb_period_array = []
    for k,i,j,shead in zip(Old_chd.k, Old_chd.i, Old_chd.j, Old_chd.shead):
        chb_period_array.append(((k, i, j), shead))
    chd_period[per] = chb_period_array
flopy.mf6.ModflowGwfchd(flow_model,stress_period_data=chd_period)

# IMS package for MF6 Simulation
ims = flopy.mf6.ModflowIms(new_sim)
new_sim.register_ims_package(ims, [new_sim.name])

# Write MF6 Simulation
new_sim.write_simulation()

# Run Simulation
success, mfoutput = new_sim.run_simulation(silent=False, pause=False, report=True)

print("Done w Flow~")
# sys.exit()

""" Transport model """  ###########################################################################################################################################

# Old mt3d trans model
sim_trans_old = mt_old

######### Build Simulation
sim_trans = flopy.mf6.MFSimulation(sim_ws=new_trans_d, continue_=True)
# TDIS Package
flopy.mf6.ModflowTdis(sim_trans, nper=sim_trans_old.btn.nper, filename="incised_mt3d.tdis")

############ Build Groundwater Transport Model
mtrans = flopy.mf6.ModflowGwt(sim_trans, modelname="incised_mt3d", save_flows=True)

# DIS Package for Gwt
flopy.mf6.ModflowGwtdis(mtrans, nlay=mf_old.dis.nlay,
                        nrow=mf_old.dis.nrow, ncol=mf_old.dis.ncol, delr=mf_old.dis.delr, delc=mf_old.dis.delc,
                        top=mf_old.dis.top.array, botm=mf_old.dis.botm.array, idomain=mf_old.modelgrid.idomain)

# IC Package
ic_r_mt3d = sim_trans_old.btn.sconc[0].array
flopy.mf6.ModflowGwfic(mtrans, strt=ic_r_mt3d)

# Adv Package
flopy.mf6.ModflowGwtadv(mtrans)  # ,scheme="UPSTREAM")

# DSP Package
dsp_r_mt3d = sim_trans_old.dsp.al.array
flopy.mf6.ModflowGwtdsp(mtrans)  # ,xt3d_off=True,alh=dsp_r_mt3d,ath1=dsp_r_mt3d*0.1,ath2 =dsp_r_mt3d*0.01,diffc=07.8600e-06)

# SSM package
flopy.mf6.ModflowGwtssm(mtrans, save_flows=True)

#MST6 Package
flopy.mf6.ModflowGwtmst(mtrans)

# OC package
flopy.mf6.ModflowGwtoc(mtrans, budget_filerecord="incised_mt3d.cbc", concentration_filerecord="incised_mt3d.ucn",
                       saverecord=[("CONCENTRATION", "last"), ("BUDGET", "last")], printrecord=[("BUDGET", "last")])

# IMS package
# ims = flopy.mf6.ModflowIms(sim_trans, pname="ims",filename="mf6_{0}.ims".format(org_trans_d),complexity="complex",outer_dvclose=1e-04,
#                           print_option ="all" )
ims = flopy.mf6.ModflowIms(sim_trans)
sim_trans.register_ims_package(ims, [mtrans.name])
sim_trans.write_simulation()

print("Done w Transport~")