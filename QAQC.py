
import os
import pandas as pd
import numpy as np
import geopandas as gpd
pd.set_option('display.max_columns', 10)
import flopy
from collections import OrderedDict
import shutil
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as ticker
from scipy import stats
import statistics

new_d = 'mf6_sce0_refined'
org_d = 'test1_newBCs_Sce0'
    
def compare_heads():

    old_hds = flopy.utils.HeadFile(os.path.join(org_d,'incised.hds'))
    old_hds=old_hds.get_data(kstpkper=(old_hds.get_kstpkper()[0]))

    new_hds = flopy.utils.HeadFile((os.path.join(new_d,'incised.hds')))
    new_hds=new_hds.get_data(kstpkper=(new_hds.get_kstpkper()[0]))

    new_wt=np.zeros(shape=(m.dis.nrow.array,m.dis.ncol.array))
    old_wt=np.zeros(shape=(m.dis.nrow.array,m.dis.ncol.array))

    old_id=nwt_m.modelgrid.idomain
    old_bot=nwt_m.dis.botm.array

    new_id=m.dis.idomain.array
    new_bot=m.dis.botm.array

    # for i in np.arange(nwt_m.dis.nrow):
    #         for j in np.arange(nwt_m.dis.ncol):
    #                 for k in np.arange(nwt_m.dis.nlay):
    #                     if old_id[k,i,j]==1 and old_hds[k,i,j]>old_bot[k,i,j]:
    #                         old_wt[i,j]=old_hds[k,i,j]
    #                         break
    #
    # for i in np.arange(m.dis.nrow.array):
    #         for j in np.arange(m.dis.ncol.array):
    #                 for k in np.arange(m.dis.nlay.array):
    #                     if new_id[k,i,j]==1 and new_hds[k,i,j]>new_bot[k,i,j]:
    #                         new_wt[i,j]=new_hds[k,i,j]
    #                         break
    #
    # diff_hds=new_wt-old_wt
    # # diff_hds[diff_hds<0.01] = np.nan
    # fig1, axes = plt.subplots(1, 1, figsize=(7, 7.5),dpi=300)
    # ax=axes
    # cb=ax.imshow(diff_hds)
    # plt.colorbar(cb,orientation='horizontal',ax=ax,pad=0.1)
    # plt.legend()
    # plt.show()
    #
    # fig2, axes = plt.subplots(1, 1, figsize=(7, 7.5),dpi=300)
    # ax=axes
    # levels=np.arange(-500,500,5)
    # co1=ax.contour(old_wt, colors='k', origin='upper',levels=levels,label="Old_Model")
    # co2=ax.contour(new_wt, colors='r', origin='upper',levels=levels,linestyles='dashed',label="Updated_Model")
    # # levels=np.arange(1,35,3)
    # ax.clabel(co1, inline=True, fontsize=12,levels=levels,inline_spacing=5)
    # plt.tight_layout()
    # plt.legend()
    # plt.show()

    fig3, axes = plt.subplots(1, 2, figsize=(14, 7.5), dpi=300)
    ax = axes[0]
    cb = ax.imshow(new_hds[10], cmap='jet', vmin=0, vmax=35)
    c1 = plt.colorbar(cb, orientation='horizontal', ax=ax, pad=0.05, label="Head (m)")
    ax.text(15, 15, "MF6 Heads @ Lay10", fontsize=15)
    ax = axes[1]
    cb = ax.imshow(old_hds[1], cmap='jet', vmin=0, vmax=35)
    c2 = plt.colorbar(cb, orientation='horizontal', ax=ax, pad=0.05, label="Head (m)")
    # ticks=[-1,1]
    # kk=["gain","loss"]
    # c2.ax.set_xticklabels(["Taking","Giving"])
    ax.text(15, 15, "NWT Heads @ Lay1", fontsize=15)
    plt.tight_layout()
    plt.legend()
    plt.show()

    print('Head Plots are made')

def compare_saturation():
    old_hds = flopy.utils.HeadFile(os.path.join(org_d, 'incised.hds'))
    old_hds = old_hds.get_data(kstpkper=(old_hds.get_kstpkper()[0]))

    new_hds = flopy.utils.HeadFile((os.path.join(new_d, 'incised.hds')))
    new_hds = new_hds.get_data(kstpkper=(new_hds.get_kstpkper()[0]))

    # saturation comparison
    new_sat = (new_hds - m.dis.botm.array) / thick_new * 100
    new_sat[new_sat <= 0] = np.nan
    new_sat[new_sat > 100] = 100
    new_sat[new_id == 0] = np.nan
    # new_id[k,i,j]

    old_sat = (old_hds - nwt_m.dis.botm.array) / thick_old * 100
    old_sat[old_sat <= 0] = np.nan
    old_sat[old_sat > 100] = 100
    old_sat[new_id == 0] = np.nan

    for i in range(8):
        fig, axes = plt.subplots(1, 2, figsize=(14, 7.5), dpi=300)
        ax = axes[0]
        cb = ax.imshow(new_sat[i], cmap="jet", vmin=0, vmax=100)
        c = plt.colorbar(cb, orientation='horizontal', ax=ax, pad=0.05)
        c.set_label(label="S%", size=16, weight='bold')
        c.ax.tick_params(labelsize=15)
        ax.text(20, 20, "MF6 saturation @ Lay {0}".format(i + 1), fontsize=18)
        ax = axes[1]
        cb = ax.imshow(old_sat[i], cmap="jet", vmin=0, vmax=100)
        c = plt.colorbar(cb, orientation='horizontal', ax=ax, pad=0.05)
        c.set_label(label="S%", size=16, weight='bold')
        c.ax.tick_params(labelsize=15)
        ax.text(20, 20, "NWT saturation  @ Lay {0}".format(i + 1), fontsize=18)
    plt.tight_layout()

def otherComparisons():
    # comparing layer heads
    fig, axes = plt.subplots(1, 2, figsize=(14, 7.5),dpi=300)
    ax=axes[0]
    cb=ax.imshow(new_hds[4],cmap='jet', vmin=0,vmax=35)
    c1=plt.colorbar(cb,orientation='horizontal',ax=ax,pad=0.05, label="Head (m)")
    ax.text(15,15,"MF6 Heads @ Lay5",fontsize=15)
    ax=axes[1]
    cb=ax.imshow(old_hds[4],cmap='jet', vmin=0,vmax=35)
    c2=plt.colorbar(cb,orientation='horizontal',ax=ax,pad=0.05, label="Head (m)")
    #ticks=[-1,1]
    # kk=["gain","loss"]
    # c2.ax.set_xticklabels(["Taking","Giving"])
    ax.text(15,15,"NWT Heads @ Lay5",fontsize=15)
    plt.tight_layout()

    fig, axes = plt.subplots(1, 2, figsize=(14, 7.5),dpi=300)
    ax=axes[0]
    cb=ax.imshow(new_hds[0],cmap='jet', vmin=0,vmax=35)
    c1=plt.colorbar(cb,orientation='horizontal',ax=ax,pad=0.05, label="Head (m)")
    ax.text(15,15,"MF6 Heads @ Lay1",fontsize=15)
    ax=axes[1]
    cb=ax.imshow(old_hds[0],cmap='jet', vmin=0,vmax=35)
    c2=plt.colorbar(cb,orientation='horizontal',ax=ax,pad=0.05, label="Head (m)")
    #ticks=[-1,1]
    # kk=["gain","loss"]
    # c2.ax.set_xticklabels(["Taking","Giving"])
    ax.text(15,15,"NWT Heads @ Lay1",fontsize=15)
    plt.tight_layout()

    # layer thickness comparison

    top_old=np.ones(shape=(nwt_m.dis.nlay,nwt_m.dis.nrow,nwt_m.dis.ncol))
    top_old[0]=nwt_m.dis.top.array
    top_old[1:]=nwt_m.dis.botm.array[:-1]
    thick_old=top_old-nwt_m.dis.botm.array

    top_new=np.ones(shape=(m.dis.nlay.array,m.dis.nrow.array,m.dis.ncol.array))
    top_new[0]=m.dis.top.array
    top_new[1:]=m.dis.botm.array[:-1]
    thick_new=top_new-m.dis.botm.array

    thick_diff=thick_new-thick_old
    thick_diff[thick_diff==0]=np.nan
    fig, axes = plt.subplots(2, 4, figsize=(28, 15),dpi=300)
    axe=axes.ravel()
    for i in range(8):
        ax=axe[i]
        cb=ax.imshow(thick_diff[i],cmap='seismic', vmin=-3,vmax=3)
        c=plt.colorbar(cb,orientation='horizontal',ax=ax,pad=0.05)
        c.set_label(label="m", size=16, weight='bold')
        c.ax.tick_params(labelsize=15)
        ax.text(20,20,"MF6 thick - MT3d thick, @ Lay {0}".format(i+1), fontsize=18)
    plt.tight_layout()

    fig, axes = plt.subplots(2, 4, figsize=(28, 15),dpi=300)
    axe=axes.ravel()
    for i in range(8):
        ax=axe[i]
        cb=ax.imshow(thick_new[i],cmap='seismic')
        c=plt.colorbar(cb,orientation='horizontal',ax=ax,pad=0.05)
        c.set_label(label="m", size=16, weight='bold')
        c.ax.tick_params(labelsize=15)
        ax.text(20,20,"MF6 thick, @ Lay {0}".format(i+1), fontsize=18)
    plt.tight_layout()



    #ET comparison
    MF6_cbb = bf.CellBudgetFile(os.path.join(new_d,'mf6_flow.cbb'))
    NWT_cbb = bf.CellBudgetFile(os.path.join("flow",'flow.cbb'))

    MF6_evt_flux=np.zeros_like(m.dis.botm.array)
    MF6_evt=MF6_cbb.get_record(6)
    nrow=m.dis.nrow.data
    ncol=m.dis.ncol.data
    # or try to use mf6_drn_array=MF6_cbb.get_record(6,full3D=True)
    for i in range(len(MF6_evt)):
        lay_no=np.ceil(MF6_evt[i][0]/nrow/ncol)
        row_no=np.ceil((MF6_evt[i][0]-(np.floor(MF6_evt[i][0]/nrow/ncol)*nrow*ncol))/ncol)
        col_no=MF6_evt[i][0]-((lay_no-1)*nrow*ncol)-((row_no-1)*ncol)
        MF6_evt_flux[int(lay_no-1),int(row_no-1),int(col_no-1)]=MF6_evt[i][2]
    mf6_evt_array=[MF6_evt_flux[0]+MF6_evt_flux[4]][0]
    mf6_evt_array=MF6_evt_flux.sum(axis=0)
    mf6_evt_array[mf6_evt_array==0]=np.nan

    NWT_evt_array=NWT_cbb.get_record(5)[0]+NWT_cbb.get_record(5)[4]
    NWT_evt_array=NWT_cbb.get_record(5).sum(axis=0)
    NWT_evt_array[NWT_evt_array==0]=np.nan

    fig, axes = plt.subplots(1, 2, figsize=(14, 7.5),dpi=300)
    ax=axes[0]
    cb=ax.imshow(mf6_evt_array,cmap='gist_ncar', vmin=-1.6,vmax=0)
    c1=plt.colorbar(cb,orientation='horizontal',ax=ax,pad=0.05, label="ET flux (m3/d)")
    ax.text(15,15,"MF6 ET flux",fontsize=15)
    ax=axes[1]
    cb=ax.imshow(NWT_evt_array,cmap='gist_ncar', vmin=-1.6,vmax=0)
    c2=plt.colorbar(cb,orientation='horizontal',ax=ax,pad=0.05, label="ET flux (m3/d)")
    #ticks=[-1,1]
    # kk=["gain","loss"]
    # c2.ax.set_xticklabels(["Taking","Giving"])
    ax.text(15,15,"NWT ET flux",fontsize=15)
    plt.tight_layout()

    #GHB comparison
    MF6_GHB_flux=np.zeros_like(m.dis.botm.array)
    MF6_GHB=MF6_cbb.get_record(4)
    nrow=m.dis.nrow.data
    ncol=m.dis.ncol.data
    # or try to use mf6_drn_array=MF6_cbb.get_record(4,full3D=True)
    for i in range(len(MF6_GHB)):
        lay_no=np.ceil(MF6_GHB[i][0]/nrow/ncol)
        row_no=np.ceil((MF6_GHB[i][0]-(np.floor(MF6_GHB[i][0]/nrow/ncol)*nrow*ncol))/ncol)
        col_no=MF6_GHB[i][0]-((lay_no-1)*nrow*ncol)-((row_no-1)*ncol)
        MF6_GHB_flux[int(lay_no-1),int(row_no-1),int(col_no-1)]=MF6_GHB[i][2]
    mf6_GHB_array=MF6_GHB_flux[4]
    mf6_GHB_array[mf6_GHB_array==0]=np.nan

    NWT_GHB_array=NWT_cbb.get_record(6)[4]
    NWT_GHB_array[NWT_GHB_array==0]=np.nan

    fig, axes = plt.subplots(1, 2, figsize=(14, 7.5),dpi=300)
    ax=axes[0]
    cb=ax.imshow(mf6_GHB_array,cmap='gist_ncar', vmin=-4.6,vmax=7.5)
    c1=plt.colorbar(cb,orientation='horizontal',ax=ax,pad=0.05, label="GHB flux (m3/d)")
    ax.text(15,15,"MF6 GHB flux",fontsize=15)
    ax=axes[1]
    cb=ax.imshow(NWT_GHB_array,cmap='gist_ncar', vmin=-4.6,vmax=7.5)
    c2=plt.colorbar(cb,orientation='horizontal',ax=ax,pad=0.05, label="GHB flux (m3/d)")
    #ticks=[-1,1]
    # kk=["gain","loss"]
    # c2.ax.set_xticklabels(["Taking","Giving"])
    ax.text(15,15,"NWT GHB flux",fontsize=15)
    plt.tight_layout()

    #DRN comparison
    MF6_drn_flux=np.zeros_like(m.dis.botm.array)
    MF6_drn=MF6_cbb.get_record(3)
    nrow=m.dis.nrow.data
    ncol=m.dis.ncol.data
    for i in range(len(MF6_drn)):
        lay_no=np.ceil(MF6_drn[i][0]/nrow/ncol)
        row_no=np.ceil((MF6_drn[i][0]-(np.floor(MF6_drn[i][0]/nrow/ncol)*nrow*ncol))/ncol)
        col_no=MF6_drn[i][0]-((lay_no-1)*nrow*ncol)-((row_no-1)*ncol)
        MF6_drn_flux[int(lay_no-1),int(row_no-1),int(col_no-1)]=MF6_drn[i][2]
    # or try to use mf6_drn_array=MF6_cbb.get_record(3,full3D=True)
    mf6_drn_array=[MF6_drn_flux[0]+MF6_drn_flux[3]][0]
    mf6_drn_array[mf6_drn_array==0]=np.nan

    NWT_drn_array=NWT_cbb.get_record(4)[0]+NWT_cbb.get_record(4)[3]
    NWT_drn_array[NWT_drn_array==0]=np.nan

    fig, axes = plt.subplots(1, 2, figsize=(14, 7.5),dpi=300)
    ax=axes[0]
    cb=ax.imshow(mf6_drn_array,cmap='tab20', vmin=-71,vmax=0)
    c1=plt.colorbar(cb,orientation='horizontal',ax=ax,pad=0.05, label="DRN flux (m3/d)")
    ax.text(15,15,"MF6 DRN flux",fontsize=15)
    ax=axes[1]
    cb=ax.imshow(NWT_drn_array,cmap='tab20', vmin=-71,vmax=0)
    c2=plt.colorbar(cb,orientation='horizontal',ax=ax,pad=0.05, label="DRN flux (m3/d)")
    #ticks=[-1,1]
    # kk=["gain","loss"]
    # c2.ax.set_xticklabels(["Taking","Giving"])
    ax.text(15,15,"NWT DRN flux",fontsize=15)
    plt.tight_layout()

    #vertical flow comparison
    flowja = MF6_cbb.get_record(0)[0]
    MF6_face_flows=flopy.mf6.utils.get_structured_faceflows(flowja, grb_file=os.path.join("mf6_flow","mf6_flow.dis.grb"))
    mf6_vflow_array=MF6_face_flows[0]
    mf6_vflow_array[mf6_vflow_array==0]=np.nan

    NWT_vflow_array=NWT_cbb.get_record(1)
    NWT_vflow_array[NWT_vflow_array==0]=np.nan

    fig, axes = plt.subplots(1, 2, figsize=(14, 7.5),dpi=300)
    ax=axes[0]
    cb=ax.imshow(mf6_vflow_array[10],cmap='tab20')#,vmin=-15,vmax=4)
    c1=plt.colorbar(cb,orientation='horizontal',ax=ax,pad=0.05, label="vflow flux (m3/d)")
    ax.text(15,15,"MF6 Lower face flux @ Lay1",fontsize=15)
    ax=axes[1]
    cb=ax.imshow(NWT_vflow_array[10],cmap='tab20')#,vmin=-15,vmax=4)
    c2=plt.colorbar(cb,orientation='horizontal',ax=ax,pad=0.05, label="vflow flux (m3/d)")
    #ticks=[-1,1]
    # kk=["gain","loss"]
    # c2.ax.set_xticklabels(["Taking","Giving"])
    ax.text(15,15,"NWT Lower face flux @ Lay1",fontsize=15)
    plt.tight_layout()

    v_flow_diff=mf6_vflow_array-NWT_vflow_array
    v_flow_diff[v_flow_diff==0]=np.nan

    fig, axes = plt.subplots(2, 1, figsize=(18, 10),dpi=300)
    ax=axes[0]
    xsect = flopy.plot.PlotCrossSection(ax=axes[0],model=nwt_m, line={"Column": 175},extent=(5500,8000,0,30))
    csa = xsect.plot_array(NWT_vflow_array,cmap="seismic",vmin=-10, vmax=10)
    patches = xsect.plot_ibound()
    linecollection = xsect.plot_grid()
    t = ax.set_title(
        "Column 217 Cross-Section in right face flux array - NWT/mt3d model", fontsize=18)
    #cb = plt.colorbar(csa, ax=axes[1],orientation='horizontal')
    ax=axes[1]
    xsect = flopy.plot.PlotCrossSection(ax=axes[1],model=m, line={"Column": 175},extent=(5500,8000,-0,30))
    csa = xsect.plot_array(mf6_vflow_array,cmap="seismic",vmin=-10, vmax=10)
    patches = xsect.plot_ibound()
    linecollection = xsect.plot_grid()
    t = ax.set_title(
        "Column 217 Cross-Section in right face flux array - MF6 model", fontsize=18)
    #cb = plt.colorbar(csa, ax=axes[1],orientation='horizontal')

    # plotting cross-sections of conc and other arrays
    conc_new_obj = flopy.utils.HeadFile((os.path.join(new_trans_d,'{0}.ucn'.format(new_trans_d))),text="concentration")
    conc_obj = flopy.utils.HeadFile(os.path.join(org_trans_d,'{0}.ucn'.format(org_trans_d)),text="concentration")
    conc_new=conc_new_obj.get_data(kstpkper=(conc_new_obj.get_kstpkper()[25]))
    conc=conc_obj.get_data(kstpkper=(conc_obj.get_kstpkper()[25]))
    conc_new[conc_new==1e30] = 0
    conc_new[new_sat==0]=0
    conc_new[id_arr==0]=0
    conc_new[conc_new>=0] = np.nan
    conc[conc>=0] = np.nan
    conc[id_arr==0]=0


    fig, axes = plt.subplots(2, 1, figsize=(18, 10),dpi=300)
    ax=axes[0]
    xsect = flopy.plot.PlotCrossSection(ax=axes[0],model=nwt_m, line={"Column": 217},extent=(5200,6500,-30,30))
    csa = xsect.plot_array(conc,cmap="jet")
    patches = xsect.plot_ibound()
    linecollection = xsect.plot_grid()
    t = ax.set_title(
        "Column 217 Cross-Section in conc array - NWT/mt3d model", fontsize=18)
    #cb = plt.colorbar(csa, ax=axes[1],orientation='horizontal')
    ax=axes[1]
    xsect = flopy.plot.PlotCrossSection(ax=axes[1],model=m, line={"Column": 217},extent=(5200,6500,-30,30))
    csa = xsect.plot_array(conc_new,cmap="jet")
    patches = xsect.plot_ibound()
    linecollection = xsect.plot_grid()
    t = ax.set_title(
        "Column 217 Cross-Section in conc array - MF6 model", fontsize=18)
    #cb = plt.colorbar(csa, ax=axes[1],orientation='horizontal')

    # making head cross-section in the usaturated zone
    from flopy.plot import styles
    old_k=nwt_m.upw.hk.array
    hsu_data = pd.read_csv(os.path.join(org_d,"hsu_zones.DAT"), header=None, delimiter=r"\s+")
    hsu_data_arr = hsu_data.values.reshape(43,mt3d_m.btn.nrow,mt3d_m.btn.ncol)
    layers=np.zeros_like(conc)
    for i, k in enumerate(np.arange(nwt_m.dis.nlay)):
        if i>19:
            layers[k]=19
        else:
           layers[k]=i
    with styles.USGSMap():
        fig, axes = plt.subplots(2, 1, figsize=(18, 10),dpi=300)
        ax=axes[0]
        cmap="gist_rainbow"
        # force the first color entry to be greyax=axes[0]
        xsect = flopy.plot.PlotCrossSection(ax=axes[0],model=nwt_m, line={"Row": 207},extent=(5400,7000,-20,30))
        csa = xsect.plot_array(hsu_data_arr,cmap=cmap,alpha=0.7)
        patches = xsect.plot_ibound()
        levels = np.arange(10, 35, 1)
        contour_set = xsect.contour_array(old_hds, levels=levels, colors="k")
        plt.clabel(contour_set, fmt="%.1f", colors="k", fontsize=11)
        wt = xsect.plot_surface(old_wt, color="blue", lw=3)
        linecollection = xsect.plot_grid()
        t = ax.set_title("Row 207 Cross-Section in grid and heads - NWT model", fontsize=18)
        ax=axes[1]
        xsect = flopy.plot.PlotCrossSection(ax=axes[1],model=nwt_m, line={"Column": 200},extent=(5400,6600,-20,30))
        csa = xsect.plot_array(hsu_data_arr,cmap=cmap,alpha=0.7)
        patches = xsect.plot_ibound()
        contour_set = xsect.contour_array(old_hds, levels=levels, colors="k")
        plt.clabel(contour_set, fmt="%.1f", colors="k", fontsize=11)
        wt = xsect.plot_surface(old_wt, color="blue", lw=3)
        linecollection = xsect.plot_grid()
        t = ax.set_title(
            "Column 200 Cross-Section in grid and heads - NWT model", fontsize=18)
        #cb = plt.colorbar(csa, ax=axes[1],orientation='horizontal')

if __name__ == '__main__':
    sim_mf6 = flopy.mf6.MFSimulation.load(sim_ws=new_d)
    m = sim_mf6.get_model()

    nwt_m = flopy.modflow.Modflow.load("incised.nam", model_ws=org_d,
                                       verbose=False,
                                       check=False,
                                       exe_name="mfnwt.exe", )
    compare_heads()
    # compare_saturation()