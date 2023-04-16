import incised,cross_valley,incised_vertical

import pandas as pd
import os,platform
import geopandas as gpd
from shapely.geometry import Point
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

def pathline_to_df(fn='mppth', xll=0, yll=0,write_shp=True):
    model_ws = './'
    columns = ['PID','Group','TimePointIndex','CumulativeTimeStep','TrackingTime','GlobX','GlobY','GlobZ','Layer',
               'Row','Column','Grid','LocalX','LocalY','LocalZ','LineSegmentIndex']
    df = pd.read_csv(os.path.join(model_ws,fn),skiprows=3,delim_whitespace=True, names=columns)

    df['UTMx'] = xll + df['GlobX']
    df['UTMy'] = (yll) + df['GlobY']
    df['Year'] = df['TrackingTime'] / (365.25)

    df['geometry'] = df.apply(lambda i: Point(i['UTMx'],i['UTMy']),axis=1)
    gdf = gpd.GeoDataFrame(df,geometry='geometry')
    if write_shp:
        if not os.path.exists(os.path.join(model_ws,'shp')):
            os.mkdir(os.path.join(model_ws,'shp'))
        gdf.to_file(os.path.join(model_ws,'shp','pth_pts_shp.shp'))
    return gdf

# from numpy import *
import numpy as np
from random import uniform
import flopy
import copy

def get_mf(mf, dis, upw, ibound, strt, chd_spd, well_cells, base_dir = '.'):
    """

    :param mf:
    :param dis:
    :param upw:
    :param ibound:
    :param strt:
    :param chd_spd:
    :param well_cells:
    :param base_dir:
    :return:
    """

    # Basic package object:
    bas = flopy.modflow.ModflowBas(mf,              # MODFLOW object
                                   ibound=ibound,   # Boundary flags (0 = inactive)
                                   strt=strt)       # Initial heads
    
    # CHD Object:
    chd = flopy.modflow.ModflowChd(mf,
                                   stress_period_data=chd_spd)
    
    
    # NWT solver object:
    nwt = flopy.modflow.ModflowNwt(mf)              # Chooses NWT matrix solver
    
    # OUTPUT CONTROL
    
    stress_period_data = {}
    stress_period_data[(0,0)] = ['save head',
                                 'save drawdown',
                                 'save budget',
                                 'print head',
                                 'print budget']
    #
    # Output control object
    oc = flopy.modflow.ModflowOc(mf,                # MODFLOW object
                                 stress_period_data=stress_period_data,
                                 compact=True)

    # LMT Linkage with MT3DMS for multi-species mass transport modeling
    output_file_name = 'mt3d_link.ftl'
    lmt = flopy.modflow.ModflowLmt(mf, output_file_name=output_file_name)

    return mf,bas,chd,nwt,oc

def get_mt(modelname,mf,dis,upw,bas,pond_cells,base_dir='.'):

    # Build MT3DMS object:
    mt = flopy.mt3d.Mt3dms(modelname=modelname,
                           version='mt3dms', 
                           exe_name='mt3dms.exe', 
                           modflowmodel=mf,
                           model_ws=base_dir)
    hdsFile = mf.name+'.hds'
    hds = flopy.utils.binaryfile.HeadFile(os.path.join(base_dir, hdsFile)).get_data(kstpkper=(0, 0))
    # Basic Transport Package (BTN):
    nrow = dis.nrow
    ncol = dis.ncol
    nlay = dis.nlay
    # wet_pond_cells = np.where(abs(hds[0,:,:] - mf.hdry) > 0.1)
    ifx = np.nonzero(abs(hds[pond_cells] - mf.hdry) > 0.1)
    ibound = bas.ibound
    icbund = np.ones((nlay, nrow, ncol), dtype=np.int32)
    icbund[pond_cells] = -1
    # icbund[0,:,:] = -1
    # ibnd_nz = np.nonzero(ibound[0,:,:])
    #print(nrow)
    #print(ibnd_nz[0].tolist())
    #print(ncol)
    #print(ibnd_nz[1].tolist())
    #icbund[0, ibnd_nz] = -1
    # icbund[0, ibnd_nz[0], ibnd_nz[1]] = -1  # Constant concentration
    sconc = np.zeros((nlay, nrow, ncol))
    #sconc[0, ibnd_nz[0], ibnd_nz[1]] = 1.0
    # sconc[0, ibnd_nz] = 1.0
    # sconc[0,:,:] = 1.0
    # sconc[0, wet_pond_cells] = 1.0
    # sconc[0, wet_pond_cells[0], wet_pond_cells[1]] = 1.0
    sconc[pond_cells] = 1.0
    btn = flopy.mt3d.Mt3dBtn(mt,sconc=0.0, prsity=upw.sy, thkmin=0.1, munit='g',
                                icbund=icbund, ttsmult=2.0, ttsmax=1.0e7)

    # Advection Package (ADV)
    mixelm = 0  # Standard finite-difference method
    percel = 0.75  # Courant number PERCEL is also a stability constraint
    adv = flopy.mt3d.Mt3dAdv(mt,mixelm=mixelm,percel=percel)

    # GCG Solver
    mxiter = 1  # Maximum number of outer iterations
    iter1 = 500  # Maximum number of inner iterations
    isolve = 3  # Preconditioner = Modified Incomplete Cholesky
    gcg = flopy.mt3d.Mt3dGcg(mt, mxiter=mxiter, iter1=iter1, isolve=isolve, cclose=1e-04)

    # DSP file
    al = 100.0  # longitudinal dispersivity
    dmcoef = 0  # effective molecular diffusion coefficient
    trpt = 0.1  # ratio of the horizontal transverse dispersivity to the longitudinal dispersivity
    trpv = 0.0  # ratio of the vertical transverse dispersivity to the longitudinal dispersivity

    dsp = flopy.mt3d.Mt3dDsp(mt, al=al, dmcoef=dmcoef, trpt=trpt, trpv=trpv)

    #ssm = None
    # Source and Sink Matching Package (SSM)
    itype = flopy.mt3d.Mt3dSsm.itype_dict()

    #[K,I,J,CSS,iSSType] = layer, row, column, source concentration, type of sink/source: well-constant concentration cell
    # num_ssm = len(np.nonzero(ibound[0,:,:])[0])
    # num_ssm = len(wet_pond_cells[0])
    num_ssm = len(ifx[0])
    ssm_data = np.recarray(num_ssm,formats="i4,i4,i4,f4,i4",names="k,i,j,css,itype")

    # ssm_data["k"] = 0                              # Choosing the pond layer
    # ssm_data["i"] = wet_pond_cells[0]              # All active rows in pond
    # ssm_data["j"] = wet_pond_cells[1]              # All active columns in pond
    ssm_data["k"] = pond_cells[0][ifx]               # Choosing the pond layer
    ssm_data["i"] = pond_cells[1][ifx]               # All active rows in pond
    ssm_data["j"] = pond_cells[2][ifx]               # All active columns in pond
    ssm_data["css"] = 1.0                            # Set concentration value
    ssm_data["itype"] = -1                           # Make concentration value fixed
    ssm_spd = {0:ssm_data}

    ssm = flopy.mt3d.Mt3dSsm(mt,stress_period_data=ssm_spd)

    return mt, btn, adv, gcg, ssm

def mp_write_loc(modelname,Np_tot,loc_data):

    # Starting locations file:
    slf=open(modelname+".loc", 'w')
    slf.write("# starting locations file\n")
    slf.write("1\n")
    slf.write("1\n") #one group
    slf.write("group1\n")
    slf.write(str(Np_tot)+"\n")
    slf.write(loc_data)
    slf.close()

    return None

def get_mp(modelname, mf, dis, upw, bas, flows, C_i = 1.0, drain_rate = 1.0, N_p = 10000, NXp = 2, NYp = 2, num_pV = 10.0, base_dir = '.'):
    """

    :param modelname:
    :param mf:
    :param dis:
    :param upw:
    :param bas:
    :param flows: # Flow rates from pond and through well
    :param C_i: # Initial concentration (mg/L)
    :param drain_rate: # Number of ponds to empty per particle drop
    :param N_p:  # Number of particles to drop each time (approximate)
    :param NXp: # Number of particles per cell in x direction
    :param NYp: # Number of particles per cell in y direction
    :param num_pV: # Number of pore volumes to flush before terminating simulation
    :param base_dir: # Save all files in prescribed directory
    :return:
    """
    prsity = upw.sy.array
    disFile = modelname+".dis"# os.path.join(base_dir,modelname+".dis")
    hdsFile = modelname+".hds"# os.path.join(base_dir,modelname+".hds")
    cbbFile = modelname+".cbc"#os.path.join(base_dir,modelname+".cbc")

    # Starting locations:
    # Import head values:
    hds = flopy.utils.binaryfile.HeadFile(os.path.join(base_dir,hdsFile)).get_data(kstpkper=(0,0))

    # Import botm values:
    botm = dis.botm.array; top = dis.top.array
    # Calculate particle starting locations:
    loc_data = ""         # Output string (to be written in .loc file)
    dx_arr = (1.0 + np.arange(NXp))/(NXp + 1.0)
    dy_arr = (1.0 + np.arange(NYp))/(NYp + 1.0)

    pond_rows = np.unique(np.nonzero(bas.ibound[0,:,:])[0])
    wet_check = hds[-2,pond_rows, 0] - botm[-2,pond_rows, 0] > 0.0 * (botm[-3,pond_rows,0] - botm[-2,pond_rows,0])
    not_wet_check = np.logical_not(wet_check)
    hw = np.zeros(len(pond_rows))
    hw[wet_check] = (hds[-2,:,0] - botm[-1,:, 0])[pond_rows][wet_check]
    hw[not_wet_check] = (hds[-1,:,0] - botm[-1,:,0])[pond_rows][not_wet_check]
    zstep = NXp*NYp*hw/np.round(N_p*hw/sum(hw))
    zbot = botm[-2, pond_rows, 0] - botm[-1, pond_rows, 0]
    ztop = botm[-3, pond_rows, 0] - botm[-2, pond_rows, 0]

    particleID = 0      # Particle counter

    for i_pond, pond_row in enumerate(pond_rows):
        # Place particles only in non-dry cells:
        # if i_pond % 2:
        #     # Ignore every other row
        #     continue
        for dx in dx_arr:
            for dy in dy_arr:
                # Uniform distribution in the vertical direction:
                z_arr = (np.arange(0.5*zstep[i_pond], hw[i_pond], zstep[i_pond]) +
                         0.5*zstep[i_pond]*uniform(-1.0,1.0))
                for z in z_arr:
                    ilay = dis.nlay if z <= zbot[i_pond] else dis.nlay-1
                    dz = z/zbot[i_pond] if z <= zbot[i_pond] else (z - zbot[i_pond])/ztop[i_pond]
                    if z == z_arr[-1]:
                        dz = (hw[i_pond] - zbot[i_pond])/ztop[i_pond]
                    particleID += 1
                    loc_data+=("{} {} {} {} {} {} {} {} {} {} {}\n"
                               .format(
                               particleID, # Particle ID (index)
                               1,          # Group number
                               1,          # Grid
                               ilay,       # Layer
                               pond_row+1, # Row
                               1,          # Column
                               dx,         # LocalX
                               dy,         # LocalY
                               dz,         # LocalZ
                               0.0,        # Release time
                               1))         # Label

    slf=open(os.path.join(base_dir,modelname+".loc"), 'w')
    slf.write("# starting locations file\n")
    slf.write("1\n")
    slf.write("1\n") #one group
    slf.write("group1\n")
    slf.write(str(particleID)+"\n")
    slf.write(loc_data)
    slf.close()

    # MODPATH (Particle Tracking):
    mp = flopy.modpath.Modpath6(modelname=modelname,
                                  exe_name='mp6',
                                  modflowmodel=mf,
                                  dis_file=disFile,
                                  head_file=hdsFile,
                                  budget_file=cbbFile,
                                  model_ws=base_dir,
                                  verbose=True,
                                  load=True,
                                  listunit=7)

    mpbas = flopy.modpath.Modpath6Bas(mp,               # MODPATH Object
                                     laytyp=upw.laytyp,# Layer type (1 - Convertible)
                                     prsity=prsity,    # Porosity
                                     ibound=bas.ibound,# ibound array
                                     hdry=upw.hdry)    # Default head value for dry cell

    mpsim = flopy.modpath.Modpath6Sim(mp,                      # MODPATH Object
                                     option_flags=[2,         # Simulation Type (1 - Endpoint simulation, 2 - Pathline simulation)
                                                   1,         # Tracking Direction
                                                   1,         # Weak Sink Option
                                                   1,         # Weak Source Option
                                                   1,         # Reference Time Option
                                                   2,         # Stop Option (2 - Track particles through to termination, 3 - Specify a time to stop particle tracking)
                                                   2,         # Particle Generation Option (2 - Start locations written in file)
                                                   2,         # Time Point Option (A specified number of time points are calculated for a fixed time increment)
                                                   1,         # Budget Output Option
                                                   1,         # Zone Array Option
                                                   1,         # Retardation Option
                                                   1])        # Advective Observations Option

    return mp, mpbas, mpsim

def get_final_concentrations(well_cells, base_dir='.'):
    filename = os.path.join(base_dir,'MT3D001.UCN')
    ucn = flopy.utils.binaryfile.UcnFile(filename)
    ucn_data = ucn.get_alldata(-2)[0]
    c_well = ucn_data[well_cells[0],well_cells[1]]
    if c_well > 1.0:
        # "Well cell" is dry, so calculate concentration from cell in lowest layer:
        c_well = ucn.get_alldata(-1)[0][well_cells[0],well_cells[1]]
    c_riv = ucn_data[well_cells[0], -1]
    if c_riv > 1.0:
        # "River cell" is dry, so calculate concentration from cell in lowest layer:
        c_riv = ucn.get_alldata(-1)[0][well_cells[0],-1]
    return c_well, c_riv

def write_inp(filename, input_data): # Name of file # Recarray of input data
    # Writes input file
    outStr=""
    for i in np.arange(len(input_data)):
        outStr+="{} {} {}\n".format(input_data["label"][i],
                                    input_data["value"][i],
                                    input_data["unit"] [i])
    slf=open(filename, 'w')
    slf.write(outStr)
    slf.close()
    return None

def mf_model_input(modelname, base_dir='.'):

    # This routine reads input parameters from a .inp file (with sensitivity
    # parameters) and a .dft file (with non-sensitivity parameters)
    # Import "default" parameters (those not in sensitivity study):
    print(f"Looking for {os.path.join(base_dir,modelname+'.dft')}")
    _,dft,_ = np.loadtxt(os.path.join(base_dir,modelname+'.dft'),
                         dtype = {'names': ('label','value','unit'),
                                  'formats':('S8','f8','S8')},
                         unpack=True,skiprows=1)

    # Import sensitivity parameters (K's and n's from different materials):
    _,inp,_ = np.loadtxt(os.path.join(base_dir,modelname+'.inp'),
                         dtype = {'names': ('label','value','unit'),
                                  'formats':('S8','f8','S8')},unpack=True)

    return dft,inp

def mf_model_output(mf, modelname, bas, well_cells, pond_cells, base_dir='.'):

    # Write model input files:
    mf.write_input()

    # Run the model:
    success, mfoutput = mf.run_model(silent=True, pause=False, report=True)

    # Check if it ran correctly:
    # if not success:
    #    raise Exception('MODFLOW did not terminate normally')

    # Calculate flows (Zone budget)
    cbb_file = os.path.join(base_dir,modelname + ".cbc")
    zon = np.ones(bas.ibound.shape,dtype=int)
    #upgradien zone
    zon[:, 38:77, 1] = 4
    # River zone:
    zon[:,38:77,-1] = 2
    # Pond zone:
    zon_rows,zon_cols = np.nonzero(bas.ibound[0,:,:])
    # zon[0,np.nonzero(bas.ibound[0,:,:])] = 3
    # zon[0,zon_rows,zon_cols] = 3
    zon[pond_cells] = 3
    # Collector well:
    # zon[-1,well_cells] = 4
    zon[-3,well_cells[0],well_cells[1]] = 4
    zb = flopy.utils.ZoneBudget(cbb_file, zon, kstpkper=(0,0))
    df = zb.get_dataframes().reset_index()
    flow_riv = df['ZONE_2'][1]*28.3168  # SS rate into river (L/day)
    flow_UPG = df['ZONE_4'][1] * 28.3168  # SS rate into river (L/day)
    flow_pond = df['ZONE_1'][3]*28.3168 # SS effluent rate from pond (L/day)
    flow_well = 0.5*(df['ZONE_1'][5] + df['ZONE_4'][2])*28.3168 # SS flow through well (L/day)
    print(df.to_string())
    print('Flow into river: ' + str(flow_riv) + ' L/day.')
    print('Flow from pond: ' + str(flow_pond) + ' L/day.')
    print('Flow through UPG: ' + str(flow_UPG) + ' L/day.')
    print('Percent: ' + str(flow_pond/flow_riv*100.0))
    exit()

    flows = (flow_riv, flow_pond, flow_well)

    return success, mfoutput, flows

def mt3_model_output(modelname, mf, dis, upw, bas, pond_cells, base_dir='.'):

    # MT3DMS:
    mt,btn,adv,gcg,ssm = get_mt(modelname+'_mt3', mf, dis, upw, bas, pond_cells, base_dir=base_dir)

    # Write model input files:
    mt.write_input()

    # Run the model:
    success, mtoutput = mt.run_model(silent=True)

    return success, mtoutput

def mp_model_output(modelname, mf, dis, upw, bas, flows, base_dir='.'):

    # Import "default" parameters (those not in sensitivity study):
    _,mpdft,_ = np.loadtxt(os.path.join(base_dir,modelname+'.mpdft'),
                           dtype = {'names': ('label','value','unit'),
                                    'formats':('S8','f8','S8')},
                           unpack=True,skiprows=1)

    N_p = int(mpdft[0])
    NXp = int(mpdft[1])
    NYp = int(mpdft[2])
    C_i = mpdft[3]
    drain_rate = mpdft[4]
    num_pV = mpdft[5]
    num_c = int(mpdft[6])
    
    # MODPATH:
    mp,mpbas,mpsim = get_mp(modelname, mf, dis, upw, bas, flows, C_i, drain_rate, N_p, NXp=NXp, NYp=NYp,
                            num_pV=num_pV, base_dir=base_dir)

    # Write model input files:
    mp.write_input()

    # Run the model:
    mp.run_model(silent=True)

    return None # num_c,mp_array

def model_output(modelname, mf, dis, upw, bas, well_cells, pond_cells, base_dir = '.', run_mp = 'True'):
    """

    :param modelname:
    :param mf:
    :param dis:
    :param upw:
    :param bas:
    :param well_cells:
    :param base_dir:
    :param run_mp:
    :return: conc_riv, conc_well
    """
    # Generate MODFLOW outputs:
    success, mfoutput, flows = mf_model_output(mf, modelname, bas, well_cells, pond_cells, base_dir=base_dir)

    # Generate MT3DMS outputs:
    success, mtoutput = mt3_model_output(modelname, mf, dis, upw, bas, pond_cells, base_dir=base_dir)

    conc_well, conc_riv = get_final_concentrations(well_cells, base_dir)

    # # Generate MODPATH outputs:
    if run_mp:
        mp_model_output(modelname, mf, dis, upw, bas, flows, base_dir=base_dir)

    return conc_riv, conc_well, flows

def run_sensitivity_case(model_output, modelname, modelname_i, input_data, run_mp = True):
    #
    if platform.system() == 'Windows':
        cp_str = 'copy '; sl_str = '\\'
    else:
        cp_str = 'cp '; sl_str = '/'
    #
    os.system('mkdir ' + modelname_i)
    os.system(cp_str + modelname+'.dft ' + modelname_i + sl_str + modelname_i + '.dft')
    os.system(cp_str + modelname+'.mpdft ' + modelname_i + sl_str + modelname_i + '.mpdft')
    write_inp(modelname_i + sl_str + modelname_i + '.inp',input_data)
    #
    c_riv, c_well = model_output(modelname_i, base_dir=modelname_i, run_mp = run_mp)
    return c_riv, c_well

def sensitivity_study(model_output, modelname, run_mp = True):
    dtype = {'names': ('label','value','unit'),
             'formats':('U10','f8','U8')}
    input_data = np.loadtxt(modelname+'.inp',dtype=dtype)

    # filename = modelname + '.par'
    # par_data = np.loadtxt(filename,skiprows=2)
    # num_par = len(par_data[0][1:])
    # suffix = '0'*num_par
    # par_data = par_data[:,1:]
    # par_array = par_data[0,:]
    # input_data["value"][0:num_par] = par_array
    # suffix = '0'
    # modelname_i = modelname+'_'+suffix
    # c_riv_0, c_well_0 = run_sensitivity_case(model_output,
    #                                          modelname,
    #                                          modelname_i,
    #                                          input_data,
    #                                          run_mp = run_mp)
    # c_riv_0 *= 100.0; c_well_0 *= 100.0
    # print('Base Case:')
    # print('Well Relative Concentration: ' + str(c_well_0) + '%')
    # print('River Relative Concentration: ' + str(c_riv_0) + '%')
    # print('\n')
    strt_log = -1
    end_log = 2
    num_log = end_log - strt_log + 1

    for i_mult, mult in enumerate(np.logspace(strt_log, end_log, num_log), start=1):
        input_data_i = copy.deepcopy(input_data)
        input_data_i['value'][np.where(input_data['label'] == 'kh_all')[0][0]] = mult*(
                input_data['value'][np.where(input_data['label'] == 'kh_all')[0][0]])
        input_data_i['value'][np.where(input_data['label'] == 'kv_all')[0][0]] = mult*(
                input_data['value'][np.where(input_data['label'] == 'kv_all')[0][0]])
        suffix = str(i_mult)
        modelname_i = modelname+'_'+suffix
        c_riv,c_well = run_sensitivity_case(model_output,
                                            modelname,
                                            modelname_i,
                                            input_data_i,
                                            run_mp = run_mp)

        _, _, flow_pond, flow_pond1 = generate_outputs(modelname_i+'.nam',  # MODFLOW NAM filename
                                           modelname_i+'.cbc',  # Budget filename
                                           'MT3D001.UCN',  # Filename of UCN output from MT3DMS
                                           base_dir=modelname_i)  # Base directory containing all input and output files

        print('Multiplier = ' + str(mult))
        print('Well Relative Concentration: ' + str(c_well*100.0) + ' %')
        print('River Relative Concentration: ' + str(c_riv*100.0) + ' %')
        print('Flow through pond: ' + str(flow_pond) + ' % of entire alluvium')
        print('Flow through pond: ' + str(flow_pond1) + ' % of top alluvium sublayer')
        print('\n')

        del input_data_i

    return None

def plot_pathline(modelname, base_dir='.'):
    mfnamFile = modelname+'.nam'
    mf = flopy.modflow.Modflow.load(mfnamFile,
                                    version = 'mfnwt',
                                    exe_name = 'mfnwt.exe',
                                    model_ws = base_dir)

    dis = mf.get_package('dis')

    grd = flopy.discretization.structuredgrid.StructuredGrid(delc = dis.delc.array,
                                                             delr = dis.delr.array,
                                                             top = dis.top.array,
                                                             botm = dis.botm.array,
                                                             nlay = dis.nlay,
                                                             nrow = dis.nrow,
                                                             ncol = dis.ncol)

    pmv = flopy.plot.map.PlotMapView(model=mf,
                                     modelgrid=grd,
                                     extent=(3500,5500,4300,5100))

    pmv.plot_grid(linewidths=0.1)
    pthFile = os.path.join(base_dir, modelname+'.mppth')
    pth = flopy.utils.PathlineFile(pthFile)
    print('collecting pathline data...')
    pth_data = pth.get_alldata()
    print('done!')
    pmv.plot_pathline(pth_data, layer ='all', colors='blue', linewidths=0.25)
    plt.xlabel('X (ft)')
    plt.ylabel('Y (ft)')
    plt.title('Particle pathlines')
    # print(dir(line_segments))
    plt.show()

    return None

def generate_outputs(NamFile, BudgetFile, UcnFile, base_dir='.'):
    """

    :param NamFile:  # MODFLOW NAM filename
    :param BudgetFile:   # Budget filename
    :param UcnFile: # Filename of UCN output from MT3DMS
    :param base_dir: # Base directory containing all input and output files
    :return:
    """
    mf = flopy.modflow.Modflow.load(NamFile,
                                    version='mfnwt',
                                    exe_name='mfnwt.exe',
                                    model_ws=base_dir)
    dis = mf.get_package('dis')
    bas = mf.get_package('bas6')

    # Calculate row and column indices of well cell:
    well_row = int(dis.nrow/2)
    well_col = np.nonzero(bas.ibound[-3,well_row,:])[0][-1]+1

    # Calculate concentrations in well and river cells:
    ucn = flopy.utils.binaryfile.UcnFile(os.path.join(base_dir,UcnFile))
    ucn_data = ucn.get_alldata(-2)[0]
    c_well = ucn_data[well_row, well_col]
    if c_well > 1.0:
        # "Well cell" is dry, so calculate concentration from cell in lowest layer:
        c_well = ucn.get_alldata(-1)[0][well_row, well_col]
    c_riv = ucn_data[well_row, -1]
    if c_riv > 1.0:
        # "River cell" is dry, so calculate concentration from cell in lowest layer:
        c_riv = ucn.get_alldata(-1)[0][well_row, -1]

    # Calculate flows (Zone budget)
    zon = np.ones(bas.ibound.shape,dtype=int)
    # Pond zone:
    zon_rows, zon_cols = np.nonzero(bas.ibound[0,:,:])
    zon[0,zon_rows,zon_cols] = 2
    # Upstream Regions:
    zon[-2,zon_rows,0] = 3
    zon[-1,zon_rows, 0] = 4
    zb = flopy.utils.ZoneBudget(os.path.join(base_dir, BudgetFile), zon, kstpkper=(0,0))
    df = zb.get_dataframes().reset_index()
    df.set_index("name", inplace=True)
    flow_ups1 = df.loc['FROM_CONSTANT_HEAD']['ZONE_3']  # SS influent rate to top alluvium sublayer
    flow_ups2 = df.loc['FROM_CONSTANT_HEAD']['ZONE_4'] # SS influent rate to bottom alluvium sublayer subset
    flow_pond = df.loc['TOTAL_IN']['ZONE_2']           # SS influent rate to pond (L/day)
    return c_well, c_riv, 100.0*flow_pond/(flow_ups1 + flow_ups2), 100.0*flow_pond/flow_ups1