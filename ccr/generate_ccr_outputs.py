import numpy as np
import flopy
import os
import matplotlib.pyplot as plt


def generate_outputs(NamFile,      # MODFLOW NAM filename
                     BudgetFile,   # Budget filename
                     UcnFile,      # Filename of UCN output from MT3DMS
                     base_dir='.'):# Base directory containing all input and output files

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
    flow_pond = df.loc['TOTAL_IN']['ZONE_2']  # SS influent rate to pond (L/day)
    flow_ups1 = df.loc['FROM_CONSTANT_HEAD']['ZONE_3']  # SS influent rate to top alluvium sublayer
    flow_ups2 = df.loc['FROM_CONSTANT_HEAD']['ZONE_4']  # SS influent rate to bottom alluvium sublayer subset

    return c_well, c_riv, 100.0*flow_pond/(flow_ups1 + flow_ups2), 100.0*flow_pond/flow_ups1

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
