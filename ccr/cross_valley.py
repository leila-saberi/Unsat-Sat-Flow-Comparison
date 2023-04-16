import ccr
import flopy
import numpy as np


def get_mf(A = 80.0,                 # Pond footprint area (acres)
           D = 100.0,                # Pond depth (feet)
           theta_x = 2.85,           # Slope in NS direction (degrees)
           theta_y = 4.0,            # Slope in EW direction (degrees)
           h_d = 10.0,               # Height top of dam reaches above pond surface (feet)
           b_d = 10.0,               # Width of top of dam (feet)
           NXL_u = 1.0,                # Number of pond lengths upstream of pond in flow domain
           NXL_d = 0.5,              # Number of pond lengths downstream of pond in flow domain
           NXW = 1.0,                # Number of pond widths before and after pond in flow domain
           num_el = 100000,          # Maximum number of flow cells
           ztop = 200.0,             # Elevation at top of pond (feet)
           up_depth = 0.0,           # Water-table depth below upstream ground surface (feet)
           down_depth = 0.0,         # Water-table depth below downstream ground surface (feet)
           modelname='cross_valley', # Prefix for all output files
           hk_ash = 1.0,             # Hydraulic conducitivity along rows of ash
           vka_ash = 0.01,           # Ratio of vertical to horizontal (along row) conductivity of ash
           sy_ash  = 0.1,            # Specific yield of ash
           ss_ash  = 1.e-4,          # Specific storage of ash
           hk_dam = 1.0,             # Hydraulic conducitivity along rows of dam
           vka_dam = 0.01,           # Ratio of vertical to horizontal (along row) conductivity of dam
           sy_dam  = 0.1,            # Specific yield of dam
           ss_dam  = 1.e-4,          # Specific storage of dam
           hk_all = 10.0,            # Hydraulic conducitivity along rows of alluvium
           vka_all = 0.1,            # Ratio of vertical to horizontal (along row) conductivity of alluvium
           sy_all = 0.1,             # Specific yield of alluvium
           ss_all = 1.e-4,           # Specific storage of alluvium
           hdry = -1000.0,           # Default output value for dry cells
           d_screen = 25.0,          # Depth the screened portion of well penetrates ground surface (feet)
           base_dir = '.'):          # Base directory to save all output files

    A *= 43560.0        # Convert to square feet
    
    tan_x = np.tan(theta_x*np.pi/180.0); cot_x = 1.0/tan_x
    tan_y = np.tan(theta_y*np.pi/180.0); cot_y = 1.0/tan_y
    
    cot_d = -cot_x + np.sqrt(A*tan_y/tan_x)/D; tan_d = 1.0/cot_d # Slope of dam
    L = D*(cot_d + cot_x) # Max length of dam (ft)
    W = 2.0*L*tan_x/tan_y # Max width of dam (ft)
    
    x_0 = NXL_u*L
    y_0 = (NXW + 0.5)*W
    
    # DEFINE MODEL EXTENT, GRID RESOLUTION, AND CHARACTERISTICS
    # Define flow domain and discretization:
    Lx = (NXL_u + NXL_d + 1.0)*L  # E-W flow-domain extent
    Ly = (2.0*NXW + 1.0)*W  # N-S flow-domain extent
    nlay = 4                     # 0: pond, 1: dam (pass through), 2: aquifer (w/well), 3: aquifer (w/o well)
    chi_T = Ly/Lx                # Aspect ratio of flow domain
    
    ncol = int(np.floor(np.sqrt(num_el/(chi_T*nlay))))
    nrow = int(np.floor(chi_T*ncol))
    
    if (nrow % 2) == 0:
        # Even number of rows, force to be odd
        nrow += 1
        ncol += 1 # Add additional column to maintain aspect ratio

    delr = Lx / ncol
    delc = Ly / nrow
    
    x = delr*(0.5 + np.arange(ncol)) # EW coordinate at center of element
    y = delc*(0.5 + np.arange(nrow)) # NS coordinate at center of element
    
    # Layers:
    # 1. Ash pond
    # 2. Dam (and pass through below pond)
    # 3. Alluvium (aquifer w/well)
    # 4. Alluvium (aquifer w/o well)
    
    top = ztop*np.ones((nrow,ncol))
    botm = np.zeros((nlay,nrow,ncol))
    ibound = np.zeros((nlay,nrow,ncol),dtype=np.int32) # Defaults as inactive
    ibound[nlay-2:,:,:] = 1 # All elements in alluvium layers are active
    
    x_3 = x_0 + L + h_d*cot_d + b_d
    # x_4 = (h_d - x_0*tan_x + x_3*tan_d)/(tan_d - tan_x)
    x_4 = x_3 + (D + h_d)*cot_d
    dy_4 = 4.0*(x_4-x_0)*tan_y; dy_4 = y_0
    y_min = y_0 - dy_4; y_max = y_0 + dy_4
    
    laytyp = 1       # Layer type (=1 means 'convertible', =0 means 'confined')

    # Default all hydraulic values to those of alluvium
    # Properties of pond and dam to be changed below
    hk = hk_all*np.ones((nlay,nrow,ncol))
    hk[0,:,:] = hk_ash
    hk[1,:,:] = hk_dam
    vka = vka_all*np.ones((nlay,nrow,ncol))
    vka[0,:,:] = vka_ash
    vka[1,:,:] = vka_dam
    sy = sy_all*np.ones((nlay,nrow,ncol))
    sy[0,:,:] = sy_ash
    sy[1,:,:] = sy_dam
    ss = ss_all*np.ones((nlay,nrow,ncol))
    ss[0,:,:] = ss_ash
    ss[1,:,:] = ss_dam
    
    tmin = 0.5 # Minimum cell thickness (ft)
    pond_x = []; pond_y = [] # Array of pond cell indices
    x_1 = x_0 + L
    z_1 = ztop - (x_1 - x_0)*tan_x + abs(y - y_0)*tan_y
    z_1 = ztop - D + abs(y - y_0)*tan_y
    z_end = ztop - D
    m_1 = (z_end - z_1)/(x[-1] - x_1)  # EW slope of ground surface after dam
    # Loop through rows:
    
    x_d0 = x_0 + D*cot_x 
    well_cells = [int(nrow/2),np.where(x > x_4)[0][0]]
    
    for i in np.arange(nrow):
        # Loop through columns:
        for j in np.arange(ncol):
            #
            # Default top array as ground surface:
            top[i,j] = ztop - (x[j] - x_0)*tan_x + abs(y[i] - y_0)*tan_y
            botm[nlay-2,i,j] = ztop - (x[j] - x_0)*tan_x - d_screen
            botm[nlay-1,i,j] = ztop - (x[j] - x_0)*tan_x - D_g
            if x[j] > x_0 + D*cot_x:
                # Ground surface no longer slopes along EW direction:
                top[i,j] = ztop - D + abs(y[i] - y_0)*tan_y
                botm[nlay-2,i,j] = ztop - D - d_screen
                botm[nlay-1,i,j] = ztop - D - D_g
            if x[j] > x_1:
                # Taper down to constant elevation at river boundary:
                top[i,j] = z_1[i] + m_1[i]*(x[j] - x_1)
                botm[nlay-2,i,j] = ztop - D - d_screen
                botm[nlay-1,i,j] = ztop - D - D_g
            botm[0,i,j] = top[i,j] # Set bottom of (inactive) pond layer
            botm[1,i,j] = top[i,j] # Set bottom of (inactive) dam layer
            # Find rectangle of cells that can contain pond or dam:
            if (((x[j] > x_0) and (x[j] < x_4)) and
                ((y[i] > y_min) and (y[i] < y_max))):
                # Inside footprint of pond + dam
                pond_surf = ztop # Pond surface elevation
                # Dam surface elevation:
                up_dam_surf = (x[j] - (x_0 + L))*tan_d
                dam_surf = ztop + min([up_dam_surf,
                                       h_d,
                                       h_d - (x[j] - x_3)*tan_d])
                if ((pond_surf > top[i,j]+tmin) and 
                    (pond_surf > ztop + up_dam_surf+tmin)):
                    # Inside of pond footprint:
                    ibound[0,i,j] = 1 # Pond cell is active
                    ibound[1,i,j] = 1 # Dam cell is active
                    if dam_surf < top[i,j]:
                        # Pond is above ground surface:
                        # Bottom of pond surface is ground surface:
                        botm[0,i,j] = top[i,j] 
                        # Bottom of "dam" layer is just below pond:
                        botm[1,i,j] = botm[0,i,j] - tmin
                        # Change properties of pass-through layers:
                        hk[1,i,j] = hk_all
                        vka[1,i,j] = vka_all
                        sy[1,i,j] = sy_all
                        ss[1,i,j] = ss_all
                    else:
                        # Pond is above dam surface
                        # Bottom of pond surface is dam surface:
                        botm[0,i,j] = dam_surf
                        botm[1,i,j] = top[i,j]
                    # "Top" surface is pond surface:
                    top[i,j] = pond_surf
                    # Collect pond indices:
                    pond_x = np.append(pond_x,j); pond_y = np.append(pond_y,i)
                elif dam_surf > top[i,j] + tmin:
                    # Inside footprint where dam is highest surface:
                    ibound[1,i,j] = 1 # Dam cell is active
                    top[i,j] = dam_surf    # New top surface is replaced with dam surface
                    botm[0,i,j] = top[i,j] # Bottom of (inactive) pond cell is dam surface

    # Initiate all cells as saturated:
    strt = np.zeros((nlay,nrow,ncol))
    for ilay in np.arange(nlay):
        strt[ilay,:,:] = top

    # Calculate upstream and downstream head values
    up_head = min(top[:,0]) - up_depth       # Inputted water-table depth below "v"
    down_head = min(top[:,-1]) - down_depth  # Inputted water-table depth below constant ground surf

    # Upstream boundary recarray:
    up_bound = np.recarray(2*nrow,formats="i8,i8,i8,f8,f8",names="lay,row,col,shead,ehead")
    up_bound["lay"] = np.append(np.tile(nlay-2,nrow),np.tile(nlay-1,nrow))  # Choosing alluvium layers to place CHD
    up_bound["col"] = 0             # First column (upstream boundary)
    up_bound["row"] = np.tile(np.arange(nrow),2) # All rows on this boundary (repeated for both alluvium layers)
    up_bound["shead"] = up_bound["ehead"] = up_head # Set boundary head values

    # Downstream boundary recarray
    down_bound = up_bound.copy() # Make copy of upstream boundary recarray
    down_bound["col"] = ncol - 1 # Replace column with last column
    down_bound["shead"] = down_bound["ehead"] = down_head # Replace boundary head value

    # Create boundary dictionary:
    chd_spd = {0:np.append(up_bound,down_bound).tolist()}

    mf = flopy.modflow.Modflow(modelname,            # Provide a name for the model
                               version = 'mfnwt',    # Pick version of MODFLOW (MODFLOW-NWT)
                               exe_name = 'mfnwt',   # Pick the executable (MODFLOW-NWT)
                               model_ws = base_dir)  # Install all files in specified directory

    # Discretization package object:
    dis = flopy.modflow.ModflowDis(mf,              # MODFLOW object
                                   nlay,            # Number of layers
                                   nrow,            # Number of rows
                                   ncol,            # Number of columns
                                   delr=delr,       # Row step size
                                   delc=delc,       # Column step size
                                   top=top,         # Top elevation of model
                                   botm=botm,       # Bottom elevations of each layer
                                   nper=1,          # Number of stress periods
                                   perlen=1.0e8,    # Length of each stress period
                                   nstp=1,          # Number of time steps
                                   steady=True)     # Array of steady/transient flags per stress period


    # UPW Object:
    upw = flopy.modflow.ModflowUpw(mf,                   # MODFLOW obje4t
                                   hk=hk,                # Hydraulic conductivity along rows (and columns, in this case)
                                   vka=vka,              # Vertical-to-horizontal conductivity ratio
                                   sy=sy,                # Specific yield
                                   ss=ss,                # Specific storage
                                   laytyp=laytyp,        # Layer type (Confined/convertible)
                                   hdry=hdry,            # Default head value for dry cells
                                   ipakcb=53)            # Flag (nonzero) dictating cell by cell budget data should
                                                         # be saved

    mf,bas,chd,nwt,oc = ccr.get_mf(mf,dis,upw,ibound,strt,chd_spd,well_cells,base_dir=base_dir)

    return mf,dis,bas,chd,upw,nwt,oc,well_cells

def mf_model_output(modelname='cross_valley', base_dir='.', run_model=False):

    dft,inp = ccr.mf_model_input(modelname, base_dir=base_dir)

    mf,dis,bas,chd,upw,nwt,oc,w_cells = get_mf(A = dft[0],                 
                                               D = dft[1],                
                                               theta_x = dft[2],           
                                               theta_y = dft[3],            
                                               D_g = dft[4],               
                                               h_d = dft[5],              
                                               b_d = dft[6],              
                                               NXL_u = dft[7],
                                               NXL_d = dft[8],
                                               NXW = dft[9],
                                               num_el = int(dft[10]),
                                               ztop = dft[11],
                                               up_depth = dft[12],
                                               down_depth = dft[13],
                                               modelname=modelname,
                                               hk_all = inp[0],
                                               vka_all = inp[1]/inp[0],
                                               sy_all = inp[2],
                                               ss_all = dft[16],
                                               hk_ash = inp[3],
                                               vka_ash = inp[4]/inp[3],
                                               sy_ash = inp[5],
                                               ss_ash = dft[14],
                                               hk_dam = inp[6],
                                               vka_dam = inp[7]/inp[6],
                                               sy_dam = inp[8],
                                               ss_dam = dft[15],
                                               hdry = dft[17],
                                               d_screen = dft[18],
                                               base_dir=base_dir)

    if run_model:
        success, mfoutput, flows = ccr.mf_model_output(mf,modelname,bas,w_cells,base_dir=base_dir)
    
    return mf, dis, bas, chd, upw, nwt, oc, w_cells

def model_output(modelname='cross_valley', base_dir='.', run_mp = True):

    # modelname=os.path.join(base_dir,modelname)
    mf, dis, bas, chd, upw, nwt, oc, w_cells = mf_model_output(modelname,base_dir=base_dir)
    conc_riv, conc_well = ccr.model_output(modelname, mf, dis, upw, bas, w_cells, base_dir=base_dir,
                                           run_mp = run_mp)
    #print("Well Relative Concentration: " + str(conc_well))
    #print("River Relative Concentration: " + str(conc_riv))

    return conc_riv, conc_well

if __name__ == 'main':
    model_output()
    print("done")