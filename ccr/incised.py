import numpy as np
import flopy
import copy
import ccr


def get_mf(A = 80.0,             # Pond footprint area (acres)
           D = 100.0,            # Pond depth (feet)
           chi_b = 0.5,          # Aspect ratio (W/L) of pond
           theta_p = 45.0,       # Pond taper slope (degrees)
           D_g = 400.0,          # Depth of alluvium (feet)
           NXL_u = 1.0,          # Number of pond lengths before pond in flow domain
           NXL_d = 1.0,          # Number of pond lengths after pond in flow domain
           NXW = 1.0,            # Number of pond widths before and after pond in flow domain
           num_el = 100000,      # Maximum number of flow cells
           ztop = 400.0,         # Top ground elevation (feet)
           up_depth = 15.0,      # Head value at upstream boundary (feet)
           down_depth = 25.0,    # Head value at downstream boundary (feet)
           modelname='incised',  # Prefix for all output files
           hk_ash = 1.0,         # Hydraulic conductivity along rows of ash
           vka_ash = 0.01,       # Ratio of vertical to horizontal (along row) conductivity of ash
           sy_ash  = 0.1,        # Specific yield of ash
           ss_ash  = 1.e-4,      # Specific storage of ash
           hk_all = 1000.0,      # Hydraulic conductivity along rows of alluvium
           vka_all = 0.1,        # Ratio of vertical to horizontal (along row) conductivity of alluvium
           sy_all = 0.1,         # Specific yield of alluvium
           ss_all = 1.e-4,       # Specific storage of alluvium
           hdry = -1000.0,       # Default output value for dry cells
           d_screen = 150.0,     # Depth of well screen
           nlay = 3,             # Number of layers in model
           base_dir = '.'):

    A *= 43560.0        # Convert to square feet
    L_p = np.sqrt(A/chi_b) # Length of pond (in feet)
    W_p = chi_b*L_p     # Width of pond (in feet)
    
    tan_p = np.tan(theta_p*np.pi/180.0)
    
    L_s = D/tan_p       # Length of the sloped portion of pond (meters)
    
    # DEFINE MODEL EXTENT, GRID RESOLUTION, AND CHARACTERISTICS
    # Define flow domain and discretization:
    Lx = (NXL_u + NXL_d + 1.0)*L_p # E-W flow-domain extent
    Ly = (2.0*NXW + 1.0)*W_p     # N-S flow-domain extent 
    # nlay = 3                     # 0: pond, 1: aquifer (w/well), 2: aquifer (w/o well)
    chi_T = Ly/Lx                # Aspect ratio of flow domain
    
    ncol = int(np.floor(np.sqrt(num_el/(chi_T*nlay))))
    nrow = int(np.floor(chi_T*ncol))

    if (nrow % 2) == 0:
        # Even number of rows, force to be odd
        nrow += 1
        ncol += 1 # Add additional column to maintain aspect ratio

    # Discretization:
    delr = Lx / ncol
    delc = Ly / nrow
    
    # Location arrays:
    x = delr*(0.5 + np.arange(ncol)) # EW coordinate at center of element
    y = delc*(0.5 + np.arange(nrow)) # NS coordinate at center of element
    
    NXL1 = NXL_u+1
    NXW1 = NXW+1
    x_0 = NXL_u*L_p; x_1 = NXL1*L_p
    y_0 = NXW*W_p; y_1 = NXW1*W_p
    
    # Layers:
    # 1. Ash pond
    # 2. Alluvium (aquifer w/well)
    # 3. Alluvium (aquifer w/o well)
    
    top = ztop*np.ones((nrow,ncol))
    botm = np.zeros((nlay,nrow,ncol))
    botm[0,:,:] = top
    botm[1,:,:] = ztop - d_screen
    botm[2,:,:] = ztop - D_g
    ibound = np.zeros((nlay,nrow,ncol),dtype=int) # Defaults as inactive
    ibound[nlay-2:,:,:] = 1 # All elements in bottom layers are active
    min_t = 0.5

    pond_cells = (np.where((x > x_0) & (x < x_1))[0], # Columns in pond footprint
                  np.where((y > y_0) & (y < y_1))[0]) # Rows in pond footprint

    well_cells = [int(nrow/2),np.where(x > x_1)[0][0]]
    # Loop through rows in pond footprint:
    for j in pond_cells[0]:
        # Loop through columns in pond footprint:
        for i in pond_cells[1]:
            ibound[0,i,j] = 1 # Pond cells are active
            botm[0,i,j] = max([ztop - D,                    # Bottom surface of pond
                               ztop - tan_p*(x[j] - x_0),   # Western slope
                               ztop + tan_p*(x[j] - x_1),   # Eastern slope
                               ztop - tan_p*(y[i] - y_0),   # Northern slope
                               ztop + tan_p*(y[i] - y_1)])  # Southern slope
            if top[i,j] - botm[0,i,j] < min_t:
                # Cell is too thin
                ibound[0,i,j] = 0       # Inactivate cell
                botm[0,i,j] = top[i,j]  # Reset bottom of this cell to be ground surface
                # min_t = top[i,j] - botm[0,i,j]

    # Set initial heads to make domain fully saturated:
    strt = ztop*np.ones((nlay,nrow,ncol))

    # Constant head boundaries:
    # Calculate upstream and downstream head values
    up_head = ztop - up_depth      # Inputted water-table depth below "v"
    down_head = ztop - down_depth  # Inputted water-table depth below constant ground surf

    # Upstream boundary recarray:
    up_bound = np.recarray(2*nrow,formats="i8,i8,i8,f8,f8",names="lay,row,col,shead,ehead")
    up_bound["lay"] = np.append(np.tile(nlay-2,nrow),
                                np.tile(nlay-1,nrow)) # Choosing the alluvium layer to place CHD
    up_bound["col"] = 0 # First column (upstream boundary)
    up_bound["row"] = np.tile(np.arange(nrow),2) # All rows on this boundary
    up_bound["shead"] = up_bound["ehead"] = up_head # Set boundary head values

    # Downstream boundary recarray
    down_bound = up_bound.copy() # Make copy of upstream boundary recarray
    down_bound["col"] = ncol - 1 # Replace column with last column
    down_bound["shead"] = down_bound["ehead"] = down_head # Replace boundary head value
    # Create boundary dictionary:
    chd_spd = {0:np.append(up_bound,down_bound).tolist()}

    # Layer Property Module:
    laytyp = 1
    hk = hk_ash*np.ones((nlay,nrow,ncol))
    hk[1:nlay-1,:,:] = hk_all
    vka = vka_ash*np.ones((nlay,nrow,ncol))
    vka[1:nlay-1,:,:] = vka_all
    sy = sy_ash*np.ones((nlay,nrow,ncol))
    sy[1:nlay-1,:,:] = sy_all
    ss = ss_ash*np.ones((nlay,nrow,ncol))
    ss[1:nlay-1,:,:] = ss_all

    # MODFLOW object:
    mf = flopy.modflow.Modflow(modelname,           # Provide a name for the model
                               exe_name = 'mfnwt',  # Pick the executable (MODFLOW-NWT)
                               version  = 'mfnwt',  # Choose model version as NWT
                               model_ws = base_dir) # Model working directory

    # Discretization package object:
    dis = flopy.modflow.ModflowDis(mf,              # MODFLOW object
                                   nlay,            # Number of layers
                                   nrow,            # Number of rows
                                   ncol,            # Number of columns
                                   delr=delr,       # Row step size
                                   delc=delc,       # Column step size
                                   top=ztop,        # Top elevation of model
                                   botm=botm,       # Bottom elevations of each layer
                                   nper=1,          # Number of stress periods
                                   perlen=1.0e8,    # Length of each stress period
                                   nstp=1,          # Number of time steps
                                   steady=True)     # Array of steady/transient flags per stress period

    # UPW Object:
    upw = flopy.modflow.ModflowUpw(mf,                   # MODFLOW object
                                   hk=hk,                # Hydraulic conductivity along rows (and columns, in this case)
                                   vka=vka,              # Vertical-to-horizontal conductivity ratio
                                   sy=sy,                # Specific yield
                                   ss=ss,                # Specific storage
                                   laytyp=laytyp,        # Layer type (Confined/convertible)
                                   iphdry=1,             # Flag to set hdry to default value (when > 0) 
                                   hdry=hdry,            # Default head value for dry cells
                                   ipakcb=53)            # Flag (nonzero) dictating cell by cell budget data should
                                                         # be saved

    mf, bas, chd, nwt, oc = ccr.get_mf(mf,dis,upw,ibound,strt,chd_spd,well_cells,base_dir=base_dir)

    return mf, dis, bas, chd, upw, nwt, oc, well_cells


def mf_model_output(modelname='incised', base_dir='.', run_model=False):

    dft,inp = ccr.mf_model_input(modelname,base_dir=base_dir)

    mf,dis,bas,chd,upw,nwt,oc,w_cells = get_mf(modelname=modelname,
                                               A = dft[0],
                                               D = dft[1],
                                               chi_b = dft[2],
                                               theta_p = dft[3],
                                               D_g = dft[4],
                                               NXL_u = dft[5],
                                               NXL_d = dft[6],
                                               NXW = dft[7],
                                               num_el = int(dft[8]),
                                               ztop = dft[9],
                                               up_depth = dft[10],
                                               down_depth = dft[11],
                                               hk_all = inp[0],
                                               vka_all = inp[1]/inp[0],
                                               sy_all = inp[2],
                                               ss_all = dft[13],
                                               hk_ash = inp[3],
                                               vka_ash = inp[4]/inp[3],
                                               sy_ash = inp[5],
                                               ss_ash = dft[12],
                                               hdry = dft[14],
                                               d_screen = dft[15],
                                               nlay = int(dft[16]),
                                               base_dir = base_dir)
    
    if run_model:
        success,mfoutput,flows = ccr.mf_model_output(mf,modelname,bas,w_cells,base_dir=base_dir)
    
    return mf, dis, bas, chd, upw, nwt, oc, w_cells


def model_output(modelname='incised', base_dir='.', run_mp = True):

    mf, dis, bas, chd, upw, nwt, oc, w_cells = mf_model_output(modelname, base_dir=base_dir)
    conc_riv, conc_well = ccr.model_output(modelname, mf, dis, upw, bas, w_cells, base_dir=base_dir,
                                           run_mp=run_mp)

    return conc_riv, conc_well