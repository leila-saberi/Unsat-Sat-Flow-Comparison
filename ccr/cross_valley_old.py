import numpy as np
import ccr,flopy,copy

def get_mf(A = 80.0,                 # Pond footprint area (acres)
           D = 100.0,                # Pond depth (feet)
           theta_x = 4.0,            # Slope in NS direction (degrees)
           theta_y = 4.0,            # Slope in EW direction (degrees)
           D_g = 500.0,              # Depth of alluvium (feet)
           h_d = 10.0,               # Height top of dam reaches above pond surface (feet)
           b_d = 10.0,               # Width of top of dam (feet)
           NXL = 1.0,                # Number of pond lengths before and after pond in flow domain
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
           prsity_ash = 0.05,        # Porosity of ash (-)
           hk_dam = 1.0,             # Hydraulic conducitivity along rows of dam
           vka_dam = 0.01,           # Ratio of vertical to horizontal (along row) conductivity of dam
           sy_dam  = 0.1,            # Specific yield of dam
           ss_dam  = 1.e-4,          # Specific storage of dam
           prsity_dam = 0.01,        # Porosity of dam (-)
           hk_all = 10.0,            # Hydraulic conducitivity along rows of alluvium
           vka_all = 0.1,            # Ratio of vertical to horizontal (along row) conductivity of alluvium
           sy_all = 0.1,             # Specific yield of alluvium
           ss_all = 1.e-4,           # Specific storage of alluvium
           prsity_all = 0.15,        # Porosity of alluvium (-)
           chi_pond = 0.2,           # Fraction of pond thickness above lowest point that remains saturated
           hdry = -1000.0):          # Default output value for dry cells

    A *= 43560.0        # Convert to square feet
    
    tan_x = np.tan(theta_x*np.pi/180.0); cot_x = 1.0/tan_x
    tan_y = np.tan(theta_y*np.pi/180.0); cot_y = 1.0/tan_y
    
    cot_d = -cot_x + np.sqrt(A*tan_y/tan_x)/D; tan_d = 1.0/cot_d # Slope of dam
    L = D*(cot_d + cot_x) # Max length of dam (ft)
    W = 2.0*L*tan_x/tan_y # Max width of dam (ft)
    
    x_0 = NXL*L
    y_0 = (NXW + 0.5)*W
    
    # DEFINE MODEL EXTENT, GRID RESOLUTION, AND CHARACTERISTICS
    # Define flow domain and discretization:
    Lx = (2.0*NXL + 1.0)*L  # E-W flow-domain extent
    Ly = (2.0*NXW + 1.0)*W  # N-S flow-domain extent 
    nlay = 3                     # 0: pond, 1: dam (pass through), 2: aquifer
    chi_T = Ly/Lx                # Aspect ratio of flow domain
    
    ncol = int(np.floor(np.sqrt(num_el/(chi_T*nlay))))
    nrow = int(np.floor(chi_T*ncol))
    
    delr = Lx / ncol
    delc = Ly / nrow
    
    x = delr*(0.5 + np.arange(ncol)) # EW coordinate at center of element
    y = delc*(0.5 + np.arange(nrow)) # NS coordinate at center of element
    
    # Layers:
    # 1. Ash pond
    # 2. Dam (and pass through below pond)
    # 3. Alluvium (aquifer)
    
    top = ztop*np.ones((nrow,ncol))
    botm = np.zeros((nlay,nrow,ncol))
    ibound = np.zeros((nlay,nrow,ncol),dtype=np.int32) # Defaults as inactive
    ibound[2,:,:] = 1 # All elements in bottom layer are active
    
    x_3 = x_0 + L + h_d*cot_d + b_d
    # x_4 = (h_d - x_0*tan_x + x_3*tan_d)/(tan_d - tan_x)
    x_4 = x_3 + (D + h_d)*cot_d
    dy_4 = 4.0*(x_4-x_0)*tan_y; dy_4 = y_0
    y_min = y_0 - dy_4; y_max = y_0 + dy_4
    
    laytyp  =     1       # Layer type (=1 means 'convertible', =0 means 'confined')
    
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
    prsity = prsity_all*np.ones((nlay,nrow,ncol))
    prsity[0,:,:] = prsity_ash
    prsity[1,:,:] = prsity_dam
    
    tmin = 0.5 # Minimum cell thickness (ft)
    pond_x = []; pond_y = [] # Array of pond cell indices
    pond_head = ztop - (1.0 - chi_pond)*D 
    x_1 = x_0 + L
    z_1 = ztop - (x_1 - x_0)*tan_x + abs(y - y_0)*tan_y
    z_1 = ztop - D + abs(y - y_0)*tan_y
    z_end = ztop - D
    m_1 = (z_end - z_1)/(x[-1] - x_1)  # EW slope of ground surface after dam
#   print(z_4 + m_4*(x[-1] - x_4)); exit()
    # Loop through rows:
    for i in np.arange(nrow):
        # Loop through columns:
        for j in np.arange(ncol):
            #
            # Default top array as ground surface:
            top[i,j] = ztop - (x[j] - x_0)*tan_x + abs(y[i] - y_0)*tan_y
            if x[j] > x_0 + D*cot_x:
                # Ground surface no longer slopes along EW direction:
                top[i,j] = ztop - D + abs(y[i] - y_0)*tan_y
            if x[j] > x_1:
                # Taper down to constant elevation at river boundary:
                top[i,j] = z_1[i] + m_1[i]*(x[j] - x_1)
            botm[0,i,j] = top[i,j] # Set bottom of (inactive) pond layer
            botm[1,i,j] = top[i,j] # Set bottom of (inactive) dam layer
            # Set bottom of alluvium surface (never changes):
            botm[2,i,j] = top[i,j] - D_g         
            #
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
                        # Bottom of "dam" layer is just below pond bottom:
                        botm[1,i,j] = botm[0,i,j] - 1.0 
                        # Change properties of pass-through layer:
                        hk[1,i,j] = hk_all
                        vka[1,i,j] = vka_all
                        sy[1,i,j] = sy_all
                        ss[1,i,j] = ss_all
                        prsity[1,i,j] = prsity_all
                    else:
                        # Pond is above dam surface
                        # Bottom of pond surface is dam surface:
                        botm[0,i,j] = dam_surf
                        # Bottom of dam is ground surface:
                        botm[1,i,j] = top[i,j]
                    # "Top" surface is pond surface:
                    top[i,j] = pond_surf
                    # Collect pond indices:
                    pond_x = np.append(pond_x,j); pond_y = np.append(pond_y,i)
                elif dam_surf > top[i,j]+tmin:
                    # Inside footprint where dam is highest surface:
                    ibound[1,i,j] = 1 # Dam cell is active
                    botm[1,i,j] = top[i,j] # Bottom of dam cell is ground surface
                    top[i,j] = dam_surf # New top surface is replaced with dam surface
                    botm[0,i,j] = top[i,j] # Bottom of (inactive) pond cell is dam surface

    pond_cells = (np.unique(pond_x).astype(int),np.unique(pond_y).astype(int))

    # Initiate all cells as saturated:
    strt = np.zeros((nlay,nrow,ncol))
    strt[0,:,:] = top
    strt[1,:,:] = top
    strt[2,:,:] = top

    # Calculate upstream and downstream head values
    up_head = top[:,0] - up_depth            # Inputted water-table depth below "v"
    down_head = min(top[:,-1]) - down_depth  # Inputted water-table depth below constant ground surf

    # Upstream boundary recarray:
    up_bound = np.recarray(nrow,formats="i8,i8,i8,f8,f8",names="lay,row,col,shead,ehead")
    up_bound["lay"] = 2 # Choosing the alluvium layer to place CHD
    up_bound["col"] = 0 # First column (upstream boundary)
    up_bound["row"] = np.arange(nrow) # All rows on this boundary
    up_bound["shead"] = up_bound["ehead"] = up_head # Set boundary head values
    # Downstream boundary recarray
    down_bound = up_bound.copy() # Make copy of upstream boundary recarray
    down_bound["col"] = ncol - 1 # Replace column with last column
    down_bound["shead"] = down_bound["ehead"] = down_head # Replace boundary head value

    # Create boundary dictionary:
    chd_spd = {0:np.append(up_bound,down_bound).tolist()}

    mf = flopy.modflow.Modflow(modelname,           # Provide a name for the model
                               version = 'mfnwt',   # Pick version of MODFLOW (MODFLOW-NWT)
                               exe_name = 'mfnwt')  # Pick the executable (MODFLOW-NWT)

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
                                   perlen=1.0,      # Length of each stress period
                                   nstp=1,          # Number of time steps
                                   steady=True)     # Array of steady/transient flags per stress period


    # UPW Object:
    upw = flopy.modflow.ModflowUpw(mf,                   # MODFLOW object
                                   hk=hk,                # Hydraulic conductivity along rows (and columns, in this case)
                                   vka=vka,              # Vertical-to-horizontal conductivity ratio
                                   sy=sy,                # Specific yield
                                   ss=ss,                # Specific storage
                                   laytyp=laytyp,        # Layer type (Confined/convertible)
                                   hdry=hdry,            # Default head value for dry cells
                                   ipakcb=53)            # Flag (nonzero) dictating cell by cell budget data should
                                                         # be saved

    mf,bas,chd,nwt,oc = ccr.get_mf(mf,dis,upw,ibound,strt,chd_spd,pond_cells)

    return mf,dis,bas,chd,upw,nwt,oc,pond_cells

def model_output(modelname):

    _,inp,_ = np.loadtxt(modelname+'.inp', 
                         dtype = {'names': ('label','value','unit'),
                                  'formats':('S8','f8','S8')},unpack=True)
    
    ztop = 200.0                
    hk_all = inp[0]
    vka_all = inp[1]/inp[0]
    prsity_all = inp[2]
    hk_ash = inp[3]
    vka_ash = inp[4]/inp[3]
    prsity_ash = inp[5]
    hk_dam = inp[6]
    vka_dam = inp[7]/inp[6]
    prsity_dam = inp[8]
    up_depth = inp[9]
    down_depth = inp[10]
    num_el = int(inp[11])
    
    mf,dis,bas,chd,upw,nwt,oc,pond_cells = get_mf(modelname=modelname,
                                                  ztop=ztop,
                                                  hk_all=hk_all,
                                                  vka_all=vka_all,
                                                  prsity_all=prsity_all,
                                                  sy_all=prsity_all,
                                                  hk_ash=hk_ash,
                                                  vka_ash=vka_ash,
                                                  prsity_ash=prsity_ash,
                                                  sy_ash=prsity_ash,
                                                  hk_dam=hk_dam,
                                                  vka_dam=vka_dam,
                                                  sy_dam=prsity_dam,
                                                  prsity_dam=prsity_dam,
                                                  up_depth=up_depth,
                                                  down_depth=down_depth,
                                                  num_el=num_el,
                                                  b_d = 1000.0)
    
    t,num_particles = ccr.model_output(modelname,mf,dis,bas,chd,upw,nwt,oc,pond_cells)

    return t,num_particles
