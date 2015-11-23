#!/usr/bin/env python
__version__="0.1"
import sys
import os
help_desc="""
This program performs a QA procedure to identify voxels in which the different echo time series have different
slow drifts. In order to do this, the program first performs a fit to the multi-echo data to obtain two time series
(one is the T2* time series and the other one is the S0 time series). It then looks at the residual of this fit to
come with a single QA value per voxel. More specifically, the QA number is the standard deviation of a fit of 10th 
order legendre polynomials to the residual. We compute the standard deviation of this fit, instead of the stardard
deviation of the residual itself to be less senstitive to voxels with large residuals that do not necessarily signal
a slow drift in the residuals (e.g., a slow difference in drift in the echoes).

Static S0 and T2* Maps
----------------------
In addition, this program can also be used to generate static T2* and S0 maps using a non-linear optimizer. To request
the computation of such maps you must use the --get_static_maps option. Other parameters associated with this particular
procdure are:
    
   --So_init:   Initial Guess for the S0 value.  Default=2500
   --T2s_init:  Initial Guess for the T2s value. Default=40ms
   --So_min:    Lowest admissible S0 value.      Default=100
   --So_max:    Highest admissible S0 value.     Default=10000
   --T2s_min:   Lowest admissible T2s value.     Default=10ms
   --T2s_max:   Highest admissible T2s value.    Default=300ms
    
If none of these are provided, the software will use the default values. The addtional outputs, when static maps are requested, 
include:
    
    * <prefix>.sTE.t2s.nii: Voxel-wise static T2* map.
    * <prefix>.sTE.S0.nii:  Voxel-wise static So map.
    * <prefix>.sTE.SSE.nii: Voxel-wise Squared Standard Errors associated with the non-linear static fit.
    * <prefix>.sTE.mask.nii:Mask with voxels for which the non-linear optimizer was able to find values for t2s and So within
      the provided ranges. I usually find a few bad voxels in the ventricles and CSF around the brain. Other than that, the 
      mask should contain the whole intracranial volume. (1=good voxels | 0=bad voxels).
    
Dependences:
------------
This program was coded and tested with Python 2.7.10.
   
This program requires:

   * Numpy:           Python package for manipulation of numeric arrays. To install, please run: pip install numpy.
   * Scipy:           Scientific Library for Python. To install, please run: pip install scipy.
   * Nibabel:         Library to load NeuroImaging datasets. To install, please run: pip install nibabel.
   * SKLearn:         Python modules for machine learning and data mining. To install, please run: pip install sklearn
   * multiprocessing: Backport for multiprocesing package. To install, please run: pip install multiprocessing.
   * argparse:        Python argument parser. To install, please run: pip install argparse
   
How to run:
-----------
This program will take require to a minimum the following inputs:
    
   * ME Dataset: This will be a 4D dataset in which the different echoes have been concatenated in the Z-direction. This
                 means that for a dataset with Ne echoes and dimensions (Nx,Ny,Nz,Nt) the input to this QA program should
                 be a dataset with dimensions (Nx,Ny,Nz*Ne,Nt). Such datasets can be easily created using AFNI program
                 3dZcat (e.g., 3dZcat -prefix MEdataset.nii E01.nii E02.nii E03.nii).
                 To pass the ME dataset please use -d or --orig_data

   * Echo times: The program needs to know the echo times (in milisenconds) used during data acquisition. There are two
                 ways to provide the echo times. You can provide them on the command line using -e and a list of echo times
                 separated by commas. You can also provide them via a text file with --tes_file. The file should contain
                 a single line of text with the echo times separated by commas (e.g., 12,24,35). 

   * Debug: By default the program will write a single output (<prefix>.QAC.ResidualFit_stdv.nii) that contains a map of the standard deviation of the legengre
            fit to the residuals of the fit to the TE-dependence model. In this map, voxels with a very high value should be
            explored as most likely they corresponds to voxels where the echoes drifted differently across time. Values in the
            thousands usually mean problematic voxels.
            If you pass the --debug option to the program, it will also write the following additional outputs:
                
                * <prefix>.QAC.Residual.nii:         timeseries of the residuals from the TE-dependence fit.
                * <prefix>.QAC.ResidualFit.nii:      legendre polynomial fit to <prefix>.QAC.Residual.nii
                * <prefix>.QAC.ResidualFit_stdv.nii: voxel-wise standard deviation across time of <prefix>.QAC.ResidualFit.nii
                * <prefix>.QAC.SSE.nii: voxel-wise Sum of Squared Errors for the TE-fit. This is a less robust QA metric (I believe).
                * <prefix>.QAC.SSE.rank.nii: a ranked version of <prefix>.QAC.SSE.nii ranging from 0 to 100.
                             
A sample command line would be:
-------------------------------

    CASE (1) Do QA Only, asking for all intermediate files
    
    $ python meqa.py -d MEdataset.nii -e 12,24,35 --prefix MEdatataset --debug
    
    CASE (2) Do QA and Static Maps using default settings.
    
    $ python meqa.py -d MEdataset.nii -e 12,24,35 --prefix MEdatataset --get_static_maps
    
    Case (3) Do QA and Static Maps, but using a different initial guess for voxel-wise T2*.
    
    $ python meqa.py -d MEdataset.nii -e 12,24,35 --prefix MEdatataset --get_static_maps --T2s_init 30
    
One way to explore the output of this program is:
-------------------------------------------------    
    1.  Open AFNI
    2.  Select as underlay the first echo E01.nii
    3.  Open a graph window (press Graph in AFNI).
    4.  On the graph window, select Detrend=0 (Opt --> Detrend = 0)
    5.  On the graph window, overlay the time series for the additional echoes (Opt --> Tran1D = Dataset #N --> Select all necessary datasets --> Press "Save + Close")
    6.  You may want to have thinner lines in the Graph window for better visualization. For this (Opt --> Colors, Etc. --> Uncheck all the "Use Thick Lines" entries)
    7.  Open a second AFNI viewer (Press New in AFNI)
    8.  Select <prefix>.QAC.Residual.nii as the underlay
    9.  Select <prefix>.QAC.ResidualFit_stdv.nii as overlay
    10. Press "Define Overlay -->" to open additional overlay controls in AFNI.
    11. Show only positives (Check Pos? On the lower bottom of the new panel). This will give more color range for visualization purposes.
    12. Deselect "autoRange:"
    13. Set the maximum value to 1000 (box right below autoRange:"
    14. Go to the voxels with the highest values, and the explore their timeseries via the graph window associated with AFNI viewer [A].
    
"""
# =================================================================================================================
# =================================         FUNCTION DEFINITIONS            =======================================
# =================================================================================================================   

# === FUNCTION: dep_check
def dep_check():
    print("++ INFO [Main]: Checking for dependencies....")
    fails                = 0
    modules = set(["numpy","argparse","scipy","sklearn","multiprocessing","nibabel"])
    
    for m in modules:
        try:
            __import__(m)
        except ImportError:
            fails += 1
            print("++ ERROR [Main]: Can't import Module %s. Please install." % m)
                
    if fails == 0:
        print(" +              All Dependencies are OK.")
    else:
        print(" +               All dependencies not available. Please install according to above error messages.")
        print(" +               Program exited.")
        sys.exit()
        
# === FUNCTION: niiLoad
def niiLoad(path):
    """
    This function reads nifti datasets
    
    Parameters:
    -----------
    path: string containing the path to the NIFTI file you want to read.
    
    Returns:
    --------
    data: a numpy array with the data
    aff:  the affine transformation associated with the dataset. Needed in order to write to disk again.
    head: the header of the nifti dataset.
    """
    mepi_dset       = nib.load(path)
    data            = mepi_dset.get_data()
    aff             = mepi_dset.get_affine()
    head            = mepi_dset.get_header()
    head.extensions = []
    head.set_sform(head.get_sform(),code=1)
    return data,aff,head

# === FUNCTION: mask4MEdata
def mask4MEdata(data):
    """
    This function will create a mask for ME data taking into
    account the time series of all echoes. It expects datasets already masked in AFNI. What this function does
    is to compose a new mask that only includes voxels that were included by 3dAutomask for all echoes.
	
    Paramters:
    ----------
    data: this is a (Nx,Ny,Nz,Ne,Nt) array with ME data.

    Returns:
    --------
    mask: this is a (Nx,Ny,Nz) marking all voxels that have valid time-series for all echo times.
    """
    # Create mask taking into account all echoes
    Nx,Ny,Nz,Ne,Nt = data.shape
    mask           = np.ones((Nx,Ny,Nz),dtype=np.bool)
    for i in range(Ne):
        tmpmask = (data[:,:,:,i,:] != 0).prod(axis=-1,dtype=np.bool)
	mask    = mask & tmpmask
    return mask

# === FUNCTION: niiwrite_nv
def niiwrite_nv(data,mask,temp_path,aff,temp_header):
        """
        This function will write NIFTI datasets

        Parameters:
        ----------
        data: this is (Nv, Nt) or (Nv,) array. No z-cat ME datasets allowed.
        mask: this is (Nx,Ny,Nz) array with Nv entries equal to True. This is an intracranial voxel mask.
        temp_path: this is the output directory.
        aff: affine transformation associated with this dataset.
        temp_header: header for the dataset.  
        
        Returns:
        --------
        None.
        """
        Nx,Ny,Nz   = mask.shape
        if (data.ndim ==1):
                temp       = np.zeros((Nx,Ny,Nz),order='F')
                temp[mask] = data
        if (data.ndim ==2):
                        _,Nt = data.shape
                        temp         = np.zeros((Nx,Ny,Nz,Nt),order='F')
                        temp[mask,:] = data
        if (data.ndim ==3):
                        Nv, Ne, Nt   = data.shape
                        temp         = np.zeros((Nx,Ny,Nz,Nt),order='F')
                        temp[mask,:] = np.squeeze(data[:,0,:])
                        for e in range(1,Ne):
                                aux       = np.zeros((Nx,Ny,Nz,Nt),order='F')
                                aux[mask,:] = np.squeeze(data[:,e,:])
                                temp = np.concatenate((temp,aux),axis=2)

        outni      = nib.Nifti1Image(temp,aff,header=temp_header)
        outni.to_filename(temp_path)
        print(" +              Dataset %s written to disk" % (temp_path))

# === FUNCTION: linearFit
def linearFit_perVoxel(item):
  Ne     = item['Ne']
  Nt     = item['Nt']
  DeltaS = item['DeltaS']
  Smean  = item['Smean']
  tes    = item['tes']
  
  datahat  = np.zeros((Ne,Nt)) 
  B        = np.zeros((Ne,2))
  B[:,0]   = Smean
  B[:,1]   = -(1.0/tes.mean())*Smean*tes
  rcond    = np.linalg.cond(B)
  res      = np.linalg.lstsq(B,DeltaS)

  drho         = res[0][0,:]
  dkappa       = res[0][1,:]
  fit_residual = res[1]
  datahat      = np.dot(B,res[0]) + np.tile(Smean.reshape((Ne,1)),Nt)
  db           = np.zeros(res[0].shape)
  db[1,:]      = dkappa
  return {'dkappa':dkappa,'drho':drho, 'fit_residual':fit_residual,'rcond':rcond,'datahat':datahat}

def linearFit(data, tes,Ncpus, dataMean=None):
    """
    This function will compute the fit of ME data using a least square approach
    
    Parameters:
    -----------
    data: ME dataset masked. This is a (Nv,Ne,Nt) numpy array
    tes: numpy array with the echo times used during acquisition.
    Ncpus: number of CPUs.
    dataMean: mean across time of the time series
	
    Returns:
    -------
    dkappa: Time series of TE dependent fluctuations.
    drho:   Time series of non-TE dependent fluctuations.
    fit_residual: voxel-wise and time point wise residuals
    rcond:         condition number for each voxel.
    datahat: Estimate of the data after fitting
    """
    Nv,Ne,Nt     = data.shape
    rcond        = np.zeros((Nv,))
    drho         = np.zeros((Nv,Nt))
    dkappa       = np.zeros((Nv,Nt))
    fit_residual = np.zeros((Nv,Nt))
    datahat      = np.zeros((Nv,Ne,Nt))
    if dataMean is  None:
	dataMean     = data.mean(axis=2)
    pool   = Pool(processes=Ncpus)
    result = pool.map(linearFit_perVoxel, [{'DeltaS':data[v,:,:] - np.tile(dataMean[v,:].reshape((Ne,1)),Nt),'Smean':dataMean[v,:],'tes':tes,'Ne':int(Ne),'Nt':int(Nt)} for v in np.arange(Nv)]) 
    
    for v in range(Nv):
	rcond[v]          = result[v]['rcond'] 
	drho[v,:]         = result[v]['drho'] 
	dkappa[v,:]       = result[v]['dkappa'] 
	fit_residual[v,:] = result[v]['fit_residual'] 
	datahat[v,:,:]    = result[v]['datahat'] 
    return dkappa,drho,fit_residual,rcond,datahat

# === FUNCTION: getMeanByPolyFit	
def getMeanByPolyFit(data,Ncpus,polort=4):
   """ 
   This function computes the mean across time for all voxels and echoes using
   legendre polynomial fitting. More robust against slow dirfts

   Parameters
   ----------
   data:   ME dataset (Nv,Ne,Nt)
   polort: order for the legendre polynomials to fit.
   Returns
   -------
   mean:   (Nv,Ne)
   """
   Nv,Ne,Nt = data.shape
   aux = np.reshape(data.copy(),(Nv*Ne,Nt))
   drift = np.zeros((Nt,polort))
   x     = np.linspace(-1,1-2/Nt,Nt)
   for n in range(polort):
      drift[:,n]=np.polyval(legendre(n),x)
   # 2. Fit polynomials to residuals
   linRegObj= LinearRegression(normalize=False,fit_intercept=False, n_jobs=Ncpus)
   linRegObj.fit(drift,aux.T)
   mean    = np.reshape(linRegObj.coef_[:,0],(Nv,Ne))
   return mean

# === FUNCTION: computeQA
def computeQA(data,tes,Ncpus,data_mean=None):
    """
    Simple function to compute the amount of variance in the data that is explained by
    the ME fit
    
    Parameters:
    -----------
    data: ME dataset (Nv,Ne,Nt)
    tes:  Echo times used during acquisition
    Ncpus: number of CPUs for parallelization across voxels
    data_mean: mean across time of the input data.
    Returns:
    --------
    SSE:  Voxelwise mean across time of the Sum of Squared Errors for the TE fit. (Nv,)  
    rankSSE: ranked version from 0 to 100 of SSE (good for kappa computation)     (Nv,)
    """
    Nv,Ne,Nt        = data.shape
    if data_mean is None:
       data_mean = getMeanByPolyFit(data,Ncpus,polort=4)
    data_demean            = data - data_mean[:,:,np.newaxis]
    dkappa, drho,residual,rcond,data_hat = linearFit(data,tes,Ncpus,data_mean)
    data_hat_mean                        = getMeanByPolyFit(data_hat,Ncpus,polort=4)
    data_demean_hat                      = data_hat - data_hat_mean[:,:,np.newaxis]

    SSE             = ((data_demean - data_demean_hat)**2).sum(axis=-1).max(axis=-1)
    rankSSE         = 100.*rankdata(1./SSE)/Nv
    
    Npoly = np.int(10)
    print("++ INFO: Generating legendre polynomials of order [%i]." % Npoly)
    x = np.linspace(-1,1-2/Nt,Nt)
    drift = np.zeros((Nt, Npoly))
    for n in range(Npoly):
        drift[:,n] = np.polyval(legendre(n),x)
    # 2. Fit polynomials to residuals
    linRegObj        = LinearRegression(normalize=False,fit_intercept=True, n_jobs=Ncpus)
    linRegObj.fit(drift,residual.T)
    fit2residual      = np.transpose(linRegObj.predict(drift))
    fit2residual_std  = fit2residual.std(axis=1)
    fit2residual_norm = ( fit2residual - fit2residual.mean(axis=1)[:,np.newaxis] ) / fit2residual.mean(axis=1)[:,np.newaxis]
    meqa_norm         = fit2residual_norm.std(axis=1)
    return SSE,rankSSE,residual,fit2residual,fit2residual_std,meqa_norm

# === FUNCTION: make_static_maps_opt
def objective(x,Sv,aux_tes):
    return np.sqrt(sum((Sv-(x[0]*np.exp(-aux_tes/x[1])))*(Sv-(x[0]*np.exp(-aux_tes/x[1])))))

def make_static_opt_perVoxel(item):
   data_mean = item['data_mean']
   tes       = item['tes']
   So_init   = item['So_init']
   T2s_init  = item['T2s_init']
   So_min    = item['So_min']
   So_max    = item['So_max']
   T2s_min   = item['T2s_min']
   T2s_max   = item['T2s_max']
   Optimizer = item['Optimizer']
   result      = opt.minimize(objective, [So_init,T2s_init],args=(data_mean,tes), method=Optimizer,bounds=((So_min, So_max), (T2s_min, T2s_max)))
   v_fiterror  = result.fun
   v_S0        = result.x[0]
   v_t2s       = result.x[1]
   if (~result.success) or (v_t2s >= T2s_max-.05) or (v_t2s <= T2s_min+.05) or (v_S0 >= So_max-.05) or (v_S0 <= So_min+.05): 
     v_badFit  = 1
   else:
     v_badFit  = 0
   return {'v_S0':v_S0, 'v_t2s':v_t2s, 'v_fiterror':v_fiterror, 'v_badFit':v_badFit}

def make_static_maps_opt(data_mean,tes,Ncpus,So_init=2500,T2s_init=40,So_min=100,So_max=10000,T2s_min=10, T2s_max=300,Optimizer='SLSQP'):
   """
   This function computes static maps of S0 and T2s using scipy optimization

   Parameters:
   -----------
   data_mean: Mean across time of the ME dataset (Nv,Ne)
   tes:       Echo times used to acquire the data
   So_init:   Initial Guess for the S0 value.  Default=2500
   T2s_init:  Initial Guess for the T2s value. Default=40ms
   So_min:    Lowest admissible S0 value.      Default=100
   So_max:    Highest admissible S0 value.     Default=10000
   T2s_min:   Lowest admissible T2s value.     Default=10ms
   T2s_max:   Highest admissible T2s value.    Default=300ms
   Optimizer: Optimization Algorithm.          Default=SLSQP
   Returns:
   --------
   S0:        Static S0 map  (Nv,)
   t2s:       Static T2s map (Nv,)
   SSE:       Sum of Squared Errors (Nv,)
   BadFits:   Voxels marked as bad fits by the optimizer (Nv,)
   """
   print(" + INFO [make_static_maps_opt]: Initial conditions [So=%i, T2s=%i]" % (So_init, T2s_init))
   print(" + INFO [make_static_maps_opt]: Bounds So=[%i,%i] & T2s=[%i,%i]" % (So_min, So_max, T2s_min, T2s_max))
   print(" + INFO [make_static_maps_opt]: Optimizer = %s" % Optimizer)
   
   Nv,Ne    = data_mean.shape
   S0       = np.zeros(Nv,)
   t2s      = np.zeros(Nv,)
   badFits  = np.zeros(Nv,)
   fiterror = np.zeros(Nv,)

   print(" +              Multi-process Static Map Fit -> Ncpus = %d" % Ncpus)
   pool   = Pool(processes=Ncpus)
   result = pool.map(make_static_opt_perVoxel, [{'data_mean':data_mean[v,:],'tes':tes,'So_init':So_init,'So_max':So_max,'So_min':So_min,'T2s_init':T2s_init,'T2s_max':T2s_max,'T2s_min':T2s_min,'Optimizer':Optimizer} for v in np.arange(Nv)]) 
   for v in np.arange(Nv):
     S0[v]  = result[v]['v_S0']
     t2s[v] = result[v]['v_t2s']
     fiterror[v] = result[v]['v_fiterror']
     badFits[v]  = result[v]['v_badFit']
   print(" + INFO [make_static_maps_opt]: Number of Voxels with errors: %i" % badFits.sum())
   return S0, t2s, fiterror, badFits

# =================================================================================================================
# =================================             MAIN PROGRAM                =======================================
# =================================================================================================================

if __name__=='__main__':
    print("------------------------------------")
    print("-- SFIM ME-QA version %s         --" % __version__)
    print("------------------------------------")
    dep_check()
    import numpy              as np
    import nibabel            as nib
    import scipy.optimize     as opt
    from argparse             import ArgumentParser,RawTextHelpFormatter
    from multiprocessing      import cpu_count,Pool
    from scipy.special        import legendre
    from sklearn.linear_model import LinearRegression
    from scipy.stats          import rankdata
    
    # =================================================================================================================
    # =================================         PARSING OF INPUT PARAMETERS     =======================================
    # =================================================================================================================
    parser = ArgumentParser(description=help_desc,formatter_class=RawTextHelpFormatter)
    statFitGrp = parser.add_argument_group('Arguments for static T2* and S0 fits')
    parser.add_argument("-d","--orig_data",        dest='data_file',      help="Spatially concatenated Multi-Echo Dataset",type=str,   default=None)
    parser.add_argument("-e","--TEs",              dest='tes',            help="Echo times (in ms) ex: 15,39,63",          type=str,   default=None)
    parser.add_argument(     "--tes_file",         dest='tes_file',       help="Path to file with Echo time information",  type=str,   default=None)
    parser.add_argument(     "--out_dir",          dest='out_dir',        help="Output Directory. Default: current directory",type=str,default='./')
    parser.add_argument(     "--prefix",           dest='prefix',         help="Output File Prefix. Default = sbj",type=str,default='sbj')
    parser.add_argument(     "--ncpus",            dest='Ncpus',          help='Number of cpus available. Default will be #available/2', default=None, type=int)
    parser.add_argument(     "--mask",             dest='mask_file',      help='Path to the mask to use during the analysis. If not provided, one will be computed automatically',type=str, default=None)
    parser.add_argument(     "--debug",            dest='debug',          help='Flag to write out additional files',action='store_true')
    
    statFitGrp.add_argument("--get_static_maps", dest='do_static_fit', help='Flag to ask for static T2* and S0 maps', action='store_true')
    statFitGrp.add_argument("--So_init", dest='So_init',  help='Initial Guess for the S0 value.  Default=2500', type=float, default=2500)
    statFitGrp.add_argument("--T2s_init",dest='T2s_init', help='Initial Guess for the T2s value. Default=40ms', type=float, default=40)
    statFitGrp.add_argument("--So_min",  dest='So_min',   help='Lowest admissible S0 value.      Default=100', type=float,  default=100)
    statFitGrp.add_argument("--So_max",  dest='So_max',   help='Highest admissible S0 value.     Default=10000', type=float,default=10000)
    statFitGrp.add_argument("--T2s_min", dest='T2s_min',  help='Lowest admissible T2s value.     Default=10ms', type=float, default=10)
    statFitGrp.add_argument("--T2s_max", dest='T2s_max',  help='Highest admissible T2s value.    Default=300ms', type=float,default=300)
    options  = parser.parse_args()
    So_init  = float(options.So_init)
    So_min   = float(options.So_min)
    So_max   = float(options.So_max)
    T2s_init = float(options.T2s_init)
    T2s_min  = float(options.T2s_min)
    T2s_max  = float(options.T2s_max)
    
    # If no number of CPUs is provided, we will use half the available number of CPUs or 1 if only 1 is available
    if cpu_count()==1:
        Ncpus = 1;
    else:
        if (options.Ncpus is None) or (options.Ncpus > cpu_count()):
            Ncpus = int(cpu_count()/2)
        else:
            Ncpus = int(options.Ncpus) 
    print("++ INFO [Main]: Number of CPUs to use: %d" % (Ncpus))
    
    # Control all necessary inputs are available
    # ------------------------------------------
    if options.tes is None and options.tes_file is None:
        print("++ Error: No information about echo times provided. Please do so via --TEs or --tes_file")
        sys.exit()
    if (not (options.tes is None)) and (not (options.tes_file is None)):
        print("++ Error: Echo times provided in two different ways (--TEs and --tes_file). Please select only one input mode.")
        sys.exit()
    if options.data_file is None:
        print("++ Error: No ME input dataset provided. Please provide one using the -d or --orig_data parameter.")
        sys.exit()
    
    # Control for existence of files and directories
    # ----------------------------------------------
    if not os.path.exists(options.data_file):
        print("++ Error: Datafile [%s] does not exists." % options.data_file)
        sys.exit()
    if options.tes_file!=None and (not os.path.exists(options.tes_file)):
        print("++ Error: Echo Times file [%s] does not exists." % options.tes_file)
        sys.exit()
    if (not os.path.exists(options.out_dir)) or (not os.path.isdir(options.out_dir)):
        print("++ Error: Output directory [%s] does not exists." % options.out_dir)
    
    # Set all output paths
    # --------------------
    outputDir = os.path.abspath(options.out_dir)
    
    # =================================================================================================================
    # =================================         LOADING DATA INTO MEMORY        =======================================
    # =================================================================================================================
    
    # Load echo times information
    # ---------------------------
    if options.tes!=None:
        print("++ INFO [Main]: Reading echo times from input parameters.")
        tes      = np.fromstring(options.tes,sep=',',dtype=np.float32)
    if options.tes_file!=None:
        print("++ INFO [Main]: Reading echo times from input echo time file.")
        tes         = np.loadtxt(options.tes_file,delimiter=',')
    Ne = tes.shape[0]
    print(" +              Echo times: %s" % (str(tes)))
    
    # Load ME-EPI data
    # ----------------
    print("++ INFO [Main]: Loading ME dataset....")
    mepi_data,mepi_aff,mepi_head  = niiLoad(options.data_file)
    Nx,Ny,Nz,Nt                   = mepi_data.shape
    Nz                            = Nz/Ne # Because the input was the Z-concatenated dataset
    mepi_data                     = mepi_data.reshape((Nx,Ny,Nz,Ne,Nt),order='F')
    print(" +              Dataset dimensions: [Nx=%i,Ny=%i,Nz=%i,Ne=%i,Nt=%i]" % (Nx,Ny,Nz,Ne,Nt))
    
    # =================================================================================================================
    # =================================        ORIGINAL MASK / RESHAPE          =======================================
    # =================================================================================================================
    origMask_path   = os.path.join(outputDir,options.prefix+'.mask.orig.nii')     
    if options.mask_file==None:
        print("++ INFO [Main]: Generating initial mask from data.")
        mask           = mask4MEdata(mepi_data)
        if options.debug:
            niiwrite_nv(mask[mask],mask,options.out_dir+options.prefix+'.mask.orig.nii',mepi_aff ,mepi_head)
    else:
        print("++ INFO [Main]: Using user-provided mask.")
        if not os.path.exists(options.data_file):
            print("++ Error: Provided mask [%s] does not exists." % options.mask_file)
            sys.exit()
        mask,_,_       = niiLoad(options.mask_file)
        mask = (mask>0)
        
    Nv             = np.sum(mask)
    print(" +              Number of Voxels in mask [Nv=%i]" % Nv)
    # Put the data into a 3Dx format (Nvoxels, Nechoes, Ndatapoints)
    SME      = mepi_data[mask,:,:].astype(float) #(Nv,Ne,Nt)
    
    # =================================================================================================================
    # =================================        COMPUTE MEAN ACROSS TIME         =======================================
    # =================================================================================================================
    print("++ INFO [Main]: Computing Mean across time for all echoes....")
    # There are two options here:
    #   (1) Simply use the mean command.
    #       Smean_case01 = SME.mean(axis=-1)
    #   (2) Compute the mean at the same time you fit some legendre polynomials, to remove the influence
    #   of small drits in the computation. It should make almost no difference.
    Smean_case02 = getMeanByPolyFit(SME,Ncpus,polort=4)
    SME_mean     = Smean_case02
    
    # =================================================================================================================
    # =================================                PERFORM QA               =======================================
    # =================================================================================================================
    
    print("++ INFO [Main]: Quality Assurance....")
    QA_SSE_path                   = os.path.join(outputDir,options.prefix+'.QAC.SSE.nii')
    QA_SSE_rank_path              = os.path.join(outputDir,options.prefix+'.QAC.SSE_rank.nii')
    QA_Residual_path              = os.path.join(outputDir,options.prefix+'.QAC.Residual.nii')
    QA_ResidualFit_path           = os.path.join(outputDir,options.prefix+'.QAC.ResidualFit.nii')
    QA_ResidualFit_stdv_path      = os.path.join(outputDir,options.prefix+'.QAC.ResidualFit_stdv.nii')
    QA_MEQAnorm_path              = os.path.join(outputDir,options.prefix+'.QAC.MEQAnorm.nii')

    print(" +              Computing QA metrics from data.")
    QA_SSE,QA_SSE_Rank,QA_Residual, QA_ResidualFit, QA_ResidualFit_stdv, QA_MEQAnorm = computeQA(SME,tes,Ncpus,data_mean=SME_mean)
    if options.debug:
        niiwrite_nv(QA_SSE_Rank,         mask, QA_SSE_rank_path, mepi_aff ,mepi_head)
        niiwrite_nv(QA_SSE,              mask, QA_SSE_path     , mepi_aff ,mepi_head)
        niiwrite_nv(QA_Residual,         mask, QA_Residual_path, mepi_aff ,mepi_head)
        niiwrite_nv(QA_ResidualFit,      mask, QA_ResidualFit_path, mepi_aff ,mepi_head)
    niiwrite_nv(QA_ResidualFit_stdv, mask, QA_ResidualFit_stdv_path, mepi_aff ,mepi_head)
    niiwrite_nv(QA_MEQAnorm, mask, QA_MEQAnorm_path, mepi_aff ,mepi_head)
 
    # =================================================================================================================
    # =================================                STATIC FIT               =======================================
    # =================================================================================================================
    if options.do_static_fit:
        print("++ INFO [Main]: Static T2* and S0 maps requested...")
        stFit_S0_path  = os.path.join(outputDir,options.prefix+'.sTE.S0.nii')
        stFit_t2s_path = os.path.join(outputDir,options.prefix+'.sTE.t2s.nii')
        stFit_SSE_path = os.path.join(outputDir,options.prefix+'.sTE.SSE.nii')
        stFit_bVx_path = os.path.join(outputDir,options.prefix+'.sTE.mask.nii')
        mask_bad_staticFit = np.zeros((Nv,), dtype=bool)
        
        S0, t2s, SSE, mask_bad_staticFit = make_static_maps_opt(SME_mean,tes,Ncpus,So_init=So_init,T2s_init=T2s_init,So_min=So_min,So_max=So_max,T2s_min=T2s_min, T2s_max=T2s_max)
       
        mask_bad_staticFit = np.logical_not(mask_bad_staticFit)
        niiwrite_nv(S0                ,mask,stFit_S0_path, mepi_aff ,mepi_head)
        niiwrite_nv(t2s               ,mask,stFit_t2s_path,mepi_aff ,mepi_head)
        niiwrite_nv(SSE               ,mask,stFit_SSE_path,mepi_aff ,mepi_head)
        niiwrite_nv(mask_bad_staticFit,mask,stFit_bVx_path,mepi_aff ,mepi_head) 
