import sys
import numpy.ma as ma
import numpy   as np
import nibabel as nib
import scipy.optimize as opt
from scipy.stats import scoreatpercentile,zscore,skew,rankdata
import matplotlib.pyplot as plt
import time
import datetime
from scipy.special        import legendre
from sklearn.linear_model import LinearRegression
from multiprocessing      import cpu_count,Pool

def computeFFT(ts,TR):
    """
    This function computes the fft of the ICA representative timeseries and
    writes the output to disk
    
    Parameters
    ----------
    ts: input timseries (Nc,Nt)
    TR: repetition time in seconds (float)
 
    Returns
    -------
    FFT:  Power spectrum for all the signals (Nc,Nt)
    freq: Frequencies for the x-axis (Nt)
    """
    Nc,Nt = ts.shape
    FFT   = np.zeros(ts[:,Nt/2:].shape)
    for c in range(Nc):
        its      = ts[c,:]
        iFFT     = np.fft.fftshift(abs(np.fft.fft(its,axis=0)))
        FFT[c,:] = iFFT[Nt/2:]
        freq     = np.fft.fftshift(np.fft.fftfreq(Nt,float(TR)))
        freq     = freq[Nt/2:]
    return FFT, freq

def getelbow(ks):
   nc = ks.shape[0]
   coords = np.array([np.arange(nc),ks])
   p  = coords - np.tile(np.reshape(coords[:,0],(2,1)),(1,nc))
   b  = p[:,-1]
   b_hat = np.reshape(b/np.sqrt((b**2).sum()),(2,1))
   proj_p_b = p - np.dot(b_hat.T,p)*np.tile(b_hat,(1,nc))
   d = np.sqrt((proj_p_b**2).sum(axis=0))
   k_min_ind = d.argmax()
   k_min  = ks[k_min_ind]
   return k_min_ind

def getMeanByPolyFit(data,polort=4):
   """ 
   This function computes the mean across time for all voxels and echoes using
   legendre polynomial fitting. More robust against slow dirfts

   Parameters
   ----------
   data:   ME dataset (Nv,Ne,Nt)
   polort: order for the legendre polynomials to fit. Default=4
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
   linRegObj= LinearRegression(normalize=False,fit_intercept=False)
   linRegObj.fit(drift,aux.T)
   mean    = np.reshape(linRegObj.coef_[:,0],(Nv,Ne))
   return mean

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
			temp[mask,:] = data[:,0,:]
			for e in range(1,Ne):
				aux       = np.zeros((Nx,Ny,Nz,Nt),order='F')
				aux[mask,:] = data[:,e,:]
				temp = np.concatenate((temp,aux),axis=2)

	outni      = nib.Nifti1Image(temp,aff,header=temp_header)
	outni.to_filename(temp_path)
	print " +              Dataset %s written to disk" % (temp_path)

def make_optcom(data,t2s,tes):
	"""
	Generates the optimally combined time series. 

	Parameters:
	-----------
	data: this is the original ME dataset with the mean in a (Nv,Ne,Nt) array.
	t2s:  this is the static T2s map in a (Nv,) array.
	tes:  echo times in a (Ne,) array.

	Returns:
	--------
	octs: optimally combined time series in a (Nv,Nt) array.	
	"""
	Nv,Ne,Nt = data.shape
	ft2s  = t2s[:,np.newaxis]
	alpha = tes * np.exp(-tes /ft2s)
	alpha = np.tile(alpha[:,:,np.newaxis],(1,1,Nt))
	octs  = np.average(data,axis = 1,weights=alpha)
	return octs

def mask4MEdata(data):
	"""
	This function will create a mask for ME data taking into
  	account the time series of all echoes
	
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

def linearFit_perVoxel(item):
  Ne     = item['Ne']
  Nt     = item['Nt']
  DeltaS = item['DeltaS']
  Smean  = item['Smean']
  tes    = item['tes']
  
  datahat  = np.zeros((Ne,Nt))
  datahat2 = np.zeros((Ne,Nt)) 
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
  datahat2     = np.dot(B,db) + np.tile(Smean.reshape((Ne,1)),Nt) 
  return {'dkappa':dkappa,'drho':drho, 'fit_residual':fit_residual,'rcond':rcond,'datahat':datahat,'datahat2':datahat2}

def linearFit(data, tes,Ncpu, dataMean=None):
	"""
	This function will compute the fit of ME data using a least square approach
	(although not sure why this is not based on the log)

	Parameters:
	-----------
	data: ME dataset masked. This is a (Nv,Ne,Nt) numpy array
	tes: numpy array with the echo times used during acquisition.

	Returns:
	-------
	dkappa: Time series of TE dependent fluctuations.
	drho:   Time series of non-TE dependent fluctuations.
	fit_residual: voxel-wise and time point wise residuals
	rcond:         condition number for each voxel.
	datahat: Estimate of the data after fitting
   datahat2: Estimate of the data after fitting and removing S0 changes
	"""
	Nv,Ne,Nt     = data.shape
	rcond        = np.zeros((Nv,))
	drho         = np.zeros((Nv,Nt))
	dkappa       = np.zeros((Nv,Nt))
	fit_residual = np.zeros((Nv,Nt))
	datahat      = np.zeros((Nv,Ne,Nt))
	datahat2     = np.zeros((Nv,Ne,Nt))
	if dataMean is  None:
		dataMean     = data.mean(axis=2)
	pool   = Pool(processes=Ncpu)
	result = pool.map(linearFit_perVoxel, [{'DeltaS':data[v,:,:] - np.tile(dataMean[v,:].reshape((Ne,1)),Nt),'Smean':dataMean[v,:],'tes':tes,'Ne':int(Ne),'Nt':int(Nt)} for v in np.arange(Nv)]) 
    
	for v in range(Nv):
		rcond[v]          = result[v]['rcond'] 
		drho[v,:]         = result[v]['drho'] 
		dkappa[v,:]       = result[v]['dkappa'] 
		fit_residual[v,:] = result[v]['fit_residual'] 
		datahat[v,:,:]    = result[v]['datahat'] 
		datahat2[v,:,:]   = result[v]['datahat2'] 
	return dkappa,drho,fit_residual,rcond,datahat,datahat2

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

def make_static_maps_opt(data_mean,tes,Ncpu,So_init=2500,T2s_init=40,So_min=100,So_max=10000,T2s_min=10, T2s_max=300,Optimizer='SLSQP'):
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
   print " + INFO [make_static_maps_opt]: Initial conditions [So=%i, T2s=%i]" % (So_init, T2s_init)
   print " + INFO [make_static_maps_opt]: Bounds So=[%i,%i] & T2s=[%i,%i]" % (So_min, So_max, T2s_min, T2s_max)
   print " + INFO [make_static_maps_opt]: Optimizer = %s" % Optimizer
   
   Nv,Ne    = data_mean.shape
   S0       = np.zeros(Nv,)
   t2s      = np.zeros(Nv,)
   badFits  = np.zeros(Nv,)
   fiterror = np.zeros(Nv,)

   print " +              Multi-process Static Map Fit -> Ncpu = %d" % Ncpu
   pool   = Pool(processes=Ncpu)
   result = pool.map(make_static_opt_perVoxel, [{'data_mean':data_mean[v,:],'tes':tes,'So_init':So_init,'So_max':So_max,'So_min':So_min,'T2s_init':T2s_init,'T2s_max':T2s_max,'T2s_min':T2s_min,'Optimizer':Optimizer} for v in np.arange(Nv)]) 
   for v in np.arange(Nv):
     S0[v]  = result[v]['v_S0']
     t2s[v] = result[v]['v_t2s']
     fiterror[v] = result[v]['v_fiterror']
     badFits[v]  = result[v]['v_badFit']
   print " + INFO [make_static_maps_opt]: Number of Voxels with errors: %i" % badFits.sum()
   return S0, t2s, fiterror, badFits

def _characterize_this_component(item):
    B          = item['B']              #(Ne,Nv)
    X1         = item['X1']             #(Ne,Nv)
    X2         = item['X2']             #(Ne,Nv)
    aff        = item['aff']
    head       = item['head']
    mask       = item['mask']
    c          = item['c']
    outDir     = item['outDir']
    outPrefix  = item['outPrefix']
    F_MAX      = item['F_MAX']
    weight_map = item['weight_map']
    c_mask     = item['c_mask']
    writeOuts  = item['writeOuts']
    Ne, Nv     = B.shape
    # Write the fits for the differente echoes into file (very useful for debugging purposes)
    if writeOuts:
           niiwrite_nv(B.T,mask,outDir+outPrefix+'.chComp.EXTRA.Beta'+str(c).zfill(3)+'.nii',aff,head)
    alpha = (B**2).sum(axis=0)                     #(Nv,)

    # S0 Model
    coeffs_S0        = (B*X1).sum(axis=0)/(X1**2).sum(axis=0)        #(Nv,)
    
    ##### DAN VERSION
    #print "NEW THING SHAPE %s" % str((np.atleast_2d(X1[:,1])).shape)
    #coeffs_S0 = np.linalg.lstsq(np.atleast_2d(X1[:,1]).T,B)
    ### END OF DAN VERSION
    if writeOuts:
        niiwrite_nv((X1*np.tile(coeffs_S0,(Ne,1))).T,mask,outDir+outPrefix+'.chComp.EXTRA.S0Fit'+str(c).zfill(3)+'.nii',aff,head)
    SSE_S0           = (B - X1*np.tile(coeffs_S0,(Ne,1)))**2         #(Ne,Nv)
    SSE_S0           = SSE_S0.sum(axis=0)                            #(Nv,)
    F_S0             = ((alpha - SSE_S0)/(SSE_S0))/((Ne-1)/(1)) #(Nv,)
    F_S0_AllValues   = F_S0.copy()
    #F_S0[F_S0>F_MAX] = F_MAX
    
    # R2 Model
    #### DAN VERSION
    #coeffs_R2 = np.linalg.lstsq(np.atleast_2d(X2[:,1]).T,B)
    ### END OF DAN VERSION
    coeffs_R2        = (B*X2).sum(axis=0)/(X2**2).sum(axis=0)        #(Nv,)
    if writeOuts:
        niiwrite_nv((X2*np.tile(coeffs_R2,(Ne,1))).T,mask,outDir+outPrefix+'.chComp.EXTRA.R2Fit'+str(c).zfill(3)+'.nii',aff,head)
    SSE_R2           = (B - X2*np.tile(coeffs_R2,(Ne,1)))**2         #(Ne,Nv)
    SSE_R2           = SSE_R2.sum(axis=0)                            #(Nv,)
    F_R2             = ((alpha - SSE_R2)/(SSE_R2))/((Ne-1)/(1)) #(Nv,)
    F_R2_AllValues   = F_R2.copy()
    #F_R2[F_R2>F_MAX] = F_MAX
    
    # Kappa Computation
    Kappa_map        = F_R2 * weight_map
    Kappa_mask_arr   = ma.masked_array(Kappa_map,  mask=c_mask)
    KappW_mask_arr   = ma.masked_array(weight_map, mask=c_mask)
    Kappa            = Kappa_mask_arr.mean() / KappW_mask_arr.mean()
    #Kappa            = np.median(Kappa_mask_arr)/np.median(KappW_mask_arr)
    #Rho Computation
    Rho_map          = F_S0 * weight_map
    Rho_mask_arr     = ma.masked_array(Rho_map,    mask=c_mask)
    RhoW_mask_arr    = ma.masked_array(weight_map, mask=c_mask)
    Rho              = Rho_mask_arr.mean() / RhoW_mask_arr.mean()
    #Rho              = np.median(Rho_mask_arr)/np.median(RhoW_mask_arr)
    return {'Kappa':Kappa, 'Rho':Rho, 'Kappa_map':Kappa_map, 'Rho_map':Rho_map, 'FS0':F_S0_AllValues, 'FR2':F_R2_AllValues, 'cS0':coeffs_S0, 'cR2':coeffs_R2}
    
def characterize_components(origTS_pc, data_mean, tes, t2s, S0, mmix, ICA_maps, voxelwiseQA, 
                            Ncpus, ICA_maps_thr=1.0, discard_mask=None,writeOuts=False,outDir=None, outPrefix=None, mask=None, aff=None, head=None, Z_MAX=8, F_MAX=500):
    """
    This function computes kappa, rho and variance for each ICA component.
    
    Parameters:
    -----------
    origTS_pc:    original ME Timeseries in signal percent units for intra-cranial voxels (Nv,Ne,Nt)
    data_mean:    voxel-wise mean across time of the ME timeseries. (Nv,Ne)
    tes:          echo times (Ne,)
    t2s:          static T2* voxel-wise map (Nv,)
    S0:           static S0  voxel-wise map (Nv,)
    mmix:         ICA mixing matrix (Nc,Nt)
    ICA_maps:     ICA spatial maps  (Nv,Nc)
    voxelwiseQA:  voxel-wise QA map obtained from attempting a static TE-dependence fit to the means. It
                  can be used as part of the wiegths used duing the averaging. (Nv,)
    Ncpus:        number of available CPUs for multi-processing.
    ICA_maps_thr: threshold for selection of voxels entering the averaging for kappa and rho computation.
    discard_mask: voxels to be discarded becuase the static field generated erroneous S0 or T2* values.
    writeOuts:    flag to instruct the function to save additional files.
    outDir:       path where the program should write NIFTI and other datasets.
    outPrefix:    prefix for datasets written to disk.
    mask:         binary map indicating intra-cranial voxels (Nv,)
    aff:          affine needed by nibabel to write NIFTI datasets.
    head:         header needed by nibabel to write NIFTI datasets.
    Z_MAX:        maximum allowed Z-score in ICA maps.
    F_MAX:        maximum allowed F-stat  in R2 and S0 fits of the TE-dependence model.
    
    Results:
    --------
    features:     a 7xNc numpy array with features for each component. Column 0 is component ID,
                  column 1 is kappa, column 2 is rho, column 3 is explained variance, column 4 is 
                  max Fstat in the R2 fit, column 5 is max Fstat in the S0 fit, column 6 is the 
                  kappa/rho ratio. Components are sorted by variance (e.g., the way they came out of the ICA)
    """
    
    Nv,Ne,Nt = origTS_pc.shape
    Nc,_     = mmix.shape
    # If no discard_mask is provided, create one in which no voxels will be discarded during the
    # averaging.
    if discard_mask is None:
       discard_mask = np.zeros((Nv,))
    
    # Get ICA-component masks based on threshold
    ICA_maps_mask = ICA_maps.copy()
    ICA_maps_mask = (np.abs(ICA_maps)<ICA_maps_thr)  #(Nv,Nc)
    
    # Compute overall variance in the ICA data
    totalvar     = (ICA_maps**2).sum()               # Single Value
    
    # Compute beta maps (fits to each echo)
    # This is equivalent to running 3dDeconvolve with the mmix being the stim_files.
    # I tried this and gives exactly the same result
    # The 100 factor is there so that the b are kind of signal percent change
    beta       = 100*np.linalg.lstsq(mmix.T, (np.reshape(origTS_pc,(Nv*Ne,Nt))).T)[0].T 
    beta       = np.reshape(beta,(Nv,Ne,Nc))                         # (Nv,Ne,Nc)
    
    # Initialize results holder
    F_S0_maps  = np.zeros((Nv,Nc))
    F_R2_maps  = np.zeros((Nv,Nc))
    c_S0_maps  = np.zeros((Nv,Nc))
    c_R2_maps  = np.zeros((Nv,Nc))
    Kappa_maps = np.zeros((Nv,Nc))
    Kappa_masks= np.zeros((Nv,Nc))
    Rho_maps   = np.zeros((Nv,Nc))
    Rho_masks  = np.zeros((Nv,Nc))
    Weight_maps= np.zeros((Nv,Nc))
    varexp     = np.zeros(Nc)
    kappas     = np.zeros(Nc)
    rhos       = np.zeros(Nc)
    
    # Compute metrics per component
    X1 = np.ones((Ne,Nv)) 
    X2 = np.repeat((-tes/tes.mean())[:,np.newaxis].T,Nv,axis=0).T   #(Ne,Nv)
    print " +              Multi-process Characterize Components -> Ncpu = %d" % Ncpus
    pool   = Pool(processes=Ncpus)
    result = pool.map(_characterize_this_component, [{
                'c':          c,   'F_MAX':F_MAX,
                'X1':         X1,  'X2':X2, 
                'aff':        aff, 'head':head, 'outDir':outDir, 'outPrefix':outPrefix, 
                'mask':       mask,
                'B':          np.atleast_3d(beta)[:,:,c].transpose(),
                'weight_map': (ICA_maps[:,c]**2.)*voxelwiseQA,
                'c_mask':     ((discard_mask + ICA_maps_mask[:,c]) > 0.5),
                'writeOuts':  writeOuts
                } for c in np.arange(Nc)]) 
    
    features = np.zeros((Nc,7))
                            
    for c in range(Nc):
        Weight_maps[:,c] = (ICA_maps[:,c]**2.)*voxelwiseQA
        Kappa_maps[:,c]  = result[c]['Kappa_map']
        Rho_maps[:,c]    = result[c]['Rho_map']
        F_S0_maps[:,c]    = result[c]['FS0']
        F_R2_maps[:,c]    = result[c]['FR2']
        c_S0_maps[:,c]    = result[c]['cS0']
        c_R2_maps[:,c]    = result[c]['cR2']
        Kappa_masks[:,c]  = ((discard_mask + ICA_maps_mask[:,c]) > 0.5)
        Rho_masks[:,c]    = ((discard_mask + ICA_maps_mask[:,c]) > 0.5)        
        features[c,0]    = c
        features[c,1]    = result[c]['Kappa']
        features[c,2]    = result[c]['Rho']
        features[c,3]    = 100*((ICA_maps[:,c]**2).sum()/totalvar)
        FR2_mask_arr     = ma.masked_array(F_R2_maps[:,c], mask=Kappa_masks[:,c])
        features[c,4]    = FR2_mask_arr.max()
        FS0_mask_arr     = ma.masked_array(F_S0_maps[:,c], mask=Rho_masks[:,c])
        features[c,5]    = FS0_mask_arr.max()
        features[c,6]    = features[c,1] / features[c,2]
        
    niiwrite_nv(beta      , mask,outDir+outPrefix+'.chComp.Beta.nii',aff ,head)
    niiwrite_nv(F_S0_maps , mask,outDir+outPrefix+'.chComp.FS0.nii',aff ,head)
    niiwrite_nv(F_R2_maps , mask,outDir+outPrefix+'.chComp.FR2.nii',aff ,head)
    niiwrite_nv(c_S0_maps , mask,outDir+outPrefix+'.chComp.cS0.nii',aff ,head)
    niiwrite_nv(c_R2_maps , mask,outDir+outPrefix+'.chComp.cR2.nii',aff ,head)
    niiwrite_nv(Kappa_maps, mask,outDir+outPrefix+'.chComp.Kappa.nii',aff ,head)
    niiwrite_nv(abs(Kappa_masks-1),mask,outDir+outPrefix+'.chComp.Kappa_mask.nii',aff ,head)
    niiwrite_nv(Rho_maps,   mask,outDir+outPrefix+'.chComp.Rho.nii',aff ,head)
    niiwrite_nv(Weight_maps,mask,outDir+outPrefix+'.chComp.weightMaps.nii',aff ,head)   
    return features
    


def niiLoad(path):
	mepi_dset       = nib.load(path)
	data            = mepi_dset.get_data()
	aff             = mepi_dset.get_affine()
	head            = mepi_dset.get_header()
	head.extensions = []
	head.set_sform(head.get_sform(),code=1)
	return data,aff,head

def getSVDThreshold(svd_input,u,s,v,verb=False):
	""" 
	This function estimates the number of non-noise components in a PCA decomposition
	The algorithm used is from: Gavish and Donoho 2014 "The optiomal hard threshold for singular values is 4/sqrt(3)"
	Parameters:
   -----------
	svd_input: Input data to the SVD decomposition in an array (Nv.Nt)
	u:   U matrix from the SVD decomposition
	s:   Eigen-values from the SVD decomposition
	v:   VT matrix from the SVD decomposition. 
		  These three inputs are the outcomes of np.linalg.svd as they come out of that function.
	verb: If true, the function will print out some metrics about the results.
		  By default is set to False

   Results:
   --------
   Nc:  Number of non-gaussian noise components.
	"""
	m = svd_input.shape[1]
	n = svd_input.shape[0]
	beta = np.float(m)/np.float(n)
	omega = 0.56*beta**3 - 0.95*beta**2 + 1.82*beta + 1.43
	ymed = np.median(s)
	tau = omega*ymed
	shat = s * (s>tau)
	svd_inputhat = np.dot(u,np.dot(np.diag(shat),v))
	err = svd_input  - svd_inputhat
	relerr = np.linalg.norm(err,'fro')**2/np.linalg.norm(svd_input,'fro')**2
	Nc = np.where(shat>0)[0][-1] + 1
	if verb:
		print(" +              M = " + str(m) + ", N = " + str(n))
		print(" +              Beta = " + str(beta) + ", omega = " + str(omega))
		print(" +              Median y = " + str(ymed) + ", Tau = "+str(tau))
		print(" +              Number of components = " + str(Nc))
		print(" +              Signal power: " + str(1.0 - relerr))
		print(" +              Noise power: " + str(relerr))
	return Nc

def computeQA(data,tes,Ncpu,data_mean=None):
    """
    Simple function to compute the amount of variance in the data that is explained by
    the ME fit
    
    Parameters:
    -----------
    data: ME dataset (Nv,Ne,Nt)
    tes:  Echo times used during acquisition
    
    Returns:
    --------
    SSE:  Voxelwise mean across time of the Sum of Squared Errors for the TE fit. (Nv,)  
    rankSSE: ranked version from 0 to 100 of SSE (good for kappa computation)     (Nv,)
    """
    Nv,Ne,Nt        = data.shape
    if data_mean is None:
       data_mean = data.mean(axis=-1)
    data_demean     = data     - data_mean[:,:,np.newaxis]
    _, _, _, _, data_hat, _  = linearFit(data,tes,Ncpu,data_mean)
    data_hat_mean            = getMeanByPolyFit(data_hat,polort=7)
    data_demean_hat          = data_hat - data_hat_mean[:,:,np.newaxis]

    SSE             = ((data_demean - data_demean_hat)**2).sum(axis=-1).max(axis=-1)
    rankSSE         = 100.*rankdata(1./SSE)/Nv
    return SSE,rankSSE

def writeCompTable(out_dir,data_file, features, varexp, psel, Nt, sort_col):
    Nc,_ = features.shape
    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    #all_sorted_idxs   = (features[np.argsort(features[:,sort_col]),0][::-1]).astype('int')
    Ngood = (psel==1).sum()
    Nbad  = (psel==0).sum()
    midk  = []
    rej                 = np.where(psel==0)[0]
    #rej_sorting_feature = features[rej,:]
    #rej_sorted_idxs     = np.atleast_2d(rej[np.argsort(rej_sorting_feature)[::-1]])
    np.savetxt(out_dir+'rejected.txt',rej,fmt='%d',delimiter=',')
    acc                 = np.where(psel==1)[0]
    #acc_sorting_feature = features[acc,:]
    #acc_sorted_idxs     = np.atleast_2d(acc[np.argsort(acc_sorting_feature)[::-1]])
    np.savetxt(out_dir+'accepted.txt',acc,fmt='%d',delimiter=',')
    open(out_dir+'midk_rejected.txt','w').write(','.join([str(int(cc)) for cc in midk]))
    with open(out_dir+'comp_table.txt','w') as f:
        f.write("#\n#ME-ICA Component statistics table for: %s \n" % (data_file))
        f.write("#Run on %s \n" % (st)) 
        f.write("#\n")
        f.write("#Dataset variance explained by ICA (VEx): %.02f \n" %  ( varexp ) )
        f.write("#Total components generated by decomposition (TCo): %i \n" %  ( Nc ) )
        f.write("#No. accepted BOLD-like components, i.e. effective degrees of freedom for correlation (lower bound; DFe): %i\n" %  ( Ngood ) )
        f.write("#Total number of rejected components (RJn): %i\n" %  (Nbad) )
        f.write("#Nominal degress of freedom in denoised time series (..._medn.nii.gz; DFn): %i \n" %  (Nt-Nbad) )
        f.write("#ACC %s \t#Accepted BOLD-like components\n" % ','.join([str(int(cc)) for cc in acc]) )
        f.write("#REJ %s \t#Rejected non-BOLD components\n" % ','.join([str(int(cc)) for cc in rej]) )
        f.write("#MID    \t#Rejected R2*-weighted artifacts\n")
        f.write("#IGN    \t#Ignored components (kept in denoised time series)\n")
        f.write("#VEx  TCo   DFe   RJn   DFn   \n")
        f.write("##%.02f  %i %i %i %i \n" % (varexp,Nc,Ngood,Nbad,Nt-Nbad))
        f.write("#  comp  Kappa Rho   %%Var %%VarN	MaxR2	MaxS0	Ratio\n")
        idx = 0
        for i in range(Nc):
            f.write('%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n'%(features[i,0],features[i,1],features[i,2],features[i,3],features[i,3],features[i,4],features[i,5],features[i,6]))
            idx=idx+1
