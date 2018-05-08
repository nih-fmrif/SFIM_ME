import sys
import numpy.ma as ma
import numpy   as np
import nibabel as nib
import scipy.optimize as opt
import pandas         as pd
from scipy.stats import scoreatpercentile,zscore,skew,rankdata
import matplotlib.pyplot as plt
import time
import datetime
from scipy.special        import legendre
from sklearn.linear_model import LinearRegression
from multiprocessing      import cpu_count,Pool
from scipy.stats          import f,t, linregress
from math                 import ceil
from meUtils import niiLoad, niiwrite_nv, mask4MEdata


def numpy_weighted_median(data, weights=None):
    """Calculate the weighted median of an array/list using numpy."""
    import numpy as np
    if weights is None:
        return np.median(np.array(data).flatten())
    data, weights = np.array(data).flatten(), np.array(weights).flatten()
    if any(weights > 0):
        sorted_data, sorted_weights = map(np.array, zip(*sorted(zip(data, weights))))
        midpoint = 0.5 * sum(sorted_weights)
        if any(weights > midpoint):
            return (data[weights == np.max(weights)])[0]
        cumulative_weight = np.cumsum(sorted_weights)
        below_midpoint_index = np.where(cumulative_weight <= midpoint)[0][-1]
        if cumulative_weight[below_midpoint_index] == midpoint:
            return np.mean(sorted_data[below_midpoint_index:below_midpoint_index+2])
        return sorted_data[below_midpoint_index+1]

def weighted_median(data, weights=None):
    """Calculate the weighted median of a list."""
    if weights is None:
        return median(data)
    midpoint = 0.5 * sum(weights)
    if any([j > midpoint for j in weights]):
        return data[weights.index(max(weights))]
    if any([j > 0 for j in weights]):
        sorted_data, sorted_weights = zip(*sorted(zip(data, weights)))
        cumulative_weight = 0
        below_midpoint_index = 0
        while cumulative_weight <= midpoint:
            below_midpoint_index += 1
            cumulative_weight += sorted_weights[below_midpoint_index-1]
        cumulative_weight -= sorted_weights[below_midpoint_index-1]
        if cumulative_weight == midpoint:
            bounds = sorted_data[below_midpoint_index-2:below_midpoint_index]
            return sum(bounds) / float(len(bounds))
        return sorted_data[below_midpoint_index-1]


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
    FFT   = np.zeros(ts[:,int(Nt/2):].shape)
    for c in range(Nc):
        its      = ts[c,:]
        iFFT     = np.fft.fftshift(abs(np.fft.fft(its,axis=0)))
        FFT[c,:] = iFFT[int(Nt/2):]
        freq     = np.fft.fftshift(np.fft.fftfreq(Nt,float(TR)))
        freq     = freq[int(Nt/2):]
    return FFT, freq


def getelbow_meanvar(ks, val=False):
    """Elbow using mean/variance method - conservative

    This function was copied from another repository on May 7, 2018 (11282c5)
    https://github.com/ME-ICA/tedana/blob/master/tedana/interfaces/tedana.py

    Parameters
    ----------
    ks : array-like (presorted)
    val : bool, optional
        Default is False
    Returns
    -------
    array-like
        Either the elbow index (if val is True) or the values at the elbow
        index (if val is False)
    """
    # ks = np.sort(ks)[::-1]
    nk = len(ks)
    temp1 = [(ks[nk - 5 - ii - 1] > ks[nk - 5 - ii:nk].mean()
              + 2 * ks[nk - 5 - ii:nk].std())
             for ii in range(nk - 5)]
    ds = np.array(temp1[::-1], dtype=np.int)
    dsum = []
    c_ = 0
    for d_ in ds:
        c_ = (c_ + d_) * d_
        dsum.append(c_)
    e2 = np.argmax(np.array(dsum))
    elind = np.max([getelbow_linproj(ks), e2])

    if val:
        return ks[elind]
    else:
        return elind


def getelbow_linproj(ks, val=False):
    """Elbow using linear projection method - moderate

    This function was copied from another repository on May 7, 2018 (11282c5)
    https://github.com/ME-ICA/tedana/blob/master/tedana/interfaces/tedana.py

    Parameters
    ----------
    ks : array-like (presorted)
    val : bool, optional
        Default is False
    Returns
    -------
    array-like
        Either the elbow index (if val is True) or the values at the elbow
        index (if val is False)
    """
    # ks = np.sort(ks)[::-1]
    n_components = ks.shape[0]
    coords = np.array([np.arange(n_components), ks])
    p = coords - np.tile(np.reshape(coords[:, 0], (2, 1)), (1, n_components))
    b = p[:, -1]
    b_hat = np.reshape(b / np.sqrt((b ** 2).sum()), (2, 1))
    proj_p_b = p - np.dot(b_hat.T, p) * np.tile(b_hat, (1, n_components))
    d = np.sqrt((proj_p_b ** 2).sum(axis=0))
    k_min_ind = d.argmax()

    if val:
        return ks[k_min_ind]
    else:
        return k_min_ind


def getelbow_curvature(ks, val=False):
    """Elbow using curvature - aggressive

    This function was copied from another repository on May 7, 2018 (11282c5)
    https://github.com/ME-ICA/tedana/blob/master/tedana/interfaces/tedana.py

    Parameters
    ----------
    ks : array-like (presorted)
    val : bool, optional
        Default is False
    Returns
    -------
    array-like
        Either the elbow index (if val is True) or the values at the elbow
        index (if val is False)
    """
    # ks = np.sort(ks)[::-1]
    dKdt = ks[:-1] - ks[1:]
    dKdt2 = dKdt[:-1] - dKdt[1:]
    curv = np.abs((dKdt2 / (1 + dKdt[:-1]**2.) ** (3. / 2.)))
    curv[np.isnan(curv)] = -1 * 10**6
    maxcurv = np.argmax(curv) + 2

    if val:
        return(ks[maxcurv])
    else:
        return maxcurv


def getelbow(ks, Method='linproj', val=False):
    """
    Take a sorted series of values, stored in ks and find an inflection
    point in those values (the elbow). There are several methods for
    finding the elbow with some more agressive than others. More aggressive
    methods are more likely to set the elbow that fewer values are larger
    than the elbow

    Parameters
    ----------
    ks : array-like list of Values
    Method: Options for the elbow finding method.
        'meanvar'   Mean/Variance method - conservative
        'linproj'   Linear projection method - moderate [DEFAULT]
        'curvature' Curvature method - agressive
    val : bool, optional
          Used to define if the elbow value or index is returned
          Default is False (return the index)

    Returns
    -------
    elbowval is the index for the elbow from the shorted ks array
    (if val is True) or the values at the elbow index (if val is False)
    """"
    ks = np.sort(ks)[::-1]

    if Method.lower() == 'meanvar':
        elbowval = getelbow_meanvar(ks, val)
    elif Method.lower() == 'linproj':
        elbowval = getelbow_linproj(ks, val)
    elif Method.lower() == 'curvature':
        elbowval = getelbow_curvature(ks, val)
    else:
        print("++ Error: Method provided to getelbow not found: %s" % Method)

    return elbowval

# def getelbow(ks):
#    nc = ks.shape[0]
#    coords = np.array([np.arange(nc),ks])
#    p  = coords - np.tile(np.reshape(coords[:,0],(2,1)),(1,nc))
#    b  = p[:,-1]
#    b_hat = np.reshape(b/np.sqrt((b**2).sum()),(2,1))
#    proj_p_b = p - np.dot(b_hat.T,p)*np.tile(b_hat,(1,nc))
#    d = np.sqrt((proj_p_b**2).sum(axis=0))
#    k_min_ind = d.argmax()
#    k_min  = ks[k_min_ind]
#    return k_min_ind


def SelectGoodComponents(fica_feats, SelectionCriteria='method1',
                         krRatio=None):
    """
    This function uses the ICA component features to decide which components
    are more BOLD-like and should be make_optcom

    Parameters
    ----------
    fica_feats: Pandas data structure with information for each component
                It is calculated in characterize_components()
    SelectionCriteria: There are several algorithms for selecting criteria.
        This variable should be assigned the name of the desired algorithms
        The options are:
        'Method1': A basic elbow-based criteria
            (Method1 is a placeholder name until I think of a logical naming
            structure for the various selection options)

    Returns
    -------
    fica_psel: The idices in arrays for the accepted components
    accepted: The component IDs of the accepted Components
    rejected: The component IDs of the rejected Components
    """
    if SelectionCriteria.lower() == 'method1':
        KappaCutOff = getelbow(fica_feats['Kappa'], Method='linproj', 'True')
        RhoCutOff = getelbow(fica_feats['Rho'], Method='curvature', 'True')
        fica_psel = ((fica_feats['Kappa'] >= KappaCutOff) or
                     (fica_feats['Rho'] < RhoCutOff)).get_values()
        accepted = fica_feats['cID'][fica_psel].get_values().astype(int)
        rejected = np.setdiff1d(fica_feats['cID'].get_values().astype(int),
                                accepted)

    return fica_psel, accepted, rejected


def getMeanByPolyFit(data, polort=4):
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
    Nv, Ne, Nt = data.shape
    aux = np.reshape(data.copy(), (Nv*Ne, Nt))
    drift = np.zeros((Nt, polort))
    x = np.linspace(-1, 1-2/Nt, Nt)
    for n in range(polort):
        drift[:, n] = np.polyval(legendre(n), x)
    # 2. Fit polynomials to residuals
    linRegObj = LinearRegression(normalize=False, fit_intercept=False)
    linRegObj.fit(drift, aux.T)
    mean = np.reshape(linRegObj.coef_[:, 0], (Nv, Ne))
    return mean

def make_optcom(data, t2s, tes):
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
   print(" + INFO [make_static_maps_opt]: Initial conditions [So=%i, T2s=%i]" % (So_init, T2s_init))
   print(" + INFO [make_static_maps_opt]: Bounds So=[%i,%i] & T2s=[%i,%i]" % (So_min, So_max, T2s_min, T2s_max))
   print(" + INFO [make_static_maps_opt]: Optimizer = %s" % Optimizer)

   Nv,Ne    = data_mean.shape
   S0       = np.zeros(Nv,)
   t2s      = np.zeros(Nv,)
   badFits  = np.zeros(Nv,)
   fiterror = np.zeros(Nv,)

   print(" +              Multi-process Static Map Fit -> Ncpu = %d" % Ncpu)
   pool   = Pool(processes=Ncpu)
   result = pool.map(make_static_opt_perVoxel, [{'data_mean':data_mean[v,:],'tes':tes,'So_init':So_init,'So_max':So_max,'So_min':So_min,'T2s_init':T2s_init,'T2s_max':T2s_max,'T2s_min':T2s_min,'Optimizer':Optimizer} for v in np.arange(Nv)])
   for v in np.arange(Nv):
     S0[v]  = result[v]['v_S0']
     t2s[v] = result[v]['v_t2s']
     fiterror[v] = result[v]['v_fiterror']
     badFits[v]  = result[v]['v_badFit']
   print(" + INFO [make_static_maps_opt]: Number of Voxels with errors: %i" % badFits.sum())
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
    Z_MAX      = item['Z_MAX']
    weight_map = item['weight_map']
    c_mask     = item['c_mask']
    writeOuts  = item['writeOuts']
    doFM       = item['doFM']
    Ne, Nv     = B.shape
    Kappa_mask = c_mask.copy()
    Rho_mask   = c_mask.copy()
    # Write the fits for the differente echoes into file (very useful for debugging purposes)
    if writeOuts:
           niiwrite_nv(B.T,mask,outDir+outPrefix+'.chComp.EXTRA.Beta'+str(c).zfill(3)+'.nii',aff,head)
    # S0 Model
    coeffs_S0        = (B*X1).sum(axis=0)/(X1**2).sum(axis=0)        #(Nv,)
    estima_S0        = X1*np.tile(coeffs_S0,(Ne,1))
    SSR_S0           = (estima_S0**2).sum(axis=0)
    SSE_S0           = ((B-estima_S0)**2).sum(axis=0)
    dofSSR_S0        = 1
    dofSSE_S0        = Ne -1
    F_S0             = (SSR_S0/dofSSR_S0) / (SSE_S0/dofSSE_S0) #(Nv,)
    p_S0             = 1-f.cdf(F_S0,1,Ne-1)
    F_S0_mask        = (p_S0 < 0.05)
    F_S0[F_S0>F_MAX] = F_MAX
    F_S0_AllValues   = F_S0.copy()
    if writeOuts:
        niiwrite_nv((X1*np.tile(coeffs_S0,(Ne,1))).T,mask,outDir+outPrefix+'.chComp.EXTRA.S0Fit'+str(c).zfill(3)+'.nii',aff,head)

   # R2 Model
    # If beta_i = alpha_i*TE/mean(TE) + error --> To solve this RTO model (Regression Through the Origin), it is possible to obtain
    # the alpha_i that provides the best fit (in a least-squares fashion) by simply using the equation below
    # Please look at: Eishenhouer JG "Regression throught the origin" Technical Statistics (25):3, 2003
    coeffs_R2        = (B*X2).sum(axis=0)/(X2**2).sum(axis=0)        #(Nv,)
    estima_R2        = X2*np.tile(coeffs_R2,(Ne,1))
    SSR_R2           = (estima_R2**2).sum(axis=0)
    SSE_R2           = ((B-estima_R2)**2).sum(axis=0)
    dofSSR_R2        = 1
    dofSSE_R2        = Ne -1
    F_R2             = (SSR_R2/dofSSR_R2) / (SSE_R2/dofSSE_R2) #(Nv,)
    p_R2             = 1-f.cdf(F_R2,1,Ne-1)
    F_R2_mask        = (p_R2 < 0.05)
    F_R2[F_R2>F_MAX] = F_MAX
    F_R2_AllValues   = F_R2.copy()

    if writeOuts:
        niiwrite_nv((X2*np.tile(coeffs_R2,(Ne,1))).T,mask,outDir+outPrefix+'.chComp.EXTRA.R2Fit'+str(c).zfill(3)+'.nii',aff,head)

    # Mask for spatial averaging
    Kappa_mask   = np.logical_and(c_mask, np.logical_or(F_R2_mask, F_S0_mask)) #(Nv,)
    Rho_mask     = Kappa_mask.copy() #(Nv,)

    Kappa_map  = F_R2 * weight_map
    Kappa_map  = Kappa_map/(weight_map[Kappa_mask].mean())
    Kappa_map  = Kappa_map * Kappa_mask           # weigths from the voxels entering the computation
    Kappa = np.mean(Kappa_map[Kappa_mask])

    # Rho Computation
    Rho_map  = F_S0 * weight_map
    Rho_map  = Rho_map/(weight_map[Rho_mask].mean())
    Rho_map  = Rho_map * Rho_mask           # weigths from the voxels entering the computation
    Rho = np.mean(Rho_map[Rho_mask])


    # EXTRA CODE
    if doFM==True:
        FullModel_Slope = np.zeros((Nv,))
        FullModel_Inter = np.zeros((Nv,))
        FullModel_p     = np.zeros((Nv,))
        FullModel_r     = np.zeros((Nv,))
        FullModel_Slope_err = np.zeros((Nv,))
        FullModel_Slope_T   = np.zeros((Nv,))
        FullModel_Slope_p   = np.zeros((Nv,))
        for v in range(Nv):
            FullModel_Slope[v], FullModel_Inter[v], FullModel_r[v], FullModel_p[v], FullModel_Slope_err[v] = linregress(X2[:,v],B[:,v])
            FullModel_Slope_T[v] = FullModel_Slope[v]/FullModel_Slope_err[v]
            FullModel_Slope_p[v] = 2*(1-t.cdf(np.abs(FullModel_Slope_T[v]),Ne-2))
        return {'Kappa':Kappa, 'Rho':Rho, 'Kappa_map':Kappa_map, 'Rho_map':Rho_map, 'FS0':F_S0_AllValues, 'FR2':F_R2_AllValues, 'cS0':coeffs_S0, 'cR2':coeffs_R2,
            'Kappa_mask':Kappa_mask,'Rho_mask':Rho_mask, 'F_R2_mask':F_R2_mask, 'F_S0_mask':F_S0_mask, 'pR2':p_R2, 'pS0':p_S0,
            'FM_Slope':FullModel_Slope, 'FM_Inter':FullModel_Inter, 'FM_p':FullModel_p, 'FM_r':FullModel_r,
            'FM_Slope_err':FullModel_Slope_err, 'FM_Slope_T':FullModel_Slope_T, 'FM_Slope_p':FullModel_Slope_p}
    else:
        return {'Kappa':Kappa, 'Rho':Rho, 'Kappa_map':Kappa_map, 'Rho_map':Rho_map, 'FS0':F_S0_AllValues, 'FR2':F_R2_AllValues, 'cS0':coeffs_S0, 'cR2':coeffs_R2,
            'Kappa_mask':Kappa_mask,'Rho_mask':Rho_mask, 'F_R2_mask':F_R2_mask, 'F_S0_mask':F_S0_mask, 'pR2':p_R2, 'pS0':p_S0}

def characterize_components(origTS_pc, data_mean, tes, t2s, S0, mmix, ICA_maps, voxelwiseQA,
                            Ncpus, ICA_maps_thr=0.0, discard_mask=None,writeOuts=False,outDir=None, outPrefix=None, mask=None, aff=None, head=None, Z_MAX=8, F_MAX=500, doFM=False):
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
       discard_mask = np.zeros((Nv,), dtype=bool)

    # Get ICA-component masks based on threshold
    ICA_maps_mask = np.zeros(ICA_maps.shape, dtype=bool)   #### THIS WAS SIZE BEFORE..... SHOULD DEFINITELY CHECK
    ICA_maps_mask = np.logical_and((np.abs(ICA_maps)>ICA_maps_thr), discard_mask[:,np.newaxis] )
    niiwrite_nv(ICA_maps_mask, mask,outDir+outPrefix+'.ICA.Zmaps.mask.nii',aff ,head)

    # Compute overall variance in the ICA data
    totalvar     = (ICA_maps**2).sum()               # Single Value

    # Compute beta maps (fits to each echo)
    # This is equivalent to running 3dDeconvolve with the mmix being the stim_files.
    # I tried this and gives exactly the same result
    # The 100 factor is there so that the b are kind of signal percent change
    beta       = 100*np.linalg.lstsq(mmix.T, (np.reshape(origTS_pc,(Nv*Ne,Nt))).T)[0].T
    beta       = np.reshape(beta,(Nv,Ne,Nc))   ### <----------------------------------------  MAYBE PUT HERE AN ABS                      # (Nv,Ne,Nc)

    # Initialize results holder
    F_S0_maps  = np.zeros((Nv,Nc))
    F_R2_maps  = np.zeros((Nv,Nc))
    F_S0_masks = np.zeros((Nv,Nc), dtype=bool)
    F_R2_masks = np.zeros((Nv,Nc), dtype=bool)
    c_S0_maps  = np.zeros((Nv,Nc))
    c_R2_maps  = np.zeros((Nv,Nc))
    p_S0_maps  = np.zeros((Nv,Nc))
    p_R2_maps  = np.zeros((Nv,Nc))
    Kappa_maps = np.zeros((Nv,Nc))
    Kappa_masks= np.zeros((Nv,Nc), dtype=bool)
    Rho_maps   = np.zeros((Nv,Nc))
    Rho_masks  = np.zeros((Nv,Nc), dtype=bool)
    Weight_maps= np.zeros((Nv,Nc))
    varexp     = np.zeros(Nc)
    kappas     = np.zeros(Nc)
    rhos       = np.zeros(Nc)
    if doFM:
        FM_Slope_map = np.zeros((Nv,Nc))
        FM_Inter_map = np.zeros((Nv,Nc))
        FM_p_map     = np.zeros((Nv,Nc))
        FM_r_map     = np.zeros((Nv,Nc))
        FM_Slope_err_map = np.zeros((Nv,Nc))
        FM_Slope_T_map   = np.zeros((Nv,Nc))
        FM_Slope_p_map   = np.zeros((Nv,Nc))
    # Compute metrics per component
    print("Hello -------------->")
    print("beta.%s" % str(beta.shape))
    print("ICA_maps.%s" % str(ICA_maps.shape))
    print("ICA_maps_mask.%s" % str(ICA_maps_mask.shape))
    X1 = np.ones((Ne,Nv))
    X2 = np.repeat((tes/tes.mean())[:,np.newaxis].T,Nv,axis=0).T   #<--------------------   MAYBE NEEDS A MINUS SIGN   (Ne,Nv)
    print("X1.%s" % str(X1.shape))
    print("X2.%s" % str(X2.shape))
    print(" +              Multi-process Characterize Components -> Ncpu = %d" % Ncpus)
    pool   = Pool(processes=Ncpus)
    result = pool.map(_characterize_this_component, [{
                'c':          c,   'F_MAX':F_MAX, 'Z_MAX':Z_MAX,
                'X1':         X1,  'X2':X2,
                'aff':        aff, 'head':head, 'outDir':outDir, 'outPrefix':outPrefix,
                'mask':       mask,
                'B':          np.atleast_3d(beta)[:,:,c].transpose(),
                'weight_map': (ICA_maps[:,c]**2.)*voxelwiseQA,
                'c_mask':     ICA_maps_mask[:,c],
                'writeOuts':  writeOuts,
                'doFM': doFM
                } for c in np.arange(Nc)])
    feat_names = ['cID','Kappa','Rho','Var','maxFR2','maxFS0','K/R', 'maxZICA','NvZmask','NvFR2mask','NvFS0mask','NvKapMask','NvRhoMask','Dan']
    feat_vals  = np.zeros((Nc,len(feat_names)))

    for c in range(Nc):
        Weight_maps[:,c] = (ICA_maps[:,c]**2.)*voxelwiseQA
        Kappa_maps[:,c]  = result[c]['Kappa_map']
        Rho_maps[:,c]    = result[c]['Rho_map']
        F_S0_maps[:,c]    = result[c]['FS0']
        F_R2_maps[:,c]    = result[c]['FR2']
        c_S0_maps[:,c]    = result[c]['cS0']
        c_R2_maps[:,c]    = result[c]['cR2']
        p_S0_maps[:,c]    = result[c]['pS0']
        p_R2_maps[:,c]    = result[c]['pR2']
        F_S0_masks[:,c]   = result[c]['F_S0_mask']
        F_R2_masks[:,c]   = result[c]['F_R2_mask']
        Kappa_masks[:,c]  = result[c]['Kappa_mask']
        Rho_masks[:,c]    = result[c]['Rho_mask']
        if doFM:
            FM_Slope_map[:,c] = result[c]['FM_Slope']
            FM_Inter_map[:,c] = result[c]['FM_Inter']
            FM_p_map[:,c]     = result[c]['FM_p']
            FM_r_map[:,c]     = result[c]['FM_r']
            FM_Slope_err_map[:,c] = result[c]['FM_Slope_err']
            FM_Slope_T_map[:,c]   = result[c]['FM_Slope_T']
            FM_Slope_p_map[:,c]   = result[c]['FM_Slope_p']
        feat_vals[c,0]    = c
        feat_vals[c,1]    = result[c]['Kappa']
        feat_vals[c,2]    = result[c]['Rho']
        feat_vals[c,3]    = 100*((ICA_maps[:,c]**2).sum()/totalvar)


        FR2_mask_arr     = ma.masked_array(F_R2_maps[:,c], mask=np.logical_not(Kappa_masks[:,c])) #abs(Kappa_masks[:,c]-1))
        feat_vals[c,4]    = FR2_mask_arr.max()

        FS0_mask_arr     = ma.masked_array(F_S0_maps[:,c], mask=np.logical_not(Rho_masks[:,c])) #abs(Rho_masks[:,c]-1))
        feat_vals[c,5]    = FS0_mask_arr.max()

        feat_vals[c,6]    = feat_vals[c,1] / feat_vals[c,2]

        ZICA_mask_arr    = ma.masked_array(ICA_maps[:,c], mask=np.logical_not(ICA_maps_mask[:,c])) #abs(ICA_maps_mask[:,c]-1))
        feat_vals[c,7]    = ZICA_mask_arr.max()

        feat_vals[c,8]    = ICA_maps_mask[:,c].sum()
        feat_vals[c,9]    = F_R2_masks[:,c].sum()
        feat_vals[c,10]   = F_S0_masks[:,c].sum()
        feat_vals[c,11]   = Kappa_masks[:,c].sum()
        feat_vals[c,12]   = Rho_masks[:,c].sum()

        # DAN METRIC
        if doFM:
            aux_mask      = np.logical_and((FM_p_map[:,c]<0.05),(FM_Slope_p_map[:,c]<0.05))
            aux_numerator = Weight_maps[aux_mask,c].sum()
            aux_denominat = Weight_maps[:,c].sum()
            aux_metric = aux_numerator/aux_denominat
            print("[%d] -> aux_mask%s | aux_numerator%s | DF=%f" % (c,str(aux_mask.shape),str(Weight_maps[aux_mask,c].shape),aux_metric))

    df_feats = pd.DataFrame(data=feat_vals,columns=feat_names)
    df_feats.to_csv(outDir+outPrefix+'.DF.csv')
    df_feats['cID'] = df_feats['cID'].astype(int)

    niiwrite_nv(beta      , mask,outDir+outPrefix+'.chComp.Beta.nii',aff ,head)
    niiwrite_nv(F_S0_maps , mask,outDir+outPrefix+'.chComp.FS0.nii',aff ,head)
    niiwrite_nv(F_R2_maps , mask,outDir+outPrefix+'.chComp.FR2.nii',aff ,head)
    niiwrite_nv(F_S0_masks, mask,outDir+outPrefix+'.chComp.FS0.mask.nii',aff ,head)
    niiwrite_nv(F_R2_masks, mask,outDir+outPrefix+'.chComp.FR2.mask.nii',aff ,head)
    niiwrite_nv(c_S0_maps , mask,outDir+outPrefix+'.chComp.cS0.nii',aff ,head)
    niiwrite_nv(c_R2_maps , mask,outDir+outPrefix+'.chComp.cR2.nii',aff ,head)
    niiwrite_nv(p_S0_maps , mask,outDir+outPrefix+'.chComp.pS0.nii',aff ,head)
    niiwrite_nv(p_R2_maps , mask,outDir+outPrefix+'.chComp.pR2.nii',aff ,head)

    niiwrite_nv(Kappa_maps,  mask,outDir+outPrefix+'.chComp.Kappa.nii',aff ,head)
    niiwrite_nv(Kappa_masks, mask,outDir+outPrefix+'.chComp.Kappa_mask.nii',aff ,head)
    niiwrite_nv(Rho_masks,   mask,outDir+outPrefix+'.chComp.Rho_mask.nii',aff ,head)
    niiwrite_nv(Rho_maps,    mask,outDir+outPrefix+'.chComp.Rho.nii',aff ,head)
    niiwrite_nv(Weight_maps, mask,outDir+outPrefix+'.chComp.weightMaps.nii',aff ,head)
    if doFM:
        niiwrite_nv(FM_Slope_map , mask,outDir+outPrefix+'.chComp.FM.Slope.nii',aff ,head)
        niiwrite_nv(FM_Inter_map , mask,outDir+outPrefix+'.chComp.FM.Inter.nii',aff ,head)
        niiwrite_nv(FM_p_map , mask,outDir+outPrefix+'.chComp.FM.p.nii',aff ,head)
        niiwrite_nv(FM_r_map , mask,outDir+outPrefix+'.chComp.FM.r.nii',aff ,head)
        niiwrite_nv(FM_Slope_err_map , mask,outDir+outPrefix+'.chComp.FM.Slope.err.nii',aff ,head)
        niiwrite_nv(FM_Slope_T_map , mask,outDir+outPrefix+'.chComp.FM.Slope.T.nii',aff ,head)
        niiwrite_nv(FM_Slope_p_map , mask,outDir+outPrefix+'.chComp.FM.Slope.p.nii',aff ,head)
    return df_feats


# SEPARATE ECHOES BUSINESS
def _characterize_this_component_se(item):
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
    Z_MAX      = item['Z_MAX']
    weight_map = item['weight_map']
    c_mask     = item['c_mask']
    writeOuts  = item['writeOuts']
    doFM       = item['doFM']
    doMedian  = item['doMedian']
    Ne, Nv     = B.shape
    Kappa_mask = c_mask.copy()
    Rho_mask   = c_mask.copy()
    # Write the fits for the differente echoes into file (very useful for debugging purposes)
    if writeOuts:
           niiwrite_nv(B.T,mask,outDir+outPrefix+'.chComp.EXTRA.Beta'+str(c).zfill(3)+'.nii',aff,head)
    # S0 Model
    coeffs_S0        = (B*X1).sum(axis=0)/(X1**2).sum(axis=0)        #(Nv,)
    estima_S0        = X1*np.tile(coeffs_S0,(Ne,1))
    SSR_S0           = (estima_S0**2).sum(axis=0)
    SSE_S0           = ((B-estima_S0)**2).sum(axis=0)
    dofSSR_S0        = 1
    dofSSE_S0        = Ne -1
    F_S0             = (SSR_S0/dofSSR_S0) / (SSE_S0/dofSSE_S0) #(Nv,)
    p_S0             = 1-f.cdf(F_S0,1,Ne-1)
    F_S0_mask        = (p_S0 < 0.05)
    F_S0[F_S0>F_MAX] = F_MAX
    F_S0_AllValues   = F_S0.copy()
    if writeOuts:
        niiwrite_nv((X1*np.tile(coeffs_S0,(Ne,1))).T,mask,outDir+outPrefix+'.chComp.EXTRA.S0Fit'+str(c).zfill(3)+'.nii',aff,head)

    # R2 Model
    # If beta_i = alpha_i*TE/mean(TE) + error --> To solve this RTO model (Regression Through the Origin), it is possible to obtain
    # the alpha_i that provides the best fit (in a least-squares fashion) by simply using the equation below
    # Please look at: Eishenhouer JG "Regression throught the origin" Technical Statistics (25):3, 2003
    coeffs_R2        = (B*X2).sum(axis=0)/(X2**2).sum(axis=0)        #(Nv,)
    estima_R2        = X2*np.tile(coeffs_R2,(Ne,1))
    SSR_R2           = (estima_R2**2).sum(axis=0)
    SSE_R2           = ((B-estima_R2)**2).sum(axis=0)
    dofSSR_R2        = 1
    dofSSE_R2        = Ne -1
    F_R2             = (SSR_R2/dofSSR_R2) / (SSE_R2/dofSSE_R2) #(Nv,)
    p_R2             = 1-f.cdf(F_R2,1,Ne-1)
    F_R2_mask        = (p_R2 < 0.05)
    F_R2[F_R2>F_MAX] = F_MAX
    F_R2_AllValues   = F_R2.copy()

    if writeOuts:
        niiwrite_nv((X2*np.tile(coeffs_R2,(Ne,1))).T,mask,outDir+outPrefix+'.chComp.EXTRA.R2Fit'+str(c).zfill(3)+'.nii',aff,head)

    # Mask for spatial averaging
    Kappa_mask   = np.logical_and(c_mask, np.logical_or(F_R2_mask, F_S0_mask)) #(Nv,)
    Rho_mask     = Kappa_mask.copy() #(Nv,)

    # Kappa Computation
    Kappa_map  = F_R2 * weight_map
    Kappa_map  = Kappa_map/(weight_map[Kappa_mask].mean())
    Kappa_map  = Kappa_map * Kappa_mask           # weigths from the voxels entering the computation
    if doMedian:
        Kappa = np.median(Kappa_map[Kappa_mask])
    else:
        Kappa = np.mean(Kappa_map[Kappa_mask])

    # Rho Computation
    Rho_map  = F_S0 * weight_map
    Rho_map  = Rho_map/(weight_map[Rho_mask].mean())
    Rho_map  = Rho_map * Rho_mask           # weigths from the voxels entering the computation
    if doMedian:
        Rho = np.median(Rho_map[Rho_mask])
    else:
        Rho = np.mean(Rho_map[Rho_mask])

    # EXTRA CODE
    if doFM==True:
        FullModel_Slope = np.zeros((Nv,))
        FullModel_Inter = np.zeros((Nv,))
        FullModel_p     = np.zeros((Nv,))
        FullModel_r     = np.zeros((Nv,))
        FullModel_Slope_err = np.zeros((Nv,))
        FullModel_Slope_T   = np.zeros((Nv,))
        FullModel_Slope_p   = np.zeros((Nv,))
        for v in range(Nv):
            FullModel_Slope[v], FullModel_Inter[v], FullModel_r[v], FullModel_p[v], FullModel_Slope_err[v] = linregress(X2[:,v],B[:,v])
            FullModel_Slope_T[v] = FullModel_Slope[v]/FullModel_Slope_err[v]
            FullModel_Slope_p[v] = 2*(1-t.cdf(np.abs(FullModel_Slope_T[v]),Ne-2))
        return {'Kappa':Kappa, 'Rho':Rho, 'Kappa_map':Kappa_map, 'Rho_map':Rho_map, 'FS0':F_S0_AllValues, 'FR2':F_R2_AllValues, 'cS0':coeffs_S0, 'cR2':coeffs_R2,
            'Kappa_mask':Kappa_mask,'Rho_mask':Rho_mask, 'F_R2_mask':F_R2_mask, 'F_S0_mask':F_S0_mask, 'pR2':p_R2, 'pS0':p_S0,
            'FM_Slope':FullModel_Slope, 'FM_Inter':FullModel_Inter, 'FM_p':FullModel_p, 'FM_r':FullModel_r,
            'FM_Slope_err':FullModel_Slope_err, 'FM_Slope_T':FullModel_Slope_T, 'FM_Slope_p':FullModel_Slope_p}
    else:
        return {'Kappa':Kappa, 'Rho':Rho, 'Kappa_map':Kappa_map, 'Rho_map':Rho_map, 'FS0':F_S0_AllValues, 'FR2':F_R2_AllValues, 'cS0':coeffs_S0, 'cR2':coeffs_R2,
            'Kappa_mask':Kappa_mask,'Rho_mask':Rho_mask, 'F_R2_mask':F_R2_mask, 'F_S0_mask':F_S0_mask, 'pR2':p_R2, 'pS0':p_S0}

def characterize_components_se(origTS_pc, data_mean, tes, t2s, S0, mmix, ICA_maps, voxelwiseQA,
                            Ncpus, ICA_maps_thr=0.0, discard_mask=None,writeOuts=False,outDir=None, outPrefix=None, mask=None, aff=None, head=None, Z_MAX=8, F_MAX=500, doFM=False, doMedian=False):
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
    discard_mask: voxels to be discarded because the static field generated erroneous S0 or T2* values.
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
       discard_mask = np.zeros((Nv*Ne,), dtype=bool)
    else:
        discard_mask = np.tile(discard_mask,Ne)
    # Get ICA-component masks based on threshold
    ICA_maps_mask = np.zeros(ICA_maps.shape, dtype=bool) #(Nv*Ne,Nc)
    ICA_maps_mask = np.logical_and((np.abs(ICA_maps)>ICA_maps_thr), discard_mask[:,np.newaxis] )
    niiwrite_nv(np.reshape(ICA_maps_mask,(Nv,Ne,Nc),order='F'), mask,outDir+outPrefix+'.ICA.Zmaps.mask.nii',aff ,head)

    # Compute overall variance in the ICA data
    totalvar     = (ICA_maps**2).sum()               # Single Value

    # Compute beta maps (fits to each echo)
    # This is equivalent to running 3dDeconvolve with the mmix being the stim_files.
    # I tried this and gives exactly the same result
    # The 100 factor is there so that the b are kind of signal percent change
    beta       = 100*np.linalg.lstsq(mmix.T, (np.reshape(origTS_pc,(Nv*Ne,Nt))).T)[0].T
    beta       = np.reshape(beta,(Nv,Ne,Nc))   #(Nv,Ne,Nc)## <----------------------------------------  MAYBE PUT HERE AN ABS

    #beta       = np.reshape(ICA_maps,(Nv,Ne,Nc),order='F') #(Nv*Ne,Nc)
    oc_ICA_maps= make_optcom(beta,t2s,tes)
    oc_ICA_maps_mask = np.sum(np.reshape(ICA_maps_mask,(Nv,Ne,Nc),order='F'),axis=1)>(ceil(Ne/2.)) #(Nv,Nc)
    niiwrite_nv(oc_ICA_maps_mask, mask,outDir+outPrefix+'.ICA.Zmaps.maskOC.nii',aff ,head) #(Nv,Nc)
    # Initialize results holder
    F_S0_maps  = np.zeros((Nv,Nc))
    F_R2_maps  = np.zeros((Nv,Nc))
    F_S0_masks = np.zeros((Nv,Nc), dtype=bool)
    F_R2_masks = np.zeros((Nv,Nc), dtype=bool)
    c_S0_maps  = np.zeros((Nv,Nc))
    c_R2_maps  = np.zeros((Nv,Nc))
    p_S0_maps  = np.zeros((Nv,Nc))
    p_R2_maps  = np.zeros((Nv,Nc))
    Kappa_maps = np.zeros((Nv,Nc))
    Kappa_masks= np.zeros((Nv,Nc), dtype=bool)
    Rho_maps   = np.zeros((Nv,Nc))
    Rho_masks  = np.zeros((Nv,Nc), dtype=bool)
    Weight_maps= np.zeros((Nv,Nc))
    varexp     = np.zeros(Nc)
    kappas     = np.zeros(Nc)
    rhos       = np.zeros(Nc)
    if doFM:
        FM_Slope_map = np.zeros((Nv,Nc))
        FM_Inter_map = np.zeros((Nv,Nc))
        FM_p_map     = np.zeros((Nv,Nc))
        FM_r_map     = np.zeros((Nv,Nc))
        FM_Slope_err_map = np.zeros((Nv,Nc))
        FM_Slope_T_map   = np.zeros((Nv,Nc))
        FM_Slope_p_map   = np.zeros((Nv,Nc))
    # Compute metrics per component
    X1 = np.ones((Ne,Nv)) #(Ne,Nv)
    X2 = np.repeat((tes/tes.mean())[:,np.newaxis].T,Nv,axis=0).T   #<--------------------   MAYBE NEEDS A MINUS SIGN   (Ne,Nv)
    print(" +              Multi-process Characterize Components -> Ncpu = %d" % Ncpus)
    pool   = Pool(processes=Ncpus)
    result = pool.map(_characterize_this_component_se, [{
                'c':          c,   'F_MAX':F_MAX, 'Z_MAX':Z_MAX,
                'X1':         X1,  'X2':X2,
                'aff':        aff, 'head':head, 'outDir':outDir, 'outPrefix':outPrefix,
                'mask':       mask,
                'B':          np.atleast_3d(beta)[:,:,c].transpose(),
                'weight_map': (oc_ICA_maps[:,c]**2.)*voxelwiseQA,
                'c_mask':     oc_ICA_maps_mask[:,c],
                'writeOuts':  writeOuts,
                'doFM': doFM, 'doMedian':doMedian
                } for c in np.arange(Nc)])
    feat_names = ['cID','Kappa','Rho','Var','maxFR2','maxFS0','K/R', 'maxZICA','NvZmask','NvFR2mask','NvFS0mask','NvKapMask','NvRhoMask','Dan']
    feat_vals  = np.zeros((Nc,len(feat_names)))

    for c in range(Nc):
        Weight_maps[:,c] = (oc_ICA_maps[:,c]**2.)*voxelwiseQA
        Kappa_maps[:,c]  = result[c]['Kappa_map']
        Rho_maps[:,c]    = result[c]['Rho_map']
        F_S0_maps[:,c]    = result[c]['FS0']
        F_R2_maps[:,c]    = result[c]['FR2']
        c_S0_maps[:,c]    = result[c]['cS0']
        c_R2_maps[:,c]    = result[c]['cR2']
        p_S0_maps[:,c]    = result[c]['pS0']
        p_R2_maps[:,c]    = result[c]['pR2']
        F_S0_masks[:,c]   = result[c]['F_S0_mask']
        F_R2_masks[:,c]   = result[c]['F_R2_mask']
        Kappa_masks[:,c]  = result[c]['Kappa_mask']
        Rho_masks[:,c]    = result[c]['Rho_mask']
        if doFM:
            FM_Slope_map[:,c] = result[c]['FM_Slope']
            FM_Inter_map[:,c] = result[c]['FM_Inter']
            FM_p_map[:,c]     = result[c]['FM_p']
            FM_r_map[:,c]     = result[c]['FM_r']
            FM_Slope_err_map[:,c] = result[c]['FM_Slope_err']
            FM_Slope_T_map[:,c]   = result[c]['FM_Slope_T']
            FM_Slope_p_map[:,c]   = result[c]['FM_Slope_p']
        feat_vals[c,0]    = c
        feat_vals[c,1]    = result[c]['Kappa']
        feat_vals[c,2]    = result[c]['Rho']
        feat_vals[c,3]    = 100*((ICA_maps[:,c]**2).sum()/totalvar)


        FR2_mask_arr     = ma.masked_array(F_R2_maps[:,c], mask=np.logical_not(Kappa_masks[:,c])) #abs(Kappa_masks[:,c]-1))
        feat_vals[c,4]    = FR2_mask_arr.max()

        FS0_mask_arr     = ma.masked_array(F_S0_maps[:,c], mask=np.logical_not(Rho_masks[:,c])) #abs(Rho_masks[:,c]-1))
        feat_vals[c,5]    = FS0_mask_arr.max()

        feat_vals[c,6]    = feat_vals[c,1] / feat_vals[c,2]

        ZICA_mask_arr    = ma.masked_array(oc_ICA_maps[:,c], mask=np.logical_not(oc_ICA_maps_mask[:,c])) #abs(ICA_maps_mask[:,c]-1))
        feat_vals[c,7]    = ZICA_mask_arr.max()

        feat_vals[c,8]    = ICA_maps_mask[:,c].sum()
        feat_vals[c,9]    = F_R2_masks[:,c].sum()
        feat_vals[c,10]   = F_S0_masks[:,c].sum()
        feat_vals[c,11]   = Kappa_masks[:,c].sum()
        feat_vals[c,12]   = Rho_masks[:,c].sum()

        # DAN METRIC
        if doFM:
            aux_mask      = np.logical_and((FM_p_map[:,c]<0.05),(FM_Slope_p_map[:,c]<0.05))
            aux_numerator = Weight_maps[aux_mask,c].sum()
            aux_denominat = Weight_maps[:,c].sum()
            aux_metric = aux_numerator/aux_denominat
            print("[%d] -> aux_mask%s | aux_numerator%s | DF=%f" % (c,str(aux_mask.shape),str(Weight_maps[aux_mask,c].shape),aux_metric))

    df_feats = pd.DataFrame(data=feat_vals,columns=feat_names)
    df_feats.to_csv(outDir+outPrefix+'.DF.csv')
    df_feats['cID'] = df_feats['cID'].astype(int)

    niiwrite_nv(beta      , mask,outDir+outPrefix+'.chComp.Beta.nii',aff ,head)
    niiwrite_nv(F_S0_maps , mask,outDir+outPrefix+'.chComp.FS0.nii',aff ,head)
    niiwrite_nv(F_R2_maps , mask,outDir+outPrefix+'.chComp.FR2.nii',aff ,head)
    niiwrite_nv(F_S0_masks, mask,outDir+outPrefix+'.chComp.FS0.mask.nii',aff ,head)
    niiwrite_nv(F_R2_masks, mask,outDir+outPrefix+'.chComp.FR2.mask.nii',aff ,head)
    niiwrite_nv(c_S0_maps , mask,outDir+outPrefix+'.chComp.cS0.nii',aff ,head)
    niiwrite_nv(c_R2_maps , mask,outDir+outPrefix+'.chComp.cR2.nii',aff ,head)
    niiwrite_nv(p_S0_maps , mask,outDir+outPrefix+'.chComp.pS0.nii',aff ,head)
    niiwrite_nv(p_R2_maps , mask,outDir+outPrefix+'.chComp.pR2.nii',aff ,head)

    niiwrite_nv(Kappa_maps,  mask,outDir+outPrefix+'.chComp.Kappa.nii',aff ,head)
    niiwrite_nv(Kappa_masks, mask,outDir+outPrefix+'.chComp.Kappa_mask.nii',aff ,head)
    niiwrite_nv(Rho_masks,   mask,outDir+outPrefix+'.chComp.Rho_mask.nii',aff ,head)
    niiwrite_nv(Rho_maps,    mask,outDir+outPrefix+'.chComp.Rho.nii',aff ,head)
    niiwrite_nv(Weight_maps, mask,outDir+outPrefix+'.chComp.weightMaps.nii',aff ,head)
    if doFM:
        niiwrite_nv(FM_Slope_map , mask,outDir+outPrefix+'.chComp.FM.Slope.nii',aff ,head)
        niiwrite_nv(FM_Inter_map , mask,outDir+outPrefix+'.chComp.FM.Inter.nii',aff ,head)
        niiwrite_nv(FM_p_map , mask,outDir+outPrefix+'.chComp.FM.p.nii',aff ,head)
        niiwrite_nv(FM_r_map , mask,outDir+outPrefix+'.chComp.FM.r.nii',aff ,head)
        niiwrite_nv(FM_Slope_err_map , mask,outDir+outPrefix+'.chComp.FM.Slope.err.nii',aff ,head)
        niiwrite_nv(FM_Slope_T_map , mask,outDir+outPrefix+'.chComp.FM.Slope.T.nii',aff ,head)
        niiwrite_nv(FM_Slope_p_map , mask,outDir+outPrefix+'.chComp.FM.Slope.p.nii',aff ,head)
    return df_feats

def getSVDThreshold(svd_input,u,s,v,verb=False,SEprogram=False,Ne=None):
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
    verb: If true, the function will print(out some metrics about the results.)
          By default is set to False

   Results:
   --------
   Nc:  Number of non-gaussian noise components.
    """
    m = svd_input.shape[1]
    n = svd_input.shape[0]
    if SEprogram==True and not (Ne is None):
        n = n / Ne
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

def writeCompTable(origCommandLine, out_dir,data_file, features, kvar_afterPCA, kvar_afterICA, kvar_FINAL, psel, Nt, sort_col):
    Nc,_ = features.shape
    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    Ngood = (psel==1).sum()
    Nbad  = (psel==0).sum()
    midk  = []
    rej                 = np.where(psel==0)[0]
    np.savetxt(out_dir+'rejected.txt',rej,fmt='%d',delimiter=',')
    acc                 = np.where(psel==1)[0]
    np.savetxt(out_dir+'accepted.txt',acc,fmt='%d',delimiter=',')
    open(out_dir+'midk_rejected.txt','w').write(','.join([str(int(cc)) for cc in midk]))
    with open(out_dir+'comp_table.txt','w') as f:
        f.write("#Original Command line:\n")
        f.write("# %s \n" % (origCommandLine))
        f.write("#\n#ME-ICA Component statistics table for: %s \n" % (data_file))
        f.write("#Run on %s \n" % (st))
        f.write("#\n")
        f.write("#Dataset variance kept after PCA: %.02f \n" %  ( kvar_afterPCA ) )
        f.write("#Dataset variance kept after ICA: %.02f \n" %  ( kvar_afterICA ) )
        f.write("#Dataset variance explained by ICA (VEx): %.02f \n" %  ( kvar_afterICA ) ) # Repeat for compatibility with report
        f.write("#Dataset variance kept in denoised TS: %.02f \n" %  ( kvar_FINAL ) )
        f.write("#Total components generated by decomposition (TCo): %i \n" %  ( Nc ) )
        f.write("#No. accepted BOLD-like components, i.e. effective degrees of freedom for correlation (lower bound; DFe): %i\n" %  ( Ngood ) )
        f.write("#Total number of rejected components (RJn): %i\n" %  (Nbad) )
        f.write("#Nominal degress of freedom in denoised time series (..._medn.nii.gz; DFn): %i \n" %  (Nt-Nbad) )
        f.write("#ACC %s \t#Accepted BOLD-like components\n" % ','.join([str(int(cc)) for cc in acc]) )
        f.write("#REJ %s \t#Rejected non-BOLD components\n" % ','.join([str(int(cc)) for cc in rej]) )
        f.write("#MID    \t#Rejected R2*-weighted artifacts\n")
        f.write("#IGN    \t#Ignored components (kept in denoised time series)\n")
        f.write("#VEx  TCo   DFe   RJn   DFn   \n")
        f.write("##%.02f  %i %i %i %i \n" % (kvar_afterICA,Nc,Ngood,Nbad,Nt-Nbad))
        f.write("#")
        features.to_csv(f,sep='\t', index=False)
