#!/usr/bin/env python
__version__="0.10"
import sys
import os

path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'melib'))
print("++ INFO [Main]: Using meBasic library located in: %s" % path)
if not path in sys.path:
    sys.path.insert(1, path)
del path

help_desc = """
This program performs multi-echo (ME) independent component analysis (ICA) denoising
of multi-echo fMRI time series. The program takes as input fully pre-processed multi-
echo fMRI datasets, echo times used duing acquisition, a series of parameters that
affect the computation of component characteristics and the outputs.
"""

def dep_check():
    print("++ INFO [Main]: Checking for dependencies....")
    fails                = 0
    if sys.version_info > (3,0):
        pickle = "pickle"
    else:
        pickle = "cPickle"
    modules = set(["numpy", "pandas", "matplotlib","argparse","scipy","sklearn",pickle,"bz2"])
    
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

if __name__=='__main__':
    print("--------------------------------------------------")
    print("-- SFIM ME-ICA (Separate Echoes) version %s --" % __version__)
    print("--------------------------------------------------")

    dep_check() 
    import numpy             as np
    import pandas            as pd
    import matplotlib.pyplot as plt
    import meBasics          as meb
    import bz2
    from argparse              import ArgumentParser,RawTextHelpFormatter
    from sklearn.decomposition import FastICA
    from sklearn.preprocessing import StandardScaler
    from multiprocessing       import cpu_count
    from scipy.stats           import skew, zscore

    if sys.version_info > (3,0):
        import pickle
    else:
        import cPickle as pickle
    # Parse input parameters
    # ----------------------
    origCommandLine = " ".join(sys.argv[:])
    parser = ArgumentParser(description=help_desc,formatter_class=RawTextHelpFormatter)
    prsCch = parser.add_argument_group('Component Characterization')
    prsInp = parser.add_argument_group('Input Options')
    prsOut = parser.add_argument_group('Output Options')
    statFitGrp = parser.add_argument_group('Arguments for static T2* and S0 fits')
    
    prsInp.add_argument("-d","--orig_data",        dest='data_file',      help="Spatially concatenated Multi-Echo Dataset",type=str,   default=None)
    prsInp.add_argument("-e","--TEs",              dest='tes',            help="Echo times (in ms) ex: 15,39,63",          type=str,   default=None)
    prsInp.add_argument(     "--TR",               dest='TR',             help="Repetion Time (in secs). ex: 2",           type=float, default=None)
    prsInp.add_argument(     "--tes_file",         dest='tes_file',       help="Path to file with Echo time information",  type=str,   default=None)
    prsOut.add_argument(     "--out_dir",          dest='out_dir',        help="Output Directory. Default: current directory",type=str,default='.')
    parser.add_argument(     "--polort",           dest='polort',         help="Order of legendre polinomial to fit to residuals. Default=5",type=int,default=5)
    prsOut.add_argument(     "--prefix",           dest='prefix',         help="Output File Prefix. Default = sbj",type=str,default='sbj')
    parser.add_argument("-v","--verbose",          dest='writeExtra',     help="Write extra files. Default = 0. Values = 1 (write all)",action='store_true')
    parser.add_argument(     "--noQA",             dest='doQA',           help='Trigger to switch QA OFF.',action='store_false')
    parser.add_argument(     "--use_preTEfit",     dest='doPreTEdn',      help='Trigger to switch TE-dependence fit prior to PCA/ICA ON.',action='store_true')
    prsCch.add_argument(     "--use_QAweight",     dest='doQAw',          help='Trigger to switch QA-weigths in Kappa/Rho computations ON.',action='store_true')
    prsCch.add_argument("-r","--krRatio",          dest='krRatio',        help='Selection criteria. Default 1.25',type=float, default=1.25)
    parser.add_argument(     "--mask",             dest='mask_file',      help='Path to the mask to use during the analysis. If not provided, one will be computed automatically',type=str, default=None)
    parser.add_argument(     "--reuse",            dest='reuse',          help='Use this option if you want to omit recomputation of any outputs that already exist in the directory', action='store_true')
    parser.add_argument(     "--ncpus",            dest='Ncpus',          help='Number of cpus available. Default will be #available/2', default=None, type=int)
    prsOut.add_argument(     "--save_extra",       dest='save_extra',     help='Trigger to write to disk additional files.',action='store_true')
    prsCch.add_argument("-z","--ICAZThr",          dest='ica_zthr',       help='Threshold for the Z-score ICA maps', default=0, type=float)
    prsCch.add_argument(     "--Fmax",             dest='Fmax',           help='Maximum value for F-stats in TE-dependent fit. Default= No Max', default=None, type=float)
    prsCch.add_argument(     "--Zmax",             dest='Zmax',           help='Maximum value for Z-stats in ICA Maps. Default= No Max', default=None, type=float)
    prsCch.add_argument(     "--do_FullModel",     dest='doFM',           help='Trigger to request Full Model Computation ON.',action='store_true')
    prsCch.add_argument(     "--doMedian",         dest='doMedian',       help='When computing Kappa and Rho, use median instead of mean',action='store_true')
    
    statFitGrp.add_argument("--So_init", dest='So_init',  help='Initial Guess for the S0 value.  Default=2500', type=float, default=2500)
    statFitGrp.add_argument("--T2s_init",dest='T2s_init', help='Initial Guess for the T2s value. Default=40ms', type=float, default=40)
    statFitGrp.add_argument("--So_min",  dest='So_min',   help='Lowest admissible S0 value.      Default=100', type=float,  default=100)
    statFitGrp.add_argument("--So_max",  dest='So_max',   help='Highest admissible S0 value.     Default=10000', type=float,default=10000)
    statFitGrp.add_argument("--T2s_min", dest='T2s_min',  help='Lowest admissible T2s value.     Default=10ms', type=float, default=10)
    statFitGrp.add_argument("--T2s_max", dest='T2s_max',  help='Highest admissible T2s value.    Default=300ms', type=float,default=300)
    options  = parser.parse_args()
    krRatio  = float(options.krRatio)
    ica_zthr = float(options.ica_zthr)
    if options.Fmax==None:
        Fmax = np.inf
    else:
        Fmax     = float(options.Fmax)
    if options.Zmax==None:
        Zmax = np.inf
    else:
        Zmax     = float(options.Zmax)
    if (options.Ncpus is None) or (options.Ncpus > cpu_count()):
        Ncpu = int(cpu_count()/2)
    else:
        Ncpu = int(options.Ncpus) 

    So_init  = float(options.So_init)
    So_min   = float(options.So_min)
    So_max   = float(options.So_max)
    T2s_init = float(options.T2s_init)
    T2s_min  = float(options.T2s_min)
    T2s_max  = float(options.T2s_max)

    print("++ INFO [Main]: K/R Ratio = %f" % krRatio)
    print("++ INFO [Main]: ICA Z Threshold = %f" % ica_zthr)
    print("++ INFO [Main]: F-stat Max Value = %f" % Fmax)
    print("++ INFO [Main]: Z-stat Max Value = %f" % Zmax)
    print("++ INFO [Main]: Reuse exiting results? %s" % options.reuse)
    print("++ INFO [Main]: Number of CPUs to use: %d" % (Ncpu))
    print("++ INFO [Main]: Use median for Kappa/Rho computations?: %s" % (options.doMedian))
    # Control all necessary inputs are available
    # ------------------------------------------
    if options.tes is None and options.tes_file is None:
        print("++ Error: No information about echo times provided. Please do so via --TEs or --tes_file")
        sys.exit()
    if (not (options.tes is None)) and (not (options.tes_file is None)):
        print("++ Error: Echo times provided in two different ways (--TEs and --tes_file). Please select only one input mode")
        sys.exit()
    if options.data_file is None:
        print("++ Error: No input dataset provided")
        sys.exit()
    if options.TR is None:
        print("++ Error: No Repetition Time provided")
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
    outputDir = os.path.abspath(options.out_dir) + "/"
    
    
    # =========  ME-ICa starts   =========
    # ====================================
    
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
    
    # Get the repetition time
    # -----------------------
    TR = float(options.TR)
    print(" +              Repetition Time: %s" % (str(TR)))
    
    # Load ME-EPI data
    # ----------------
    print("++ INFO [Main]: Loading ME dataset....")
    mepi_data,mepi_aff,mepi_head  = meb.niiLoad(options.data_file)
    Nx,Ny,Nz,Nt                   = mepi_data.shape
    Nz                            = Nz/Ne # Because the input was the Z-concatenated dataset
    mepi_data                     = mepi_data.reshape((Nx,Ny,Nz,Ne,Nt),order='F')
    print(" +              Dataset dimensions: [Nx=%i,Ny=%i,Nz=%i,Ne=%i,Nt=%i]" % (Nx,Ny,Nz,Ne,Nt))
    
    # Generate Working mask
    # ---------------------
    origMask_needed = True
    origMask_path   = os.path.join(outputDir,options.prefix+'.mask.orig.nii')
    if options.reuse:
       if os.path.exists(origMask_path):
           print("++ INFO [Main]: Re-using existing original mask [%s]" % (origMask_path))
           mask,_,_ = meb.niiLoad(origMask_path)
           mask = (mask>0)
           origMask_needed = False
    if origMask_needed:       
       if options.mask_file==None:
          print("++ INFO [Main]: Generating initial mask from data.")
          mask           = meb.mask4MEdata(mepi_data)
          meb.niiwrite_nv(mask[mask],mask,options.out_dir+options.prefix+'.mask.orig.nii',mepi_aff ,mepi_head)
       else:
          print("++ INFO [Main]: Using user-provided mask.")
          if not os.path.exists(options.data_file):
             print("++ Error: Provided mask [%s] does not exists." % options.mask_file)
             sys.exit()
          mask,_,_       = meb.niiLoad(options.mask_file)
          mask = (mask>0)
    Nv             = np.sum(mask)
    print(" +              Number of Voxels in mask [Nv=%i]" % Nv)
    
    # Reshape the input data
    # ----------------------
    SME      = mepi_data[mask,:,:].astype(float) #(Nv,Ne,Nt)
    
    # Compute mean across time
    # ------------------------
    print("++ INFO [Main]: Computing Mean across time for all echoes....")
    # There are two options here:
    #   (1) Simply use the mean command.
    Smean_case01 = SME.mean(axis=-1)
    #   (2) Compute the mean at the same time you fit some legendre polynomials, to remove the influence
    #   of small drits in the computation. It should make almost no difference.
    Smean_case02 = meb.getMeanByPolyFit(SME,polort=7)
    SME_mean     = Smean_case02
    meb.niiwrite_nv(SME_mean,mask,options.out_dir+options.prefix+'.SME.mean.nii',mepi_aff ,mepi_head)
    
    # Perform QA (New updated version)
    # --------------------------------
    print("++ INFO [Main]: Quality Assurance....")
    QA_SSE_path      = os.path.join(outputDir,options.prefix+'.QAC.SSE.nii')
    QA_SSE_rank_path = os.path.join(outputDir,options.prefix+'.QAC.SSE_rank.nii')
    if options.doQA==1:
        if options.reuse and os.path.exists(QA_SSE_path) and os.path.exists(QA_SSE_rank_path):
           print(" +              Loading pre-existing QA_SSE [%s]" % (QA_SSE_path))
           print(" +              Loading pre-existing QA_SSE_Rank [%s]" % (QA_SSE_rank_path))
           QA_SSE,_,_      = meb.niiLoad(QA_SSE_path);QA_SSE=QA_SSE[mask]
           QA_SSE_Rank,_,_ = meb.niiLoad(QA_SSE_rank_path);QA_SSE_Rank=QA_SSE_Rank[mask]
        else: 
           print(" +              Computing QA metrics from data.")
           QA_SSE,QA_SSE_Rank=meb.computeQA(SME,tes,Ncpu,data_mean=SME_mean)
           #meb.niiwrite_nv(QA_SSE_Rank, mask, QA_SSE_rank_path, mepi_aff ,mepi_head)
           #meb.niiwrite_nv(QA_SSE,      mask, QA_SSE_path     , mepi_aff ,mepi_head)
        
    # Compute static S0 and R2* maps
    # ------------------------------
    print("++ INFO [Main]: Compute Static S0 and R2* maps...")
    stFit_S0_path  = os.path.join(outputDir,options.prefix+'.sTE.S0.nii')
    stFit_t2s_path = os.path.join(outputDir,options.prefix+'.sTE.t2s.nii')
    stFit_SSE_path = os.path.join(outputDir,options.prefix+'.sTE.SSE.nii')
    stFit_bVx_path = os.path.join(outputDir,options.prefix+'.sTE.mask.bad.nii')
    mask_bad_staticFit = np.zeros((Nv,), dtype=bool)
    if options.reuse and os.path.exists(stFit_S0_path) and os.path.exists(stFit_t2s_path) and os.path.exists(stFit_SSE_path) and os.path.exists(stFit_bVx_path):
       print(" +              Loading pre-existing static S0 [%s] map." % (stFit_S0_path))
       print(" +              Loading pre-existing static t2s [%s] map." % (stFit_t2s_path))
       print(" +              Loading pre-existing bad static fit voxels [%s] map." % (stFit_bVx_path))
       S0,_,_  = meb.niiLoad(stFit_S0_path); S0 = S0[mask] 
       t2s,_,_ = meb.niiLoad(stFit_t2s_path); t2s     = t2s[mask]
       mask_bad_staticFit,_,_ = meb.niiLoad(stFit_bVx_path); mask_bad_staticFit     = mask_bad_staticFit[mask]
    else:
       print(" +              Computing static S0 and t2s maps from the data.")
       # Do a non-linear optimization to fit the original curve (no log-linear transformation)
       # using an optimization algorithm that takes boundaries.
       S0, t2s, SSE, mask_bad_staticFit = meb.make_static_maps_opt(SME_mean,tes,Ncpu,So_init=So_init,T2s_init=T2s_init,So_min=So_min,So_max=So_max,T2s_min=T2s_min, T2s_max=T2s_max)
       mask_bad_staticFit = np.logical_not(mask_bad_staticFit)
       meb.niiwrite_nv(S0                ,mask,stFit_S0_path, mepi_aff ,mepi_head)
       meb.niiwrite_nv(t2s               ,mask,stFit_t2s_path,mepi_aff ,mepi_head)
       meb.niiwrite_nv(SSE               ,mask,stFit_SSE_path,mepi_aff ,mepi_head)
       meb.niiwrite_nv(mask_bad_staticFit,mask,stFit_bVx_path,mepi_aff ,mepi_head) # In mask_bad_staticFit, voxels that should be discarded are = 1 both in Python and in the associated AFNI brik
       # There is a simpler, yet less accurate way to do this. 
       # Do a log-linear fit. The lack of boundaries and the log-linear transformation leads to
       # some voxels having completely wrong values for TE and S0. The benefit of this is that 
       # this method is really fast.
       # S0,t2s,_,_ = meb.make_static_maps(S,tes)
       
    # Compute Linear Fit to Remove non-S0/T2* signals
    # -----------------------------------------------
    print("++ INFO [Main]: Remove non-S0/T2* signals....")
    dTE_TSb4Fit_path = os.path.join(outputDir,options.prefix+'.dTE.TSb4Fit.nii')
    dTE_TSa4Fit_path = os.path.join(outputDir,options.prefix+'.dTE.TSa4Fit.nii')
    dTE_FitRes_path  = os.path.join(outputDir,options.prefix+'.dTE.FitRes.nii')
    if (options.doPreTEdn==1):
       if options.reuse and os.path.exists(dTE_TSb4Fit_path) and os.path.exists(dTE_TSa4Fit_path) and os.path.exists(dTE_FitRes_path):
          print(" +              Load pre-existing ME data without non-S0/T2* signals.")
          datahat,_,_ = meb.niiLoad(dTE_TSa4Fit_path)
          datahat     = np.reshape(datahat,(Nx,Ny,Nz,Ne,Nt),order='F')
          datahat     = datahat[mask,:,:]
       else: 
          print(" +              Computing dynamic S0/T2* fits.")
          dkappa,drho,fitres,rcond,datahat,_ = meb.linearFit(SME,tes,Ncpu)  
          meb.niiwrite_nv(SME,    mask,dTE_TSb4Fit_path, mepi_aff ,mepi_head)
          meb.niiwrite_nv(datahat,mask,dTE_TSa4Fit_path, mepi_aff ,mepi_head)
          meb.niiwrite_nv(fitres ,mask,dTE_FitRes_path,  mepi_aff ,mepi_head)
       SME = datahat
       SME_mean = meb.getMeanByPolyFit(SME,polort=7)
    else:
       print(" +              OFF")
 
    # Create all necessary versions of the signal
    # -------------------------------------------
    SME_nm  = SME - SME_mean[:,:,np.newaxis]     #(Nv,Ne,Nt) SME_mn = Multi-echo data with no mean (deltaS)
    SME_pc  = SME_nm / SME_mean[:,:,np.newaxis]  #(Nv,Ne,Nt) SME_pc   = Multi-echo data with no mean / mean (deltaS/Smean)

    # Compute Optimally Combined Time-series
    # --------------------------------------
    print("++ INFO [Main]: Compute Optimally Combined Time series...")
    octs_path = os.path.join(outputDir,'ts_OC.nii')
    if options.reuse and os.path.exists(octs_path):
        print(" +              Load pre-existing optimally combined time series [%s]" % octs_path)
        octs,_,_ = meb.niiLoad(octs_path)
        octs = octs[mask,:]
    else:
        print(" +              Compute Optimally Combined Time series from data.")
        octs = meb.make_optcom(SME,t2s,tes)                                                             #(Nv,Nt)
        meb.niiwrite_nv(octs,mask,octs_path,mepi_aff ,mepi_head)

    # SVD/PCA Decomposition
    # ---------------------
    print("++ INFO [Main]: Removal of Gaussian Noise...")
    pca_octs_path   = os.path.join(outputDir,options.prefix+'.PCA.sets.nii')
    pca_out_path    = os.path.join(outputDir,options.prefix+'.PCA.outts.nii')
    pca_pickle_path = os.path.join(outputDir,options.prefix+'.PCA.pklbz')
    
    # Remove the mean from the OCTS prior to PCA
    pca_input = SME_nm.copy()
    pca_input = np.reshape(pca_input,(Nv*Ne,Nt),order='F')
    
    # THIS IS HERE FOR REFERNECE. DO NOT USE.
    #input_scaler  = StandardScaler(copy=True, with_mean=True, with_std=False)
    #pca_input   = input_scaler.fit_transform(np.reshape(SME,(Nv*Ne,Nt),order='F').T).T        # Z-score in the time dimension (Nv,Nt)
    #meb.niiwrite_nv(np.reshape(pca_input,(Nv,Ne,Nt),order='F'),    mask,os.path.join(outputDir,options.prefix+'.PCA.input.nii'), mepi_aff ,mepi_head)
    #pca_input = np.reshape(pca_input,(Nv*Ne,Nt),order='F')
    
    if options.reuse and os.path.exists(pca_octs_path) and os.path.exists(pca_out_path) and os.path.exists(pca_pickle_path):
        print(" +              Loading previous results from PCA computations: OCTS   =%s" % (pca_octs_path))
        print(" +              Loading previous results from PCA computations: PCA_OUT=%s" % (pca_out_path))
        print(" +              Loading previous results from PCA computations: OTHERS =%s" % (pca_pickle_path))
        
        pcastate_f = bz2.BZ2File(pca_pickle_path,'rb')
        pcastate = pickle.load(pcastate_f)
        for key,val in pcastate.items(): exec(key + '=val')
        print(" +              Number of PCA Components = %d" % Nc)
        pca_out,_,_       = meb.niiLoad(pca_out_path)
        pca_out           = pca_out.reshape((Nx,Ny,Nz,Ne,Nt),order='F')
        pca_out           = pca_out[mask,:,:].astype(float)
        pca_out           = np.reshape(pca_out,(Nv*Ne,Nt),order='F')
        se_afterPCA,_,_       = meb.niiLoad(pca_octs_path)
        se_afterPCA           = se_afterPCA.reshape((Nx,Ny,Nz,Ne,Nt),order='F')
        se_afterPCA           = se_afterPCA[mask,:,:].astype(float)
        se_afterPCA           = np.reshape(se_afterPCA,(Nv*Ne,Nt),order='F')
    else:
        print(" +              Performing PCA on optimally combined timeseries.")
        u,s,v       = np.linalg.svd(pca_input, full_matrices=0) # u(Nv,Nt),s(Nt,),v(Nc_pca,Nt)
        pca_mmix     = v #(Npca,Nt) --> N_pca = Nt as of now
        print(" +              Computing number of PCA components in dataset.")
        Nc = meb.getSVDThreshold(pca_input,u,s,v, verb=True) # Add paper reference here.
        print(" +              Writing PCA variance plot.")
        f,ax = plt.subplots(1,1)
        pca_var = 100*s**2/((s**2).sum())
        ax.plot(pca_var,'k')
        ax2 = ax.twinx()
        ax2.plot(np.cumsum(pca_var),'b--')
        ax2.plot([Nc,Nc],[0,100],'r--')
        ax.set_xlabel('PCA Components')
        ax.set_ylabel('Variance Explained')
        ax2.set_ylabel('Cumulative Variance Explained',color='b')    
        f.savefig(options.out_dir+options.prefix+'.PCA_VariancePlot.pdf')

        print(" +              Writing OCTS after PCA.")
        pca_out       = np.dot(np.dot(v[0:Nc,:].T,v[0:Nc,:]),pca_input.T).T    #(Nv*Ne,Nt)
        se_afterPCA   = np.reshape(pca_out.copy(),(Nv,Ne,Nt),order='F')
        se_afterPCA   = se_afterPCA + SME_mean[:,:,np.newaxis] 
        meb.niiwrite_nv(np.reshape(se_afterPCA,(Nv,Ne,Nt),order='F'),mask,pca_octs_path,mepi_aff ,mepi_head)
        meb.niiwrite_nv(np.reshape(pca_out,(Nv,Ne,Nt),order='F'),    mask,pca_out_path, mepi_aff ,mepi_head)
    
        print(" +              Serializing, compressing and writting to disk of PCA output objects.")
        pcastate   = {'Nc':Nc}
        pcastate_f = bz2.BZ2File(pca_pickle_path,'wb')
        pickle.dump(pcastate,pcastate_f)
        pcastate_f.close()
    
    origME_varexp = (np.reshape(SME,(Nv*Ne,Nt),order='F').std(axis=-1)**2).sum()
    ppcaME_varexp = (np.reshape(se_afterPCA,(Nv*Ne,Nt),order='F').std(axis=-1)**2).sum()
    kvar_afterPCA = 100*ppcaME_varexp/origME_varexp
    print(" +              Var Explained after PCA: %0.2f" % (kvar_afterPCA))
    
    # ICA Decomposition
    # -----------------
    print("++ INFO [Main]: ICA Decomposition...")
    ica_mix_path       = os.path.join(outputDir,'meica_mix.1D')
    ica_mix_zsc_path   = os.path.join(outputDir,'meica_mix_zsc.1D')
    ica_mix_fft_path   = os.path.join(outputDir,'meica_mix_fft.1D')
    ica_mix_freqs_path = os.path.join(outputDir,'meica_mix_fft_freqs.1D')
    ica_maps_path      = os.path.join(outputDir,options.prefix+'.ICA.Zmaps.nii')
    ica_octs_path      = os.path.join(outputDir,options.prefix+'.ICA.sets.nii')
    fica_input         = pca_out.copy()              #(Nv*Ne,Nt)
    if options.reuse and os.path.exists(ica_mix_path) and os.path.exists(ica_mix_zsc_path) and os.path.exists(ica_mix_fft_path) and os.path.exists(ica_mix_freqs_path) and os.path.exists(ica_maps_path) and os.path.exists(ica_octs_path):
        print(" +              Loading pre-existing ICA Mixing matrix [%s]." % (ica_mix_path))
        print(" +              Loading pre-existing Normalized ICA Mixing matrix [%s]." % (ica_mix_zsc_path))
        print(" +              Loading pre-existing ICA Maps [%s]." % (ica_maps_path))
        fica_mmix      = np.loadtxt(ica_mix_path).T
        fica_mmix_zsc  = np.loadtxt(ica_mix_zsc_path).T
        fica_out,_,_   = meb.niiLoad(ica_maps_path)
        fica_out       = fica_out.reshape((Nx,Ny,Nz,Ne,Nc),order='F')
        fica_out       = fica_out[mask,:,:].astype(float)
        fica_out       = np.reshape(fica_out,(Nv*Ne,Nc),order='F')
        se_afterICA,_,_  = meb.niiLoad(ica_octs_path)
        se_afterICA      = se_afterICA.reshape((Nx,Ny,Nz,Ne,Nt),order='F')
        se_afterICA      = se_afterICA[mask,:,:].astype(float)
        se_afterICA      = np.reshape(se_afterICA,(Nv*Ne,Nt),order='F')
    else:
        print(" +              Perform ICA....")
        fica       = FastICA(n_components=Nc, max_iter=1000)
        fica_out   = fica.fit_transform(fica_input).T #(Nc,Nv*Ne)
        fica_out -= fica_out.mean(axis=0)
        fica_out /= fica_out.std(axis=0)
        fica_out  = fica_out.T        #(Nv*Ne,Nc)
        fica_mmix  = fica.mixing_.T   #(Nc,Nt)
        # Correct the sign of components
        # ------------------------------ 
        print(" +              Correct the sign of the ICA components....")
        fica_signs          = skew(fica_out,axis=0) #(Nc,)
        fica_signs         /= np.abs(fica_signs)                        #(Nc,)
        fica_out            = (fica_out.T*fica_signs[:,np.newaxis]).T   #(Nv*Ne,Nc)
        fica_mmix           = (fica_mmix*fica_signs[:,np.newaxis])      #(Nc,Nt)
        fica_mmix_zsc       = zscore(fica_mmix,axis=-1)                 #(Nc,Nt)
        # Save ICA Mixing matrix, its normalized version and the ICA maps
        print(" +              Saving ICA maps as %s" % ica_maps_path)
        print(" +              Saving ICA representative timeseries as %s" % ica_mix_path)
        print(" +              Saving Z-scores ICA representative timeseries as %s" % ica_mix_zsc_path)
        np.savetxt(ica_mix_path,fica_mmix.T)
        np.savetxt(ica_mix_zsc_path,fica_mmix_zsc.T)
        meb.niiwrite_nv(np.reshape(fica_out,(Nv,Ne,Nc),order='F'),mask,ica_maps_path,mepi_aff ,mepi_head)
        # Compute the FFT of the ICA representative time series for report  
        # ----------------------------------------------------------------
        print(" +              Saving the Fourier Transform of the ICA timeseries as %s" % ica_mix_fft_path)
        FFT, freq = meb.computeFFT(fica_mmix,TR)
        np.savetxt(ica_mix_fft_path,FFT.T)
        np.savetxt(ica_mix_freqs_path,freq)
        # Compute a few additional outputs for the ICA
        # --------------------------------------------
        print(" +              Writing OCTS after ICA [%s]" % ica_octs_path)
        se_beta     = np.linalg.lstsq(fica_mmix_zsc.T,np.reshape(SME_nm,(Nv*Ne,Nt),order='F').T)[0].T  #(Nv*Ne,Nt)
        se_afterICA = np.dot(se_beta,fica_mmix_zsc)                  #(Nv*Ne,Nt)
        se_afterICA = np.reshape(se_afterICA,(Nv,Ne,Nt),order='F')   #(Nv,Ne,Nt)
        se_afterICA = se_afterICA + SME_mean[:,:,np.newaxis]         #(Nv,Ne,Nt)
        meb.niiwrite_nv(se_afterICA,mask,ica_octs_path,mepi_aff ,mepi_head)

    picaME_varexp = (np.reshape(se_afterICA,(Nv*Ne,Nt),order='F').std(axis=-1)**2).sum()
    kvar_afterICA = 100*picaME_varexp/origME_varexp
    print(" +              Var Explained after ICA: %0.2f" % (kvar_afterICA))
    
    # Characterize ICA Components
    # ---------------------------
    print("++ INFO [Main]: Computing Kappa, Rho and Variance Explained...")
    if (options.doQAw==0): 
       voxelwiseQA=np.ones((Nv,))
       print("++ INFO [Main]: QA weigths not used during kappa/rho computation")
    else: 
       voxelwiseQA=QA_SSE_Rank
       print("++ INFO [Main]: QA weigths were used during kappa/rho computation")
    fica_feats = meb.characterize_components_se(SME_pc, SME_mean, tes, t2s, S0, fica_mmix_zsc, fica_out, voxelwiseQA, Ncpu, ICA_maps_thr=ica_zthr, F_MAX=Fmax, Z_MAX=Zmax,
                 outDir=options.out_dir,
                 outPrefix=options.prefix, 
                 mask=mask,writeOuts=options.save_extra,
                 aff=mepi_aff, head=mepi_head, discard_mask=mask_bad_staticFit,
                 doFM=options.doFM, doMedian=options.doMedian)
    pd.options.display.float_format = '{:,.2f}'.format
    print(fica_feats)
    # Selection good and bad components
    # ---------------------------------
    print("++ INFO [Main]: Component Selection...")
    
    # (1) Selection based on kappa/rho ratio
    fica_psel      = (fica_feats['K/R']>krRatio).get_values()
    accepted       = fica_feats['cID'][fica_psel].get_values().astype(int)
    print(" +              Unsorted List of accepted components=%s" % str(accepted))
    
    # Generating  output time-series
    # ------------------------------
    print("++ INFO [Main]: Generating output time series...")
    list_goodComp = np.where(fica_psel==1)[0].tolist()
    list_badComp  = np.where(fica_psel==0)[0].tolist()
    octs_nomean   = octs - octs.mean(axis=-1)[:,np.newaxis]
    beta          = np.linalg.lstsq(fica_mmix.T, octs_nomean.T)[0].T # MAY WANT TO CHANGE THIS TO MMIX_ZSC, BUT THEN CHANGE THINGS BELOW TOO 
    goodTS        = np.dot(beta[:,list_goodComp], fica_mmix[list_goodComp,:])
    badTS         = np.dot(beta[:,list_badComp], fica_mmix[list_badComp,:])
    denoisedTS    = octs - badTS
    goodTS        = goodTS + octs.mean(axis=-1)[:,np.newaxis]
    badTS         = badTS  + octs.mean(axis=-1)[:,np.newaxis]
    meb.niiwrite_nv(goodTS,mask,options.out_dir+options.prefix+'.OUT.gTS.nii',mepi_aff,mepi_head)  
    meb.niiwrite_nv(badTS, mask,options.out_dir+options.prefix+'.OUT.bTS.nii',mepi_aff,mepi_head)
    meb.niiwrite_nv(beta ,mask ,options.out_dir+'betas_OC.nii',mepi_aff ,mepi_head)        
    meb.niiwrite_nv(denoisedTS       ,mask ,options.out_dir+'dn_ts_OC.nii',mepi_aff ,mepi_head)
    
    octsORIG_varexp  = (octs.std(axis=-1)**2).sum()
    octsFINAL_varexp = (denoisedTS.std(axis=-1)**2).sum()
    kvar_FINAL = 100*octsFINAL_varexp/octsORIG_varexp
    print(" + INFO [Main]: Var Kept in the data: %0.2f" % (kvar_FINAL))
    
    ## Output some informative metrics
    ## -------------------------------
    meb.writeCompTable(origCommandLine, options.out_dir, options.data_file, fica_feats, kvar_afterPCA, kvar_afterICA, kvar_FINAL, fica_psel, Nt, 6)
    print("++ INFO [Main]: Successfull Completion of the Analysis.")
    print("++ INFO [Main]: =======================================")
