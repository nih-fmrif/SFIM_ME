#!/usr/bin/env bash

# This script is hard coded to analyze two data sets stored on felix.nimh.nih.gov at
# /data/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SFIM_ME_TestingData
# Both data sets are from SBJ01 of the multiecho 100 runs study. The first has 3 echos
# and the second has 5 echos. This script has three required options that are just listed
# after the command call:
# OutputPrefix: This is the label to identify the output. If you use a new name here, 
#  it will save the output into a new directory. If you use an existing name, it will
#  overwrite what's already there. This option allows one to compare the results across
#  versions of the code
# NumCPUS: The number of parallel threads that can be used
# Python27Env Python3Env: This program runs data with both Python 2.7 and Python 3.5.
#   For each version "source activate $Python??Env" will be called and then the script will be run
# RootPath: If this is being run on felix/biowulf, this should be:
#    '/data/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SFIM_ME_TestingData' If this us run on a 
#    mounted drive on a Mac, you can use '//Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SFIM_ME_TestingData'

#
# An example call is: 
#  ./ProcessingADataSample.sh Original2015_11_25 /Users/handwerkerd/Documents/Code/SFIM_ME 12 3 p27werker p35werker  /data/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SFIM_ME_TestingData 12
# Running line-by-line D Handwerker's options: OutputPrefix=Original2015_11_25; SFIMMEloc=/Users/handwerkerd/Documents/Code/SFIM_ME; NumCPUS=4; NumIter=3; Python2Env=p27; Python3Env=p3; RootPath='//Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SFIM_ME_TestingData'
#
# D.A. Handwerker 11/25/2015

OutputPrefix=$1
SFIMMEloc=$2
NumCPUS=$3
NumIter=$4
Python2Env=$5
Python3Env=$6
RootPath=$7



###################################################################
# Setting up the original directory with the volume files to test
# This only needs to be done once so I'm keeping the code here, but setting
# a variable so that it never actually runs.
OrganizeInitialFiles=0
if [ ${OrganizeInitialFiles} -gt 0 ]; then
  mkdir /data/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SFIM_ME_TestingData
  mkdir /data/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SFIM_ME_TestingData/Echos3
  mkdir /data/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SFIM_ME_TestingData/Echos5
  mkdir /data/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SFIM_ME_TestingData/Outputs
  mkdir /data/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SFIM_ME_TestingData/Outputs/Echos3
  mkdir /data/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SFIM_ME_TestingData/Outputs/Echos5
  cd mkdir /data/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SFIM_ME_TestingData
  3dcalc -a ../CrossRunAnalyses.AnatAlign/SBJ01/RegistrationComparisons/SBJ01_Anatomy+orig \
    -prefix SBJ01_Anatomy.nii.gz -expr 'a'


  cd /data/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SFIM_ME_TestingData/Echos3
  ln -s ../../SBJ01_S01/D01_Version02.AlignByAnat.Cubic/Task01/p06.SBJ01_S01_Task01_e1.sm.nii.gz ./BlockDesign_Echos3_e1.nii.gz
  ln -s ../../SBJ01_S01/D01_Version02.AlignByAnat.Cubic/Task01/p06.SBJ01_S01_Task01_e2.sm.nii.gz ./BlockDesign_Echos3_e2.nii.gz
  ln -s ../../SBJ01_S01/D01_Version02.AlignByAnat.Cubic/Task01/p06.SBJ01_S01_Task01_e3.sm.nii.gz ./BlockDesign_Echos3_e3.nii.gz
  ln -s ../../SBJ01_S01/D01_Version02.AlignByAnat.Cubic/Task01/zcat_ffd_SBJ01_S01_Task01.nii.gz ./zcat_BlockDesign_Echos3.nii.gz
  ln -s ../SBJ01_Anatomy.nii.gz ./
  # At some step in the processing pipeline, these NIFTI images get labelled TLRC, but their actually ORIG
  3drefit -view 'orig' -space ORIG p06.SBJ01_S01_Task01_e1.sm.nii.gz
  3drefit -view 'orig' -space ORIG p06.SBJ01_S01_Task01_e2.sm.nii.gz
  3drefit -view 'orig' -space ORIG p06.SBJ01_S01_Task01_e3.sm.nii.gz
  3drefit -view 'orig' -space ORIG zcat_ffd_SBJ01_S01_Task01.nii.gz 

  cd /data/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SFIM_ME_TestingData/Echos5
  ln -s ../../SBJ01_S09/D01_Version02.AlignByAnat.Cubic/Task11/p06.SBJ01_S09_Task11_e1.sm.nii.gz ./BlockDesign_Echos5_e1.nii.gz
  ln -s ../../SBJ01_S09/D01_Version02.AlignByAnat.Cubic/Task11/p06.SBJ01_S09_Task11_e2.sm.nii.gz ./BlockDesign_Echos5_e2.nii.gz
  ln -s ../../SBJ01_S09/D01_Version02.AlignByAnat.Cubic/Task11/p06.SBJ01_S09_Task11_e3.sm.nii.gz ./BlockDesign_Echos5_e3.nii.gz
  ln -s ../../SBJ01_S09/D01_Version02.AlignByAnat.Cubic/Task11/p06.SBJ01_S09_Task11_e4.sm.nii.gz ./BlockDesign_Echos5_e4.nii.gz
  ln -s ../../SBJ01_S09/D01_Version02.AlignByAnat.Cubic/Task11/p06.SBJ01_S09_Task11_e5.sm.nii.gz ./BlockDesign_Echos5_e5.nii.gz
  ln -s ../../SBJ01_S09/D01_Version02.AlignByAnat.Cubic/Task11/zcat_ffd_SBJ01_S09_Task11.nii.gz ./zcat_BlockDesign_Echos5.nii.gz
  ln -s ../SBJ01_Anatomy.nii.gz ./
  # At some step in the processing pipeline, these NIFTI images get labelled TLRC, but their actually ORIG
  3drefit -view 'orig' -space ORIG p06.SBJ01_S09_Task11_e1.sm.nii.gz
  3drefit -view 'orig' -space ORIG p06.SBJ01_S09_Task11_e2.sm.nii.gz
  3drefit -view 'orig' -space ORIG p06.SBJ01_S09_Task11_e3.sm.nii.gz
  3drefit -view 'orig' -space ORIG p06.SBJ01_S09_Task11_e4.sm.nii.gz
  3drefit -view 'orig' -space ORIG p06.SBJ01_S09_Task11_e5.sm.nii.gz
fi



export OMP_NUM_THREADS=${NumCPUS}

for echotype in Echos3 Echos5; do
  cd $RootPath/Outputs	
  mkdir ${echotype}/${OutputPrefix}
  mkdir ${echotype}/${OutputPrefix}/Python2
  mkdir ${echotype}/${OutputPrefix}/Python3
  for pythonversion in Python2 Python3; do
    cd $RootPath/Outputs/${echotype}/${OutputPrefix}/${pythonversion}
    for meicatype in meica_OC meica_SE; do
        mkdir $meicatype
        for ((iter=1;iter<=NumIter;++iter)); do
          mkdir ${meicatype}/iter${iter}
          mkdir ${meicatype}/iter${iter}
        done
    done
  done  
done


for ((iter=1;iter<=NumIter;++iter)); do
  for pythonversion in 2 3; do
    if [ $pythonversion -eq 2 ]; then
      source activate ${Python2Env}
    else
      source activate ${Python3Env}
    fi
    for echotype in Echos3 Echos5; do
      cd $RootPath/Outputs/${echotype}/${OutputPrefix}/Python${pythonversion}/meica_OC/iter${iter}
      pwd
      # WORKING HERE. THERE IS A BUG IN THE MAIN CALL:
	  ++ INFO [Main]: Computing Kappa, Rho and Variance Explained...
	  ++ INFO [Main]: QA weigths were used during kappa/rho computation
	   +              Dataset ./Original2015_11_25.ICA.Zmaps.mask.nii written to disk
	  Hello -------------->
	  beta.(24610, 5, 45)
	  ICA_maps.(24610, 45)
	  ICA_maps_mask.(24610, 45)
	  X1.(5, 24610)
	  X2.(5, 24610)
	   +              Multi-process Characterize Components -> Ncpu = 4
	   +              Dataset ./Original2015_11_25.chComp.EXTRA.Beta000.nii written to disk
	   +              Dataset ./Original2015_11_25.chComp.EXTRA.Beta003.nii written to disk
	   +              Dataset ./Original2015_11_25.chComp.EXTRA.Beta006.nii written to disk
	   +              Dataset ./Original2015_11_25.chComp.EXTRA.Beta009.nii written to disk
	   +              Dataset ./Original2015_11_25.chComp.EXTRA.S0Fit000.nii written to disk
	   +              Dataset ./Original2015_11_25.chComp.EXTRA.S0Fit003.nii written to disk
	   +              Dataset ./Original2015_11_25.chComp.EXTRA.S0Fit006.nii written to disk
	   +              Dataset ./Original2015_11_25.chComp.EXTRA.S0Fit009.nii written to disk
	   +              Dataset ./Original2015_11_25.chComp.EXTRA.R2Fit000.nii written to disk
	  Traceback (most recent call last):
	    File "/Users/handwerkerd/Documents/Code/SFIM_ME/sfim_meica_OC.py", line 423, in <module>
	      doFM=options.doFM)
	    File "/Users/handwerkerd/Documents/Code/SFIM_ME/melib/meBasics.py", line 530, in characterize_components
	      } for c in np.arange(Nc)])
	    File "/Users/handwerkerd/anaconda/envs/p27/lib/python2.7/multiprocessing/pool.py", line 251, in map
	       +              Dataset ./Original2015_11_25.chComp.EXTRA.R2Fit003.nii written to disk
	  return self.map_async(func, iterable, chunksize).get()
	    File "/Users/handwerkerd/anaconda/envs/p27/lib/python2.7/multiprocessing/pool.py", line 567, in get
	      raise self._value
	  NameError: global name 'doMedian' is not defined
	   +              Dataset ./Original2015_11_25.chComp.EXTRA.R2Fit006.nii written to disk
	   
	   
	   
	   # ALSO NEED TO ADD THE MOTION .1D file to the orig data directories
	   # NEED TO MAKE SURE THE REPORTS RUN AND ADD THE META REPORT
	   # NEED TO MAKE SURE THE meica_SE HAS THE CORRECT COMMAND CALL AND ALSO RUNS
	  
	  
	  python $SFIMMEloc/sfim_meica_OC.py --tes_file ${RootPath}/${echotype}/TES_${echotype}.1D \
		  -d ${RootPath}/${echotype}/zcat_BlockDesign_${echotype}.nii.gz --prefix ${OutputPrefix} \
			  --krRatio=2 --TR=2.0 --out_dir=./ --use_QAweight --ncpus ${NumCPUS} --save_extra --Fmax 500 --Zmax 8 -z 0.25
	  python $SFIMMEloc/report/meica_report.py -TED ./ -motion BlockDesign_${echotype}_motion.1D --ncpus ${NumCPUS} --overwrite
	  python $SFIMMEloc/report_bokeh/meica_report.py -meicaDir ./
	  
      cd $RootPath/Outputs/${echotype}/${OutputPrefix}/Python${pythonversion}/meica_SE/iter${iter}
      pwd
	  python $SFIMMEloc/sfim_meica_SE.py --tes_file ${RootPath}/${echotype}/TES_${echotype}.1D \
		  -d ${RootPath}/${echotype}/zcat_BlockDesign_${echotype}.nii.gz --prefix ${OutputPrefix} \
			  --krRatio=2 --TR=2.0 --out_dir=./ --use_QAweight --ncpus ${NumCPUS} --save_extra --Fmax 500 --Zmax 8 -z 0.25
	  python $SFIMMEloc/report/meica_report.py -TED ./ -motion BlockDesign_${echotype}_motion.1D --ncpus ${NumCPUS} --overwrite
	  python $SFIMMEloc/report_bokeh/meica_report.py -meicaDir ./	  
	  

    done
  done
done


    mkdir ${OutputPrefix}
    mkdir
    cd $OutputPrefix
    
  done
done


cd /data/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/DXX_TestingJavierMeica1

ln -s ../D01_Version02.AlignByAnat.Cubic/Task01/zcat_ffd_SBJ01_S01_Task01.nii.gz ./
ln -s ../D01_Version02.AlignByAnat.Cubic/Task01/SBJ01_Anatomy.nii.gz ./
ln -s ../D01_Version02.AlignByAnat.Cubic/Task01/SBJ01_S01_Task01_e2_Motion.1D ./

On felix
 export OMP_NUM_THREADS=16
python /data/handwerkerd/SFIM_ME/jmeica_optcom9.py --tes_file Echos_100runs.1D -d zcat_ffd_SBJ01_S01_Task01.nii.gz --prefix SBJ01_Sess01_Task01 --krRatio=5 --TR=2.0 --out_dir=./ --use_QAweight --ncpus 16 --save_extra --Fmax 500 --Zmax 8 -z 0.8 --reuse


python /data/handwerkerd/SFIM_ME/report/meica_report.py -TED ./ -motion SBJ01_S01_Task01_e2_Motion.1D -dir ./ -ax -sag -title SBJ01_Sess01_Task01 -setname ./ -label SBJ01_Sess01_Task01 -sort_col 7 -overwrite

python /data/handwerkerd/SFIM_ME/report_bokeh/meica_report.py -meicaDir ./ -runID SBJ01_Sess01_Task01 -picDir ./meica.SBJ01_Sess01_Task01/Report_Figures/

On Mac
 export OMP_NUM_THREADS=4

cd //Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/D01_Version02.AlignByAnat.Cubic/Task01
3dZcat -overwrite -prefix zcat_noaddmean.nii p04.SBJ01_S01_Task01_e1.align_clp+orig p04.SBJ01_S01_Task01_e2.align_clp+orig p04.SBJ01_S01_Task01_e3.align_clp+orig
3dcalc -overwrite -a zcat_mask_SBJ01_S01_Task01.nii.gz -b zcat_noaddmean.nii -prefix zcat_noaddmean.nii -expr 'b*ispositive(a)'
mv zcat_noaddmean.nii ../../DXX_TestingJavierMeica1/

# Using optimally combined as input to ICA
python /Users/handwerkerd/Documents/Code/SFIM_ME/sfim_meica_OC.py --tes_file Echos_100runs.1D -d zcat_noaddmean.nii --prefix SBJ01_Sess01_Task01 --krRatio=2 --TR=2.0 --out_dir=./ --use_QAweight --ncpus 4 --save_extra --Fmax 500 --Zmax 8 -z 0.25


python /Users/handwerkerd/Documents/Code/SFIM_ME/report/meica_report.py -TED ./ -motion SBJ01_S01_Task01_e2_Motion.1D -dir ./ -ax -sag -title SBJ01_Sess01_Task01 -setname ./ -label SBJ01_Sess01_Task01 -sort_col 7 -overwrite

# using concatinated echos as input to ICA
python /Users/handwerkerd/Documents/Code/SFIM_ME/sfim_meica_SE.py --tes_file Echos_100runs.1D -d zcat_noaddmean.nii --prefix SE_SBJ01_Sess01_Task01 --krRatio=2 --TR=2.0 --out_dir=./ --use_QAweight --ncpus 4 --save_extra --Fmax 500 --Zmax 8 -z 0.25

python /Users/handwerkerd/Documents/Code/SFIM_ME/report/meica_report.py -TED ./ -motion SBJ01_S01_Task01_e2_Motion.1D -dir ./ -ax -sag -title SE_SBJ01_Sess01_Task01 -setname ./ -label SE_SBJ01_Sess01_Task01 -sort_col 7 -overwrite






 -anat SBJ01_Anatomy.nii.gz

python /Users/handwerkerd/Documents/Code/SFIM_ME/report_bokeh/meica_report.py -meicaDir ./ -runID SBJ01_Sess01_Task01 -picDir ./meica.SBJ01_Sess01_Task01/Report_Figures/





# Checking smoothness estimations with the two preprocessing streams
cd //Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/Data4Dan/SBJ01FCH/D02_Preproc
3dFWHMx -automask -detrend -input pc03.SBJ01FCH_S01Run01_E01.volreg.AX.nii.gz
 3.08354  2.84319  2.58934

cd //Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/D01_Version02.AlignByAnat.Cubic/Task01
3dFWHMx -automask -detrend -input p06.SBJ01_S01_Task01_e1.sm.nii.gz
 3.61379  3.85934  3.70338

# normalizing first results in noticably more smoothing


# My preprocessing:
python /Users/handwerkerd/Documents/Code/SFIM_ME/jmeica_optcom9.py --tes_file Echos_100runs.1D -d zcat_ffd_SBJ01_S01_Task01.nii.gz --prefix SBJ01_Sess01_Task01 --krRatio=2 --TR=2.0 --out_dir=./ --use_QAweight --ncpus 4 --save_extra --Fmax 500 --Zmax 8 -z 0.25
++ INFO [Main]: Using meBasic library located in: /Users/handwerkerd/Documents/Code/SFIM_ME/melib
--------------------------------------------------
-- SFIM ME-ICA (Optimally Combined) version 0.10 --
--------------------------------------------------
++ INFO [Main]: Checking for dependencies....
 +              All Dependencies are OK.
++ INFO [Main]: K/R Ratio = 2.000000
++ INFO [Main]: ICA Z Threshold = 0.250000
++ INFO [Main]: F-stat Max Value = 500.000000
++ INFO [Main]: Z-stat Max Value = 8.000000
++ INFO [Main]: Reuse exiting results? False
++ INFO [Main]: Number of CPUs to use: 4
++ INFO [Main]: Reading echo times from input echo time file.
 +              Echo times: [ 15.4  29.7  44. ]
 +              Repetition Time: 2.0
++ INFO [Main]: Loading ME dataset....
 +              Dataset dimensions: [Nx=41,Ny=52,Nz=28,Ne=3,Nt=151]
++ INFO [Main]: Generating initial mask from data.
 +              Dataset ./SBJ01_Sess01_Task01.mask.orig.nii written to disk
 +              Number of Voxels in mask [Nv=27105]
++ INFO [Main]: Computing Mean across time for all echoes....
 +              Dataset ./SBJ01_Sess01_Task01.SME.mean.nii written to disk
++ INFO [Main]: Quality Assurance....
 +              Computing QA metrics from data.
 +              Dataset /Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/DXX_TestingJavierMeica1/SBJ01_Sess01_Task01.QAC.SSE_rank.nii written to disk
 +              Dataset /Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/DXX_TestingJavierMeica1/SBJ01_Sess01_Task01.QAC.SSE.nii written to disk
++ INFO [Main]: Compute Static S0 and R2* maps...
 +              Computing static S0 and t2s maps from the data.
 + INFO [make_static_maps_opt]: Initial conditions [So=2500, T2s=40]
 + INFO [make_static_maps_opt]: Bounds So=[100,10000] & T2s=[10,300]
 + INFO [make_static_maps_opt]: Optimizer = SLSQP
 +              Multi-process Static Map Fit -> Ncpu = 4
 + INFO [make_static_maps_opt]: Number of Voxels with errors: 27060
 +              Dataset /Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/DXX_TestingJavierMeica1/SBJ01_Sess01_Task01.sTE.S0.nii written to disk
 +              Dataset /Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/DXX_TestingJavierMeica1/SBJ01_Sess01_Task01.sTE.t2s.nii written to disk
 +              Dataset /Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/DXX_TestingJavierMeica1/SBJ01_Sess01_Task01.sTE.SSE.nii written to disk
 +              Dataset /Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/DXX_TestingJavierMeica1/SBJ01_Sess01_Task01.sTE.mask.bad.nii written to disk
++ INFO [Main]: Remove non-S0/T2* signals....
 +              OFF
++ INFO [Main]: Compute Optimally Combined Time series...
 +              Compute Optimally Combined Time series from data.
 +              Dataset /Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/DXX_TestingJavierMeica1/ts_OC.nii written to disk
++ INFO [Main]: Removal of Gaussian Noise...
 +              Performing PCA on optimally combined timeseries.
 +              Computing number of PCA components in dataset.
 +              M = 151, N = 27105
 +              Beta = 0.00557092787309, omega = 1.44010970207
 +              Median y = 5128.58221278, Tau = 7385.72100251
 +              Number of components = 40
 +              Signal power: 0.927520821061
 +              Noise power: 0.0724791789394
 +              Writing OCTS after PCA.
 +              Dataset /Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/DXX_TestingJavierMeica1/SBJ01_Sess01_Task01.PCA.octs.nii written to disk
 +              Dataset /Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/DXX_TestingJavierMeica1/SBJ01_Sess01_Task01.PCA.outts.nii written to disk
 +              Serializing, compressing and writting to disk of PCA output objects.
 +              Var Explained after PCA: 92.75
++ INFO [Main]: ICA Decomposition...
 +              Perform ICA....
 +              Correct the sign of the ICA components....
 +              Saving ICA maps as /Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/DXX_TestingJavierMeica1/SBJ01_Sess01_Task01.ICA.Zmaps.nii
 +              Saving ICA representative timeseries as /Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/DXX_TestingJavierMeica1/meica_mix.1D
 +              Saving Z-scores ICA representative timeseries as /Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/DXX_TestingJavierMeica1/meica_mix_zsc.1D
 +              Dataset /Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/DXX_TestingJavierMeica1/SBJ01_Sess01_Task01.ICA.Zmaps.nii written to disk
 +              Saving the Fourier Transform of the ICA timeseries as /Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/DXX_TestingJavierMeica1/meica_mix_fft.1D
 +              Writing OCTS after ICA [/Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/DXX_TestingJavierMeica1/SBJ01_Sess01_Task01.ICA.octs.nii]
 +              Dataset /Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/DXX_TestingJavierMeica1/SBJ01_Sess01_Task01.ICA.octs.nii written to disk
 +              Var Explained after ICA: 92.75
++ INFO [Main]: Computing Kappa, Rho and Variance Explained...
++ INFO [Main]: QA weigths were used during kappa/rho computation
 +              Dataset ./SBJ01_Sess01_Task01.ICA.Zmaps.mask.nii written to disk
 +              Multi-process Characterize Components -> Ncpu = 4
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta000.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta003.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta006.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta009.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit000.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit003.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit006.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit009.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit000.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit003.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit006.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit009.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta001.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta004.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta007.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta010.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit001.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit004.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit007.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit010.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit001.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit004.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit007.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit010.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta002.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta005.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta008.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta011.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit002.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit005.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit008.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit011.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit002.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit005.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit008.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit011.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta012.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta015.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta018.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta021.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit012.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit015.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit018.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit021.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit012.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit015.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit018.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit021.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta013.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta016.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta019.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta022.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit013.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit016.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit019.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit022.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit013.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit016.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit019.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit022.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta014.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta017.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta020.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta023.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit014.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit017.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit020.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit023.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit014.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit017.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit020.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit023.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta024.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta027.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta030.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta033.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit024.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit027.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit030.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit033.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit024.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit027.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit030.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit033.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta025.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta028.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta031.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta034.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit025.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit028.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit031.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit034.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit025.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit028.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit031.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit034.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta026.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta029.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta032.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta035.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit026.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit029.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit035.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit032.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit026.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit029.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit032.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit035.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta036.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta039.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit036.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit039.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit036.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit039.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta037.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit037.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit037.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.Beta038.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.S0Fit038.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.EXTRA.R2Fit038.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.Beta.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.FS0.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.FR2.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.FS0.mask.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.FR2.mask.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.cS0.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.cR2.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.pS0.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.pR2.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.Kappa.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.Kappa_mask.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.Rho_mask.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.Rho.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.chComp.weightMaps.nii written to disk
    compID  Kappa    Rho  varExp  maxR2  maxS0  Ratio  maxZICA  NvZmask  \
0     0.00 130.57  87.03    2.07 500.00 122.49   1.50     3.43    35.00
1     1.00  53.47 113.95    2.69 500.00 216.68   0.47     1.35    35.00
2     2.00 227.07  24.10    2.49 315.64  77.43   9.42     2.55    32.00
3     3.00  74.03  46.89    3.01 143.15 102.53   1.58     2.21    34.00
4     4.00  54.60  82.80    2.64 165.02 323.70   0.66     2.65    37.00
5     5.00 163.33 166.74    2.89 500.00 500.00   0.98     3.28    35.00
6     6.00  72.70  30.08    2.19 500.00  68.90   2.42     0.81    35.00
7     7.00  71.60 128.75    3.32 500.00 262.68   0.56     1.90    31.00
8     8.00  97.41  65.38    2.57 500.00 500.00   1.49     2.03    33.00
9     9.00  82.88  27.30    1.34 500.00  48.87   3.04     0.84    31.00
10   10.00  54.79 394.05    2.17 175.01 500.00   0.14     1.20    34.00
11   11.00  92.62  60.77    1.92 500.00 500.00   1.52     3.89    32.00
12   12.00 144.35  25.19    2.46 494.33 500.00   5.73     3.75    33.00
13   13.00 167.70 265.82    2.44 500.00 500.00   0.63     5.12    40.00
14   14.00 235.74  76.79    2.17 383.23  89.53   3.07     0.91    33.00
15   15.00  65.24  54.35    2.02 120.46 500.00   1.20     4.09    38.00
16   16.00  32.03 192.12    2.45 322.03 500.00   0.17     1.66    36.00
17   17.00 135.64  57.83    2.67 500.00  88.76   2.35     4.19    37.00
18   18.00 228.25  47.42    2.10 500.00 176.92   4.81     1.90    33.00
19   19.00  32.94  70.79    3.04  40.23  71.29   0.47     3.70    33.00
20   20.00  43.43  37.07    2.61 137.64  74.99   1.17     1.02    35.00
21   21.00  82.73  87.06    2.08 239.54 125.61   0.95     4.55    28.00
22   22.00  58.27  23.22    2.26  63.58  23.22   2.51     0.88    34.00
23   23.00 221.14  52.54    1.43 500.00  83.60   4.21     1.26    26.00
24   24.00 214.38  29.24    3.17 500.00  35.31   7.33     2.41    34.00
25   25.00  76.49  26.55    2.58 249.75  30.85   2.88     2.58    34.00
26   26.00 122.37  63.44    2.59 155.75 226.78   1.93     1.78    33.00
27   27.00  44.52  33.42    2.73  75.56  83.44   1.33     2.15    33.00
28   28.00 131.60  37.70    2.57 500.00 100.40   3.49     0.93    30.00
29   29.00  32.03  27.75    3.21  60.45  28.70   1.15     2.15    36.00
30   30.00  66.05  47.64    3.05 159.54 260.01   1.39     2.72    41.00
31   31.00  97.27  76.75    2.47 169.08 500.00   1.27     1.97    37.00
32   32.00 102.43 145.95    2.42 500.00 242.34   0.70     2.23    34.00
33   33.00  37.66 192.18    2.38  81.57 500.00   0.20     2.73    31.00
34   34.00  92.73 118.45    2.82 151.74 500.00   0.78     2.44    34.00
35   35.00  24.31  79.29    2.72  36.60 351.58   0.31     1.63    30.00
36   36.00 148.71 158.11    2.23 500.00 500.00   0.94     1.52    35.00
37   37.00  83.22  60.62    2.95 319.23 500.00   1.37     2.21    35.00
38   38.00  41.85 106.68    2.43 159.45 160.79   0.39     1.66    28.00
39   39.00  43.07 141.56    2.66 281.69 381.48   0.30     1.72    32.00

    NvFR2mask  NvFS0mask  NvKapMask  NvRhoMask
0    4,149.00   2,713.00       6.00       6.00
1    5,055.00   2,763.00       7.00       6.00
2    6,721.00   3,403.00       7.00       4.00
3    4,366.00   2,744.00       5.00       3.00
4    4,816.00   2,761.00      10.00       7.00
5    5,836.00   3,178.00      10.00       7.00
6    5,001.00   2,809.00       6.00       2.00
7    6,453.00   3,465.00       8.00       3.00
8    5,283.00   2,880.00       7.00       3.00
9    4,984.00   2,899.00      11.00       7.00
10   5,262.00   2,947.00       8.00       5.00
11   4,099.00   2,357.00       8.00       6.00
12   5,372.00   2,877.00       6.00       3.00
13   6,111.00   4,023.00      12.00      22.00
14   5,235.00   3,000.00       6.00       2.00
15   4,285.00   2,615.00       7.00       9.00
16   4,993.00   2,699.00       6.00       3.00
17   4,780.00   2,836.00       9.00       5.00
18   4,466.00   2,670.00       5.00       7.00
19   3,664.00   1,999.00       3.00       2.00
20   5,076.00   2,788.00       3.00       4.00
21   4,194.00   2,671.00       6.00       3.00
22   3,693.00   2,414.00       3.00       1.00
23   4,455.00   2,480.00       6.00       4.00
24   4,600.00   2,793.00       5.00       3.00
25   4,063.00   2,597.00      12.00       5.00
26   5,057.00   2,792.00       4.00       6.00
27   4,486.00   2,761.00       3.00       5.00
28   6,615.00   3,231.00      10.00       4.00
29   5,185.00   3,004.00       4.00       4.00
30   4,492.00   2,769.00      10.00      10.00
31   4,963.00   3,045.00      10.00      12.00
32   4,580.00   3,147.00       8.00       6.00
33   5,269.00   4,032.00       6.00      10.00
34   4,913.00   2,731.00       5.00       5.00
35   5,933.00   3,266.00       9.00       4.00
36   5,479.00   2,894.00      14.00      15.00
37   5,076.00   2,838.00       5.00       7.00
38   4,425.00   3,036.00       6.00       4.00
39   5,455.00   3,241.00      12.00      13.00
++ INFO [Main]: Component Selection...
 +              Unsorted List of accepted components=[ 2  6  9 12 14 17 18 22 23 24 25 28]
 +              Variance criteria not active [OFF].
++ INFO [Main]: Generating output time series...
 +              Dataset ./SBJ01_Sess01_Task01.OUT.gTS.nii written to disk
 +              Dataset ./SBJ01_Sess01_Task01.OUT.bTS.nii written to disk
 +              Dataset ./betas_OC.nii written to disk
 +              Dataset ./dn_ts_OC.nii written to disk
++ INFO [Main]: Successfull Completion of the Analysis.
++ INFO [Main]: =======================================




 python /Users/handwerkerd/Documents/Code/SFIM_ME/report/meica_report.py -TED ./ -motion SBJ01_S01_Task01_e2_Motion.1D -dir ./ -ax -sag -title SBJ01_Sess01_Task01 -setname ./ -label SBJ01_Sess01_Task01 -sort_col 7 -overwrite
++ INFO: Checking system for dependencies...
++ INFO: Sphinx version: 1.3.1
++ INFO: Numpy version: 1.10.1
++ INFO: Dependencies OK.
++ 3daxialize: AFNI version=AFNI_2011_12_21_1014 (Mar 23 2015) [64-bit]
++ 3daxialize: AFNI version=AFNI_2011_12_21_1014 (Mar 23 2015) [64-bit]
++ 3daxialize: AFNI version=AFNI_2011_12_21_1014 (Mar 23 2015) [64-bit]
/Users/handwerkerd/anaconda/lib/python2.7/site-packages/numpy/lib/npyio.py:891: UserWarning: loadtxt: Empty input file: "/Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/DXX_TestingJavierMeica1/midk_rejected.txt"
  warnings.warn('loadtxt: Empty input file: "%s"' % fname)
    compID    Kappa     Rho  varExp  varExpN    maxR2    maxS0  Ratio  \
0        2  227.072  24.098   2.489    2.489  315.643   77.432  9.423
1       24  214.382  29.238   3.168    3.168  500.000   35.310  7.332
2       12  144.347  25.191   2.461    2.461  494.326  500.000  5.730
3       18  228.245  47.422   2.102    2.102  500.000  176.918  4.813
4       23  221.142  52.536   1.433    1.433  500.000   83.602  4.209
5       28  131.600  37.705   2.567    2.567  500.000  100.404  3.490
6       14  235.735  76.787   2.167    2.167  383.228   89.532  3.070
7        9   82.880  27.295   1.336    1.336  500.000   48.874  3.036
8       25   76.495  26.554   2.576    2.576  249.746   30.852  2.881
9       22   58.270  23.217   2.257    2.257   63.580   23.217  2.510
10       6   72.704  30.076   2.185    2.185  500.000   68.903  2.417
11      17  135.641  57.829   2.666    2.666  500.000   88.757  2.346

    maxZICA  NvZmask  NvFR2mask  NvFS0mask  NvKapMask  NvRhoMask
0     2.553       32       6721       3403          7          4
1     2.411       34       4600       2793          5          3
2     3.746       33       5372       2877          6          3
3     1.896       33       4466       2670          5          7
4     1.256       26       4455       2480          6          4
5     0.925       30       6615       3231         10          4
6     0.908       33       5235       3000          6          2
7     0.837       31       4984       2899         11          7
8     2.584       34       4063       2597         12          5
9     0.883       34       3693       2414          3          1
10    0.813       35       5001       2809          6          2
11    4.194       37       4780       2836          9          5
++ INFO: Making figures
   Making Kappa vs component number plot
[14 18  2 23 24 13  5 36 12 17 28  0 26 32  8 31 34 11 37  9 21 25  3  6  7
 30 15 22 10  4  1 27 20 39 38 33 19 16 29 35]
   Making Kappa vs Rho plot
/Users/handwerkerd/anaconda/lib/python2.7/site-packages/matplotlib/collections.py:590: FutureWarning: elementwise comparison failed; returning scalar instead, but in the future will perform elementwise comparison
  if self._edgecolors == str('face'):
   Making TSNR plots
   Making motion plots
   Making component images.  This set of figures may take awhile...
    Figures created for Component 00
    Figures created for Component 02
    Figures created for Component 04
    Figures created for Component 06
    Figures created for Component 08
    Figures created for Component 10
    Figures created for Component 12
    Figures created for Component 14
    Figures created for Component 07
    Figures created for Component 01
    Figures created for Component 03
    Figures created for Component 15
    Figures created for Component 11
    Figures created for Component 05
    Figures created for Component 13
    Figures created for Component 09
    Figures created for Component 16
    Figures created for Component 18
    Figures created for Component 20
    Figures created for Component 22
    Figures created for Component 24
    Figures created for Component 26
    Figures created for Component 28
    Figures created for Component 30
    Figures created for Component 19
    Figures created for Component 21
    Figures created for Component 17
    Figures created for Component 23
    Figures created for Component 25
    Figures created for Component 27
    Figures created for Component 29
    Figures created for Component 31
    Figures created for Component 32
    Figures created for Component 34
    Figures created for Component 36
    Figures created for Component 38
    Figures created for Component 35
    Figures created for Component 33
    Figures created for Component 37
    Figures created for Component 39
++ Error: No anatomical specified, cannot create coregistration or correlation maps
++ INFO: Occupying sphinx directory with .rst files
sphinx-build -b html -d _build/doctrees   . _build/html
Running Sphinx v1.3.1
making output directory...
loading pickled environment... not yet created
WARNING: 'default' html theme has been renamed to 'classic'. Please change your html_theme setting either to the new 'alabaster' default theme, or to 'classic' to keep using the old default.
building [mo]: targets for 0 po files that are out of date
building [html]: targets for 4 source files that are out of date
updating environment: 4 added, 0 changed, 0 removed
reading sources... [100%] intro
/Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/DXX_TestingJavierMeica1/meica.SBJ01_Sess01_Task01/analysis.rst:883: ERROR: Unexpected indentation.
/Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/DXX_TestingJavierMeica1/meica.SBJ01_Sess01_Task01/analysis.rst:7: ERROR: Undefined substitution referenced: "warning".
/Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/DXX_TestingJavierMeica1/meica.SBJ01_Sess01_Task01/analysis.rst:881: ERROR: Undefined substitution referenced: "warning".
looking for now-outdated files... none found
pickling environment... done
checking consistency... done
preparing documents... done
writing output... [100%] intro
generating indices... genindex
writing additional pages... search
copying images... [100%] Report_Figures/Axial_Component_01.png
copying static files... WARNING: html_static_path entry u'/Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/DXX_TestingJavierMeica1/meica.SBJ01_Sess01_Task01/_static' does not exist
done
copying extra files... done
dumping search index in English (code: en) ... done
dumping object inventory... done
build succeeded, 5 warnings.

Build finished. The HTML pages are in _build/html.
sphinx-build -b latex -d _build/doctrees   . _build/latex
Running Sphinx v1.3.1
making output directory...
loading pickled environment... done
building [mo]: targets for 0 po files that are out of date
building [latex]: all documents
updating environment: 0 added, 0 changed, 0 removed
looking for now-outdated files... none found
processing MeicaReport.tex... index intro diagnostics analysis
resolving references...
writing... /Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/DXX_TestingJavierMeica1/meica.SBJ01_Sess01_Task01/analysis.rst:: WARNING: unusable reference target found: ../Report_Figures/kappa_vs_rho.html
done
copying images... Report_Figures/Axial_Component_32.png Report_Figures/FFT_Component_34.png Report_Figures/Axial_Component_22.png Report_Figures/FFT_Component_31.png Report_Figures/TimeSeries_28.png Report_Figures/TimeSeries_01.png Report_Figures/TimeSeries_00.png Report_Figures/Sagittal_Component_32.png Report_Figures/Sagittal_Component_13.png Report_Figures/Axial_Component_36.png Report_Figures/Sagittal_Component_02.png Report_Figures/FFT_Component_24.png Report_Figures/FFT_Component_03.png Report_Figures/TimeSeries_23.png Report_Figures/Axial_Component_30.png Report_Figures/Sagittal_Component_12.png Report_Figures/Axial_Component_11.png Report_Figures/Axial_Component_34.png Report_Figures/FFT_Component_32.png Report_Figures/TimeSeries_39.png Report_Figures/TimeSeries_07.png Report_Figures/Axial_Component_07.png Report_Figures/Axial_Component_19.png Report_Figures/TimeSeries_15.png Report_Figures/Sagittal_Component_18.png Report_Figures/Sagittal_Component_24.png Report_Figures/FFT_Component_16.png Report_Figures/FFT_Component_26.png Report_Figures/Axial_Component_39.png Report_Figures/Axial_Component_02.png Report_Figures/TimeSeries_11.png Report_Figures/Axial_Component_09.png Report_Figures/Sagittal_Component_34.png Report_Figures/Sagittal_Component_19.png Report_Figures/TimeSeries_33.png Report_Figures/FFT_Component_04.png Report_Figures/TimeSeries_32.png Report_Figures/motion_rate.png Report_Figures/kappa_rho_vs_components.png Report_Figures/FFT_Component_33.png Report_Figures/TimeSeries_12.png Report_Figures/Axial_Component_38.png Report_Figures/Axial_Component_18.png Report_Figures/FFT_Component_38.png Report_Figures/Axial_Component_04.png Report_Figures/Sagittal_Component_21.png Report_Figures/Sagittal_Component_30.png Report_Figures/Sagittal_Component_33.png Report_Figures/Axial_Component_14.png Report_Figures/TimeSeries_22.png Report_Figures/Axial_Component_08.png Report_Figures/Axial_Component_10.png Report_Figures/TimeSeries_26.png Report_Figures/TimeSeries_36.png Report_Figures/TimeSeries_06.png Report_Figures/Sagittal_Component_08.png Report_Figures/Sagittal_Component_09.png Report_Figures/Sagittal_Component_20.png Report_Figures/Axial_Component_16.png Report_Figures/FFT_Component_17.png Report_Figures/Sagittal_Component_03.png Report_Figures/FFT_Component_25.png Report_Figures/TimeSeries_21.png Report_Figures/TimeSeries_20.png Report_Figures/Axial_Component_12.png Report_Figures/TimeSeries_24.png Report_Figures/TimeSeries_05.png Report_Figures/TimeSeries_13.png Report_Figures/FFT_Component_27.png Report_Figures/FFT_Component_29.png Report_Figures/Sagittal_Component_26.png Report_Figures/Sagittal_Component_01.png Report_Figures/Sagittal_Component_04.png Report_Figures/Sagittal_Component_05.png Report_Figures/Axial_Component_21.png Report_Figures/Axial_Component_33.png Report_Figures/Axial_Component_25.png Report_Figures/Axial_Component_29.png Report_Figures/TimeSeries_16.png Report_Figures/Sagittal_Component_37.png Report_Figures/Axial_Component_05.png Report_Figures/FFT_Component_19.png Report_Figures/Sagittal_Component_25.png Report_Figures/Sagittal_Component_00.png Report_Figures/kappa_vs_rho.png Report_Figures/Sagittal_Component_10.png Report_Figures/Sagittal_Component_14.png Report_Figures/Axial_Component_23.png Report_Figures/TimeSeries_19.png Report_Figures/Axial_Component_26.png Report_Figures/TimeSeries_04.png Report_Figures/TimeSeries_10.png Report_Figures/Axial_Component_35.png Report_Figures/Sagittal_Component_31.png Report_Figures/FFT_Component_06.png Report_Figures/FFT_Component_07.png Report_Figures/Sagittal_Component_17.png Report_Figures/Sagittal_Component_16.png Report_Figures/TimeSeries_34.png Report_Figures/TimeSeries_08.png Report_Figures/TimeSeries_03.png Report_Figures/Sagittal_Component_22.png Report_Figures/Sagittal_Component_23.png Report_Figures/FFT_Component_28.png Report_Figures/FFT_Component_36.png Report_Figures/Sagittal_Component_29.png Report_Figures/Sagittal_Component_39.png Report_Figures/motion_plot.png Report_Figures/Sagittal_Component_27.png Report_Figures/FFT_Component_10.png Report_Figures/Sagittal_Component_06.png Report_Figures/TimeSeries_31.png Report_Figures/Axial_Component_31.png Report_Figures/TimeSeries_27.png Report_Figures/TimeSeries_35.png Report_Figures/Axial_Component_03.png Report_Figures/Axial_Component_00.png Report_Figures/FFT_Component_39.png Report_Figures/Sagittal_Component_35.png Report_Figures/tsoc_tsnr.png Report_Figures/tsoc_tsnr_hist.png Report_Figures/FFT_Component_22.png Report_Figures/medn_tsnr_hist.png Report_Figures/FFT_Component_12.png Report_Figures/FFT_Component_11.png Report_Figures/FFT_Component_37.png Report_Figures/Axial_Component_24.png Report_Figures/TimeSeries_37.png Report_Figures/TimeSeries_02.png Report_Figures/FFT_Component_08.png Report_Figures/FFT_Component_05.png Report_Figures/FFT_Component_20.png Report_Figures/FFT_Component_21.png Report_Figures/FFT_Component_13.png Report_Figures/Sagittal_Component_38.png Report_Figures/Sagittal_Component_15.png Report_Figures/Axial_Component_20.png Report_Figures/FFT_Component_35.png Report_Figures/Axial_Component_15.png Report_Figures/FFT_Component_09.png Report_Figures/tsnr_ratio_hist.png Report_Figures/Sagittal_Component_11.png Report_Figures/FFT_Component_02.png Report_Figures/Sagittal_Component_07.png Report_Figures/tsnr_ratio.png Report_Figures/FFT_Component_30.png Report_Figures/TimeSeries_29.png Report_Figures/FFT_Component_01.png Report_Figures/FFT_Component_15.png Report_Figures/FFT_Component_14.png Report_Figures/FFT_Component_00.png Report_Figures/Axial_Component_17.png Report_Figures/TimeSeries_09.png Report_Figures/TimeSeries_25.png Report_Figures/Axial_Component_28.png Report_Figures/TimeSeries_17.png Report_Figures/TimeSeries_14.png Report_Figures/FFT_Component_18.png Report_Figures/Sagittal_Component_28.png Report_Figures/Sagittal_Component_36.png Report_Figures/Axial_Component_37.png Report_Figures/TimeSeries_30.png Report_Figures/TimeSeries_18.png Report_Figures/medn_tsnr.png Report_Figures/Axial_Component_27.png Report_Figures/TimeSeries_38.png Report_Figures/Axial_Component_13.png Report_Figures/Axial_Component_06.png Report_Figures/FFT_Component_23.png Report_Figures/Axial_Component_01.png
copying TeX support files...
done
build succeeded, 1 warning.

Build finished; the LaTeX files are in _build/latex.
Run `make' in that directory to run these through (pdf)latex (use `make latexpdf' here to do that automatically).
mv: rename /Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/DXX_TestingJavierMeica1/meica.SBJ01_Sess01_Task01/_static/* to /Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/DXX_TestingJavierMeica1/meica.SBJ01_Sess01_Task01/*: No such file or directory



python /Users/handwerkerd/Documents/Code/SFIM_ME/report_bokeh/meica_report.py -meicaDir ./ -runID SBJ01_Sess01_Task01 -picDir ./meica.SBJ01_Sess01_Task01/Report_Figures/
++ INFO [Main]: Using meBasic library located in: /Users/handwerkerd/Documents/Code/SFIM_ME/melib

++ INFO [Main]: ME-ICA Output Directory is /Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/DXX_TestingJavierMeica1
++ INFO [Main]: Path to PNG ICA files is /Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/DXX_TestingJavierMeica1/meica.SBJ01_Sess01_Task01/Report_Figures
++ INFO [Main]: Number of bins for histograms 200
++ INFO [Main]: FR2 Maps [./SBJ01_Sess01_Task01.chComp.FR2.nii] loaded successfully. [(41, 52, 28, 40)]
++ INFO [Main]: FR2 Maps [./SBJ01_Sess01_Task01.chComp.cR2.nii] loaded successfully. [(41, 52, 28, 40)]
++ INFO [Main]: FS0 Maps [./SBJ01_Sess01_Task01.chComp.FS0.nii] loaded successfully. [(41, 52, 28, 40)]
++ INFO [Main]: FS0 Maps [./SBJ01_Sess01_Task01.chComp.cS0.nii] loaded successfully. [(41, 52, 28, 40)]
++ INFO [Main]: Kappa Masks [./SBJ01_Sess01_Task01.chComp.Kappa_mask.nii] loaded successfully.
++ INFO [Main]: Rho Masks [./SBJ01_Sess01_Task01.chComp.Rho_mask.nii] loaded successfully.
++ INFO [Main]: Kappa Maps [./SBJ01_Sess01_Task01.chComp.Kappa.nii] loaded successfully.
++ INFO [Main]: Rho Maps [./SBJ01_Sess01_Task01.chComp.Rho.nii] loaded successfully.
++ INFO [Main]: ICA Maps [./SBJ01_Sess01_Task01.ICA.Zmaps.nii] loaded successfully.
++ INFO [Main]: ICA Masks [./SBJ01_Sess01_Task01.ICA.Zmaps.mask.nii] loaded successfully.
++ INFO [Main]: ICA Masks [./SBJ01_Sess01_Task01.mask.orig.nii] loaded successfully.
++ INFO [Main]: ICA Masks [./SBJ01_Sess01_Task01.chComp.weightMaps.nii] loaded successfully.
/Users/handwerkerd/Documents/Code/SFIM_ME/report_bokeh/AxSag_Component_XX.png
++ INFO: Number of ICA Maps: 40
No handlers could be found for logger "/Users/handwerkerd/anaconda/lib/python2.7/site-packages/bokeh/validation/check.pyc"




# Javier's pre-processed version
cd //Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/Data4Dan/SBJ01FCH/D03_MEICA

python /Users/handwerkerd/Documents/Code/SFIM_ME/report_bokeh/meica_report.py -meicaDir ./ -runID SBJ01FCH_S01Run01 -picDir ./meica.SBJ01FCH_S01Run01/Report_Figures/

++ INFO [Main]: Using meBasic library located in: /Users/handwerkerd/Documents/Code/SFIM_ME/melib
++ INFO [Main]: ME-ICA Output Directory is /Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/Data4Dan/SBJ01FCH/D03_MEICA
++ INFO [Main]: Path to PNG ICA files is /Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/Data4Dan/SBJ01FCH/D03_MEICA/meica.SBJ01FCH_S01Run01/Report_Figures
++ INFO [Main]: Number of bins for histograms 200
++ INFO [Main]: FR2 Maps [./SBJ01FCH_S01Run01.chComp.FR2.nii] loaded successfully. [(64, 64, 33, 39)]
++ INFO [Main]: FR2 Maps [./SBJ01FCH_S01Run01.chComp.cR2.nii] loaded successfully. [(64, 64, 33, 39)]
++ INFO [Main]: FS0 Maps [./SBJ01FCH_S01Run01.chComp.FS0.nii] loaded successfully. [(64, 64, 33, 39)]
++ INFO [Main]: FS0 Maps [./SBJ01FCH_S01Run01.chComp.cS0.nii] loaded successfully. [(64, 64, 33, 39)]
++ INFO [Main]: Kappa Masks [./SBJ01FCH_S01Run01.chComp.Kappa_mask.nii] loaded successfully.
++ INFO [Main]: Rho Masks [./SBJ01FCH_S01Run01.chComp.Rho_mask.nii] loaded successfully.
++ INFO [Main]: Kappa Maps [./SBJ01FCH_S01Run01.chComp.Kappa.nii] loaded successfully.
++ INFO [Main]: Rho Maps [./SBJ01FCH_S01Run01.chComp.Rho.nii] loaded successfully.
++ INFO [Main]: ICA Maps [./SBJ01FCH_S01Run01.ICA.Zmaps.nii] loaded successfully.
++ INFO [Main]: ICA Masks [./SBJ01FCH_S01Run01.ICA.Zmaps.mask.nii] loaded successfully.
++ INFO [Main]: ICA Masks [./SBJ01FCH_S01Run01.mask.orig.nii] loaded successfully.
++ INFO [Main]: ICA Masks [./SBJ01FCH_S01Run01.chComp.weightMaps.nii] loaded successfully.
/Users/handwerkerd/Documents/Code/SFIM_ME/report_bokeh/AxSag_Component_XX.png
++ INFO: Number of ICA Maps: 39
No handlers could be found for logger "/Users/handwerkerd/anaconda/lib/python2.7/site-packages/bokeh/validation/check.pyc"


cd //Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SBJ01_S01/DXX_TestingJavierMeica1

python /Users/handwerkerd/Documents/Code/SFIM_ME/report_bokeh/meica_report.py -meicaDir ./ -runID SBJ01_Sess01_Task01 -picDir ./meica.SBJ01_Sess01_Task01/Report_Figures/


