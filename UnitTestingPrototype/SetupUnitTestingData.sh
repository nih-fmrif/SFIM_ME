###################################################################
# Setting up the original directory with the volume files to test
# This moves one block design data set with 3 echos and another with 5 echos along with all necessary
#   files to run SFIM_ME programs into a common location
# This only needs to be done once on felix so there is there's no need to re-run this code

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
  3dcalc -a ../../SBJ01_S01/D01_Version02.AlignByAnat.Cubic/Task01/p04.SBJ01_S01_Task01_e1.align_clp+orig \
	     -b ../../SBJ01_S01/D01_Version02.AlignByAnat.Cubic/Task01/p04.SBJ01_S01_Task01_e2.align.mask+orig \
		 -prefix BlockDesign_Echos3_e1.nii.gz -expr 'ispositive(b)*a'
  3dcalc -a ../../SBJ01_S01/D01_Version02.AlignByAnat.Cubic/Task01/p04.SBJ01_S01_Task01_e2.align_clp+orig \
	     -b ../../SBJ01_S01/D01_Version02.AlignByAnat.Cubic/Task01/p04.SBJ01_S01_Task01_e2.align.mask+orig \
		 -prefix BlockDesign_Echos3_e2.nii.gz -expr 'ispositive(b)*a'
  3dcalc -a ../../SBJ01_S01/D01_Version02.AlignByAnat.Cubic/Task01/p04.SBJ01_S01_Task01_e3.align_clp+orig \
	     -b ../../SBJ01_S01/D01_Version02.AlignByAnat.Cubic/Task01/p04.SBJ01_S01_Task01_e2.align.mask+orig \
		 -prefix BlockDesign_Echos3_e3.nii.gz -expr 'ispositive(b)*a'
  3dZcat -prefix zcat_BlockDesign_Echos3.nii.gz BlockDesign_Echos3_e1.nii.gz BlockDesign_Echos3_e2.nii.gz BlockDesign_Echos3_e3.nii.gz
  
  ln -s ../../SBJ01_S01/D01_Version02.AlignByAnat.Cubic/Task01/SBJ01_S01_Task01_e2_Motion.1D ./BlockDesign_Echos3_Motion.1D
  ln -s ../SBJ01_Anatomy.nii.gz ./
  # At some step in the processing pipeline, these NIFTI images get labelled TLRC, but their actually ORIG
  3drefit -view 'orig' -space ORIG BlockDesign_Echos3_e1.nii.gz
  3drefit -view 'orig' -space ORIG BlockDesign_Echos3_e2.nii.gz
  3drefit -view 'orig' -space ORIG BlockDesign_Echos3_e3.nii.gz
  3drefit -view 'orig' -space ORIG zcat_BlockDesign_Echos3.nii.gz

  cd /data/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SFIM_ME_TestingData/Echos5
  3dcalc -a ../../SBJ01_S09/D01_Version02.AlignByAnat.Cubic/Task11/p04.SBJ01_S09_Task11_e1.align_clp+orig \
	     -b ../../SBJ01_S09/D01_Version02.AlignByAnat.Cubic/Task11/p04.SBJ01_S09_Task11_e2.align.mask+orig \
		 -prefix BlockDesign_Echos5_e1.nii.gz -expr 'ispositive(b)*a'
  3dcalc -a ../../SBJ01_S09/D01_Version02.AlignByAnat.Cubic/Task11/p04.SBJ01_S09_Task11_e2.align_clp+orig \
	     -b ../../SBJ01_S09/D01_Version02.AlignByAnat.Cubic/Task11/p04.SBJ01_S09_Task11_e2.align.mask+orig \
		 -prefix BlockDesign_Echos5_e2.nii.gz -expr 'ispositive(b)*a'
  3dcalc -a ../../SBJ01_S09/D01_Version02.AlignByAnat.Cubic/Task11/p04.SBJ01_S09_Task11_e3.align_clp+orig \
	     -b ../../SBJ01_S09/D01_Version02.AlignByAnat.Cubic/Task11/p04.SBJ01_S09_Task11_e2.align.mask+orig \
		 -prefix BlockDesign_Echos5_e3.nii.gz -expr 'ispositive(b)*a'
  3dcalc -a ../../SBJ01_S09/D01_Version02.AlignByAnat.Cubic/Task11/p04.SBJ01_S09_Task11_e4.align_clp+orig \
	     -b ../../SBJ01_S09/D01_Version02.AlignByAnat.Cubic/Task11/p04.SBJ01_S09_Task11_e2.align.mask+orig \
		 -prefix BlockDesign_Echos5_e4.nii.gz -expr 'ispositive(b)*a'
  3dcalc -a ../../SBJ01_S09/D01_Version02.AlignByAnat.Cubic/Task11/p04.SBJ01_S09_Task11_e5.align_clp+orig \
	     -b ../../SBJ01_S09/D01_Version02.AlignByAnat.Cubic/Task11/p04.SBJ01_S09_Task11_e2.align.mask+orig \
		 -prefix BlockDesign_Echos5_e5.nii.gz -expr 'ispositive(b)*a'
  3dZcat -prefix zcat_BlockDesign_Echos5.nii.gz BlockDesign_Echos5_e1.nii.gz BlockDesign_Echos5_e2.nii.gz \
	  BlockDesign_Echos5_e3.nii.gz  BlockDesign_Echos5_e4.nii.gz  BlockDesign_Echos5_e5.nii.gz 
  
  
  ln -s ../../SBJ01_S09/D01_Version02.AlignByAnat.Cubic/Task11/SBJ01_S09_Task11_e2_Motion.1D ./BlockDesign_Echos5_Motion.1D
  ln -s ../SBJ01_Anatomy.nii.gz ./
  # At some step in the processing pipeline, these NIFTI images get labelled TLRC, but their actually ORIG
  3drefit -view 'orig' -space ORIG BlockDesign_Echos5_e1.nii.gz
  3drefit -view 'orig' -space ORIG BlockDesign_Echos5_e2.nii.gz
  3drefit -view 'orig' -space ORIG BlockDesign_Echos5_e3.nii.gz
  3drefit -view 'orig' -space ORIG BlockDesign_Echos5_e4.nii.gz
  3drefit -view 'orig' -space ORIG BlockDesign_Echos5_e5.nii.gz  
  3drefit -view 'orig' -space ORIG zcat_BlockDesign_Echos5.nii.gz

