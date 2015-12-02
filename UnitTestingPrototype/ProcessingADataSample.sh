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
# SFIMMEloc: The location of the SFIM_ME directory. (i.e. "python ${SFIMMEloc}/sfim_meica_oc.py" should run a program)
# NumCPUS: The number of parallel threads that can be used
# NumIter: The number of iterations to run for each task condition 
#    (i.e. 3 means run sfim_meica_??.py with Python verison X on the # echo data 3 times)
# Python27Env Python3Env: This program runs data with both Python 2.7 and Python 3.5.
#   For each version "source activate $Python??Env" will be called and then the script will be run
# RootPath: If this is being run on felix/biowulf, this should be:
#    '/data/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SFIM_ME_TestingData' If this is run on a 
#    mounted drive on a Mac, you can use '//Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SFIM_ME_TestingData'

#
# An example call is: 
#  ./ProcessingADataSample.sh Original2015_11_25 /Users/handwerkerd/Documents/Code/SFIM_ME 12 3 p27werker p35werker  /data/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SFIM_ME_TestingData 12
# Running line-by-line D Handwerker's options: OutputPrefix=Original2015_11_25; SFIMMEloc=/Users/handwerkerd/Documents/Code/SFIM_ME; NumCPUS=4; NumIter=3; Python2Env=p27; Python3Env=p3; RootPath='//Volumes/NIMH_SFIM/100RUNS_3Tmultiecho/PrcsData/SFIM_ME_TestingData'
#
# D.A. Handwerker 12/02/2015
#
# Future work: Add more flexability to this code so that users can run different combinations of unit tests
#   Make it possible to change the options used for sfim_meica_OC.py and sfim_meica_SE.py
#   Make some automatic measures of accuracy
#   Learn more about the proper way to set up a unit testing script so that this can be done right!


OutputPrefix=$1
SFIMMEloc=$2
NumCPUS=$3
NumIter=$4
Python2Env=$5
Python3Env=$6
RootPath=$7


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
      python $SFIMMEloc/sfim_meica_OC.py --tes_file ${RootPath}/${echotype}/TES_${echotype}.1D \
           -d ${RootPath}/${echotype}/zcat_BlockDesign_${echotype}.nii.gz --prefix ${OutputPrefix} \
           --krRatio=3 --TR=2.0 --out_dir=./ --use_QAweight --ncpus ${NumCPUS} --save_extra --Fmax 500 --Zmax 8 -z 0.8
      python $SFIMMEloc/report/meica_report.py --TED_dir ./ --motion BlockDesign_${echotype}_motion.1D --ncpus ${NumCPUS} --overwrite
      python $SFIMMEloc/report_bokeh/meica_report.py -meicaDir ./ -picDir meica.Report/Report_Figures/ -runID ${OutputPrefix} -saveReport
  
      cd $RootPath/Outputs/${echotype}/${OutputPrefix}/Python${pythonversion}/meica_SE/iter${iter}
      pwd
      python $SFIMMEloc/sfim_meica_SE.py --tes_file ${RootPath}/${echotype}/TES_${echotype}.1D \
        -d ${RootPath}/${echotype}/zcat_BlockDesign_${echotype}.nii.gz --prefix ${OutputPrefix} \
        --krRatio=3 --TR=2.0 --out_dir=./ --use_QAweight --ncpus ${NumCPUS} --save_extra --Fmax 500 --Zmax 8 -z 0.8
      python $SFIMMEloc/report/meica_report.py --TED_dir ./ --motion BlockDesign_${echotype}_motion.1D --ncpus ${NumCPUS} --overwrite
      python $SFIMMEloc/report_bokeh/meica_report.py  -meicaDir ./ -picDir meica.Report/Report_Figures/ -runID Original2015_11_25 -saveReport	  
    done
  done
done

cd $RootPath/Outputs/
mkdir ${OutputPrefix} 
python $SFIMMEloc/report/meta_report.py \
  -pattern_1 "./Echos3/${OutputPrefix}/Python2/meica_OC/iter?/meica.Report/meica_report.txt" \
  -pattern_2 "./Echos3/${OutputPrefix}/Python3/meica_OC/iter?/meica.Report/meica_report.txt" \
  -pattern_3 "./Echos3/${OutputPrefix}/Python2/meica_SE/iter?/meica.Report/meica_report.txt" \
  -pattern_4 "./Echos3/${OutputPrefix}/Python3/meica_SE/iter?/meica.Report/meica_report.txt" \
  -pattern_5 "./Echos5/${OutputPrefix}/Python2/meica_OC/iter?/meica.Report/meica_report.txt" \
  -pattern_6 "./Echos5/${OutputPrefix}/Python3/meica_OC/iter?/meica.Report/meica_report.txt" \
  -pattern_7 "./Echos5/${OutputPrefix}/Python2/meica_SE/iter?/meica.Report/meica_report.txt" \
  -pattern_8 "./Echos5/${OutputPrefix}/Python3/meica_SE/iter?/meica.Report/meica_report.txt" \
  -dest '.' -label meta.Report.${OutputPrefix} --ncpus 4 -var_component -high_var

