python /myOpt/METools/V03/jmeica_optcom3.py --tes_file ../D00_OriginalData/SBJ03T01_RestE5A3_Echoes.1D -d ../D02_Preproc/pc04.SBJ03T01_RestE5A3.zcat.data.nii.gz --prefix SBJ03T01_RestE5A3 --krRatio=5 --TR=2.0 --out_dir=./ --use_QAweight --reuse
python /myOpt/METools/V03/report/meica_report.py -TED ./ -motion ../D02_Preproc/SBJ03T01_RestE5A3_Motion.1D -dir ./ -ax -sag -title SBJ03T01_RestE5A3 -setname ./ -label SBJ03T01_RestE5A3 -sort_col 7 -overwrite
python /myOpt/METools/V03/report_bokeh/meica_report.py -meicaDir ./ -runID SBJ03T01_RestE5A3 -picDir ./meica.SBJ03T01_RestE5A3/Report_Figures/
