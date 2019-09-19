#!/bin/bash

libpep=/home/data/wang558/NistSpectralLibrary/human_synthetic_hcd_selected.peptides
for charge in 2 3;
do
  psm=/home/data/wang558/PXD001197-full/msgf_results/tsv/raw/psm/${charge}_full_PXD001197.uniprot.tsv.psm.i2l
  bf=/home/wang558/research/msms_analysis/msms_analysis/src/search/dataframes_for_search/df_search_pxd1197c${charge}_against_nist_bruteforce 

  spst=/home/data/wang558/PXD001197-full/raw_spectra/${charge}_full_PXD001197_nospace_mztol_0.05da.xls 

  # Run SpectraST with -s_SP4 useSp4Scoring, http://tools.proteomecenter.org/wiki/index.php?title=Software:SpectraST
  #spst=/home/data/wang558/PXD001197-full/raw_spectra/${charge}_full_PXD001197_nospace_mztol_0.05da_SP4.xls 

  #sed -i '1 s/### Query/Query/' $spst
  for it in 10 20 50 100 125 150 200;
  do
    lsh=/home/wang558/research/msms_analysis/msms_analysis/src/search/dataframes_for_search/df_search_pxd1197c${charge}_against_nist_lsh_${it}iter

    out=analysis/analysis_pxd1197c${charge}_against_nist_bruteforce_n_spectrast_n_lsh${it}iter 
    
    # Align with the results of -s_SP4 with SpectraST.
    #out=analysis/analysis_pxd1197c${charge}_against_nist_bruteforce_n_spectrast_SP4_n_lsh${it}iter 

    nohup python analyze.py  --library_peptides=${libpep} --psm_res=${psm} --bruteforce_res=${bf} --spectrast_res=${spst} --lsh_res=${lsh} --charge_state=${charge} &> $out &
  done
done
