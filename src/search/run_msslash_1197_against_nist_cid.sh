#!/bin/bash

# CHARGE 2
#bash run_msslash.sh \
#  -t "/home/data/wang558/SpectralLibrary/NistSpectralLibrary/cid/human_consensus_final_true_lib.mgf" \
#  -d "/home/data/wang558/SpectralLibrary/NistSpectralLibrary/cid/decoy_library_by_precursor_swap/human_cid_consensus_decoy.mgf" \
#  -i "log" \
#  -c 2 \
#  -e "/home/data/wang558/PXD001197-full/msslash/2_full_PXD001197.mgf" \
#  -p "/home/data/wang558/PXD001197-full/msgf_results/tsv/raw/psm/2_full_PXD001197.uniprot.tsv.psm.i2l" \
#  -m 0

# CHARGE 3
bash run_msslash.sh \
  -t "/home/data/wang558/SpectralLibrary/NistSpectralLibrary/cid/human_consensus_final_true_lib.mgf" \
  -d "/home/data/wang558/SpectralLibrary/NistSpectralLibrary/cid/decoy_library_by_precursor_swap/human_cid_consensus_decoy.mgf" \
  -i "log" \
  -c 3 \
  -e "/home/data/wang558/PXD001197-full/msslash/charge3/3_full_PXD001197.mgf" \
  -p "/home/data/wang558/PXD001197-full/msgf_results/tsv/raw/psm/3_full_PXD001197.uniprot.tsv.psm.i2l" \
  -m 0
