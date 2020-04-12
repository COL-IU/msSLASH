""" Count number of PSMs within FDR for SpectraST spectral library searching results
    Example:
        python count_psms_within_fdr_spectrast.py --spectrast-filename /home/data/wang558/PXD001197-full/raw_spectra/spectrast_related/tol0.05_against_nist_cid_by_precursor_swap_06302019/3_full_PXD001197_nospace.xls --msgf-psm-filename /home/data/wang558/PXD001197-full/msgf_results/tsv/raw/psm/3_full_PXD001197.uniprot.tsv.psm.i2l --spectral-library-peptide-filename /home/data/wang558/SpectralLibrary/NistSpectralLibrary/cid/peptides/human_consensus_final_true_lib.i2l.peptides_n_charge --transpose-underscore --output-psms --output-df --output-decoys 2>&1 | tee /home/data/wang558/PXD001197-full/raw_spectra/spectrast_related/tol0.05_against_nist_cid_by_precursor_swap_06302019/3_full_PXD001197_nospace.calcfdr.log
"""

__author__ = "Lei Wang"
__description__ = "count psms within fdr for spectrast searching results"
__email__ = "wang558@indiana.edu"

import logging
import pandas as pd
import time

from executor.cli import get_parameters

logger = logging.getLogger(__name__)


def modify_pep(peptide):
    # For spectrast
    peptide = peptide.replace("I", "L")
    peptide = peptide.replace("C[160]", "C(C)")
    peptide = peptide.replace("M[147]", "M(O)")

    # For MSGF identification.
    peptide = peptide.replace("C+57.021", "C(C)")
    peptide = peptide.replace("M+15.995", "M(O)")
    return peptide


# Calculate number of decoy and target, using preset FDR, for SpectraST result.
def count_psm_within_fdr(
    spectrast_filename,
    msgf_psm_filename,
    spectral_library_peptide_filename,
    fdr_cutoff,
    use_second_approach_for_fdr,
    output_df,
    transpose_underscore,
    output_decoys,
    output_psms,
    **kwargs,
):

    logger.info(
        "use FDR approach: {}".format(
            "[Second Approach] find last FDR that is smaller than preset FDR"
            if use_second_approach_for_fdr
            else "[First Approach] find first FDR that is greater than preset"
        )
    )

    df = pd.read_csv(spectrast_filename, sep="\t", index_col=False)
    logger.info(f"read dataframe from {spectrast_filename}, length is {len(df)}")

    try:
        df.rename(columns={"### Query": "Query"}, inplace=True)
        logger.info("replace '### Query' with 'Query' in dataframe")
    except:
        pass

    logger.info("sort and reset index inplace by Dot descendingly.")
    start = time.time()
    df.sort_values(by=["Dot"], ascending=False, inplace=True)
    df.reset_index(drop=True, inplace=True)
    end = time.time()
    logger.info(f"finish sorting and reseting in {(end-start):.4f}s")

    if output_df:
        start = time.time()
        df.to_csv(spectrast_filename + ".sorted", sep="\t", na_rep="NA", index=False)
        end = time.time()
        logger.info(f"write sorted df to file in {(end-start):.4f}s")

    # READ MSGF PSM
    logger.info(f"read MSGF PSM from {msgf_psm_filename}")
    dict_psm = {}
    with open(msgf_psm_filename, "r") as fd:
        for line in fd:
            line = line.strip("\r\n")
            title, peptide = line.split("\t")
            dict_psm[title] = peptide
    logger.info(f"finish reading MSGF PSM file in {(end-start):.4f}s, size of dict is {len(dict_psm)}")

    # READ TARGET LIBRARY *MODIFIED* PEPTIDE/CHARGE
    logger.info(f"read library modified(I2L and mods) peptides w/ charge from {spectral_library_peptide_filename}")
    set_lib_modified_peptides_w_charge = set()
    start = time.time()
    with open(spectral_library_peptide_filename, "r") as fd:
        for line in fd:
            peptide_n_charge = line.strip("\r\n")
            if peptide_n_charge not in set_lib_modified_peptides_w_charge:
                set_lib_modified_peptides_w_charge.add(peptide_n_charge)
            else:
                pass
    end = time.time()
    logger.info(f"finish reading {len(set_lib_modified_peptides_w_charge)} unique modified peptides in {(end-start):.4f}s")

    # First approach, i.e. break when fdr is greater than given FDR.
    num_decoy, num_target = 0, 0
    num_msgf_in_lib = 0
    num_match_n_msgf_in_lib = 0
    top_peptide, top_score = None, None

    # Second approach to calculate FDR, i.e. last FDR thats smaller than given.
    num_reverse, num_forward = 0, 0
    num_unique_peptides = 0
    cut_off_sim, cut_off_fdr = 0, 0
    cut_off_pos = 0
    cut_off_num_msgf_in_lib, cut_off_num_match = 0, 0

    psms_within_fdr, decoys = [], []
    peptides_withini_fdr = set()
    psms_within_fdr_size, decoys_size = None, None

    cnt = 0
    for i, row in enumerate(df.itertuples()):
        top_peptide_w_charge = row.ID
        # I2L, C[160] -> C(C),  M[147] -> M(O)
        top_peptide_w_charge = modify_pep(top_peptide_w_charge)
        top_peptide, charge = top_peptide_w_charge.split("/")

        top_score = row.Dot

        title = row.Query

        # Enable the following replace clause for PXD001197.
        # Disable the following replace clause for PXD000561.
        if transpose_underscore:
            title = title.replace("_", " ")

        msgf_modified_peptide = modify_pep(dict_psm.get(title, ""))
        msgf_modified_peptide_w_charge = "{}/{}".format(msgf_modified_peptide, charge)
        match = 0
        if msgf_modified_peptide_w_charge in set_lib_modified_peptides_w_charge:
            num_msgf_in_lib += 1
            if msgf_modified_peptide == top_peptide:
                num_match_n_msgf_in_lib += 1
                match = 1
        else:
            # logger.info("{} has MSGF ID: {} which does not exist in Spectral Lib".format(title, msgf_modified_peptide))
            pass

        if "DECOY_" in row.Proteins:
            num_decoy += 1
            decoys.append("{}\t{}\t{}".format(title, top_peptide, top_score))
        else:
            num_target += 1
            psms_within_fdr.append("{}\t{}\t{}\t{}".format(title, top_peptide, top_score, match))
            peptides_withini_fdr.add(top_peptide)

        fdr = num_decoy * 1.0 / num_target

        # Second approach to calculate FDR, i.e. last FDR smaller than given.
        if use_second_approach_for_fdr:
            if fdr < fdr_cutoff:
                num_reverse, num_forward = num_decoy, num_target
                num_unique_peptides = len(peptides_withini_fdr)
                cut_off_sim, cut_off_fdr = top_score, fdr
                cut_off_pos = i + 1
                cut_off_num_msgf_in_lib = num_msgf_in_lib
                cut_off_num_match = num_match_n_msgf_in_lib
                psms_within_fdr_size = len(psms_within_fdr)
                decoys_size = len(decoys)
        else:
            # First approach to calculate FDR.
            if fdr > fdr_cutoff:
                break

    if output_decoys:
        out_file = spectrast_filename[: spectrast_filename.rfind(".")] + ".decoys"
        logger.info(f"write decoys to {out_file}")
        with open(out_file, "w") as writer:
            writer.write("top_title\ttop_peptide\ttop_score\n")
            writer.write("\n".join(decoys if decoys_size is None else decoys[:decoys_size]))

    if output_psms:
        out_file = spectrast_filename[: spectrast_filename.rfind(".")] + ".psms"
        logger.info(f"write psms to {out_file}")
        with open(out_file, "w") as writer:
            writer.write("top_title\ttop_peptide\ttop_score\tmatch\n")
            writer.write(
                "\n".join(psms_within_fdr if psms_within_fdr_size is None else psms_within_fdr[:psms_within_fdr_size])
            )

    if use_second_approach_for_fdr:
        logger.info(f"{fdr_cutoff} FDR, {cut_off_sim} similarity, cut off @{cut_off_pos}L")
        logger.info(f"{num_forward} targets vs. {num_reverse} decoys")
        logger.info(f"{num_unique_peptides} unique peptides")
        logger.info(
            f"{cut_off_num_msgf_in_lib} msgf identified psms also exist in library, {cut_off_num_match} among those match with library searching"
        )
    else:
        logger.info(f"{fdr_cutoff} FDR, {top_score} similarity, cut off @{i+1}L")
        logger.info(f"{num_target} targets vs. {num_decoy} decoys")
        logger.info(f"{len(peptides_withini_fdr)} unique peptides")
        logger.info(
            f"{num_msgf_in_lib} msgf identified psms also exist in library, {num_match_n_msgf_in_lib} among those match with library searching"
        )


def run(**params):
    count_psm_within_fdr(**params)
    # calcTargetNDecoyByThreshold(df)


def parser_add(parser):
    add = parser.add_argument

    add("--spectrast-filename", required=True, help="Spectrast file holding the searching result")
    add("--msgf-psm-filename", required=True, help="File holding MSGF PSM results of 1percent FDR")
    add("--spectral-library-peptide-filename", required=True, help="File holding library peptides.")
    add("--fdr-cutoff", type=float, default=0.01, help="FDR for target decoy search.")
    add(
        "--transpose-underscore", action="store_true", help="whether to transpose underscore in title to space",
    )
    add(
        "--use-second-approach-for-fdr", action="store_true", help="whether to use second approach for FDR calc.",
    )
    add("--output-df", action="store_true", help="whether to save sorted data frame.")
    add("--output-decoys", action="store_true", help="whether to print decoys within FDR")
    add("--output-psms", action="store_true", help="whether to print PSMs within FDR")

    return parser


def cli(args=None):
    kwargs = {"description": __description__}
    params = get_parameters(args, parser_add=parser_add, **kwargs)
    run(**params)


if __name__ == "__main__":
    cli()
