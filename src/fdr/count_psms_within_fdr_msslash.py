"""Count number of PSMs within FDR for msSLASH spectral library searching results"""

__author__ = "Lei Wang"
__description__ = "count psms within fdr for msslash searching results"
__email__ = "wang558@indiana.edu"

import logging
import numpy as np
import pandas as pd
import time

from executor.cli import get_parameters

logger = logging.getLogger(__name__)


def count_psm_within_fdr(
    msslash_filename,
    fdr_cutoff,
    use_second_approach_for_fdr,
    output_df,
    output_decoys,
    output_psms,
    **kwargs
):

    logger.info(
        "use FDR approach: {}".format(
            "[Second Approach] find last FDR that is smaller then preset FDR"
            if use_second_approach_for_fdr
            else "[First Approach] find first FDR that's greater than preset"
        )
    )

    start = time.time()

    df = pd.read_csv(msslash_filename, sep="\t")
    df.sort_values(by=["TopScore"], ascending=False, inplace=True)
    df.reset_index(drop=True, inplace=True)

    end = time.time()
    logger.info(f"finish sorting and reseting in {(end-start):.4f}s")

    if output_df:
        start = time.time()
        df.to_csv(msslash_filename + ".sorted", sep="\t", na_rep="NA", index=False)
        end = time.time()
        logger.info(f"write sorted df to file in {(end-start):.4f}s")

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
    cut_off_title = None
    cut_off_num_msgf_in_lib, cut_off_num_match = 0, 0

    psms_within_fdr = []
    peptides_within_fdr = set()
    psms_within_fdr_size = None
    decoys = []
    decoys_size = None
    pepmass_shift = []
    pepmass_shift_size = None

    for i, row in enumerate(df.itertuples()):
        top_peptide = row.TopPep
        top_score = row.TopScore
        top_title = row.Title

        match = 0
        if 1 == int(row.MsgfPepExistInLib):
            num_msgf_in_lib += 1
            if "YES" == row.Match:
                num_match_n_msgf_in_lib += 1
                pepmass_shift.append(row.MsgfPepMassDelta)
                match = 1

        if pd.isna(top_peptide):
            num_decoy += 1
            decoys.append("{}\t{}\t{}".format(top_title, top_peptide, top_score))
        else:
            num_target += 1
            psms_within_fdr.append(
                "{}\t{}\t{}\t{}".format(top_title, top_peptide, top_score, match)
            )
            peptides_within_fdr.add(top_peptide)

        fdr = num_decoy * 1.0 / num_target

        # Second approach for FDR, i.e. last FDR thats smaller than given.
        if use_second_approach_for_fdr:
            if fdr < fdr_cutoff:
                num_reverse, num_forward = num_decoy, num_target
                num_unique_peptides = len(peptides_within_fdr)
                cut_off_sim, cut_off_fdr = top_score, fdr
                cut_off_pos = i + 1
                cut_off_num_msgf_in_lib = num_msgf_in_lib
                cut_off_num_match = num_match_n_msgf_in_lib
                cut_off_title = row.Title
                psms_within_fdr_size = len(psms_within_fdr)
                decoys_size = len(decoys)
                pepmass_shift_size = len(pepmass_shift)
        else:
            # First approach to calculate FDR.
            if fdr > fdr_cutoff:
                break

    if output_decoys:
        out_file = msslash_filename[: msslash_filename.rfind(".")] + ".decoys"
        logger.info(f"write decoys to {out_file}")
        with open(out_file, "w") as writer:
            writer.write("top_title\ttop_peptide\ttop_score\n")
            writer.write("\n".join(decoys if decoys_size is None else decoys[:decoys_size]))

    if output_psms:
        out_file = msslash_filename[: msslash_filename.rfind(".")] + ".psms"
        logger.info(f"write decoys to {out_file}")
        with open(out_file, "w") as writer:
            writer.write("top_title\ttop_peptide\ttop_score\tmatch\n")
            writer.write(
                "\n".join(
                    psms_within_fdr
                    if psms_within_fdr_size is None
                    else psms_within_fdr[:psms_within_fdr_size]
                )
            )

    if use_second_approach_for_fdr:
        logger.info(f"{fdr_cutoff} FDR, {cut_off_sim} similarity, cut off @{cut_off_pos}L")
        logger.info(f"title at cutoff: {cut_off_title}")
        logger.info(f"{num_forward} targets vs {num_reverse} decoys, fdr: {cut_off_fdr}")
        logger.info(f"{num_unique_peptides} unique peptides within fdr")
        logger.info(
            f"{cut_off_num_msgf_in_lib} msgf identified psms also exist in library, {cut_off_num_match} among those match with library searching")
    else:
        logger.info(f"{fdr_cutoff} FDR, {top_score} similarity, cut off @{i+1}L")
        logger.info(f"{num_target} targets vs {num_decoy} decoys, fdr: {fdr}")
        logger.info(f"{len(peptides_within_fdr)} unique peptides within fdr")
        logger.info(
            f"{num_msgf_in_lib} msgf identified psms also exist in library, {num_match_n_msgf_in_lib} among those match with library searching"
        )

    logger.info("len of pepmass_shift: {}".format(len(pepmass_shift)))
    logger.info("min of pepmass_shift: {}".format(np.min(pepmass_shift)))
    logger.info("max of pepmass_shift: {}".format(np.max(pepmass_shift)))
    logger.info("std of pepmass_shift: {}".format(np.std(pepmass_shift)))


def run(**params):
    count_psm_within_fdr(**params)


def parser_add(parser):
    add = parser.add_argument

    add("--msslash-filename", required=True, help="file holding the msSLASH searching result")
    add("--fdr-cutoff", type=float, default=0.01, help="FDR for target decoy search")
    add(
        "--use-second-approach-for-fdr",
        action="store_true",
        help="whether to use second approach for FDR calc",
    )
    add("--output-df", action="store_true", help="whether to save sorted data frame")
    add("--output-decoys", action="store_true", help="whether to print decoys within FDR")
    add("--output-psms", action="store_true", help="whether to print PSMs within FDR")

    return parser


def cli(args=None):
    kwargs = {"description": __description__}
    params = get_parameters(args, parser_add=parser_add, **kwargs)
    run(**params)


if __name__ == "__main__":
    cli()
