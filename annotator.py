#!/usr/bin/env python
"""
annotator.py -i INPUT -o OUTPUT -d DEPTH_FIELD -v VARIANT_DEPTH_FIELD [-g <True | False>]

Runs annotator on a specified input VCF, and writes the result to a specified output file.

Args:
    -i  INPUT (str):               path to the input VCF
    -o  OUTPUT (str):              path to the output tab-separated file
    -d  DEPTH_COLUMN (str):        the name of the FORMAT field that specifies total read coverage 
                                   of the variant locus
    -v  VARIANT_DEPTH_FIELD (str): the name of the FORMAT field that specifies variant read 
                                   coverage of the variant locus
    -g  PER_GENE (boolean):        if enabled, annotator will return one gene annotation per 
                                   variant. Otherwise (by default) will return one annotation per
                                   transcript
"""

import argparse

## Import annotator modules
import annotator.vcf
import annotator.vep
import annotator.main


def main():
    "runs the annotator program on the input data if this module is main"
    parser = argparse.ArgumentParser(
        prog = "annotator.py",
        description = (
            "Runs annotator on a specified input VCF, and writes the result to a specified output "
            "file."
        )
    )
    parser.add_argument(
        "-i", 
        "--input",
        help = "Path to input file. It must be a VCF with a properly formatted header.",
        type = str,
        required = True)
    parser.add_argument(
        "-o",
        "--output", 
        help = "Path to the output tab-separated file. If it does not exist, it will be created.",
        type = str,
        required = True)
    parser.add_argument(
        "-d",
        "--total_depth",
        help = "Name of the FORMAT field that describes the total read depth of the variant locus.",
        type = str,
        required = True)
    parser.add_argument(
        "-v",
        "--variant_depth",
        help = "Name of the FORMAT field that describes the variant read depth of the locus.",
        type = str,
        required = True)
    parser.add_argument(
        "-g",
        "--per-gene",
        help = "Should annotator output annotations per gene instead of per transcript?",
        default = False,
        action = argparse.BooleanOptionalAction
    )
    args = parser.parse_args()

    annotator.main.run_annotator(
        args.input,
        args.output,
        args.total_depth,
        args.variant_depth,
        args.per_gene
    )

if __name__ == "__main__":
    main()
