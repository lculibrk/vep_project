#!/usr/bin/env python
"""
annotator.py -i INPUT -o OUTPUT -d DEPTH_FIELD -v VARIANT_DEPTH_FIELD

Runs annotator on a specified input VCF, and writes the result to a specified output file.

Args:
    INPUT (str):               path to the input VCF
    OUTPUT (str):              path to the output tab-separated file
    DEPTH_COLUMN (str):        the name of the FORMAT field that specifies total read coverage of
                               the variant locus
    VARIANT_DEPTH_FIELD (str): the name of the FORMAT field that specifies variant read coverage of
                               the variant locus
"""

import argparse

## Import sys modules
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
        type = str)
    parser.add_argument(
        "-o",
        "--output", 
        help = "Path to the output tab-separated file. If it does not exist, it will be created.",
        type = str)
    parser.add_argument(
        "-d",
        "--total_depth",
        help = "Name of the FORMAT field that describes the total read depth of the variant locus.",
        type = str)
    parser.add_argument(
        "-v",
        "--variant_depth",
        help = "Name of the FORMAT field that describes the variant read depth of the locus.",
        type = str)
    args = parser.parse_args()
    if not args.input:
        raise ValueError("You have not specified an input")
    elif not args.output:
        raise ValueError("You have not specified an output")
    elif not args.total_depth:
        raise ValueError("You have not specified the total depth field")
    elif not args.variant_depth:
        raise ValueError("You have not specified the variant read depth field")
    annotator.main.run_annotator(
        args.input,
        args.output,
        args.total_depth,
        args.variant_depth
    )

if __name__ == "__main__":
    main()
