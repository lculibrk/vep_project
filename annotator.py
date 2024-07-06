#!/usr/bin/env python

## Import sys modules
import sys
import annotator.vcf


def main():
    ## Begin by parsing VCF and returning a table of the desired columns
    vcf = annotator.vcf.read_vcf("data/test_vcf_data.txt")
    annotator.vcf.parse_vcf(vcf, "NR", "NV")
    return()

if __name__ == "__main__":
    main()