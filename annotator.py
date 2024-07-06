#!/usr/bin/env python

## Import sys modules
import sys
import annotator.vcf
import annotator.vep
import pprint


def main():
    vcf = annotator.vcf.read_vcf("data/test_vcf_data.txt")
    vcf_lines = annotator.vcf.parse_vcf(vcf, "NR", "NV")
    pprint.pprint(annotator.vep.annotate_chunks(vcf_lines[:2], 2))
    #pprint.pprint(annotator.vep.annotate_chunks(vcf_lines[0:10]).json(), 1)
    return()

if __name__ == "__main__":
    main()