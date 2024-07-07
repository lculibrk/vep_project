""" This module contains the main runner function, run_annotator"""
import os

import annotator.vcf
import annotator.vep

def run_annotator(vcf_path, output, total_cov_field, var_cov_field, sample_name = None,
                  chunk_size = 200):
    """
    Reads and parses the input VCF, makes an API call to Ensembl VEP, annotates variants and 
    returns the resulting JSON
    
    Args:
        vcf_path (str):        path to the input VCF
        output (str):          path to the output tab-separated variants file
        total_cov_field (str): the name of the FORMAT field that contains TOTAL coverage.
        var_cov_field (str):   the name of the FORMAT field that contains VARIANT (non-reference) 
                               coverage.
        sample_name (str):     the name of the sample to be processed in the VCF. If not specified,
                               the program will take the first column after FORMAT. Default: None
        by_gene (boolean):     whether to query for impact by transcript (default, False) or the 
                               highest impact per gene (True). Default: False
        chunk_size (int):      the maximum number of variants to submit per VEP API call. 
                               Default: 200

    Returns:
        a list of lists, each corresponding to an variant:genetic feature combination. For example,
        a variant affecting five transcripts/genes would have five corresponding lines.
    """
    ## Call the VCF reading function
    vcf = annotator.vcf.read_vcf(vcf_path)

    ## Call the parsing function
    vcf_lines = annotator.vcf.parse_vcf(
        vcf,
        total_cov_field,
        var_cov_field,
        sample_name = sample_name
        )

    ## Annotate the VCF lines using the received annotations
    out_list = annotator.vep.annotate_variants(vcf_lines, chunk_size = chunk_size)

    ## Create output path if it does not exist
    os.makedirs(os.path.dirname(output), exist_ok = True)

    ## Write the tsv to output
    with open(output, "w", encoding = "UTF-8") as f:
        f.write(
            "CHROM\tPOS\tREF\tALT\tFILTER\tN_VARIANT_READS\tTOTAL_READS\tVARIANT_READ_FRACTION\t"
            "VARIANT_TYPE\tGENE_ID\tGENE_SYMBOL\tTRANSCRIPT_ID\tIMPACT\tCONSEQUENCE_TERMS\tHGVSP\t"
            "DBSNP_ID\tCOSMIC_ID\tPOP_AF\n"
            )
        for line in out_list:
            f.write("%s\n" % "\t".join(line))

    return out_list
