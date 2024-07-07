"""
This module handles VCF loading and parsing. It implements functions for reading and parsing VCFs
"""
import annotator.exceptions

def read_vcf(vcf_path) -> list:
    """
    Read a VCF specified by the path to a list and strip the header. Enforce format compliance. 

    Args:
        vcf_path (str): path to input VCF.
    
    Returns:
        A list of lists, each with the contents of a VCF line. Line 1 contains column names.

    Excepts:
        MalformedDataError: If the VCF is completely empty
        MalformedDataError: If the VCF is missing fields in the header
        MalformedDataError: If the VCF has an inconsistent number of columns
    """
    ## Read the VCF, skipping header lines. Strip whitespace and split by tab.
    with open(vcf_path, "r", encoding = "utf-8") as f:
        vcf = [l for l in f.readlines() if not l.startswith("##")]
    vcf = [l.strip().split("\t") for l in vcf]

    ## Throw a useful error if VCF is completely empty with no header line
    if len(vcf) == 0:
        raise annotator.exceptions.MalformedDataError(
            f"File {vcf_path} is completely empty (it must have a header at least)"
            )

    ## Enforce that all the VCF ver.4 fixed fields are present in line 1
    vcf_fields = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    for field in vcf_fields:
        if field not in vcf[0]:
            raise annotator.exceptions.MalformedDataError(
                f"Field {field} not found in the VCF! Are you sure this is a format compliant VCF?"
                )
    ## If the file is format-compliant, but empty, it should run without errors. So we won't check
    ## for length at this stage.

    ## Verify that the VCF is not truncated.
    expected_length = len(vcf[0])
    for sublist in vcf:
        if len(sublist) != expected_length:
            raise annotator.exceptions.MalformedDataError(
                "One of your VCF rows has fewer columns than expected. Possible truncation."
                f"Problematic column:\n{sublist}"
                )
    return vcf

def parse_vcf(vcf, total_cov_field, var_cov_field, sample_name = None) -> list:
    """
    Takes in a list of VCF lines, parses and outputs as a list of lists, one per line of data.
    Multiallelic sites are split into one list per alt allele.

    Args:
        vcf (list):            list of VCF lines, split by tab. 
        total_cov_field (str): the name of the FORMAT field that contains TOTAL coverage.
        var_cov_field (str):   the name of the FORMAT field that contains VARIANT (non-reference) 
                               coverage.
        sample_name (str):     the name of the sample to be processed in the VCF. If not specified, 
                               the program will take the first column after FORMAT. Default: None

    Returns:
        A list of lists, in this format: 
        [CHROM, POS, REF, ALT, FILTER, N_VARIANT_READS, TOTAL_READS]
    
    Excepts:
        ValueError: If a sample that is not in the VCF is specified by sample_name
        MalformedDataError: If the VCF has no genotype (sample) column
    """
    ## Extract VCF column names
    vcf_colnames = vcf[0]
    vcf = vcf[1:]

    ## Extract indices of the fixed VCF fields
    chrom_ind = vcf_colnames.index("#CHROM")
    pos_ind = vcf_colnames.index("POS")
    ref_ind = vcf_colnames.index("REF")
    alt_ind = vcf_colnames.index("ALT")
    filt_ind = vcf_colnames.index("FILTER")
    fmt_ind = vcf_colnames.index("FORMAT")

    ## indices of the values in the output file
    out_chrom = 0
    out_pos = 1
    out_ref = 2
    out_alt = 3
    out_filt = 4
    out_var_cov = 5
    out_tot_cov = 6


    ## If the sample name is specified, then check that it exists and use it. Otherwise, use the
    ## first column after FORMAT
    if sample_name:
        if sample_name not in vcf_colnames:
            raise ValueError(
                (
                f"Sample {sample_name} was specified, but not found in the VCF! VCF columns:"
                f"{vcf_colnames}"
                )
            )
        samp_ind = vcf_colnames.index(sample_name)
    else:
        samp_ind = fmt_ind + 1
        ## Throw an error if the field after FORMAT doesn't exist
        if len(vcf_colnames) < samp_ind:
            raise annotator.exceptions.MalformedDataError(
                "Input VCF does not contain any genotype fields!"
                )

    ## Split the genotype and format columns by colons
    genotype_column = [l[samp_ind].split(":") for l in vcf]
    fmt_column = [l[fmt_ind].split(":") for l in vcf]

    ## Zip the vcf, genotype, and format lists. Get the required fields from vcf by indexing
    ## using the prior extracted indices.
    ## We get the variant/total coverage fields by getting the matching index for the
    ## var/total_cov_field and indexing the genotype fields by that index.
    tbl = [
        [
            v[chrom_ind],
            v[pos_ind],
            v[ref_ind],
            v[alt_ind],
            v[filt_ind],
            geno[fmt.index(var_cov_field)],
            geno[fmt.index(total_cov_field)]
        ]
            for (v, geno, fmt) in zip(vcf, genotype_column, fmt_column)
    ]
    ## Separate multiallelic variants from non-multiallelic variants
    tbl_nonmultiallelic = [element for element in tbl if "," not in element[out_alt]]
    tbl_multiallelic = [element for element in tbl if "," in element[out_alt]]

    ## Loop over multiallelic sites, split by comma and create a new sublist for each allele
    tbl_multiallelic = [
        [
            [
                element[out_chrom],
                element[out_pos],
                element[out_ref],
                allele,
                element[out_filt],
                var,
                total
            ]
            for allele, var, total in zip(element[out_alt].split(","),
                                          element[out_var_cov].split(","),
                                          element[out_tot_cov].split(","))
        ] for element in tbl_multiallelic
    ]

    ## Flatten the list of lists of lists
    tbl_multiallelic = [element for sublist in tbl_multiallelic for element in sublist]

    ## Concatenate the list back together and sort by chrom, pos
    tbl = tbl_nonmultiallelic + tbl_multiallelic
    tbl.sort(key = lambda x: (x[0], int(x[1])))
    return tbl
