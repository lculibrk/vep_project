import annotator.exceptions

def read_vcf(vcf_path) -> list:
    """
    Read a VCF specified by the path to a list and strip the header. Enforce format compliance. 
    """
    ## Read the VCF, skipping header lines. Strip whitespace and split by tab.
    with open(vcf_path, "r") as f:
        vcf = [l for l in f.readlines() if not l.startswith("##")]
    vcf = [l.strip().split("\t") for l in vcf]

    ## Throw a useful error if VCF is completely empty with no header line
    if len(vcf) == 0:
        raise ValueError(f"File {vcf_path} has no VCF header line (beginning with #CHROM)")
    
    ## Enforce that all the VCF ver.4 fixed fields are present in line 1
    vcf_fields = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    for field in vcf_fields:
        if field not in vcf[0]:
            raise ValueError(f"Field {field} not found in the VCF! Are you sure this is a format compliant VCF?")
    ## If the file is format-compliant, but empty, it should run without errors. So we won't check for length at this stage.
    return(vcf)

def parse_vcf(vcf, total_cov_field, var_cov_field, sample_name = None) -> list:
    """
    Takes in a list of VCF lines, parses and outputs as a list of lists, one per line of data

    Args:
        vcf (list):            list of VCF lines, split by tab. 
        total_cov_field (str): the name of the FORMAT field that contains TOTAL coverage.
        var_cov_field (str):   the name of the FORMAT field that contains VARIANT (non-reference) coverage.
        sample_name (str):     the name of the sample to be processed in the VCF. If not specified, the program will take the first column after FORMAT. Default: None

    Returns:
        A list of lists, in this format:
        [CHROM, POS, REF, ALT, FILTER, N_VARIANT_READS, TOTAL_READS]
    """
    vcf_colnames = vcf[0]
    vcf = vcf[1:]

    ## Extract indices of the fixed VCF fields
    chrom_ind = vcf_colnames.index("#CHROM")
    pos_ind = vcf_colnames.index("POS")
    ref_ind = vcf_colnames.index("REF")
    alt_ind = vcf_colnames.index("ALT")
    filt_ind = vcf_colnames.index("FILTER")
    fmt_ind = vcf_colnames.index("FORMAT")

    ## If the sample name is specified, then check that it exists and use it. Otherwise, use the first column after FORMAT
    if sample_name:
        if sample_name not in vcf_colnames:
            raise ValueError(f"Sample {sample_name} was specified, but not found in the VCF! VCF columns: {vcf_colnames}")
        else:
            samp_ind = vcf_colnames.index(sample_name)
    else:
        samp_ind = fmt_ind + 1
        ## Throw an error if the field after FORMAT doesn't exist
        if len(vcf_colnames) < samp_ind:
            raise annotator.exceptions.MalformedDataError("Input VCF does not contain any genotype fields!")
    
    ## Split the genotype and format columns by colons
    genotype_column = [l[samp_ind].split(":") for l in vcf]
    fmt_column = [l[fmt_ind].split(":") for l in vcf]

    ## Zip the vcf, genotype, and format lists. Get the required fields from vcf by indexing using the prior extracted indices.
    ## We get the variant/total coverage fields by getting the matching index for var/total_cov_field and indexing the genotype 
    ## fields by that index. 
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
    return(tbl)