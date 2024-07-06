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

