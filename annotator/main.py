import annotator.vcf
import annotator.vep

def annotator(vcf_path, total_cov_field, var_cov_field, sample_name = None):
    """
    Reads and parses the input VCF, makes an API call to Ensembl VEP, annotates variants and returns the resulting JSON
    
    Args:
        vcf_path (str):        path to the input VCF
        output (str):          path to the output tab-separated variants file
        total_cov_field (str): the name of the FORMAT field that contains TOTAL coverage.
        var_cov_field (str):   the name of the FORMAT field that contains VARIANT (non-reference) coverage.
        sample_name (str):     the name of the sample to be processed in the VCF. If not specified, the program will take the first column after FORMAT. Default: None

    Returns:
        a dictionary 
    """
    ## Call the VCF reading function
    vcf = annotator.vcf.read_vcf(vcf_path)

    ## Call the parsing function
    vcf_lines = annotator.vcf.parse_vcf(vcf, total_cov_field, var_cov_field, sample_name = sample_name)
    
    ## Annotate variants
    annotations = annotator.vep.annotate_variants(vcf_lines).json()
    return(annotations)
