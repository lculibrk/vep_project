import requests
import annotator.exceptions

def annotate_variants(variant_lines) -> list:
    """
    Takes in a list of parsed variant lines in the format returned by vcf.parse_vcf(), performs a VEP API call and returns the resulting JSON for each variant. 

    Args:
        variant_lines (list):  list of variant lines (lists), as returned from vcf.parse_vcf() 
        total_cov_field (str): the name of the FORMAT field that contains TOTAL coverage.
        var_cov_field (str):   the name of the FORMAT field that contains VARIANT (non-reference) coverage.
        sample_name (str):     the name of the sample to be processed in the VCF. If not specified, the program will take the first column after FORMAT. Default: None

    Returns:
        A list of API returns, one per input element. 
    """
    ## Error handling of empty VCFs is not this function's job, so we return an empty list for empty variants
    if len(variant_lines) == 0:
        return([])
    ## hardcoded REST API URL/headers
    url = "https://rest.ensembl.org"
    ext = "/vep/homo_sapiens/region"
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
    ## Build the query string. Per VEP documentation, format is same as VCF 4.0: CHROM POS ID REF ALT QUAL FILTER INFO. 
    ## Query string is JSON formatted and space-separated. Sub a period in for ID, QUAL, FILTER, INFO. Concatenate the
    ## resulting list elements with a ", ", and then add the head and tail of the json with string addition.
    query_string = '", "'.join([f"{line[0]} {line[1]} . {line[2]} {line[3]} . . ." for line in variant_lines])
    query_string = '{ "variants" : ["' + query_string + '" ] }'

    ## Post the request
    r = requests.post(url+ext, headers=headers, data=query_string)
    ## Throw an error if the request is an error
    if r.status_code != 200:
        raise annotator.exceptions.RequestError(f"API request failed with status {r.status_code}! See below for the data given to the API:\n{query_string}")
    return(r.json())

def annotate_chunks(variant_lines, chunk_size = 200):
    """
    Runs the annotate_variants() function on sublists of variant_lines, chunk_size elements at a time.
    Ostensibly this exists because the VEP API has a limit of 200 items per request. 

    Args:
        variant_lines (list):   list of variant lines (lists), as returned from vcf.parse_vcf()
        chunk_size (int):       number of elements in each batch API request
    
    Returns:
        A list of API returns, one per input element.
    """
    ## Verify that the provided chunk size is an integer and error if not
    if not isinstance(chunk_size, int):
        raise ValueError("Provided chunk size must be an integer > 0")
    ## Verify that the chunk size is a positive integer
    if chunk_size < 1:
        raise ValueError("Provided chunk size must be an integer > 0")
    
    ## Run annotate_variants on chunk_size length sublists of the variants
    ret = [annotate_variants(variant_lines[sublist:sublist+chunk_size]) for sublist in range(0, len(variant_lines), chunk_size)]

    ## Flatten the resulting list of chunk-size lists
    ret = [element for sublist in ret for element in sublist]
    return(ret)

def parse_vep(vep_return) -> list:
    """
    Parse the output of annotate_chunks/annotate_variants
    """
    return(None)