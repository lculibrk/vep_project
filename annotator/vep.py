""" 
This module handles VEP API calls. It implements functions for querying the Ensembl VEP API with
variants processed by the vcf module.
"""

import requests
import annotator.exceptions

def get_variant_annotations(vcf_lines, by_gene = False) -> list:
    """
    Takes in a list of parsed variant lines in the format returned by vcf.parse_vcf(), performs a
    VEP API call and returns the resulting JSON for each variant. 

    Args:
        vcf_lines (list):  a list of variants, as output by vcf.parse_vcf(), formatted: 
                           [CHROM, POS, REF, ALT, FILTER, N_VARIANT_READS, TOTAL_READS].
                           We explicitly trust this object to be correct by being created by
                           vcf.parse_vcf(). Use other data sources at your own risk.
        by_gene (boolean): whether to query for impact by transcript (default, False) or the
                           highest impact per gene (True). Default: False

    Returns:
        A list of API returns, one per input element. 
    """
    ## Error handling of empty VCFs is not this function's job, so we return an empty list for empty
    ## variants
    if len(vcf_lines) == 0:
        return []
    ## hardcoded REST API URL/headers
    url = "https://rest.ensembl.org"
    ext = "/vep/homo_sapiens/region"
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
    ## Build the query string. Per VEP documentation, format is same as VCF 4.0: CHROM POS ID REF
    ## ALT QUAL FILTER INFO.
    ## Query string is JSON formatted and space-separated. Sub a period in for ID, QUAL, FILTER,
    ## INFO. Concatenate the resulting list elements with a ", ", and then add the head and tail
    ## of the json with string addition.
    query_string = '", "'.join(
        [
            f"{line[0]} {line[1]} . {line[2]} {line[3]} . . ." 
            for line in vcf_lines
        ]
    )
    query_string = '{ "variants" : ["' + query_string + '" ] }'

    if by_gene:
        per_gene = "true"
    else:
        per_gene = "false"

    ## Parameters for the API call
    params = {"variant_class": "true", "hgvs": "true", "per_gene":per_gene}

    ## Post the request
    r = requests.post(url+ext,
                      headers = headers,
                      data = query_string,
                      params = params,
                      timeout = 60
                      )

    ## Throw an error if the request is an error
    if r.status_code != 200:
        raise annotator.exceptions.RequestError(
            (
                f"API request failed with status {r.status_code}!"
                f"See below for the data given to the API:\n{query_string}"
            )
            )
    return r.json()

def get_chunked_annotations(variant_lines, chunk_size = 200, by_gene = False):
    """
    Runs the get_variant_annotations() function on sublists of variant_lines, 
    chunk_size elements at a time.
    Ostensibly this exists because the VEP API has a limit of 200 items per request. 

    Args:
        variant_lines (list): list of variant lines (lists), as returned from vcf.parse_vcf()
        chunk_size (int):     number of elements in each batch API request
        by_gene (boolean):    whether to query for impact by transcript (default, False) or the 
                              highest impact per gene (True). Default: False
    
    Returns:
        A list of API returns, one per input element.
    """
    ## Verify that the provided chunk size is an integer and error if not
    if not isinstance(chunk_size, int):
        raise ValueError("Provided chunk size must be an integer > 0")
    ## Verify that the chunk size is a positive integer
    if chunk_size < 1:
        raise ValueError("Provided chunk size must be an integer > 0")

    ## Run get_variant_annotations on chunk_size length sublists of the variants
    ret = [
            get_variant_annotations(variant_lines[sublist:sublist+chunk_size], by_gene = by_gene)
            for sublist in range(0, len(variant_lines), chunk_size)
          ]

    ## Flatten the resulting list of chunk-size lists
    ret = [element for sublist in ret for element in sublist]
    return ret

def merge_variant_annotation(variant, annotation, consequences_list) -> list:
    """
    Parse a single variant and its associated annotations, extracting annotations according to a 
    list of possible consequences keys.
    This function is intended to be used ie. by looping on a zip of the full list of 
    variants/annotations.

    Args:
        variant (list):          a list corresponding to a single variant from the list returned 
                                 from vcf.parse_vcf()
        annotation (dict):       a dict of VEP annotation, corresponding to 
                                 get_variant_annotations([variant])
        consequence_list (list): a list of strings corresponding to the consequences dicts in VEP 
                                 annotations. 
    
    Returns:
        a list of variant annotations for the provided variant. One element per annotated 
        gene/transcript. 
    """
    ## Compute allele frequency in the variant calls
    var_af = str(round(int(variant[5])/int(variant[6]), ndigits = 5))

    ## Extract variant class
    variant_class = annotation.get("variant_class", "NA")

    ## Get each consequence entry if it exists
    features = [annotation.get(consequence) for consequence in consequences_list]
    features = [feature for feature in features if feature is not None]

    ## If we have transcript_consequences or intergenic_consequences, loop through and extract
    ## gene_id, gene_symbol, transcript_id, impact, and consequence_terms
    if features:
        ## Flatten the list
        features = [element for sublist in features for element in sublist]

        ## Extract the desired annotations, returning "NA" if they do not exist
        extracted = [
            [
                feature.get("gene_id", "NA"),
                feature.get("gene_symbol", "NA"),
                feature.get("transcript_id", "NA"),
                feature.get("impact", "NA"),
                feature.get("hgvsp", "NA"),
                ";".join(feature.get("consequence_terms", "NA"))
            ] for feature in features
        ]
    else:
        ## If there's nothing, then we simply return NAs
        extracted = ["NA"] * 6

    ## Allele frequency data is contained within colocated_variants.
    ## colocated_variants is a list of dicts, as multiple variant databases are used
    ## (ie. HGMD, COSMIC, dbSNP).
    ## if population allele frequencies are available, they are a key named "frequencies"
    colocated_data = annotation.get("colocated_variants")

    ## If there's annotated variants at this position (colocated variants), pull COSMIC
    ## IDs and minor allele frequencies
    if colocated_data:
        ## Check for COSMIC IDs. VEP doesn't match COSMIC ID with the actual allele of interest,
        ## so in this case we return the list of all COSMIC variants at this locus.
        cosmic_strings = [
                            element.get("id")
                            for element in colocated_data if
                                element.get("allele_string") == "COSMIC_MUTATION"
                        ]
        if cosmic_strings:
            ## Concatenate with semicolons for a single output field
            cosmic_strings = ";".join(cosmic_strings)
        else:
            ## If there's no COSMIC IDs, then set it to None
            cosmic_strings = "NA"

        ## Check for dbSNP (rs) IDs. There should only be one dbSNP ID per variant.
        ## Need to do two loops here, first to get all the "id" values, second to select those
        ## starting with "rs"
        rs_strings = [element.get("id") for element in colocated_data if element.get("id")]
        rs_strings = [rs for rs in rs_strings if rs.startswith("rs")]

        ## Extract rs IDs
        if rs_strings:
            ## Concatenate with semicolons for a single output field
            rs_string = ";".join(rs_strings)
        else:
            ## If no dbSNP IDs, set to None
            rs_string = "NA"

        ## Check for allele frequency data. Frequency is a dict, each key is an alt allele, which
        ## is another dict of allele frequencies from various sources, such as gnomAD
        af_data = [element for element in colocated_data if element.get("frequencies")]

        if af_data:
            ## Get the frequencies dict for the variant matching our alt allele
            ## We use the first instance of frequencies because population allele frequencies
            ## should be the same for multiple colocated_data entries
            af_dict = af_data[0].get("frequencies", {}).get(variant[3], {})
            ## Get the allele frequency. We first try the "af" key, and if it isn't found, we fall
            ## back to "gnomade", finally returning "NA" if neither
            af = str(af_dict.get("af", af_dict.get("gnomade", "NA")))
        else:
            af = "NA"
    else: ## There's no annotated variant here
        af = "NA"
        cosmic_strings = "NA"
        rs_string = "NA"

    result = [
        variant + [var_af, variant_class] + e + [rs_string, cosmic_strings, af] for e in extracted
        ]

    return result

def annotate_variants(vcf_lines, chunk_size = 200, by_gene = False) -> list:
    """
    Annotates variants and outputs a list of lists containing variant information, variant
    annotation, and variant impact, one list per genetic feature that each variant impacts.

    Args:
        vcf_lines (list): a list of variants, as output by vcf.parse_vcf(), formatted: 
                          [CHROM, POS, REF, ALT, FILTER, N_VARIANT_READS, TOTAL_READS].
                          We explicitly trust this object to be correct by being created by
                          vcf.parse_vcf(). Use other data sources at your own risk.
        chunk_size (int): number of elements in each batch API request
    
    Returns:
        A list of lists, each in this format, with "NA" for missing values: 
        [CHROM, POS, REF, ALT, FILTER, N_VARIANT_READS, TOTAL_READS, VARIANT_READ_FRACTION, 
        VARIANT_TYPE, GENE_ID, GENE_SYMBOL, TRANSCRIPT_ID, IMPACT, CONSEQUENCE_TERMS, HGVSP, 
        DBSNP_ID, COSMIC_ID, POP_AF]
    """
    ## Get annotations for the variants
    annotations = get_chunked_annotations(vcf_lines, chunk_size = chunk_size, by_gene = by_gene)

    ## Set which VEP consequences we want
    consequences_list = ["transcript_consequences", "intergenic_consequences"]

    ## Use a list comprehension to run merge_variant_annotation on every variant and its associated
    ## VEP annotations
    out_list =  [
                    annotator.vep.merge_variant_annotation(variant, annotation, consequences_list)
                    for variant, annotation in zip(vcf_lines, annotations)
                ]

    ## Flatten the out_list
    out_list = [element for sublist in out_list for element in sublist]

    return out_list
