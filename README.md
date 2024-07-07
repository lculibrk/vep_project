# VEP annotator project

Website: https://github.com/lculibrk/vep_project

Author: [Luka Culibrk](https://github.com/lculibrk)

This program annotates variants using [Ensembl's](https://www.ensembl.org/) [Variant Effect Predictor](https://www.ensembl.org/info/docs/tools/vep/index.html) [REST API](https://rest.ensembl.org/#VEP). 
It implements a pip-installable package, `annotator` and a script, `annotator.py`, which functions as a wrapper for `annotator`. It takes input as a [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) file, 
and outputs a tab-separated file, described below. 

The tool has been tested using some example data, with the results available at [output/main/variants.tsv](https://github.com/lculibrk/vep_project/blob/main/output/main/variants.tsv)

## Installation and dependencies
The package requires python, version 3.6-3.11, and `requests`. Two convenient ways of satisfying these dependencies are provided below:

### Option 1: conda

This option requires [conda](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html) to be installed
```
conda env create -n env/ -f environment.yml
conda activate env/
pip install ./
```

### Option 2: docker

This option requires [docker](https://www.docker.com/) to be installed. 
Alternatively if you use a shared HPC that forbids docker, you can use [singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html)

```
docker pull lculibrk/vep_project
```

## How to run this tool
### Execution:
```
./annotator.py -i INPUT -o OUTPUT -d DEPTH_FIELD -v VARIANT_DEPTH_FIELD [-g | --per-gene | --no-per-gene]
```
### Parameters:

| Parameter name | Description |
|----------------|-------------|
| INPUT | Path to the input VCF |
| OUTPUT | Path to the desired output file |
| DEPTH_FIELD | the FORMAT field in your VCF that corresponds to the total depth of the locus |
| VARIANT_DEPTH_FIELD | the FORMAT field in your VCF that corresponds to the variant allele depth |
| per-gene | whether to annotate by transcript (default, --no-per-gene) or per gene (-g, --per-gene) |


### Example using provided data

The data at [output/main/variants.tsv](https://github.com/lculibrk/vep_project/blob/main/output/main/variants.tsv) was generated using the below command:

```
./annotator.py -i data/test_vcf_data.txt -o output/main/variants.tsv -d NR -v NV
```

## Description of the output

`annotator` will produce a tab-separated file with 18 columns, described as follows:
	
| Column | Type | Description |
|---|---|---|
| CHROM | string | The chromosome where the variant is located |
| POS | int | The [1-based position](https://www.biostars.org/p/84686/) of the variant |
| REF | string | The reference sequence (allele) at this position |
| ALT | string | The variant sequence (allele) at this position |
| FILTER | string | The FILTER field carried over from the VCF. PASS indicates the variant passed all filters. |
| N_VARIANT_READS | int | The number of variant-supporting reads aligned at the variant position. |
| TOTAL_READS | int | The number of total (variant + reference) reads aligned at the variant position |
| VARIANT_READ_FRACTION | float | The fraction of variant-supporting reads aligned at the variant position, rounded to five digits. |
| VARIANT_TYPE | string | The type of variant, e.g. SNV, substitution (for multi-nucleotide variants), insertion, deletion |
| GENE_ID | string | Ensembl gene ID at this locus |
| GENE_SYMBOL | string | HUGO gene ID at this locus |
| TRANSCRIPT_ID | string | Ensembl transcript ID at this locus |
| IMPACT | string | VEP-estimated variant impact |
| CONSEQUENCE_TERMS | string | The consequence terms provided by VEP for this variant/gene/transcript, separated by a semicolon ";" if multiple were found |
| HGVSP | string | Ensembl protein ID and amino acid change, if applicable |
| DBSNP_ID | string | [dbSNP](https://www.ncbi.nlm.nih.gov/snp/) (rs) IDs for the variant, separated by a semicolon ";" if multiple were found |
| COSMIC_ID | string | [COSMIC](https://cancer.sanger.ac.uk/cosmic) variant IDs for the variant position, separated by a semicolon ";" if multiple were found |
| POP_AF | float | minor (population) allele frequency for the variant |

