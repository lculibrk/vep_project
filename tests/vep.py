import os
import unittest

from annotator import vcf, vep, main

class test_get_variant_annotations(unittest.TestCase):
    def test_get_variant_annotations_platypus(self):
        """
        Test annotating a platypus VCF
        """
        vcf_path = os.path.join(os.path.dirname(__file__), "data", "platypus.vcf")
        vcf_file = vcf.read_vcf(vcf_path)
        vcf_parsed = vcf.parse_vcf(vcf_file, "NR", "NV")
        annotations = vep.get_variant_annotations(vcf_parsed)
        self.assertEqual(len(vcf_parsed), len(annotations))

    def test_get_variant_annotations_empty_vcf(self):
        """
        Test annotating a VCF with no variants
        """
        vcf_path = os.path.join(os.path.dirname(__file__), "data", "empty.vcf")
        vcf_file = vcf.read_vcf(vcf_path)
        vcf_parsed = vcf.parse_vcf(vcf_file, "NR", "NV")
        annotations = vep.get_variant_annotations(vcf_parsed)
        self.assertEqual(len(annotations), 0)
    
class test_chunked_annotations(unittest.TestCase):
    def test_chunked_vep_size_1(self):
        """
        Test submitting 2 variants, one at a time to the VEP API
        """
        vcf_path = os.path.join(os.path.dirname(__file__), "data", "platypus.vcf")
        vcf_file = vcf.read_vcf(vcf_path)
        vcf_parsed = vcf.parse_vcf(vcf_file, "NR", "NV")
        annotations = vep.get_chunked_annotations(vcf_parsed[0:2], 1)
        self.assertEqual(len(annotations), 2)
    
    def test_chunked_vep_same_size(self):
        """
        Test submitting data whose size is == chunk size
        """
        vcf_path = os.path.join(os.path.dirname(__file__), "data", "platypus.vcf")
        vcf_file = vcf.read_vcf(vcf_path)
        vcf_parsed = vcf.parse_vcf(vcf_file, "NR", "NV")
        annotations = vep.get_chunked_annotations(vcf_parsed[0:2], 2)
        self.assertEqual(len(annotations), 2)

    def test_chunked_vep_gr_size(self):
        """
        Test submitting data whose size is > chunk size
        """
        vcf_path = os.path.join(os.path.dirname(__file__), "data", "platypus.vcf")
        vcf_file = vcf.read_vcf(vcf_path)
        vcf_parsed = vcf.parse_vcf(vcf_file, "NR", "NV")
        annotations = vep.get_chunked_annotations(vcf_parsed[0:2], 3)
        self.assertEqual(len(annotations), 2)

    def test_chunked_vep_emptydata(self):
        """
        Test submitting empty data to get_chunked_annotations
        """
        vcf_path = os.path.join(os.path.dirname(__file__), "data", "empty.vcf")
        vcf_file = vcf.read_vcf(vcf_path)
        vcf_parsed = vcf.parse_vcf(vcf_file, "NR", "NV")
        annotations = vep.get_chunked_annotations(vcf_parsed)
        self.assertEqual(len(annotations), 0)

class test_merge_variant_annotation(unittest.TestCase):
    def test_intergenic_variant(self):
        """
        Test annotation of a variant that is intergenic. 
        """
        variant_list = [["9", "82929050", "A", "T", "PASS", 100, 150]]
        expected_output = [['9', '82929050', 'A', 'T', 'PASS', 100, 150, '0.66667', 'SNV', 'NA', 'NA', 'NA', 'MODIFIER', 'NA', 'intergenic_variant', 'NA', 'NA', 'NA']]
        annotations = vep.annotate_variants(variant_list)
        self.assertEqual(len(annotations), 1)
        self.assertEqual(annotations, expected_output)


    def test_mixed_genic_intergenic_variants(self):
        """
        Test annotation of a variant that is intergenic, followed by a variant that is genic but unannotated.
        """
        variant_list = [["9", "82929050", "A", "T", "PASS", 100, 150],
                        ["9", "83004550", "T", "G", "PASS", 58, 92]]
        expected_output = [['9', '82929050', 'A', 'T', 'PASS', 100, 150, '0.66667', 'SNV', 'NA', 'NA', 'NA', 'MODIFIER', 'NA', 'intergenic_variant', 'NA', 'NA', 'NA'], 
                           ['9', '83004550', 'T', 'G', 'PASS', 58, 92, '0.63043', 'SNV', 'ENSG00000165105', 'RASEF', 'ENST00000376447', 'LOW', 'ENSP00000365630.3:p.Arg384=', 'synonymous_variant', 'NA', 'NA', 'NA']]
        annotations = vep.annotate_variants(variant_list)
        self.assertEqual(annotations, expected_output)

    def test_known_pathogenic_variants(self):
        """
        Test some known pathogenic cancer variants: (KRAS G12C, 12:25245351:G>T), (BRAF V600E, 7:140753336:A>T)
        Ensure that they return COSMIC IDs so they aren't missed in a real analysis
        """
        variant_list = [["12", "25245351", "G", "T", "PASS", 100, 150],
                        ["7", "140753336", "A", "T", "PASS", 58, 92]]
        annotations = vep.annotate_variants(variant_list)
        cosmic_ids = [line[16] for line in annotations]
        self.assertTrue(all([id.startswith("COS") for id in cosmic_ids]))

        



class test_annotate_variants(unittest.TestCase):
    def test_line_column_counts(self):
        """
        Test that parse_variants returns a constant number of columns
        """
        vcf_path = os.path.join(os.path.dirname(__file__), "data", "platypus.vcf")
        vcf_file = vcf.read_vcf(vcf_path)
        vcf_parsed = vcf.parse_vcf(vcf_file, "NR", "NV")
        out_list = vep.annotate_variants(vcf_parsed)
        it = iter(out_list)
        the_len = len(next(it))
        self.assertTrue(all(len(l) == the_len for l in it))

    def test_full_annotation_empty(self):
        """
        Test that an empty input VCF results in an empty annotated mutation table
        """
        vcf_path = os.path.join(os.path.dirname(__file__), "data", "empty.vcf")
        vcf_file = vcf.read_vcf(vcf_path)
        vcf_parsed = vcf.parse_vcf(vcf_file, "NR", "NV")
        annotations = vep.annotate_variants(vcf_parsed)
        self.assertEqual(len(annotations), 0)

    def test_by_gene_arg(self):
        """
        Test that the by_gene argument to annotate_variants correctly results in one gene per annotation
        """
        vcf_path = os.path.join(os.path.dirname(__file__), "data", "platypus.vcf")
        vcf_file = vcf.read_vcf(vcf_path)
        vcf_parsed = vcf.parse_vcf(vcf_file, "NR", "NV")
        annotations_gene = vep.annotate_variants(vcf_parsed, by_gene = True)
        ## Get the gene IDs and show that the number of gene IDs == the length of the table
        gene_ids = set([str(line[1]) + "_" + line[3] + "_" + line[9] for line in annotations_gene])
        self.assertEqual(len(gene_ids), len(annotations_gene))

class test_main(unittest.TestCase):
    def test_correct_rows(self):
