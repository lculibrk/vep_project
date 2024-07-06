import os
import unittest

from annotator import vcf, vep

class test_annotate_vcf(unittest.TestCase):
    def test_annotate_platypus_vcf(self):
        """
        Test annotating a platypus VCF
        """
        vcf_path = os.path.join(os.path.dirname(__file__), "data", "platypus.vcf")
        vcf_file = vcf.read_vcf(vcf_path)
        vcf_parsed = vcf.parse_vcf(vcf_file, "NR", "NV")
        annotations = vep.annotate_variants(vcf_parsed)
        self.assertEqual(len(vcf_parsed), len(annotations))

    def test_annotate_empty_vcf(self):
        """
        Test annotating a VCF with no variants
        """
        vcf_path = os.path.join(os.path.dirname(__file__), "data", "empty.vcf")
        vcf_file = vcf.read_vcf(vcf_path)
        vcf_parsed = vcf.parse_vcf(vcf_file, "NR", "NV")
        annotations = vep.annotate_variants(vcf_parsed)
        self.assertEqual(len(annotations), 0)
    
    def test_chunked_vep_size_1(self):
        """
        Test submitting 2 variants, one at a time to the VEP API
        """
        vcf_path = os.path.join(os.path.dirname(__file__), "data", "platypus.vcf")
        vcf_file = vcf.read_vcf(vcf_path)
        vcf_parsed = vcf.parse_vcf(vcf_file, "NR", "NV")
        annotations = vep.annotate_chunks(vcf_parsed[0:2], 1)
        self.assertEqual(len(annotations), 2)
    
    def test_chunked_vep_same_size(self):
        """
        Test submitting data whose size is == chunk size
        """
        vcf_path = os.path.join(os.path.dirname(__file__), "data", "platypus.vcf")
        vcf_file = vcf.read_vcf(vcf_path)
        vcf_parsed = vcf.parse_vcf(vcf_file, "NR", "NV")
        annotations = vep.annotate_chunks(vcf_parsed[0:2], 2)
        self.assertEqual(len(annotations), 2)

    def test_chunked_vep_gr_size(self):
        """
        Test submitting data whose size is > chunk size
        """
        vcf_path = os.path.join(os.path.dirname(__file__), "data", "platypus.vcf")
        vcf_file = vcf.read_vcf(vcf_path)
        vcf_parsed = vcf.parse_vcf(vcf_file, "NR", "NV")
        annotations = vep.annotate_chunks(vcf_parsed[0:2], 3)
        self.assertEqual(len(annotations), 2)