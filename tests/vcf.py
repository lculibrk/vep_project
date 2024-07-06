import os
import unittest

from annotator import vcf

class test_load_vcf(unittest.TestCase):
    def test_load_platypus_vcf(self):
        """
        Test loading a platypus-formatted VCF
        """
        vcf_path = os.path.join(os.path.dirname(__file__), "data", "platypus.vcf")
        vcf_file = vcf.read_vcf(vcf_path)
        self.assertEqual(len(vcf_file), 16)

    def test_load_empty_vcf(self):
        """
        Test loading a VCF with no variants
        """
        vcf_path = os.path.join(os.path.dirname(__file__), "data", "empty.vcf")
        vcf_file = vcf.read_vcf(vcf_path)
        self.assertEqual(len(vcf_file), 1)

class test_parse_vcf(unittest.TestCase):
    def test_parse_platypus_vcf(self):
        """
        Test parsing a platypus-formatted VCF
        """
        vcf_path = os.path.join(os.path.dirname(__file__), "data", "platypus.vcf")
        vcf_file = vcf.read_vcf(vcf_path)
        vcf_parsed = vcf.parse_vcf(vcf_file, "NR", "NV")
        self.assertEqual(len(vcf_parsed), 15)
    def test_parse_empty_vcf(self):
        """
        Test parsing an empty VCF
        """
        vcf_path = os.path.join(os.path.dirname(__file__), "data", "empty.vcf")
        vcf_file = vcf.read_vcf(vcf_path)
        vcf_parsed = vcf.parse_vcf(vcf_file,  "NR", "NV")
        self.assertEqual(len(vcf_parsed), 0)