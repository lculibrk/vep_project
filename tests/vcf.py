""" This module implements testing for the vcf module. """
import os
import unittest
from unittest.mock import mock_open, patch

import annotator.exceptions
from annotator import vcf


class test_load_vcf(unittest.TestCase):
    """ Unit tests for the load_vcf function """
    def test_load_platypus_vcf(self):
        """
        Test loading a platypus-formatted VCF
        """
        vcf_path = os.path.join(os.path.dirname(__file__), "data", "platypus.vcf")
        vcf_file = vcf.read_vcf(vcf_path)
        self.assertEqual(len(vcf_file), 18)

    def test_load_empty_vcf(self):
        """
        Test loading a VCF with no variants
        """
        vcf_path = os.path.join(os.path.dirname(__file__), "data", "empty.vcf")
        vcf_file = vcf.read_vcf(vcf_path)
        self.assertEqual(len(vcf_file), 1)

    def test_error_malformed_vcf(self):
        """
        Test error handling for VCFs without correct column headers 
        """
        m = mock_open(read_data = '#CHROM\tBLAHBLAH\tI_LAUGH\tAT_FILE_FORMAT\tCONFORMITY\n')
        with patch('builtins.open', m):
            with self.assertRaises(annotator.exceptions.MalformedDataError):
                vcf.read_vcf("dummy_path.vcf")

    def test_error_empty_vcf(self):
        """
        Test error handling when loading an empty file
        """
        m = mock_open(read_data = '')
        with patch('builtins.open', m):
            with self.assertRaises(annotator.exceptions.MalformedDataError):
                vcf.read_vcf("dummy_path.vcf")

class test_parse_vcf(unittest.TestCase):
    """ Unit tests for the parse_vcf function """
    def test_parse_platypus_vcf(self):
        """
        Test parsing a platypus-formatted VCF and verify that it has a consistent number of columns.
        """
        vcf_path = os.path.join(os.path.dirname(__file__), "data", "platypus.vcf")
        vcf_file = vcf.read_vcf(vcf_path)
        vcf_parsed = vcf.parse_vcf(vcf_file, "NR", "NV")
        it = iter(vcf_parsed)
        the_len = len(next(it))
        self.assertTrue(all(len(l) == the_len for l in it))

    def test_parse_empty_vcf(self):
        """
        Test parsing an empty VCF and verify that it produces an empty output
        """
        vcf_path = os.path.join(os.path.dirname(__file__), "data", "empty.vcf")
        vcf_file = vcf.read_vcf(vcf_path)
        vcf_parsed = vcf.parse_vcf(vcf_file,  "NR", "NV")
        self.assertEqual(len(vcf_parsed), 0)

    def test_multiallelic_vcf(self):
        """
        Test splitting multiallelic sites
        """
        vcf_lines = [
            [
                "#CHROM" ,
                "POS" ,
                "ID" ,
                "REF" ,
                "ALT" ,
                "QUAL" ,
                "FILTER",
                "INFO" ,
                "FORMAT",
                "sample"
            ],
            [
                "3", 
                "64527465", 
                ".", 
                "C" ,
                "A,T" ,
                "2906" ,
                "PASS" ,
                "",
                "GT:GL:GOF:GQ:NR:NV" ,
                "1/2:-1,-1,-1:8:99:169,169:75,92"
            ]
        ]
        expected = [['3', '64527465', 'C', 'A', 'PASS', '75', '169'],
                    ['3', '64527465', 'C', 'T', 'PASS', '92', '169']]
        vcf_parsed = vcf.parse_vcf(vcf_lines, "NR", "NV")
        self.assertEqual(vcf_parsed, expected)
