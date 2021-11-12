#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 09:41:56 2018

@author: beccalove
"""

import os
import subprocess
import sys
import unittest

import allel
import compkaryo
import numpy as np

WORKING_DIR = os.path.dirname(os.path.abspath(__file__))

class TestImportData(unittest.TestCase):

    def setUp(self):

        path = os.path.join(WORKING_DIR, "test.vcf")

        self.callset = compkaryo.import_data(path)

    def test_samples(self):

        samples = np.array(["S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8",
                            "S9", "S10"])

        np.testing.assert_array_equal(samples, self.callset["samples"])

    def test_length(self):

        self.assertEqual(len(self.callset["variants/POS"]), 103)

    def test_fails(self):

        bad_path = "foo"

        self.assertRaises(FileNotFoundError, compkaryo.import_data, bad_path)

class TestImportInversion(unittest.TestCase):

    def setUp(self):

        inversion = "test"

        self.targets = compkaryo.import_inversion(inversion)

    def test_targets_are_correct(self):

        targets = np.array([17, 25, 33, 39, 52, 69, 112, 125, 141, 166,
                            175, 182, 191, 207, 225, 232, 242, 255, 267, 278,
                            299, 304, 310, 318, 322, 337, 345, 352, 380, 393,
                            401, 434, 445, 456, 476, 485, 499, 515, 525, 534,
                            543, 555, 570, 585, 592, 613, 627, 639, 642, 678])

        np.testing.assert_array_equal(targets, self.targets)

class TestExtractVtblIndices(unittest.TestCase):

    def setUp(self):

        ##need to mock up a callset and targets
        self.callset = allel.VariantTable(
                {"variants/POS" : np.array([5, 10, 15, 20, 25, 30, 35, 40]),
                 "variants/CHROM" : np.array(["foo", "foo", "foo", "foo",
                                              "foo", "foo", "foo", "foo"])})

        self.targets = np.array([10, 15, 35, 40])

    def test_basic_extraction(self):

        indices = compkaryo.extract_vtbl_indices(self.targets, self.callset,
                                                 chrom="foo")

        np.testing.assert_array_equal(np.array([1, 2, 6, 7]), indices)

    def test_genotypes_out_of_order(self):

        callset = allel.VariantTable(
                {"variants/POS" : np.array([10, 5, 20, 15, 30, 40, 35, 25]),
                 "variants/CHROM" : np.array(["foo", "foo", "foo", "foo",
                                              "foo", "foo", "foo", "foo"])})

        indices = compkaryo.extract_vtbl_indices(self.targets, callset,
                                                 chrom="foo")

        np.testing.assert_array_equal(np.array([0, 3, 6, 5]), indices)

class TestGenotypeAtConcordantSites(unittest.TestCase):

    def setUp(self):

        gt = np.array([[[0, 1], [-1, -1], [0, 0], [1, 1]],
                       [[0, 0], [0, 0], [0, 0], [1, 1]],
                       [[0, 1], [0, 0], [0, 0], [0, 1]],
                       [[1, 1], [0, 0], [0, 0], [0, 0]],
                       [[-1, -1], [0, 1], [0, 0], [1, 1]],
                       [[0, 1], [0, 0], [0, 1], [1, 1]],
                       [[0, 1], [0, 0], [0, 0], [1, 1]],
                       [[0, 1], [-1, -1], [0, 1], [1, 1]],
                       [[0, 1], [0, 0], [0, 0], [1, 1]],
                       [[0, 1], [0, 0], [0, 1], [1, 1]]])

        alt = np.array(["A", "T", "C", "G", "T", "A", "C", "G", "A", "T"])

        ref = np.array(["T", "C", "G", "T", "A", "C", "G", "A", "T", "A"])

        chrom = np.array(["2", "2", "2", "2", "2", "2", "2", "2", "2", "2"])

        filter_pass = np.array([False] * 10)

        pos = np.array([5, 14, 16, 28, 45, 67, 81, 105, 148, 177])

        self.callset = {"samples" : np.array(["bar1", "bar2", "bar6", "bar23"]),
                        "calldata/GT" : gt,
                        "variants/ALT" : alt,
                        "variants/CHROM" : chrom,
                        "variants/FILTER_PASS" : filter_pass,
                        "variants/POS" : pos,
                        "variants/REF" : ref}

    def test_without_samples_bool(self):

        indices = [0, 1, 2, 3, 4, 5, 6]

        results =\
        compkaryo.calculate_genotype_at_concordant_sites(self.callset,
                                                         indices)

        test_results = [1, (1/6), (1/7), (11/7)]

        total_called = [6, 6, 7, 7]

        np.testing.assert_array_equal(results,
                                      (test_results, total_called))

    def test_with_samples_bool(self):

        indices = [4, 5, 6, 7, 8, 9]

        samples_bool = [True, True, False, False]

        results =\
        compkaryo.calculate_genotype_at_concordant_sites(self.callset,
                                                         indices,
                                                         samples_bool)

        test_results = [1, 0.2]

        total_called = [5, 5]

        np.testing.assert_array_equal(results, (test_results, total_called))

    def test_with_totals(self):

        indices = [0, 2, 4, 6, 8]

        results =\
        compkaryo.calculate_genotype_at_concordant_sites(self.callset,
                                                         indices, totals=True)

        test_results = [1, 1/4, 0, 9/5]

        total_called = [4, 4, 5, 5,]

        num_0 = [0, 3, 5, 0]

        num_1 = [4, 1, 0, 1]

        num_2 = [0, 0, 0, 4]

        np.testing.assert_array_equal(results, (test_results, total_called,
                                                num_0, num_1, num_2))

class IntegrationTest(unittest.TestCase):

    def setUp(self):

        self.pkg_dir =\
        os.path.dirname(sys.modules["compkaryo"].__file__).replace(" ", "\ ")  

    def test_basic_output(self):

        command = (sys.executable +
                   " " + self.pkg_dir +
                   "/compkaryo.py" + " " + WORKING_DIR.replace(" ", "\ ") +\
                   "/test.vcf test")

        output = subprocess.getoutput(command)

        clean = np.array([[float(item) for item in line.split("\t")]\
                          for line in output.split("\n")])

        truth_set = np.array([[75/42, 42],
                              [7/43, 43],
                              [43/43, 43],
                              [41/43, 43],
                              [5/42, 42],
                              [75/41, 41],
                              [3/42, 42],
                              [72/39, 39],
                              [42/43, 43],
                              [77/41, 41]])

        np.testing.assert_array_equal(clean, truth_set)

    def test_with_samples_bool_file(self):

        command = (sys.executable +
                   " " + self.pkg_dir +
                   "/compkaryo.py" + " " + WORKING_DIR.replace(" ", "\ ") +\
                   "/test.vcf test -s " + WORKING_DIR.replace(" ", "\ ") +\
                   "/test_samples_bool.txt")

        output = subprocess.getoutput(command)

        clean = np.array([[float(item) for item in line.split("\t")]\
                          for line in output.split("\n")])

        truth_set = np.array([[75/41, 41],
                              [3/42, 42],
                              [72/39, 39],
                              [77/41, 41]])

        np.testing.assert_array_equal(clean, truth_set)

    def test_with_samples_bool_command_line(self):

        command = (sys.executable +\
                   " " + self.pkg_dir +\
                   "/compkaryo.py" + " " + WORKING_DIR.replace(" ", "\ ") +\
                   "/test.vcf test -s S6 S7 S8 S10")

        output = subprocess.getoutput(command)

        clean = np.array([[float(item) for item in line.split("\t")]\
                          for line in output.split("\n")])

        truth_set = np.array([[75/41, 41],
                              [3/42, 42],
                              [72/39, 39],
                              [77/41, 41]])

        np.testing.assert_array_equal(clean, truth_set)

    def test_with_totals(self):

        command = (sys.executable +\
                   " " + self.pkg_dir +\
                   "/compkaryo.py" + " " + WORKING_DIR.replace(" ", "\ ") +\
                   "/test.vcf test -t")

        output = subprocess.getoutput(command)

        clean = np.array([[float(item) for item in line.split("\t")]\
                          for line in output.split("\n")])

        truth_set = np.array([[75/42, 42, 2, 5, 35],
                              [7/43, 43, 37, 5, 1],
                              [43/43, 43, 3, 37, 3],
                              [41/43, 43, 5, 35, 3],
                              [5/42, 42, 38, 3, 1],
                              [75/41, 41, 2, 3, 36],
                              [3/42, 42, 39, 3, 0],
                              [72/39, 39, 1, 4, 34],
                              [42/43, 43, 2, 40, 1],
                              [77/41, 41, 1, 3, 37]])

        np.testing.assert_array_equal(clean, truth_set)
