# This file is part of pyplink.
#
# The MIT License (MIT)
#
# Copyright (c) 2014 Louis-Philippe Lemieux Perreault
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


from __future__ import print_function

import os
import random
import shutil
import unittest
from tempfile import mkdtemp

try:
    from itertools import zip_longest as zip
except ImportError:
    from itertools import izip_longest as zip

from pkg_resources import resource_filename

import numpy as np
import pandas as pd

from six.moves import range

from ..pyplink import PyPlink


class TestPyPlink(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Creating a temporary directory
        cls.tmp_dir = mkdtemp(prefix="beelinetools_test_")

        # Getting the BED/BIM/FAM files
        cls.bed = resource_filename(
            __name__,
            os.path.join("data", "test_data.bed"),
        )
        cls.bim = resource_filename(
            __name__,
            os.path.join("data", "test_data.bim"),
        )
        cls.fam = resource_filename(
            __name__,
            os.path.join("data", "test_data.fam"),
        )

        # Getting the prefix of the files
        cls.prefix = os.path.splitext(cls.bed)[0]

        # The list of markers
        cls.markers = ["rs10399749", "rs2949420", "rs2949421", "rs2691310",
                       "rs4030303", "rs4030300", "rs3855952", "rs940550",
                       "rs13328714", "rs11490937"]

        # The genotypes
        cls.genotypes = [[0, 0, 1], [0, 1, 0], [-1, -1, -1], [-1, -1, 1],
                         [0, 0, 0], [0, 0, 0], [0, 1, 2], [0, 0, 0], [1, 0, 0],
                         [0, 1, 0]]
        cls.acgt_genotypes = [["CC", "CC", "GC"], ["TT", "CT", "TT"],
                              ["00", "00", "00"], ["00", "00", "AT"],
                              ["GG", "GG", "GG"], ["CC", "CC", "CC"],
                              ["AA", "GA", "GG"], ["TT", "TT", "TT"],
                              ["GC", "CC", "CC"], ["GG", "AG", "GG"]]

        # Reading the files
        cls.pedfile = PyPlink(cls.prefix)

    @classmethod
    def tearDownClass(cls):
        # Cleaning the temporary directory
        shutil.rmtree(cls.tmp_dir)

    def test_pyplink_object_integrity(self):
        """Checks the integrity of the PyPlink object."""
        # Checking the name of the BED file
        self.assertTrue(hasattr(self.pedfile, "bed_filename"))
        self.assertEqual(self.bed, self.pedfile.bed_filename)

        # Checking the BED object
        self.assertTrue(hasattr(self.pedfile, "_bed"))
        self.assertTrue(isinstance(self.pedfile._bed, np.ndarray))

        # Checking the name of the BIM file
        self.assertTrue(hasattr(self.pedfile, "bim_filename"))
        self.assertEqual(self.bim, self.pedfile.bim_filename)

        # Checking the BIM object
        self.assertTrue(hasattr(self.pedfile, "_bim"))
        self.assertTrue(isinstance(self.pedfile._bim, pd.DataFrame))

        # Checking the name of the FAM file
        self.assertTrue(hasattr(self.pedfile, "fam_filename"))
        self.assertEqual(self.fam, self.pedfile.fam_filename)

        # Checking the FAM object
        self.assertTrue(hasattr(self.pedfile, "_fam"))
        self.assertTrue(isinstance(self.pedfile._fam, pd.DataFrame))

    def test_pyplink_object_error(self):
        """Checks what happens when we play with the PyPlink object."""
        # Creating a new object
        data = PyPlink(self.prefix)

        # Changing the BIM to None
        ori = data._bim
        data._bim = None
        with self.assertRaises(RuntimeError) as cm:
            data._read_bed()
        self.assertEqual("no BIM or FAM file were read", str(cm.exception))
        data._bim = ori

        # Changing the FAM to None
        ori = data._fam
        data._fam = None
        with self.assertRaises(RuntimeError) as cm:
            data._read_bed()
        self.assertEqual("no BIM or FAM file were read", str(cm.exception))
        data._fam = ori

    def test_pyplink_bad_bed(self):
        """Checks what happens when we read a bad BED file."""
        # The new file prefix
        new_prefix = os.path.join(self.tmp_dir, "bad_data")

        # Copying the FAM file
        new_fam = new_prefix + ".fam"
        with open(new_fam, "w") as o_file, open(self.fam, "r") as i_file:
            o_file.write(i_file.read())

        # Copying the BIM file
        new_fam = new_prefix + ".bim"
        with open(new_fam, "w") as o_file, open(self.fam, "r") as i_file:
            o_file.write(i_file.read())

        # Creating a new BED file with invalid number of bytes
        new_bed = new_prefix + ".bed"
        with open(new_bed, "wb") as o_file:
            o_file.write(bytearray([108, 27, 1, 1, 2, 3, 4]))

        # This should raise an exception
        with self.assertRaises(ValueError) as cm:
            PyPlink(new_prefix)
        self.assertEqual("invalid number of entries: {}".format(new_bed),
                         str(cm.exception))

        # Creating a new BED file with invalid first byte
        new_bed = new_prefix + ".bed"
        with open(new_bed, "wb") as o_file:
            o_file.write(bytearray([107, 27, 1, 1, 2, 3, 4]))

        # This should raise an exception
        with self.assertRaises(ValueError) as cm:
            PyPlink(new_prefix)
        self.assertEqual("not a valid BED file: {}".format(new_bed),
                         str(cm.exception))

        # Creating a new BED file with invalid second byte
        new_bed = new_prefix + ".bed"
        with open(new_bed, "wb") as o_file:
            o_file.write(bytearray([108, 28, 1, 1, 2, 3, 4]))

        # This should raise an exception
        with self.assertRaises(ValueError) as cm:
            PyPlink(new_prefix)
        self.assertEqual("not a valid BED file: {}".format(new_bed),
                         str(cm.exception))

        # Creating a new BED file not in SNP-major format
        new_bed = new_prefix + ".bed"
        with open(new_bed, "wb") as o_file:
            o_file.write(bytearray([108, 27, 0, 1, 2, 3, 4]))

        # This should raise an exception
        with self.assertRaises(ValueError) as cm:
            PyPlink(new_prefix)
        self.assertEqual("not in SNP-major format: {}".format(new_bed),
                         str(cm.exception))

    def test_missing_files(self):
        """Checks that an exception is raised when an input file is missing."""
        # Creating dummy BED/BIM/FAM files
        prefix = os.path.join(self.tmp_dir, "test_missing")
        for extension in (".bed", ".bim", ".fam"):
            with open(prefix + extension, "w") as o_file:
                pass

        # Removing the files (one by one) and checking the exception is raised
        for extension in (".bed", ".bim", ".fam"):
            os.remove(prefix + extension)
            with self.assertRaises(IOError) as cm:
                PyPlink(prefix)
            self.assertEqual("No such file: '{}'".format(prefix + extension),
                             str(cm.exception))
            with open(prefix + extension, "w") as o_file:
                pass

    def test_get_nb_markers(self):
        """Tests that the correct number of markers is returned."""
        self.assertEqual(self.pedfile.get_nb_markers(), 10)

    def test_get_nb_samples(self):
        """Tests that the correct number of samples is returned."""
        self.assertEqual(self.pedfile.get_nb_samples(), 3)

    def test_get_bim(self):
        """Tests the 'get_bim' function."""
        # The original BIM file (with the 'i' column)
        ori_bim = self.pedfile._bim

        # The expected values
        chromosomes = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        positions = [45162, 45257, 45413, 46844, 72434, 72515, 77689, 78032,
                     81468, 222077]
        cms = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        a1s = ["G", "C", "0", "A", "0", "0", "G", "0", "G", "A"]
        a2s = ["C", "T", "0", "T", "G", "C", "A", "T", "C", "G"]

        # Getting the BIM file
        bim = self.pedfile.get_bim()

        # Checking the columns
        self.assertTrue(
            set(bim.columns.values) == {"chrom", "pos", "cm", "a1", "a2"}
        )

        # Checking the indexes
        self.assertTrue(set(bim.index.values) == set(self.markers))

        # Checking the values for the markers
        zipped = zip(self.markers, chromosomes, positions, cms, a1s, a2s)
        for marker, chrom, pos, cm, a1, a2 in zipped:
            self.assertEqual(chrom, bim.loc[marker, "chrom"])
            self.assertEqual(pos, bim.loc[marker, "pos"])
            self.assertEqual(cm, bim.loc[marker, "cm"])
            self.assertEqual(a1, bim.loc[marker, "a1"])
            self.assertEqual(a2, bim.loc[marker, "a2"])

        # Comparing with the original values
        comparison = ori_bim.loc[:, ["chrom", "pos", "cm", "a1", "a2"]] == bim
        self.assertTrue(comparison.all().all())

        # Testing that changing a values in the BIM, doesn't change the value
        # in the original BIM
        bim.loc["rs4030300", "chrom"] = 2
        bim.loc["rs2949420", "cm"] = 0.1
        comparison = ori_bim.loc[:, ["chrom", "pos", "cm", "a1", "a2"]] == bim
        self.assertFalse(comparison.all().chrom)
        self.assertFalse(comparison.all().cm)
        self.assertTrue(comparison.all()[["pos", "a1", "a2"]].all())

    def test_get_fam(self):
        """Tests the 'get_fam' function."""
        # The original FAM file (with the 'byte' and 'bit' columns)
        ori_fam = self.pedfile._fam

        # The expected values
        fids = ["Sample_1", "Sample_2", "Sample_3"]
        iids = ["Sample_1", "Sample_2", "Sample_3"]
        fathers = ["0", "0", "Sample_1"]
        mothers = ["0", "0", "Sample_2"]
        genders = [1, 2, 2]
        status = [-9, -9, -9]

        # Getting the FAM file
        fam = self.pedfile.get_fam()

        # Checking the columns
        self.assertTrue(
            set(fam.columns.values) == {"fid", "iid", "father", "mother",
                                        "gender", "status"}
        )

        # Checking the values
        zipped = zip(fids, iids, fathers, mothers, genders, status)
        for i, (fid, iid, father, mother, gender, s) in enumerate(zipped):
            self.assertEqual(fid, fam.loc[i, "fid"])
            self.assertEqual(iid, fam.loc[i, "iid"])
            self.assertEqual(father, fam.loc[i, "father"])
            self.assertEqual(mother, fam.loc[i, "mother"])
            self.assertEqual(gender, fam.loc[i, "gender"])
            self.assertEqual(s, fam.loc[i, "status"])

        # Comparing with the original values
        comparison = ori_fam.loc[:, ["fid", "iid", "father", "mother",
                                     "gender", "status"]] == fam
        self.assertTrue(comparison.all().all())

        # Testing that changing a values in the FAM, doesn't change the value
        # in the original FAM
        fam.loc[2, "father"] = "0"
        fam.loc[0, "status"] = 2
        comparison = ori_fam.loc[:, ["fid", "iid", "father", "mother",
                                     "gender", "status"]] == fam
        self.assertFalse(comparison.all().father)
        self.assertFalse(comparison.all().status)
        self.assertTrue(
            comparison.all()[["fid", "iid", "mother", "gender"]].all()
        )

    def test_generator(self):
        """Testing the class as a generator."""
        # Zipping and checking
        zipped = zip(
            [i for i in zip(self.markers, self.genotypes)],
            self.pedfile,
        )
        for (e_marker, e_geno), (marker, geno) in zipped:
            self.assertEqual(e_marker, marker)
            self.assertTrue((e_geno == geno).all())

        # The generator should be empty
        remaining = [(marker, geno) for marker, geno in self.pedfile]
        self.assertEqual(0, len(remaining))

        # Just to be sure, we seek at the beginning of the file
        self.pedfile.seek(0)

    def test_seek(self):
        """Testing the seeking (for the generator)."""
        for marker, geno in self.pedfile:
            pass

        # The generator should be empty
        remaining = [(marker, geno) for marker, geno in self.pedfile]
        self.assertEqual(0, len(remaining))

        # Seeking at the second position
        zipped = zip(
            [i for i in zip(self.markers[1:], self.genotypes[1:])],
            self.pedfile,
        )
        self.pedfile.seek(1)
        for (e_marker, e_geno), (marker, geno) in zipped:
            self.assertEqual(e_marker, marker)
            self.assertTrue((e_geno == geno).all())

        # Seeking at the fourth position
        zipped = zip(
            [i for i in zip(self.markers[3:], self.genotypes[3:])],
            self.pedfile,
        )
        self.pedfile.seek(3)
        for (e_marker, e_geno), (marker, geno) in zipped:
            self.assertEqual(e_marker, marker)
            self.assertTrue((e_geno == geno).all())

        # Seeking at an invalid position
        with self.assertRaises(ValueError) as cm:
            self.pedfile.seek(-1)
        self.assertEqual("invalid position in BED: -1", str(cm.exception))

        # Seeking at an invalid position
        with self.assertRaises(ValueError) as cm:
            self.pedfile.seek(100)
        self.assertEqual("invalid position in BED: 100", str(cm.exception))

        # Just to be sure, we seek at the beginning of the file
        self.pedfile.seek(0)

    def test_iter_geno(self):
        """Tests the 'iter_geno' function."""
        zipped = zip(
            [i for i in zip(self.markers, self.genotypes)],
            self.pedfile.iter_geno(),
        )
        for (e_marker, e_geno), (marker, geno) in zipped:
            self.assertEqual(e_marker, marker)
            self.assertTrue((e_geno == geno).all())

    def test_iter_acgt_geno(self):
        """Tests the 'iter_acgt_geno" function."""
        zipped = zip(
            [i for i in zip(self.markers, self.acgt_genotypes)],
            self.pedfile.iter_acgt_geno(),
        )
        for (e_marker, e_geno), (marker, geno) in zipped:
            self.assertEqual(e_marker, marker)
            self.assertTrue((e_geno == geno).all())

    def test_iter_geno_marker(self):
        """Tests the 'iter_geno_marker' function."""
        # Getting a subset of indexes
        indexes = random.sample(range(len(self.markers)), 4)

        # Getting the markers and genotypes
        markers = [self.markers[i] for i in indexes]
        genotypes = [self.genotypes[i] for i in indexes]

        # Zipping and comparing
        zipped = zip(
            [i for i in zip(markers, genotypes)],
            self.pedfile.iter_geno_marker(markers),
        )
        for (e_marker, e_geno), (marker, geno) in zipped:
            self.assertEqual(e_marker, marker)
            self.assertTrue((e_geno == geno).all())

        # Testing a single marker
        index = random.randint(0, len(self.markers) - 1)
        e_marker = self.markers[index]
        e_geno = self.genotypes[index]
        for marker, geno in self.pedfile.iter_geno_marker(e_marker):
            self.assertEqual(e_marker, marker)
            self.assertTrue((e_geno == geno).all())

        # Adding a marker that doesn't exist
        markers.extend(["unknown_1", "unknown_2"])
        with self.assertRaises(ValueError) as cm:
            [i for i in self.pedfile.iter_geno_marker(markers)]
        self.assertEqual("['unknown_1', 'unknown_2']: marker not in BIM",
                         str(cm.exception))

    def test_iter_acgt_geno_marker(self):
        """Tests the 'iter_acgt_geno_marker' function."""
        # Getting a subset of indexes
        indexes = random.sample(range(len(self.markers)), 4)

        # Getting the markers and genotypes
        markers = [self.markers[i] for i in indexes]
        genotypes = [self.acgt_genotypes[i] for i in indexes]

        # Zipping and comparing
        zipped = zip(
            [i for i in zip(markers, genotypes)],
            self.pedfile.iter_acgt_geno_marker(markers),
        )
        for (e_marker, e_geno), (marker, geno) in zipped:
            self.assertEqual(e_marker, marker)
            self.assertTrue((e_geno == geno).all())

        # Testing a single marker
        index = random.randint(0, len(self.markers) - 1)
        e_marker = self.markers[index]
        e_geno = self.acgt_genotypes[index]
        for marker, geno in self.pedfile.iter_acgt_geno_marker(e_marker):
            self.assertEqual(e_marker, marker)
            self.assertTrue((e_geno == geno).all())

        # Adding a marker that doesn't exist
        markers.extend(["unknown_3", "unknown_4"])
        with self.assertRaises(ValueError) as cm:
            [i for i in self.pedfile.iter_acgt_geno_marker(markers)]
        self.assertEqual("['unknown_3', 'unknown_4']: marker not in BIM",
                         str(cm.exception))

    def test_repr(self):
        """Tests the object representation of the string."""
        # Counting the number of samples
        nb_samples = None
        with open(self.fam, "r") as i_file:
            nb_samples = len(i_file.read().splitlines())

        # Counting the number of markers
        nb_markers = None
        with open(self.bim, "r") as i_file:
            nb_markers = len(i_file.read().splitlines())

        # Creating the expected string representation
        e_repr = "PyPlink({:,d} samples; {:,d} markers)".format(nb_samples,
                                                                nb_markers)

        # Getting the observed string representation
        o_repr = str(self.pedfile)

        # Comparing
        self.assertEqual(e_repr, o_repr)

    def test_get_geno_marker(self):
        """Tests the 'get_geno_marker' function."""
        # Getting a random marker to test
        i = random.choice(range(len(self.markers)))
        marker = self.markers[i]
        e_geno = self.genotypes[i]

        # Getting the genotype
        o_geno = self.pedfile.get_geno_marker(marker)
        self.assertTrue((o_geno == e_geno).all())

        # Asking for an unknown marker should raise an ValueError
        with self.assertRaises(ValueError) as cm:
            self.pedfile.get_geno_marker("dummy_marker")
        self.assertEqual(
            "dummy_marker: marker not in BIM",
            str(cm.exception),
        )

    def test_get_acgt_geno_marker(self):
        """Tests the 'get_acgt_geno_marker' function."""
        # Getting a random marker to test
        i = random.choice(range(len(self.markers)))
        marker = self.markers[i]
        e_geno = self.acgt_genotypes[i]

        # Getting the genotype
        o_geno = self.pedfile.get_acgt_geno_marker(marker)
        self.assertTrue((o_geno == e_geno).all())

        # Asking for an unknown marker should raise an ValueError
        with self.assertRaises(ValueError) as cm:
            self.pedfile.get_acgt_geno_marker("dummy_marker")
        self.assertEqual("dummy_marker: marker not in BIM", str(cm.exception))

    def test_get_context_read_mode(self):
        """Tests the PyPlink object as context manager."""
        with PyPlink(self.prefix) as genotypes:
            self.assertEqual(3, len(genotypes.get_fam().head(n=3)))

    def test_invalid_mode(self):
        """Tests invalid mode when PyPlink as context manager."""
        with self.assertRaises(ValueError) as cm:
            with PyPlink(self.prefix, "u") as genotypes:
                pass
        self.assertEqual("invalid mode: 'u'", str(cm.exception))

    def test_write_binary(self):
        """Tests writing a Plink binary file."""
        # The expected genotypes
        expected_genotypes = [
            np.array([0, 0, 0, 1, 0, 1, 2], dtype=int),
            np.array([0, 0, 0, 0, -1, 0, 1], dtype=int),
            np.array([0, -1, -1, 2, 0, 0, 0], dtype=int),
        ]

        # The prefix
        test_prefix = os.path.join(self.tmp_dir, "test_write")

        # Writing the binary file
        with PyPlink(test_prefix, "w") as plink:
            for genotypes in expected_genotypes:
                plink.write_marker(genotypes)

        # Checking the file exists
        self.assertTrue(os.path.isfile(test_prefix + ".bed"))

        # Writing the FAM file
        with open(test_prefix + ".fam", "w") as o_file:
            for i in range(7):
                print("f{}".format(i+1), "s{}".format(i+1), "0", "0",
                      random.choice((1, 2)), "-9", sep=" ", file=o_file)

        # Writing the BIM file
        with open(test_prefix + ".bim", "w") as o_file:
            for i in range(3):
                print(i+1, "m{}".format(i+1), "0", i+1, "T", "A", sep="\t",
                      file=o_file)

        # Reading the written binary file
        plink = PyPlink(test_prefix)
        for i, (marker, genotypes) in enumerate(plink):
            self.assertEqual("m{}".format(i+1), marker)
            self.assertTrue((expected_genotypes[i] == genotypes).all())

    def test_write_binary_error(self):
        """Tests writing a binary file, with different number of sample."""
        # The expected genotypes
        expected_genotypes = [
            np.array([0, 0, 0, 1, 0, 1, 2], dtype=int),
            np.array([0, 0, 0, 0, -1, 0], dtype=int),
            np.array([0, -1, -1, 2, 0, 0, 0], dtype=int),
        ]

        # The prefix
        test_prefix = os.path.join(self.tmp_dir, "test_write")

        # Writing the binary file
        with self.assertRaises(ValueError) as cm:
            with PyPlink(test_prefix, "w") as plink:
                for genotypes in expected_genotypes:
                    plink.write_marker(genotypes)
        self.assertEqual("7 samples expected, got 6", str(cm.exception))

    def test_grouper_padding(self):
        """Tests the _grouper function (when padding is required)."""
        expected_chunks = [
            (0, 1, 2),
            (3, 4, 5),
            (6, 7, 8),
            (9, 0, 0),
        ]
        observed_chunks = PyPlink._grouper(range(10), 3)
        for expected, observed in zip(expected_chunks, observed_chunks):
            self.assertEqual(expected, observed)

    def test_grouper_no_padding(self):
        """Tests the _grouper function (when padding is not required)."""
        expected_chunks = [
            (0, 1, 2, 3, 4),
            (5, 6, 7, 8, 9),
        ]
        observed_chunks = PyPlink._grouper(range(10), 5)
        for expected, observed in zip(expected_chunks, observed_chunks):
            self.assertEqual(expected, observed)
