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
import sys
import stat
import random
import shutil
import zipfile
import platform
import unittest
from tempfile import mkdtemp
from io import UnsupportedOperation
from distutils.spawn import find_executable
from subprocess import check_call, PIPE, CalledProcessError

try:
    from itertools import zip_longest as zip
except ImportError:
    from itertools import izip_longest as zip

try:
    from urllib.request import urlretrieve
except ImportError:
    from urllib import urlretrieve

try:
    from unittest import mock
except ImportError:
    import mock

from pkg_resources import resource_filename

import numpy as np
import pandas as pd

from six.moves import range

from .. import pyplink


def get_plink(tmp_dir):
    """Gets the Plink binary, if required."""
    # Checking if Plink is in the path
    plink_path = "plink"
    if platform.system() == "Windows":
        plink_path += ".exe"

    if find_executable(plink_path) is None:
        print("Downloading Plink", file=sys.stderr)

        # The url for each platform
        url = ("http://statgen.org/wp-content/uploads/Softwares/"
               "plink-1.0.7/{filename}")

        # Getting the name of the file
        filename = ""
        if platform.system() == "Windows":
            filename = "plink-1.07-dos.zip"
        elif platform.system() == "Darwin":
            filename = "plink-1.07-mac-intel.zip"
        elif platform.system() == "Linux":
            if platform.architecture()[0].startswith("32"):
                filename = "plink-1.07-i686.zip"
            elif platform.architecture()[0].startswith("64"):
                filename = "plink-1.07-x86_64.zip"
            else:
                return None, "System not compatible for Plink"
        else:
            return None, "System not compatible for Plink"

        # Downloading Plink
        zip_path = os.path.join(tmp_dir, filename)
        try:
            urlretrieve(
                url.format(filename=filename),
                zip_path,
            )
        except:
            return None, "Plink's URL is not available"

        # Unzipping Plink
        with zipfile.ZipFile(zip_path, "r") as z:
            z.extractall(tmp_dir)
        plink_path = os.path.join(tmp_dir, os.path.splitext(filename)[0],
                                  plink_path)
        if not os.path.isfile(plink_path):
            return None, "Cannot use Plink"

        # Making the script executable
        if platform.system() in {"Darwin", "Linux"}:
            os.chmod(plink_path, stat.S_IRWXU)

    # Testing Plink works
    try:
        check_call([
            plink_path,
            "--noweb",
            "--help",
            "--out", os.path.join(tmp_dir, "execution_test")
        ], stdout=PIPE, stderr=PIPE)
    except CalledProcessError:
        return None, "Plink cannot be properly used"
    except IOError:
        return None, "Plink was not properly installed"

    return plink_path, "OK"


class TestPyPlink(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Creating a temporary directory
        cls.tmp_dir = mkdtemp(prefix="pyplink_test_")

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

        # Getting Plink
        cls.plink_path, cls.plink_message = get_plink(cls.tmp_dir)

    def setUp(self):
        # Reading the plink binary file
        self.pedfile = pyplink.PyPlink(self.prefix)

    @classmethod
    def tearDownClass(cls):
        # Cleaning the temporary directory
        shutil.rmtree(cls.tmp_dir)

    def tearDown(self):
        # Closing the PyPlink object
        self.pedfile.close()

    def test_pyplink_object_integrity(self):
        """Checks the integrity of the PyPlink object."""
        # Checking the name of the BED file
        self.assertTrue(hasattr(self.pedfile, "bed_filename"))
        self.assertEqual(self.bed, self.pedfile.bed_filename)

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
        # Changing the BIM to None
        ori = self.pedfile._bim
        self.pedfile._bim = None
        with self.assertRaises(RuntimeError) as cm:
            self.pedfile._read_bed()
        self.assertEqual("no BIM or FAM file were read", str(cm.exception))
        self.pedfile._bim = ori

        # Changing the FAM to None
        ori = self.pedfile._fam
        self.pedfile._fam = None
        with self.assertRaises(RuntimeError) as cm:
            self.pedfile._read_bed()
        self.assertEqual("no BIM or FAM file were read", str(cm.exception))
        self.pedfile._fam = ori

    def test_pyplink_bad_bed(self):
        """Checks what happens when we read a bad BED file."""
        # The new file prefix
        new_prefix = os.path.join(self.tmp_dir, "bad_data")

        # Copying the FAM file
        new_fam = new_prefix + ".fam"
        with open(new_fam, "w") as o_file, open(self.fam, "r") as i_file:
            o_file.write(i_file.read())

        # Copying the BIM file
        new_bim = new_prefix + ".bim"
        with open(new_bim, "w") as o_file, open(self.bim, "r") as i_file:
            o_file.write(i_file.read())

        # Creating a new BED file with invalid number of bytes
        new_bed = new_prefix + ".bed"
        with open(new_bed, "wb") as o_file:
            o_file.write(bytearray([108, 27, 1, 1, 2, 3, 4]))

        # This should raise an exception
        with self.assertRaises(ValueError) as cm:
            pyplink.PyPlink(new_prefix)
        self.assertEqual("invalid number of entries: corrupted BED?",
                         str(cm.exception))

        # Creating a new BED file with invalid first byte
        new_bed = new_prefix + ".bed"
        with open(new_bed, "wb") as o_file:
            o_file.write(bytearray([107, 27, 1, 1, 2, 3, 4]))

        # This should raise an exception
        with self.assertRaises(ValueError) as cm:
            pyplink.PyPlink(new_prefix)
        self.assertEqual("not a valid BED file: {}".format(new_bed),
                         str(cm.exception))

        # Creating a new BED file with invalid second byte
        new_bed = new_prefix + ".bed"
        with open(new_bed, "wb") as o_file:
            o_file.write(bytearray([108, 28, 1, 1, 2, 3, 4]))

        # This should raise an exception
        with self.assertRaises(ValueError) as cm:
            pyplink.PyPlink(new_prefix)
        self.assertEqual("not a valid BED file: {}".format(new_bed),
                         str(cm.exception))

        # Creating a new BED file not in SNP-major format
        new_bed = new_prefix + ".bed"
        with open(new_bed, "wb") as o_file:
            o_file.write(bytearray([108, 27, 0, 1, 2, 3, 4]))

        # This should raise an exception
        with self.assertRaises(ValueError) as cm:
            pyplink.PyPlink(new_prefix)
        self.assertEqual(
            "not in SNP-major format (please recode): {}".format(new_bed),
            str(cm.exception),
        )

    def test_missing_files(self):
        """Checks that an exception is raised when an input file is missing."""
        # Creating dummy BED/BIM/FAM files
        prefix = os.path.join(self.tmp_dir, "test_missing")
        for extension in (".bed", ".bim", ".fam"):
            with open(prefix + extension, "w"):
                pass

        # Removing the files (one by one) and checking the exception is raised
        for extension in (".bed", ".bim", ".fam"):
            os.remove(prefix + extension)
            with self.assertRaises(IOError) as cm:
                pyplink.PyPlink(prefix)
            self.assertEqual("No such file: '{}'".format(prefix + extension),
                             str(cm.exception))
            with open(prefix + extension, "w"):
                pass

    def test_get_nb_markers(self):
        """Tests that the correct number of markers is returned."""
        self.assertEqual(self.pedfile.get_nb_markers(), 10)

    def test_get_nb_markers_w_mode(self):
        """Tests that an exception is raised if in write mode."""
        with self.assertRaises(UnsupportedOperation) as cm:
            # Creating the dummy PyPlink object
            prefix = os.path.join(self.tmp_dir, "test_error")
            with pyplink.PyPlink(prefix, "w") as p:
                p.get_nb_markers()
        self.assertEqual("not available in 'w' mode", str(cm.exception))

    def test_get_nb_samples(self):
        """Tests that the correct number of samples is returned."""
        self.assertEqual(self.pedfile.get_nb_samples(), 3)

    def test_get_nb_samples_w_mode(self):
        """Tests that an exception is raised if in write mode."""
        with self.assertRaises(UnsupportedOperation) as cm:
            # Creating the dummy PyPlink object
            prefix = os.path.join(self.tmp_dir, "test_error")
            with pyplink.PyPlink(prefix, "w") as p:
                p.get_nb_samples()
        self.assertEqual("not available in 'w' mode", str(cm.exception))

    def test_get_bim(self):
        """Tests the 'get_bim' function."""
        # The original BIM file (with the 'i' column)
        ori_bim = self.pedfile._bim

        # The expected values
        chromosomes = [1, 2, 3, 4, 4, 5, 6, 6, 6, 8]
        positions = [45162, 45257, 45413, 46844, 72434, 72515, 77689, 78032,
                     81468, 222077]
        cms = [0, 1, 1, 2, 2, 3, 4, 4, 5, 6]
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

    def test_get_bim_w_mode(self):
        """Tests that an exception is raised if in write mode."""
        with self.assertRaises(UnsupportedOperation) as cm:
            # Creating the dummy PyPlink object
            prefix = os.path.join(self.tmp_dir, "test_error")
            with pyplink.PyPlink(prefix, "w") as p:
                p.get_bim()
        self.assertEqual("not available in 'w' mode", str(cm.exception))

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

    def test_get_fam_w_mode(self):
        """Tests that an exception is raised if in write mode."""
        with self.assertRaises(UnsupportedOperation) as cm:
            # Creating the dummy PyPlink object
            prefix = os.path.join(self.tmp_dir, "test_error")
            with pyplink.PyPlink(prefix, "w") as p:
                p.get_fam()
        self.assertEqual("not available in 'w' mode", str(cm.exception))

    def test_generator(self):
        """Testing the class as a generator."""
        # Zipping and checking
        zipped = zip(
            [i for i in zip(self.markers, self.genotypes)],
            self.pedfile,
        )
        for (e_marker, e_geno), (marker, geno) in zipped:
            self.assertEqual(e_marker, marker)
            np.testing.assert_array_equal(e_geno, geno)

        # The generator should be empty
        remaining = [(marker, geno) for marker, geno in self.pedfile]
        self.assertEqual(0, len(remaining))

    def test_generator_w_mode(self):
        """Tests that an exception is raised if in write mode."""
        with self.assertRaises(UnsupportedOperation) as cm:
            # Creating the dummy PyPlink object
            prefix = os.path.join(self.tmp_dir, "test_error")
            with pyplink.PyPlink(prefix, "w") as p:
                marker, genotypes = next(p)
        self.assertEqual("not available in 'w' mode", str(cm.exception))

    def test_next(self):
        """Tests that an exception is raised when calling next in w mode."""
        marker, genotypes = self.pedfile.next()

        # Comparing
        self.assertEqual(self.markers[0], marker)
        np.testing.assert_array_equal(self.genotypes[0], genotypes)

    def test_next_w_mode(self):
        """Tests that an exception is raised when calling next in w mode."""
        with self.assertRaises(UnsupportedOperation) as cm:
            # Creating the dummy PyPlink object
            prefix = os.path.join(self.tmp_dir, "test_error")
            with pyplink.PyPlink(prefix, "w") as p:
                p.next()
        self.assertEqual("not available in 'w' mode", str(cm.exception))

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
            np.testing.assert_array_equal(e_geno, geno)

        # Seeking at the fourth position
        zipped = zip(
            [i for i in zip(self.markers[3:], self.genotypes[3:])],
            self.pedfile,
        )
        self.pedfile.seek(3)
        for (e_marker, e_geno), (marker, geno) in zipped:
            self.assertEqual(e_marker, marker)
            np.testing.assert_array_equal(e_geno, geno)

        # Seeking at the tenth position
        zipped = zip(
            [i for i in zip(self.markers[9:], self.genotypes[9:])],
            self.pedfile,
        )
        self.pedfile.seek(9)
        for (e_marker, e_geno), (marker, geno) in zipped:
            self.assertEqual(e_marker, marker)
            np.testing.assert_array_equal(e_geno, geno)

        # Seeking at an invalid position
        with self.assertRaises(ValueError) as cm:
            self.pedfile.seek(-1)
        self.assertEqual("invalid position in BED: -1", str(cm.exception))

        # Seeking at an invalid position
        with self.assertRaises(ValueError) as cm:
            self.pedfile.seek(100)
        self.assertEqual("invalid position in BED: 100", str(cm.exception))

        # Seeking at an invalid position
        with self.assertRaises(ValueError) as cm:
            self.pedfile.seek(10)
        self.assertEqual("invalid position in BED: 10", str(cm.exception))

    def test_seek_w_mode(self):
        """Tests that an exception is raised if in write mode."""
        with self.assertRaises(UnsupportedOperation) as cm:
            # Creating the dummy PyPlink object
            prefix = os.path.join(self.tmp_dir, "test_error")
            with pyplink.PyPlink(prefix, "w") as p:
                p.seek(100)
        self.assertEqual("not available in 'w' mode", str(cm.exception))

    def test_iter_geno(self):
        """Tests the 'iter_geno' function."""
        zipped = zip(
            [i for i in zip(self.markers, self.genotypes)],
            self.pedfile.iter_geno(),
        )
        for (e_marker, e_geno), (marker, geno) in zipped:
            self.assertEqual(e_marker, marker)
            np.testing.assert_array_equal(e_geno, geno)

    def test_iter_geno_w_mode(self):
        """Tests that an exception is raised if in write mode."""
        with self.assertRaises(UnsupportedOperation) as cm:
            # Creating the dummy PyPlink object
            prefix = os.path.join(self.tmp_dir, "test_error")
            with pyplink.PyPlink(prefix, "w") as p:
                marker, genotypes = next(p.iter_geno())
        self.assertEqual("not available in 'w' mode", str(cm.exception))

    def test_iter_acgt_geno(self):
        """Tests the 'iter_acgt_geno" function."""
        zipped = zip(
            [i for i in zip(self.markers, self.acgt_genotypes)],
            self.pedfile.iter_acgt_geno(),
        )
        for (e_marker, e_geno), (marker, geno) in zipped:
            self.assertEqual(e_marker, marker)
            np.testing.assert_array_equal(e_geno, geno)

    def test_iter_acgt_geno_w_mode(self):
        """Tests that an exception is raised if in write mode."""
        with self.assertRaises(UnsupportedOperation) as cm:
            # Creating the dummy PyPlink object
            prefix = os.path.join(self.tmp_dir, "test_error")
            with pyplink.PyPlink(prefix, "w") as p:
                marker, genotypes = next(p.iter_acgt_geno())
        self.assertEqual("not available in 'w' mode", str(cm.exception))

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
            np.testing.assert_array_equal(e_geno, geno)

        # Testing a single marker
        index = random.randint(0, len(self.markers) - 1)
        e_marker = self.markers[index]
        e_geno = self.genotypes[index]
        for marker, geno in self.pedfile.iter_geno_marker(e_marker):
            self.assertEqual(e_marker, marker)
            np.testing.assert_array_equal(e_geno, geno)

        # Adding a marker that doesn't exist
        markers.extend(["unknown_1", "unknown_2"])
        with self.assertRaises(ValueError) as cm:
            [i for i in self.pedfile.iter_geno_marker(markers)]
        self.assertEqual("unknown_1: marker not in BIM", str(cm.exception))

    def test_iter_geno_marker_w_mode(self):
        """Tests that an exception is raised if in write mode."""
        with self.assertRaises(UnsupportedOperation) as cm:
            # Creating the dummy PyPlink object
            prefix = os.path.join(self.tmp_dir, "test_error")
            with pyplink.PyPlink(prefix, "w") as p:
                marker, genotypes = next(p.iter_geno_marker(["M1", "M2"]))
        self.assertEqual("not available in 'w' mode", str(cm.exception))

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
            np.testing.assert_array_equal(e_geno, geno)

        # Testing a single marker
        index = random.randint(0, len(self.markers) - 1)
        e_marker = self.markers[index]
        e_geno = self.acgt_genotypes[index]
        for marker, geno in self.pedfile.iter_acgt_geno_marker(e_marker):
            self.assertEqual(e_marker, marker)
            np.testing.assert_array_equal(e_geno, geno)

        # Adding a marker that doesn't exist
        markers.extend(["unknown_3", "unknown_4"])
        with self.assertRaises(ValueError) as cm:
            [i for i in self.pedfile.iter_acgt_geno_marker(markers)]
        self.assertEqual("unknown_3: marker not in BIM", str(cm.exception))

    def test_iter_acgt_geno_marker_w_mode(self):
        """Tests that an exception is raised if in write mode."""
        with self.assertRaises(UnsupportedOperation) as cm:
            # Creating the dummy PyPlink object
            prefix = os.path.join(self.tmp_dir, "test_error")
            with pyplink.PyPlink(prefix, "w") as p:
                marker, genotypes = next(p.iter_acgt_geno_marker(["M1", "M2"]))
        self.assertEqual("not available in 'w' mode", str(cm.exception))

    def test_repr_r_mode(self):
        """Tests the object representation of the string (r mode)."""
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

    def test_repr_w_mode(self):
        """Tests the object representation of the string (w mode)."""
        # The expected representation
        e_repr = 'PyPlink(mode="w")'

        # Creating the dummy PyPlink object
        prefix = os.path.join(self.tmp_dir, "test_repr")
        with pyplink.PyPlink(prefix, "w") as pedfile:
            # Comparing the expected with the observed representation
            o_repr = str(pedfile)
            self.assertEqual(e_repr, o_repr)

    def test_get_geno_marker(self):
        """Tests the 'get_geno_marker' function."""
        # Getting a random marker to test
        i = random.choice(range(len(self.markers)))
        marker = self.markers[i]
        e_geno = self.genotypes[i]

        # Getting the genotype
        o_geno = self.pedfile.get_geno_marker(marker)
        np.testing.assert_array_equal(o_geno, e_geno)

        # Asking for an unknown marker should raise an ValueError
        with self.assertRaises(ValueError) as cm:
            self.pedfile.get_geno_marker("dummy_marker")
        self.assertEqual(
            "dummy_marker: marker not in BIM",
            str(cm.exception),
        )

    def test_get_geno_marker_w_mode(self):
        """Tests that an exception is raised if in write mode."""
        with self.assertRaises(UnsupportedOperation) as cm:
            # Creating the dummy PyPlink object
            prefix = os.path.join(self.tmp_dir, "test_error")
            with pyplink.PyPlink(prefix, "w") as p:
                p.get_geno_marker("M1")
        self.assertEqual("not available in 'w' mode", str(cm.exception))

    def test_get_iter_w_mode(self):
        """Tests that an exception is raised if in write mode."""
        with self.assertRaises(UnsupportedOperation) as cm:
            # Creating the dummy PyPlink object
            prefix = os.path.join(self.tmp_dir, "test_error")
            with pyplink.PyPlink(prefix, "w") as p:
                iter(p)
        self.assertEqual("not available in 'w' mode", str(cm.exception))

    def test_get_acgt_geno_marker(self):
        """Tests the 'get_acgt_geno_marker' function."""
        # Getting a random marker to test
        i = random.choice(range(len(self.markers)))
        marker = self.markers[i]
        e_geno = self.acgt_genotypes[i]

        # Getting the genotype
        o_geno = self.pedfile.get_acgt_geno_marker(marker)
        np.testing.assert_array_equal(o_geno, e_geno)

        # Asking for an unknown marker should raise an ValueError
        with self.assertRaises(ValueError) as cm:
            self.pedfile.get_acgt_geno_marker("dummy_marker")
        self.assertEqual("dummy_marker: marker not in BIM", str(cm.exception))

    def test_get_acgt_geno_marker_w_mode(self):
        """Tests that an exception is raised if in write mode."""
        with self.assertRaises(UnsupportedOperation) as cm:
            # Creating the dummy PyPlink object
            prefix = os.path.join(self.tmp_dir, "test_error")
            with pyplink.PyPlink(prefix, "w") as p:
                p.get_acgt_geno_marker("M1")
        self.assertEqual("not available in 'w' mode", str(cm.exception))

    def test_get_context_read_mode(self):
        """Tests the PyPlink object as context manager."""
        with pyplink.PyPlink(self.prefix) as genotypes:
            self.assertEqual(3, len(genotypes.get_fam().head(n=3)))

    def test_invalid_mode(self):
        """Tests invalid mode when PyPlink as context manager."""
        with self.assertRaises(ValueError) as cm:
            pyplink.PyPlink(self.prefix, "u")
        self.assertEqual("invalid mode: 'u'", str(cm.exception))

    def test_write_binary(self):
        """Tests writing a Plink binary file."""
        # The expected genotypes
        expected_genotypes = [
            np.array([0,  0,  0, 1,  0, 1, 2], dtype=int),
            np.array([0,  0,  0, 0, -1, 0, 1], dtype=int),
            np.array([0, -1, -1, 2,  0, 0, 0], dtype=int),
        ]

        # The prefix
        test_prefix = os.path.join(self.tmp_dir, "test_write")

        # Writing the binary file
        with pyplink.PyPlink(test_prefix, "w") as pedfile:
            for genotypes in expected_genotypes:
                pedfile.write_genotypes(genotypes)

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
        with pyplink.PyPlink(test_prefix) as pedfile:
            for i, (marker, genotypes) in enumerate(pedfile):
                self.assertEqual("m{}".format(i+1), marker)
                np.testing.assert_array_equal(expected_genotypes[i], genotypes)

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
            with pyplink.PyPlink(test_prefix, "w") as pedfile:
                pedfile.write_genotypes(expected_genotypes[0])  # 7 genotypes
                pedfile.write_genotypes(expected_genotypes[1])  # 6 genotypes
        self.assertEqual("7 samples expected, got 6", str(cm.exception))

    def test_grouper_padding(self):
        """Tests the _grouper function (when padding is required)."""
        expected_chunks = [
            (0, 1, 2),
            (3, 4, 5),
            (6, 7, 8),
            (9, 0, 0),
        ]
        observed_chunks = pyplink.PyPlink._grouper(range(10), 3)
        for expected, observed in zip(expected_chunks, observed_chunks):
            self.assertEqual(expected, observed)

    def test_grouper_no_padding(self):
        """Tests the _grouper function (when padding is not required)."""
        expected_chunks = [
            (0, 1, 2, 3, 4),
            (5, 6, 7, 8, 9),
        ]
        observed_chunks = pyplink.PyPlink._grouper(range(10), 5)
        for expected, observed in zip(expected_chunks, observed_chunks):
            self.assertEqual(expected, observed)

    @unittest.skipIf(platform.system() not in {"Darwin", "Linux", "Windows"},
                     "Plink not available for {}".format(platform.system()))
    def test_with_plink(self):
        """Tests to read a binary file using Plink."""
        # Checking if we need to skip
        if self.plink_path is None:
            self.skipTest(self.plink_message)

        # Creating the BED file
        all_genotypes = [
            [0, 1, 0, 0, -1, 0, 1, 0, 0, 2],
            [2, 1, 2, 2,  2, 2, 2, 1, 0, 1],
            [0, 0, 0, 0,  0, 1, 0, 0, 0, 0],
        ]
        prefix = os.path.join(self.tmp_dir, "test_output")
        with pyplink.PyPlink(prefix, "w") as pedfile:
            for genotypes in all_genotypes:
                pedfile.write_genotypes(genotypes)

        # Creating the FAM file
        fam_content = [
            ["F0", "S0", "0", "0", "1", "-9"],
            ["F1", "S1", "0", "0", "1", "-9"],
            ["F2", "S2", "0", "0", "2", "-9"],
            ["F3", "S3", "0", "0", "1", "-9"],
            ["F4", "S4", "0", "0", "1", "-9"],
            ["F5", "S5", "0", "0", "2", "-9"],
            ["F6", "S6", "0", "0", "1", "-9"],
            ["F7", "S7", "0", "0", "0", "-9"],
            ["F8", "S8", "0", "0", "1", "-9"],
            ["F9", "S9", "0", "0", "2", "-9"],
        ]
        with open(prefix + ".fam", "w") as o_file:
            for sample in fam_content:
                print(*sample, sep=" ", file=o_file)

        # Creating the BIM file
        bim_content = [
            ["1", "M0", "0", "123", "A", "G"],
            ["1", "M1", "0", "124", "C", "T"],
            ["2", "M2", "0", "117", "G", "C"],
        ]
        with open(prefix + ".bim", "w") as o_file:
            for marker in bim_content:
                print(*marker, sep="\t", file=o_file)

        # Creating a transposed pedfile using Plink
        out_prefix = prefix + "_transposed"
        try:
            check_call([
                self.plink_path,
                "--noweb",
                "--bfile", prefix,
                "--recode", "--transpose", "--tab",
                "--out", out_prefix,
            ], stdout=PIPE, stderr=PIPE)
        except CalledProcessError:
            self.fail("Plink could not recode file")

        # Checking the two files exists
        self.assertTrue(os.path.isfile(out_prefix + ".tped"))
        self.assertTrue(os.path.isfile(out_prefix + ".tfam"))

        # Checking the content of the TFAM file
        expected = "\n".join("\t".join(sample) for sample in fam_content)
        with open(out_prefix + ".tfam", "r") as i_file:
            self.assertEqual(expected + "\n", i_file.read())

        # Checking the content of the TPED file
        with open(out_prefix + ".tped", "r") as i_file:
            # Checking the first marker
            marker_1 = i_file.readline().rstrip("\r\n").split("\t")
            self.assertEqual(["1", "M0", "0", "123"], marker_1[:4])
            self.assertEqual(["G G", "A G", "G G", "G G", "0 0", "G G", "A G",
                              "G G", "G G", "A A"],
                             marker_1[4:])

            # Checking the second marker
            marker_2 = i_file.readline().rstrip("\r\n").split("\t")
            self.assertEqual(["1", "M1", "0", "124"], marker_2[:4])
            self.assertEqual(["C C", "T C", "C C", "C C", "C C", "C C", "C C",
                              "T C", "T T", "T C"],
                             marker_2[4:])

            # Checking the third marker
            marker_3 = i_file.readline().rstrip("\r\n").split("\t")
            self.assertEqual(["2", "M2", "0", "117"], marker_3[:4])
            self.assertEqual(["C C", "C C", "C C", "C C", "C C", "G C", "C C",
                              "C C", "C C", "C C"],
                             marker_3[4:])

            # Checking this is the end of the file
            self.assertEqual("", i_file.readline())

    @unittest.skipIf(platform.system() not in {"Darwin", "Linux", "Windows"},
                     "Plink not available for {}".format(platform.system()))
    def test_with_plink_individual_major(self):
        """Tests to read a binary file (INDIVIDUAL-major) using Plink."""
        # Checking if we need to skip
        if self.plink_path is None:
            self.skipTest(self.plink_message)

        # The genotypes
        all_genotypes = [
            [0, 1, 0, 0, -1, 0, 1, 0, 0, 2],
            [2, 1, 2, 2,  2, 2, 2, 1, 0, 1],
            [0, 0, 0, 0,  0, 1, 0, 0, 0, 0],
        ]
        transposed_genotypes = [
            [row[i] for row in all_genotypes]
            for i in range(len(all_genotypes[0]))
        ]

        # Creating the BED file (INDIVIDUAL-major)
        prefix = os.path.join(self.tmp_dir, "test_output")
        with pyplink.PyPlink(prefix, "w", "INDIVIDUAL-major") as pedfile:
            for genotypes in transposed_genotypes:
                pedfile.write_genotypes(genotypes)

        # Creating the FAM file
        fam_content = [
            ["F0", "S0", "0", "0", "1", "-9"],
            ["F1", "S1", "0", "0", "1", "-9"],
            ["F2", "S2", "0", "0", "2", "-9"],
            ["F3", "S3", "0", "0", "1", "-9"],
            ["F4", "S4", "0", "0", "1", "-9"],
            ["F5", "S5", "0", "0", "2", "-9"],
            ["F6", "S6", "0", "0", "1", "-9"],
            ["F7", "S7", "0", "0", "0", "-9"],
            ["F8", "S8", "0", "0", "1", "-9"],
            ["F9", "S9", "0", "0", "2", "-9"],
        ]
        with open(prefix + ".fam", "w") as o_file:
            for sample in fam_content:
                print(*sample, sep=" ", file=o_file)

        # Creating the BIM file
        bim_content = [
            ["1", "M0", "0", "123", "A", "G"],
            ["1", "M1", "0", "124", "C", "T"],
            ["2", "M2", "0", "117", "G", "C"],
        ]
        with open(prefix + ".bim", "w") as o_file:
            for marker in bim_content:
                print(*marker, sep="\t", file=o_file)

        # Creating a transposed pedfile using Plink
        out_prefix = prefix + "_transposed"
        try:
            check_call([
                self.plink_path,
                "--noweb",
                "--bfile", prefix,
                "--recode", "--transpose", "--tab",
                "--out", out_prefix,
            ], stdout=PIPE, stderr=PIPE)
        except CalledProcessError:
            self.fail("Plink could not recode file")

        # Checking the two files exists
        self.assertTrue(os.path.isfile(out_prefix + ".tped"))
        self.assertTrue(os.path.isfile(out_prefix + ".tfam"))

        # Checking the content of the TFAM file
        expected = "\n".join("\t".join(sample) for sample in fam_content)
        with open(out_prefix + ".tfam", "r") as i_file:
            self.assertEqual(expected + "\n", i_file.read())

        # Checking the content of the TPED file
        with open(out_prefix + ".tped", "r") as i_file:
            # Checking the first marker
            marker_1 = i_file.readline().rstrip("\r\n").split("\t")
            self.assertEqual(["1", "M0", "0", "123"], marker_1[:4])
            self.assertEqual(["G G", "A G", "G G", "G G", "0 0", "G G", "A G",
                              "G G", "G G", "A A"],
                             marker_1[4:])

            # Checking the second marker
            marker_2 = i_file.readline().rstrip("\r\n").split("\t")
            self.assertEqual(["1", "M1", "0", "124"], marker_2[:4])
            self.assertEqual(["C C", "T C", "C C", "C C", "C C", "C C", "C C",
                              "T C", "T T", "T C"],
                             marker_2[4:])

            # Checking the third marker
            marker_3 = i_file.readline().rstrip("\r\n").split("\t")
            self.assertEqual(["2", "M2", "0", "117"], marker_3[:4])
            self.assertEqual(["C C", "C C", "C C", "C C", "C C", "G C", "C C",
                              "C C", "C C", "C C"],
                             marker_3[4:])

            # Checking this is the end of the file
            self.assertEqual("", i_file.readline())

    def test_wrong_bed_format(self):
        """Tests opening a BED file with unknown format."""
        with self.assertRaises(ValueError) as cm:
            pyplink.PyPlink(self.prefix, bed_format="UNKNOWN-major")
        self.assertEqual(
            "invalid bed format: UNKNOWN-major",
            str(cm.exception),
        )

    def test_invalid_bed_format_with_r_mode(self):
        """Tests an invalid BED format with r mode."""
        with self.assertRaises(ValueError) as cm:
            pyplink.PyPlink(self.prefix, bed_format="INDIVIDUAL-major")
        self.assertEqual(
            "only SNP-major format is supported with mode 'r'",
            str(cm.exception),
        )

    def test_write_genotypes_in_r_mode(self):
        """Tests to use the 'write_genotypes' function in read mode."""
        with self.assertRaises(UnsupportedOperation) as cm:
            self.pedfile.write_genotypes([0, 0, 0])
        self.assertEqual("not available in 'r' mode", str(cm.exception))

    @mock.patch.object(pyplink, "logger")
    def test_dup_markers(self, pyplink_logger):
        """Tests when there are duplicated markers."""
        # Checking the original one has no duplicates
        self.assertEqual(len(self.pedfile.get_duplicated_markers()), 0)

        # Copying the BED and the FAM to the temporary directory
        new_prefix = os.path.join(self.tmp_dir, "with_dup")
        for suffix in (".bed", ".fam"):
            shutil.copyfile(self.prefix + suffix, new_prefix + suffix)

        # Copying the BIM file to the temporary directory
        shutil.copyfile(self.prefix + "_with_dup.bim", new_prefix + ".bim")

        # Reading the new files
        pedfile = pyplink.PyPlink(new_prefix)

        # Checking a warning was called
        self.assertTrue(pyplink_logger.warning.called)

        # Checking the BIM
        chromosomes = [1, 2, 3, 4, 4, 5, 6, 6, 6, 8]
        markers = ["rs10399749", "rs2949420:dup1", "rs2949421", "rs2691310",
                   "rs4030303:dup1", "rs4030303:dup2", "rs4030303:dup3",
                   "rs940550:dup1", "rs940550:dup2", "rs2949420:dup2"]
        positions = [45162, 45257, 45413, 46844, 72434, 72515, 77689, 78032,
                     81468, 222077]
        cms = [0, 1, 1, 2, 2, 3, 4, 4, 5, 6]
        a1s = ["G", "C", "0", "A", "0", "0", "G", "0", "G", "A"]
        a2s = ["C", "T", "0", "T", "G", "C", "A", "T", "C", "G"]

        # Getting the BIM file
        bim = pedfile.get_bim()

        # Checking the columns
        self.assertTrue(
            set(bim.columns.values) == {"chrom", "pos", "cm", "a1", "a2"}
        )

        # Checking the indexes
        self.assertTrue(set(bim.index.values) == set(markers))

        # Checking the values for the markers
        zipped = zip(markers, chromosomes, positions, cms, a1s, a2s)
        for marker, chrom, pos, cm, a1, a2 in zipped:
            self.assertEqual(chrom, bim.loc[marker, "chrom"])
            self.assertEqual(pos, bim.loc[marker, "pos"])
            self.assertEqual(cm, bim.loc[marker, "cm"])
            self.assertEqual(a1, bim.loc[marker, "a1"])
            self.assertEqual(a2, bim.loc[marker, "a2"])

        # Checking only one duplicated markers
        for i, marker in enumerate(markers):
            geno = pedfile.get_geno_marker(marker)
            np.testing.assert_array_equal(geno, self.genotypes[i])

        # Checking the list of duplicated markers
        self.assertEqual(
            set(m.split(":")[0] for m in markers if ":" in m),
            pedfile.get_duplicated_markers(),
        )

        # Closing the file
        pedfile.close()
