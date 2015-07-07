"""Module that reads binary Plink files."""

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


import os

import numpy as np
import pandas as pd


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014 Louis-Philippe Lemieux Perreault"
__license__ = "MIT"


__all__ = ["PyPlink"]


# The recoding values
_geno_recode = {1: -1,  # Unknown genotype
                2: 1,   # Heterozygous genotype
                0: 2,   # Homozygous A1
                3: 0}   # Homozygous A2


class PyPlink(object):
    """Reads and store a set of binary Plink files."""

    # The genotypes values
    _geno_values = np.array(
        [
            [_geno_recode[(i >> j) & 3] for j in range(0, 7, 2)]
            for i in range(256)
        ],
        dtype=np.int8,
    )

    def __init__(self, i_prefix):
        """Initializes a new PyPlink object.

        :param i_prefix: the prefix of the binary Plink files
        :type i_prefix: str

        Reads binary Plink files (BED, BIM and FAM) and create the objects.

        """
        # These are the name of the input files
        self.bed_filename = "{}.bed".format(i_prefix)
        self.bim_filename = "{}.bim".format(i_prefix)
        self.fam_filename = "{}.fam".format(i_prefix)

        # Checking that all the files exists (otherwise, there's nothing to do)
        for filename in (self.bed_filename, self.bim_filename,
                         self.fam_filename):
            if not os.path.isfile(filename):
                raise IOError("No such file: '{}'".format(filename))

        # Setting BIM and FAM to None
        self._bim = None
        self._fam = None

        # Reading the input files
        self._read_bim()
        self._read_fam()
        self._read_bed()

        # Where we're at
        self._n = 0

    def __iter__(self):
        """The __iter__ function."""
        return self

    def __next__(self):
        """The __next__ function."""
        return self.next()

    def next(self):
        """The next function."""
        if self._n < self._nb_markers:
            self._n += 1

            # We want to return information about the marker and the genotypes
            geno = self._geno_values[self._bed[self._n - 1]].flatten(order="C")
            return self._markers[self._n - 1], geno[:self._nb_samples]
        else:
            raise StopIteration()

    def seek(self, n):
        """Gets to a certain position in the BED file when iterating."""
        if 0 <= n < len(self._bed):
            self._n = n

        else:
            # Invalid seek value
            raise ValueError("invalid position in BED: {}".format(n))

    def _read_bim(self):
        """Reads the BIM file."""
        # The original BIM columns
        original_bim_cols = ["chrom", "snp", "cm", "pos", "a1", "a2"]

        # Reading the BIM file and setting the values
        bim = pd.read_csv(self.bim_filename, sep="\t",
                          names=original_bim_cols)
        bim = bim.set_index("snp", drop=False)
        bim["i"] = range(len(bim))
        bim[2] = bim.a1 * 2           # Original '0'
        bim[1] = bim.a1 + bim.a2      # Original '2'
        bim[0] = bim.a2 * 2           # Original '3'
        bim[-1] = "00"                # Original 1

        # Testing something
        allele_encoding = np.array(
            [bim[0], bim[1], bim[2], bim[-1]],
            dtype="U2",
        )
        self._allele_encoding = allele_encoding.T

        # Saving the data in the object
        self._bim = bim[["chrom", "pos", "cm", "a1", "a2", "i"]]
        self._nb_markers = len(self._bim)
        self._markers = np.array(list(self._bim.index))

    def get_bim(self):
        """Returns the BIM file."""
        return self._bim.drop("i", axis=1)

    def get_nb_markers(self):
        """Returns the number of markers."""
        return self._nb_markers

    def _read_fam(self):
        """Reads the FAM file."""
        # The original FAM columns
        self.original_fam_cols = ["fid", "iid", "father", "mother", "gender",
                                  "status"]
        # Reading the FAM file and setting the values
        fam = pd.read_csv(self.fam_filename, sep=" ",
                          names=self.original_fam_cols)
        fam["byte"] = [
            int(np.ceil((1 + 1) / 4.0)) - 1 for i in range(len(fam))
        ]
        fam["bit"] = [(i % 4) * 2 for i in range(len(fam))]

        # Saving the data in the object
        self._fam = fam
        self._nb_samples = len(self._fam)

    def get_fam(self):
        """Returns the FAM file."""
        return self._fam.drop(["byte", "bit"], axis=1)

    def get_nb_samples(self):
        """Returns the number of samples."""
        return self._nb_samples

    def _read_bed(self):
        """Reads the BED file."""
        # Checking if BIM and BAM files were both read
        if (self._bim is None) or (self._fam is None):
            raise RuntimeError("no BIM or FAM file were read")

        data = None
        with open(self.bed_filename, "rb") as bed_file:
            # Checking that the first two bytes are OK
            if (ord(bed_file.read(1)) != 108) or (ord(bed_file.read(1)) != 27):
                raise ValueError("not a valid BED file: "
                                 "{}".format(self.bed_filename))

            # Checking that the format is SNP-major
            if ord(bed_file.read(1)) != 1:
                raise ValueError("not in SNP-major format: "
                                 "{}".format(self.bed_filename))

            # The number of bytes per marker
            nb_bytes = int(np.ceil(self._nb_samples / 4.0))

            # Reading the data
            data = np.fromfile(bed_file, dtype=np.uint8)
            data.shape = (self._nb_markers, nb_bytes)

        # Saving the data in the object
        self._bed = data

    def iter_geno(self):
        """Iterates over genotypes."""
        for i in range(len(self._bed)):
            geno = self._geno_values[self._bed[i]].flatten(order="C")
            yield self._markers[i], geno[:self._nb_samples]

    def iter_acgt_geno(self):
        """Iterates over genotypes (ACGT format)."""
        for i in range(len(self._bed)):
            geno = self._geno_values[self._bed[i]].flatten(order="C")
            yield (self._markers[i],
                   self._allele_encoding[i][geno[:self._nb_samples]])

    def iter_geno_marker(self, markers):
        """Iterates over genotypes for given markers."""
        # If string, we change to list
        if isinstance(markers, str):
            markers = [markers]

        # Checking the list of required markers
        unknown_markers = set(markers) - set(self._bim.index)
        if len(unknown_markers) > 0:
            raise KeyError("marker not in BIM: {}".format(
                sorted(unknown_markers)
            ))

        # Getting the required markers
        required_markers = self._bim.loc[markers, :]

        # Then, we iterate
        for snp, i in required_markers.i.iteritems():
            geno = self._geno_values[self._bed[i]].flatten(order="C")
            yield snp, geno[:self._nb_samples]

    def iter_acgt_geno_marker(self, markers):
        """Iterates over genotypes for given markers (ACGT format)."""
        # If string, we change to list
        if isinstance(markers, str):
            markers = [markers]

        # Checking the list of required markers
        unknown_markers = set(markers) - set(self._bim.index)
        if len(unknown_markers) > 0:
            raise KeyError("marker not in BIM: {}".format(
                sorted(unknown_markers)
            ))

        # Getting the required markers
        required_markers = self._bim.loc[markers, :]

        # Then, we iterate
        for snp, i in required_markers.i.iteritems():
            geno = self._geno_values[self._bed[i]].flatten(order="C")
            yield snp, self._allele_encoding[i][geno[:self._nb_samples]]
