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
import logging
from itertools import repeat
from collections import Counter
from io import UnsupportedOperation

try:
    from itertools import zip_longest
except ImportError:
    from itertools import izip_longest as zip_longest

import numpy as np
import pandas as pd

from six.moves import range


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014 Louis-Philippe Lemieux Perreault"
__license__ = "MIT"


__all__ = ["PyPlink"]


# The logger
logger = logging.getLogger(__name__)


# The recoding values
_geno_recode = {1: -1,  # Unknown genotype
                2:  1,  # Heterozygous genotype
                0:  2,  # Homozygous A1
                3:  0}  # Homozygous A2
_byte_recode = dict(value[::-1] for value in _geno_recode.items())


class PyPlink(object):
    """Reads and store a set of binary Plink files.

    Args:
        prefix (str): The prefix of the binary Plink files.
        mode (str): The open mode for the binary Plink file.
        bed_format (str): The type of bed (SNP-major or INDIVIDUAL-major).

    Reads or write binary Plink files (BED, BIM and FAM).

    .. code-block:: python

        from pyplink import PyPlink

        # Reading BED files
        with PyPlink("plink_file_prefix") as bed:
            pass

        # Writing BED files
        with PyPlink("plink_file_prefix", "w") as bed:
            pass

    """

    # The genotypes values
    _geno_values = np.array(
        [
            [_geno_recode[(i >> j) & 3] for j in range(0, 7, 2)]
            for i in range(256)
        ],
        dtype=np.int8,
    )

    def __init__(self, prefix, mode="r", bed_format="SNP-major"):
        """Initializes a new PyPlink instance."""
        # The mode
        self._mode = mode

        # The bed format
        if bed_format not in {"SNP-major", "INDIVIDUAL-major"}:
            raise ValueError("invalid bed format: {}".format(bed_format))
        self._bed_format = bed_format

        # These are the name of the files
        self.bed_filename = "{}.bed".format(prefix)
        self.bim_filename = "{}.bim".format(prefix)
        self.fam_filename = "{}.fam".format(prefix)

        if self._mode == "r":
            if self._bed_format != "SNP-major":
                raise ValueError("only SNP-major format is supported "
                                 "with mode 'r'")

            # Checking that all the files exists (otherwise, error...)
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

        elif self._mode == "w":
            # The dummy number of samples and bytes
            self._nb_values = None

            # Opening the output BED file
            self._bed = open(self.bed_filename, "wb")
            self._write_bed_header()

        else:
            raise ValueError("invalid mode: '{}'".format(self._mode))

    def __repr__(self):
        """The representation of the PyPlink object."""
        if self._mode == "r":
            return "PyPlink({:,d} samples; {:,d} markers)".format(
                self.get_nb_samples(),
                self.get_nb_markers(),
            )

        return 'PyPlink(mode="w")'

    def __iter__(self):
        """The __iter__ function."""
        if self._mode != "r":
            raise UnsupportedOperation("not available in 'w' mode")

        return self

    def __next__(self):
        """The __next__ function."""
        return self.next()

    def __enter__(self):
        """Entering the context manager."""
        return self

    def __exit__(self, *args):
        """Exiting the context manager."""
        self.close()

    def close(self):
        """Closes the BED file."""
        # Closing the BED file
        self._bed.close()

    def next(self):
        """Returns the next marker.

        Returns:
            tuple: The marker name as a string and its genotypes as a
            :py:class:`numpy.ndarray`.

        """
        if self._mode != "r":
            raise UnsupportedOperation("not available in 'w' mode")

        self._n += 1
        if self._n > self._nb_markers:
            raise StopIteration()

        return self._bim.index[self._n - 1], self._read_current_marker()

    def _read_current_marker(self):
        """Reads the current marker and returns its genotypes."""
        return self._geno_values[
            np.frombuffer(self._bed.read(self._nb_bytes), dtype=np.uint8)
        ].flatten(order="C")[:self._nb_samples]

    def seek(self, n):
        """Gets to a certain marker position in the BED file.

        Args:
            n (int): The index of the marker to seek to.

        """
        if self._mode != "r":
            raise UnsupportedOperation("not available in 'w' mode")

        if 0 <= n < self._nb_markers:
            self._n = n
            self._bed.seek(self._get_seek_position(n))

        else:
            # Invalid seek value
            raise ValueError("invalid position in BED: {}".format(n))

    def _get_seek_position(self, n):
        """Gets the seek position in the file (including special bytes).

        Args:
            n (int): The index of the marker to seek to.

        """
        return 3 + self._nb_bytes * n

    def _read_bim(self):
        """Reads the BIM file."""
        # Reading the BIM file and setting the values
        bim = pd.read_csv(self.bim_filename, sep=r"\s+",
                          names=["chrom", "snp", "cm", "pos", "a1", "a2"],
                          dtype=dict(snp=str, a1=str, a2=str))

        # Saving the index as integer
        bim["i"] = bim.index

        # Checking for duplicated markers
        try:
            bim = bim.set_index("snp", verify_integrity=True)
            self._has_duplicated = False

        except ValueError:
            # Setting this flag to true
            self._has_duplicated = True

            # Finding the duplicated markers
            duplicated = bim.snp.duplicated(keep=False)
            duplicated_markers = bim.loc[duplicated, "snp"]
            duplicated_marker_counts = duplicated_markers.value_counts()

            # The dictionary that will contain information about the duplicated
            # markers
            self._dup_markers = {
                m: [] for m in duplicated_marker_counts.index
            }

            # Logging a warning
            logger.warning("Duplicated markers found")
            for marker, count in duplicated_marker_counts.items():
                logger.warning("  - {}: {:,d} times".format(marker, count))
            logger.warning("Appending ':dupX' to the duplicated markers "
                           "according to their location in the BIM file")

            # Renaming the markers
            counter = Counter()
            for i, marker in duplicated_markers.items():
                counter[marker] += 1
                new_name = "{}:dup{}".format(marker, counter[marker])
                bim.loc[i, "snp"] = new_name

                # Updating the dictionary containing the duplicated markers
                self._dup_markers[marker].append(new_name)

            # Resetting the index
            bim = bim.set_index("snp", verify_integrity=True)

        # Encoding the allele
        #   - The original 0 is the actual 2 (a1/a1)
        #   - The original 2 is the actual 1 (a1/a2)
        #   - The original 3 is the actual 0 (a2/a2)
        #   - The original 1 is the actual -1 (no call)
        allele_encoding = np.array(
            [bim.a2 * 2, bim.a1 + bim.a2, bim.a1 * 2,
             list(repeat("00", bim.shape[0]))],
            dtype="U2",
        )
        self._allele_encoding = allele_encoding.T

        # Saving the data in the object
        self._bim = bim[["chrom", "pos", "cm", "a1", "a2", "i"]]
        self._nb_markers = self._bim.shape[0]

    def get_bim(self):
        """Returns the BIM file.

        Returns:
            pandas.DataFrame: The BIM file.

        """
        if self._mode != "r":
            raise UnsupportedOperation("not available in 'w' mode")

        return self._bim.drop("i", axis=1)

    def get_nb_markers(self):
        """Returns the number of markers.

        Returns:
            int: The number of markers in the dataset.

        """
        if self._mode != "r":
            raise UnsupportedOperation("not available in 'w' mode")

        return self._nb_markers

    def get_duplicated_markers(self):
        """Returns the duplicated markers, if any.

        Returns:
            dict: The set of duplicated marker (might be empty).

        """
        if self._has_duplicated:
            return self._dup_markers
        else:
            return {}

    def _read_fam(self):
        """Reads the FAM file."""
        # Reading the FAM file and setting the values
        fam = pd.read_csv(self.fam_filename, sep=r"\s+",
                          names=["fid", "iid", "father", "mother", "gender",
                                 "status"],
                          dtype=dict(fid=str, iid=str, father=str, mother=str))

        # Saving the data in the object
        self._fam = fam
        self._nb_samples = self._fam.shape[0]

    def get_fam(self):
        """Returns the FAM file.

        Returns:
            pandas.DataFrame: The FAM file.

        """
        if self._mode != "r":
            raise UnsupportedOperation("not available in 'w' mode")

        return self._fam.copy()

    def get_nb_samples(self):
        """Returns the number of samples.

        Returns:
            int: The number of samples in the dataset.

        """
        if self._mode != "r":
            raise UnsupportedOperation("not available in 'w' mode")

        return self._nb_samples

    def _read_bed(self):
        """Reads the BED file."""
        # Checking if BIM and BAM files were both read
        if (self._bim is None) or (self._fam is None):
            raise RuntimeError("no BIM or FAM file were read")

        # The number of bytes per marker
        self._nb_bytes = int(np.ceil(self._nb_samples / 4.0))

        # Checking the file is valid by looking at the first 3 bytes and the
        # last entry (correct size)
        with open(self.bed_filename, "rb") as bed_file:
            # Checking that the first two bytes are OK
            if (ord(bed_file.read(1)) != 108) or (ord(bed_file.read(1)) != 27):
                raise ValueError("not a valid BED file: "
                                 "{}".format(self.bed_filename))

            # Checking that the format is SNP-major
            if ord(bed_file.read(1)) != 1:
                raise ValueError("not in SNP-major format (please recode): "
                                 "{}".format(self.bed_filename))

            # Checking the last entry (for BED corruption)
            seek_index = self._get_seek_position(self._bim.iloc[-1, :].i)
            bed_file.seek(seek_index)
            geno = self._geno_values[
                np.frombuffer(bed_file.read(self._nb_bytes), dtype=np.uint8)
            ].flatten(order="C")[:self._nb_samples]
            if geno.shape[0] != self._nb_samples:
                raise ValueError("invalid number of entries: corrupted BED?")

        # Opening the file for the rest of the operations (reading 3 bytes)
        self._bed = open(self.bed_filename, "rb")
        self._bed.read(3)

    def _write_bed_header(self):
        """Writes the BED first 3 bytes."""
        # Writing the first three bytes
        final_byte = 1 if self._bed_format == "SNP-major" else 0
        self._bed.write(bytearray((108, 27, final_byte)))

    def iter_geno(self):
        """Iterates over genotypes from the beginning of the BED file.

        Returns:
            tuple: The name of the marker as a string, and its genotypes as a
            :py:class:`numpy.ndarray` (additive format).

        """
        if self._mode != "r":
            raise UnsupportedOperation("not available in 'w' mode")

        # Seeking back at the beginning of the file
        self.seek(0)

        # Return itself (the generator)
        return self

    def iter_acgt_geno(self):
        """Iterates over genotypes (ACGT format).

        Returns:
            tuple: The name of the marker as a string, and its genotypes as a
            :py:class:`numpy.ndarray` (ACGT format).

        """
        # Need to iterate over itself, and modify the actual genotypes
        for i, (marker, geno) in enumerate(self.iter_geno()):
            yield marker, self._allele_encoding[i][geno]

    def iter_geno_marker(self, markers, return_index=False):
        """Iterates over genotypes for a list of markers.

        Args:
            markers (list): The list of markers to iterate onto.
            return_index (bool): Wether to return the marker's index or not.

        Returns:
            tuple: The name of the marker as a string, and its genotypes as a
            :py:class:`numpy.ndarray` (additive format).

        """
        if self._mode != "r":
            raise UnsupportedOperation("not available in 'w' mode")

        # If string, we change to list
        if isinstance(markers, str):
            markers = [markers]

        # Iterating over all markers
        if return_index:
            for marker in markers:
                geno, seek = self.get_geno_marker(marker, return_index=True)
                yield marker, geno, seek
        else:
            for marker in markers:
                yield marker, self.get_geno_marker(marker)

    def iter_acgt_geno_marker(self, markers):
        """Iterates over genotypes for a list of markers (ACGT format).

        Args:
            markers (list): The list of markers to iterate onto.

        Returns:
            tuple: The name of the marker as a string, and its genotypes as a
            :py:class:`numpy.ndarray` (ACGT format).

        """
        # We iterate over the markers
        for snp, geno, s in self.iter_geno_marker(markers, return_index=True):
            # Getting the SNP position and converting to ACGT
            yield snp, self._allele_encoding[s][geno]

    def get_geno_marker(self, marker, return_index=False):
        """Gets the genotypes for a given marker.

        Args:
            marker (str): The name of the marker.
            return_index (bool): Wether to return the marker's index or not.

        Returns:
            numpy.ndarray: The genotypes of the marker (additive format).

        """
        if self._mode != "r":
            raise UnsupportedOperation("not available in 'w' mode")

        # Check if the marker exists
        if marker not in self._bim.index:
            raise ValueError("{}: marker not in BIM".format(marker))

        # Seeking to the correct position
        seek_index = self._bim.loc[marker, "i"]
        self.seek(seek_index)

        if return_index:
            return self._read_current_marker(), seek_index
        return self._read_current_marker()

    def get_acgt_geno_marker(self, marker):
        """Gets the genotypes for a given marker (ACGT format).

        Args:
            marker (str): The name of the marker.

        Returns:
            numpy.ndarray: The genotypes of the marker (ACGT format).

        """
        # Getting the marker's genotypes
        geno, snp_position = self.get_geno_marker(marker, return_index=True)

        # Returning the ACGT's format of the genotypes
        return self._allele_encoding[snp_position][geno]

    def write_genotypes(self, genotypes):
        """Write genotypes to binary file.

        Args:
            genotypes (numpy.ndarray): The genotypes to write in the BED file.

        """
        if self._mode != "w":
            raise UnsupportedOperation("not available in 'r' mode")

        # Initializing the number of samples if required
        if self._nb_values is None:
            self._nb_values = len(genotypes)

        # Checking the expected number of samples
        if self._nb_values != len(genotypes):
            raise ValueError("{:,d} samples expected, got {:,d}".format(
                self._nb_values,
                len(genotypes),
            ))

        # Writing to file
        byte_array = [
            g[0] | (g[1] << 2) | (g[2] << 4) | (g[3] << 6) for g in
            self._grouper((_byte_recode[geno] for geno in genotypes), 4)
        ]
        self._bed.write(bytearray(byte_array))

    @staticmethod
    def _grouper(iterable, n, fillvalue=0):
        """Collect data into fixed-length chunks or blocks.

        Args:
            n (int): The size of the chunk.
            fillvalue (int): The fill value.

        Returns:
            iterator: An iterator over the chunks.

        """
        args = [iter(iterable)] * n
        return zip_longest(fillvalue=fillvalue, *args)
