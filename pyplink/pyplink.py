"""Module that reads binary Plink files."""

__all__ = ["PyPlink"]

import os

import numpy as np
import pandas as pd


# The recoding values
_geno_recode = {1: -9,  # Unknown genotype
                2: 1,   # Heterozygous genotype
                0: 2,   # Homozygous A1
                3: 0}   # Homozygous A2

class PyPlink(object):
    """Reads and store a set of binary Plink files."""

    # The genotypes values
    _geno_values = np.array([[_geno_recode[(i >> j) & 3]
                                           for j in range(0, 7, 2)]
                                           for i in range(256)], dtype=np.int8)


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

        # Reading the input files
        self.read_bim()
        self.read_fam()
        self.read_bed()

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
        if self._n < self.nb_markers:
            self._n += 1

            # We want to return information about the marker and the genotypes
            geno = self._geno_values[self._bed[self._n - 1]].flatten(order="C")
            return geno[:self.nb_samples]
        else:
            raise StopIteration()


    def read_bim(self):
        """Reads the BIM file."""
        # The original BIM columns
        self.original_bim_cols = ["chr", "snp", "cm", "pos", "a1", "a2"]

        # Reading the BIM file and setting the values
        bim = pd.read_csv(self.bim_filename, sep="\t",
                          names=self.original_bim_cols)
        bim["2"] = bim.a1 * 2           # Original '0'
        bim["1"] = bim.a1 + bim.a2      # Original '2'
        bim["0"] = bim.a2 * 2           # Original '3'
        bim["-9"] = "00"                # Original 1

        # Saving the data in the object
        self._bim = bim
        self.nb_markers = len(self._bim)


    def get_bim(self):
        """Returns the BIM file."""
        return self._bim


    def read_fam(self):
        """Reads the FAM file."""
        # The original FAM columns
        self.original_fam_cols = ["fid", "iid", "father", "mother", "gender",
                                  "status"]
        # Reading the FAM file and setting the values
        fam = pd.read_csv(self.fam_filename, sep=" ",
                          names=self.original_fam_cols)
        fam["byte"] = [int(np.ceil((1 + 1) / 4.0)) - 1 for i in range(len(fam))]
        fam["bit"] = [(i % 4) * 2 for i in range(len(fam))]

        # Saving the data in the object
        self._fam = fam
        self.nb_samples = len(self._fam)
        self.samples = ["{};{}".format(fid, iid)
                            for fid, iid in zip(self._fam.fid, self._fam.iid)]


    def get_fam(self):
        """Returns the FAM file."""
        return self._fam


    def read_bed(self):
        """Reads the BED file."""
        # Checking if BIM and BAM files were both read
        if (not hasattr(self, "_bim")) or (not hasattr(self, "_fam")):
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
            nb_bytes = int(np.ceil(self.nb_samples / 4.0))

            # Reading the data
            data = np.fromfile(bed_file, dtype=np.uint8)
            data.shape = (self.nb_markers, nb_bytes)

        # Saving the data in the object
        self._bed = data


    def get_bed(self):
        """Returns the BED file."""
        return self._bed
