"""Module that reads binary Plink files."""

__all__ = ["PyPlink"]

import os

import numpy as np
import pandas as pd


# The recoding values
_geno_recode = {1: -1,  # Unknown genotype
                2: 1,   # Heterozygous genotype
                0: 2,   # Homozygous A1
                3: 0}   # Homozygous A2

class PyPlink(object):
    """Reads and store a set of binary Plink files."""

    # The genotypes values
    __geno_values = np.array([[_geno_recode[(i >> j) & 3]
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

        # Setting BIM and FAM to None
        self.__bim = None
        self.__fam = None

        # Reading the input files
        self.read_bim()
        self.read_fam()
        self.read_bed()

        # Where we're at
        self.__n = 0


    def __iter__(self):
        """The __iter__ function."""
        return self


    def __next__(self):
        """The __next__ function."""
        return self.next()


    def next(self):
        """The next function."""
        if self.__n < self.__nb_markers:
            self.__n += 1

            # We want to return information about the marker and the genotypes
            geno = self.__geno_values[self.__bed[self.__n - 1]].flatten(order="C")
            return geno[:self.__nb_samples]
        else:
            raise StopIteration()


    def seek(self, n):
        """Gets to a certain position in the BED file when iterating."""
        if 0 <= n < len(self.__bed):
            self.__n = n

        else:
            # Invalid seek value
            raise ValueError("invalid position in BED: {}".format(n))


    def read_bim(self):
        """Reads the BIM file."""
        # The original BIM columns
        original_bim_cols = ["chr", "snp", "cm", "pos", "a1", "a2"]

        # Reading the BIM file and setting the values
        bim = pd.read_csv(self.bim_filename, sep="\t",
                          names=original_bim_cols)
        bim[2] = bim.a1 * 2           # Original '0'
        bim[1] = bim.a1 + bim.a2      # Original '2'
        bim[0] = bim.a2 * 2           # Original '3'
        bim[-1] = "00"                # Original 1

        # Testing something
        allele_encoding = np.array([bim[0], bim[1], bim[2], bim[-1]], dtype="U2")
        self._allele_encoding = allele_encoding.T

        # Saving the data in the object
        self.__bim = bim[original_bim_cols]
        self.__nb_markers = len(self.__bim)


    def get_bim(self):
        """Returns the BIM file."""
        return self.__bim.copy()


    def get_nb_markers(self):
        """Returns the number of markers."""
        return self.__nb_markers


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
        self.__fam = fam
        self.__nb_samples = len(self.__fam)


    def get_fam(self):
        """Returns the FAM file."""
        return self.__fam.copy()


    def get_nb_samples(self):
        """Returns the number of samples."""
        return self.__nb_samples


    def read_bed(self):
        """Reads the BED file."""
        # Checking if BIM and BAM files were both read
        if (self.__bim is None) or (self.__fam is None):
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
            nb_bytes = int(np.ceil(self.__nb_samples / 4.0))

            # Reading the data
            data = np.fromfile(bed_file, dtype=np.uint8)
            data.shape = (self.__nb_markers, nb_bytes)

        # Saving the data in the object
        self.__bed = data


    def iter_geno(self):
        """Iterates over genotypes."""
        for i in range(len(self.__bed)):
            geno = self.__geno_values[self.__bed[i]].flatten(order="C")
            yield geno[:self.__nb_samples]


    def iter_acgt_geno(self):
        """Iterates over genotypes."""
        for i in range(len(self.__bed)):
            geno = self.__geno_values[self.__bed[i]].flatten(order="C")
            yield self._allele_encoding[i][geno[:self.__nb_samples]]
