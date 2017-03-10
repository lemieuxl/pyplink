.. pyplink documentation master file, created by
   sphinx-quickstart on Fri Mar 10 09:00:31 2017.

pyplink
========

:py:mod:`pyplink` is a Python (2 and 3) binary Plink file parser and writer
released under the MIT licence. The difference with existing parsers (and Plink
itself) is that :py:mod:`pyplink` does not load the BED file in memory, making
possible to work with extremely large files (*e.g.* 1000 Genomes Phase 3
files).

For more information on how to use :py:mod:`pyplink`, refer to the
:doc:`API documentation <pyplink>`. Below is a snippet describing the most
common usage of the module.


.. code-block:: python

    from pyplink import PyPlink

    with PyPlink("plink_file_prefix") as bed:
        # Getting the BIM and FAM
        bim = bed.get_bim()
        fam = bed.get_fam()

        # Iterating over all loci
        for loci_name, genotypes in bed:
            pass

        # Getting the genotypes of a single marker (numpy.ndarray)
        genotypes = bed.get_geno_marker("rs12345")


.. toctree::
   :hidden:

   installation
   pyplink
