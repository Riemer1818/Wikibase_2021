#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import logging

__author__ = "Riemer van der Vliet"
__copyright__ = "Copyright 2020, Laboratory of Systems and Synthetic Biology"
__credits__ = ["Riemer van der Vliet", "Jasper Koehorst"]
__license__ = "GPL"
__version__ = "2.0.0"
__maintainer__ = "Riemer van der Vliet"
__email__ = "riemer.vandervliet@wur.nl"
__status__ = "Development"

"""
functions called by MAIN.py script. performs te blast.
"""

logging.basicConfig(level=logging.INFO)


def blast(database: str, fasta_loc: str, blasted_file: str):
    """
    Performs the blast given an '../Parts/input_fasta_file' and deletes this after.

    :param database: location of root database file
    :param fasta_loc: location of fasta file
    :param blasted_blasted file: location of output file (blasted file)
    """

    # performs the DIAMOND blast command. output is set to 5 (XML format), max-target-seqs is set 1 hit.
    os.system(
        'diamond blastx -d ' + database + ' -q ' + fasta_loc + ' -o ' + blasted_file + ' --outfmt 5 --max-target-seqs 1')
    logging.info("done blast")

    # removes fasta file
    os.remove(fasta_loc)


if __name__ == '__main__':
    database = '/nvme1/riemer/uniprot/UniProt_TREMBL_2020-06.dmnd'
    fasta_loc = './fastafile.fna'
    blasted_file = '../Parts/Db_output.xml'
    blast(database, fasta_loc, blasted_file)
