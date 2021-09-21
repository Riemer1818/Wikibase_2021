#!/usr/bin/env python
# -*- coding: utf-8 -*-

from bs4 import BeautifulSoup
from Bio.Restriction import EcoRI, NheI, XbaI, SpeI, PstI, BglII, BamHI, XhoI, AgeI, AarI, RestrictionBatch
import os
import logging
from BB_parser_functions import *
from Diamondblast_functions import *
from WDI_writer import *
from WDI_writer_functions import prepare
import pickle
import sys

__author__ = "Riemer van der Vliet"
__copyright__ = "Copyright 2020, Laboratory of Systems and Synthetic Biology"
__credits__ = ["Riemer van der Vliet", "Jasper Koehorst"]
__license__ = "GPL"
__version__ = "2.0.0"
__maintainer__ = "Riemer van der Vliet"
__email__ = "riemer.vandervliet@wur.nl"
__status__ = "Development"

"""
The Main script orchestrates the entirety of the uploading. Within the script it has two possible paths
to run, the "new" pathway and the "old" pathway. main_old(), uses prior created files to upload to the Wikibase
and the main_new(), deletes old files and creates new files to upload to the Wikibase. To run the script, provide the 
following system arguments: "old/new username password"
"""

sys.setrecursionlimit(10000)


# defines paths to objects and directories
input_path = '../Parts/xml_parts.txt'
blasted_file = '../Parts/Db_output.xml'
fasta_loc = '../Parts/fastafile.fna'
T_directory = "../Parts/Temp_pickle/"
F_directory = "../Parts/Final_pickle/"

# Path to database .dmnd file
database = '/nvme1/riemer/uniprot/UniProt_TREMBL_2020-06.dmnd'

# Wikibase SPARQL endpoint
endpoint_url = "https://bioparts.wiki.opencura.com/query/sparql?"

# UniProt SPARQL endpoint
Sparql_endpoint = 'http://purl.uniprot.org/uniprot/'

# Wikibase API
mediawiki_api_url = "https://bioparts.wiki.opencura.com/w/api.php"

# iGEM HTML sequence
iGEM_sequence_url = "http://parts.iGEM.org/cgi/partsdb/composite_edit/putseq.cgi?part="

logging.basicConfig(level=logging.INFO)

# Unwanted keys to be removed
BB_unwanted = ["seq_edit_cache", "p_status_cache", "s_status_cache", "sequence_update", "review_result",
               "review_count", "review_total", "flag", "rating", "notes", "ok", "has_barcode", "temp_1", "temp_2",
               "temp_3", "temp4", "ps_string", "favorite", "works", "doc_size", "uses", "m_user_id", "m_datetime",
               "sample_status", "discontinued", "informational", "dominant", "sequence_sha1", "categories",
               "group_u_list", "deep_count", "scars", "default_scars", "owner_id", "dominant", "informational",
               "sample_status",
               "owning_group_id"]

BL_unwanted = ["Hit_num", "Hit_hsps", "Hsp_num", "Hsp_query-from", "Hsp_query-to", "Hsp_hit-from", "Hsp_hit-to",
               "Hsp_query-frame",
               "Hsp_hit-frame", "Hsp_identity", "Hsp_positive", "Hsp_gaps", "Hsp_qseq", "Hsp_hseq", "Hsp_midline",
               "Hsp_align-len",
               "Hsp_score", "Hit_len", "Hsp_score"]
# Restriction site objects from Biopython
RS = [EcoRI, XbaI, SpeI, PstI, NheI, BglII, BamHI, XhoI, PstI, AarI]

# Items used by WDI writer
items = ["UniProt", "TrEMBL", "iGEM Parts Registry", "RFC10", "RFC12", "RFC21", "RFC23", "RFC25", "EcoRI", "XbaI",
         "SpeI", "PstI", "NheI", "BglII", "BamHI", "XhoI", "AgeI", "AarI", "The DIAMOND protein aligner", "Coding",
         "Intermediate", "Regulatory", "Generator", "Plasmid", "Composite", "RNA", "RBS", "Plasmid_Backbone",
         "Biobrick", "Biopython", "Planning", "Unavailable", "Available", "Deleted", "Reporter", "DNA", "Terminator",
         "Inverter", "Project", "Measurement", "Device", "Signalling", "Translational_Unit", "Primer", "Temporary",
         "Protein_Domain", "Other"
         ]

# Opens input path
input_file = open(input_path, 'r', encoding="utf8", errors="ignore")

def BB_parser(input_file, BB_unwanted: list):
    """
    Parses the input file and runs the WDI for each of the different biobricks.

    :param input_file: initial XML file
    :param BB_unwanted: Unwanted items to be removed
    """

    # creates Bs4 object
    soup = BeautifulSoup(input_file, "xml")

    # parses Bs4 object
    BB_unparsed = soup.find_all('row')

    # iterates through object
    for children in BB_unparsed:
        # creates dictionary using BB_func()
        BB_dict = BB_func(children, BB_unwanted)

        # filters out empty dictionary
        if BB_dict is not None:
            BB_int_prepare(BB_dict)


def BB_int_prepare(BB_dict: dict):
    """
    Parses the children of the main file into dictionaries and stores these. Also makes a fastafile
    where it adds the sequences to be blasted.

    :param BB_dict: dictionary containing data extracted from XML child.
    """

    # placeholder for dictionary
    WDI_dict = {}

    # names biobrick item
    BB_name = BB_dict['part_name']

    # iterates items in dictionary and performes functions if needed\
    # creates WDI_dict from previous dictionary
    for key, value in BB_dict.items():

        # part_name
        if key == "part_name":
            WDI_dict['part name'] = BB_name

        # sequence, assembly compatability and restriction sites.
        elif key == "sequence":

            # test is a sequence is present
            if len(value) <= 1:
                WDI_dict['RS_Dict'] = None
                WDI_dict['AS_Dict'] = None
            else:
                # adds sequence to fasta file
                create_fastafile(BB_name, value, fasta_loc)

                # adds restriction sites to nested dictionary
                WDI_dict["RS_dict"] = RS_finder(value, RS)

                # iterates dictionary of assembly standards and adds to dictionary
                for key2, value2 in AS_finder(value).items():
                    WDI_dict[key2] = value2

        # long_description
        elif key == "description":
            WDI_dict['Long description'] = HTML_strip(value)

        # description
        elif key == "short_desc":
            WDI_dict['description'] = HTML_strip(value)

        # authors
        elif key == "author":
            authors = []
            # seperates author names from eachother
            authors = Author_name_sep(value)
            WDI_dict['authors'] = authors

        # part_type
        elif key == "part_type":

            # items are defined list, value must be in this list.
            if value in items:
                WDI_dict['part type'] = value
            else:
                logging.info("part type not found " + value)

        # anything else
        else:
            WDI_dict[key] = value

    # filters empty dictionary, otherwise creates pickle file
    if len(WDI_dict) <= 1:
        pass
    else:
        Part_pickle(WDI_dict=WDI_dict, directory=T_directory)


def Part_pickle(WDI_dict: dict, directory: str):
    """Makes pickle file for dictionary

    :param WDI_dict: dictionary containing information of biobrick
    :param directory: directory location where pickle needs to be stored
    """
    with open(directory + WDI_dict['part name'] + ".pickle", 'wb') as handle:
        pickle.dump(WDI_dict, handle, protocol=pickle.DEFAULT_PROTOCOL)


def WDI_dict_blast_add(WDI_dict: dict, BL_unparsed) -> dict:
    """Adds ID to WDI_dict from unparsed blast file

    :param WDI_dict: dictionary containing information of biobrick
    :param BL_unparsed: unparsed Bs4 object from BLAST XML file

    :return: WDI_dict containing information of biobrick
    """

    # placeholder for dictionary
    ID = {}

    # iterates through parsed blast file in dictionary format.
    for key, value in blast_BB_parser(WDI_dict['part name'], unparsed=BL_unparsed, unwanted=BL_unwanted).items():

        # currently sets Hit_nmr to one
        ID["Hit_nmr"] = 1
        ID[key] = value

        # retrieves additional information
        if key == "Hit_accession":

            # location is URL
            loc = Sparql_endpoint + value
            logging.info("retrieving uniprot info on " + loc)

            # iterates through SPARQL results dictionary and adds to ID dictionary
            for key2, value2 in SPARQLWrapper_IDs(loc).items():
                ID[key2] = value2

            # retrieves SPARQL results for EC number
            EC = SPARQLWrapper_EC(loc)

            # checks if EC is present and adds to ID dictionary
            if EC != None:
                ID["EC number"] = EC
            else:
                pass

        elif key == "Hit_def":
            # regular expressions on hit names
            ID['organism'] = organism_sep(value)
            ID['uniprot name'] = uniprot_name_sep(value)

    # adds dictionary to WDI_dict
    if len(ID) <= 1:
        logging.info("No blast hit found for " + WDI_dict["part name"])
    else:
        WDI_dict["UniProt dict"] = ID
    return WDI_dict


def main_old():
    """performs the uploading of pickle items using previous pickle files
    """

    # logs in and retrieves items and property IDs
    login_instance = wdi_login.WDLogin(user=username, pwd=password,
                                       mediawiki_api_url=mediawiki_api_url)
    [item_lookup, property_lookup] = prepare(items, endpoint_url)

    # iterates the files in the final directory, uploads these using WDI_writer function
    for filename in os.listdir(F_directory):
        if filename.endswith(".pickle"):
            with open(F_directory + filename, "rb") as handle:
                WDI_dict = pickle.load(handle)

                WDI_writer(WDI_dict, item_lookup, property_lookup, login_instance, endpoint_url, mediawiki_api_url)
        else:
            continue


def main_new():
    """Makes new files for Final_pickle, performing the blast, parsing the BB dictionary, querying uniprot etc.
    """
    logging.warning("Deleting old files and making new files")

    # removes old files.
    for filename in os.listdir(F_directory):
        if filename.endswith(".pickle"):
            os.remove(F_directory + filename)
        else:
            pass

    # logs in and retrieves items and property IDs
    login_instance = wdi_login.WDLogin(user=username, pwd=password,
                                       mediawiki_api_url=mediawiki_api_url)
    [item_lookup, property_lookup] = prepare(items, endpoint_url)

    # runs the parser
    BB_parser(input_file, BB_unwanted)

    # checking if blastfile is present, if present does not perform BLAST
    if os.path.isfile(blasted_file) != True:
        blast(database, fasta_loc, blasted_file)
    else:
        logging.warning("NOT performING BLAST, REMOVE BLAST FILE")

    # parses the blastfile
    blasted_file_handle = open(blasted_file, encoding="utf8", errors="ignore")
    BL_unparsed = blast_parser(blasted_file_handle)

    # iterates through the temporary pickle files. Makes final pickle files and writes updates or creates new files
    for filename in os.listdir(T_directory):
        if filename.endswith(".pickle"):
            with open(T_directory + filename, "rb") as handle:

                WDI_dict_V1 = pickle.load(handle)

                # if values are too long, they are popped
                for key, value in WDI_dict_V1.copy().items():
                    if value is None:
                        WDI_dict_V1.pop(key)
                for key, value in WDI_dict_V1.copy().items():
                    if len(str(value)) >= 250:
                        WDI_dict_V1.pop(key)

                # adds information retrieved from BLAST
                WDI_dict_V2 = WDI_dict_blast_add(WDI_dict=WDI_dict_V1, BL_unparsed=BL_unparsed)

                # pickles file as final pickle
                Part_pickle(WDI_dict=WDI_dict_V2, directory=F_directory)

                # uploads final pickle files
                WDI_writer(WDI_dict_V2, item_lookup, property_lookup, login_instance, endpoint_url, mediawiki_api_url)

                # closes file
                handle.close()

                # removes temporary pickle file
                os.remove(T_directory + filename)
        else:
            continue


if __name__ == '__main__':
    """Contains two functions and ask which one to use via sys.argv. 
    """

    method      = sys.argv[1]
    username    = sys.argv[2]
    password    = sys.argv[3]

    if len(sys.argv) < 3:
        print("""
            wrong arguments provided, please inform if running "new" files are to be used or "old" files.
            [username] [password]
            """)
        quit()

    # checks if the pipeline is supposed to create new files or use the old ones
    if method == "old":
        logging.info("running the pipeline using old files")
        main_old()

    if method == "new":
        logging.info("running the pipeline creating new files")
        main_new()

    # if method is not "old" or "new".
    else:
        logging.info("""
        no correct system argument provided, inform if running the script using "new" or "old" files, not""", method)
        quit()

