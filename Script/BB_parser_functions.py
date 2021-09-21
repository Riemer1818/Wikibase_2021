#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA as Nor_DNA
import re
from bs4 import BeautifulSoup
from Bio.Restriction import EcoRI, NheI, XbaI, SpeI, PstI, BglII, BamHI, XhoI, AgeI, AarI, RestrictionBatch
from Bio.Seq import Seq
import logging
from SPARQLWrapper import SPARQLWrapper, JSON

__author__ = "Riemer van der Vliet"
__copyright__ = "Copyright 2020, Laboratory of Systems and Synthetic Biology"
__credits__ = ["Riemer van der Vliet", "Jasper Koehorst"]
__license__ = "GPL"
__version__ = "2.0.0"
__maintainer__ = "Riemer van der Vliet"
__email__ = "riemer.vandervliet@wur.nl"
__status__ = "Development"

"""
Functions used in the MAIN script related to biobrick parsing. 
"""

logging.basicConfig(level=logging.INFO)


def RS_finder(value: str, RS: list) -> dict:
    """Finds just the restriction sites and also where these are if present. Returns a dictionary with lists

    :param value: nucleotide sequence
    :param RS: list of restriction sites objects

    :return: RS_dict with present restriction sites as key and list of sites as value.
    """
    RS_dict = {}
    my_seq = Seq((value), Nor_DNA())

    # iterates through the list of restriction sites.
    try:
        for i in RS:
            RS_dict[i] = i.search(my_seq)
    except TypeError:
        pass

    for key in RS_dict.copy():
        if len(RS_dict[key]) <= 0:
            RS_dict.pop(key)

    return RS_dict


def create_fastafile(BB_name: str, value: str, fasta_loc: str):
    """Creates a fasta file of all the different sequences. Also makes sure that the fasta format is complied with.

    :param BB_name: biobrick part name
    :param value: nucleotide sequence
    :fasta_loc: location of fasta file
    """

    # appends to fasta file in FASTA format https://en.wikipedia.org/wiki/FASTA_format
    fasta = open(fasta_loc, 'a+')
    to_be_written = '>' + BB_name + '\n' + value + '\n'
    fasta.write(to_be_written)
    fasta.close()


def BB_func(children: dict, BB_unwanted: list) -> dict or bool:
    """Does all the parsing and premature edits of the dictionary

    :param children: Bs4 object from XML file
    :param BB_unwanted: Unwanted items from XML file

    :return: BB_dict dictionary of statements extracted from XML or None
    """

    # finds all child of children statment
    child = children.find_next()
    BB = child.find_next_siblings()

    # placeholder for dictionary
    BB_dict = {}

    # fills the dictionary
    for G_child in BB:
        key = G_child.get("name")
        BB_dict.setdefault(key)
        BB_dict[key] = G_child.get_text()

    # adds BB_part_id to dictionary
    BB_part_id = children.find(attrs={"name": "part_id"}).get_text()
    BB_dict["part_id"] = BB_part_id

    # filters out the incorrectly formatted biobricks
    try:
        x = BB_dict["part_name"]
    except KeyError:
        logging.info("skipped {part_id} " + BB_dict["part_id"] + " for now")
        return None

    # filters out the unwanted aspects and the empty values
    for key, value in BB_dict.items():
        if key in BB_unwanted or value is None:
            BB_dict[key] = []
    for key in BB_dict.copy():
        if not BB_dict[key]:
            BB_dict.pop(key)

    return BB_dict


def AS_finder(value: str) -> dict:
    """Finds the restriction sites and assembly compatibilities of common restriction sites in iGEM.

    :param value: nucleotide sequence

    :return: dictionary with key compatible or incompatible and list of assembly compatibilities as value.
    """

    # set assembly standards and sources linking.
    RFC10 = RestrictionBatch([EcoRI, XbaI, SpeI, PstI])
    # http://dspace.mit.edu/handle/1721.1/45138
    RFC12 = RestrictionBatch([EcoRI, SpeI, NheI, PstI])
    # https://dspace.mit.edu/bitstream/handle/1721.1/45139/BBFRFC12.txt?sequence=1&isAllowed=y
    RFC21 = RestrictionBatch([EcoRI, BglII, BamHI, XhoI])
    # https://dspace.mit.edu/bitstream/handle/1721.1/46747/BBFRFC21.pdf?sequence=1&isAllowed=y
    RFC23 = RestrictionBatch([EcoRI, XbaI, SpeI, PstI])
    # https://dspace.mit.edu/bitstream/handle/1721.1/32535/PhillipsSilverFusion.pdf?sequence=1&isAllowed=y
    RFC25 = RestrictionBatch([EcoRI, XbaI, AgeI, SpeI, PstI])
    # https://dspace.mit.edu/bitstream/handle/1721.1/45140/BBF_RFC%2025.pdf?sequence=1&isAllowed=y

    # creates list of assemblies
    assembly_list = [RFC10, RFC12, RFC21, RFC23, RFC25]
    assembly_list_2 = ['RFC10', 'RFC12', 'RFC21', 'RFC23', 'RFC25']

    # helps iteration
    i = 0

    # creates sequence object
    my_seq = Seq(value, Nor_DNA())

    # placeholder for list objects
    lst_incom = []
    lst_comp = []

    # iterates through assembly list
    while i < len(assembly_list):

        # searches for each of the restriciton batches
        RFC = assembly_list[i]
        RFC_AS = RFC.search(my_seq)

        # adds restriction compatability to value list
        for RS_key, RS_value in RFC_AS.items():
            if len(RS_value) > 0:
                lst_incom.append(assembly_list_2[i])

        if assembly_list_2[i] not in lst_incom:
            lst_comp.append(assembly_list_2[i])
        i = i + 1

    # removes duplicates
    lst_incom = list(dict.fromkeys(lst_incom))

    if len(lst_comp) < 1:
        lst_comp = None
    if len(lst_incom) < 1:
        lst_incom = None

    # creates final assembly compatibility dictionary
    AS_dict = {"Incompatible": lst_incom, "Compatible": lst_comp}

    return AS_dict


def HTML_strip(value: str) -> str:
    """Regular expression to find anything that looks like HTML and remove this.

    :param value: text with HTML to be removed

    :return: string with HTML removed
    """

    data = str(value)
    re_html = re.compile(r'<.*?>')
    return re_html.sub('', data)


def Author_name_sep(value: str) -> list:
    """Splits the authors up into a list of different authors. checks on commas and 'and'.

    :param value: string of author text

    :return: list of authors
    """

    x = re.split(", +|  +and +", value)
    return x


def blast_parser(blasted_file):
    """Creates a soup item for the blasted_file and finds also the "iterations" in here.

    :param blasted_file: Blasted file in xml format

    return: unparsed Bs4 object of "Iteration"
    """

    # creates Bs4 item from xml file
    soup = BeautifulSoup(blasted_file, 'xml')

    # returns unparsed children of Bs4 file
    unparsed = soup.find_all('Iteration')
    return unparsed


def blast_BB_parser(BB_name: str, unparsed, unwanted: list) -> dict:
    """Creates the BL_dict for each of the biobricks and removes unwanted items

    :param BB_name: biobrick  part name
    :param unparsed: Bs4 object of "iteration"
    :param unwanted: list of unwanted items

    :return: BL_dict containing statements and items extracted from blast hit
    """
    # placeholder for dictionary
    BL_dict = {}

    # tries to iterate the Bs4 object
    try:
        for children in unparsed:
            if children.find(string=BB_name) is not None:
                # extracts hits and iterates them
                hits = children.find('Hit')
                for line in hits:
                    if line.name is not None:
                        BL_dict[line.name] = line.get_text()
                        if line.name == "Hit_hsps":

                            # extracts children of Hsps to extract text.
                            Hsps = line.children
                            for Hsp in Hsps:
                                if Hsp.name is not None:
                                    CHsp = (Hsp.children)
                                    for child in CHsp:
                                        if child.name is not None:
                                            BL_dict[child.name] = child.text
                # creates copy of dictionary and removes unwanted items in unwanted
                for key in BL_dict.copy():
                    if key in unwanted:
                        BL_dict.pop(key)
    except TypeError:
        pass

    return BL_dict


def SPARQLWrapper_EC(loc: str) -> str:
    """performs the query to the UniProt endpoint.

    :param loc: Location in UniProt (URL)

    :return EC: string
    """

    # endpoint URL
    sparql = SPARQLWrapper("https://sparql.uniprot.org/sparql/")

    # sets query
    query = """
     PREFIX core:<http://purl.uniprot.org/core/>
     SELECT ?enzyme
     WHERE {
         <""" + loc + """> core:enzyme ?enzyme .
     }"""

    # performs query
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    # parses and iterates JSON
    for result in results["results"]["bindings"]:
        EC = result["enzyme"]["value"][31:]

        return EC


def SPARQLWrapper_IDs(loc: str) -> dict:
    """Finds identifiers and values from the UniProt database.

    :param loc: location in UniProt (URL)

    :return: ID keys and ID values
    """

    # endpoint URL
    sparql = SPARQLWrapper("https://sparql.uniprot.org/sparql/")

    # placeholder for dictionary
    ID = {}

    # sets query
    query = """
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    
    SELECT ?enzyme ?ID
    WHERE {
    <""" + loc + """> rdfs:seeAlso ?ID
    }
    """

    # performs query
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    # parses and iterates JSON
    for result in results["results"]["bindings"]:
        key_value = result["ID"]["value"][24:]
        result_ls = re.split("/", key_value)

        # fills dictionary and eliminates duplicates
        for index, object in enumerate(result_ls):
            if index == 0:
                key = object
            elif index == 1:
                value = object
            else:
                pass

        ID[key] = value

    return ID


def organism_sep(value: str) -> str:
    """Separates the blast hit and isolates the organism name.

    :param value: blast hit

    :return: name of organism
    """

    # sets regular expression
    x = re.findall("OS=.*OX", value)

    # iterates list
    for string in x:
        o = string.strip(" OX")
        organism = o.strip("OS=")
        return organism


def uniprot_name_sep(value: str) -> str:
    """Separates the blast hit and isolates the UniProt protein name.

    :param value: blast hit

    :return: UniProt name
    """
    x = re.findall(".*OS=", value)
    for string in x:
        name = string.strip(" OS=")
        return name
