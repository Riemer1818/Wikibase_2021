#!/usr/bin/env python
# -*- coding: utf-8 -*-

import datetime
from wikidataintegrator import wdi_core

__author__ = "Riemer van der Vliet"
__copyright__ = "Copyright 2020, Laboratory of Systems and Synthetic Biology"
__credits__ = ["Riemer van der Vliet", "Jasper Koehorst"]
__license__ = "GPL"
__version__ = "2.0.0"
__maintainer__ = "Riemer van der Vliet"
__email__ = "riemer.vandervliet@wur.nl"
__status__ = "Development"

"""
functions called by WDI_writer script. performs uploading of dictionary items to wikibase.
"""


def create_iGEM_reference(value: str, item_lookup: dict, property_lookup: dict) -> list:
    """Creates iGEM references to the original iGEM Registry of Standard Biological parts.
    [http://parts.igem.org/Main_Page]

    :param value: Biobrick label
    :param item_lookup: Wikibase item IDs.
    :param property_lookup: Wikibase property IDs.

    return: iGEM references as WDI objects
    """

    iGEM_part_url = "http://parts.iGEM.org/Part:"

    # iGEM registry reference
    ref_stated_in = wdi_core.WDItemID(
        prop_nr=property_lookup['stated in'],
        value=item_lookup['iGEM Parts Registry'],
        is_reference=True)

    # biobrick reference
    biobrick_url = iGEM_part_url + value
    ref_reference_biobrick_url = wdi_core.WDUrl(
        prop_nr=property_lookup['reference URL'],
        value=biobrick_url,
        is_reference=True)

    return [ref_stated_in, ref_reference_biobrick_url]


def create_blast_reference(item_lookup: dict, property_lookup: dict) -> list:
    """Creates Blast references to the bioparts wikidata item of the database used and also to the Uniprot database
    wikidata item.

    :param item_lookup: Wikibase item IDs.
    :param property_lookup: Wikibase property IDs.

    return: BLAST references as WDI objects
    """

    # TrEMBL reference
    ref_stated_in_database = wdi_core.WDItemID(
        value=item_lookup['TrEMBL'],
        prop_nr=property_lookup['stated in'],
        is_reference=True)

    # DIAMOND reference
    ref_computation_inference = wdi_core.WDItemID(
        value=item_lookup['The DIAMOND protein aligner'],
        prop_nr=property_lookup['computational inference'],
        is_reference=True)

    return [ref_stated_in_database, ref_computation_inference]


def create_biopython_reference(item_lookup: dict, property_lookup: dict) -> list:
    """Creates a reference to the Biopython python module.

    :param item_lookup: Wikibase item IDs.
    :param property_lookup: Wikibase property IDs.

    return: Biopython references as WDI objects
    """

    # Biopython reference
    ref_computation_inference = wdi_core.WDItemID(
        value=item_lookup['Biopython'],
        prop_nr=property_lookup['computational inference'],
        is_reference=True)

    return [ref_computation_inference]


def create_RS_qualifier(value2: list, property_lookup: dict) -> list:
    """Creates the qualifier to annotate where the restriction sites has been located.

    :param value2: list of restriction sites strings
    :param property_lookup: dictionary containing wikibase property IDs.

    :return: RS qualifiers as WDI object
    """

    # placeholder list
    qualifier_RS = []

    # iterates through list and makes qualifiers
    for RS_hit in value2:
        qualifier_RS.append(wdi_core.WDString(
            value=str(RS_hit),
            prop_nr=property_lookup['at site'],
            is_qualifier=True))

    return qualifier_RS


def create_datetime_qualifier(property_lookup: dict):
    """Creates the date time qualifier.

    :param property_lookup: Wikibase property IDs.

    :return: WDI object
    """
    local_date_time = datetime.datetime.now().strftime("+%Y-%m-%dT00:00:00Z")

    qual_retrieved = wdi_core.WDTime(
        local_date_time,
        prop_nr=property_lookup['retrieved'],
        is_qualifier=True)

    return qual_retrieved
