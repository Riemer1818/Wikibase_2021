#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from wikidataintegrator import wdi_core, wdi_login
import logging
import pickle

__author__ = "Riemer van der Vliet"
__copyright__ = "Copyright 2020, Laboratory of Systems and Synthetic Biology"
__credits__ = ["Riemer van der Vliet", "Jasper Koehorst"]
__license__ = "GPL"
__version__ = "2.0.0"
__maintainer__ = "Riemer van der Vliet"
__email__ = "riemer.vandervliet@wur.nl"
__status__ = "Development"

"""
functions used by WDI writer file
"""


def get_properties(endpoint_url: str) -> dict:
    """Finds properties on the endpoint url and returns the IDs

    :param endpoint_url: Wikibase SPARQL endpoint

    :return: Property lookup dictionary of key property string and value property ID of Wikibase
    """

    # placeholder for dictionary
    property_lookup = {}

    # creates query
    query = """SELECT ?property ?label WHERE {
        ?property a wikibase:Property .
        ?property rdfs:label ?label .
        FILTER (LANG(?label) = "en" )}
        """

    # gets results
    results = wdi_core.WDItemEngine.execute_sparql_query(query=query, endpoint=endpoint_url)

    # iterates iterates data
    for result in results["results"]["bindings"]:
        label = result["label"]["value"].split("/")[-1]
        property_lookup[label] = result["property"]["value"].split("/")[-1]

    return property_lookup


def get_items(items: list, endpoint_url: str) -> dict:
    """Gets the IDs for each of the items in the item list. First tries to find it in the pickle file.

    :param items: list of items of which IDs need to be traced
    :param endpoint_url: Wikibase SPARQL endpoint

    :return: item_lookup dictionary of key item string and value item ID of Wikibase
    """
    if os.path.isfile("../Parts/item_lookup.pickle"):
        with open('../Parts/item_lookup.pickle', 'rb') as handle:
            item_lookup = pickle.load(handle)
    else:
        item_lookup = {}

    for item_x in items:
        logging.info("Retrieving item " + item_x)
        if item_x in item_lookup: continue
        item_lookup[item_x] = get_item_by_name(item_x, endpoint_url)

    with open('../Parts/item_lookup.pickle', 'wb') as handle:
        pickle.dump(item_lookup, handle, protocol=pickle.DEFAULT_PROTOCOL)

    return item_lookup


def get_item_by_name(label: str, endpoint_url: str) -> str or None:
    """Finds items on the endpoint url and returns the IDs

    :param label: Item label
    :param endpoint_url: Wikibase SPARQL endpoint

    :return: string of Wikibase ID or None
    """

    # set query
    query = """
    SELECT DISTINCT ?item WHERE { 
      VALUES ?label { \"""" + label + """\"@en }
      ?item rdfs:label ?label .
    }"""

    # get results
    try:
        results = wdi_core.WDItemEngine.execute_sparql_query(query, endpoint=endpoint_url)
    except:
        print("Query failed: ")
        raise Exception("Query failed")

    # parse and return results
    for result in results["results"]["bindings"]:
        return result["item"]["value"].split("/")[-1]
    return None


def prepare(items: list, endpoint_url: str) -> list:
    """Returns a list of lists of items ID and property IDs

    :param items: list of items of which IDs need to be traced
    :param endpoint_url: Wikibase SPARQL endpoint

    :return: list of item dictionary and of property dictionary
    """
    return [get_items(items, endpoint_url), get_properties(endpoint_url)]


def get_properties(endpoint_url: str) -> dict:
    """Finds properties on the endpoint url and returns the IDs

    :param endpoint_url: Wikibase SPARQL endpoint

    :return: property_lookup dictionary of key property string and value property ID of Wikibase
    """

    # placeholder for dictionary
    property_lookup = {}

    # set query
    query = """SELECT ?property ?label WHERE {
        ?property a wikibase:Property .
        ?property rdfs:label ?label .
        FILTER (LANG(?label) = "en" )}
        """

    # get results
    results = wdi_core.WDItemEngine.execute_sparql_query(query=query, endpoint=endpoint_url)

    # parse results
    for result in results["results"]["bindings"]:
        label = result["label"]["value"].split("/")[-1]
        property_lookup[label] = result["property"]["value"].split("/")[-1]

    return property_lookup


def get_items(items: list, endpoint_url: str) -> dict:
    """Gets the IDs for each of the items in the item list. First tries to find it in the pickle file.

    :param items: list of items of which IDs need to be traced
    :param endpoint_url: Wikibase SPARQL endpoint

    :return: item_lookup dictionary with item strings and value IDs
    """

    # checks if there is a pickle file under name item_lookup.pickle,
    # otherwise creates dictionary placeholder
    if os.path.isfile("../Parts/item_lookup.pickle"):
        with open('../Parts/item_lookup.pickle', 'rb') as handle:
            item_lookup = pickle.load(handle)
    else:
        item_lookup = {}

    # iterates items and gets the item ID by name
    for item_x in items:
        logging.info("Retrieving item " + item_x)
        if item_x in item_lookup: continue
        item_lookup[item_x] = get_item_by_name(item_x, endpoint_url)

    # dumps object as pickle file
    with open('../Parts/item_lookup.pickle', 'wb') as handle:
        pickle.dump(item_lookup, handle, protocol=pickle.DEFAULT_PROTOCOL)

    return item_lookup


def get_item_by_name(label: str, endpoint_url: str) -> str or bool:
    """Finds items on the endpoint url and returns the IDs

    :param label: Item label
    :param endpoint_url: Wikibase SPARQL endpoint

    :return: result string of wikibase ID or None
    """

    # sets query
    query = """
    SELECT DISTINCT ?item WHERE { 
      VALUES ?label { \"""" + label + """\"@en }
      ?item rdfs:label ?label .
    }"""

    # gets results
    try:
        results = wdi_core.WDItemEngine.execute_sparql_query(query, endpoint=endpoint_url)
    except:
        print("Query failed: ")
        raise Exception("Query failed")

    # iterates results
    for result in results["results"]["bindings"]:
        return result["item"]["value"].split("/")[-1]
    return None
