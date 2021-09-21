#!/usr/bin/env python
# -*- coding: utf-8 -*-

from wikidataintegrator import wdi_core, wdi_login
from WDI_value_functions import *
from Parser_functions_V1 import *
from Diamondblast_V2 import *
from WDI_functions_V2 import prepare, get_item_by_name
import copy
import pickle

__author__ = "Riemer van der Vliet"
__copyright__ = "Copyright 2020, Laboratory of Systems and Synthetic Biology"
__credits__ = ["Riemer van der Vliet", "Jasper Koehorst"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Riemer van der Vliet"
__email__ = "riemer.vandervliet@wur.nl"
__status__ = "Development"

"""
The Add_assembly.py file is used to create linkage between BioParts item pages. The script adds the deep_u_list assembly
to the Wikibase for each of the pickled dictionaries in the final directory. It uses the "contains" property on the 
Wikibase to connect the files. 
"""

F_directory = "../Parts/Final_pickle/"
F2_directory = "../Parts/Final_pickle2/"
endpoint_url = "https://bioparts.wiki.opencura.com/query/sparql"
mediawiki_api_url = "https://bioparts.wiki.opencura.com/w/api.php"

# list of item strings used by get_item_by_name() func.
items = ["iGEM Parts Registry"]


def func(handle: str) -> bool:
    """Uploads assembly (property: "contains") to wikibase.

    :param handle: pickled dictionary
    :return: True or False
    """

    # creates WDI_dict from pickle file
    WDI_dict = pickle.load(handle)
    try:
        string = WDI_dict['deep_u_list']
    except KeyError:
        return True

    # creates deep copy of reference
    iGEM_ref = copy.deepcopy(create_iGEM_reference(WDI_dict['part name'], item_lookup, property_lookup))

    # splits the deep_u_list into list of part_id
    x = re.split("_", string)

    # checks that no instance of self are passed
    if len(x) <= 3:
        logging.info(filename + " contains instance of self, skipping")
        return True
    else:
        pass

    statements = []

    # for part_id in deep list statement is appended
    for part_id in x:
        if len(part_id) <= 1:
            pass
        else:
            y = get_and_parse_results(part_id)
            if y is None:
                return False
            else:
                statements.append(wdi_core.WDItemID(
                    value=y,
                    prop_nr=property_lookup['Contains'],
                    references=[iGEM_ref],
                    qualifiers=[datetime_qual]))

    # finds parts page
    parts_page_identifier = get_item_by_name(WDI_dict['part name'], endpoint_url)

    # writes statements to parts page
    parts_page = wdi_core.WDItemEngine(
        wd_item_id=parts_page_identifier,
        new_item=False,
        data=statements,
        mediawiki_api_url=mediawiki_api_url,
        sparql_endpoint_url=endpoint_url)
    parts_page.write(login_instance)
    logging.info("part " + filename + " assembly is added")
    return True


def make_query(part_id: str) -> str:
    """Creates a SPARQL query

    :param part_id: biobrick identifier
    :returns: SPARQL format query

    """
    query = """PREFIX wd: <http://www.wikidata.org/entity/>
        PREFIX wdt: <http://bioparts.wiki.opencura.com/prop/direct/>
        PREFIX wikibase: <http://wikiba.se/ontology#>
        PREFIX p: <http://bioparts.wiki.opencura.com/prop/>
        PREFIX ps: <http://www.wikidata.org/prop/statement/>
        PREFIX pq: <http://www.wikidata.org/prop/qualifier/>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        PREFIX bd: <http://www.bigdata.com/rdf#>
    
    
        SELECT ?item
        Where {
          ?item wdt:P11 \"""" + part_id + """\"; 
        }"""
    return query


def get_and_parse_results(part_id: str) -> str:
    """performes the query using "make_query(param)" function and WDI package SPARQLwrapper

    :param part_id: biobrick identifier
    :return: wikibase item ID
    """
    results = wdi_core.WDItemEngine.execute_sparql_query(query=make_query(part_id), endpoint=endpoint_url)
    for result in results["results"]["bindings"]:
        return result['item']['value'][41:]


if __name__ == '__main__':

    # request username and password for log in.
    username    = sys.argv[1]
    password    = sys.argv[2]

    if len(sys.argv) < 2:
        print("""
            wrong arguments provided [username] [password] need to be provided
            """)
        quit()

    # logs in to the wikibase instance using prepare function
    # creates a datetime qualifier
    [item_lookup, property_lookup] = prepare(items, endpoint_url)
    login_instance = wdi_login.WDLogin(user=username, pwd=password,
                                       mediawiki_api_url=mediawiki_api_url)
    datetime_qual = copy.deepcopy(create_datetime_qualifier(property_lookup))

    handle = []

    # iterates over the files in F_directory and checks if the actions have been performed.
    for filename in os.listdir(F_directory):
        if filename.endswith(".pickle"):
            with open(F_directory + filename, "rb") as handle:
                completed = func(handle)
            if completed is True:
                os.rename(F_directory + filename, F2_directory + filename)
            else:
                logging.warning("part " + filename + " assembly is NOT added")
                pass
