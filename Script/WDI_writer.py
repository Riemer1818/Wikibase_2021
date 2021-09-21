#!/usr/bin/env python
# -*- coding: utf-8 -*-

import copy
from wikidataintegrator import wdi_core, wdi_login
from WDI_value_functions import *
from WDI_writer_functions import get_item_by_name
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
functions called by MAIN.py script. performs uploading of dictionary items to wikibase.
"""

logging.basicConfig(level=logging.INFO)


def WDI_writer(WDI_dict: dict, item_lookup: dict, property_lookup: dict,
               login_instance, endpoint_url: str, mediawiki_api_url: str
               ):
    """Creates statements for the iterated dictionary and creates an item page is this is not already present.
    Otherwise updates the item page.

    :param WDI_dict: Containing data to be iterated and uploaded,
    :param item_lookup: Wikibase item IDs.
    :param property_lookup: Wikibase property IDs.
    :param login_instance: login instance of the Wikibase bot.
    :param endpoint_url: SPARQL endpoint of Wikibase.
    :param mediawiki_api_url: API of Wikibase.
    """

    logging.info("-------------------------next biobrick-------------------------")

    # sets label
    label = WDI_dict['part name']
    logging.info("Parsing biobrick " + label)

    # statement holder for the biobrick page
    statements = []
    # holder for aliases
    aliases = []

    # creating deep copies of references and qualifiers.
    iGEM_ref = copy.deepcopy(create_iGEM_reference(WDI_dict['part name'], item_lookup, property_lookup))
    Blast_ref = copy.deepcopy(create_blast_reference(item_lookup, property_lookup))
    Biopython_ref = copy.deepcopy(create_biopython_reference(item_lookup, property_lookup))
    datetime_qual = copy.deepcopy(create_datetime_qualifier(property_lookup))

    # instance of biobrick is a conserved statement
    statements.append(wdi_core.WDItemID(
        value=item_lookup['Biobrick'],
        prop_nr=property_lookup['instance of'],
        references=[iGEM_ref]))

    # parsing the WDI dictionary and iterating and creating statements of each of the items.
    for key, value in WDI_dict.items():

        # description
        if key == 'description':
            logging.info("Parsing description")
            description = value

        # part name
        elif key == 'part name':
            logging.info("Parsing iGem Parts Identifier")
            statements.append(wdi_core.WDExternalID(
                value=WDI_dict[key],
                prop_nr=property_lookup['iGem Parts ID'],
                references=[iGEM_ref]))

        # long description
        elif key == 'long description':
            logging.info("Parsing long description")
            statements.append(wdi_core.WDString(
                value=value,
                prop_nr=property_lookup['long description'],
                references=[iGEM_ref]))

        # skipping sequence
        elif key == 'sequence':
            pass

        # sequence and sequence length (+URL)
        elif key == 'sequence_length':
            logging.info("Parsing sequence")
            qual_list = []
            sl_qual_list = []

            sl_qual_list.append(wdi_core.WDString(
                value=str(value),
                prop_nr=property_lookup['sequence length'],
                is_qualifier=True))

            qual_list = sl_qual_list.append(datetime_qual)

            statements.append(wdi_core.WDString(
                value="http://parts.igem.org/cgi/partsdb/composite_edit/putseq.cgi?part=" + label,
                prop_nr=property_lookup['sequence'],
                references=[iGEM_ref],
                qualifiers=qual_list))

        # authors
        elif key == 'authors':
            logging.info("Parsing author")
            for author in value:
                statements.append(wdi_core.WDString(
                    value=author.strip(),
                    prop_nr=property_lookup['author'],
                    references=[iGEM_ref]))

        # restriction sites
        elif key == 'RS_dict':
            logging.info("Parsing RS_dict")
            for key2, value2 in value.items():
                print(key2, value2)
                RS_qual_list = create_RS_qualifier(value2, property_lookup)
                print(RS_qual_list)

                statements.append(
                    wdi_core.WDItemID(
                        value=item_lookup[str(key2)],
                        prop_nr=property_lookup['restriction site'],
                        references=[iGEM_ref, Biopython_ref],
                        qualifiers=RS_qual_list
                    )
                )

        # compatibilities
        elif key == 'Compatible':
            logging.info("Parsing assembly compatibilities")
            for AS in value:
                item = str(item_lookup[AS])
                statements.append(wdi_core.WDItemID(
                    value=item, prop_nr=property_lookup['Compatible with'],
                    references=[Biopython_ref],
                    qualifiers=[datetime_qual]))

        # incompatibilities
        elif key == 'Incompatible':
            logging.info("Parsing assembly incompatibilities")
            for AS in value:
                item = str(item_lookup[AS])
                statements.append(wdi_core.WDItemID(
                    value=item,
                    prop_nr=property_lookup['Incompatible with'],
                    references=[Biopython_ref],
                    qualifiers=[datetime_qual]))

        # part type
        elif key == 'part type':
            logging.info("Parsing part type " + value)
            statements.append(wdi_core.WDItemID(
                value=item_lookup[value],
                prop_nr=property_lookup['part type'],
                references=[iGEM_ref]))

        # UniProt_hit creates statement with qualifiers
        elif key == 'UniProt dict':
            logging.info("Parsing blast hit")

            qual_list = []

            # iterates the items of the dictionary
            for key2, value2 in value.items():

                # UniProt name
                if key2 == 'uniprot name':
                    qual_list.append(wdi_core.WDExternalID(
                        value=value2,
                        prop_nr=property_lookup['UniProt name'],
                        is_qualifier=True))

                # EC number
                elif key2 == 'EC number':
                    qual_list.append(wdi_core.WDExternalID(
                        value=value2,
                        prop_nr=property_lookup['EC number'],
                        is_qualifier=True))

                # KO number
                elif key2 == 'ko':
                    qual_list.append(wdi_core.WDExternalID(
                        value=value2,
                        prop_nr=property_lookup['KO number'],
                        is_qualifier=True))

                # organism
                elif key2 == 'organism':
                    qual_list.append(wdi_core.WDString(
                        value=value2,
                        prop_nr=property_lookup['Organism'],
                        is_qualifier=True))

                # UniProt protein ID
                elif key2 == 'Hit_accession':
                    qual_list.append(wdi_core.WDExternalID(
                        value=value2,
                        prop_nr=property_lookup['UniProt protein ID'],
                        is_qualifier=True))

                # Expect Value
                elif key2 == 'Hsp_evalue':
                    qual_list.append(wdi_core.WDString(
                        value=str(value2),
                        prop_nr=property_lookup['Expect Value'],
                        is_qualifier=True))

                # Bit Score
                elif key2 == 'Hsp_bit-score':
                    qual_list.append(wdi_core.WDString(
                        value=str(value2),
                        prop_nr=property_lookup['Bit Score'],
                        is_qualifier=True))

                else:
                    pass

            # Alignment hit number
            statements.append(wdi_core.WDString(
                value=str(WDI_dict['UniProt dict']['Hit_nmr']),
                prop_nr=property_lookup['Alignment'],
                references=[Blast_ref],
                qualifiers=qual_list
            ))

        # status
        elif key == "status":
            logging.info("Parsing status")
            statements.append(wdi_core.WDItemID(
                value=item_lookup[value],
                prop_nr=property_lookup["Status"],
                references=[iGEM_ref],
                qualifiers=[datetime_qual]))

        # part ID
        elif key == 'part_id':
            aliases.append(str(value))
            statements.append(wdi_core.WDString(
                value=str(value),
                prop_nr=property_lookup["part ID"],
                references=[iGEM_ref],
                qualifiers=[datetime_qual]))

        # nickname
        elif key == 'nickname':
            aliases.append(value)

        else:
            logging.info("skipped " + key)

    # finding parts page
    parts_page_identifier = get_item_by_name(WDI_dict['part name'], endpoint_url)

    # if parts page present, update and not write.
    if parts_page_identifier is not None:
        logging.info("Part " + label + " " + parts_page_identifier.strip(
            "Q") + " already exists, writing update")

        parts_page = wdi_core.WDItemEngine(
            wd_item_id=parts_page_identifier,
            new_item=False,
            data=statements,
            mediawiki_api_url=mediawiki_api_url,
            sparql_endpoint_url=endpoint_url)

        if len(description) > 1:
            parts_page.set_description(description=description, lang="en")
        else:
            pass

        parts_page.set_aliases(aliases=aliases, lang="en", append=False)

        print(parts_page.wd_json_representation)
        parts_page.write(login_instance)
        logging.info("part " + label + " page is updated")

    else:
        parts_page = wdi_core.WDItemEngine(
            new_item=True,
            data=statements,
            mediawiki_api_url=mediawiki_api_url,
            sparql_endpoint_url=endpoint_url)

        parts_page.set_label(label=label, lang="en")

        # checks if description is present, then adds discription
        if len(description) > 1:
            parts_page.set_description(description=description, lang="en")
        else:
            pass

        parts_page.set_aliases(aliases=aliases, lang="en")

        print(parts_page.wd_json_representation)

        parts_page.write(login_instance)
        logging.info("part " + label + " page is created")
