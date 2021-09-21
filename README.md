# BioParts Wikibase

The Python pipeline is used to create the first framework of a unified semantic database for the
biobricks in the current [iGEM Registry](http://parts.igem.org/Main_Page). This Wikibase instance is stored in 
a Resource Description Framework (RDF) format that makes the data easily retrievable using a triple format query. 
The pipeline annotates the biological parts with protein functions, restriction sites from restriction
enzymes, and interconnects information of parts of the biobricks that are composed of other biological parts. 
Wikibase is a semantic database. The [BioParts Wikibase](https://bioparts.wiki.opencura.com/wiki/Main_Page) 
is an instance that helps navigate the iGEM registry and it is constructed in the iGEM 2021 competition and to be updated upon
by the iGEM community. This git repository contains a current version and previous versions used to
upload to the BioParts Wikibase. Additionally, this git repository also contains query examples to
retrieve information intelligently from the BioParts Wikibase.


## Installation of modules

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install the following packages.


[WikidataIntegrator](https://github.com/SuLab/WikidataIntegrator) 0.7.4 This package is used to read and write to the 
Wikibase. It is specifically created to upload biological
information to Wikibase instances. The package is used to create statements in the correct format
and write these to the designated item page. Additionally, it also acts as a wrapper for the queries
performed to the BioParts Wikibase to retrieve information.
```bash
pip3 install wikidataintegrator
```
[Biopython](https://biopython.org/wiki/Download) 1.77 This package is used to find restriction sites in the nucleotide 
sequences of the biological parts and
makes sure that the annotation of this is correct with the standards as indicated on
 [the iGEM registry their assembly information pages](http://parts.igem.org/Help:Assembly_Compatibility).
```bash
pip3 install biopython 
```
[XlsxWriter](https://xlsxwriter.readthedocs.io/index.html) 1.3.7 This package is used to write the output from the
 queries to the Wikibase instance to the excel file.
```bash
pip3 install xlsxwriter
```
[SPARQLWrapper](https://pypi.org/project/SPARQLWrapper/) 1.8.5 This package is used as a python SPARQL wrapper when
 querying external databases, such as the [UniProt](https://sparql.uniprot.org/sparql) database to retrieve additional information on the biological part. 
```bash
pip install SPARQLWrapper
```
[Beautiful Soup](https://www.crummy.com/software/BeautifulSoup/bs4/doc/) 4.9.1 This package is used to parse XML and HTML
 formatted data. It is used both for parsing of the point in time [datadump file](http://parts.igem.org/Registry_API) 
 and for the webpages used for storage of nucleotide sequences of the biological parts. 
```bash
pip3 install beautifulsoup4
```

## Installation of DIAMOND
[The DIAMOND protein aligner](http://www.diamondsearch.org/index.php) This external programme is used as Basic Local 
 Alignment Search Tool (BLAST). To install, follow the start guide on their website to set up DIAMOND.

## Usage
There are two methods of usage, the creation of a Wikibase instance and the query of such instance.

## To construct the Wikibase

The python scripts can be found under ../Script. For more information regarding the pipeline and functions used by the
MAIN.py file please consult to the manual under ./Manual. 

Before usage the following paths have to be provided in the MAIN.py file. 

```python3
input_path      = # path to XML file 
blasted_file    = # path to blasted file location (standard is '../Parts/Db_output.xml')
fasta_loc       = # path to fasta file location (standard is '../Parts/fastafile.fna')
T_directory     = # path to temporary pickle file directory (standard is '../Parts/Temp_pickle/') 
F_directory     = # path to final pickle file directory (standard is '../Parts/Temp_pickle/')
database        = # path to database .dmnd file
```

Before usage the following URLs have to be provided in the MAIN.py file.

For the Wikibase
```python3
endpoint_url        = # Wikibase SPARQL endpoint
mediawiki_api_url   = # Wikibase API
```

For UniProt
```python3
Sparql_endpoint     = # UniProt SPARQL endpoint (standard is 'http://purl.uniprot.org/uniprot/')
```

For iGEM
```python3
iGEM_sequence_url = # iGEM HTML sequence (standard is 'http://parts.iGEM.org/cgi/partsdb/composite_edit/putseq.cgi?part=')
```

To use the script and when new final pickle files have to be created, run the following code. The "new" parameter indicates that the parsing
and aligning of biological parts has not been done prior and therefore no dictionary objects are present (or they are to
be replaced). The username and password parameters are to authenticate the bot and make sure that changes can be made to
the Wikibase instance. Further documentation is present in the manual. 
```bash
python3 MAIN.py new username password
```

On the contrary if final pickle files can be used again, run the following code. The "old" parameter indicates that
 the previous pickled dictionary objects can be used. This has decreased runtime and can be used when changed to WDI_writer.py
is updated. Further documentation is present in the manual.
```bash
python3 MAIN.py old username password
```

Running Add_assembly.py is always done after pickled dictionary objects have been created. The script adds links between 
biobrick item pages via the "contains" statement. Recommended to run after the previous script, but running in 
tandem is also possible. Again the username and password have to be provided as arguments. 
```bash
python3 Add_assembly.py username password
```

## To query the Wikibase 

The Wikibase uses as triple formatted query. At present only the query method used by [WikidataIntegrator](https://github.com/SuLab/WikidataIntegrator)
can be used to read from the Wikibase and no other SPARQL wrappers. Additionally, the prefixes need to be set manually. 
For the [BioParts Wikibase](https://bioparts.wiki.opencura.com/wiki/Main_Page) the following prefixes have to be set. 

For this specific Wikibase
```python
PREFIX wd: <http://bioparts.wiki.opencura.com/entity/>
PREFIX wdt: <http://bioparts.wiki.opencura.com/prop/direct/>
PREFIX p: <http://bioparts.wiki.opencura.com/prop/>
PREFIX ps: <http://bioparts.wiki.opencura.com/prop/statement/>
PREFIX pq: <http://bioparts.wiki.opencura.com/prop/qualifier/>
```
Standard unchanged prefixes 
```python
PREFIX wikibase: <http://wikiba.se/ontology#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX bd: <http://www.bigdata.com/rdf#>
```

A specific example of a query incorporated in a python script can be found under ./Manual/query_example. 
Which performs a SPARQL query to the [BioParts Wikibase](https://bioparts.wiki.opencura.com/wiki/Main_Page)
using the wrapper provided by [WikidataIntegrator](https://github.com/SuLab/WikidataIntegrator). 
It searches for an item given an Enzyme Commision number (EC number) and can be used as guidelines to construct a 
personalised query script. The current query example uses the xlsxwriter package to create an .xlsx (excel) file as 
output. This can easily be changed if another output format is required.

Before using the following has to be provided.

```python
endpoint_url = # URL to Wikibase SPARQL endpoint
Out_name = # Output file name
```

## Contribution

Riemer van der Vliet,
Jasper Koehorst,
Robert Smith

## Licence

Copyright 2020, Laboratory of Systems and Synthetic Biology 
