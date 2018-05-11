# pubmed-csv

Under construction using the [biopython](https://github.com/biopython/biopython) and [pdfminer.six](https://github.com/pdfminer/pdfminer.six) packages.
Data is collected using the [NCBI Entrez E-Utilities API](https://www.ncbi.nlm.nih.gov/books/NBK25501/).

## Objectives:
- ingest articles in PDF form and output them as text.
- retrieve article data given a file containing PMIDs
- output article properties to .csv format


## Usage:

Make sure biopython is installed, as the Bio.Entrez module is required.
```
python3 -m pip install --user biopython
```

The script entrez_abstracts.py can then be called with python3:
```
python3 entrez_abstracts.py
```

On first-time use, the user will be prompted to enter an e-mail address for communication with NCBI, as well as their [API key](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/). This will allow the user to be contacted in case they are making excessive calls to NCBI services.
