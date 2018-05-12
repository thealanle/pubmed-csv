#!/usr/bin/env python3

from Bio import Entrez
import csv
import configurator


CONFIG_FILE = 'pubmed-csv.cfg'


def ids_from_file(filepath):
    """
    Given a file with a PMID on each line, returns a list of lists containing
    PMIDs as strings. Due to memory restrictions with Entrez.read, PMIDs are
    batched in groups of 200.
    """

    result = [[]]  # A list containing an empty list.
    count = 0

    with open(filepath, 'r') as f_in:
        for line in f_in:
            if count < 200:
                # Append to the last sublist in the list
                result[-1].append(line.strip())
                count += 1
            else:
                # Create a new sublist once a group of 200 is reached.
                result.append([line.strip()])
                count = 1

    return result


def get_abs_text(record):
    """
    Given a record from Entrez.read(handle), attempt to return a string
    containing the text of the article's abstract.
    """

    result = []

    try:
        for item in record['MedlineCitation']['Article']['Abstract']['AbstractText']:
            result.append(item)
    except KeyError:
        pass

    if len(result) == 0:
        return '*Abstract Unavailable*'
    else:
        return ' '.join(result)


def write_txt_file(title, abstract, f_out):
    """
    Writes article properties to a given output .txt file.
    """

    print(
        '*' * 40,
        '\n',
        '[Title] {}'.format(title),
        '\n',
        '[Abstract] {}'.format(abstract),
        file=f_out
    )


def write_csv_file(records, output_file):
    """
    Writes rows to output_file in append mode, which is required for datasets
    with more than 200 articles.
    """
    with open(output_file, 'a') as csv_out:

        csv_writer = csv.writer(csv_out)
        csv_writer.writerow(['PMID', 'Title', 'Abstract'])

        for record in records['PubmedArticle']:
            pmid_text = record['MedlineCitation']['PMID']
            title_text = record['MedlineCitation']['Article']['ArticleTitle']
            abs_text = get_abs_text(record)
            csv_writer.writerow([pmid_text, title_text, abs_text])


if __name__ == '__main__':
    config = configurator.setup(CONFIG_FILE)

    Entrez.email = config['User']['Email']
    Entrez.api_key = config['User']['API key']
    output_file = config['Paths']['Output Filename']
    pmid_list = ids_from_file(config['Paths']['PMID Input'])

    for sublist in pmid_list:
        with Entrez.efetch(db='pubmed', id=','.join(sublist),
                           retmode='xml') as handle:
            records = Entrez.read(handle)

            write_csv_file(records, output_file)
