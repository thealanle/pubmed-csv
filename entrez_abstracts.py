#!/usr/bin/env python3

from Bio import Entrez
import csv


# The config file should consist of three lines, of the format:
# <user@example.com>
# <NCBI_api_key>
# <path/to/pmids.txt>
CONFIG_PATH = 'pubmed-csv-config.txt'


def config_setup(path):
    """
    Given a path to the config file, returns a tuple:
    (user@example.com, <NCBI_api_key>)
    If the config file does not exist, user is prompted to create one.
    """

    email = None
    api_key = None
    pmid_path = None

    try:
        with open(path, 'r') as config:
            print("Config file found at {}".format(path))
            email = config.readline().strip()
            api_key = config.readline().strip()
            pmid_path = config.readline().strip()
    except FileNotFoundError:
        query = input(
            "No config file found at '{}'. Create a config file now? y/N\n >".format(path)).lower()
        if query != 'y':
            print("No config file will be created.")
            pass
        else:
            email = input("Input an e-mail address for use with NCBI:\n>")
            api_key = input("Input NCBI API key:\n>")
            pmid_path = input("Input path to file containing PMIDs:\n>")
            with open(path, 'w') as config:
                config.write(email)
                config.write(api_key)
                config.write(pmid_path)
            print("Config file created at {}".format(path))

    return (email, api_key, pmid_path)


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


def write_csv_file(pmid, title, abstract, f_out):
    """
    Writes article properties to a given output .csv file.
    """
    writer = csv.writer(f_out)
    writer.writerow([pmid, title, abstract])


if __name__ == '__main__':

    Entrez.email, Entrez.api_key, pmid_path = config_setup(CONFIG_PATH)
    pmid_list = ids_from_file(pmid_path)

    for sublist in pmid_list:

        with Entrez.efetch(db='pubmed',
                           id=','.join(sublist),
                           retmode='xml') as handle:

            records = Entrez.read(handle)

            # Writes to a .csv file in 'append' mode.
            with open('abstracts_out.csv', 'a') as csv_out:

                csv_writer = csv.writer(csv_out)
                csv_writer.writerow(['PMID', 'Title', 'Abstract'])

                for record in records['PubmedArticle']:
                    pmid_text = record['MedlineCitation']['PMID']
                    title_text = record['MedlineCitation']['Article']['ArticleTitle']
                    abs_text = get_abs_text(record)
                    write_csv_file(pmid_text, title_text,
                                   abs_text, csv_out)

            # # Writes to a .txt file. TO-DO: Implement write_txt_file.
            # with open('abstracts_out.txt', 'a') as txt_out:
            #     for record in records['PubmedArticle']:
            #         title_text = record['MedlineCitation']['Article']['ArticleTitle']
            #         abs_text = get_abs_text(record)
