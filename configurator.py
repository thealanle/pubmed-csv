#!/usr/bin/env python3

import configparser
import os


def new_cfg(filename):
    """
    Creates a new .cfg file and returns the new config object
    """
    config = configparser.ConfigParser()

    # Ingest user properties
    config['User'] = {}
    config['User']['Email'] = input(
        "E-mail address for use with NCBI:\n>").strip()
    config['User']['API Key'] = input("NCBI API key:\n>").strip()

    # Ingest file and directory locations
    config['Paths'] = {}
    config['Paths']['PMID Input'] = input(
        "Path to file containing PMIDs:\n>").strip()
    config['Paths']['Output Filename'] = input("Output filename:\n>").strip()

    with open(filename, 'w') as configfile:
        config.write(configfile)

    config.read(filename)

    return config


def setup(filename):
    result = None
    if os.path.isfile(filename):
        result = configparser.ConfigParser()
        result.read(filename)
    else:
        result = new_cfg(filename)
    return result


if __name__ == '__main__':
    setup('pubmed-csv.cfg')
