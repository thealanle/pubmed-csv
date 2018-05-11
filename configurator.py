#!/usr/bin/env python3

import configparser
import os

DEFAULT_CFG = 'pubmed-csv.cfg'


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

    # Write properties to file.
    with open(filename, 'w') as configfile:
        config.write(configfile)

    config.read(filename)

    return config


def setup(filename):
    """
    Ingests or creates a .cfg file. The result is a ConfigParser() object.
    """

    result = None

    # Read the .cfg file. If it does not exist, create one and query user for properties.
    if os.path.isfile(filename):
        result = configparser.ConfigParser()
        # .read(filename) is an in-place operation and must be called before returning result
        result.read(filename)
    else:
        result = new_cfg(filename)
    return result


if __name__ == '__main__':
    setup(DEFAULT_CFG)
