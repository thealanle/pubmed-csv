#!/usr/bin/env python3

import csv
import sys


def extract_id_pairs(path_in, path_out):
    with open(path_in, 'r') as f_in, open(path_out, 'a') as f_out:

        reader = csv.reader(f_in)
        writer = csv.writer(f_out)

        for row in reader:
            if row[8] and row[9]:
                writer.writerow([row[8], row[9]])


def get_pmcids(map_path, pmid_path, path_out):
    pairs = set()
    pmids = set()
    with open(map_path, 'r') as map_in, open(pmid_path, 'r') as pmid_in:

        map_reader = csv.reader(map_in)
        pmid_reader = csv.reader(pmid_in)

        for row in map_reader:
            if len(row) == 2:
                pairs.add((row[0], row[1]))

        for row in pmid_reader:
            pmids.add(row[0])

    inter = {x for x, y in pairs if y in pmids}

    with open(path_out, 'a') as f_out:
        for pmcid in inter:
            f_out.write(pmcid + '\n')


if __name__ == '__main__':
    args = sys.argv
    get_pmcids(args[1], args[2], args[3])
