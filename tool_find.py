#!/usr/bin/env python3

import csv
import string
import sys

"""Heuristically search for interesting names and terms.

This tool is designed to find names of interest using different heuristic
methods. Command line arguments are used to specify the input and output
files.

Usage:
    tool_find.py input.csv output.csv

"""


def parens(s):
    """Return a substring from between a pair of parentheses.

    Args:
        s: a string
    Returns:
        The return value. If no parentheses pair is found, returns an empty
        string.
    """

    if '(' in s and ')' in s:
        result = s[s.index('(') + 1: s.index(')')]
        if len(result) > 1 and len(result) < 10:
            return result
    else:
        return ''


def caps(s):
    """Return substrings that do not contain lowercase characters.

    Args:
        s: a string
    Returns:
        A string composed of substrings that meet the criteria
    """
    tokens = s.split(' ')
    result = []

    def is_caps(t):
        for char in t:
            if char in string.ascii_lowercase:
                return False
        return True

    for token in tokens:
        token = token.strip('()-[]')
        if is_caps(token):
            result.append(token)

    if len(result) > 1 and any(c.isupper() for c in result):
        return ' '.join(result)
    else:
        return ''


def initials(s):
    """


    """
    tokens = [t.strip('()-[]') for t in s.split(' ')]
    # Not yet implemented
    # freebies = ['a', 'in', 'the', 'an', 'of', 'with']

    result = []

    current_run = []

    for t in tokens:
        if len(t) and t[0] in string.ascii_uppercase:
            current_run.append(t)
        else:
            if len(current_run) > 1:
                result.append(' '.join(current_run))
            current_run = []

    if len(current_run) > 1:
        result.append(' '.join(current_run))

    if len(result) > 0:
        return result

    else:
        return ''


def tokenize(s):
    result = []
    tmp = s.split(' ()[]')
    for each in tmp:
        result.append(each.strip(' []()'))

    return result


def export_csv(filename):

    with open(filename, 'w', newline='') as f_out:
        writer = csv.writer(f_out, delimiter=',')
        writer.writerow(['Title', 'Parentheses', 'Abbreviation', 'Initials'])

        for key, value in sorted(results.items()):
            try:
                writer.writerow([key, value[0], value[1], str(value[2])])
            except IndexError:
                writer.writerow([key, value])
    pass


if __name__ == '__main__':
    titles = []
    results = {}
    with open(sys.argv[1]) as f_in:
        reader = csv.reader(f_in)
        for row in reader:
            titles.append(row[0])

    for title in titles:

        result = []

        result.append(parens(title))

        result.append(caps(title))

        result.append(', '.join(initials(title)))

        if len(result) > 0:
            results[title] = result
        else:
            results[title] = ['', '', '']

    export_csv(sys.argv[2])
