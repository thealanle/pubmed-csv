#!/usr/bin/env python3

from bs4 import BeautifulSoup
import csv
import collections


class Library:
    """
    A container for Documents. Handles writing of aggregate data to a .CSV file.
    """

    def __init__(self, filename):
        data_in = open(filename)
        contents = data_in.read()
        soup = BeautifulSoup(contents, 'xml')
        data_in.close()

        self.docs = collections.OrderedDict()

        # Creates a Document object for each article and adds it to the
        # Library's docs dict
        for article in soup.find_all(['article']):
            try:
                doc = Document(article)
                self.add_doc(doc)
                # print(doc.meta['body'])
            except AttributeError:
                pass

    # TODO: Update to use PMID or PMCID instead of title
    def add_doc(self, doc):
        self.docs[doc.meta['article-title']] = doc

    def export_csv(self, filename):  # Not yet implemented
        headers = [
            'journal-title',
            'article-title',
            'first-author',
            'year',
            'body'
        ]

        with open(filename, encoding='utf-8', mode='w', newline='\n') as f_out:
            writer = csv.DictWriter(f_out, fieldnames=headers)
            writer.writeheader()
            for doc in self.docs.values():
                writer.writerow(doc.meta)


class Document:
    """
    An article with its properties parsed out for access and modification.
    """

    def __init__(self, article):
        front = article.find('front')
        body = article.find('body')
        self.meta = {
            'journal-title': front.find('journal-title').string,
            'article-title': front.find('article-title').string,
            'first-author': front.find('contrib').find('surname').string,
            'year': front.find('year').string,
            'body': self.get_body(body),
        }

    def __str__(self):
        """
        Return select metadata elements as a printable string.
        """
        summary = [doc.meta['article-title'],
                   doc.meta['first-author'] + ' ' + doc.meta['year'],
                   doc.meta['body']
                   ]
        try:
            return '\n'.join(summary)
        except TypeError:
            return ''

    def get_body(self, body, inline=True):
        """
        Return all primary body text as a single string without linebreaks.
        """

        # Make a list of all strings under <p> tags
        temp = [p.get_text()
                for p in body.find_all('p') if p.parent.name == 'sec']

        # Use tokens to avoid whitespace errors
        text_list = []
        for p in temp:
            tokens = p.split()
            for each in tokens:
                text_list.append(each)

        return ' '.join(text_list) if inline else '\n'.join(text_list)


def print_pmids(article):
    # Get PMID from article
    article_ids = [id_tag for id_tag in article.find_all(
        'article-id') if id_tag.attrs['pub-id-type'] == 'pmid']
    for tag in article_ids:
        # id_list.append(' '.join([tag.attrs['pub-id-type'], tag.string]))
        print(tag.attrs['pub-id-type'], tag.string)
        # Output:
        # pmid 25615823
        # pmid 19686402


if __name__ == '__main__':
    library = Library("xml_input/pmcids-mini.xml")
    for doc in library.docs.values():
        print(doc, '\n')
    # library.export_csv('test.csv')