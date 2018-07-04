#!/usr/bin/env python3

from bs4 import BeautifulSoup
import csv


class Library:
    """
    TODO: Implement a container to enable access and modification of Documents.
    """

    def __init__(self, filename):
        data_in = open(filename)
        contents = data_in.read()
        soup = BeautifulSoup(contents, 'xml')
        data_in.close()
        self.article_list = soup.find_all(['article'])


class Document:
    """
    TODO: Implement a class representation of each journal article.
    """
    pass


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


def print_all(article):
    text_list = article.find_all('p')
    for each in text_list:
        print(each.get_text())


def print_metadata(article):
    front = article.find('front')
    journal_title = front.find('journal-title').string
    print(journal_title)
    article_title = front.find('article-title').string
    print(article_title)
    first_author = front.find('contrib').find('surname').string
    author_date = first_author + ' ' + article.find('year').string
    print(author_date)


def print_abstract(article):
    abstract = article.find('abstract')
    print(abstract.get_text(), end='')


def print_body(article):
    body = article.find('body')
    p_list = [p.get_text() for p in body.find_all('p') if p.parent.name ==
              'sec']
    text_list = [p.replace('\n', '') for p in p_list]
    print(' '.join(text_list))


if __name__ == '__main__':
    library = Library("xml_input/efetch-pmc.xml")
    for article in library.article_list:
        print_metadata(article)
        print_body(article)
