#!/usr/bin/env python3

from bs4 import BeautifulSoup as bs


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


def print_article_text(article):
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
    print('\n')
    # for child in front:
    #     print(front.get_text())


def print_full_text(article):
    body = article.find('body')
    p_list = body.find_all('p')
    # secs = [sec for sec in body.find_all(
    #     'sec') if 'supplementary' not in sec.attrs['sec-type']]

    print(secs.get_text())
    for p in p_list:
        print(p.get_text())


def print_abstract(article):
    abstract = article.find('abstract')
    print(abstract.get_text())


if __name__ == '__main__':
    infile = open("xml_input/efetch-pmc.xml")
    contents = infile.read()
    soup = bs(contents, 'xml')

    article_list = soup.find_all(['article'])
    # print_article_text(article_list)

    for article in article_list:
        print_metadata(article)
        print_abstract(article)
        print_full_text(article)

    infile.close()
