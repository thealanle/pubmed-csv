#!/usr/bin/env python3

from bs4 import BeautifulSoup as bs

infile = open("xml_input/efetch-pmc.xml")
contents = infile.read()
soup = bs(contents, 'lxml')


articles = soup.find_all('sec')
for article in articles:
    print(article.get_text())

infile.close()
