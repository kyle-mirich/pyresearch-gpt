import argparse
from Bio import Entrez
import os
import json
from colorama import Fore, Style, init

# Initialize colorama
init(autoreset=True)

def print_box(message, color=Fore.RESET):
    print(color + "-" * 120)
    print(message.center(120))
    print("-" * 120)


def search_pubmed(query, max_results=10, email='your_email@example.com'):
    Entrez.email = email
    handle = Entrez.esearch(db='pubmed', term=query, retmax=max_results)
    record = Entrez.read(handle)
    return record['IdList']

def get_pubmed_article(pmid, email='your_email@example.com'):
    Entrez.email = email
    handle = Entrez.efetch(db='pubmed', id=pmid, retmode='xml')
    record = Entrez.read(handle)
    return record


def get_abstract_from_pubmed_article(pmid, email='your_email@example.com'):
    article = get_pubmed_article(pmid, email)
    abstract = article['PubmedArticle'][0]['MedlineCitation']['Article']['Abstract']['AbstractText']
    return ' '.join(abstract)


def get_authors_from_pubmed_article(pmid, email='your_email@example.com'):
    article = get_pubmed_article(pmid, email)
    authors = article['PubmedArticle'][0]['MedlineCitation']['Article']['AuthorList']
    return [f"{author['ForeName']} {author['LastName']}" for author in authors]


def get_journal_from_pubmed_article(pmid, email='your_email@example.com'):
    article = get_pubmed_article(pmid, email)
    journal_info = article['PubmedArticle'][0]['MedlineCitation']['Article']['Journal']
    return {
        'title': journal_info['Title'],
        'issn': journal_info['ISSN'],
        'journal_issue': journal_info['JournalIssue']
    }


def get_publication_date_from_pubmed_article(pmid, email='your_email@example.com'):
    article = get_pubmed_article(pmid, email)
    pub_date = article['PubmedArticle'][0]['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']
    return f"{pub_date['Year']} {pub_date.get('Month', '')} {pub_date.get('Day', '')}".strip()


