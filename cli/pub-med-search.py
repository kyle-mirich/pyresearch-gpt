import argparse
from Bio import Entrez
import os
import json
from colorama import Fore, Style, init
import csv

def create_csv_file(article_details_list, query):
    csv_filename = f"{query.replace(' ', '_')}.csv"
    with open(csv_filename, 'w', newline='', encoding='utf-8') as csvfile:
        fieldnames = ['title', 'authors', 'journal', 'publication_date', 'abstract']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for article in article_details_list:
            writer.writerow({
                'title': article['Title'],
                'authors': article['Authors'],
                'journal': article['Journal'],
                'publication_date': article['Publication Date'],
                'abstract': article['Abstract']
            })
    print(f"{Fore.GREEN}CSV file created: {csv_filename}")
# Initialize colorama
init(autoreset=True)

def print_box(message, color=Fore.RESET):
    print(color + "-" * 120)
    print(message.center(120))
    print("-" * 120)

def create_directory(query):
    directory_name = query.replace(" ", "_")
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)
    return directory_name

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
    author_names = []
    for author in authors:
        forename = author.get('ForeName', '')
        lastname = author.get('LastName', '')
        author_names.append(f"{forename} {lastname}".strip())
    return author_names


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


# New function for searching and previewing PubMed articles
def search_and_preview(query, max_results, email):
    Entrez.email = email
    handle = Entrez.esearch(db='pubmed', term=query, retmax=max_results)
    record = Entrez.read(handle)
    pmid_list = record['IdList']

    print(f"\n{Fore.GREEN}Found {len(pmid_list)} results for '{query}':\n")
    for i, pmid in enumerate(pmid_list):
        print(f"{Fore.BLUE}{i + 1}: {Fore.GREEN}PMID: {pmid}")
    return pmid_list


def display_article_details(pmid, email, write_to_file=False):
    article = get_pubmed_article(pmid, email)
    abstract = get_abstract_from_pubmed_article(pmid, email)
    authors = get_authors_from_pubmed_article(pmid, email)
    journal = get_journal_from_pubmed_article(pmid, email)
    pub_date = get_publication_date_from_pubmed_article(pmid, email)

    article_details = {
        "Title": article['PubmedArticle'][0]['MedlineCitation']['Article']['ArticleTitle'],
        "Authors": ', '.join(authors),
        "Journal": journal['title'],
        "Publication Date": pub_date,
        "Abstract": abstract
    }

    for key, value in article_details.items():
        print(f"{Fore.CYAN}{key}: {value}")

    if write_to_file:
        write_article_details_to_file(article_details, pmid)

def create_combined_file(directory, query):
    combined_filename = f"{directory}/{query.replace(' ', '_')}-combined.txt"
    with open(combined_filename, 'w', encoding='utf-8') as combined_file:
        for file in os.listdir(directory):
            if file.endswith('_details.txt'):
                with open(f"{directory}/{file}", 'r', encoding='utf-8') as individual_file:
                    combined_file.write(individual_file.read())
                    combined_file.write("\n\n")
    print(f"{Fore.GREEN}Combined file created: {combined_filename}")

def write_article_details_to_file(article_details, pmid, directory):
    filename = f"{directory}/PMID_{pmid}_details.txt"
    if os.path.exists(filename):
        base, extension = os.path.splitext(filename)
        counter = 1
        while os.path.exists(f"{base}_{counter}{extension}"):
            counter += 1
        filename = f"{base}_{counter}{extension}"

    with open(filename, 'w', encoding='utf-8') as file:
        for key, value in article_details.items():
            file.write(f"{key}: {value}\n")
    print(f"{Fore.GREEN}Article details written to {filename}")


def display_article_details(pmid, email, directory, write_to_file=False):
    article = get_pubmed_article(pmid, email)
    abstract = get_abstract_from_pubmed_article(pmid, email)
    authors = get_authors_from_pubmed_article(pmid, email)
    journal = get_journal_from_pubmed_article(pmid, email)
    pub_date = get_publication_date_from_pubmed_article(pmid, email)

    article_details = {
        "Title": article['PubmedArticle'][0]['MedlineCitation']['Article']['ArticleTitle'],
        "Authors": ', '.join(authors),
        "Journal": journal['title'],
        "Publication Date": pub_date,
        "Abstract": abstract
    }

    for key, value in article_details.items():
        print(f"{Fore.CYAN}{key}: {value}")

    if write_to_file:
        write_article_details_to_file(article_details, pmid, directory)

def main():

    email = 'kylemirich20@gmail.com'
    print_box(f"{Fore.RED}Welcome to PubMed Article Fetcher!\n", Fore.RED)

    while True:
        try:
            query = input(Fore.GREEN + "Enter search query: ")
            max_results = int(input(Fore.GREEN + "Enter the number of results to preview: "))
            id_list = search_and_preview(query, max_results, email)

            query_directory = create_directory(query)
            if input(Fore.GREEN + "Do you want to save all results to files? (yes/no): ").lower() == 'yes':
                for pmid in id_list:
                    display_article_details(pmid, email, query_directory, write_to_file=True)
                create_combined_file(query_directory, query)
                break

            pmid_input = input(Fore.GREEN + "Enter the number of the article for details (or 'exit' to quit): ")
            if pmid_input.lower() == 'exit':
                break

            pmid = id_list[int(pmid_input) - 1]
            display_article_details(pmid, email, query_directory, write_to_file=True)
            create_combined_file(query_directory, query)

        except Exception as e:
            print(f"{Fore.RED}An unexpected error occurred: {e}\n")

    print_box("Thank you for using PubMed Article Fetcher!", Fore.RED)

if __name__ == "__main__":
    main()
