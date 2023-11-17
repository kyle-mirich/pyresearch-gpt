from serpapi import GoogleSearch
import requests
from bs4 import BeautifulSoup
from textwrap import wrap

def search(query, api_key):
    params = {
        "q": query,
        "api_key": api_key
    }
    search = GoogleSearch(params)
    results = search.get_dict()
    return results['organic_results']

def click_link(url):
    response = requests.get(url)
    if response.status_code == 200:
        soup = BeautifulSoup(response.content, 'html.parser')
        return soup.prettify()
    else:
        return None

def read_content(url, chunk_size=1000):
    html_content = click_link(url)
    if html_content:
        soup = BeautifulSoup(html_content, 'html.parser')
        text_content = soup.get_text()
        chunks = wrap(text_content, chunk_size)
        chunked_content = [{"part": f"{i+1}/{len(chunks)}", "content": chunk} for i, chunk in enumerate(chunks)]
        return {
            "url": url,
            "total_parts": len(chunks),
            "chunks": chunked_content
        }
    else:
        return None


from gensim.summarization import summarize

def generate_summary(url, ratio=0.2):
    content = read_content(url)
    if content:
        text = ' '.join([chunk['content'] for chunk in content['chunks']])
        summary = summarize(text, ratio=ratio)
        return summary
    else:
        return None


def highlight_keywords(content, keywords):
    for keyword in keywords:
        content = content.replace(keyword, f'**{keyword}**')
    return content


from googletrans import Translator

def translate_content(content, target_language):
    translator = Translator()
    translation = translator.translate(content, dest=target_language)
    return translation.text


from textblob import TextBlob

def analyze_sentiment(content):
    blob = TextBlob(content)
    return blob.sentiment.polarity, blob.sentiment.subjectivity


import json

def save_bookmark(url, filename='bookmarks.json'):
    try:
        with open(filename, 'r') as file:
            bookmarks = json.load(file)
    except FileNotFoundError:
        bookmarks = []

    bookmarks.append(url)

    with open(filename, 'w') as file:
        json.dump(bookmarks, file)

class HistoryTracker:
    def __init__(self):
        self.history = []

    def add_to_history(self, action, details):
        self.history.append({"action": action, "details": details})

    def get_history(self):
        return self.history


def extract_images(url):
    html_content = click_link(url)
    soup = BeautifulSoup(html_content, 'html.parser')
    images = [img['src'] for img in soup.find_all('img') if 'src' in img.attrs]
    return images

def extract_videos(url):
    html_content = click_link(url)
    soup = BeautifulSoup(html_content, 'html.parser')
    videos = [video['src'] for video in soup.find_all('video') if 'src' in video.attrs]
    return videos

def extract_news(query, api_key):
    params = {
        "q": query,
        "api_key": api_key,
        "engine": "google_news"
    }
    search = GoogleSearch(params)
    results = search.get_dict()
    return results['news_results']


def extract_citations(url):
    html_content = click_link(url)
    soup = BeautifulSoup(html_content, 'html.parser')
    citations = [citation.get_text() for citation in soup.find_all('cite')]
    return citations


import matplotlib.pyplot as plt

def plot_trends(data, title='Trends over Time', xlabel='Time', ylabel='Value'):
    plt.plot(data)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()


class UserProfile:
    def __init__(self, preferences=None):
        self.preferences = preferences or {}

    def set_preference(self, key, value):
        self.preferences[key] = value

    def get_preference(self, key):
        return self.preferences.get(key, None)


def batch_download(urls, destination_folder):
    for url in urls:
        response = requests.get(url)
        filename = url.split('/')[-1]
        with open(f'{destination_folder}/{filename}', 'wb') as file:
            file.write(response.content)



def generic_scraper(url):
    page = requests.get(url)
    soup = BeautifulSoup(page.content, "html.parser")
    # Modify the following lines to suit the specific website's HTML structure
    title = soup.find("title").text
    content = soup.find("div", {"class": "content-class"}).text # Change the class name as per the website
    return {
        "title": title,
        "content": content
    }