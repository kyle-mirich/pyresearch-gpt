import requests

def search_books(query):
    url = f"https://openlibrary.org/search.json?q={query}"
    response = requests.get(url)
    return response.json().get('docs', [])


def get_book_details(olid):
    url = f"https://openlibrary.org/api/books?bibkeys=OLID:{olid}&format=json"
    response = requests.get(url)
    return response.json().get(f'OLID:{olid}', {})


def search_authors(query):
    url = f"https://openlibrary.org/search/authors.json?q={query}"
    response = requests.get(url)
    return response.json().get('docs', [])


def get_author_details(author_id):
    url = f"https://openlibrary.org/authors/{author_id}.json"
    response = requests.get(url)
    return response.json()


def search_subjects(query):
    url = f"https://openlibrary.org/subjects/{query}.json"
    response = requests.get(url)
    return response.json()


def get_work_details(work_id):
    url = f"https://openlibrary.org/works/{work_id}.json"
    response = requests.get(url)
    return response.json()


def get_editions_by_work(work_id):
    url = f"https://openlibrary.org/works/{work_id}/editions.json"
    response = requests.get(url)
    return response.json().get('entries', [])


def get_covers_by_olid(olid, size='M'):
    covers = {
        'small': f'http://covers.openlibrary.org/b/olid/{olid}-S.jpg',
        'medium': f'http://covers.openlibrary.org/b/olid/{olid}-M.jpg',
        'large': f'http://covers.openlibrary.org/b/olid/{olid}-L.jpg'
    }
    return covers.get(size, '')
