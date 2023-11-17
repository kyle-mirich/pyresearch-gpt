import wikipediaapi

def wikipedia_search(query, language='en'):
    wiki = wikipediaapi.Wikipedia(language)
    return wiki.search(query)

def get_wikipedia_summary(title, language='en'):
    wiki = wikipediaapi.Wikipedia(language)
    page = wiki.page(title)
    return page.summary

def get_wikipedia_article(title, language='en'):
    wiki = wikipediaapi.Wikipedia(language)
    page = wiki.page(title)
    return page.text

def get_wikipedia_links(title, language='en'):
    wiki = wikipediaapi.Wikipedia(language)
    page = wiki.page(title)
    return [link.title for link in page.links.values()]

def get_wikipedia_categories(title, language='en'):
    wiki = wikipediaapi.Wikipedia(language)
    page = wiki.page(title)
    return [category.title for category in page.categories.values()]

def get_wikipedia_translations(title, language='en'):
    wiki = wikipediaapi.Wikipedia(language)
    page = wiki.page(title)
    return page.langlinks

def get_wikipedia_sections(title, language='en'):
    wiki = wikipediaapi.Wikipedia(language)
    page = wiki.page(title)
    return [section.title for section in page.sections]

def get_wikipedia_images(title, language='en'):
    wiki = wikipediaapi.Wikipedia(language)
    page = wiki.page(title)
    return [image.url for image in page.images]
