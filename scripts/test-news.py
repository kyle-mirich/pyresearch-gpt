from gnews import GNews
import json
from newspaper import Article

def gnews_toolkit(action, max_results=10, keyword=None, topic=None, location=None, language='en', country='United States', url=None):
    google_news = GNews()
    google_news.max_results = max_results
    
    if action == 'fetch_top_news':
        return json.dumps(google_news.get_top_news())
    
    if action == 'fetch_news_by_keyword':
        return json.dumps(google_news.get_news(keyword))
    
    if action == 'fetch_news_by_topic':
        return json.dumps(google_news.get_news_by_topic(topic))
    
    if action == 'fetch_news_by_location':
        location_encoded = location.replace(" ", "%20")
        return json.dumps(google_news.get_news_by_location(location_encoded))
    
    if action == 'customize_news':
        google_news.language = language
        google_news.country = country
    
    if action == 'fetch_full_article':
        article = Article(url)
        article.download()
        article.parse()
        full_article_data = {
            'title': article.title,
            'url': url,
            'authors': article.authors,
            'published_date': article.publish_date,
            'text': article.text,
            'images': list(article.images), # Convert images to a list
            # Add other attributes as needed
        }
        return full_article_data

    raise ValueError("Invalid action specified")
