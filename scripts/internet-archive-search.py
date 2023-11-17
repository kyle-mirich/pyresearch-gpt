from internetarchive import search_items, get_item, download
import os
import requests
import json

def sanitize_filename(filename):
    return "".join(c for c in filename if c.isalnum() or c in (' ', '.', '_')).rstrip()

def search_archive(query, max_results=50):
    results = []
    for item in search_items(query).iter_as_items():
        if len(results) >= max_results:
            break
        metadata = item.item_metadata["metadata"]
        result = {
            "identifier": item.identifier,
            "title": metadata.get("title"),
            "creator": metadata.get("creator"),
            "date": metadata.get("date"),
            "description": metadata.get("description"),
            "license": metadata.get("licenseurl")
        }
        results.append(result)
    return results

def download_files(identifier, file_types, download_path):
    item = get_item(identifier)
    download(item.identifier, destdir=download_path, formats=file_types)

def create_source_file(filename):
    source_path = os.path.join("sources", f"{filename}.json")
    with open(source_path, 'w') as f:
        json.dump({}, f)
    return source_path

def append_to_sources(identifier, title, download_path, source_file_path):
    item = get_item(identifier)
    metadata = item.metadata

    metadata_dict = {
        "identifier": identifier,
        "title": metadata.get('title', 'None'),
        "collection": metadata.get('collection', []) or [],
        "creator": metadata.get('creator', 'None') or 'None',
        "date": metadata.get('date', 'None') or 'None',
        "language": metadata.get('language', 'None') or 'None',
        "mediatype": metadata.get('mediatype', 'None') or 'None',
        "scanner": metadata.get('scanner', 'None') or 'None',
        "subject": metadata.get('subject', []) or [],
        "uploader": metadata.get('uploader', 'None') or 'None',
        "publicdate": metadata.get('publicdate', 'None') or 'None',
        "addeddate": metadata.get('addeddate', 'None') or 'None',
    }

    try:
        with open(source_file_path, 'r') as f:
            sources_data = json.load(f)
    except FileNotFoundError:
        sources_data = {}

    key = f"{identifier},{title}"
    sources_data[key] = metadata_dict

    with open(source_file_path, 'w') as f:
        json.dump(sources_data, f, indent=4)
