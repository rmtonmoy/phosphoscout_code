"""
MCP server with literature search and web search tools.

Available Tools:
- web_search: Search the public web for up-to-date information and sources (Tavily)
- query_pubmed: Search PubMed for biomedical literature
- query_arxiv: Search arXiv for papers
- query_scholar: Search Google Scholar for papers
- extract_url_content: Extract text content from a webpage
- extract_pdf_content: Extract text content from a PDF URL

The literature tool implementations are ported from Biomni
(https://github.com/snap-stanford/Biomni, biomni/tool/literature.py) so the
project no longer depends on the biomni package (and its ~15 GB data lake).
"""

import os
import re
import time
from io import BytesIO
from pathlib import Path
from typing import Union
from urllib.parse import urljoin

import PyPDF2
import dotenv
import requests
from bs4 import BeautifulSoup
from mcp.server.fastmcp import FastMCP
from langchain_tavily import TavilySearch

PROJECT_ROOT = Path(__file__).resolve().parents[2]

dotenv.load_dotenv(str(PROJECT_ROOT / '.env'))
mcp = FastMCP("LiteratureTools")


@mcp.tool()
def web_search(query: Union[str, dict]) -> str:
    """Search the public web for up-to-date information and sources

    Args:
        query: A concise query string or a dictionary with a key-value pair where the key is the query

    Returns:
        Summarized results with citations
    """
    if isinstance(query, dict):
        if not query:
            raise ValueError("Dictionary input cannot be empty")
        search_query = list(query.keys())[0]
    else:
        search_query = query

    tavily_search = TavilySearch(
        max_results=5,
        topic="general",
    )

    results = tavily_search.invoke(search_query)

    if isinstance(results, str):
        return results
    elif isinstance(results, dict):
        formatted_results = f"Search results for: {search_query}\n\n"
        if "results" in results:
            for i, result in enumerate(results["results"], 1):
                formatted_results += f"{i}. "
                if "title" in result:
                    formatted_results += f"{result['title']}\n"
                if "content" in result:
                    formatted_results += f"   {result['content']}\n"
                if "url" in result:
                    formatted_results += f"   Source: {result['url']}\n"
                formatted_results += "\n"
        else:
            formatted_results += str(results)
        return formatted_results
    else:
        return str(results)


@mcp.tool()
def query_pubmed(query: str, max_papers: int = 10, max_retries: int = 3) -> str:
    """Query PubMed for papers based on the provided search query.

    Args:
        query: The search query string.
        max_papers: The maximum number of papers to retrieve (default: 10).
        max_retries: Maximum number of retry attempts with simplified queries (default: 3).

    Returns:
        The formatted search results or an error message.
    """
    from pymed import PubMed

    try:
        pubmed = PubMed(tool="PhosphoScout", email="your-email@example.com")

        papers = list(pubmed.query(query, max_results=max_papers))

        retries = 0
        while not papers and retries < max_retries:
            retries += 1
            simplified_query = " ".join(query.split()[:-retries]) if len(query.split()) > retries else query
            time.sleep(1)
            papers = list(pubmed.query(simplified_query, max_results=max_papers))

        if papers:
            return "\n\n".join(
                [f"Title: {paper.title}\nAbstract: {paper.abstract}\nJournal: {paper.journal}" for paper in papers]
            )
        return "No papers found on PubMed after multiple query attempts."
    except Exception as e:
        return f"Error querying PubMed: {e}"


@mcp.tool()
def query_arxiv(query: str, max_papers: int = 10) -> str:
    """Query arXiv for papers based on the provided search query.

    Args:
        query: The search query string.
        max_papers: The maximum number of papers to retrieve (default: 10).

    Returns:
        The formatted search results or an error message.
    """
    import arxiv

    try:
        client = arxiv.Client()
        search = arxiv.Search(query=query, max_results=max_papers, sort_by=arxiv.SortCriterion.Relevance)
        results = "\n\n".join(
            [f"Title: {paper.title}\nSummary: {paper.summary}" for paper in client.results(search)]
        )
        return results if results else "No papers found on arXiv."
    except Exception as e:
        return f"Error querying arXiv: {e}"


@mcp.tool()
def query_scholar(query: str) -> str:
    """Query Google Scholar for papers based on the provided search query.

    Args:
        query: The search query string.

    Returns:
        The first search result formatted or an error message.
    """
    from scholarly import ProxyGenerator, scholarly

    pg = ProxyGenerator()
    pg.FreeProxies()
    scholarly.use_proxy(pg)
    try:
        search_query = scholarly.search_pubs(query)
        result = next(search_query, None)
        if result:
            return (
                f"Title: {result['bib']['title']}\n"
                f"Year: {result['bib']['pub_year']}\n"
                f"Venue: {result['bib']['venue']}\n"
                f"Abstract: {result['bib']['abstract']}"
            )
        return "No results found on Google Scholar."
    except Exception as e:
        return f"Error querying Google Scholar: {e}"


@mcp.tool()
def extract_url_content(url: str) -> str:
    """Extract the text content of a webpage using requests and BeautifulSoup.

    Args:
        url: Webpage URL to extract content from.

    Returns:
        Text content of the webpage.
    """
    response = requests.get(url, headers={"User-Agent": "Mozilla/5.0"})

    content_type = response.headers.get("Content-Type", "")
    if "text/plain" in content_type or "application/json" in content_type:
        return response.text.strip()

    soup = BeautifulSoup(response.text, "html.parser")
    content = soup.find("main") or soup.find("article") or soup.body

    if content is None:
        return response.text.strip()

    for element in content(["script", "style", "nav", "header", "footer", "aside", "iframe"]):
        element.decompose()

    paragraphs = content.find_all(["p", "h1", "h2", "h3", "h4", "h5", "h6"])
    cleaned_text = []
    for p in paragraphs:
        text = p.get_text().strip()
        if text:
            cleaned_text.append(text)

    return "\n\n".join(cleaned_text)


@mcp.tool()
def extract_pdf_content(url: str) -> str:
    """Extract the text content of a PDF file given its URL.

    Args:
        url: URL of the PDF file to extract text from.

    Returns:
        The extracted text content from the PDF.
    """
    try:
        if not url.lower().endswith(".pdf"):
            response = requests.get(url, timeout=30)
            if response.status_code == 200:
                pdf_links = re.findall(r'href=[\'"]([^\'"]+\.pdf)[\'"]', response.text)
                if pdf_links:
                    if not pdf_links[0].startswith("http"):
                        base_url = "/".join(url.split("/")[:3])
                        url = base_url + pdf_links[0] if pdf_links[0].startswith("/") else base_url + "/" + pdf_links[0]
                    else:
                        url = pdf_links[0]
                else:
                    return f"No PDF file found at {url}. Please provide a direct link to a PDF file."

        response = requests.get(url, timeout=30)

        content_type = response.headers.get("Content-Type", "").lower()
        if "application/pdf" not in content_type and not response.content.startswith(b"%PDF"):
            return f"The URL did not return a valid PDF file. Content type: {content_type}"

        pdf_file = BytesIO(response.content)

        text = ""
        try:
            pdf_reader = PyPDF2.PdfReader(pdf_file)
            for page_num in range(len(pdf_reader.pages)):
                page = pdf_reader.pages[page_num]
                text += page.extract_text() + "\n\n"
        except Exception as e:
            return f"Error extracting text from PDF: {str(e)}"

        text = re.sub(r"\s+", " ", text).strip()

        if not text:
            return "The PDF file did not contain any extractable text. It may be an image-based PDF requiring OCR."

        return text

    except requests.exceptions.RequestException as e:
        return f"Error downloading PDF: {str(e)}"
    except Exception as e:
        return f"Error extracting text from PDF: {str(e)}"


if __name__ == "__main__":
    mcp.run(transport="stdio")
