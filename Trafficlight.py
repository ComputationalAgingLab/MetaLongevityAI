class TrafficLight:
    def __init__(self):
        self.q_weights = {'Q1': 4, 'Q2': 3, 'Q3': 2, 'Q4': 1}
        self.art_typ_weights = {'SystematicReview': 5.00,
                                'Meta-Analysis':
                                    5.00,
                                'RandomizedControlledTrial':
                                    4.50,
                                'ClinicalTrial':
                                    4.00,
                                'ReviewArticle':
                                    3.50,
                                'CaseReports':
                                    2.50,
                                'JournalArticle':
                                    3.00,
                                'Editorial':
                                    1.50,
                                'Letter':
                                    1.20,
                                'Commentary':
                                    1.00,
                                'Review':
                                    1.00,
                                'None':
                                    0.00
                                }

    def calculate(self, j_q, art_typ):
        cur = self.q_weights[j_q] * self.art_typ_weights[art_typ]
        if cur >= 0 and cur <= 5:
          return "Red"
        elif cur > 5 and cur < 12:
          return "Yellow"
        else:
          return "Green"


import requests

def get_publication_type_from_semantic_scholar(doi):
    url = f"https://api.semanticscholar.org/graph/v1/paper/DOI:{doi}?fields=title,publicationTypes"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        return data.get("publicationTypes", ["None"]) if data.get("publicationTypes", ["None"]) else ["None"]
    return Exception

from Bio import Entrez
from Bio.Entrez import read

Entrez.email = ("example@example.com")


def get_pmid_by_doi(doi):
    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "pubmed",
        "term": f"{doi}[doi]",
        "retmode": "json",
        "retmax": 1
    }

    response = requests.get(search_url, params=params)
    response.raise_for_status()
    data = response.json()

    pmid_list = data.get("esearchresult", {}).get("idlist", [])
    return pmid_list[0] if pmid_list else ["None"]


def get_publication_type_from_pubmed(pmid):
    handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
    records = read(handle)

    try:
        article = records["PubmedArticle"][0]["MedlineCitation"]["Article"]
        pub_types = article.get("PublicationTypeList", [])
        pub_types = (','.join(pub_types)).split(',')
        for i in range(len(pub_types)):
            pub_types[i] = ''.join(pub_types[i].split(' '))
        return pub_types if pub_types else ["None"]
    except Exception as e:
        return ["None"]

def GetPublicationType(doi):
  try:
    return get_publication_type_from_semantic_scholar(doi)
  except Exception:
    try:
        return get_publication_type_from_pubmed(get_pmid_by_doi(doi))
    except Exception as e:
      return ["None"]

def GetJournalByTitle(title):
    url = f"https://api.crossref.org/works"
    params = {
        'query.title': title,
        'rows': 1  # To limit the response to the first match (adjust as needed)
    }
    response = requests.get(url, params=params)
    if response.status_code == 200:
        data = response.json()
        if data['message']['items']:
            item = data['message']['items'][0]  # Take the first result
            journal_title = item.get('container-title', ['Not Found'])[0]  # Get the journal title
            return journal_title
        else:
            return "No results found for the given title."
    else:
        return f"Error: {response.status_code}"


import pandas as pd
from fuzzywuzzy import process


def GetQ(title: str) -> str:
    title = GetJournalByTitle(title)
    df = pd.read_excel("data.xlsx")
    df["Квартиль"] = pd.qcut(df["Cites / Doc. (2years)"], q=4, labels=["Q1", "Q2", "Q3", "Q4"])
    matches = process.extract(title, df['Title'], limit=1)
    result = df[df["Title"] == matches[0][0]]
    return result.iat[0, 23]


L = TrafficLight()
print(GetQ("Hesperetin promotes longevity and delays aging via activation of Cisd2 in naturally aged mice"))
print(GetPublicationType("10.1186/s12866-019-1494-7"))
print(L.calculate(GetQ("Hesperetin promotes longevity and delays aging via activation of Cisd2 in naturally aged mice"), GetPublicationType("10.1186/s12929-022-00838-7")[0]))