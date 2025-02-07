from data import Paper
from datetime import datetime
import paperscraper.pubmed as pm
from logging import info

class PaperSource:
    def search(self, query: str, n=10) -> list[Paper]:
        raise NotImplementedError()

    def get(self, doi: str) -> Paper:
        raise NotImplementedError()

class PubMedSource(PaperSource):
    def search(self, query, n=100) -> list[Paper]:
        info(f"Searching PubMed for {query}")
        res = pm.get_pubmed_papers(query, max_results=n)

        return [
            Paper(
                p["title"],
                p["abstract"],
                p["journal"],
                datetime.fromisoformat(p["date"]),
                p["doi"],
                None,
                None,
            )
            for _, p in res.iterrows()
        ]
