import dataclasses
import itertools
import pathlib
from typing import Optional
import paperscraper
import scihub
from parse import LLamaParser
from data import Paper
from datetime import datetime
import paperscraper.pubmed as pm
import paperscraper.scholar as pscholar
import paperscraper.pdf as pspdf
import curl_cffi.requests as requests
from hashlib import sha256
import json
import jsonpickle

import scholarly
from logging import info
from pathlib import Path


class PaperSource:
    def search(self, query: str, n=100) -> list[Paper]:
        raise NotImplementedError()

    def get(self, paper: Paper) -> Optional[Paper]:
        raise NotImplementedError()


class PdfSource(PaperSource):
    def __init__(self, pdfpath: Path):
        self.pdfpath = pdfpath
        self.pdfpath.mkdir(exist_ok=True, parents=True)

        self.scihub = scihub.SciHub()

        super().__init__()

    def _resolve_doi(self, title: str) -> Optional[str]:
        url = f"https://api.crossref.org/works?query.title={title}&rows=1"
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()

            if (
                "message" in data
                and "items" in data["message"]
                and data["message"]["items"]
            ):
                return data["message"]["items"][0].get("DOI")
        return None

    def get(self, paper: Paper) -> Optional[Paper]:
        if paper.doi is None:
            paper = dataclasses.replace(paper, doi=self._resolve_doi(paper.title))
        
        if paper.doi is None:
            return None

        path = self.pdfpath / f"{paper.doi.replace('/', '--')}.pdf"

        try:
            pspdf.save_pdf({"doi": paper.doi}, str(path))
        except Exception as e:
            info(f"Failed to download PDF via PaperScraper for {paper.doi}: {e}")

        try:
            res = self.scihub.fetch(paper.doi)
            if not res:
                raise Exception("Failed to fetch PDF via SciHub")

            pdf: bytes = res["pdf"]  # type: ignore

            path.write_bytes(pdf)
        except Exception as e:
            info(f"Failed to download PDF via SciHub for {paper.doi}: {e}")

        if not path.exists():
            return None

        paper = dataclasses.replace(paper, pdf=path)

        return paper


class ScholarSource(PdfSource):
    def search(self, query, n=100) -> list[Paper]:
        info(f"Searching Google Scholar for {query}")
        (self.pdfpath / "scholar_cache").mkdir(exist_ok=True, parents=True)

        hs = sha256((query + str(n)).encode(), usedforsecurity=False).hexdigest()[:32]

        if (self.pdfpath / "scholar_cache" / hs).exists():
            with open(self.pdfpath / "scholar_cache" / hs) as f:
                return jsonpickle.decode(f.read())

        results = scholarly.scholarly.search_pubs(query)

        ans = [
            Paper(
                p["bib"]["title"],
                p["bib"]["abstract"],
                p["bib"]["venue"],
                (
                    datetime(year=int(p["bib"]["pub_year"]), month=1, day=1)
                    if "pub_year" in p["bib"] and p["bib"]["pub_year"].isdigit()
                    else None
                ),
                # None,
                None,
                None,
                None,
            )
            for p in itertools.islice(results, n)
        ]

        with open(self.pdfpath / "scholar_cache" / hs, "w") as f:
            f.write(jsonpickle.encode(ans, unpicklable=True))

        return ans


class PubMedSource(PdfSource):
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
