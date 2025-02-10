from logging import info
from os import environ
from typing import Optional
from data import Paper
from pathlib import Path
from llama_cloud_services import LlamaParse


class PdfParser:
    def parse(self, path: str) -> str:
        raise NotImplementedError()


class LLamaParser(PdfParser):
    def __init__(self):
        self.parser = LlamaParse(
            api_key=environ["LLAMAPARSE_API_KEY"], result_type="markdown", language="en"
        )

    def parse(self, path: Path) -> Optional[str]:
        info(f"Parsing {path}")
        return next((x.text for x in self.parser.load_data(path)), None)
