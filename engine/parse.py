from logging import info, warning
from os import environ
from typing import Optional
from data import Paper
from pathlib import Path
from llama_cloud_services import LlamaParse
from llama_cloud_services.parse.utils import ResultType as LLamaResult
from hashlib import sha256


class PdfParser:
    def parse(self, path: Path) -> Optional[str]:
        raise NotImplementedError()


class LLamaParser(PdfParser):
    def __init__(self, cachepath: Path):
        self.parser = LlamaParse(
            api_key=environ["LLAMAPARSE_API_KEY"], result_type=LLamaResult.MD, language="en"
        )
        self.cache = cachepath

    def parse(self, path: Path) -> Optional[str]:
        info(f"Parsing {path}")

        if not path.exists():
            warning(f"Tried to parse non-existent file {path}")
            return None

        hs = sha256(path.read_bytes(), usedforsecurity=False).hexdigest()[:32]

        if (self.cache / hs).exists():
            info(f"Cache hit for {path}")
            return (self.cache / hs).read_text()

        result = "\n\n".join(x.text for x in self.parser.load_data(str(path)))

        self.cache.mkdir(exist_ok=True, parents=True)
        (self.cache / hs).write_text(result)

        return result