from parse import LLamaParser
from source import ScholarSource
from analysis import DfsAnalysisV1
from dotenv import load_dotenv
from pathlib import Path
import logging
import warnings

warnings.filterwarnings("ignore")
logging.getLogger("httpx").setLevel(logging.ERROR)
logging.getLogger("scholarly").setLevel(logging.ERROR)
logging.getLogger("scihub").setLevel(logging.ERROR)

load_dotenv()

da1 = DfsAnalysisV1(ScholarSource(Path("pdfs/")), LLamaParser(Path("pdfs/cache/")))

res = da1.run(input("Query: "))

print(res)