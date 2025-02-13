import pandas as pd
from parse import LLamaParser
from source import ScholarSource
from analysis import DfsAnalysisV1
from dotenv import load_dotenv
from pathlib import Path
import logging
import warnings
import streamlit as st
from logging import info
import json

warnings.filterwarnings("ignore")
logging.getLogger("httpx").setLevel(logging.ERROR)
logging.getLogger("scholarly").setLevel(logging.ERROR)
logging.getLogger("scihub").setLevel(logging.ERROR)

load_dotenv()

def flatten(l):
    if isinstance(l, list):
        return [item for sublist in l for item in flatten(sublist)]
    else:
        return [l]


class BasicStreamlitLogHandler(logging.StreamHandler):
    def emit(self, record):
        if record.levelname == "DEBUG":
            return
        txt = self.format(record)
        status.write(txt)
        status.update(label=txt[:100] + ("..." if len(txt) > 100 else ""))


st.title("Automated meta-analysis")

with st.form(key="prompt"):
    prompt = st.text_area("Intervention").strip()
    n = st.number_input(
        "Number of papers to process", min_value=1, max_value=200, value=10
    )
    submit = st.form_submit_button("Submit")


if submit:
    progress = st.progress(0, "Progress")
    status = st.status("Loading...")

    logging.getLogger().addHandler(BasicStreamlitLogHandler())

    with status:
        info(f"Running basic analysis for {prompt}")

        summ = DfsAnalysisV1(ScholarSource(Path("pdfs/")), LLamaParser(Path("pdfs/cache/")))

        result = summ.run(prompt, n, progress=lambda x: progress.progress(x, "Progress"))

        status.write("Done")

    st.header("Results")

    st.write(result["report"])

    st.dataframe(pd.DataFrame(result["data"]))
