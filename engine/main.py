import streamlit as st

from analysis import SummaryAnalysis
from parse import PdfParser
from source import PaperSource, PubMedSource, ScholarSource
import logging
from logging import info
from dotenv import load_dotenv

load_dotenv()

logging.getLogger("httpx").setLevel(logging.ERROR)
logging.getLogger("scholarly").setLevel(logging.ERROR)

st.title("Automated meta-analysis")

with st.form(key="prompt"):
    prompt = st.text_area("Intervention").strip()
    n = st.number_input(
        "Number of papers to process", min_value=1, max_value=200, value=30
    )
    submit = st.form_submit_button("Submit")


if submit:
    status = st.status("Loading...")

    class BasicStreamlitLogHandler(logging.StreamHandler):
        def emit(self, record):
            txt = self.format(record)
            status.write(txt)
            status.update(label=txt[:100] + ("..." if len(txt) > 100 else ""))

    logging.getLogger().addHandler(BasicStreamlitLogHandler())

    with status:
        info(f"Running basic analysis for {prompt}")

        summ = SummaryAnalysis(ScholarSource())

        result = summ.run(prompt, n)

        status.write("Done")

    st.header("Results")

    st.write(result)
