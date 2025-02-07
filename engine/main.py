import streamlit as st

from analysis import SummaryAnalysis
from parse import PdfParser
from source import PaperSource, PubMedSource
import logging
from logging import info

logging.getLogger("httpx").setLevel(logging.ERROR)

st.title("Automated meta-analysis")

with st.form(key="prompt"):
    prompt = st.text_area("Intervention").strip()
    n = st.number_input("Number of papers to process", min_value=1, max_value=200, value=30)
    submit = st.form_submit_button("Submit")

status = st.empty()

class BasicStreamlitLogHandler(logging.StreamHandler):
    def emit(self, record):
        status.write(self.format(record))

logging.basicConfig(handlers=[BasicStreamlitLogHandler(), logging.StreamHandler()], level=logging.INFO)

if submit:
    info(f"Running basic analysis for {prompt}")

    summ = SummaryAnalysis(PubMedSource())

    result = summ.run(prompt, n)
    
    st.write(result)
    
    


