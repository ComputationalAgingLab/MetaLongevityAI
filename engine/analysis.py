from os import environ
import re
from data import Paper
from source import PaperSource
from gigachat import GigaChat
import json
from logging import info

# response = model.chat("Расскажи о себе в двух словах?")
# print(response.choices[0].message.content)


class SummaryAnalysis:
    def __init__(self, source: PaperSource):
        self.model = GigaChat(
            credentials=environ["GIGACHAT_API_KEY"],
            scope="GIGACHAT_API_CORP",
            model="GigaChat-Pro",
            verify_ssl_certs=False,
        )
        self.source = source
        self.prologue = f"You are an expert scientist. \n"
    
    def filter_one(self, query: str, paper: Paper) -> bool:
        prompt = self.prologue + "You are researching the effects of {query} on human longevity. \n"
        prompt += f"You need to filter the papers you have found. \n"
        prompt += f"Reply to the following question with one word: Yes or No. Other responses are forbidden. \n"
        prompt += f"Does this paper investigate the effects of {query} on human longevity? \n"
        prompt += f"Paper: {paper.title}\nAbstract: {paper.abstract}\n\n"

        info(f"Evaluating {paper.title}")

        response = self.model.chat(prompt)

        resp_text = response.choices[0].message.content

        result = "yes" in resp_text.strip().lower()

        if not result:
            info(f"Rejected with response: {resp_text}")
        
        return result

    def filter(self, query: str, data: list[Paper]) -> list[Paper]:
        info(f"Filtering {len(data)} papers")
        return [paper for paper in data if self.filter_one(query, paper)]

    def summarize(self, query: str, data: list[Paper]):
        prompt = self.prologue
        prompt += f"You are researching the effects of {prompt} on human longevity. "
        prompt += f"You need to summarize the findings present in the papers you have found. \n"
        prompt += f"Present a conclusion based on the data (magnitude and direction of effect, i.e. inconclusive, longevity increased or decreased) "
        prompt += f"Make sure to support your findings with evidence from the papers. "
        prompt += f"You must present your findings in the following fixed format: \n"
        prompt += f"Thought process: <clearly explain your reasoning> \n"
        prompt += f"Conclusion: {query} has <increased|decreased|no effect> human longevity. \n"
        prompt += f"Evidence: <provide numerical and other citations from the abstracts justifying the conclusion> \n"
        prompt += f"Here are the papers you have found: \n"
        prompt += "\n".join([f"{i+1}. {paper.title}\nAbstract: {paper.abstract}\n\n" for i, paper in enumerate(data)])

        prompt += "You must respond in the aforementioned format. \n"

        info(f"Summarizing papers")

        response = self.model.chat(prompt)

        return response.choices[0].message.content
    
    def run(self, query: str, n: int = 20):
        info(f"Running basic analysis for {query}")

        orig_papers = self.source.search(f"({query}) AND (longevity OR mortality OR lifespan) AND (RCT OR trial OR study) AND (human OR man OR woman)", n)
        papers = self.filter(query, orig_papers)
        info(f"Accepted {len(papers)}/{len(orig_papers)} papers")

        if len(papers) == 0:
            return "No papers found"
        
        summary = self.summarize(query, papers)
        return summary

class StudyAnalysis:
    def __init__(self, source: PaperSource):
        self.model = GigaChat(
            credentials=environ["GIGACHAT_API_KEY"],
            scope="GIGACHAT_API_CORP",
            model="GigaChat-Pro",
            verify_ssl_certs=False,
        )
        self.source = source
        self.prologue = f"You are an expert scientist. \n"
    
    def filter_one(self, query: str, paper: Paper) -> bool:
        prompt = self.prologue + "You are researching the effects of {query} on human longevity. \n"
        prompt += f"You need to filter the papers you have found. \n"
        prompt += f"Reply to the following question with one word: Yes or No. Other responses are forbidden. \n"
        prompt += f"Does this paper present a clinical study on the effects of {query} on human longevity? \n"
        prompt += f"Paper: {paper.title}\nAbstract: {paper.abstract}\n\n"

        info(f"Evaluating {paper.title}")

        response = self.model.chat(prompt)

        resp_text = response.choices[0].message.content

        result = "yes" in resp_text.strip().lower()

        if not result:
            info(f"Rejected with response: {resp_text}")
        
        return result

    def filter(self, query: str, data: list[Paper]) -> list[Paper]:
        info(f"Filtering {len(data)} papers")
        return [paper for paper in data if self.filter_one(query, paper)]

    def summarize(self, query: str, data: list[Paper]):
        prompt = self.prologue
        prompt += f"You are researching the effects of {prompt} on human longevity. "
        prompt += f"You need to summarize the findings present in the papers you have found. \n"
        prompt += f"Present a conclusion based on the data (magnitude and direction of effect, i.e. inconclusive, longevity increased or decreased) "
        prompt += f"Make sure to support your findings with evidence from the papers. "
        prompt += f"You must present your findings in the following fixed format: \n"
        prompt += f"Thought process: <clearly explain your reasoning> \n"
        prompt += f"Conclusion: {query} has <increased|decreased|no effect> human longevity. \n"
        prompt += f"Evidence: <provide numerical and other citations from the abstracts justifying the conclusion> \n"
        prompt += f"Here are the papers you have found: \n"
        prompt += "\n".join([f"{i+1}. {paper.title}\nAbstract: {paper.abstract}\n\n" for i, paper in enumerate(data)])

        prompt += "You must respond in the aforementioned format. \n"

        info(f"Summarizing papers")

        response = self.model.chat(prompt)

        info(response)

        return response.choices[0].message.content
    
    def run(self, query: str, n: int = 20):
        info(f"Running basic analysis for {query}")

        orig_papers = self.source.search(f"({query}) AND (longevity OR mortality OR lifespan) AND (RCT OR trial OR study) AND (human OR man OR woman)", n)
        papers = self.filter(query, orig_papers)
        info(f"Accepted {len(papers)}/{len(orig_papers)} papers")

        if len(papers) == 0:
            return "No papers found"
        
        summary = self.summarize(query, papers)
        return summary
