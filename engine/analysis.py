from abc import ABC
from os import environ
from re import A
from typing import Any, Optional, TypeVar

from attr import dataclass

from data import Paper
from source import PaperSource
from gigachat import GigaChat
import json
from logging import debug, info, warning

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
        prompt = (
            self.prologue
            + "You are researching the effects of {query} on human longevity. \n"
        )
        prompt += f"You need to filter the papers you have found. \n"
        prompt += f"Reply to the following question with one word: Yes or No. Other responses are forbidden. \n"
        prompt += (
            f"Does this paper investigate the effects of {query} on human longevity? \n"
        )
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
        prompt += "\n".join(
            [
                f"{i+1}. {paper.title}\nAbstract: {paper.abstract}\n\n"
                for i, paper in enumerate(data)
            ]
        )

        prompt += "You must respond in the aforementioned format. \n"

        info(f"Summarizing papers")

        response = self.model.chat(prompt)

        return response.choices[0].message.content

    def run(self, query: str, n: int = 20):
        info(f"Running basic analysis for {query}")

        orig_papers = self.source.search(
            f"({query}) AND (longevity OR mortality OR lifespan) AND (RCT OR trial OR study) AND (human OR man OR woman)",
            n,
        )
        papers = self.filter(query, orig_papers)
        info(f"Accepted {len(papers)}/{len(orig_papers)} papers")

        if len(papers) == 0:
            return "No papers found"

        summary = self.summarize(query, papers)
        return summary


class Question(ABC):
    RETRIES = 3

    R = TypeVar("R")

    def __init__(self, q_id: str, prompter):
        self.q_id = q_id
        self.prompter = prompter
        self.log_asked = None
        self.log_res = None

    def ask(self, data: dict) -> str:
        pass

    def answer(self, data: dict, result: str) -> Optional[R]:
        debug(f"Answered question {self.q_id}: {result}")

        result = result.strip()

        if "</thought>" in result:
            result = result.split("</thought>")[1].strip()

        return result

    def branch(self, data: dict, result: str) -> list[dict]:
        return [{self.q_id: result, **data}]

    def execute(self, llm, data) -> Optional[R]:
        self.log_asked = self.ask(data)
        response = llm(self.log_asked)
        for _ in range(self.RETRIES):
            try:
                self.log_res = self.answer(data, response)
                return self.log_res
            except Exception as e:
                warning("Error in answering question", e)
        raise e

    def dfs(self, llm, path: list[Question], data: dict) -> list[dict]:
        response = self.execute(llm, data)

        ans = []
        for q in self.branch(data, response):
            if not path:
                ans.append(q)
            else:
                ans.append(path[0].dfs(llm, path[1:], q))

        return ans


class Choice(Question):
    def __init__(self, q_id: str, prompter, categories: list[str]):
        super().__init__(q_id, prompter)

        assert "none" not in categories
        self.categories = categories


class SingleChoice(Choice):
    def ask(self, data) -> str:
        return (
            f"{self.prompter(data)}\n"
            "You must first provide the reasoning behind your choice in the following format:\n"
            "<thought> YOUR THOUGTHS HERE </thought>\n"
            "Then you must answer with exactly one of the following options on a single line and nothing else: \n"
            f"{', '.join(self.categories)}\n"
            "If neither of the options is applicable, you may respond with 'None'.\n"
        )

    def answer(self, data, result: str) -> Optional[str]:
        result = super().answer(data, result)

        if result == "none":
            return None

        if result not in self.categories:
            raise ValueError("Invalid response")

        return result


class AnyList(Question):
    def ask(self, data) -> str:
        return (
            f"{self.prompter(data)}\n"
            "You must first provide the reasoning behind your choice in the following format:\n"
            "<thought> YOUR THOUGTHS HERE </thought>\n"
            "Then you must answer with ONE OR MORE options on a single line and nothing else: \n"
            f"{', '.join(self.categories)}\n"
            "You must answer in JSON format, i.e. ['Option 1', 'Option 2']\n"
            "If neither of the options is applicable, you may respond with [].\n"
        )

    def branch(self, data, result):
        return [{self.q_id: x, **data} for x in result]

    def answer(self, data, result: str) -> Optional[str]:
        result = super().answer(data, result)

        try:
            l = json.loads(result)
        except json.JSONDecodeError:
            raise ValueError("Invalid response format")

        return l


class MultiChoice(Choice):
    def ask(self, data) -> str:
        return (
            f"{self.prompter(data)}\n"
            "You must first provide the reasoning behind your choice in the following format:\n"
            "<thought> YOUR THOUGTHS HERE </thought>\n"
            "Then you must answer with ONE OR MORE of the following options on a single line and nothing else: \n"
            f"{', '.join(self.categories)}\n"
            "You must answer in JSON format, i.e. ['option1', 'option2']\n"
            "If neither of the options is applicable, you may respond with [].\n"
        )

    def branch(self, data, result):
        return [{self.q_id: x, **data} for x in result]

    def answer(self, data, result: str) -> Optional[str]:
        result = super().answer(data, result)

        try:
            l = json.loads(result)
        except json.JSONDecodeError:
            raise ValueError("Invalid response format")

        if any([x not in self.categories for x in l]):
            raise ValueError("Invalid response value")

        return l


class AnyText(Question):
    def ask(self, data) -> str:
        return (
            f"{self.prompter(data)}\n"
            "You must first provide the reasoning behind your choice in the following format:\n"
            "<thought> YOUR THOUGTHS HERE </thought>\n"
            "Then you must write your answer on a single line and nothing else.\n"
            "If you are unable to provide a truthful answer, you may respond with 'None'.\n"
        )

    def answer(self, data, result: str) -> Optional[str]:
        result = super().answer(data, result)

        if result == "none":
            return None

        return result


class Numeric(Question):
    def ask(self, data) -> str:
        return (
            f"{self.prompter(data)}\n"
            "You must first provide the reasoning behind your choice in the following format:\n"
            "<thought> YOUR THOUGTHS HERE </thought>\n"
            "Then you must answer with a single number (or a decimal fraction) on a single line and nothing else.\n"
            "If you are unable to provide a number, you may respond with 'None'.\n"
        )

    def answer(self, data, result: str) -> Optional[float]:
        result = super().answer(data, result)

        result = result.replace(",", ".")

        if result == "none":
            return None

        try:
            return float(result)
        except ValueError:
            raise ValueError("Invalid response")


class DfsAnalysisV1:
    def __init__(self, source: PaperSource):
        self.model = GigaChat(
            credentials=environ["GIGACHAT_API_KEY"],
            scope="GIGACHAT_API_CORP",
            model="GigaChat-Pro",
            verify_ssl_certs=False,
        )
        self.source = source
        self.categories = {
            "species": [
                "mice",
                "rats",
                "humans",
                "non-human primates",
                "drosophila",
                "bacteria",
            ],
        }

    def ask_llm(self, query: str) -> str:
        return self.model.chat(query).choices[0].message.content
    
    def run(self, query, n: int = 10):
        papers = self.source.search("({query}) AND (longevity OR mortality OR lifespan) AND (RCT OR trial OR study)")

    def process_paper(self, query, paper: Paper, n: int = 20):
        info(f"Processing {paper.title}")

        assert paper.full_text is not None

        prologue = (
            "You are an expert scientist. \n"
            f"You are researching the impact of {query} on longevity. \n"
            "You need to analyze the following paper. \n"
            f"Paper: {paper.title}\nAbstract: {paper.abstract}\n\n"
            f"{paper.full_text}\n\n"
        )

        q = {}

        q[0] = MultiChoice(
            "species",
            lambda d: prologue
            + (
                f"Which species is the study conducted on?"
                "List all species that were subject to interventions in the study."
            ),
            self.categories["species"],
        )

        q[1] = AnyList(
            "intervention",
            lambda d: prologue
            + (
                f"What interventions are studied in groups consisting of {d['species']}?\n"
                "List one for each group, including the control group (if applicable).\n"
                "Include groups that are subjected to combinations of interventions as well (if applicable)"
            ),
        )

        q[2] = Numeric(
            "n",
            lambda d: prologue
            + f"How many subjects are in the '{d['intervention']}' group of {d['species']}?",
        )

        q[3] = AnyList(
            "target",
            lambda d: prologue
            + f"What outcomes (targets) are measured in the '{d['intervention']}' group of {d['species']}" + f" (N = {d['n']})?" if d['n'] else "",
        )

        q[4] = AnyText(
            "result",
            lambda d: prologue
            + (
                f"What is the final result (impact) on the '{d['target']}' target in the '{d['intervention']}' group of {d['species']} (N = {d['n']}) as the result of the study?\n"
                "Provide a clear and concise answer, including the direction of the effect (if applicable).\n"
                "If the results are inconclusive, please state so."
            ),
        )

        q[5] = AnyText(
            "evidence",
            lambda d: prologue
            + (
                f"What is the evidence supporting the result on the '{d['target']}' target in the '{d['intervention']}' group of {d['species']} (N = {d['n']})?\n"
                f"The result is {d['result']}.\n"
                "Provide numerical and other citations from the abstracts justifying EXACTLY THIS conclusion."
            ),
        )

        q[6] = AnyText(
            "p_value",
            lambda d: prologue
            + (
                f"What is the p-value of the result on the '{d['target']}' target in the '{d['intervention']}' group of {d['species']} (N = {d['n']})?\n"
                f"The result is {d['result']}, as justified by the following: {d['evidence']}\n"
                "Provide the p-value of THIS SPECIFIC result if it is stated in the paper"
            ),
        )

        q = list(q.values())

        result = q[0].dfs(self.ask_llm, q[1:], {})

        return result
