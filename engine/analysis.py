from datetime import datetime
import os
from abc import ABC
from os import environ, error
from re import A
from typing import Any, Optional, TypeVar

from attr import dataclass
from numpy import source

from description import PaperDescriptionSource
from trafficlight import GetQ, GetPublicationType, TrafficLightClassifier, get_doi
from parse import LLamaParser, PdfParser
from data import Paper
from source import PaperSource
from gigachat import GigaChat
import json
from logging import debug, info, warning
from pathlib import Path

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
                f"{i + 1}. {paper.title}\nAbstract: {paper.abstract}\n\n"
                for i, paper in enumerate(data)
            ]
        )

        prompt += "You must respond in the aforementioned format. \n"

        info(f"Summarizing papers")

        response = self.model.chat(prompt)

        return response.choices[0].message.content

    def classifier(self, data: list[Paper]):
        dct = {"Red": 0, "Green": 0, "Yellow": 0}
        for paper in data:
            tl = TrafficLightClassifier()
            paper.color = tl.classify(title=paper.title)
            dct[paper.color] += 1
            print(paper.doi, paper.color)
        return dct

    def run(self, query: str, n: int = 20):
        info(f"Running basic analysis for {query}")

        orig_papers = self.source.search(
            f"({query}) AND (longevity OR mortality OR lifespan) AND (RCT OR trial OR study) AND (human OR man OR woman)",
            n,
        )
        papers = self.filter(query, orig_papers)
        info(f"Accepted {len(papers)}/{len(orig_papers)} papers")
        dct = self.classifier(papers)
        info(f"Classifying success {dct}")
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
        raise NotImplementedError()

    def answer(self, data: dict, result: str) -> Any:
        info(f"Answered question {self.q_id}: {result}")

        result = result.strip()

        REMOVE_BEFORE = [
            "thought>",
            "answer**:",
            "answer:**",
            "result**:",
            "result:**",
            "answer:",
            "result:",
        ]

        for i in REMOVE_BEFORE:
            if i in result.lower():
                result = result[result.lower().rfind(i) + len(i) :]

        return result.strip()

    def branch(self, data: dict, result: Optional[str]):
        return [{self.q_id: result, **data}]

    def execute(self, llm, data) -> Optional[str]:
        for _ in range(self.RETRIES):
            try:
                self.log_asked = self.ask(data)
                response = llm(self.log_asked)
                self.log_res = self.answer(data, response)
                return self.log_res
            except Exception as e:
                warning(f"Error in answering question {e}")
        raise e  # type: ignore

    def dfs(self, llm, path: list, data: dict) -> list[dict]:
        response = self.execute(llm, data)

        ans = []
        children = self.branch(data, response)
        for q in children:
            if not path:
                ans.append(q)
            else:
                ans += path[0].dfs(llm, path[1:], q)

        if children == []:
            raise ValueError("No children while evaluating path")

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
            "You must select exactly one option of the following:\n"
            f"{', '.join(self.categories)}\n"
            'If neither of the options is applicable, you may respond with "None".\n'
            "You must answer in the following format:\n"
            "**Reasoning:** YOUR REASONING HERE\n"
            '**Answer:** YOUR ANSWER HERE (i.e "**Answer:** Option 1" and NOTHING ELSE)\n'
            # "You must first provide the reasoning behind your choice in the following format:\n"
            # "<thought> YOUR THOUGTHS HERE </thought>\n"
            # "Then you must answer with exactly one of the following options on a single line and nothing else: \n"
            # f"{', '.join(self.categories)}\n"
            # "If neither of the options is applicable, you may respond with 'None'.\n"
        )

    def answer(self, data, result: str):
        res = [x for x in super().answer(data, result).split("\n") if x.strip() != ""][
            -1
        ].strip()

        if res is None or res.lower() == "none":
            return None

        if res not in self.categories:
            raise ValueError("Invalid response")

        return res


class AnyList(Question):
    def ask(self, data) -> str:
        return (
            f"{self.prompter(data)}\n"
            'You must provide the answer as a list of double-quoted strings encased with square braces i.e. ["Option 1", "Option 2"] or ["Option 1"]\n'
            "If neither of the options is applicable, you may respond with [].\n"
            "You must answer in the following format:\n"
            "**Reasoning:** YOUR REASONING HERE\n"
            '**Answer:** YOUR ANSWER HERE (i.e ["Option 1", "Option 2"] and NOTHING ELSE)\n'
            # "You must first provide the reasoning behind your choice in the following format:\n"
            # "<thought> YOUR THOUGTHS HERE </thought>\n"
            # "Then you must answer with ONE OR MORE options on a single line and nothing else. \n"
            # 'You must provide them as a list of double-quoted strings encased with square braces i.e. ["Option 1", "Option 2"] or ["Option 1"]\n'
            # "If neither of the options is applicable, you may respond with [].\n"
        )

    def branch(self, data, result):
        return (
            [{self.q_id: x, **data} for x in result]
            if result
            else [{self.q_id: None, **data}]
        )

    def answer(self, data, result: str) -> list[str]:
        res = [x for x in super().answer(data, result).split("\n") if x.strip() != ""][
            -1
        ].strip()

        if res is None:
            return res

        try:
            info("Parsing response: " + res)
            l = json.loads(res.replace("'", '"'))
        except json.JSONDecodeError:
            raise ValueError("Invalid response format")

        return l


class MultiChoice(Choice):
    def __init__(self, q_id: str, prompter, categories: list[str]):
        super().__init__(q_id, prompter, categories)

    def ask(self, data) -> str:
        return (
            f"{self.prompter(data)}\n"
            "You must answer with ZERO OR MORE of the following options: \n"
            f"{', '.join(self.categories)}\n"
            'You must provide them as a list of double-quoted strings encased with square braces i.e. ["Option 1", "Option 2"] or ["Option 1"]\n'
            "If neither of the options is applicable, you may respond with [].\n"
            "You must answer in the following format:\n"
            "**Reasoning:** YOUR REASONING HERE\n"
            '**Answer:** YOUR ANSWER HERE (i.e ["Option 1", "Option 2"] and NOTHING ELSE)\n'
            # "You must first provide the reasoning behind your choice in the following format:\n"
            # "<thought> YOUR THOUGTHS HERE </thought>\n"
            # "Then you must answer with ONE OR MORE of the following options on a single line and nothing else: \n"
            # f"{', '.join(self.categories)}\n"
            # 'You must provide them as a list of double-quoted strings encased with square braces i.e. ["Option 1", "Option 2"] or ["Option 1"]\n'
            # "If neither of the options is applicable, you may respond with [].\n"
        )

    def branch(self, data, result):
        return (
            [{self.q_id: x, **data} for x in result]
            if result
            else [{self.q_id: None, **data}]
        )

    def answer(self, data, result: str) -> list[str]:
        res = [x for x in super().answer(data, result).split("\n") if x.strip() != ""][
            -1
        ].strip()

        if res is None:
            return res

        try:
            info("Parsing response: " + res)
            if res in self.categories:
                return [res]

            l = json.loads(res.replace("'", '"'))
        except json.JSONDecodeError:
            raise ValueError("Invalid response format")

        if any([x not in self.categories for x in l]):
            raise ValueError("Invalid response value")

        return l


class AnyText(Question):
    def ask(self, data) -> str:
        return (
            f"{self.prompter(data)}\n"
            "If you are unable to provide a truthful answer, you may respond with 'None'.\n"
            "You must answer in the following format:\n"
            "**Reasoning:** YOUR REASONING HERE\n"
            "**Answer:** YOUR ANSWER HERE\n"
            # "You must first provide the reasoning behind your choice in the following format:\n"
            # "<thought> YOUR THOUGTHS HERE </thought>\n"
            # "Then you must write your answer\n"
            # "If you are unable to provide a truthful answer, you may respond with 'None'.\n"
        )

    def answer(self, data, result: str) -> Optional[str]:
        res = super().answer(data, result).strip()

        info("Got response: " + res)

        if res.lower() == "none":
            return None

        return res


class Numeric(Question):
    def ask(self, data) -> str:
        return (
            f"{self.prompter(data)}\n"
            "You must answer with a single number (or a decimal fraction).\n"
            "If you are unable to provide a number, you may respond with 'None'.\n"
            "You must answer in the following format:\n"
            "**Reasoning:** YOUR REASONING HERE\n"
            "**Answer:** YOUR ANSWER HERE (A SINGLE NUMBER/DECIMAL FRACTION AND NOTHING ELSE)\n"
        )

    def answer(self, data, result: str) -> Optional[float]:
        res = [x for x in super().answer(data, result).split("\n") if x.strip() != ""][
            -1
        ].strip()

        if res is None:
            return None

        res = res.replace(",", "")

        if res.lower() == "none":
            return None

        try:
            info("Parsing response: " + res)
            return float(res)
        except ValueError:
            raise ValueError("Invalid response")


class DfsAnalysisV1:
    def __init__(self, source: PaperSource, parser: PdfParser):
        self.model = GigaChat(
            credentials=environ["GIGACHAT_API_KEY"],
            scope="GIGACHAT_API_CORP",
            model="GigaChat-Pro",
            verify_ssl_certs=False,
        )

        self.cache_path = "llm_cache.json"
        try:
            self.cache = json.load(open(self.cache_path))
        except:
            self.cache = {}

        self.tl_classifier = TrafficLightClassifier()
        self.desc_source = PaperDescriptionSource(self.model)
        self.source = source
        self.parser = parser
        self.categories = {
            "species": [
                "mice",
                "rats",
                "humans",
                "non-human primates",
                "drosophila",
                "bacteria",
            ],
            "targets": [
                "lifespan/mortality",
                "age-related diseases",
                "biomarkers",
            ],
        }
    
    # @cachier() # sus
    def ask_llm(self, query: str) -> str:
        info("Asking LLM: ...\n" + "\n".join(query.split("\n")[-20:]))

        if query in self.cache:
            info("Cache hit")
            return self.cache[query]

        res = self.model.chat(query).choices[0].message.content
        self.cache[query] = res

        json.dump(self.cache, open(self.cache_path, "w"))

        return res

    def run(self, query, n: int = 10, progress=None):
        if progress is None:
            progress = lambda x: None

        if not self.desc_source.fool_check(query):
            return {"report": "Invalid query", "data": []}

        desc = self.desc_source.get_description(query)

        papers = self.source.search(
            f"({query}) AND (longevity OR mortality OR lifespan) AND (RCT OR trial OR study)",
            n=n,
        )
        info(f"Found {len(papers)} papers")

        results = []
        rows = []

        left = n

        for paper in papers:
            progress(0.8 * (1.0 - left / n))

            left -= 1
            info(f"Fetching pdf for {paper.title}")
            paper = self.source.get(paper)

            if paper is None or (paper.full_text is None and paper.pdf is None):
                warning(
                    f"Failed to fetch pdf for {paper.title if paper else '<unknown>'}"
                )
                continue

            if paper.full_text is None:
                info(f"Parsing pdf for {paper.title}")
                paper.full_text = self.parser.parse(paper.pdf)  # type: ignore

            info("Classifying...")
            paper.color = self.tl_classifier.classify(paper.title)

            res = self.process_paper(query, paper)
            info(f"Found: " + "\n\n".join([json.dumps(x, indent=2) for x in res]))
            results.append({"result": res, "paper": paper})
            rows += [
                {**{f"paper_{k}": v for k, v in paper.__dict__.items() if type(v) in [str, datetime, int, float]}, **x}
                for x in res
            ]

        progress(0.8)

        report = (
            "== Query: " + query + "\n" +
            "Description: \n" + desc + "\n" + self.build_report(query, results)
        )

        conclusion = (
            "== Conclusion: \n" +
            self.make_conclusion(desc, report)
        )

        progress(1.0)

        return {
            "report": report,
            "data": rows,
            "conclusion": conclusion,
        }

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
                f"Which species is the study conducted on?\n"
                "List all species that were subject to interventions in the study.\n"
            ),
            self.categories["species"],
        )

        q[1] = AnyList(
            "intervention",
            lambda d: prologue
            + (
                f"What interventions are studied in groups consisting of {d['species']}?\n"
                "List each intervention group, NOT including the control group.\n"
                "Include groups that are subjected to combinations of interventions as well (if applicable)\n"
            ),
        )

        group_template = lambda d: prologue + (
            f"You need to provide data for the following intervention group: \n"
            f"Species: {d['species']}\n"
            f"Intervention: {d['intervention']}\n"
            + (f"Measured outcome/target: {d['target']}\n" if d.get("target") else "")
            + (f"Result: {d['result']}\n" if d.get("result") else "")
            + (f"Evidence: {d['evidence']}\n" if d.get("evidence") else "")
            + (f"p-value: {d['p_value']}\n" if d.get("p_value") else "")
        )

        q[2] = Numeric(
            "n",
            lambda d: group_template(d) + f"How many subjects are in this group?\n",
        )

        q[3] = AnyList(
            "target",
            lambda d: group_template(d)
            + (
                f"What outcomes (targets) are measured in this group?\n"
                f"You may include values like all-cause mortality, lifespan or effect on specific biomarkers or age-related diseases\n"
            ),
        )

        q[4] = AnyText(
            "result",
            lambda d: group_template(d)
            + (
                f"What is the result of the intervention on the target in this group?\n"
                "Provide a clear and concise answer, including the direction of the effect (if applicable).\n"
                "If there is no effect or the results are inconclusive, please state so.\n"
            ),
        )

        q[5] = AnyText(
            "evidence",
            lambda d: group_template(d)
            + (
                f"What is the evidence from the paper supporting this result?\n"
                "Provide numerical and other quotes/excerpts from the paper explaining or justifying this conclusion.\n"
            ),
        )

        q[6] = AnyText(
            "p_value",
            lambda d: group_template(d)
            + (
                f"What is the p-value of this result?\n"
                "Provide the p-value of THIS SPECIFIC result if it is stated in the paper\n"
            ),
        )

        q[7] = AnyText(
            "excerpt",
            lambda d: group_template(d)
            + (
                "Provide paragraph(s) from the paper that mention the intervention, conclusion and other data about this intervention group.\n"
                "If you believe this intervention group was not actually studied or the above data is not present in the paper, please reply with None or Error.\n"
            ),
        )

        q = list(q.values())

        result = []
        try:
            result = q[0].dfs(self.ask_llm, q[1:], {})
        except Exception as e:
            error(f"Error in processing {paper.title}", e)

        result_filtered = []

        for i in result:
            if not (
                i["excerpt"] in [None, ""]
                or "none" in i["excerpt"].lower()
                or "error" in i["excerpt"].lower()
            ):
                result_filtered.append(i)
            else:
                warning(f"Filtered out {i}")

        return result

    def build_report(self, query: str, results: list[dict]):
        prompt = (
            (
                "You are an expert scientist. \n"
                f"You are researching the impact of {query} on longevity. \n"
                "You did a thorough and rigorous analysis of papers related to the subject and extracted a number of key excerpts and data points. \n"
                "The papers are colored yellow or green according to their quality/relevance. \n"
                "Some of your findings might contain errors, be duplicates, or miss important data. \n"
                "These are the papers you have analyzed: \n"
            )
            + (
                "\n\n\n".join(
                    f"== Paper: {x['paper'].title}\n"
                    f"Abstract: {x['paper'].abstract}\n"
                    f"Year: {x['paper'].date.year if x['paper'].date is not None else '<unknown>'}\n"
                    f"Authors: {x['paper'].authors if x['paper'].authors is not None else '<unknown>'}\n"
                    f"Color: {x['paper'].color}\n"
                    f"Extracted study data:\n\n"
                    + "\n".join(
                        (
                            "=== Intervention group:\n"
                            + f"Species: {d['species']}\n"
                            + f"Intervention: {d['intervention']}\n"
                        )
                        + (
                            f"Measured outcome/target: {d['target']}\n"
                            if d.get("target")
                            else ""
                        )
                        + (f"Result: {d['result']}\n" if d.get("result") else "")
                        + (f"Evidence: {d['evidence']}\n" if d.get("evidence") else "")
                        + (f"p-value: {d['p_value']}\n" if d.get("p_value") else "")
                        + (
                            f"Supporting reference: {d['excerpt']}\n"
                            if d.get("excerpt")
                            else ""
                        )
                        for d in x["result"]
                    )
                    for x in results
                )
            )
            + "You need to make a report based on the data you have extracted. \n"
            + "Make sure to reference the papers\n"
            + "You can use Markdown\n\n\n"
            + "== Example: \n"
            "Analysis results:\n"
            "In humans:\n"
            "   — Effect on lifespan:\n"
            "        [green] 0.5 mg metformin administration increases human lifespan by <...> [John Doe et al. 2024]; 0.1 mg metformin doesn't affect human lifespan significantly [Jane Smith et al. 1999].\n"
            "        [yellow] metformin supplementation is associated with higher lifespan in diabetes patiens [Ref3].\n"
            "   — Effect on age-related diseases:\n"
            "        [green] 0.1% metformin decreases risks of developing colorectal cancer by <...> [Ref4].\n"
            "        [yellow] metformin doesn't change the risks of developing Alzheimer's [Ref5], 0.5% metformin mitigates insulin resistance [Ref6].\n"
            "   — Effect on biomarkers:\n"
            "        [green] 0.1% metformin increases SIRT1 expression [Ref7].\n"
            "        [yellow] 0.1% metformin upregulates AMPK [Ref8], and decreases arterial hypertension [Ref9]; 0.5% metformin doesn't affect blood cholesterol levels [Ref10].\n"
            "\n"
            "In mice:\n"
            "   — ...\n"
            "   — ...\n"
            "   — ...\n"
            "\n"
            "In Drosophila:\n"
            "   — ...\n"
            "   — ...\n"
            "   — ...\n"
            "\n"
            "Conclusions:\n"
            "[green] According to systematic reviews and meta-analyses published in high quality journals, metformin administration increases ..., decreases ..., and has no effect on ... [Ref, Ref, Ref...]. In mice, metformin was shown to ... [Refs...]\n"
            "[yellow] According to clinical trials and high-quality research articles, metformin ... [Refs...]. In mice, such articles demonstrated that ... [Refs...]\n"
            "[red] In addition, reviews and papers published in lower-quality journals claim that metformin ... [Ref], ... [Ref], and ... [Ref].\n"
            "To sum up, metformin increases human lifespan and ameliorates a variety of age-related biomarkers.\n"
        )

        response = self.model.chat(prompt)

        return response.choices[0].message.content

    def make_conclusion(self, desc, report):
        prompt = (
            f"""
                Analyze the results of the articles categorized by their color-coded trust levels. Provide structured summaries with key findings, numerical data, and statistical insights where available. Conduct a meta-analysis to identify patterns, trends, and discrepancies across sources.

                Task Breakdown:
                Summarize Results by Trust Category

                Meta-Analysis

                Compare trends across categories: Do lower-trust sources show biases or misinformation? Do higher-trust sources provide converging evidence?
                Highlight major points of agreement and disagreement.
                Quantify key insights where possible (e.g., "X% of sources agree on Y").
                Overall Conclusion

                Provide a synthesized final assessment based on the collective findings.
                Address reliability concerns and suggest confidence levels in the overall conclusions.
                Data Input:
                
                Description: {desc}
                
                Results: {report}
            """
        )

        response = self.model.chat(prompt)

        return response.choices[0].message.content
