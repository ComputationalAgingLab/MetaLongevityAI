from dataclasses import dataclass
from datetime import datetime
from enum import Enum
from typing import Optional

@dataclass
class Paper:
    title: str
    abstract: str
    journal: str
    date: datetime
    doi: Optional[str]
    pdf: Optional[bytes]
    full_text: Optional[str]

class PaperType(Enum):
    SYSTEMATIC_REVIEW = 1
    RANDOMIZED_CONTROLLED_TRIAL = 2
    OBSERVATIONAL_STUDY = 3
    META_ANALYSIS = 4

class OutcomeType(Enum):
    BINARY = 1
    CONTINUOUS = 2

class ResultType(Enum):
    DECREASED = 1
    INCONCLUSIVE = 2
    INCREASED = 3

@dataclass
class ResultL1:
    intervention: str
    comparator: str
    outcome_type: OutcomeType
    result: ResultType

@dataclass
class ReportL1:
    paper: Paper
    full_text: str # MD
    results: list[ResultL1]