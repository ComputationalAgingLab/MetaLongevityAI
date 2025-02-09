import pandas as pd
from fuzzywuzzy import process


def get_q_of_article(title: str) -> str:
    df = pd.read_excel("data.xlsx")
    print(df.head())  # Вывести первые 5 строк

    df["Квартиль"] = pd.qcut(df["Cites / Doc. (2years)"], q=4, labels=["Q1", "Q2", "Q3", "Q4"])

    matches = process.extract(title, df['Title'], limit=1)
    result = df[df["Title"] == matches[0][0]]
    print(matches[0][0])
    return result.loc[0, 'Квартиль']


print(get_q_of_article('Nature Reviews Molecular  Biology'))
