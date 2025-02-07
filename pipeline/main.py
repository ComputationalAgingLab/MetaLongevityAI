import streamlit as st
from langchain_gigachat.chat_models import GigaChat
from langchain.schema import HumanMessage
from langchain.document_loaders import DirectoryLoader, TextLoader


MAX_SYMBOLS = 90_000
PATH = r'YOUR_PATH'
KEY = "YOUR_KEY"
llm = GigaChat(credentials=KEY, scope="GIGACHAT_API_CORP", model="GigaChat-Pro", verify_ssl_certs=False)


def load_documents(path):
    text_loader_kwargs = {"encoding": "utf-8"}
    loader = DirectoryLoader(path, glob="**/*.txt", loader_cls=TextLoader, loader_kwargs=text_loader_kwargs)
    documents = loader.load()
    return documents


def qa_with_context(question, context):
    system_prompt = (
        "Ты — эксперт в области медицины. Твоя задача — очно ответить на вопрос пользователя, используя только предоставленный контекст. "
        "Ответ должен быть кратким и точным. Если в контексте нет информации для ответа на конкретный вопрос напиши в ответе на него 'Недостаточно данных'. "
        "Названия или заголовки нужно писать на английском. \n\n"
        "Контекст:\n"
        f"{context}\n\n"
        "Вопросы:\n"
        f"{question}\n\n"
        "Ответ:"
    )
    return system_prompt


st.title("Анализ медицинских статей с помощью GigaChat")
documents = load_documents(PATH)
user_question = st.text_input("Введите ваш вопрос:")

if st.button("Создать таблицу с ответами"):
    questions = [
        "Как называется статья?",
        "Организация, проводившая исследование?",
        "Какой тип у этой статьи?",
        "В каком году опубликована эта статья?",
        "Какой эффект имел объект исследования на продолжительность жизни? (положительный, нейтральный, отрицательный, неопределенный)",
        "Общее число испытуемых в исследовании?",
        "Число испытуемых в контрольной группе?",
        "Число мужчин в контрольной и тестовой группах?",
        "Число женщин в контрольной и тестовой группах?",
        "Возрастной диапазон испытуемых?",
        "Длительность исследования?",
        "Страна, в которой проводилось исследование?"
    ]

    answers = []
    with st.spinner("Обработка документов..."):
        for con in documents:
            for que in questions:
                try:
                    prompt = qa_with_context(que, con.page_content[:min(MAX_SYMBOLS, len(con.page_content))])
                    answer = llm([HumanMessage(content=prompt)]).content
                    answers.append(answer)
                    st.write(answer)
                except Exception as e:
                    st.error(f"Ошибка при обработке документа: {e}")

    with st.spinner("Создание таблицы..."):
        try:
            table_prompt = (
                "На основе следующих вопросов и ответов создай таблицу. "
                "Ответы и вопросы должны быть сокращены до пары слов, но отражать всю суть. "
                "Так же нужно переписать вопросы в форме более подходящей для таблиц. \n"
                f"Вопросы: {'; '.join(questions)}. "
                f"Ответы: {'; '.join(answers)}."
            )
            table_answer = llm([HumanMessage(content=table_prompt)]).content

            st.write(table_answer)
        except Exception as e:
            st.error(f"Ошибка при создании таблицы: {e}")


if st.button("Сделать метаанализ"):
    questions = [
        "Как называется статья?",
        "Организация, проводившая исследование?",
        "Какой тип у этой статьи?",
        "В каком году опубликована эта статья?",
        "Какой эффект имел объект исследования на продолжительность жизни? (положительный, нейтральный, отрицательный, неопределенный)",
        "Общее число испытуемых в исследовании?",
        "Число испытуемых в контрольной группе?",
        "Число мужчин в контрольной и тестовой группах?",
        "Число женщин в контрольной и тестовой группах?",
        "Возрастной диапазон испытуемых?",
        "Длительность исследования?",
        "Страна, в которой проводилось исследование?"
    ]

    answers = []
    with st.spinner("Обработка документов..."):
        for con in documents:
            for que in questions:
                try:
                    prompt = qa_with_context(que, con.page_content[:min(MAX_SYMBOLS, len(con.page_content))])
                    answer = llm([HumanMessage(content=prompt)]).content
                    answers.append(answer)
                except Exception as e:
                    st.error(f"Ошибка при обработке документа: {e}")

    with st.spinner("Создание метаанализа..."):
        try:
            meta_prompt = (
                "На основе следующих вопросов и ответов сделай подробный метаанализ статей. "
                "Сделай выводы о вреде и пользе. Используй статистику и удели внимание результатам воздействия. "
                "Очень важно отвечать по теме основного вопроса пользователя: "
                f"{user_question}"
                f"Вопросы: {'; '.join(questions)}. "
                f"Ответы: {'; '.join(answers)}."
            )
            meta_analysis = llm([HumanMessage(content=meta_prompt)]).content

            st.success("Метаанализ:")
            st.write(meta_analysis)
        except Exception as e:
            st.error(f"Ошибка при создании метаанализа: {e}")
