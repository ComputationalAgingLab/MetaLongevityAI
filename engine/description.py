from langchain_gigachat.chat_models import GigaChat
from langchain.schema import HumanMessage
from logging import info

class PaperDescriptionSource:
    def __init__(self, llm):
        self.llm = llm

    def ask_llm(self, query: str) -> str:
        info("Asking LLM: ...\n" + "\n".join(query.split("\n")[-20:]))
        return self.llm.chat(query).choices[0].message.content

    def binary_question(self, question: str):
        """Function for asking yes/no question"""

        question += " Respond with ONLY yes/no"
        response = self.ask_llm(question)
        info(f"Q: {question} A: {response}")
        return response.lower().startswith("yes")


    # Add front-end feedback to bad outcomes
    def fool_check(self, query):
        """Check if [query] is an existing word and if it is related to medicine"""

        q_typos = f"Are there any typing mistakes in any of these words: {str(query.split())[1:-1]}?"
        # ans_typos = self.llm.invoke([HumanMessage(content=q_typos)]).content
        ans_typos = self.ask_llm(q_typos)
        info(f"Q: {q_typos} A: {ans_typos}")
        typos = ans_typos.lower().startswith("yes")
        if typos:
            # The word does not exist or typo, return template answer
            info(f"Typos in '{query}'")
            return False
        health_related = self.binary_question(
            f"Can '{query}' affect health in a medical setting?"
        )
        if not health_related:
            # Is not related to medicine, return template answer
            info(f"'{query}' is not related to medicine")
            return False
        return True


    def get_description(self, query):
        """Get description of intervention [query]"""

        chem_question = f"Is {query} a chemical substance?"
        is_chem = self.binary_question(chem_question)

        sentence_limit = (3, 5)
        question = (
            "You are an expert medical scientist. Describe  "
            f"{query.lower()} using only {'-'.join(map(str, sentence_limit))} sentences. "
        )
        question += (
            (
                "Give a short chemical description and explain where it can be found naturally. "
                "Dont say anything about its influence on diseases. "
                "State in which ways it can act on a human from the perspective of biochemistry. "
            )
            if is_chem
            else ("Give a short description of types used in medical studies. ")
        )
        question += "Be concise and coherent. "
        # response = self.llm.invoke([HumanMessage(content=question)]).content
        response = self.ask_llm(question)
        # log_question(question, response)
        return response