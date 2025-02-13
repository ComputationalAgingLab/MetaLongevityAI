from langchain_gigachat.chat_models import GigaChat
from langchain.schema import HumanMessage
import giga_token


GIGA_KEY = giga_token.GIGA_KEY
LOG = True

llm = GigaChat( credentials=GIGA_KEY,
                scope="GIGACHAT_API_CORP",
                model="GigaChat-Pro",
                verify_ssl_certs=False,)


def log_info(msg):
    if LOG:
        print("[INFO] " + msg)


def log_question(question, answer=''):
    if not LOG:
        return
    print("[QUESTION]", question)
    if answer != '':
        print("[ANSWER]", answer)


def binary_question(question):
    """Function for asking yes/no question"""

    question += " Respond with ONLY yes/no"
    response = llm.invoke([HumanMessage(content=question)]).content
    log_question(question, response)
    return response.lower().startswith('yes')


# Add front-end feedback to bad outcomes
def fool_check(query):
    """Check if [query] is an existing word and if it is related to medicine"""

    q_typos = f"Are there any typing mistakes in any of these words: {str(query.split())[1:-1]}?"
    ans_typos = llm.invoke([HumanMessage(content=q_typos)]).content
    log_question(q_typos, ans_typos)
    typos = ans_typos.lower().startswith('yes')
    if typos:
        # The word does not exist or typo, return template answer
        log_info(f"Typos in '{query}'")
        return False
    health_related = binary_question(f"Can '{query}' affect health in a medical setting?")
    if not health_related:
        # Is not related to medicine, return template answer
        log_info(f"'{query}' is not related to medicine")
        return False
    log_info("Fool check passed")
    return True
    

def get_description(query):
    """Get description of intervention [query]"""

    chem_question = f"Is {query} a chemical substance?"
    is_chem = binary_question(chem_question)
    
    sentence_limit = (3, 5)
    question = ("Imagine you are a medical scientist, describe intervention "
                f"{query.lower()} using only {'-'.join(map(str, sentence_limit))} sentences. ")
    question += ("Give a short chemical description and write where it can be found naturally. "
                 "Dont say anything about it influence on diseases. "
                 "Also write in wich molecular ways it can act on an organism. ")\
                if is_chem else \
                ("Give a short description of types used in medical studies. ")
    question += "Be concise and coherent. "
    response = llm.invoke([HumanMessage(content=question)]).content
    log_question(question, response)
    return response


if __name__ == '__main__':
    query = 'vitamin K2+D3'
    if fool_check(query):
        desc = get_description(query)
