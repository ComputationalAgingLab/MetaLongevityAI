import json
import time

import requests
from PyPDF2 import PdfReader
from bs4 import BeautifulSoup

from gigachat import GigaChat


def extract_abstract(article_text: str) -> str:
    model = GigaChat(
        credentials="YTUwN2I4ZTQtM2Q1ZC00NzY3LTkzYWQtOTI0MTlkMGQ3NGI4OmI3YzFkMTM4LTEyNTgtNDNkYS1hZTNkLTNiYmY4ZTc5NjFiMA==",
        scope="GIGACHAT_API_CORP", model="GigaChat-Pro",
        verify_ssl_certs=False, )
    response = model.chat("Извлеки раздел 'Abstract' (аннотацию) из следующей статьи. "
                          "Текст ответа должен состоять только из текста запроса. "
                          "Ответ должен содержать только текст раздела 'Abstract', без дополнительных комментариев:\n\n"
                          f"{article_text}")
    return response.choices[0].message.content


def from_title_to_filename(title):
    return "".join([x for x in title if x.isalpha() or x.isdigit()]) + '.pdf'


def pdf_to_json(filename):
    reader = PdfReader(filename)
    number_of_pages = len(reader.pages)
    text = ''
    for i in range(number_of_pages):
        page = reader.pages[i]
        text += page.extract_text()

    abstract = extract_abstract(text)

    return abstract, text


def extract_download_link(html):
    soup = BeautifulSoup(html, 'html.parser')
    button_tag = soup.find('button', onclick=True)
    if button_tag:
        link = button_tag['onclick'].split("'")[1]
        link = 'https:' + link

        return link
    return None


def get_doi_year(title):
    """Ищет DOI статьи по её названию через CrossRef API"""
    url = f"https://api.crossref.org/works?query.title={title}&rows=1"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()

        if "message" in data and "items" in data["message"] and data["message"]["items"]:
            try:
                year = str(data["message"]["items"][0].get('license')[0].get('start').get('date-parts')[0][0])
            except Exception as e:
                year = None
            return data["message"]["items"][0].get("DOI"), year
    return None, None


def download_paper(doi, output_filename):
    time.sleep(1)
    sci_hub_url = "https://sci-hub.se"  # Рабочее зеркало может измениться
    search_url = f"{sci_hub_url}/{doi}"

    headers = {"User-Agent": "Mozilla/5.0"}  # Имитация браузера
    session = requests.Session()

    # Получаем HTML страницы
    response = session.get(search_url, headers=headers)
    if response.status_code != 200:
        print("Ошибка доступа к Sci-Hub!")
        print(sci_hub_url)
        return

    # Извлекаем ссылку на PDF
    download_url = extract_download_link(response.text)
    # print('not formated', download_url)

    if download_url is None:
        print('Cтатья отсутствует в sci hub')
        return False
    else:
        download_url = download_url.replace('https:/downloads', 'https://sci-hub.ru/downloads')
        download_url = download_url.replace('https:/tree', 'https://sci-hub.ru/tree')
        print('formated', download_url)

        pdf_response = session.get(download_url, headers=headers, stream=True)
        if pdf_response.status_code == 200:
            with open(output_filename, "wb") as f:
                for chunk in pdf_response.iter_content(1024):
                    f.write(chunk)
            print(f"Файл сохранен как {output_filename}")
        else:
            print("Ошибка при загрузке PDF!")
        return True


# Функция для поиска статей
def search_google_scholar(query, num_results=10):
    pages = num_results // 10

    # Находим все результаты
    results = []

    for i in range(pages):
        time.sleep(5)
        search_url = f'https://scholar.google.com/scholar?start={i * 10}&q={query}'
        headers = {
            'User-Agent': 'AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36'
        }
        # Выполнение запроса

        response = requests.get(search_url, headers=headers)

        soup = BeautifulSoup(response.text, "html.parser")

        for item in soup.find_all("h3"):  # Заголовки статей обычно находятся в <h3>
            title = item.get_text()
            title = ''.join([x for x in title if x.isalpha() or x == ' '])

            link = item.a['href'] if item.a else None
            results.append({
                'title': title,
                'link': link
            })

    return results


def main(obj=None):
    # searching at google scholar
    query = ('smoking longevity interventions')
    articles = search_google_scholar(query, num_results=20)

    cnt = 0
    papers = []
    for a in articles:
        title = a['title']

        doi, year = get_doi_year(title)
        # stroing papers locally
        filename = from_title_to_filename(title)
        if download_paper(doi, output_filename=f"{filename}.pdf"):
            cnt += 1
            papers += [(title, doi, year, filename)]

    print(f"{cnt}/{len(articles)} papers are saved")
    # forming json file
    general = []
    for p in papers:
        title, doi, year, filename = p
        abstract, fulltext = pdf_to_json(f"{title}.pdf")
        dct = {
            'title': title,
            'year': year,
            'doi': doi,
            'abstract': abstract,
            'fulltext': fulltext,
            'filename': filename
        }
        general += [dct]
    with open("data.json", "w") as file:
        json.dump(general, file)


if __name__ == "__main__":
    main()
