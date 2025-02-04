from pathlib import Path
import pubget


def main():
    # Аналог запроса 'coffee longevity effect' в
    # https://www.ncbi.nlm.nih.gov/pmc/?term=coffee+longevity+effect
    query = 'coffee AND longevity AND effect'
    pubget_files_dir = Path("C:/Users/Егор Швецов/PycharmProjects/metalongai/data").expanduser()

    print('download_query_results')
    created_dir, err = pubget.download_query_results(query, pubget_files_dir, n_docs=3, retmax=3)
    print(created_dir)

    print('extract_articles')
    articlesets_dir = str(created_dir)
    articles_dir = f'{pubget_files_dir}/articles'
    created_dir, err = pubget.extract_articles(articlesets_dir, articles_dir, n_jobs=-1)
    print(created_dir)

    print('extract_data_to_csv')
    extracted_data_dir = f'{pubget_files_dir}/extracted_data'
    created_dir, err = pubget.extract_data_to_csv(articles_dir, extracted_data_dir, n_jobs=-1)

    print(created_dir)

    print('extract_vocabulary_to_csv')
    vocabulary_dir = f'{pubget_files_dir}/vocabulary'
    created_dir, err = pubget.extract_vocabulary_to_csv(extracted_data_dir, vocabulary_dir)
    print(created_dir)


if __name__ == "__main__":
    main()
