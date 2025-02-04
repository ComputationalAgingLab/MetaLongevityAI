#!/usr/bin/env python3

# %%writefile tojson.py
import os
import json


# This function converts brat format to json.
def ann2json(txt_link, ann_link, output_folder, text_id="name"):
    """ convert txt and ann files to JSON;
        text_id: default = "name" - take name from file basename,
        or any integer for manual name assigning"""

    if text_id == "name":
        text_id = os.path.basename(txt_link).split(".txt")[0]

    with open(txt_link, "r") as f:
        text = f.read()
    with open(ann_link, "r") as f:
        raw_lines = f.readlines()

    jout = {"text": {}, "id": {}, "entities": {}, "relations": {}}
    jout["text"] = text
    jout["id"] = text_id  # <- need text IDx

    if raw_lines:
        lines_T = [l.strip().split("\t") for l in raw_lines if l[0] == "T"]
        lines_R = [l.strip().split("\t") for l in raw_lines if l[0] == "R"]

        if lines_T:
            ids, lbl_spans, entities = list(zip(*lines_T))
            types, starts, ends = list(zip(*[f.split() for f in lbl_spans]))
            for i, t in enumerate(ids):
                jout["entities"].update(
                    {
                        t: {
                            "end": str(int(starts[i]) + len(entities[i])),
                            "text": entities[i],
                            "entity": types[i],
                            "start": starts[i],
                        }
                    }
                )
        if lines_R:
            rs, relations = list(zip(*lines_R))
            rels = [
                [rel[0], rel[1].split(":")[1], rel[2].split(":")[1]]
                for rel in [r.split() for r in relations]
            ]
            for i, r in enumerate(rs):
                jout["relations"].update(
                    {r: {"Arg1": rels[i][1], "Arg2": rels[i][2], "type": rels[i][0]}}
                )

    with open("{0}/{1}.json".format(output_folder, text_id), "w") as f:
        json.dump(jout, f)

    return None
    # return jout


from pdfminer.high_level import extract_pages
from pdfminer.layout import LTTextContainer


def pdf2txt(path, path2):
    with open(path2, "a+") as fout:
        for page_layout in extract_pages(path):
            for element in page_layout:
                if isinstance(element, LTTextContainer):
                    fout.write(element.get_text())


# %%writefile query_parser.py
import json
import pubmed_parser as pp
from Bio import Entrez


def do_query(
        name, search_query, papers_count, folder_path, email="firepaladiner@gmail.com"
):
    all_names = list(
        map(
            lambda x: x if os.path.isdir(os.path.join(folder_path, x)) else "",
            os.listdir(folder_path),
        )
    )
    if not name:
        all_untitled = filter(lambda x: x.startswith("Untitled"), all_names)
        max_num = 0
        for element in all_untitled:
            element = element.replace("Untitled", "")
            if "-" in element:
                element = element.replace("-", "")
            if not element:
                element = 1
            if int(element) > max_num:
                max_num = int(element)
        name = "Untitled-%d" % (max_num + 1)
    elif name in all_names:
        all_clones = filter(lambda x: x.startswith(name), all_names)
        max_num = 0
        for element in all_clones:
            element = element.replace(name, "")
            if "-" in element:
                element = element.replace("-", "")
            if not element:
                element = 1
            if int(element) > max_num:
                max_num = int(element)
        name = name + "-%d" % (max_num + 1)
    with open(os.path.join(folder_path, name + ".query"), "w+") as query_file:
        query_file.write(search_query)
    folder_path = os.path.join(folder_path, name)
    Entrez.email = email

    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

    def make_query(search_query, papers_count=10):
        """Return list with top PMC IDs by query"""
        # getting search results for the query
        handle = Entrez.esearch(
            db="pmc",
            term=search_query,
            retmax=papers_count,
            usehistory="y",
            prefix="xlink",
        )
        search_results = Entrez.read(handle)
        handle.close()
        return search_results["IdList"]  # list with PMC IDs

    def load_papers(idlist, folder_path):
        """Download papers from Id list"""
        for pmc in idlist:
            ### take a paper xml
            handle = Entrez.efetch(db="pmc", id=pmc, rettype="full", retmode="xml")
            ###
            xml_path = os.path.join(folder_path, pmc + ".xml")
            # save xml
            with open(xml_path, "w") as f:
                f.write(handle.read().decode("utf-8"))


            dict_meta = pp.parse_pubmed_xml(xml_path)
            dict_par = pp.parse_pubmed_paragraph(xml_path)

            meta_path = os.path.join(folder_path, pmc + ".meta.json")
            text_path = os.path.join(folder_path, pmc + ".txt")

            ### save a paper meta info
            with open(meta_path, "w") as fp:
                json.dump(dict_meta, fp)

            ### save a paper text
            full_text = "\n".join([p["text"] for p in dict_par])
            with open(text_path, "w",encoding="utf-8") as fp:
                fp.write(full_text)

    idlist = make_query(search_query=search_query, papers_count=papers_count)
    print(
        "Query finished: %d papers from %d will be downloaded."
        % (len(idlist), papers_count)
    )
    load_papers(idlist=idlist, folder_path=folder_path)
    print("FINISHED")


import pandas as pd
import numpy as np
from glob import glob


def get_snippet(text, coo1, coo2, expand=100):
    """Return text snippet around two entities"""
    tlen = len(text)
    cmin = min(coo1[0], coo2[0])
    cmax = max(coo1[1], coo2[1])
    cmin = cmin - expand if cmin - expand >= 0 else 0
    cmax = cmax + expand if cmax + expand < tlen else tlen
    return text[cmin:cmax]


def create_dataframe_from_jsons(folder_path):
    all_filenames = os.listdir(folder_path)
    paths = []
    for filename in all_filenames:
        name, extension = os.path.splitext(filename)
        if extension == ".json" and ".meta.json" not in filename:
            paths.append(os.path.join(folder_path, filename))

    if len(paths) == 0:
        print("No files found.")
        return None
    # load all jsons
    docs = []
    for p in paths:
        with open(p) as f:
            docs.append(json.load(f))

    # create template for dataFrame
    df = pd.DataFrame(
        {
            "Document": [],
            "Entity A": [],
            "Relation type": [],
            "Entity B": [],
            "Entity A type": [],
            "Entity B type": [],
            "Text snippet": [],
        }
    )

    # fill dataFrame
    for p, j in zip(paths, docs):
        if len(j["relations"]) != 0:
            for R, obj in j["relations"].items():
                # extract relation and corresponding entities
                r_arg1, r_arg2 = obj["Arg1"], obj["Arg2"]
                r_type = obj["type"]
                e_arg1 = j["entities"][r_arg1]["text"].lower()
                e_arg2 = j["entities"][r_arg2]["text"].lower()
                e_type1 = j["entities"][r_arg1]["entity"]
                e_type2 = j["entities"][r_arg2]["entity"]
                e_coo1 = (
                    int(j["entities"][r_arg1]["start"]),
                    int(j["entities"][r_arg1]["end"]),
                )
                e_coo2 = (
                    int(j["entities"][r_arg2]["start"]),
                    int(j["entities"][r_arg2]["end"]),
                )
                text_snippet = (
                        "..." + get_snippet(j["text"], e_coo1, e_coo2, expand=50) + "..."
                )
                # construct DataFrame
                df = df.append(
                    {
                        "Document": p,
                        "Entity A": e_arg1,
                        "Relation type": r_type,
                        "Entity B": e_arg2,
                        "Entity A type": e_type1,
                        "Entity B type": e_type2,
                        "Text snippet": text_snippet,
                    },
                    ignore_index=True,
                )
    return df


def make_relation_dataframe(df):
    cols = ["Entity A", "Relation type", "Entity B", "Text snippet"]
    rdf = df[cols].copy()
    rdf["Count"] = rdf["Entity A"] + rdf["Relation type"] + rdf["Entity B"]
    mapper = rdf.groupby("Count", as_index=False)
    rdf = mapper.agg(
        {
            "Entity A": "first",
            "Relation type": "first",
            "Entity B": "first",
            "Count": "size",
            "Text snippet": "first",
        }
    )
    rdf["In documents"] = [
        list(set(df.loc[mapper.indices[key]]["Document"]))
        for key in mapper.indices.keys()
    ]
    rdf = rdf.sort_values("Count", ascending=False)
    return rdf


def find_keys_to_remove(doc, EA, R, EB):
    """Check if the doc contains the relation R"""
    de = doc["entities"]
    dr = doc["relations"]
    _keys = [k for k in dr.keys() if dr[k]["type"] == R]
    keys2remove = []
    for k in _keys:
        e1 = de[dr[k]["Arg1"]]["text"].lower()
        e2 = de[dr[k]["Arg2"]]["text"].lower()
        if (e1 == EA) and (e2 == EB):
            keys2remove.append(k)
    return keys2remove


def delete_relation(rdf, idx):
    # rdf = pd.DataFrame(rdf)
    paths = rdf.loc[idx]["In documents"]
    EA, R, EB = rdf.loc[idx][["Entity A", "Relation type", "Entity B"]]

    # load corresponding docs
    docs = []
    for p in paths:
        with open(p) as f:
            docs.append(json.load(f))

    # delete keys from corresponding docs
    for doc in docs:
        rkeys = find_keys_to_remove(doc, EA, R, EB)
        for r in rkeys:
            doc["relations"].pop(r)

    # rewrite corresponding docs
    for doc, p in zip(docs, paths):
        with open(p, "w") as f:
            json.dump(doc, f)
    rdf = rdf.drop(idx)
    return rdf


def save_relation(rdf, idx, knowledge_base_path):
    # rdf = pd.DataFrame(rdf)
    new = str(rdf.loc[idx].to_dict())
    if not os.path.exists(knowledge_base_path):
        with open(knowledge_base_path, "w+") as tmp:
            tmp.write("")

    with open(knowledge_base_path, "r") as f:
        if new not in f.read():
            with open(knowledge_base_path, "a") as f:
                f.write(new + "\n")
