import requests
import pandas as pd
from io import StringIO


def read_input_file(fp):
    with open(fp, "r") as file:
        file_contents = file.readline().split("\t")
    return file_contents

def get_compound_info(compound):
    res = requests.get(f"https://rest.kegg.jp/find/compound/{compound}")
    return res.content.decode('UTF-8')


def parse_res_to_df(info):
    res_df = pd.read_csv(StringIO(info), sep="\t", header=None)
    res_df.columns = ["CO", "Names"]
    return res_df


def main():
    compounds = read_input_file("selected_metabolites.txt")
    print(compounds)

    info = get_compound_info("phenyllactate")

    if info is None:
        return

    info = parse_res_to_df(info)

    print(info)



if __name__=="__main__":
    main()