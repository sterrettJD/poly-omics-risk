import requests
import pandas as pd
from io import StringIO

def get_compound_info(compound):
    res = requests.get(f"https://rest.kegg.jp/find/compound/{compound}")
    return res.content.decode('UTF-8')


def parse_res_to_df(info):
    res_df = pd.read_csv(StringIO(info), sep="\t", header=None)
    res_df.columns = ["CO", "Names"]
    return res_df


def main():
    info = get_compound_info("phenyllactate")
    info = parse_res_to_df(info)

    print(info)



if __name__=="__main__":
    main()