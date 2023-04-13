import requests
import pandas as pd
def get_compound_info(compound):
    res = requests.get(f"https://rest.kegg.jp/find/compound/{compound}")
    return res.content
def main():
    info = get_compound_info("phenyllactate")
    print(info)



if __name__=="__main__":
    main()