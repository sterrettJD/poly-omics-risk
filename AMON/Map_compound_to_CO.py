import requests
import pandas as pd
from io import StringIO


def read_input_file(fp):
    with open(fp, "r") as file:
        file_contents = [x.strip() for x in file.readline().split("\t")]
    return file_contents

def get_compound_info(compound):
    res = requests.get(f"https://rest.kegg.jp/find/compound/{compound}")
    return res.content.decode('UTF-8')


def parse_res_to_df(info):
    try:
        res_df = pd.read_csv(StringIO(info), sep="\t", header=None)
        res_df.columns = ["CO", "Names"]
    except pd.errors.EmptyDataError:
        return None
    return res_df

def find_matching_compound(compound, info):
    compound_names = [compound.lower(), compound.lower().replace("oate", "oic acid")]

    info["Names"] = info["Names"].apply(lambda x: [name.strip() for name in x.split(";")])
    info["Perfect Match"] = info["Names"].apply(lambda x: any([compound_name in [name.lower()
                                                                                 for name in x]
                                                               for compound_name in compound_names]))

    match_row = info.loc[info["Perfect Match"] == True]
    if len(match_row) == 0:
        info["Imperfect Match"] = info["Names"].apply(lambda x: any([compound_name in name.lower()
                                                                     for name in x for compound_name in compound_names]))
        match_row = info.loc[info["Imperfect Match"] == True]

    return match_row


def main():
    compounds = read_input_file("selected_metabolites.txt")
    print(compounds)

    COs = []
    imperfect_COs = []
    unknown_COs = []

    for compound in compounds:
        info = get_compound_info(compound)

        info = parse_res_to_df(info)

        # If not getting any results, try changing the naming convention from -oate to -oic acid
        if info is None:
            info = get_compound_info(compound.replace("oate", "oic"))
            info = parse_res_to_df(info)

        if info is not None:

            match_row = find_matching_compound(compound, info)

            if any(match_row["Perfect Match"] == True):
                COs = COs + [x.split(":")[1] for x in match_row["CO"].to_list()]

            elif any(match_row["Imperfect Match"] == True):
                imperfect_COs.append((compound,
                                      match_row["CO"].to_list(),
                                      match_row["Names"].to_list()))

            else:
                unknown_COs.append(compound)

        else:
            unknown_COs.append(compound)
            print(f"No info found for {compound}")

    print("The accepted COs are as follows: ", COs)
    with open("selected_COs_draft.txt", "w") as f:
        for CO in COs:
            f.write(CO)
            f.write("\t")

    print("Accepted COs have been written to selected_COs_draft.txt")

    print("The following compounds had imperfect matches: ", imperfect_COs)
    print("Please select the COs you would like to use, and manually add them to the file.")

    print("Information could not be found for the following compounds. Please manually curate COs", unknown_COs)

if __name__=="__main__":
    main()