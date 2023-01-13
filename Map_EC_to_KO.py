from Bio.KEGG import REST
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas

def get_mapper():
    # Pull the linker from KEGG's REST API
    result = REST.kegg_link("ko", "enzyme").read()

    #convert to dictionary for mapping
    splitted = result.split()
    mapper = {splitted[i]: splitted[i+1] for i in range(0, len(splitted), 2)}
    return mapper

def map_ecs_to_ko(ecs, mapper):
    return [mapper[ec] if ec in mapper.keys() else None for ec in ecs]

if __name__=="__main__":
    mapper = get_mapper()
    with open("selected_taxa_ECs.txt") as f:
        ec_list = f.read().split()
    ec_list = [f"ec:{ec}" for ec in ec_list]
    mapped = map_ecs_to_ko(ec_list, mapper)

    unmapped = sum([m is None for m in mapped])/len(mapped) * 100
    print(f"Percent unmapped ECs: {round(unmapped,2)}")

    mapped_only = [m for m in mapped if m is not None]
    print(mapped_only)

    with open(r'selected_taxa_KOs.txt', 'w') as fp:
        fp.write('\t'.join(mapped_only))