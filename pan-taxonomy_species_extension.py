import os
import subprocess
import json
from Bio import Phylo
from collections import defaultdict


def read_json(json_dir):
    file = open(json_dir, 'r', encoding='UTF-8')
    dic = json.loads(file.read())
    file.close()
    return dic


def pan_homology_extension():
    print('start')
    pan_species_list = ["_".join(i.split("_")[0:2]).capitalize() for i in
                        os.listdir(work_dir + "homology_json/pan_homology")]
    pan_homology = read_json(work_dir + "homology_json/pan_homology" + species)

    # make total species list
    for gene in pan_homology[species]:
        mid = dict(zip(total_species_list, [[] for i in total_species_list]))
        mid.update(pan_homology[species][gene])
        for i in mid:
            mid[i] = set(mid[i])
        pan_homology[species][gene] = mid

    # the pan-taxonomy entension species' own domain homology info
    for division in os.listdir(work_dir + "homology_json"):
        print(division)
        if division == 'pan_homology' or division == species2division[species]:
            continue
        division_homology = read_json(work_dir + "division_json/" + division)

        print('have read the division json')
        for gene in pan_homology[species]:
            for species_in_division_homology in set(division2species[division]) & set(pan_species_list):
                for gene_in_division_homology in pan_homology[species][gene][species_in_division_homology]:
                    try:
                        for species_to_ex in set(
                                division_homology[species_in_division_homology][gene_in_division_homology]) & set(
                            division2species[division]):
                            pan_homology[species][gene][species_to_ex] = pan_homology[species][gene][species_to_ex] | \
                                                                         set(division_homology[
                                                                                 species_in_division_homology][
                                                                                 gene_in_division_homology][
                                                                                 species_to_ex])
                    except:
                        pass

    # add the out_group species' own domain homology info
    for species_in_out_group in out_group:
        for division_in_out_group in out_group[species_in_out_group]:
            if species2division[species] != division_in_out_group:
                continue
            print(species_in_out_group, division_in_out_group)
            out_group_homology = {species_in_out_group: read_json(
                work_dir + "homology_json/" + division_in_out_group + '/' + species2species[
                    species_in_out_group])[species2species[species_in_out_group]]}
            for gene in pan_homology[species]:
                for homology_gene in pan_homology[species][gene][species_in_out_group]:
                    try:
                        for species_add in out_group_homology[species_in_out_group][homology_gene]:
                            pan_homology[species][gene][species_add] = pan_homology[species][gene][species_add] | \
                                                                       set(out_group_homology[species_in_out_group][
                                                                               homology_gene][species_add])
                    except:
                        pass

    # trans the set to list
    for gene in pan_homology[species]:
        for species_homology in pan_homology[species][gene]:
            pan_homology[species][gene][species_homology] = list(pan_homology[species][gene][species_homology])
    function_dir = work_dir + "pan_extension/" + species
    file_content = """%s""" % json.dumps(pan_homology)
    outfile = open(function_dir, "w")
    outfile.write(file_content)
    outfile.close()


work_dir = 'GenOrigin/tree/master/'

nwk_dir = work_dir + "nwk/total.nwk"

out_group = {"Ustilago_maydis": ["bacteria"], "Saccharomyces_cerevisiae": ["plants", "vertebrates"],
             "Homo_sapiens": ["plants"], "Drosophila_melanogaster": ["plants", "vertebrates"],
             "Caenorhabditis_elegans": ["plants", "vertebrates"]}

species2division = {}
division2species = defaultdict(lambda: [])
species2species = {}
for division in os.listdir(work_dir + "homology_species/"):
    if division == "pan_homology":
        continue
    for species in os.listdir(work_dir + 'homology_species/' + division):
        if "_".join(species.split("_")[0:2]).capitalize() in out_group:
            if division in out_group["_".join(species.split("_")[0:2]).capitalize()]:
                continue
        division2species[division].append("_".join(species.split("_")[0:2]).capitalize())
        species2division[species] = division
        species2species['_'.join(species.split('_')[0:2]).capitalize()] = species


total_species_list = list(set(["_".join(i.split("_")[0:2]).capitalize() for j in
                               os.listdir(work_dir + "homology_species/") for i
                               in os.listdir(work_dir + "homology_species/" + j)]))

for species in os.listdir(work_dir + 'homology_species/pan_homology'):
    pan_homology_extension()
