import os
import json
from Bio import Phylo
from collections import defaultdict


# species_to_pan_homology
# species:gene:map_gene:map_species
# make a list for pan_homology_species
# the output= {species: {gene id: {identity query: {ortholog gene id: [ortholog species]}}}}
def make_species_to_pan_homology_json():
    judge_gene = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: []))))

    with open('/home/lchen/genorigin/homology_species/' + division + "/" + species) as f:
        next(f)
        for line in f.readlines():
            line = line.strip()
            line = line.split("\t")
            if line[2] != line[7]:
                sim_species_name = "_".join(line[2].split("_")[0:2]).capitalize()
                sim_homology_name = "_".join(line[7].split("_")[0:2]).capitalize()

                if line[2] == species and sim_homology_name in pan_homology_species_list:
                    judge_gene[sim_species_name][line[0]][sim_homology_name][float(line[3])].append(line[5])
                elif line[7] == species and sim_species_name in pan_homology_species_list:
                    judge_gene[sim_homology_name][line[5]][sim_species_name][float(line[8])].append(line[0])

    species_to_pan_output = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: []))))

    for gene in judge_gene["_".join(species.split("_")[0:2]).capitalize()]:

        select_pool = defaultdict(lambda: [])
        for distance in sorted(species_to_pan_homology_distance["_".join(species.split("_")[0:2]).capitalize()]):
            species_list_choose_by_distance = species_to_pan_homology_distance["_".join(species.split("_")[0:2]).capitalize()][distance]

            for species_choose_by_distance in set(species_list_choose_by_distance) & set(judge_gene["_".join(species.split("_")[0:2]).capitalize()][gene]):
                for score in judge_gene["_".join(species.split("_")[0:2]).capitalize()][gene][species_choose_by_distance]:
                    select_pool[score].append([judge_gene["_".join(species.split("_")[0:2]).capitalize()][gene][
                                     species_choose_by_distance][score], species_choose_by_distance])

            if select_pool:

                select_species = select_pool[sorted(select_pool, reverse=True)[0]][0][1]
                select_gene = select_pool[sorted(select_pool, reverse=True)[0]][0][0][0]
                select_score = str(sorted(select_pool, reverse=True)[0])
                species_to_pan_output["_".join(species.split("_")[0:2]).capitalize()][gene][select_score][select_gene] = [select_species]
                break

    homologyB_dir = work_dir + "species_to_pan_species/" + species
    file_content = """%s""" % json.dumps(species_to_pan_output)
    outfile = open(homologyB_dir, "w")
    outfile.write(file_content)
    outfile.close()


out_group = {"Ustilago_maydis": ["bacteria"], "Saccharomyces_cerevisiae": ["plants", "vertebrates"],
             "Homo_sapiens": ["plants"], "Drosophila_melanogaster": ["plants", "vertebrates"],
             "Caenorhabditis_elegans": ["plants", "vertebrates"]}

work_dir = 'GenOrigin/tree/master/'

with open(work_dir + 'genorigin_species_list/pan_homology') as f:
    pan_homology_species_list = [line.strip() for line in f.readlines()]

species_list = {}

for division in os.listdir(work_dir+'homology_split_by_species/'):
    if division != "pan_homology":
        for species_to_make_species_list in os.listdir(work_dir+'homology_split_by_species/' + division):
            if "_".join(species_to_make_species_list.split("_")[0:2]).capitalize() in out_group:
                if division in out_group["_".join(species_to_make_species_list.split("_")[0:2]).capitalize()]:
                    continue
            species_list[
                "_".join(species_to_make_species_list.split("_")[0:2]).capitalize()] = species_to_make_species_list

# species to pan homology distance
species_to_pan_homology_distance = defaultdict(lambda: defaultdict(lambda: []))
tree = Phylo.read(work_dir+'nwk/total.nwk', "newick")

for all_species in species_list:
    for pan_homology_species in pan_homology_species_list:
        if all_species == pan_homology_species:
            continue
        common_ancestor = tree.common_ancestor(pan_homology_species, all_species)

        if 'Clade()' in repr(common_ancestor):
            distance = 4290

        else:
            up = tree.get_path(all_species)[tree.get_path(all_species).index(tree.common_ancestor(all_species,pan_homology_species)):]
            up = sum([float(i.split('branch_length=')[1].split(',')[0]) for i in repr(up).split('), Clade(')])
            down = tree.get_path(all_species)[tree.get_path(all_species).index(tree.common_ancestor(all_species,pan_homology_species))+1:]
            down = sum([float(i.split('branch_length=')[1].split(',')[0]) for i in repr(down).split('), Clade(')])
            distance = (up+down)/2
        species_to_pan_homology_distance[all_species][distance].append(pan_homology_species)

for division in os.listdir(work_dir+'homology_split_by_species/'):
    if division != "pan_homology":
        for species in os.listdir(work_dir+'homology_split_by_species/' + division):
            if "_".join(species.split("_")[0:2]).capitalize() in out_group:
                if division in out_group["_".join(species.split("_")[0:2]).capitalize()]:
                    continue
            make_species_to_pan_homology_json()

