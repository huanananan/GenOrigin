import json
from collections import defaultdict


work_dir = 'GenOrigin/tree/master/'
# homologyB structure
# homologyB:species:gene:homology_species:gene
def pre(homology_dir, species_list, species):
    # homology B:species:gene ID:"_".join(homology_species.split("_")[0:2]).capitalize():homology species gene ID
    # gene ID：[species]
    homologyB = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: [])))
    # paralogues pairs
    paralogues_list = []
    # paralogues gene
    collcetion_list = []
    print(species, "pre")
    with open(homology_dir, "r") as f:
        next(f)
        for line in f.readlines():
            line = line.strip("\n")
            line = line.split("\t")
            # print(line[0],"1",species)
            name = "_".join(line[2].split("_")[0:2]).capitalize()
            homology_name = "_".join(line[7].split("_")[0:2]).capitalize()
            # paralogues

            # if the gene is paralog
            if line[2] == line[7] and line[2] == species:
                # creat keys for each species
                if line[0] not in homologyB[line[2]]:
                    homologyB[line[2]][line[0]] = dict(
                        zip(["_".join(i.split("_")[0:2]).capitalize() for i in species_list],
                            [[] for i in range(0, len(species_list))]))

                # creat keys for each species
                if line[5] not in homologyB[line[7]]:
                    homologyB[line[7]][line[5]] = dict(
                        zip(["_".join(i.split("_")[0:2]).capitalize() for i in species_list],
                            [[] for i in range(0, len(species_list))]))

                homologyB[line[2]][line[0]][name].append(line[0])
                homologyB[line[2]][line[0]][name].append(line[5])

                homologyB[line[2]][line[5]][name].append(line[0])
                homologyB[line[2]][line[5]][name].append(line[5])

            # ortholog
            elif name != homology_name:
                if line[2] == species:
                    if len(homologyB[line[2]][line[0]]) == 0:
                        for i in species_list:
                            homologyB[line[2]][line[0]]["_".join(i.split("_")[0:2]).capitalize()] = []

                    if line[0] not in homologyB[line[2]][line[0]][name]:
                        homologyB[line[2]][line[0]][name].append(line[0])

                    if line[5] not in homologyB[line[2]][line[0]][homology_name]:
                        homologyB[line[2]][line[0]][homology_name].append(line[5])

                elif line[7] == species:
                    if len(homologyB[line[7]][line[5]]) == 0:
                        for i in species_list:
                            homologyB[line[7]][line[5]]["_".join(i.split("_")[0:2]).capitalize()] = []

                    if line[0] not in homologyB[line[7]][line[5]][name]:
                        homologyB[line[7]][line[5]][name].append(line[0])

                    if line[5] not in homologyB[line[7]][line[5]][homology_name]:
                        homologyB[line[7]][line[5]][homology_name].append(line[5])
    return homologyB


# Supplement the homology list for pan-taxon compara species.
# Not contain the paralog.
def pan_homology(homologyB, species_list, species):
    # See the file in https://github.com/huanananan/GenOrigin/tree/master/genorigin_species_list
    with open(work_dir + 'genorigin_species_list/pan_species') as f:
        pan_list = [line.strip() for line in f.readlines()]

    # See the file in https://github.com/huanananan/GenOrigin/blob/master/homology_split_by_species.rar
    homology_dir = work_dir + 'homology_split_by_species/pan_homology/' + species
    if species in pan_list:
        with open(homology_dir, "r") as f:
            next(f)
            for line in f.readlines():
                line = line.strip("\n")
                line = line.split("\t")
                # print(line[0],"1",species)
                name = "_".join(line[2].split("_")[0:2]).capitalize()
                homology_name = "_".join(line[7].split("_")[0:2]).capitalize()
                # prologues
                if line[2] == line[7] and line[2] == species:
                    if len(homologyB[line[2]][line[0]]) == 0:
                        homologyB[line[2]][line[0]] = dict(
                            zip(["_".join(i.split("_")[0:2]).capitalize() for i in species_list],
                                [[] for i in range(0, len(species_list))]))

                    if len(homologyB[line[7]][line[5]]) == 0:
                        homologyB[line[7]][line[5]] = dict(
                            zip(["_".join(i.split("_")[0:2]).capitalize() for i in species_list],
                                [[] for i in range(0, len(species_list))]))

                    if line[0] not in homologyB[line[2]][line[0]][name]:
                        homologyB[line[2]][line[0]][name].append(line[0])

                    if line[5] not in homologyB[line[2]][line[5]][name]:
                        homologyB[line[2]][line[5]][name].append(line[5])
                # homology
                elif name != homology_name:
                    if line[2] == species:
                        if len(homologyB[line[2]][line[0]]) == 0:
                            for i in species_list:
                                homologyB[line[2]][line[0]]["_".join(i.split("_")[0:2]).capitalize()] = []

                        if line[0] not in homologyB[line[2]][line[0]][name]:
                            homologyB[line[2]][line[0]][name].append(line[0])

                        if line[5] not in homologyB[line[2]][line[0]][homology_name]:
                            homologyB[line[2]][line[0]][homology_name].append(line[5])

                    elif line[7] == species:
                        if len(homologyB[line[7]][line[5]]) == 0:
                            for i in species_list:
                                homologyB[line[7]][line[5]]["_".join(i.split("_")[0:2]).capitalize()] = []

                        if line[0] not in homologyB[line[7]][line[5]][name]:
                            homologyB[line[7]][line[5]][name].append(line[0])

                        if line[5] not in homologyB[line[7]][line[5]][homology_name]:
                            homologyB[line[7]][line[5]][homology_name].append(line[5])
        f.close()

    return homologyB


def run_above(species, division):
    # See the file in https://github.com/huanananan/GenOrigin/blob/master/homology_split_by_species.rar
    homology_dir = work_dir + 'homology_split_by_species'+division+'/' + species

    # See the file in https://github.com/huanananan/GenOrigin/tree/master/genorigin_species_list
    with open(work_dir + 'genorigin_species_list/'+division) as f:
        species_list = [line.strip() for line in f.readlines()]

    homologyB = pre(homology_dir, species_list, species)
    if division != 'pan_homology':
        homologyB = pan_homology(homologyB, species_list, species)

    # write the json file
    # See the file in https://github.com/huanananan/GenOrigin/tree/master/homology_json.rar
    function_dir = work_dir + 'homology_json/'+division+'/' + species
    file_content = """%s""" % json.dumps(homologyB)
    outfile = open(function_dir, "w")
    outfile.write(file_content)
    outfile.close()
