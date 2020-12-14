import os
import json


def read_json(json_dir):
    file = open(json_dir, 'r', encoding='UTF-8')
    dic = json.loads(file.read())
    file.close()
    return dic


# not contain the outgroup
def make_division_extension(division):
    division_homology = {}
    for species in os.listdir(work_dir + 'homology_json/' + division):
        if species in out_group:
            if division in out_group:
                continue

        division_homology["_".join(species.split('_')[0:2])] = read_json(
            work_dir + 'homology_json/' + division + '/' + species)[species]

    function_dir = work_dir + "division_json/" + division
    file_content = """%s""" % json.dumps(division_homology)
    outfile = open(function_dir, "w")
    outfile.write(file_content)
    outfile.close()


out_group = {"Ustilago_maydis": ["bacteria"], "Saccharomyces_cerevisiae": ["plants", "vertebrates"],
             "Homo_sapiens": ["plants"], "Drosophila_melanogaster": ["plants", "vertebrates"],
             "Caenorhabditis_elegans": ["plants", "vertebrates"]}


work_dir = 'GenOrigin/tree/master/'

for division in os.listdir(work_dir + 'homology_json/'):
    if division == "pan_homology":
        continue
    make_division_extension(division)

