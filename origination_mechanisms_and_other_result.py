import os
import subprocess
import json
from copy import deepcopy
from Bio import Phylo
from collections import defaultdict


def homology(homologyB_dir):
    file = open(homologyB_dir, 'r', encoding='UTF-8')
    homologyB = json.loads(file.read())
    file.close()
    return homologyB


# GO_term / 100mys / the range of analyse group / the drop homolog
def main(species, division):
    tree = Phylo.read("/home/lchen/genorigin/total.nwk", "newick")
    length_list = repr(tree.get_path("_".join(species.split("_")[0:2]).capitalize())).split("Clade(branch_length=")[1:]
    length_dic = {}
    for i in length_list:
        if "name" in i:
            length_dic[i.split("name='")[1].split("'")[0]] = i.split(",")[0]
        else:
            length_dic[i.split(" ")[1].split(")")[0]] = i.split(",")[0]
    age_dic = {}
    count = 0
    for i in length_dic:
        age_dic[i] = sum([float(l) for l in list(length_dic.values())[-len(length_dic) + count:]])
        count += 1
    age_dic["aa"] = 0
    less_count = 0
    Trigeminal_tree = defaultdict(lambda: defaultdict(lambda: ""))
    for i in range(0, len(length_dic)):
        if round(float(length_dic[list(length_dic.keys())[i]])) == 0:
            if "_".join(species.split("_")[0:2]).capitalize() == 'Citrobacter_freundii' or "_".join(
                    species.split("_")[0:2]).capitalize() == 'Escherichia_coli':
                Trigeminal_tree[list(length_dic.keys())[i]]["gene_age"] = int(
                    round(float(age_dic[list(age_dic.keys())[i + 1]]) + float(age_dic[list(age_dic.keys())[i]])) / 2)
                Trigeminal_tree[list(length_dic.keys())[i]]["gene_Interval"] = str(
                    round(age_dic[list(age_dic.keys())[i + 1]])) + "-" + str(round(age_dic[list(age_dic.keys())[i]]))
                Trigeminal_tree[list(length_dic.keys())[i]]["ori_branch"] = list(length_dic.keys())[i]
                Trigeminal_tree[list(length_dic.keys())[i]]["gene_branch"] = i + 1 - less_count
                continue
            less_count += 1
            Trigeminal_tree[list(length_dic.keys())[i]]["gene_age"] = int(
                round(float(age_dic[list(age_dic.keys())[i - 1]]) + float(age_dic[list(age_dic.keys())[i]])) / 2)
            Trigeminal_tree[list(length_dic.keys())[i]]["gene_Interval"] = str(
                round(age_dic[list(age_dic.keys())[i]])) + "-" + str(round(age_dic[list(age_dic.keys())[i - 1]]))
            Trigeminal_tree[list(length_dic.keys())[i]]["ori_branch"] = list(length_dic.keys())[i - 1]
            Trigeminal_tree[list(length_dic.keys())[i]]["gene_branch"] = i + 1 - less_count
        else:
            Trigeminal_tree[list(length_dic.keys())[i]]["gene_age"] = int(
                round(float(age_dic[list(age_dic.keys())[i + 1]]) + float(age_dic[list(age_dic.keys())[i]])) / 2)
            Trigeminal_tree[list(length_dic.keys())[i]]["gene_Interval"] = str(
                round(age_dic[list(age_dic.keys())[i + 1]])) + "-" + str(round(age_dic[list(age_dic.keys())[i]]))
            Trigeminal_tree[list(length_dic.keys())[i]]["ori_branch"] = list(length_dic.keys())[i]
            Trigeminal_tree[list(length_dic.keys())[i]]["gene_branch"] = i + 1 - less_count

    function_old = homology("/home/lchen/genorigin/function/" + species)
    function = defaultdict(lambda: defaultdict(lambda: []))
    function_origin_mechanism = defaultdict(lambda: defaultdict(str))

    function_GO_term = defaultdict(lambda: defaultdict(lambda: []))
    function_paralogue = defaultdict(lambda: defaultdict(lambda: []))
    function_gene_detail = defaultdict(lambda: defaultdict(lambda: []))
    function_gene_structures = defaultdict(lambda: defaultdict(lambda: []))
    function_protein = defaultdict(lambda: defaultdict(lambda: []))
    function_homolog = defaultdict(lambda: defaultdict(lambda: []))
    function_kegg = defaultdict(lambda: defaultdict(lambda: []))

    homology_web = defaultdict(lambda: [])
    homology_web_paralog = defaultdict(lambda: [])
    protein_stable_id = defaultdict(lambda: [])
    with open('/home/lchen/genorigin/homology_species/' + division + '/' + species) as f:
        for line in f.readlines():
            line = line.strip()
            line = line.split("\\t")
            if line[2] == species and line[7] != species:
                homology_web[line[0]].append(
                    {"ortholog": line[5], "homolog_perc_id": line[3], "homolog_wga_coverage": line[12],
                     "homolog_orthology_confidence": line[13], "homolog_perc_id_r1": line[8],
                     "Scientific_name": line[7], 'dn': line[9], 'ds': line[10], 'goc_score': line[11]
                     })
                protein_stable_id[line[0]] = line[1]

            if line[7] == species and line[2] != species:
                homology_web[line[5]].append(
                    {"ortholog": line[0], "homolog_perc_id": line[8], "homolog_wga_coverage": line[12],
                     "homolog_orthology_confidence": line[13], "homolog_perc_id_r1": line[3],
                     "Scientific_name": line[2], 'dn': line[9], 'ds': line[10], 'goc_score': line[11]
                     })
                protein_stable_id[line[5]] = line[6]

            if line[7] == line[2]:
                homology_web_paralog[line[0]].append(
                    {"ortholog": line[5], "homolog_perc_id": line[3], "homolog_wga_coverage": line[12],
                     "homolog_orthology_confidence": line[13], "homolog_perc_id_r1": line[8],
                     "Scientific_name": line[7], 'dn': line[9], 'ds': line[10], 'goc_score': line[11]
                     })
                homology_web_paralog[line[5]].append(
                    {"ortholog": line[0], "homolog_perc_id": line[8], "homolog_wga_coverage": line[12],
                     "homolog_orthology_confidence": line[13], "homolog_perc_id_r1": line[3],
                     "Scientific_name": line[2], 'dn': line[9], 'ds': line[10], 'goc_score': line[11]
                     })
                protein_stable_id[line[5]] = line[6]
                protein_stable_id[line[0]] = line[1]

    homologyB = homology("/home/lchen/genorigin/homology_new_json/" + species)

    annotation = homology('/home/lchen/genorigin/ensembl_json/' + species+'.json')
    organism = annotation['organism']
    genes_annotation = {}
    for i in annotation['genes']:
        genes_annotation[i['id']] = i
    annotation = ''
    GO_list = defaultdict(lambda: [])
    for gene_id in genes_annotation:
        genes = genes_annotation[gene_id]
        if 'GO' in genes:
            if genes['GO'] != []:
                for go in genes['GO']:
                    if go["evidence"] != [None]:
                        for i in go["evidence"]:
                            if i != None:
                                GO_list[go['term']].append(i)
    for i in GO_list:
        GO_list[i] = " / ".join(list(set(GO_list[i])))

    for gene_id in function_old:
        genes = genes_annotation[gene_id]
        if function_old[gene_id]['gene_branch'] == 0:
            function_old[gene_id]['gene_age'] = ">4290"
            function_old[gene_id]['gene_Interval'] = ">4290"
            function_old[gene_id]["ori_branch"] = "confidence=None"
        else:
            function_old[gene_id]['gene_age'] = Trigeminal_tree[function_old[gene_id]["ori_branch"]]['gene_age']
            function_old[gene_id]['gene_Interval'] = Trigeminal_tree[function_old[gene_id]["ori_branch"]][
                'gene_Interval']
            function_old[gene_id]['gene_branch'] = Trigeminal_tree[function_old[gene_id]["ori_branch"]]['gene_branch']
            function_old[gene_id]["ori_branch"] = Trigeminal_tree[function_old[gene_id]["ori_branch"]]["ori_branch"]

        if type(function_old[gene_id]['gene_age']) == float or type(function_old[gene_id]['gene_age']) == int:
            if function_old[gene_id]['gene_age'] <= 100:
                origin_mechanism = origin(homologyB, gene_id, species, genes, genes_annotation)
                function_origin_mechanism[gene_id]["origin_mechanism"] = origin_mechanism
                function_origin_mechanism[gene_id]["ensembl_gene_id"] = gene_id.upper()

        function[gene_id]["ensembl_gene_id"] = gene_id.upper()
        if 'name' in genes:
            function[gene_id]["external_gene_name"] = genes['name'].upper()
        else:
            function[gene_id]["external_gene_name"] = ""

        function[gene_id]["chromosome_name"] = genes['seq_region_name']
        function[gene_id]["gene_biotype"] = genes['biotype']
        function[gene_id]["start_position"] = genes['start']
        function[gene_id]["end_position"] = genes['end']
        function[gene_id]["Taxonomy_Id"] = organism['taxonomy_id']
        function[gene_id]["common_name"] = organism['display_name'].lower().capitalize()
        function[gene_id]["scientific_name"] = species.lower()
        function[gene_id]["gene_age"] = function_old[gene_id]["gene_age"]
        function[gene_id]["gene_branch"] = function_old[gene_id]["gene_branch"]
        function[gene_id]["gene_homolog"] = function_old[gene_id]["gene_homolog"]
        function[gene_id]["ori_branch"] = function_old[gene_id]["ori_branch"]
        function[gene_id]["trace"] = function_old[gene_id]["trace"]
        function[gene_id]["gene_Interval"] = function_old[gene_id]["gene_Interval"]
        function[gene_id]["protein"] = protein_stable_id[gene_id]
        function[gene_id]["division"] = function_old[gene_id]['division']
        for i in homologyB[species][gene_id]:
            if len(homologyB[species][gene_id][i]) != 0:
                function[gene_id]["gene_homolog"].append(i)
        function[gene_id]["gene_homolog"] = list(set(function[gene_id]["gene_homolog"]))

        if 'description' in genes:
            function_gene_detail[gene_id]["description"] = genes['description']
        else:
            function_gene_detail[gene_id]["description"] = ""

        function_gene_detail[gene_id]["end_position"] = genes['end']
        function_gene_detail[gene_id]["gene_biotype"] = genes['biotype']
        function_gene_detail[gene_id]["chromosome_name"] = genes['seq_region_name']
        function_gene_detail[gene_id]["Taxonomy_Id"] = genes['taxon_id']
        function_gene_detail[gene_id]["start_position"] = genes['start']
        function_gene_detail[gene_id]["ensembl_gene_id"] = gene_id.upper()
        if 'name' in genes:
            function_gene_detail[gene_id]["external_gene_name"] = genes['name'].upper()
        else:
            function_gene_detail[gene_id]["external_gene_name"] = ""

        function_gene_structures[gene_id]['ensembl_gene_id'] = gene_id.upper()
        function_gene_structures[gene_id]['Transcripts'] = []

        function_protein[gene_id]['Proteins'] = []
        function_protein[gene_id]['ensembl_gene_id'] = gene_id.upper()

        for transcript in genes['transcripts']:
            Transcripts = {"transcript_biotype": transcript['biotype'],
                           "transcript_start": transcript['start'],
                           "transcript_end": transcript['end'],
                           "transcript_length": int(transcript['end']) - int(transcript['start']) + 1,
                           "ensembl_transcript_id": transcript['id'],
                           'exons': transcript['exons'],
                           "ensembl_peptide_id": protein_stable_id[gene_id],
                           "cds_length": ""
                           }

            function_gene_structures[gene_id][
                'Transcripts'].append(Transcripts)

            if 'translations' in transcript:
                for translations in transcript['translations']:
                    if 'protein_features' in translations:
                        if translations['protein_features'] != []:
                            for fetures in translations['protein_features']:
                                if 'Pfam' not in translations:
                                    continue
                                a = {"ensembl_peptide_id": fetures['translation_id'],
                                     "pfam": translations['Pfam'][0],
                                     "pfam_end": fetures['end'],
                                     "pfam_start": fetures['start']
                                     }
                                if a not in function_protein[gene_id]['Proteins']:
                                    function_protein[gene_id]['Proteins'].append(a)

        function_homolog[gene_id]["ensembl_gene_id"] = gene_id.upper()

        function_homolog[gene_id]["homologs"] = homology_web[gene_id]

        function_paralogue[gene_id]["ensembl_gene_id"] = gene_id.upper()

        function_paralogue[gene_id]["paralog_ensembl_gene"] = homology_web_paralog[gene_id]

        function[gene_id]["GO_term"] = []
        if 'GO' in genes:
            if genes['GO'] != []:
                function_GO_term[gene_id]["ensembl_gene_id"] = gene_id.upper()
                for go in genes['GO']:
                    function[gene_id]["GO_term"].append(go["term"])

                    if GO_list[go["term"]] != [None] and GO_list[go["term"]] != []:
                        function[gene_id]["GO_term"].append(GO_list[go["term"]])
                        go['evidence'] = GO_list[go["term"]]
                    else:
                        go['evidence'] = ''

                    function_GO_term[gene_id]["GO"].append(go)

                mid = {}
                for i in function_GO_term[gene_id]["GO"]:
                    mid[i['term']] = i
                function_GO_term[gene_id]["GO"] = list(mid.values())
                function[gene_id]["GO_term"] = list(set(function[gene_id]["GO_term"]))
        if 'KEGG' in genes:
            function_kegg[gene_id]["ensembl_gene_id"] = gene_id.upper()
            function_kegg[gene_id]['kegg'] = genes['KEGG']

    function_dir = "/home/lchen/genorigin/out/" + "_".join(
        species.split("_")[0:2]).capitalize()
    file_content = '''%s''' % json.dumps(function)
    outfile = open(function_dir, "w")
    outfile.write(file_content)
    outfile.close()

    function_dir = "/home/lchen/genorigin/out/" + "_".join(
        species.split("_")[0:2]).capitalize() + "_origin_mechanism"
    file_content = '''%s''' % json.dumps(function_origin_mechanism)
    outfile = open(function_dir, "w")
    outfile.write(file_content)
    outfile.close()

    function_dir = "/home/lchen/genorigin/out/" + "_".join(
        species.split("_")[0:2]).capitalize() + "_paralogue"
    file_content = '''%s''' % json.dumps(function_paralogue)
    outfile = open(function_dir, "w")
    outfile.write(file_content)
    outfile.close()

    function_dir = "/home/lchen/genorigin/out/" + "_".join(
        species.split("_")[0:2]).capitalize() + "_protein"
    file_content = '''%s''' % json.dumps(function_protein)
    outfile = open(function_dir, "w")
    outfile.write(file_content)
    outfile.close()

    function_dir = "/home/lchen/genorigin/out/" + "_".join(
        species.split("_")[0:2]).capitalize() + "_GO_term"
    file_content = '''%s''' % json.dumps(function_GO_term)
    outfile = open(function_dir, "w")
    outfile.write(file_content)
    outfile.close()

    if species != 'homo_sapiens' and species != 'mus_musculus':
        function_dir = "/home/lchen/genorigin/out/" + "_".join(
            species.split("_")[0:2]).capitalize() + "_kegg"
        file_content = '''%s''' % json.dumps(function_kegg)
        outfile = open(function_dir, "w")
        outfile.write(file_content)
        outfile.close()

    function_dir = "/home/lchen/genorigin/out/" + "_".join(
        species.split("_")[0:2]).capitalize() + "_gene_structures"
    file_content = '''%s''' % json.dumps(function_gene_structures)
    outfile = open(function_dir, "w")
    outfile.write(file_content)
    outfile.close()

    function_dir = "/home/lchen/genorigin/out/" + "_".join(
        species.split("_")[0:2]).capitalize() + "_gene_detail"
    file_content = '''%s''' % json.dumps(function_gene_detail)
    outfile = open(function_dir, "w")
    outfile.write(file_content)
    outfile.close()
    function_dir = "/home/lchen/genorigin/out/" + "_".join(
        species.split("_")[0:2]).capitalize() + "_homolog"
    file_content = '''%s''' % json.dumps(function_homolog)
    outfile = open(function_dir, "w")
    outfile.write(file_content)
    outfile.close()

    GO_tool = defaultdict(lambda: defaultdict(lambda: []))
    for gene_id in function_GO_term:
        age = function[gene_id]["gene_age"]
        for item in function_GO_term[gene_id]["GO"]:
            if type(GO_tool[item["term"]]["age"]) == list:
                GO_tool[item["term"]]["age"] = defaultdict(lambda: [])

            if item["evidence"] not in GO_tool[item["term"]]["evidence"]:
                if item["evidence"] != "":
                    GO_tool[item["term"]]["evidence"].append(item["evidence"])

            GO_tool[item["term"]]["age"][age].append(gene_id)

    sort_list = defaultdict(lambda: [])
    for term in GO_tool:
        total = sum([len(GO_tool[term]["age"][age]) for age in GO_tool[term]["age"]]) / 2
        mid = []
        for i in GO_tool[term]["age"].keys():
            if type(i) == str:
                pass
            else:
                mid.append(i)

        mid = sorted(mid)
        mid.append('>4290')
        count = 0
        for i in mid:
            count += len(GO_tool[term]["age"][i])
            if count > total:

                if i == ">4290":
                    i = 4290

                sort_list[i].append(term)
                break
            else:
                pass

    out_list = []
    for i in sorted(sort_list.keys()):
        for term in sort_list[i]:
            y_axis = []
            for age in GO_tool[term]["age"]:
                if type(age) == str:
                    age_use = 4290
                else:
                    age_use = age
                y_axis.append({
                    "gene_age": age_use,
                    "size": len(list(set(GO_tool[term]["age"][age])))
                    # "ensembl_gene_id": [0 for i in range(0,len(list(set(GO_tool[term]["age"][age]))))],
                }
                )
            GO_name = []
            for i in GO_tool[term]["evidence"]:

                if i == [None]:
                    continue

                GO_name.append(i)
            hhh = []
            for i in GO_name:
                if i == None:
                    continue
                else:
                    hhh.append(i)
            GO_name = hhh
            GO_name = list(set(GO_name))
            print(GO_name)
            if len(GO_name) > 1:
                GO_name = " / ".join(GO_name)
            elif len(GO_name) == 0:
                GO_name = ""
            else:
                GO_name = GO_name[0]
            out_list.append({
                "evidence": term,
                "species": species,
                "GO_name": GO_name,
                "y_axis": y_axis
            })

    function_dir = "/home/lchen/genorigin/out/" + "_".join(
        species.split("_")[0:2]).capitalize() + "_GO_tool"
    outfile = open(function_dir, "w")
    outfile.write(json.dumps(out_list))
    outfile.close()


def origin(homologyB, gene_id, species, genes, genes_annotation):
    paralogues = homologyB[species][gene_id]["_".join(species.split("_")[0:2]).capitalize()]
    if len(paralogues) == 1:
        origin_mechanism = "de novo origination"
        paralogues = ""
        return [{"parental_gene": gene_id, "origin_mechanism": origin_mechanism}]
    else:

        child_gene_exon = []
        try:
            transcripts = genes["transcripts"]
        except:
            return [{"parental_gene": gene_id, "origin_mechanism": "NA"}]
        for transcript_id in transcripts:
            child_gene_exon.append(len(transcript_id['exons']))

        origin_list = []

        for gene_paralogues in paralogues:
            if gene_paralogues != gene_id:
                try:
                    parental_gene_exon = []
                    transcripts = genes_annotation[gene_paralogues]
                    if "transcripts" in transcripts:
                        for transcript_id in transcripts["transcripts"]:
                            parental_gene_exon.append(len(transcript_id['exons']))
                    if max(child_gene_exon) == 1 and max(parental_gene_exon) > 1:
                        origin_mechanism = "retroposition"
                        origin_list.append({"parental_gene": gene_paralogues, "origin_mechanism": origin_mechanism})
                    elif max(child_gene_exon) > 1 and max(parental_gene_exon) > 1:
                        origin_mechanism = "DNA-based duplication"
                        origin_list.append({"parental_gene": gene_paralogues, "origin_mechanism": origin_mechanism})
                except:
                    pass
        if len(origin_list) == 0:
            return [{"parental_gene": gene_id, "origin_mechanism": "NA"}]
        return origin_list
main(division,species)
