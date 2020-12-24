import os
import json


def get_json(json_dir):
    file = open(json_dir, 'r', encoding='UTF-8')
    get_json = json.loads(file.read())
    file.close()
    return get_json
    
    
work_dir = 'GenOrigin/tree/master/'


for species in os.listdir(work_dir+'ensembl_json/'):
    ensembl_json = get_json(work_dir+'ensembl_json/'+species):
    for gene in ensembl_json['genes']:
            file_content = """%s""" % json.dumps(gene)
            outfile = open(work_dir+'ensembl_json_split/'+species+'/'+gene['gene_id'], "w")
            outfile.write(file_content)
            outfile.close()
