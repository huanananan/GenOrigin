import subprocess
import os
import requests
import time
from threading import Thread
from concurrent.futures import ThreadPoolExecutor


def vet(hh):
    species = hh[0]

    web = 'http://ftp.ensembl.org/pub/release-98/gff3/'
    page = requests.Session().get(web + species).text
    # time.sleep(1)
    # print(page)
    for line in page.split('\n'):
        if species.capitalize() in line and '.98.gff3' in line:
            ftp = web + species + '/' + line.split('<a href="')[1].split('">')[0]
            ftp = ftp.replace('http', 'ftp')
            ftp = ftp.replace('vol1/', '')
    print(species, '\t', 'vertebrates', '\t', ftp)


def collection11(hh):
    division, species = hh
    web = 'http://ftp.ensemblgenomes.org/vol1/pub/release-45/' + division + '/gff3/'
    if species in species_collection:
        page = requests.Session().get(web + species_collection[species] + species).text
        # print(page)
        # time.sleep(1)
        for line in page.split('\n'):
            if (species.capitalize() in line and '.37.gff3' in line) or (
                    species.capitalize() in line and '.45.gff3' in line):
                ftp = web + species_collection[species] + species + '/' + line.split('<a href="')[1].split('">')[0]
                ftp = ftp.replace('http', 'ftp')
                ftp = ftp.replace('vol1/', '')

        print(species, '\t', division, '\t', ftp)


    else:
        page = requests.Session().get(web + species).text
        # print(page)
        # time.sleep(1)
        for line in page.split('\n'):
            if (species.capitalize() in line and '.37.gff3' in line) or (
                    species.capitalize() in line and '.45.gff3' in line):
                ftp = web + species + '/' + line.split('<a href="')[1].split('">')[0]
                ftp = ftp.replace('http', 'ftp')
                ftp = ftp.replace('vol1/', '')
        print(species, '\t', division, '\t', ftp)


def other(hh):
    division, species = hh
    web = 'http://ftp.ensemblgenomes.org/vol1/pub/release-45/' + division + '/gff3/'
    page = requests.Session().get(web + species).text
    time.sleep(1)
    for line in page.split('\n'):
        if (species.capitalize() in line and '.37.gff3' in line) or (
                species.capitalize() in line and '.45.gff3' in line):
            ftp = web + species + '/' + line.split('<a href="')[1].split('">')[0]
            ftp = ftp.replace('http', 'ftp')
            ftp = ftp.replace('vol1/', '')
    print(species, '\t', division, '\t', ftp)


out_group = {"Ustilago_maydis": ["bacteria"], "Saccharomyces_cerevisiae": ["plants", "vertebrates"],
             "Homo_sapiens": ["plants"], "Drosophila_melanogaster": ["plants", "vertebrates"],
             "Caenorhabditis_elegans": ["plants", "vertebrates"]}
os.chdir('/home/lchen/genorigin/ensembl_gff3/')
for division in os.listdir('/home/lchen/genorigin/homology_species/'):
    if division == 'pan_homology':
        continue

    if division == 'vertebrates':
        # threadPool = ThreadPoolExecutor(max_workers=12, thread_name_prefix="test_")
        for species in os.listdir('/home/lchen/genorigin/homology_species/' + division + '/'):
            if "_".join(species.split("_")[0:2]).capitalize() in out_group:
                if division in out_group["_".join(species.split("_")[0:2]).capitalize()]:
                    continue

            vet([species])
            # future = threadPool.submit(vet, [species])
            # subprocess.getoutput('wget ' + ftp)

    elif division == 'bacteria' or division == 'fungi' or division == 'protists':
        # if division == 'fungi':
        species_list = os.listdir('/home/lchen/genorigin/homology_species/' + division + '/')
        web = 'http://ftp.ensemblgenomes.org/vol1/pub/release-45/' + division + '/gff3/'
        page = requests.Session().get(web).text
        # time.sleep(1)
        species_collection = {}
        for line in page.split('\n'):
            if division in line and 'collection' in line:
                collection = line.split('<a href="')[1].split('">')[0]
                collection_page = requests.Session().get(web + collection).text.split('\n')
                # time.sleep(1)
                for line in collection_page:
                    if '<a href="' in line and line.split('<a href="')[1].split('/">')[0] in species_list:
                        species_collection[line.split('<a href="')[1].split('/">')[0]] = collection

        # threadPool = ThreadPoolExecutor(max_workers=12, thread_name_prefix="test_")
        for species in species_list:
            if "_".join(species.split("_")[0:2]).capitalize() in out_group:
                if division in out_group["_".join(species.split("_")[0:2]).capitalize()]:
                    continue
            collection11([division, species])
            # future = threadPool.submit(collection11, [division,species])

    else:
        # threadPool = ThreadPoolExecutor(max_workers=12, thread_name_prefix="test_")
        for species in os.listdir('/home/lchen/genorigin/homology_species/' + division + '/'):
            if "_".join(species.split("_")[0:2]).capitalize() in out_group:
                if division in out_group["_".join(species.split("_")[0:2]).capitalize()]:
                    continue

            other([division, species])
            # future = threadPool.submit(other, [division,species])
