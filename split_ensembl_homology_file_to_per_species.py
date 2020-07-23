import os

# The output see https://github.com/huanananan/GenOrigin/blob/master/homology_split_by_species.rar
work_dir = 'GenOrigin/tree/master/'

def split_homology(homology_dir, division, species_list):
    with open(homology_dir, "r") as f:
        # skip the headline and make file with head line
        headline = f.readline().strip("\n")
        for species in species_list:
            file_content = """%s""" % headline + "\n"
            outfile = open(work_dir+"homology_split_by_species/" + division + "/" + species, "w")
            outfile.write(file_content)
            outfile.close()

        for line in f:
            line = line.strip("\n")
            line = line.split("\t")
            # Only used species
            # The output see https://github.com/huanananan/GenOrigin/blob/master/homology_split_by_species.rar
            if line[2] in species_list and line[7] in species_list:
                file_content = """%s""" % "\t".join(line) + "\n"
                outfile = open(work_dir+"homology_split_by_species/" + division + "/" + line[2], "a+")
                outfile.write(file_content)
                outfile.close()

                file_content = """%s""" % "\t".join(line) + "\n"
                outfile = open(work_dir+"homology_split_by_species/" + division + "/" + line[7], "a+")
                outfile.write(file_content)
                outfile.close()




for division in ['bacteria', 'fungi', 'metazoa', 'pan_homology', 'plants', 'protists', 'vertebrates']:
    # See the file in https://github.com/huanananan/GenOrigin/tree/master/genorigin_species_list
    # Read the species list. Filter the species not used.
    species_list = []
    with open(work_dir+'genorigin_species_list/' + division) as f:
        for line in f.readlines():
            species_list.append(line.strip())

    # Creat split file dir and file

    # See the file in https://github.com/huanananan/GenOrigin/tree/master/homology_download_from_ensembl
    homology_dir = work_dir+'homology_download_from_ensembl/' + division

    split_homology(homology_dir, division, species_list)
