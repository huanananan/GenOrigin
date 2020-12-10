# [GenOrigin](http://genorigin.chenzxlab.cn/#!/)
web: http://genorigin.chenzxlab.cn/#!/

#Collect data 
## Step 1 Ensembl homology data
Download the homology data from Ensembl FTP. release 98/45.  
ftp://ftp.ensemblgenomes.org/pub/plants/release-45/tsv/ensembl-compara/homologies/  
ftp://ftp.ensemblgenomes.org/pub/fungi/release-45/tsv/ensembl-compara/homologies/  
ftp://ftp.ensemblgenomes.org/pub/protists/release-45/tsv/ensembl-compara/homologies/  
ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/tsv/ensembl-compara/homologies/  
ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/tsv/ensembl-compara/homologies/  
ftp://ftp.ensembl.org/pub/release-98/tsv/ensembl-compara/homologies/  

## Step 2 Ensembl annotation
Download the annotation json file from Ensembl FTP  
eg:  
ftp://ftp.ensemblgenomes.org/pub/plants/release-45/json/  
ftp://ftp.ensembl.org/pub/release-98/json  

## Step 3 time-tree
Use the total species list from both [Ensembl](www.ensembl.org), [EnsemblGenomes](www.ensemblgenomes.org) (GenOrigin only provide the used species).   
Then upload species list to [TimeTree](www.timetree.org).   
Finally, we get the [timetree file](https://github.com/huanananan/GenOrigin/tree/master/nwk).   

## Step 4 Count
[Count download](https://www.iro.umontreal.ca/~csuros/gene_content/count.html)  
param -gain 1.4

## Step 5 uniform time-tree scientific name
The scientific name in [TimeTree](timetree.org) and [Ensembl](ensembl.org) might be different.  
But it show the same taxonomy id in [UniProt](www.uniprot.org).  
So we uniform it to Ensembl scientific name first two word.  
See [change_time_tree_nwk_file_to_ensembl_species_name.txt](https://github.com/huanananan/GenOrigin/blob/master/change_time_tree_nwk_file_to_ensembl_species_name.txt)  

## Step 6 split the Ensembl homology data & make it to json file
Using the python script [split_ensembl_homology_file_to_per_species.py](https://github.com/huanananan/GenOrigin/blob/master/split_ensembl_homology_file_to_per_species.py) splits the homology data by species.  
The output file, see [homology_split_by_species.rar](https://github.com/huanananan/GenOrigin/blob/master/homology_split_by_species.rar)  
Also [homology_split_by_species](https://github.com/huanananan/GenOrigin/tree/master/homology_split_by_species)
Split it could save memory and improve parallel efficiency.  
Some species even split it as 100 genes to run pipeline.  
Then, make it to the json file. [trans_homology_tsv_to_json.py](https://github.com/huanananan/GenOrigin/blob/master/trans_homology_tsv_to_json.py)

## Step 7 split the annotation json file (Optional)
Using the python script [split_ensembl_homology_file_to_per_species.py](https://github.com/huanananan/GenOrigin/blob/master/split_ensembl_homology_file_to_per_species.py) split it by gene id, organism or assembly.  
The output file, see [split_ensembl_homology_file_to_per_species.py](https://github.com/huanananan/GenOrigin/blob/master/split_ensembl_homology_file_to_per_species.py).  
Split it could save memory, but not improve parallel efficiency.  
Because reading and writing a large number of files at the same time will seriously slow down the speed.  
But, it use disk place replace memory.  
The biggest json file is about 44GB memory.  
So, we provide and use the split annotation json in this project.  

#Infer origin
## Step 8 species connect to pan-taxonomy compara
Using

## Step 9 pan-taxonomy species to pan-taxonomy species

## Step 10 pan-taxonomy species extension

## step 11 other species extension
