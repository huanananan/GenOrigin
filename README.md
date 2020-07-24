# [GenOrigin](http://genorigin.chenzxlab.cn/#!/)
web: http://genorigin.chenzxlab.cn/#!/


## Step 1 ensembl homology data
Download the homology data from ensembl FTP. release 98/45.  
[plants](ftp://ftp.ensemblgenomes.org/pub/plants/release-45/tsv/ensembl-compara/homologies/)  
[fungi](ftp://ftp.ensemblgenomes.org/pub/fungi/release-45/tsv/ensembl-compara/homologies/)  
[protists](ftp://ftp.ensemblgenomes.org/pub/protists/release-45/tsv/ensembl-compara/homologies/)  
[bacteria](ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/tsv/ensembl-compara/homologies/)  
[metazoa](ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/tsv/ensembl-compara/homologies/)  
[vertebrate](ftp://ftp.ensembl.org/pub/release-98/tsv/ensembl-compara/homologies/)  

## Step 2 ensembl annotation
Download the annotation json file from ensembl FTP  
eg:  
ftp://ftp.ensemblgenomes.org/pub/plants/release-45/json/  
ftp://ftp.ensembl.org/pub/release-98/json  

## Step 3 time-tree
Use the total species list from both [ensembl](www.ensembl.org), [ensemblgenomes](www.ensemblgenomes.org) (GenOrigin only provide the used species).   
Then upload species list to [timetree](www.timetree.org).   
Finally,we get the timetree file.   

## Step 4 uniform time-tree scientific name
The scientific name in timetree.org and ensembl.org might be different.  
But it show the same taxonomy id in [uniprot](www.uniprot.org).  
So we uniform it to ensembl scientific name first two word.  
See [change_time_tree_nwk_file_to_ensembl_species_name.txt](https://github.com/huanananan/GenOrigin/blob/master/change_time_tree_nwk_file_to_ensembl_species_name.txt)  

## Step 5 split the ensembl homology data
Using the python script [split_ensembl_homology_file_to_per_species.py](https://github.com/huanananan/GenOrigin/blob/master/split_ensembl_homology_file_to_per_species.py) splits the homology data by species.  
The output file, see [homology_split_by_species.rar](https://github.com/huanananan/GenOrigin/blob/master/homology_split_by_species.rar)  
Split it could save memory and improve parallel efficiency.  
Some species even split it as 100 genes to run pipline.  

## Step 6 split the annotation json file (Optional)
Using the python script [split_ensembl_homology_file_to_per_species.py](https://github.com/huanananan/GenOrigin/blob/master/split_ensembl_homology_file_to_per_species.py) split it by gene_id, organism or assembly.  
The output file, see  [split_ensembl_homology_file_to_per_species.py](https://github.com/huanananan/GenOrigin/blob/master/split_ensembl_homology_file_to_per_species.py).  
Split it could save memory, but not improve parallel efficiency.  
Because reading and writing a large number of files at the same time will seriously slow down the speed.  
But, it use disk place replace memory.  
The biggest json file is about 44GB memory.  
So, we provide and use the split annotation json in this project.  
You can easy to change it, because i will provide the convert method in the python script as comment.  
