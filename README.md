# [GenOrigin v2.01](http://genorigin.chenzxlab.cn/#!/)
web: http://genorigin.chenzxlab.cn/#!/  
  * [Collect data](#collect-data)
    + [Step 1 Ensembl homology data](#step-1-ensembl-homology-data)
    + [Step 2 Ensembl annotation](#step-2-ensembl-annotation)
    + [Step 3 time-tree](#step-3-time-tree)
    + [Step 4 Count](#step-4-count)
    + [Step 5 uniform time-tree scientific name](#step-5-uniform-time-tree-scientific-name)
    + [Step 6 split the Ensembl homology data & make it to json file](#step-6-split-the-ensembl-homology-data-&-make-it-to-json-file)
    + [Step 7 split the annotation json file](#step-7-split-the-annotation-json-file)
  * [Infer origin](#infer-origin)
    + [Step 8 species connect to pan-taxonomy compara](#step-8-species-connect-to-pan-taxonomy-compara)
    + [Step 9 make division homology json file](#step-9-make-division-homology-json-file)
    + [Step 10 pan-taxonomy species extension](#step-10-pan-taxonomy-species-extension)
    + [Step 11 other species extension](#step-11-other-species-extension)
    + [Step 12 origination mechanisms and other result](#step-12-origination-mechanisms-and-other-result)

## Collect data 
### Step 1 Ensembl homology data
Download the homology data from Ensembl FTP. release 98/45.  
ftp://ftp.ensemblgenomes.org/pub/plants/release-45/tsv/ensembl-compara/homologies/  
ftp://ftp.ensemblgenomes.org/pub/fungi/release-45/tsv/ensembl-compara/homologies/  
ftp://ftp.ensemblgenomes.org/pub/protists/release-45/tsv/ensembl-compara/homologies/  
ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/tsv/ensembl-compara/homologies/  
ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/tsv/ensembl-compara/homologies/  
ftp://ftp.ensembl.org/pub/release-98/tsv/ensembl-compara/homologies/  

### Step 2 Ensembl annotation
1)[from_ensembl_FTP_get_the_gff3_dir.py](https://github.com/huanananan/GenOrigin/blob/master/from_ensembl_FTP_get_the_gff3_dir.py)  

Download the annotation json file from Ensembl FTP:  
eg:  
ftp://ftp.ensemblgenomes.org/pub/plants/release-45/json/  
ftp://ftp.ensembl.org/pub/release-98/json  

Download the gff3 file from Ensembl FTP:  
[from_ensembl_FTP_get_the_gff3_dir.py](https://github.com/huanananan/GenOrigin/blob/master/from_ensembl_FTP_get_the_gff3_dir.py)  

### Step 3 time-tree
Use the total species list from both [Ensembl](http://www.ensembl.org/), [EnsemblGenomes](http://www.ensemblgenomes.org/) (GenOrigin only provide the used species).   
Then upload species list to [TimeTree](http://www.timetree.org/).   
Finally, we get the [timetree file](https://github.com/huanananan/GenOrigin/tree/master/nwk).   

### Step 4 Count
1)[Count download](https://www.iro.umontreal.ca/~csuros/gene_content/count.html)  
param -gain 1.4

### Step 5 uniform time-tree scientific name
1)[change_time_tree_nwk_file_to_ensembl_species_name.txt](https://github.com/huanananan/GenOrigin/blob/master/change_time_tree_nwk_file_to_ensembl_species_name.txt)  

The scientific name in [TimeTree](timetree.org) and [Ensembl](ensembl.org) might be different.  
But it show the same taxonomy id in [UniProt](www.uniprot.org).  
So we uniform it to Ensembl scientific name first two word.  
See [change_time_tree_nwk_file_to_ensembl_species_name.txt](https://github.com/huanananan/GenOrigin/blob/master/change_time_tree_nwk_file_to_ensembl_species_name.txt)  

### Step 6 split the Ensembl homology data & make it to json file
1)[split_ensembl_homology_file_to_per_species.py](https://github.com/huanananan/GenOrigin/blob/master/split_ensembl_homology_file_to_per_species.py)  

2)[homology_split_by_species.rar](https://github.com/huanananan/GenOrigin/blob/master/homology_split_by_species.rar)   

3)[homology_split_by_species](https://github.com/huanananan/GenOrigin/tree/master/homology_split_by_species)  

4)[trans_homology_tsv_to_json.py](https://github.com/huanananan/GenOrigin/blob/master/trans_homology_tsv_to_json.py)  

Using the python script [split_ensembl_homology_file_to_per_species.py](https://github.com/huanananan/GenOrigin/blob/master/split_ensembl_homology_file_to_per_species.py) splits the homology data by species.  
The output file, see [homology_split_by_species.rar](https://github.com/huanananan/GenOrigin/blob/master/homology_split_by_species.rar)  
Also [homology_split_by_species](https://github.com/huanananan/GenOrigin/tree/master/homology_split_by_species)  
Split it could save memory and improve parallel efficiency.    
Some species even split it as 100 genes to run pipeline.    
Then, make it to the json file. [trans_homology_tsv_to_json.py](https://github.com/huanananan/GenOrigin/blob/master/trans_homology_tsv_to_json.py)  

### Step 7 split the annotation json file
1)[split_ensembl_annotation_file_to_per_gene.py](https://github.com/huanananan/GenOrigin/blob/master/split_ensembl_annotation_file_to_per_gene.py)   

Using the python script [split_ensembl_annotation_file_to_per_gene.py](https://github.com/huanananan/GenOrigin/blob/master/split_ensembl_annotation_file_to_per_gene.py)   split it by gene id, organism or assembly.      
Split it could save memory, but not improve parallel efficiency.  
Because reading and writing a large number of files at the same time will seriously slow down the speed.  
But, it use disk place replace memory.  
The biggest json file is about 44GB memory.  
So, we provide and use the split annotation json in this project.  

## Infer origin
### Step 8 species connect to pan-taxonomy compara
1)[species2pan-taxonomy_compara.py](https://github.com/huanananan/GenOrigin/blob/master/species2pan-taxonomy_compara.py)  

Using [species2pan-taxonomy_compara.py](https://github.com/huanananan/GenOrigin/blob/master/species2pan-taxonomy_compara.py) to make the json file, which can lead the gene to its represent gene in the pan-taxonomy compara with the closest distance species on the time-tree.  

### Step 9 make division homology json file
1)[make_division_homology_json.py](https://github.com/huanananan/GenOrigin/blob/master/make_division_homology_json.py)  

Using [make_division_homology_json.py](https://github.com/huanananan/GenOrigin/blob/master/make_division_homology_json.py) to make the division homology json file for next step.
This script merge the homology json file, which, the species, in the same division, to one file.  
Only the pan-taxonomy species have been merged.  
It also can use the split one, but the json file saving time.  

### Step 10 pan-taxonomy species extension
1)[pan-taxonomy_species_extension.py](https://github.com/huanananan/GenOrigin/blob/master/pan-taxonomy_species_extension.py)  

Using the pan-taxonomy gene (gene A) to find the representative pan-taxonomy homology genes (for example gene B, gene C, gene D).  
Then using the representative pan-taxonomy homology genes (gene B, gene C, gene D) to find the representative own domain's (as the gene B from fungi) homology genes (gene E, gene F, gene G), but not pan-taxonomy gene (gene A)'s domain (as gene A from vetebrates).  
Finally, as above, adding the outgroup for some domain'species.

### Step 11 other species extension
1)[All_species-ortholog-new_homologyB.py](https://github.com/huanananan/GenOrigin/blob/master/All_species-ortholog-new_homologyB.py)  

Using the pan-taxonomy_species_extension result and species2pan-taxonomy_compara result to build a huge homology tree.  
Frist, all species should be infer the age in their own analyse group, generally which is its domain, but the plants have two different analyse group because the fungi split them on the Phylogenetic tree.  
Second, those gene break out the analyse group show use the species2pan-taxonomy_compara result connect to the representative gene, then use the pan-taxonomy_species_extension result to extension its own homology info.  

### Step 12 origination mechanisms and other result
1)[origination_mechanisms_and_other_result.py](https://github.com/huanananan/GenOrigin/blob/master/origination_mechanisms_and_other_result.py)  

Using the CDS exon and two protein sequence identity to judge the potential origination mechanism.  

