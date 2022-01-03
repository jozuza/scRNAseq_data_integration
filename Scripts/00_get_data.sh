#!/bin/bash
mkdir -p ../DataForIntegration

#Petropoulos data
mkdir -p ../RawData/EMTAB3929
wget -O ../RawData/EMTAB3929/EMTAB3929.rds http://imlspenticton.uzh.ch/robinson_lab/conquer/data-mae/EMTAB3929.rds
mkdir -p ../external_metadata/EMTAB3929
wget -O ../external_metadata/EMTAB3929/stirparo2018_tableS4.xlsx http://www.biologists.com/DEV_Movies/DEV158501/TableS4.xlsx

#Xiang data
mkdir -p ../RawData/GSE136447
wget -O ../RawData/GSE136447/GSE136447_555sample_gene_count_matrix.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136447/suppl/GSE136447%5F555sample%5Fgene%5Fcount%5Fmatrix%2Etxt%2Egz
mkdir -p ../external_metadata/GSE136447
#NOTE: Embryo_Info.txt comes from supplementary table8 supplied with the paper

#Zhou data, I got permission to download 233 cells of the full lenght dataset HRR000128

