system('mkdir -p data/raw/expression')
system('mkdir -p data/processed/expression')
system('mkdir -p results')
system('mkdir -p scripts')
library(tidyverse)
system('wget -nv https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx -P data/raw/')
system('wget -nv https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDD.xlsx -P data/raw/')
system('wget -nv https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt -P data/raw/')
system('wget -nv https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt -P data/raw/')
system('wget -nv https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz -P data/raw/expression')
system('cd data/raw/expression; gunzip *')


