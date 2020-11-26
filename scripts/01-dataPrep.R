library(data.table)
library(tidyverse)

att = read_tsv('./data/raw/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt')
phe = read_tsv('./data/raw/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt')
att$SUBJID = sapply(strsplit(att$SAMPID,'-'),function(x)paste(x[1],x[2],sep='-'))
phe$SEX = c('male','female')[phe$SEX]
phe = phe %>%
  set_names(c('id','sex','age','death')) %>%
  mutate(sex = as.factor(sex),
         age = as.factor(age))
allattr = att %>%
  select(SAMPID,SMTS,SMTSD,SMNABTCH,SMGEBTCH,SUBJID,SMRIN,SMTSISCH) %>%
  set_names(c('sample_id','major_tissue','minor_tissue','batch1','batch2','id','rin','ischemic_time')) %>%
  unique() %>%
  full_join(phe)

allattr = allattr %>%
  filter(death %in% c(1,2))

dat=fread('./data/raw/expression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct',select = c('Name',allattr$sample_id))
samplesx=colnames(dat)[-1]

dat = as.data.frame(dat)
rownames(dat) = dat$Name
dat$Name = NULL
dat = as.matrix(dat)
samplesx = intersect(allattr$sample_id, colnames(dat))
samples_by_tissues = tapply(as.character(allattr$sample_id),INDEX=allattr$minor_tissue,FUN = function(x)unique(c(x)))
dat = lapply(samples_by_tissues,function(samps){
  samps = intersect(samps,samplesx)
  dat[,samps]
})

system('mkdir -p data/processed/expression/tpm/')
names(dat) = sapply(strsplit(gsub(' ','',names(dat)),'[(]'),function(x)x[[1]])

sapply(names(dat),function(nm){
  saveRDS(dat[[nm]],paste('data/processed/expression/tpm/',nm,'.rds',sep=''))
})
saveRDS(allattr,'./data/processed/attr.rds')
