library(tidyverse)
library(ggpubr)
qc_check = function(mytis, outliers = NA,round = 0) {
  att = readRDS('./data/processed/attr.rds')
  print(mytis)
  tpm = readRDS(paste('./data/processed/expression/tpm/',mytis,sep=''))
  sampnames=setdiff(intersect(colnames(tpm),att$sample_id),outliers)
  if(length(sampnames)>=10){
    tisname = gsub('.rds','',mytis)
  rownames(tpm) = sapply(strsplit(rownames(tpm),'[.]'),function(x)x[1])
  
  att = att %>% filter(sample_id %in% sampnames)
  tpm = tpm[,sampnames]
    sex = setNames(att$sex,att$sample_id)
    sex = sex[colnames(tpm)]
    death = setNames(att$death, att$sample_id)
    death = death[colnames(tpm)]
    medtpm = apply(tpm,1,median)
    tpm = tpm[names(which(medtpm>1)),]
    tpm = log2(tpm+1)
    if(length(unique(sex)) == 1 & length(unique(death)) == 1) {
      tpm = tpm
    } else if(length(unique(sex)) == 1){
      tpm = t(apply(tpm,1,function(x){
        lm(x ~ death)$resid
      }))
    } else if(length(unique(death)) == 1){
      tpm = t(apply(tpm,1,function(x){
        lm(x ~ sex )$resid
      }))
    } else {
      tpm = t(apply(tpm,1,function(x){
        lm(x ~ sex + death)$resid
      }))
    }
    tpm_qn = preprocessCore::normalize.quantiles(tpm)
    colnames(tpm_qn) = colnames(tpm)
    rownames(tpm_qn) = rownames(tpm)
    pcxx = prcomp(t(tpm_qn),scale = T)
    pcx = data.frame(pcxx$x)
    pcx$sample_id = rownames(pcx)
    pcx = left_join(pcx,att) %>%
      mutate(death = as.factor(death))
    
    varexplained = 100*summary(pcxx)$imp[2,]
    p1 = pcx %>%
      ggplot(aes(x = PC1, y = PC2, color = age, shape = sex)) +
      geom_point(size = 3) +
      scale_color_brewer(type='seq',palette = 1) +
      scale_shape_manual(values =c(19,15)) +
      theme_pubr(legend = 'right') +
      ggtitle(tisname) +
      xlab(paste('PC1 (',round(varexplained['PC1'],2),'%)',sep=''))+
      ylab(paste('PC2 (',round(varexplained['PC2'],2),'%)',sep='')) +
      coord_fixed()
    p2 = pcx %>%
      ggplot(aes(x = PC3, y = PC4, color = age, shape = sex)) +
      geom_point(size = 3) +
      scale_color_brewer(type='seq',palette = 1) +
      scale_shape_manual(values = c(19,15)) +
      theme_pubr(legend = 'right') +
      ggtitle(tisname) +
      xlab(paste('PC3 (',round(varexplained['PC3'],2),'%)',sep=''))+
      ylab(paste('PC4 (',round(varexplained['PC4'],2),'%)',sep='')) +
      coord_fixed()
    
    p3 = pcx %>%
      ggplot(aes(x = PC5, y = PC6, color = age, shape = sex)) +
      geom_point(size = 3) +
      scale_color_brewer(type='seq',palette = 1) +
      scale_shape_manual(values = c(19,15)) +
      theme_pubr(legend = 'right') +
      ggtitle(tisname) +
      xlab(paste('PC5 (',round(varexplained['PC5'],2),'%)',sep=''))+
      ylab(paste('PC6 (',round(varexplained['PC6'],2),'%)',sep='')) +
      coord_fixed()
    
    p4 = pcx %>%
      select(PC1,PC2,PC3,PC4,PC5,PC6,age) %>%
      gather(key = 'PC', value = 'value',-age) %>%
      ggplot(aes(x = age, y = value, fill = age)) +
      scale_fill_brewer(type='seq',palette = 1) +
      geom_boxplot() +
      facet_wrap(~PC,scales = 'free_y') +
      theme_pubr(x.text.angle = 90)
    
    system(paste('mkdir -p results/qc/',tisname,sep=''))
    ggsave(paste('results/qc/',tisname,'/',round,'pc1-pc2.pdf',sep=''),p1,units = 'cm',width = 18,height = 18,useDingbats = F)
    ggsave(paste('results/qc/',tisname,'/',round,'pc1-pc2.png',sep=''),p1,units = 'cm',width = 18,height = 18)
    
    ggsave(paste('results/qc/',tisname,'/',round,'pc3-pc4.pdf',sep=''),p2,units = 'cm',width = 18,height = 18,useDingbats = F)
    ggsave(paste('results/qc/',tisname,'/',round,'pc3-pc4.png',sep=''),p2,units = 'cm',width = 18,height = 18)
    
    
    ggsave(paste('results/qc/',tisname,'/',round,'pc5-pc6.pdf',sep=''),p3,units = 'cm',width = 18,height = 18,useDingbats = F)
    ggsave(paste('results/qc/',tisname,'/',round,'pc5-pc6.png',sep=''),p3,units = 'cm',width = 18,height = 18)
    
    
    ggsave(paste('results/qc/',tisname,'/',round,'pc-age.pdf',sep=''),p4,units = 'cm',width = 18,height = 18,useDingbats = F)
    ggsave(paste('results/qc/',tisname,'/',round,'pc-age.png',sep=''),p4,units = 'cm',width = 18,height = 18)
    
    pcsdx = apply(select(pcx,PC1,PC2,PC3,PC4),2,sd)
    x = t(apply(select(pcx,PC1,PC2,PC3,PC4),1,function(x) any(abs(x)>=(3*pcsdx))))
    if(any(x)) {
      round = round + 1
      outx =c(outliers,pcx$sample_id[which(t(apply(select(pcx,PC1,PC2,PC3,PC4),1,function(x) any(abs(x)>=(3*pcsdx)))))])
      qc_check(mytis,outliers = outx , round = round)
    } else { 
      saveRDS(tpm_qn,paste('./data/processed/expression/tpm_lm_l2_qn/',tisname,'.rds',sep='')) 
    }
  }
}
system('mkdir -p ./data/processed/expression/tpm_lm_l2_qn')
alltis = list.files('data/processed/expression/tpm/')
sapply(alltis,qc_check)
