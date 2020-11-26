library(tidyverse)
library(ggpubr)
library(ggthemes)
alltis = list.files('./data/processed/expression/tpm_lm_l2_qn/')
att = readRDS('./data/processed/attr.rds')

for(mytis in alltis){
  tisname = gsub('.rds','',mytis)
  expr = readRDS(paste('./data/processed/expression/tpm_lm_l2_qn/',mytis,sep=''))
  # if(ncol(expr)>=20){
  # if(!file.exists(paste('./data/processed/expression/change_w_age/',mytis,sep=''))){
    print(mytis)
    age = setNames(att$age,att$sample_id)[colnames(expr)]
    age = c(seq(20,70,by=10) + 4.5)[age]
    ch = t(apply(expr,1,function(x){
      co = cor.test(x,age,method='s')
      c(co$est,co$p.val)
    }))
    ch = as.data.frame(ch)
    colnames(ch) = c('rho','pval')
    ch$padj = p.adjust(ch$pval,method = 'fdr')
    # system('mkdir -p ./data/processed/expression/change_w_age/')
    saveRDS(ch, paste('./data/processed/expression/change_w_age/',mytis,sep=''))
    ch$gene = rownames(ch)
    ch$sign = c('negative','positive')[1+(ch$rho>=0)]
    ch$signif = ch$padj<=0.1
    pp = ggplot(ch,aes(x=sign,fill=signif)) +
      geom_bar(position = 'stack') +
      scale_y_log10() +
      scale_fill_few(guide = guide_legend('Significant\n(q<=0.1)')) +
      xlab('') + ylab('') +
      ggtitle(tisname)+
      theme_pubr(legend = 'right', x.text.angle = 90) 
    system(paste('mkdir -p ./results/change_with_age/',tisname,sep=''))
    ggsave(paste('./results/change_with_age/',tisname,'/change_in_numbers.pdf',sep=''),pp,units = 'cm',width = 8, height = 8, useDingbats = F)
    ggsave(paste('./results/change_with_age/',tisname,'/change_in_numbers.png',sep=''),pp,units = 'cm',width = 8, height = 8)
    expr_sc = t(apply(expr,1,scale))
    colnames(expr_sc) = colnames(expr)
    # system('mkdir -p ./data/processed/expression/tpm_lm_l2_qn_sc/')
    saveRDS(expr_sc, paste('./data/processed/expression/tpm_lm_l2_qn_sc/',mytis,sep=''))
    clx = kmeans(expr_sc,16)$clust
    clx = reshape2::melt(clx) %>%
      rename(cluster = value) %>%
      mutate(geneid = names(clx))
    allcluster = as.data.frame(expr_sc) %>%
      mutate(geneid = rownames(expr_sc)) %>%
      gather(key = 'sample_id', value = 'expression', -geneid) %>%
      left_join(clx) %>%
      group_by(sample_id,cluster) %>%
      summarise(meanx = mean(expression),
                medianx = median(expression),
                sdx = sd(expression),
                fq = quantile(expression,probs = 0.25),
                tq = quantile(expression,probs = 0.75),
                clsize = length(unique(geneid))) %>%
      left_join(att) %>%
      mutate(age = c(seq(20,70,by=10) + 4.5)[age]) %>%
      mutate(clname = paste('cluster',cluster,' #=',clsize,sep=''))
    p1 = ggplot(allcluster,aes(x = age)) +
      geom_point(aes(y=medianx),width = 0.3)+
      geom_smooth(aes(y= medianx), color = 'darkred', fill ='pink') +
      facet_wrap(~clname, scales = 'free_y')+
      theme_pubr()
    ggsave(paste('./results/change_with_age/',tisname,'/allcluster16.pdf',sep=''),p1,units = 'cm',width = 18, height = 18, useDingbats = F)
    ggsave(paste('./results/change_with_age/',tisname,'/allcluster16.png',sep=''),p1,units = 'cm',width = 18, height = 18)
    signifgenes = rownames(ch)[which(ch$padj<=0.1)]
    if(length(signifgenes)>=6){
      clx = kmeans(expr_sc[signifgenes,],6)$clust
      clx = reshape2::melt(clx) %>%
        rename(cluster = value) %>%
        mutate(geneid = names(clx))
      allcluster = as.data.frame(expr_sc[signifgenes,]) %>%
        mutate(geneid = rownames(expr_sc[signifgenes,])) %>%
        gather(key = 'sample_id', value = 'expression', -geneid) %>%
        left_join(clx) %>%
        group_by(sample_id,cluster) %>%
        summarise(meanx = mean(expression),
                  medianx = median(expression),
                  sdx = sd(expression),
                  fq = quantile(expression,probs = 0.25),
                  tq = quantile(expression,probs = 0.75),
                  clsize = length(unique(geneid))) %>%
        left_join(att) %>%
        mutate(age = c(seq(20,70,by=10) + 4.5)[age]) %>%
        mutate(clname = paste('cluster',cluster,' #=',clsize,sep=''))
      p2 = ggplot(allcluster,aes(x = age)) +
        geom_jitter(aes(y=medianx),width = 0.3)+
        geom_smooth(aes(y= medianx), color = 'darkred', fill ='pink') +
        facet_wrap(~clname, scales = 'free_y')+
        theme_pubr()
      ggsave(paste('./results/change_with_age/',tisname,'/signifcluster6.pdf',sep=''),p2,units = 'cm',width = 18, height = 12, useDingbats = F)
      ggsave(paste('./results/change_with_age/',tisname,'/signifcluster6.png',sep=''),p2,units = 'cm',width = 18, height = 12)
    }
    signifgenes = names(sort(abs(setNames(ch$rho,rownames(ch))),dec=T)[1:2000])
    clx = kmeans(expr_sc[signifgenes,],8)$clust
    clx = reshape2::melt(clx) %>%
      rename(cluster = value) %>%
      mutate(geneid = names(clx))
    allcluster = as.data.frame(expr_sc[signifgenes,]) %>%
      mutate(geneid = rownames(expr_sc[signifgenes,])) %>%
      gather(key = 'sample_id', value = 'expression', -geneid) %>%
      left_join(clx) %>%
      group_by(sample_id,cluster) %>%
      summarise(meanx = mean(expression),
                medianx = median(expression),
                sdx = sd(expression),
                fq = quantile(expression,probs = 0.25),
                tq = quantile(expression,probs = 0.75),
                clsize = length(unique(geneid))) %>%
      left_join(att) %>%
      mutate(age = c(seq(20,70,by=10) + 4.5)[age]) %>%
      mutate(clname = paste('cluster',cluster,' #=',clsize,sep=''))
    p3 = ggplot(allcluster,aes(x = age)) +
      geom_jitter(aes(y=medianx),width = 0.3)+
      geom_smooth(aes(y= medianx), color = 'darkred', fill ='pink') +
      facet_wrap(~clname, scales = 'free_y')+
      theme_pubr()
    ggsave(paste('./results/change_with_age/',tisname,'/top2000cluster8.pdf',sep=''),p3,units = 'cm',width = 18, height = 15, useDingbats = F)
    ggsave(paste('./results/change_with_age/',tisname,'/top2000cluster8.png',sep=''),p3,units = 'cm',width = 18, height = 15)
  # }
  rm(list=setdiff(ls(),c('alltis','att')))
}
