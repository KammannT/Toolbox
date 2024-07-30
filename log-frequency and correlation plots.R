# Script for analysis of frequency and correlation
# Setup
required_packages <- c("readxl", "reshape2", "see", "purrr", "ggprism", "rstatix", "scales", "stats", "pheatmap", "RColorBrewer", "Hmisc", "corrplot", "vtable", "ggpubr", "tidyverse")
for(package in required_packages){
  if(!require(package, character.only=TRUE)){install.packages(package, dependencies=TRUE)
    library(package, character.only = TRUE) } }

script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_dir)


# Frequency analysis (on log-scale plot)
# dummy data
donors <- paste0('Donor', 1:10)
variable2 <- c('Celltype1', 'Celltype2') # could be extended for more complex analysis

data <- data.frame(
  Donor = rep(donors, times=length(variable2), replace=FALSE),
  variable2 = rep(variable2, each=length(donors)),
  Frequency = c(rnorm(10, mean=2, sd=0.8),
                rnorm(10, mean=12, sd=2)) )
print(data)

# frequency plot
data_stats <- data %>% wilcox_test(Frequency ~ variable2, p.adjust.method='none') %>% add_significance(p.col='p', cutpoints=c(0,0.001,0.01,0.05,1), symbols=c('***','**','*','ns')) %>% add_xy_position(x='variable2')
data_mean  <- data %>% group_by(variable2) %>% summarise(Mean_Frequency=mean(Frequency)) # mean calculation before log transformation, otherwise wrong on ggplot

ggplot(data, aes(x=variable2, y=Frequency))+
  geom_boxplot(aes(fill=variable2), color='black', show.legend=FALSE, outlier.shape=NA)+
  scale_fill_manual(values=c('lightgrey', 'lightblue'))+
  geom_point(position=position_jitter(width=0.1, height=0), shape=21, size=2, color='black', fill='white', show.legend=FALSE)+
  stat_summary(data=data_mean, aes(x=variable2, y=Mean_Frequency), fun='mean', geom='point', shape=15, size=1, color='black')+
  scale_y_continuous(limits=c(0.1,20), breaks=c(0.1,1,10,20), labels=scales::label_number(big.mark=',', accuracy=0.01, drop0trailing=TRUE), expand=c(0,0), oob=squish)+
  coord_trans(y='log10', clip='off')+
  annotation_logticks(base=10, sides='l', scaled=FALSE, outside=TRUE, size=0.35, long=unit(0.15, 'cm'), mid=unit(0.1, 'cm'), short=unit(0.1, 'cm'))+
  labs(x='', y='Frequency\nof celltype (%)')+
  add_pvalue(data_stats, tip.length=0, bracket.size=0.35, size=5)+
  theme_classic()+
  theme(axis.text.x=element_text(angle=45, hjust=1))

# Correlation "heatmap"-like matrix
# adopted from hmisc and corrplot
donor_age <- data.frame(Age=round(rnorm(mean=40, sd=15, length(donors)), digits=0), Donor=donors)
data_correlation <- data %>% left_join(., donor_age) %>% select(Donor, Age, variable2, Frequency) %>% pivot_wider(names_from='variable2', names_prefix='', values_from='Frequency') %>% select(-Donor)

t1cor      <- rcorr(as.matrix(data_correlation), type='spearman') #rcorr from Hmisc package
print(as.data.frame(t1cor[['P']]))
print(as.data.frame(t1cor[['r']]))

corrplot(t1cor$r, mar=c(0,0,0,0), type="lower", diag=FALSE, tl.col='black', #corrplot from corrplot package
         tl.srt=45, tl.cex=1, pch.cex=0.5, cl.cex=0.5, method='square',
         col=rev(brewer.pal(n=8, name="RdBu")),
         p.mat=t1cor$P, insig='p-value', sig.level=-1,
         na.label=' ', na.label.col='white')

# 1vs1 correlation graph
ggplot(data_correlation, aes(x=Celltype1, y=Celltype2)) +
  geom_smooth(method="lm", se=FALSE, fullrange=TRUE, show.legend=FALSE, color='darkgrey') +
  geom_point(size=2, show.legend=FALSE, shape=21, color='black')+
  scale_x_continuous(limits=c(0,17.6), breaks=c(0,4,8,12,16), expand=c(0, 0.02), guide='prism_offset')+
  scale_y_continuous(limits=c(0,17.6), breaks=c(0,4,8,12,16), expand=c(0, 0.02), guide='prism_offset')+
  coord_cartesian()+
  theme_classic()+
  theme(axis.text=element_text(size=8))



