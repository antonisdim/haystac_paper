library(ggplot2)
library(reshape2)
library(wesanderson)
library(ggforce)
library(dplyr)
library(scales)
library(gridExtra)

pal_final <- c('#d7191c', '#fdae61', '#abd9e9', '#2c7bb6')

lm_eqn <- function(df){
  m <- lm(Mean_Posterior_Abundance_Genus ~ Mean_Posterior_Abundance_RefSeq, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

lmp <- function (df) {
  modelobject <- lm(Mean_Posterior_Abundance_Genus ~ Mean_Posterior_Abundance_RefSeq, df);
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


theme_Publication_simple <- function(base_size=14, base_family="Helvetica", xangle=0, 
                                     position="bottom", direction="horizontal", 
                                     legend.title=element_text(face="italic"),
                                     legend.text=element_text()) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.text.x = element_text(angle=xangle,hjust=0.5),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = position,
            legend.direction = direction,
            legend.key.size= unit(1, "cm"),
            legend.key.width = unit(1.5, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = legend.title,
            legend.text = legend.text,
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

##### ORAL MICRO DATASET #####

summary_stats1 <- read.csv(paste("IrinaSims/", "IrinaSims_FP_FN_Rates_All_001_new.tsv", sep=''), header=TRUE, sep='\t')

summary_stats1$Status <- ifelse(grepl("anc", summary_stats1$Sample), 
                                "Oral microbiome aDNA damage", "Oral microbiome modern")

summary_stats1$label <- ifelse(grepl("anc", summary_stats1$Sample), 
                               "B", "A")

sum_stat_counts1 <- summary_stats1[, c('Sample', 'False_Positive', 'False_Negative', 'True_Detected_Species', 'Method', 'Status', 'label')]

sum_stas_counts_melt1 <- melt(sum_stat_counts1, variable.names = "variable", value.names = "value", 
                              id.vars = c("Sample", "Method", "Status", 'label'))

avg_counts_per_sample_summary_1 <- sum_stas_counts_melt1 %>% group_by(Method, Status, label, variable) %>% 
  summarise_all(funs(mean, sd))

levels(avg_counts_per_sample_summary_1$Method) <- c(levels(avg_counts_per_sample_summary_1$Method), "Kraken2/Bracken")
avg_counts_per_sample_summary_1$Method[avg_counts_per_sample_summary_1$Method == 'kraken'] <- 'Kraken2/Bracken'

levels(avg_counts_per_sample_summary_1$Method) <- c(levels(avg_counts_per_sample_summary_1$Method), "Krakenuniq")
avg_counts_per_sample_summary_1$Method[avg_counts_per_sample_summary_1$Method == 'krakenuniq'] <- 'Krakenuniq'

levels(avg_counts_per_sample_summary_1$Method) <- c(levels(avg_counts_per_sample_summary_1$Method), "MALT")
avg_counts_per_sample_summary_1$Method[avg_counts_per_sample_summary_1$Method == 'malt'] <- 'MALT'

levels(avg_counts_per_sample_summary_1$Method) <- c(levels(avg_counts_per_sample_summary_1$Method), "HAYSTAC")
avg_counts_per_sample_summary_1$Method[avg_counts_per_sample_summary_1$Method == 'haystac'] <- 'HAYSTAC'

levels(avg_counts_per_sample_summary_1$Method) <- c(levels(avg_counts_per_sample_summary_1$Method), "Sigma")
avg_counts_per_sample_summary_1$Method[avg_counts_per_sample_summary_1$Method == 'sigma'] <- 'Sigma'

avg_counts_per_sample_summary_1$Method <- factor(avg_counts_per_sample_summary_1$Method, levels = c('HAYSTAC', 'Sigma', 'Kraken2/Bracken', 'Krakenuniq','MALT'))

avg_counts_per_sample_summary_1$Status <- factor(avg_counts_per_sample_summary_1$Status, levels = c('Oral microbiome modern', 'Oral microbiome aDNA damage'))

# png(paste("Oral_Micro_Counts_Per_Sample.png", sep=''), height = 650, width = 1300)

ggsave("Oral_Micro_Counts_Per_Sample.pdf", height = 10, width = 20, units = "in")

ggplot(avg_counts_per_sample_summary_1, aes(fill=variable, y=value_mean, x=Method)) + 
  geom_bar(position="dodge", stat="identity") + 
  geom_text(aes(label=format(round(value_mean, 1),nsmall=1), 
                y=value_mean + avg_counts_per_sample_summary_1$value_sd/sqrt(nrow(avg_counts_per_sample_summary_1))), 
            position=position_dodge(width=0.9), vjust=-0.5, size=7) +
  geom_errorbar(aes(ymin=value_mean - avg_counts_per_sample_summary_1$value_sd/sqrt(nrow(avg_counts_per_sample_summary_1)), 
                    ymax=value_mean + avg_counts_per_sample_summary_1$value_sd/sqrt(nrow(avg_counts_per_sample_summary_1))), 
                width=.2, position=position_dodge(.9)) +
  scale_fill_manual(values=pal_final, name = "", labels = c("False positive", "False negative", "True positive")) +
  facet_wrap(~Status, ncol=2) +
  geom_hline(data=filter(avg_counts_per_sample_summary_1, Status=="Oral microbiome aDNA damage"), aes(yintercept=178), linetype="dashed") +
  geom_hline(data=filter(avg_counts_per_sample_summary_1, Status=="Oral microbiome modern"), aes(yintercept=178), linetype="dashed") +
  geom_hline(data=filter(avg_counts_per_sample_summary_1, Status=="Oral microbiome aDNA damage"), aes(yintercept=190), alpha=0, linetype="dashed") +
  geom_hline(data=filter(avg_counts_per_sample_summary_1, Status=="Oral microbiome modern"), aes(yintercept=350), alpha=0, linetype="dashed") +
  geom_vline(xintercept=seq(1.5, length(unique(avg_counts_per_sample_summary_1$Method))-0.5, 1), 
             lwd=1, colour="#f0f0f0") +
  theme(text = element_text(size=20, face="bold"), axis.title=element_text(size=22, face="bold")) +
  ylab("Mean species count") + xlab('') + 
  theme_Publication_simple(base_size = 22, base_family = "Helvetica") +
  geom_text(aes(x = -Inf, y = Inf, label = label, group = label), size= theme_get()$text[["size"]], hjust = -0.5, vjust = 1.5)

dev.off()


sum_stats_percent1 <- summary_stats1[, c('Sample','False_Positive_Rate', 
                                         'False_Negative_Rate', 'Method', 'Status', 'label')]

sum_stats_percent_melt_1 <- melt(sum_stats_percent1, variable.names = "variable", value.names = "value", 
                                 id.vars = c("Sample", "Method", "Status", "label"))

avg_rates_per_sample_summary_1 <- sum_stats_percent_melt_1 %>% group_by(Method, Status, label, variable) %>% 
  summarise_all(funs(mean, sd))

levels(avg_rates_per_sample_summary_1$Method) <- c(levels(avg_rates_per_sample_summary_1$Method), "Kraken2/Bracken")
avg_rates_per_sample_summary_1$Method[avg_rates_per_sample_summary_1$Method == 'kraken'] <- 'Kraken2/Bracken'

levels(avg_rates_per_sample_summary_1$Method) <- c(levels(avg_rates_per_sample_summary_1$Method), "Krakenuniq")
avg_rates_per_sample_summary_1$Method[avg_rates_per_sample_summary_1$Method == 'krakenuniq'] <- 'Krakenuniq'

levels(avg_rates_per_sample_summary_1$Method) <- c(levels(avg_rates_per_sample_summary_1$Method), "MALT")
avg_rates_per_sample_summary_1$Method[avg_rates_per_sample_summary_1$Method == 'malt'] <- 'MALT'

levels(avg_rates_per_sample_summary_1$Method) <- c(levels(avg_rates_per_sample_summary_1$Method), "HAYSTAC")
avg_rates_per_sample_summary_1$Method[avg_rates_per_sample_summary_1$Method == 'haystac'] <- 'HAYSTAC'

levels(avg_rates_per_sample_summary_1$Method) <- c(levels(avg_rates_per_sample_summary_1$Method), "Sigma")
avg_rates_per_sample_summary_1$Method[avg_rates_per_sample_summary_1$Method == 'sigma'] <- 'Sigma'

avg_rates_per_sample_summary_1$Method <- factor(avg_rates_per_sample_summary_1$Method, levels = c('HAYSTAC', 'Sigma', 'Kraken2/Bracken', 'Krakenuniq', 'MALT'))
avg_rates_per_sample_summary_1$Status <- factor(avg_rates_per_sample_summary_1$Status, levels = c('Oral microbiome modern', 'Oral microbiome aDNA damage'))

# png(paste("Oral_Micro_Per_Sample_Rates.png", sep=''), height = 650, width = 1300)

ggsave("Oral_Micro_Per_Sample_Rates.pdf", height = 10, width = 20, units = "in")

ggplot(avg_rates_per_sample_summary_1, aes(fill=variable, y=value_mean, x=Method)) + 
  geom_bar(position="dodge", stat="identity") + 
  geom_text(aes(label=format(round(value_mean, 1),nsmall=1), 
                y=value_mean + avg_rates_per_sample_summary_1$value_sd/sqrt(nrow(avg_rates_per_sample_summary_1))), 
            position=position_dodge(width=0.9), vjust=-0.5, size=7) +
  geom_errorbar(aes(ymin=value_mean - avg_rates_per_sample_summary_1$value_sd/sqrt(nrow(avg_rates_per_sample_summary_1)), 
                    ymax=value_mean + avg_rates_per_sample_summary_1$value_sd/sqrt(nrow(avg_rates_per_sample_summary_1))), 
                width=.2, position=position_dodge(.9)) +
  scale_fill_manual(values=pal_final, name = "", labels = c("False positive", "False negative")) +
  facet_wrap(~Status, ncol=2) +
  geom_hline(data=filter(avg_rates_per_sample_summary_1, Status=="Oral microbiome aDNA damage"), aes(yintercept=50), linetype="dashed", alpha=0) +
  geom_hline(data=filter(avg_rates_per_sample_summary_1, Status=="Oral microbiome modern"), aes(yintercept=50), linetype="dashed", alpha=0) +
  theme(text = element_text(size=20, face="bold"), axis.title=element_text(size=22, face="bold")) +
  ylab("Rate %") + xlab('') +
  geom_vline(xintercept=seq(1.5, length(unique(avg_rates_per_sample_summary_1$Method))-0.5, 1), 
             lwd=1, colour="#f0f0f0") +
  theme_Publication_simple(base_size = 25, base_family = "Helvetica") +
  geom_text(aes(x = -Inf, y = Inf, label = label, group = label), size= theme_get()$text[["size"]], hjust = -0.5, vjust = 1.5)

dev.off()



##### GENERAL MICRO DATASET #####

summary_stats2 <- read.csv(paste("NewSims/", "NewSims_FP_FN_Rates_All_001_new.tsv", sep=''), header=TRUE, sep='\t')

summary_stats2$Status <- ifelse(grepl("100_anc", summary_stats2$Sample), 
                                "Microbiome 100 species ancient", "")
summary_stats2$Status <- ifelse(grepl("100_mod", summary_stats2$Sample), 
                                "Microbiome 100 species modern", summary_stats2$Status)
summary_stats2$Status <- ifelse(grepl("500_anc", summary_stats2$Sample), 
                                "Microbiome 500 species ancient", summary_stats2$Status)
summary_stats2$Status <- ifelse(grepl("500_mod", summary_stats2$Sample), 
                                "Microbiome 500 species modern", summary_stats2$Status)

summary_stats2$label <- ifelse(grepl("100_anc", summary_stats2$Sample), 
                               "A", "")
summary_stats2$label <- ifelse(grepl("100_mod", summary_stats2$Sample), 
                               "B", summary_stats2$label)
summary_stats2$label <- ifelse(grepl("500_anc", summary_stats2$Sample), 
                               "C", summary_stats2$label)
summary_stats2$label <- ifelse(grepl("500_mod", summary_stats2$Sample), 
                               "D", summary_stats2$label)

sum_stat_counts2 <- summary_stats2[, c('Sample', 'False_Positive', 'False_Negative', 'True_Detected_Species', 'Method', 'Status', 'label')]

sum_stas_counts_melt2 <- melt(sum_stat_counts2, variable.names = "variable", value.names = "value", 
                              id.vars = c("Sample", "Method", "Status", 'label'))

avg_counts_per_sample_summary_2 <- sum_stas_counts_melt2 %>% group_by(Method, Status, variable, label) %>% 
  summarise_all(funs(mean, sd))

levels(avg_counts_per_sample_summary_2$Method) <- c(levels(avg_counts_per_sample_summary_2$Method), "Kraken2/Bracken")
avg_counts_per_sample_summary_2$Method[avg_counts_per_sample_summary_2$Method == 'kraken'] <- 'Kraken2/Bracken'

levels(avg_counts_per_sample_summary_2$Method) <- c(levels(avg_counts_per_sample_summary_2$Method), "Krakenuniq")
avg_counts_per_sample_summary_2$Method[avg_counts_per_sample_summary_2$Method == 'krakenuniq'] <- 'Krakenuniq'

levels(avg_counts_per_sample_summary_2$Method) <- c(levels(avg_counts_per_sample_summary_2$Method), "MALT")
avg_counts_per_sample_summary_2$Method[avg_counts_per_sample_summary_2$Method == 'malt'] <- 'MALT'

levels(avg_counts_per_sample_summary_2$Method) <- c(levels(avg_counts_per_sample_summary_2$Method), "HAYSTAC")
avg_counts_per_sample_summary_2$Method[avg_counts_per_sample_summary_2$Method == 'haystac'] <- 'HAYSTAC'

levels(avg_counts_per_sample_summary_2$Method) <- c(levels(avg_counts_per_sample_summary_2$Method), "Sigma")
avg_counts_per_sample_summary_2$Method[avg_counts_per_sample_summary_2$Method == 'sigma'] <- 'Sigma'


avg_counts_per_sample_summary_2$Method <- factor(avg_counts_per_sample_summary_2$Method, levels = c('HAYSTAC', 'Sigma', 'Kraken2/Bracken', 'Krakenuniq', 'MALT'))

# png(paste("General_Micro_Counts_Per_Sample.png", sep=''), height = 800, width = 1400)

ggsave("General_Micro_Counts_Per_Sample.pdf", height = 12, width = 22, units = "in")

ggplot(avg_counts_per_sample_summary_2, aes(fill=variable, y=value_mean, x=Method)) + 
  geom_bar(position="dodge", stat="identity") + 
  geom_text(aes(label=format(round(value_mean, 1),nsmall=1), 
                y=value_mean + avg_counts_per_sample_summary_2$value_sd/sqrt(nrow(avg_counts_per_sample_summary_2))), 
            position=position_dodge(width=0.9), vjust=-0.5, size=6) +
  geom_errorbar(aes(ymin=value_mean - avg_counts_per_sample_summary_2$value_sd/sqrt(nrow(avg_counts_per_sample_summary_2)), 
                    ymax=value_mean + avg_counts_per_sample_summary_2$value_sd/sqrt(nrow(avg_counts_per_sample_summary_2))), 
                width=.2, position=position_dodge(.9)) +
  scale_fill_manual(values=pal_final, name = "", labels = c("False positive", "False negative", "True positive")) +
  facet_wrap(~Status, ncol=2, scales="free") +
  geom_hline(data=filter(avg_counts_per_sample_summary_2, Status=="Microbiome 100 species ancient"), aes(yintercept=100), linetype="dashed") +
  geom_hline(data=filter(avg_counts_per_sample_summary_2, Status=="Microbiome 100 species modern"), aes(yintercept=100), linetype="dashed") +
  geom_hline(data=filter(avg_counts_per_sample_summary_2, Status=="Microbiome 500 species ancient"), aes(yintercept=500), linetype="dashed") +
  geom_hline(data=filter(avg_counts_per_sample_summary_2, Status=="Microbiome 500 species modern"), aes(yintercept=500), linetype="dashed") +
  geom_hline(data=filter(avg_counts_per_sample_summary_2, Status=="Microbiome 100 species ancient"), aes(yintercept=105), linetype="dashed", alpha=0) +
  geom_hline(data=filter(avg_counts_per_sample_summary_2, Status=="Microbiome 100 species modern"), aes(yintercept=105), linetype="dashed", alpha=0) +
  geom_hline(data=filter(avg_counts_per_sample_summary_2, Status=="Microbiome 500 species ancient"), aes(yintercept=510), linetype="dashed", alpha=0) +
  geom_hline(data=filter(avg_counts_per_sample_summary_2, Status=="Microbiome 500 species modern"), aes(yintercept=510), linetype="dashed", alpha=0) +
  geom_vline(xintercept=seq(1.5, length(unique(avg_counts_per_sample_summary_2$Method))-0.5, 1), 
             lwd=1, colour="#f0f0f0") +
  theme(text = element_text(size=20, face="bold"), axis.title=element_text(size=22, face="bold")) +
  ylab("Mean species count") + xlab('') + 
  theme_Publication_simple(base_size = 25, base_family = "Helvetica") +
  geom_text(aes(x = -Inf, y = Inf, label = label, group = label), size= theme_get()$text[["size"]], hjust = -0.5, vjust = 1.5)

dev.off()


sum_stats_percent2 <- summary_stats2[, c('Sample','False_Positive_Rate', 
                                         'False_Negative_Rate', 'Method', 'Status', 'label')]

sum_stats_percent_melt_2 <- melt(sum_stats_percent2, variable.names = "variable", value.names = "value", 
                                 id.vars = c("Sample", "Method", "Status", "label"))

avg_rates_per_sample_summary_2 <- sum_stats_percent_melt_2 %>% group_by(Method, Status, variable, label) %>% 
  summarise_all(funs(mean, sd))

levels(avg_rates_per_sample_summary_2$Method) <- c(levels(avg_rates_per_sample_summary_2$Method), "Kraken2/Bracken")
avg_rates_per_sample_summary_2$Method[avg_rates_per_sample_summary_2$Method == 'kraken'] <- 'Kraken2/Bracken'

levels(avg_rates_per_sample_summary_2$Method) <- c(levels(avg_rates_per_sample_summary_2$Method), "Krakenuniq")
avg_rates_per_sample_summary_2$Method[avg_rates_per_sample_summary_2$Method == 'krakenuniq'] <- 'Krakenuniq'

levels(avg_rates_per_sample_summary_2$Method) <- c(levels(avg_rates_per_sample_summary_2$Method), "MALT")
avg_rates_per_sample_summary_2$Method[avg_rates_per_sample_summary_2$Method == 'malt'] <- 'MALT'

levels(avg_rates_per_sample_summary_2$Method) <- c(levels(avg_rates_per_sample_summary_2$Method), "HAYSTAC")
avg_rates_per_sample_summary_2$Method[avg_rates_per_sample_summary_2$Method == 'haystac'] <- 'HAYSTAC'

levels(avg_rates_per_sample_summary_2$Method) <- c(levels(avg_rates_per_sample_summary_2$Method), "Sigma")
avg_rates_per_sample_summary_2$Method[avg_rates_per_sample_summary_2$Method == 'sigma'] <- 'Sigma'

avg_rates_per_sample_summary_2$Method <- factor(avg_rates_per_sample_summary_2$Method, levels = c('HAYSTAC', 'Sigma', 'Kraken2/Bracken', 'Krakenuniq', 'MALT'))

#png(paste("General_Micro_Per_Sample_Rates.png", sep=''), height = 600, width = 1200)

ggsave("General_Micro_Per_Sample_Rates.pdf", height = 18, width = 22, units = "in")

ggplot(avg_rates_per_sample_summary_2, aes(fill=variable, y=value_mean, x=Method)) + 
  geom_bar(position="dodge", stat="identity") + 
  geom_text(aes(label=format(round(value_mean, 1),nsmall=1), 
                y=value_mean + avg_rates_per_sample_summary_2$value_sd/sqrt(nrow(avg_rates_per_sample_summary_2))), 
            position=position_dodge(width=0.9), vjust=-0.5, size=7) +
  geom_errorbar(aes(ymin=value_mean - avg_rates_per_sample_summary_2$value_sd/sqrt(nrow(avg_rates_per_sample_summary_2)), 
                    ymax=value_mean + avg_rates_per_sample_summary_2$value_sd/sqrt(nrow(avg_rates_per_sample_summary_2))), 
                width=.2, position=position_dodge(.9)) +
  scale_fill_manual(values=pal_final, name = "", labels = c("False positive", "False negative")) +
  facet_wrap(~Status, ncol=2, scales="free_y") +
  geom_hline(data=filter(avg_rates_per_sample_summary_2, Status=="Microbiome 100 species ancient"), aes(yintercept=10), linetype="dashed", alpha=0) +
  geom_hline(data=filter(avg_rates_per_sample_summary_2, Status=="Microbiome 100 species modern"), aes(yintercept=10), linetype="dashed", alpha=0) +
  geom_hline(data=filter(avg_rates_per_sample_summary_2, Status=="Microbiome 500 species ancient"), aes(yintercept=10), linetype="dashed", alpha=0) +
  geom_hline(data=filter(avg_rates_per_sample_summary_2, Status=="Microbiome 500 species modern"), aes(yintercept=10), linetype="dashed", alpha=0) +
  theme(text = element_text(size=20, face="bold"), axis.title=element_text(size=22, face="bold")) +
  ylab("Rate %") + xlab('') +
  geom_vline(xintercept=seq(1.5, length(unique(avg_rates_per_sample_summary_2$Method))-0.5, 1), 
             lwd=1, colour="#f0f0f0") +
  theme_Publication_simple(base_size = 25, base_family = "Helvetica") +
  geom_text(aes(x = -Inf, y = Inf, label = label, group = label), size= theme_get()$text[["size"]], hjust = -0.5, vjust = 1.5)

dev.off()



##### SIMPLE SIMS DATASET #####

summary_stats3 <- read.csv(paste("SimpleSims/", "SimpleSims_FP_FN_Rates_All_001_new.tsv", sep=''), header=TRUE, sep='\t')

summary_stats3$Status <- ifelse(grepl("nocont", summary_stats3$Sample), 
                                "Non Contaminated", "Contaminated")

summary_stats3$label <- ifelse(grepl("nocont", summary_stats3$Sample), 
                                "A", "B")

sum_stats_percent3 <- summary_stats3[, c('Sample','False_Positive', 
                                         'False_Negative', 'True_Detected_Species', 'Method', 'Status', 'label')]

sum_stats_melt3 <- melt(sum_stats_percent3, variable.names = "variable", value.names = "value", 
                        id.vars = c("Sample", "Method", "Status", "label"))

avg_rates_per_cont <- sum_stats_melt3 %>% group_by(Method, Status, variable, label) %>% 
  summarise_all(funs(mean, sd))

levels(avg_rates_per_cont$Method) <- c(levels(avg_rates_per_cont$Method), "Kraken2/Bracken")
avg_rates_per_cont$Method[avg_rates_per_cont$Method == 'kraken'] <- 'Kraken2/Bracken'

levels(avg_rates_per_cont$Method) <- c(levels(avg_rates_per_cont$Method), "Krakenuniq")
avg_rates_per_cont$Method[avg_rates_per_cont$Method == 'krakenuniq'] <- 'Krakenuniq'

levels(avg_rates_per_cont$Method) <- c(levels(avg_rates_per_cont$Method), "MALT")
avg_rates_per_cont$Method[avg_rates_per_cont$Method == 'malt'] <- 'MALT'

levels(avg_rates_per_cont$Method) <- c(levels(avg_rates_per_cont$Method), "HAYSTAC")
avg_rates_per_cont$Method[avg_rates_per_cont$Method == 'haystac'] <- 'HAYSTAC'

levels(avg_rates_per_cont$Method) <- c(levels(avg_rates_per_cont$Method), "Sigma")
avg_rates_per_cont$Method[avg_rates_per_cont$Method == 'sigma'] <- 'Sigma'

avg_rates_per_cont$Method <- factor(avg_rates_per_cont$Method, levels = c('HAYSTAC', 'Sigma', 'Kraken2/Bracken', 'Krakenuniq', 'MALT'))

facet_names <- list(
  "Non Contaminated"='Simple simulations without human DNA',
  "Contaminated"='Simple simulations with human DNA'
)
facet_labeller <- function(variable,value){
  return(facet_names[value])
}

avg_rates_per_cont$Status = factor(avg_rates_per_cont$Status, levels=c("Non Contaminated", "Contaminated"))

# png(paste("SimpleSims", "_FN_FP_True_Detected.png", sep=''), height = 700, width = 1400)

ggsave(paste("SimpleSims", "_FN_FP_True_Detected.pdf", sep=''), height = 11, width = 23, units = "in")

ggplot(avg_rates_per_cont, aes(fill=variable, 
                               y=sapply(value_mean, FUN=function(x) ifelse(x==0, 0.05,x) ), 
                               x=Method)) + 
  geom_bar(position="dodge", stat="identity") + xlab('') +
  geom_text(aes(label=format(round(value_mean, 1),nsmall=1), 
                y=value_mean + avg_rates_per_cont$value_sd/sqrt(nrow(avg_rates_per_cont))), 
            position=position_dodge(width=0.9), vjust=-0.5, size=7) +
  geom_errorbar(aes(ymin=value_mean - avg_rates_per_cont$value_sd/sqrt(nrow(avg_rates_per_cont)), 
                    ymax=value_mean + avg_rates_per_cont$value_sd/sqrt(nrow(avg_rates_per_cont))), 
                width=.2, position=position_dodge(.9)) +
  geom_hline(aes(yintercept=10), linetype="dashed") +
  scale_fill_manual(values=pal_final, name = "", 
                    labels = c("False positive", "False negative", "True positive")) +
  facet_wrap(~Status, ncol=2, labeller=facet_labeller) + 
  #facet_grid(Status~Method, switch='x') +
  ylab('Mean species count') +
  geom_vline(xintercept=seq(1.5, length(unique(avg_rates_per_cont$Method))-0.5, 1), 
             lwd=1, colour="#f0f0f0") +
  theme(text = element_text(size=30, face="bold"), axis.title=element_text(size=32, face="bold")) + 
  theme_Publication_simple(base_size = 28, base_family = "Helvetica") +
  geom_text(aes(x = -Inf, y = Inf, label = label, group = label), size= theme_get()$text[["size"]], hjust = -0.5, vjust = 1.5)

dev.off()



# GENUS AVERAGE RATES FOR ALL SAMPLES FOR ALL THE METHODS

pal_example <- c('#fe9929', '#4575b4')

comparison_abund <- read.csv("Genus_Analysis/Genus_vs_RefSeq_Abundance.tsv", header=TRUE, sep='\t')
comparison_abund$Mean_Posterior_Abundance_Genus <- comparison_abund$Mean_Posterior_Abundance_Genus/100
comparison_abund$Mean_Posterior_Abundance_RefSeq <- comparison_abund$Mean_Posterior_Abundance_RefSeq/100

comparison_abund <- comparison_abund[!comparison_abund$Genus=='Lactobacillus',]

c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black",
  "maroon",
  "gray70", "khaki2",
  "deeppink1", "blue1",
  "darkturquoise",
  "brown"
)

comparison_abund <- comparison_abund %>% arrange(desc(Genus))

# png(paste("Genus_vs_RefSeq_Abund.png", sep=''), height = 650, width = 800)

ggsave("Genus_vs_RefSeq_Abund.pdf", height = 9, width = 10, units = "in")

ggplot(filter(comparison_abund, Mean_Posterior_Abundance_Genus < 0.006), aes(y=Mean_Posterior_Abundance_Genus, x=Mean_Posterior_Abundance_RefSeq, 
                                                                             color=Genus, shape=Genus)) + 
  geom_abline() + geom_point(size=7) +
  scale_shape_manual(values=1:9) +
  scale_color_manual(values=c25) +
  scale_fill_manual(values=c25) +
  theme(text = element_text(size=20, face="bold"), axis.title=element_text(size=22, face="bold")) +
  #  scale_y_continuous(trans = log_trans()) + scale_x_continuous(trans = log_trans()) +
  scale_y_continuous(breaks =c(0, 0.001, 0.002, 0.003, 0.004, 0.005), limits = c(0,0.0055)) +
  scale_x_continuous(breaks =c(0, 0.001, 0.002, 0.003, 0.004, 0.005), limits = c(0,0.0055)) +
  coord_fixed() +
  ylab("Abundance Genus") + xlab("Abundance RefSeq") +
  theme_Publication_simple(base_size = 25, position="right", direction="vertical", 
                           legend.title=element_text(), 
                           legend.text=element_text(face='italic'),
                           base_family = "Helvetica")

dev.off()

ggsave("Genus_vs_RefSeq_Abund_nontrunc.pdf", height = 9, width = 10, units = "in")

ggplot(comparison_abund, aes(y=Mean_Posterior_Abundance_Genus, x=Mean_Posterior_Abundance_RefSeq, 
                                                                             color=Genus, shape=Genus)) + 
  geom_abline() + geom_point(size=7) +
  scale_shape_manual(values=1:9) +
  scale_color_manual(values=c25) +
  scale_fill_manual(values=c25) +
  theme(text = element_text(size=20, face="bold"), axis.title=element_text(size=22, face="bold")) +
  #  scale_y_continuous(trans = log_trans()) + scale_x_continuous(trans = log_trans()) +
  scale_y_continuous(breaks =c(0, 0.001, 0.002, 0.003, 0.004, 0.005), limits = c(0,0.0055)) +
  scale_x_continuous(breaks =c(0, 0.001, 0.002, 0.003, 0.004, 0.005), limits = c(0,0.0055)) +
  coord_fixed() +
  ylab("Abundance Genus") + xlab("Abundance RefSeq") +
  theme_Publication_simple(base_size = 25, position="right", direction="vertical", 
                           legend.title=element_text(), 
                           legend.text=element_text(face='italic'),
                           base_family = "Helvetica")

dev.off()

t.test(comparison_abund$Mean_Posterior_Abundance_Genus, comparison_abund$Mean_Posterior_Abundance_RefSeq, paired=TRUE)

average_rates <- read.csv("Genus_Analysis/Genus_Mean_Count.tsv", header=TRUE, sep='\t')
average_rates <- average_rates[!average_rates$Genus=='Lactobacillus',]

avg_rate_summary <- average_rates %>% group_by(Method) %>% summarise_all(funs(mean, sd))

avg_rate_summary_melt <- melt(as.data.frame(avg_rate_summary[,c("Method", "False_Positive_mean", 
                                                                "False_Negative_mean", 
                                                                "True_Detected_Species_mean")]), 
                              variable.names = "mean", value.names = "value", 
                              id.vars = c("Method"))

sd_rate_sumamry_melt <- melt(as.data.frame(avg_rate_summary[,c("Method", "False_Positive_sd", 
                                                               "False_Negative_sd",
                                                               "True_Detected_Species_sd")]), 
                             variable.names = "sd", value.names = "value", 
                             id.vars = c("Method"))

base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}


levels(avg_rate_summary_melt$Method) <- c(levels(avg_rate_summary_melt$Method), "Kraken2/Bracken")
avg_rate_summary_melt$Method[avg_rate_summary_melt$Method == 'kraken'] <- 'Kraken2/Bracken'

levels(avg_rate_summary_melt$Method) <- c(levels(avg_rate_summary_melt$Method), "Krakenuniq")
avg_rate_summary_melt$Method[avg_rate_summary_melt$Method == 'krakenuniq'] <- 'Krakenuniq'

levels(avg_rate_summary_melt$Method) <- c(levels(avg_rate_summary_melt$Method), "MALT")
avg_rate_summary_melt$Method[avg_rate_summary_melt$Method == 'malt'] <- 'MALT'

levels(avg_rate_summary_melt$Method) <- c(levels(avg_rate_summary_melt$Method), "HAYSTAC")
avg_rate_summary_melt$Method[avg_rate_summary_melt$Method == 'haystac'] <- 'HAYSTAC'

levels(avg_rate_summary_melt$Method) <- c(levels(avg_rate_summary_melt$Method), "Sigma")
avg_rate_summary_melt$Method[avg_rate_summary_melt$Method == 'sigma'] <- 'Sigma'

avg_rate_summary_melt$Method <- factor(avg_rate_summary_melt$Method, levels = c('HAYSTAC', 'Sigma', 'Kraken2/Bracken', 'Krakenuniq', 'MALT'))

# png("Average_Mean_Genus_FP_FN_Rates.png", height = 700, width = 1100)

ggsave("Average_Mean_Genus_FP_FN_Rates.pdf", height = 11, width = 12, units = "in")

avg_rate_summary_melt$title <- "Genus specific database analysis"

ggplot(avg_rate_summary_melt, aes(fill=variable, y=value, x=Method)) + 
  geom_text(aes(label=format(round(value, 1),nsmall=1), 
                y=value + sd_rate_sumamry_melt$value/sqrt(nrow(average_rates))), 
            position=position_dodge(width=0.9), vjust=-0.5, size=6.5) +
  geom_bar(position="dodge", stat="identity") + xlab('') + 
  ylab('Mean species count') +
  facet_wrap(~title) +
  geom_errorbar(aes(ymin=value - sd_rate_sumamry_melt$value/sqrt(nrow(average_rates)), 
                    ymax=value + sd_rate_sumamry_melt$value/sqrt(nrow(average_rates))), 
                width=.2, position=position_dodge(.9)) +
  geom_hline(aes(yintercept=2.01), linetype="dashed") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), breaks =c(0, 1, 2.5, 5, 10, 15, 30)) +
  scale_fill_manual(values=pal_final, name = "", labels = c("False positive", "False negative", "True positive")) + 
  theme(text = element_text(size=25, face="bold"), axis.title=element_text(size=22, face="bold")) +
  theme_Publication_simple(base_size = 25, base_family = "Helvetica")

dev.off()

