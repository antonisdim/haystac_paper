setwd("~/haystac_results/results_elife/logs/")

library(ggplot2)
library(reshape2)
library(wesanderson)
library(ggforce)
library(dplyr)
library(scales)
library(gridExtra)
library(ggpubr)

pal_final <- c('#d7191c', '#fdae61', '#abd9e9', '#2c7bb6')

theme_Publication <- function(base_size=14, base_family="Helvetica") {
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
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(1, "cm"),
            legend.key.width = unit(1.5, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}


prefix1 <- c("DB")
prefix2 <- c("Haystack_conda")
prefix3 <- c("Sample")

suffix <- c("_benchmark_new.txt")

db <- read.csv(paste(prefix1, suffix, sep=''), header=TRUE, sep='\t')
db <- db[!db[, "Method"] == 'haystac_conda',]

#db <- db %>% group_by(Database, Method) %>% 
#  summarize(Runtime_m = mean(Runtime_m), Max.RSS..KB. = mean(Max.RSS..KB.))

#db$Max.RSS..KB. <- db$Max.RSS..MB. * 1000

sample <- read.csv(paste(prefix3, suffix, sep=''), header=TRUE, sep='\t')
sample <- sample[!sample[, "Method"] == 'haystac_conda',]

#sample <- sample %>% group_by(Database, Method, Number.of.Reads) %>% 
#  summarize(Runtime_m = mean(Runtime_m), Max.RSS..KB. = mean(Max.RSS..KB.))

#sample$Max.RSS..KB. <- sample$Max.RSS..MB. * 1000 

db <- db[, c('Database', 'Method', 'RSS', 'Runtime_m')]
# https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=649402
db$RSS <- db$RSS/4

db_melt <- melt(as.data.frame(db), variable.names = "variable", value.names = "value", 
                id.vars = c("Database", "Method"))

db_melt$variable <- ifelse(grepl("RSS", db_melt$variable), 
                           "Memory", 'Runtime')

sample <- as.data.frame(sample[, c('Number.of.Reads', 'Database', 'Method', 
                     'RSS', 'Runtime_m')])
sample$RSS <- sample$RSS/4

sample_500 <- sample[sample[, "Database"] == 500,]

sample_500 <- sample_500[, c('Number.of.Reads', 'Method', 
                             'RSS', 'Runtime_m')]

sample_500_melt <-melt(sample_500, variable.names = "variable", 
                       value.names = "value", 
                       id.vars = c("Number.of.Reads", "Method"))

sample_500_melt$variable <- ifelse(grepl("RSS", sample_500_melt$variable), 
                                   "Memory", 'Runtime')


sample_1M <- sample[sample[, "Number.of.Reads"] == 1000000,]

sample_1M <- sample_1M[, c('Database', 'Method', 
                           'RSS', 'Runtime_m')]

sample_1M_melt <-melt(sample_1M, variable.names = "variable", 
                      value.names = "value", 
                      id.vars = c("Database", "Method"))

sample_1M_melt$variable <- ifelse(grepl("RSS", sample_1M_melt$variable), 
                                  "Memory", 'Runtime')

db_melt <- db_melt %>% rename(xaxis=Database)
sample_500_melt <- sample_500_melt %>% rename(xaxis=Number.of.Reads)
sample_1M_melt <- sample_1M_melt %>% rename(xaxis=Database)

db_melt$category <- 'Database'
sample_500_melt$category <- 'Analysis'
sample_1M_melt$category <- 'Analysis'

total <- rbind(db_melt, sample_1M_melt)

levels(total$Method) <- c(levels(total$Method), "Kraken2/Bracken")
total$Method[total$Method == 'kraken'] <- 'Kraken2/Bracken'

levels(total$Method) <- c(levels(total$Method), "Krakenuniq")
total$Method[total$Method == 'krakenuniq'] <- 'Krakenuniq'

levels(total$Method) <- c(levels(total$Method), "MALT")
total$Method[total$Method == 'malt'] <- 'MALT'

levels(total$Method) <- c(levels(total$Method), "Haystac")
total$Method[total$Method == 'haystac'] <- 'Haystac'

levels(total$Method) <- c(levels(total$Method), "Sigma")
total$Method[total$Method == 'sigma'] <- 'Sigma'     

total$Method <- factor(total$Method, levels = c('Haystac', 'Sigma', 
                                                'Kraken2/Bracken', 'Krakenuniq', 'MALT'))
total$category <- factor(total$category, levels = c('Database', 
                                                    'Analysis'))

# png(paste("total.png", sep=''), height = 900, width = 1000)

ggsave("benchmark_total.pdf", height = 10, width = 14, units = "in")

ggplot(total, aes(x=xaxis, y=value, color=Method)) + 
  stat_smooth(method='lm', se=F, size = 1.5, 
              fullrange = T) + geom_point(shape=21, size=2) + 
  ylab('') + 
  scale_y_continuous(trans = "log10",
                     labels=c("0.167" = "10 sec", "1.0" = "1 min", "10.0" = "10 min", 
                     "100.0" = "100 min", "1e+05" = "100 MB", 
                     "1e+06" = "1 GB", "1e+07" = "10 GB", "1e+08" = "100 GB"), 
                     breaks=c(0.167, 1.0, 10.0, 100.0, 1e+05, 1e+06, 1e+07, 1e+8)) +
  scale_x_continuous(trans = "log10", 
                     breaks=c(10,100, 500, 5000)) + 
  facet_grid(variable~category, scales='free', labeller = function (labels) {
    labels <- lapply(labels, as.character)
    list(do.call(paste, c(labels, list(sep = "\n"))))}) +
  xlab("Database size (species)") + 
  theme_Publication(base_size = 22, base_family = "Helvetica") +
  theme(panel.border = element_rect(colour = "grey", fill = NA))
      

dev.off()


# png(paste("sample_500.png", sep=''), height = 900, width = 1000)

levels(sample_500_melt$Method) <- c(levels(sample_500_melt$Method), "Kraken2/Bracken")
sample_500_melt$Method[sample_500_melt$Method == 'kraken'] <- 'Kraken2/Bracken'

levels(sample_500_melt$Method) <- c(levels(sample_500_melt$Method), "Krakenuniq")
sample_500_melt$Method[sample_500_melt$Method == 'krakenuniq'] <- 'Krakenuniq'

levels(sample_500_melt$Method) <- c(levels(sample_500_melt$Method), "MALT")
sample_500_melt$Method[sample_500_melt$Method == 'malt'] <- 'MALT'

levels(sample_500_melt$Method) <- c(levels(sample_500_melt$Method), "Haystac")
sample_500_melt$Method[sample_500_melt$Method == 'haystac'] <- 'Haystac'

levels(sample_500_melt$Method) <- c(levels(sample_500_melt$Method), "Sigma")
sample_500_melt$Method[sample_500_melt$Method == 'sigma'] <- 'Sigma'     

sample_500_melt$Method <- factor(sample_500_melt$Method, levels = c('Haystac', 'Sigma', 
                                                                    'Kraken2/Bracken', 'Krakenuniq', 'MALT'))
sample_500_melt$category <- factor(sample_500_melt$category, levels = c('Database', 
                                                    'Analysis'))

ggsave("benachmark_sample_500sp.pdf", height = 9, width = 11, units = "in")

ggplot(sample_500_melt, aes(x=xaxis, y=value, color=Method)) + 
  stat_smooth(method='lm', se=F, size = 2, fullrange = T, size =1.5) + geom_point(size = 2, shape = 21) +
  ylab('') + 
  scale_y_continuous(trans = "log10",
                     labels=c("0.167" = "10 sec", "1.0" = "1 min", "10.0" = "10 min", 
                              "100.0" = "100 min", "1e+05" = "100 MB", 
                              "1e+06" = "1 GB", "1e+07" = "10 GB", "1e+08" = "100 GB"), 
                     breaks=c(0.167, 1.0, 10.0, 100.0, 1e+05, 1e+06, 1e+07, 1e+8)) +
  scale_x_continuous(trans = "log10", breaks=c(10000, 100000, 1000000)) +
  facet_grid(variable~category, scales='free', labeller = function (labels) {
    labels <- lapply(labels, as.character)
    list(do.call(paste, c(labels, list(sep = "\n"))))}) +
  xlab("Sample size") + 
  theme_Publication(base_size = 20, base_family = "Helvetica") +
  theme(panel.border = element_rect(colour = "grey", fill = NA))

dev.off()