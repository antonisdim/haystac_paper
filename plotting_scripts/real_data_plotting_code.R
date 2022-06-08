library(ggplot2)
library(wesanderson)
library(viridis)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(cowplot)

theme_Publication_heat <- function(base_size=14, base_family="helvetica") {
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
            axis.text.x = element_text(angle=90,hjust=0),
            axis.text.y = element_text(face = 'italic'),
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
            strip.text = element_text(face="bold"), legend.box.just = "center"
    ))
  
}

pal <- wes_palette("Zissou1", 100, type = "continuous")

###### DENTAL CALCULUS ######

HACEK <- read.csv(file="Dent_Calc/DentCalc_HACEK_Pathogens_Abund_Filter.tsv", header=T, sep='\t')

HACEK <- HACEK[(HACEK$Dirichlet_Read_Num >= 20),]

HACEK <- HACEK %>% complete(Taxon, Sample, fill = list(Mean_Posterior_Abundance = NA))

HACEK <- filter(HACEK, Sample != "CSL")

HACEK$Taxon <- gsub(pattern = "_", replacement = " ", HACEK$Taxon)

hacek_order_sample <- HACEK %>% filter(!is.na(Mean_Posterior_Abundance))  %>% count(Sample, wt = Mean_Posterior_Abundance) %>% arrange(desc(n), Sample)
hacek_order_taxon <- HACEK %>% filter(!is.na(Mean_Posterior_Abundance))  %>% count(Taxon, wt = Mean_Posterior_Abundance) %>% arrange(n, Taxon)

Heatmap_HACEK <- ggplot(data=HACEK, aes(x=factor(Sample, level = hacek_order_sample$Sample), y=factor(Taxon, level = hacek_order_taxon$Taxon), fill=Mean_Posterior_Abundance)) + 
  geom_tile() +
  scale_fill_gradientn(colors = pal, na.value = "#f0f0f0",  values = scales::rescale(c(0, 0.1,4)), breaks=c(0.12,0.57,1)) +
  labs(fill = "Mean posterior abundance") + ylab("") + xlab("Sample") +
  scale_x_discrete(expand = c(0, 0), position = 'top') +
  scale_y_discrete(expand = c(0, 0)) +
  theme_classic() +
  theme_Publication_heat(base_size = 20, base_family = "Helvetica") +
  coord_equal(ratio = 1)

# png(paste("HACEK_Taxon_Vs_Sample_Heatmap.png", sep=''), height = 800, width = 1400)
ggsave("HACEK_Taxon_Vs_Sample_Heatmap.pdf", height = 9, width = 18, units = "in")
Heatmap_HACEK
dev.off()


Resp <- read.csv(file="Dent_Calc/DentCalc_Respiratory_Pathogens_Abund_Filter.tsv", header=T, sep='\t')

Resp <- Resp[(Resp$Dirichlet_Read_Num >= 20),]

Resp <- Resp %>% complete(Taxon, Sample, fill = list(Mean_Posterior_Abundance = NA))

Resp$Taxon <- gsub(pattern = "_", replacement = " ", Resp$Taxon)

Respiratory_order_sample <- Resp %>% filter(!is.na(Mean_Posterior_Abundance))  %>% count(Sample, wt = Mean_Posterior_Abundance) %>% arrange(desc(n), Sample)
Respiratory_order_taxon <- Resp %>% filter(!is.na(Mean_Posterior_Abundance))  %>% count(Taxon, wt = Mean_Posterior_Abundance) %>% arrange(n, Taxon)

Heatmap_Respiratory <- ggplot(data=Resp, aes(x=factor(Sample, level = Respiratory_order_sample$Sample), y=factor(Taxon, level = Respiratory_order_taxon$Taxon), fill=Mean_Posterior_Abundance)) + 
  geom_tile() +
  scale_fill_gradientn(colors = pal, na.value = "#f0f0f0", values = scales::rescale(c(0, 0.4,3)), breaks=c(0.3,1.15,2.0)) +
  labs(fill = "Mean posterior abundance") + ylab("") + xlab("Sample") +
  scale_x_discrete(expand = c(0, 0), position = 'top') +
  scale_y_discrete(expand = c(0, 0)) +
  theme_Publication_heat(base_size = 20, base_family = "Helvetica") +
  coord_equal(ratio = 1) 

#png(paste("Respiratory_Taxon_Vs_Sample_Heatmap.png", sep=''), height = 800, width = 1150)
ggsave("Respiratory_Taxon_Vs_Sample_Heatmap.pdf", height = 9, width = 14, units = "in")
Heatmap_Respiratory
dev.off()



###### BRONZE AGE PLAGUE ######


RISE <- read.csv(file="Plague/RISE_Yersinia_Abund_Filter.tsv", header=T, sep='\t')

RISE$Mean_Posterior_Abundance <- RISE$Mean_Posterior_Abundance/100

RISE <- RISE %>% complete(Taxon, Sample, fill = list(Mean_Posterior_Abundance = NA, Dirichlet_Read_Num = NA))

RISE <- filter(RISE, Taxon != "Grey_Matter")
RISE <- filter(RISE, Taxon != "Dark_Matter")

RISE$Taxon <- gsub(pattern = "_", replacement = " ", RISE$Taxon)

RISE_order_sample <- RISE %>% filter(!is.na(Mean_Posterior_Abundance))  %>%
  group_by(Sample) %>% 
  summarise(n = sum(Mean_Posterior_Abundance)) %>% 
  arrange(desc(n), Sample)
RISE_order_taxon <- RISE %>% filter(!is.na(Mean_Posterior_Abundance))  %>%
  group_by(Taxon) %>% 
  summarise(n = sum(Mean_Posterior_Abundance)) %>% 
  arrange(n, Taxon)

Heatmap_RISE <- ggplot(data=RISE, aes(x=factor(Sample, level = RISE_order_sample$Sample), y=factor(Taxon, level = RISE_order_taxon$Taxon), fill=Mean_Posterior_Abundance)) + 
  geom_tile() +
  scale_fill_gradientn(colors = pal, na.value = "#f0f0f0", values = scales::rescale(c(0, 0.01,1)), breaks=c(0.001, 0.005, 0.009)) +
  labs(fill = "Mean posterior \nabundance") + ylab("") + xlab("Sample") +
  scale_x_discrete(expand = c(0, 0), position = 'top') +
  scale_y_discrete(expand = c(0, 0)) +
  theme_Publication_heat(base_size = 22, base_family = "Helvetica") +
  coord_equal(ratio = 1)

#png(paste("RISE_Taxon_Vs_Sample_Heatmap.png", sep=''), height = 800, width = 1800)
ggsave("RISE_Taxon_Vs_Sample_Heatmap.pdf", height = 9, width = 10, units = "in")
Heatmap_RISE
dev.off()


