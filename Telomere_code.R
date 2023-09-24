"R version 4.1.0 (2021-05-18)"

library(ggplot2)
library(dplyr)
library(ggtext)


#set working directory#
setwd("C:/Users/miche/OneDrive - The Pennsylvania State University/Documents/Michelle/PhD/PhD project/Telomere/qPCR VS WGS")


#### Reads quality control descriptive analyses####
#Fasta_QC Before trimming#
reads_telo <- read.csv("multiqc_report_beforeTRIM/Generalstat_before trimming.csv")

#compute mean, standard deviation of number of reads#
mean(reads_telo$Total.Sequences)#44739631
sd(reads_telo$Total.Sequences) #8731701

#compute mean, standard deviation of GC content#
mean(reads_telo$X.GC)#36.98
sd(reads_telo$X.GC) #0.8203443


##Fasta_QC After trimming#
reads_telo_QC <- read.csv("multiqc_report_afterTRIM/GeneralStatAfterTrim.csv")

#compute mean, standard deviation of number of reads#
mean(reads_telo_QC$Total.Sequences)#37333533
sd(reads_telo_QC$Total.Sequences) #6591673

#compute mean, standard deviation of GC content#
mean(reads_telo_QC$X.GC)#35.655
sd(reads_telo_QC$X.GC) #0.5361083

#### Telomere lenght estimates from qPCR and computational approaches####

#open telomere data frame#
df<-read.csv("Poplar_telomere_data.csv")
str(df)

#Transforming copy numbers of k_seek and TRIP in mean telomere length#
df$k_seek_bp <- ((df$k_seek)/(19*2))
df$trip_bp <- ((df$TRIP)/(19*2))

#Assesing data distribution#
hist(df$T_S_Ratio)
hist(df$Computel)#365 outlier
hist(df$k_seek_bp)
hist(df$trip_bp)#365 outlier


#eliminate outlier because it has low genome coverage#

df2 <- df %>% filter(PLANT_ID != 365)

#Assesing data distribution#
hist(df2$T_S_Ratio)#365 outlier
hist(df2$Computel)#365 outlier
hist(df2$k_seek_bp)
hist(df2$trip_bp)#365 outlier

#Shapiro-Wilk normality test#
shapiro.test(df2$T_S_Ratio) 
shapiro.test(df2$Computel) 
shapiro.test(df2$k_seek_bp)
shapiro.test(df2$trip_bp)


#transforming data to meet normality and stabilize variation within groups (log10)#
df2$LogTS <- (log10(df2$T_S_Ratio))
df2$LogComputel <- (log10(df2$Computel))
df2$Logk_seek <- (log10(df2$k_seek_bp))
df2$Log_trip <- (log10(df2$trip_bp))

#Shapiro-Wilk normality test in log10#
shapiro.test(df2$LogTS) 
shapiro.test(df2$LogComputel) 
shapiro.test(df2$Logk_seek)
shapiro.test(df2$Log_trip)

##### Descriptive results from qPCR, Computel, k-seek, and trip####

## descriptive results - QPCR##
max(df2$T_S_Ratio)#New values 2.056228##
min(df2$T_S_Ratio)#New values 0.6736168##
mean(df2$T_S_Ratio)#New values 1.184656##
sd(df2$T_S_Ratio) #0.283018

max(df2$LogTS)#New values 0.3130712##
min(df2$LogTS)#New values -0.1715871##
mean(df2$LogTS)#New values 0.06184798##
sd(df2$LogTS) #New values 0.1009953

## descriptive results - computel##
max(df2$Computel) #6184.88
min(df2$Computel) #2160.66
mean(df2$Computel) # 4143.746
sd(df2$Computel) # 794.5442

max(df2$LogComputel) #3.791331
min(df2$LogComputel) #3.334586
mean(df2$LogComputel) # 3.609281
sd(df2$LogComputel) # 0.08536622



## descriptive results - k seek##

max(df2$k_seek_bp) #8807.684
min(df2$k_seek_bp) #1969
mean(df2$k_seek_bp) #5069.59
sd(df2$k_seek_bp) # 1617.943

max(df2$Logk_seek) #3.944862
min(df2$Logk_seek) #3.294246
mean(df2$Logk_seek) #3.681459
sd(df2$Logk_seek) # 0.1474401

#descriptive results - Trip#

max(df2$trip_bp) #14132
min(df2$trip_bp) #3195.632
mean(df2$trip_bp) #8229.17
sd(df2$trip_bp) #2484.83

max(df2$Log_trip) #4.150204
min(df2$Log_trip) #3.504557
mean(df2$Log_trip) #3.894341
sd(df2$Log_trip) # 0.1393044


####Comparing telomere estimates calculated from whole genome sequence data####
#using different bioinformatic approaches#

#Organizing data frame#
Computel <- df2 %>% select(LogComputel, Computel) %>% mutate(Program = "Computel")
names(Computel) <- c("Telomere_length", "TLbp", "Program")
k_seek <- df2 %>% select(Logk_seek, k_seek_bp) %>% mutate(Program = "K seek")
names(k_seek) <- c("Telomere_length", "TLbp", "Program")
trip <- df2 %>% select(Log_trip, trip_bp) %>% mutate(Program = "TRIP")
names(trip) <- c("Telomere_length", "TLbp", "Program")
qpcr <- df2 %>% select(LogTS, T_S_Ratio) %>% mutate(Program = "qpcr")
names(qpcr) <- c("Telomere_length", "TLbp", "Program")

WGS <- rbind(Computel,k_seek, trip)
WGS$Program <- as.factor(WGS$Program)


#visualizing data #
boxplot_WGS <- ggplot(data= WGS, aes(x = Program, y= Telomere_length))+ 
  geom_boxplot(aes(color=Program), alpha=1)+
  geom_jitter(aes(color = Program), shape=16, alpha=0.6,  position=position_jitter(0.2)) +
  scale_color_manual(values = c("#0F52BA", "#999999", "#E69F00"))+
  theme_classic()+ 
  theme(legend.position= "none", axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12))+
  labs( x = '', y = (log[10]~'('~Telomere~Repeat~Copy~Number~')'), size = 12 )+
  annotate("text", x=1, y=4, label= "A")+
  annotate("text", x=2, y=4.2, label= "B")+
  annotate("text", x=3, y=4.3, label= "C"); boxplot_WGS

ggsave("Telomere_articleVs/Figures/Revision1/boxplot_3_bioprograms.pdf", width = 11,  height = 11,  units = "cm")
ggsave("Telomere_articleVs/Figures/Revision1/boxplot_3_bioprograms.png", dpi = 600, width = 11,  height = 11,  units = "cm")


WGS_qpcr <- rbind(Computel,k_seek, trip, qpcr)
WGS_qpcr$Program <- as.factor(WGS_qpcr$Program)


#data:  Telomere_length by Program

#Bartlett Test to test Homogeneity of Variance #
variance<-bartlett.test(Telomere_length ~ Program, data = WGS)
variance

#Bartlett's K-squared = 30.616, df = 2, p-value = 2.248e-07

#Comparing telomere estimates from the three bioinformatic 
#approaches using Kruskal–Wallis 

kruskal.test(Telomere_length ~ Program, data = WGS)
#Kruskal-Wallis chi-squared = 140.95, df = 2, p-value < 2.2e-16

#Dunn test#
library(dunn.test)
dunn.test(WGS$Telomere_length, WGS$Program,method = "bonferroni")


#pairwise correlations no genome correction, figure for supplementary material, Figure S3#
library(GGally)
parwise_cor <- df2 %>% select(PLANT_ID, LogComputel, Logk_seek, Log_trip)
ggpairs(parwise_cor, columns = 2:4,
        lower=list(continuous="points", wrap=c(colour="blue")))+theme_bw()
ggsave("Telomere_articleVs/Figures/Revision1/Parwise_cor.pdf", width = 12.5,  height = 12.5,  units = "cm")
ggsave("Telomere_articleVs/Figures/Revision1/Parwise_cor.png", dpi = 600, width = 12.5,  height = 12.5,  units = "cm")

####Assessing correlation of telomere length assesed by qPCR and bioinformatic approaches#### 

###Computel####

cor.test(df2$LogComputel, df2$LogTS) ##new value 0.6637789 

graph_computel <- ggplot(data= df2, aes(x = LogComputel, y = LogTS))+
  geom_point(fill = '#0F52BA', size= 3,alpha=0.9, shape =21, color ="black")+
  geom_smooth(method = "lm", se = FALSE, color ="#333333")+
  theme_classic()+
  geom_richtext(x = 3.47, y = 0.3, label = "<i>r</i> = 0.66, <i>p</i> < 0.001", label.color = NA, size = 4)+
  theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12))+
  labs( x = (log[10]~'('~Telomere~Repeat~Copy~Number~')'), y = (log[10]~'('~rTL~')') ); graph_computel


####k_seek results###

cor.test(df2$Logk_seek, df2$LogTS) ##new value 0.4989876 

graph_kseek <- ggplot(data= df2, aes(x = Logk_seek, y = LogTS ))+
  geom_point(size= 3, alpha =1, fill ="#999999", shape =21, color ="black")+
  geom_smooth(method = "lm", se = FALSE, color ="#333333")+
  #scale_x_continuous(breaks = seq(from = 4.1, to = 4.8, by = 0.1))+#
  geom_richtext(x = 3.5, y = 0.3, label = "<i>r</i> = 0.50, <i>p</i> < 0.001", label.color = NA, size = 4)+
  theme_classic()+
  theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12))+
  labs( x = (log[10]~'('~Telomere~Repeat~Copy~Number~')'), y = (log[10]~'('~rTL~')'));graph_kseek


###TRIP Results#

cor.test(df2$Log_trip, df2$LogTS) ##new value 0.5208103

graph_trip <- ggplot(data= df2, aes(x = Log_trip, y = LogTS ))+
  geom_point(size= 3, alpha = 0.9, fill ="#E69F00", shape =21, color ="black")+
  geom_smooth(method = "lm", se = FALSE, color ="#333333")+
  theme_classic()+
  #scale_x_continuous(breaks = seq(from = 4.3, to = 5.0, by = 0.1)) +#
  geom_richtext(x = 3.7, y = 0.3, label = "<i>r</i> = 0.52, <i>p</i> < 0.001", label.color = NA, size = 4)+
  theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12))+
  labs( x = (log[10]~'('~Telomere~Repeat~Copy~Number~')'), y = (log[10]~'('~rTL~')')); graph_trip

####graph trip&k_seek&computel ####
#join boxplot for computel,  k_eek and trip kmers relationship with T/S#

library(ggpubr)
cor_qPcr_WGS <- ggarrange(graph_computel, graph_trip + rremove("ylab"), graph_kseek + rremove("ylab"), 
          hjust = -0.5,
          vjust = 1,
          align = "v",
          labels = 'AUTO',
          ncol = 3, nrow = 1)

ggsave("Telomere_articleVs/Figures/Revision1/cor_qPcr_WGS.pdf", width = 24.5,  height = 8.5,  units = "cm")
ggsave("Telomere_articleVs/Figures/Revision1/cor_qPcr_WGS.png", dpi = 600, width = 24.5,  height = 8.5,  units = "cm")

####Equation 2 - genome coverage correction####

#Adding genome coverage in k_seek and TRIP#

#K-seek#
df2$k_seek_Cpb <- df2$k_seek_bp/df2$Coverage #results in pb
max(df2$k_seek_Cpb) #302.7736
min(df2$k_seek_Cpb) #105.3505
mean(df2$k_seek_Cpb) #203.2711
sd(df2$k_seek_Cpb) #42.95003


df2$k_seek_C <- log10(df2$k_seek_bp/df2$Coverage)

cor.test(df2$k_seek_C, df2$LogTS) #  0.6175783 
str(df2)
graph_k_seek_C <- ggplot(data= df2, aes(x = k_seek_C, y = LogTS))+
  geom_point(size= 3, alpha = 0.9, color ="black", fill='#999999',shape =21 )+
  geom_smooth(method = "lm", se = FALSE, color ="#333333")+
  scale_x_continuous(breaks = seq(from = 2, to = 2.5, by = 0.1))+
  geom_richtext(x = 2.17, y = 0.3, label = "<i>r</i> = 0.62, <i>p</i> < 0.001", label.color = NA, size = 4)+
  theme_classic()+
  theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12))+
  labs( x = (log[10]~'('~Telomere~Repeat~Copy~Number~')'), y = (log[10]~'('~rTL~')')); graph_k_seek_C

#TRIP#
df2$TRIP_Cpb <- df2$trip_bp/df2$Coverage #results in pb
max(df2$TRIP_Cpb) #485.8027
min(df2$TRIP_Cpb) #170.9808
mean(df2$TRIP_Cpb) #330.7605
sd(df2$TRIP_Cpb) #64.52415

df2$TRIP_C <- log10(df2$trip_bp/df2$Coverage)
cor.test(df2$TRIP_C, df2$LogTS) #  0.6562003 

graph_trip_C <- ggplot(data= df2, aes(x = TRIP_C, y = LogTS ))+
  geom_point(size= 3, alpha = 0.9, fill ="#E69F00", shape =21, color ="black")+
  geom_smooth(method = "lm", se = FALSE, color ="#333333")+
  scale_x_continuous(breaks = seq(from = 2.2, to = 2.7, by = 0.1))+
  geom_richtext(x = 2.37, y = 0.3, label = "<i>r</i> = 0.66, <i>p</i> < 0.001", label.color = NA, size = 4)+
  theme_classic()+
  theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12))+
  labs( x = (log[10]~'('~Telomere~Repeat~Copy~Number~')'), y = (log[10]~'('~rTL~')')); graph_trip_C

####graph trip&k_seek coverage####

corrected_TL <- ggarrange(graph_computel,graph_trip_C+rremove('ylab'), graph_k_seek_C +rremove('ylab'),
          hjust = -0.5,
          vjust = 1,
          align = "v",
          labels = 'AUTO', ncol = 3, nrow = 1)
ggsave("Telomere_articleVs/Figures/Revision1/cor_qPcr_WGS_GC.pdf", width = 30.5,  height = 8.5,  units = "cm")
ggsave("Telomere_articleVs/Figures/Revision1/cor_qPcr_WGS_GC.png", dpi = 600, width = 30.5,  height = 8.5,  units = "cm")


##Boxplot of corrected telomere estimates (eq 2)#

Computel_C <- df2 %>% select(LogComputel, Computel) %>% mutate(Program = "Computel")
names(Computel_C) <- c("Telomere_length", "TLbp", "Program")
k_seek_C <- df2 %>% select(k_seek_C, k_seek_bp) %>% mutate(Program = "K seek")
names(k_seek_C) <- c("Telomere_length", "TLbp", "Program")
trip_C <- df2 %>% select(TRIP_C, trip_bp) %>% mutate(Program = "TRIP")
names(trip_C) <- c("Telomere_length", "TLbp", "Program")

WGS_C <- rbind(Computel_C,k_seek_C, trip_C)
WGS_C$Program <- as.factor(WGS_C$Program)


#visualizing data #
boxplot_WGS_C <- ggplot(data= WGS_C, aes(x = Program, y= Telomere_length))+ 
  geom_boxplot(aes(color=Program), alpha=1)+
  geom_jitter(aes(color = Program), shape=16, alpha=0.6,  position=position_jitter(0.2)) +
  scale_color_manual(values = c("#0F52BA", "#999999", "#E69F00"))+
  theme_classic()+ 
  theme(legend.position= "none", axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12))+
  labs( x = '', y = (log[10]~'('~Telomere~Repeat~Copy~Number~')'), size = 12 )+
  annotate("text", x=1, y=4, label= "A")+
  annotate("text", x=2, y=2.6, label= "B")+
  annotate("text", x=3, y=2.8, label= "C"); boxplot_WGS_C

ggsave("Telomere_articleVs/Figures/Revision1/boxplot_3_bioprograms_C.pdf", width = 12.5,  height = 12.5,  units = "cm")
ggsave("Telomere_articleVs/Figures/Revision1/boxplot_3_bioprograms_C.png", dpi = 600, width = 12.5,  height = 12.5,  units = "cm")

#Comparing telomere estimates from the three bioinformatic 
#approaches using Kruskal–Wallis 

kruskal.test(Telomere_length ~ Program, data = WGS_C)
#Kruskal-Wallis chi-squared = 251.69, df = 2, p-value < 2.2e-16

#Dunn test#
dunn.test(WGS_C$Telomere_length, WGS_C$Program,method = "bonferroni")


##parwaise correlations across approaches###
parwise_cor_C <- df2 %>% select(PLANT_ID, LogComputel, k_seek_C, TRIP_C)
cor_three_programs<- ggpairs(parwise_cor_C, columns = 2:4, 
        columnLabels = c("Computel", "K-seek", "TRIP"),
        lower=list(continuous="points", wrap=c(colour="blue")))+theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x = (log[10]~'('~Telomere~Repeat~Copy~Number~')'), y = (log[10]~'('~Telomere~Repeat~Copy~Number~')'))
cor_1 <- ggally_text( "r = 0.96*", color = I("black"))
cor_2 <- ggally_text( "r = 0.99*", color = I("black"))
cor_3 <- ggally_text( "r = 0.97*", color = I("black"))
cor_three_programs[1, 2]<-cor_1
cor_three_programs[1, 3]<-cor_2
cor_three_programs[2, 3]<-cor_3


ggsave("Telomere_articleVs/Figures/Revision1/pairwise_cor_C.pdf", width = 12.5,  height = 12.5,  units = "cm")
ggsave("Telomere_articleVs/Figures/Revision1/pairwise_cor_C.png", dpi = 600, width = 12.5,  height = 12.5,  units = "cm")


#correlation between telomere length and genome coverage#

cor(df2$Computel, df2$Coverage) #0.4487291

#Equation1
cor(df2$k_seek_bp, df2$Coverage) #0.8082736
cor(df2$trip_bp, df2$Coverage) #0.8096098

#equation2 
cor(df2$k_seek_Cpb, df2$Coverage) # 0.4691928
cor(df2$TRIP_Cpb, df2$Coverage)#0.4329931

