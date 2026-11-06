library(GenomicRanges)
library(karyoploteR)
library(dplyr)
library(ggplot2)




# Input files

df <- read.table("Hypanus-sabinus-GN6226_mapq60_sorted_coverage_clean.txt", header = F)
df2 <- read.table("Hypanus-sabinus-GN21844_mapq60_sorted_coverage_clean.txt", header = F)
df <- mutate(df, pool = (df$V5 + df2$V5))

dm <- read.table("Hypanus-sabina-GN21920_mapq60_sorted_coverage_clean.txt", header = F)
dm2 <- read.table("Hypanus-sabinus-GN21826_mapq60_sorted_coverage_clean.txt", header = F)
dm <- mutate(dm, pool = (dm$V5 + dm2$V5))

# Pooled sequences from each sex normalized using the number of reads mapped in each alignment file
ratio <- data.frame(chr=df[,1], chr_name=df[,2], start=df[,3],
                    female = df$pool, male = dm$pool, wd_size = df[,7],
                    ratio = ((df$pool*(400530412+297371103)) / (dm$pool*(362415393+457772565))))

# Remove all where male read depth == 0
ratio <- ratio[ratio$male!=0,]

# Source file is the sequence report from the NCBI genome page of sHypSab1.hap1
Chr <- read.table(
  file = "HypSab_chr.tsv",
  sep = "\t")

X1<-sort(Chr[Chr[,4] == "X1",10])
X2<-sort(Chr[Chr[,4] == "X2",10])
Y<-sort(Chr[Chr[,4] == "Y",10])

x1par <- ratio[ratio$chr_name %in% X1,]
x1par_md <- x1par %>%
  group_by(chr_name) %>%
  summarise(
    p50 = quantile(ratio, 0.5),
  ) %>%
  mutate(
    Type = case_when(
      p50 < 1.8 & p50 > 1.2 ~ "Uncertain",
      TRUE ~ ifelse(abs(p50 - 1) < abs(p50 - 2), "PAR", "Non-PAR")
    )
  ) %>% print(n=35)

x1par <- left_join(x1par, x1par_md, by = "chr_name")

#X2
x2par <- ratio[ratio$chr_name %in% X2,]
x2par_md <- x2par %>%
  group_by(chr_name) %>%
  summarise(
    p50 = quantile(ratio, 0.5),
  ) %>%
  mutate(
    Type = case_when(
      p50 < 1.8 & p50 > 1.2 ~ "Uncertain",
      TRUE ~ ifelse(abs(p50 - 1) < abs(p50 - 2), "PAR", "Non-PAR")
    )
  ) %>% print(n=35)

x2par <- left_join(x2par, x2par_md, by = "chr_name")

ypar <- ratio[ratio$chr_name %in% Y,]
ypar_md <- ypar %>%
  group_by(chr_name) %>%
  summarise(
    p50 = quantile(ratio, 0.5),
  ) %>%
  mutate(
    Type = case_when(
      p50 < 0.8 & p50 > 0.2 ~ "Uncertain",
      TRUE ~ ifelse(abs(p50 - 1) < abs(p50 - 0), "PAR", "Non-PAR")
    )
  ) %>% print(n=35)

ypar <- left_join(ypar, ypar_md, by = "chr_name")


X1plot <- ggplot(x1par, aes(x = chr_name, y = ratio, fill = Type))+
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red")+
  geom_hline(yintercept = 2, linetype = "dashed", color = "blue")+
  coord_cartesian(ylim = c(0, 5)) + 
  scale_fill_manual(
    values = c("PAR" = "red", "Non-PAR" = "blue", Uncertain = "grey80")
  )+
  labs(x = "Scaffold", y = "F-to-M Ratio", 
  ) +
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
  )

X2plot <- ggplot(x2par, aes(x = chr_name, y = ratio, fill = Type))+
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red")+
  geom_hline(yintercept = 2, linetype = "dashed", color = "blue")+
  coord_cartesian(ylim = c(0, 5)) + 
  scale_fill_manual(
    values = c("PAR" = "red", "Non-PAR" = "blue", Uncertain = "grey80")
  )+
  labs(x = "Scaffold", y = "F-to-M Ratio", 
  ) +
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
  )

Yplot <- ggplot(ypar, aes(x = chr_name, y = ratio, fill = Type))+
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue")+
  coord_cartesian(ylim = c(0, 5)) + 
  scale_fill_manual(
    values = c("PAR" = "red", "Non-PAR" = "blue", Uncertain = "grey80")
  )+
  labs(x = "Scaffold", y = "F-to-M Ratio", 
  )+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
  )

svg("Hsab_X1_box.svg", width = 6, height = 3)
X1plot
dev.off()

svg("Hsab_X2_box.svg", width = 6, height = 3)
X2plot
dev.off()

svg("Hsab_Y_box.svg", width = 6, height = 3)
Yplot
dev.off()



pp <- getDefaultPlotParams(plot.type=3)
pp$leftmargin <- 0.2
pp$rightmargin <- 0.1
pp$topmargin <- 30
pp$bottommargin <- 10
pp$data2height <- 20
pp$ideogramheight <- 2
pp$ideogramlateralmargin <- 0.005


chromlen <- read.table("GCF_030144855.1_sHypSab1.hap1_genomic.fna.fai", header=FALSE)[,c(1,2)] %>%
  arrange(V1) %>% 
  cbind(seq(1:nrow(.)),.)
colnames(chromlen)=c("chr", "chr_name","length")

custom.genome <- GRanges(seqnames=chromlen$chr_name, ranges=IRanges(start=1, end=chromlen$length))
seqlengths(custom.genome) <- chromlen$length

mychr <- c(X1[1],X2[1],Y[1])

svg("Hsab_XY_ratio.svg", width = 12, height = 2)
kp <- plotKaryotype(genome=custom.genome,
                    chr= mychr,
                    plot.type=3,
                    plot.params = pp,
                    labels.plotter=NULL)
kpAddBaseNumbers(kp, tick.dist = 10000000, add.units= TRUE, cex = 0.75)

kpAxis(kp, ymax=3, tick.pos = c(0,1,2,3), 
       data.panel = 1,side = 2)
kpAddLabels(kp, 
            labels = "F-to-M depth", 
            data.panel = 1)

kpPoints(kp,
         chr=ratio$chr_name,
         x=ratio$start,
         y=ratio$ratio,
         ymax=3,
         col="purple",
         data.panel = 1)
kpAbline(kp, h=c(0,1,2,3), ymax=3, lty=2, 
         col="#666666")

dev.off()



### unloc X1PAR length 
sum(chromlen[chromlen$chr_name %in% 
               as.character(unlist(x1par_md[x1par_md[,3]=="PAR",1])),3])

### unloc length / total length
sum(chromlen[chromlen$chr_name %in% X1[c(2:length(X1))],3]) /
  sum(chromlen[chromlen$chr_name %in% X1,3])

### Total estimated X1 length
sum(chromlen[chromlen$chr_name %in% X1,3]) +
  sum(chromlen[chromlen$chr_name %in%
                 as.character(unlist(ypar_md[ypar_md[,3]=="PAR",1])),3])


### unloc X2PAR length 
sum(chromlen[chromlen$chr_name %in% 
               as.character(unlist(x2par_md[x2par_md[,3]=="PAR",1])),3])

### unloc length / total length
sum(chromlen[chromlen$chr_name %in% X2[c(2:length(X2))],3]) /
  sum(chromlen[chromlen$chr_name %in% X2,3])

### Total estimated X2 length
sum(chromlen[chromlen$chr_name %in% X2,3]) +
  sum(chromlen[chromlen$chr_name %in%
                 c(as.character(unlist(ypar_md[ypar_md[,3]=="PAR",1]))),3])


### unloc YPAR length
sum(chromlen[chromlen$chr_name %in%
               c(as.character(unlist(ypar_md[ypar_md[,3]=="PAR",1]))),3])
sum(chromlen[chromlen$chr_name %in%
               c(as.character(unlist(ypar_md[ypar_md[,3]=="Non-PAR",1]))),3])
### unloc length / total length
sum(chromlen[chromlen$chr_name %in% Y[c(2:length(Y))],3]) /
  sum(chromlen[chromlen$chr_name %in% Y,3])

### Total estimated Y length
sum(chromlen[chromlen$chr_name %in% Y,3]) +
  sum(chromlen[chromlen$chr_name %in% 
                 c(as.character(unlist(x2par_md[x2par_md[,3]=="PAR",1]))),3]) +
  sum(chromlen[chromlen$chr_name %in% 
                 c(as.character(unlist(x1par_md[x1par_md[,3]=="PAR",1]))),3])
