library(GenomicRanges)
library(karyoploteR)
library(dplyr)
library(ggplot2)


## FM ratio

df <- read.table("GN18983_remapped_mapq60_sorted_coverage_clean.txt", header = F)
df2 <- read.table("GN18986_remapped_mapq60_sorted_coverage_clean.txt", header = F)
df <- mutate(df, pool = (df$V5 + df2$V5))

dm <- read.table("GN18973_remapped_mapq60_sorted_coverage_clean.txt", header = F)
dm2 <- read.table("GN18972_remapped_mapq60_sorted_coverage_clean.txt", header = F)
dm <- mutate(dm, pool = (dm$V5 + dm2$V5))

# Pooled sequences from each sex normalized using the number of reads mapped in each alignment file
ratio <- data.frame(chr=df[,1], chr_name=df[,2], start=df[,3],
                    female = df$pool, male = dm$pool, wd_size = df[,7],
                    ratio = ((df$pool*(454374336+475772020)) / (dm$pool*(519724485+492814442))))

# Remove all where male read depth == 0
ratio <- ratio[ratio$male!=0,]

chromlen <- read.table("GCF_017639515.1_sCarCar2.pri_genomic.fasta.fai", 
                       header=FALSE)[,c(1,2)] %>%
  arrange(V1) %>% 
  cbind(seq(1:nrow(.)),.)
colnames(chromlen)=c("chr", "chr_name","length")

custom.genome <- GRanges(seqnames=chromlen$chr_name,
                         ranges=IRanges(start=1, end=chromlen$length))
seqlengths(custom.genome) <- chromlen$length


# Source file is the sequence report from the NCBI genome page of sCarCar2.pri
Chr <- read.table(
  file = "CarCar_chr.tsv",
  sep = "\t")
X<-sort(Chr[Chr[,4] == "X",10])
Y<-sort(Chr[Chr[,4] == "Y",10])

xpar <- ratio[ratio$chr_name %in% X,]
xpar_md <- xpar %>%
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

xpar <- left_join(xpar, xpar_md, by = "chr_name")


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

Xplot <- ggplot(xpar, aes(x = chr_name, y = ratio, fill = Type)) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 2, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
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

Yplot <- ggplot(ypar, aes(x = chr_name, y = ratio, fill = Type)) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
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

svg("Ccar_X_box.svg", width = 9, height = 3)
Xplot
dev.off()


svg("Ccar_Y_box.svg", width = 3, height = 3)
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


svg("Ccar_XY_ratio.svg", width = 12, height = 2)
mychr<-c(X[1],Y[1])
kp <- plotKaryotype(genome=custom.genome,
                    chr= mychr,
                    plot.type=3,
                    plot.params = pp,
                    labels.plotter=NULL)

kpAddBaseNumbers(kp, tick.dist = 5000000, add.units= TRUE, cex = 0.75)
kpAxis(kp, ymax=3, tick.pos = c(0,1,2,3), 
       data.panel = 1,side = 2)
kpAddLabels(kp, labels = "F-to-M Depth", 
            data.panel = 1)

kpAbline(kp, h=c(0,1,2,3), ymax=3, lty=2, 
         col="#666666")
kpPoints(kp,
         chr=ratio$chr_name,
         x=ratio$start,
         y=ratio$ratio,
         ymax=3,
         col="purple",
         data.panel = 1)
dev.off()


### unloc XPAR length 
sum(chromlen[chromlen$chr_name %in% as.character(unlist(xpar_md[xpar_md[,3]=="PAR",1])),3])

### unloc length / total length
sum(chromlen[chromlen$chr_name %in% X[c(2:length(X))],3]) / sum(chromlen[chromlen$chr_name %in% X,3])

### Total estimated X length. Additional 2 Mb from the XPAR sequences not present on Y
sum(chromlen[chromlen$chr_name %in% X,3]) + sum(chromlen[chromlen$chr_name %in%
                 as.character(unlist(ypar_md[ypar_md[,3]=="PAR",1])),3]) + 2000000


### unloc YPAR length
sum(chromlen[chromlen$chr_name %in%
               as.character(unlist(ypar_md[ypar_md[,3]=="PAR",1])),3])
sum(chromlen[chromlen$chr_name %in%
               as.character(unlist(ypar_md[ypar_md[,3]=="Non-PAR",1])),3])
### unloc length / total length
sum(chromlen[chromlen$chr_name %in% Y[c(2:length(Y))],3]) / sum(chromlen[chromlen$chr_name %in% Y,3])

### Total estimated Y length. Additional 2 Mb from the XPAR sequences not present on Y
sum(chromlen[chromlen$chr_name %in% Y,3]) + sum(chromlen[chromlen$chr_name %in% 
                 as.character(unlist(xpar_md[xpar_md[,3]=="PAR",1])),3]) + 2000000

### Total number of PAR scaffolds
length(c(as.character(unlist(ypar_md[ypar_md[,3]=="PAR",1])),as.character(unlist(xpar_md[xpar_md[,3]=="PAR",1]))))
