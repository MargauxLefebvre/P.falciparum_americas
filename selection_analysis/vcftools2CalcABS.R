# Script R : vcftools2CalcABS.R
# Margaux LEFEBVRE
# 10/01/2022
# Goal : creating ABS input from .frq.count (vcftools)

args=(commandArgs(TRUE))
print(args)

# Packages
library(readr)
library(data.table)

#inputs
pop1 <- read_table2(args[1]) #vcftools allele count output for pop1
pop2 <- read_table2(args[2]) #vcftools allele count output for pop2
pop3 <- read_table2(args[3]) #vcftools allele count output for pop3
pop4 <- read_table2(args[4]) #vcftools allele count output for pop4

# fonction for formatting the inputs
ABS_formatting <- function(df) {
  colnames(df)<-c("CHROM","POS","N_ALLELE","SAMPLED","N_ANC")
  df$N_DERIV<- df$SAMPLED-df$N_ANC
  as.data.frame(df)
}

# formatting the inputs for each population
cnt_1<-ABS_formatting(pop1)
cnt_2<-ABS_formatting(pop2)
cnt_3<-ABS_formatting(pop3)
cnt_4<-ABS_formatting(pop4)

# creating the output
finl<-cbind(cnt_1$CHROM, cnt_1$POS, cnt_1$N_DERIV, cnt_1$SAMPLED, cnt_2$N_DERIV, cnt_2$SAMPLED, cnt_3$N_DERIV, cnt_3$SAMPLED, cnt_4$N_DERIV, cnt_4$SAMPLED)
colnames(finl)<-c("CHROM", "POS", "DERIV_1", "SAMPLED_1","DERIV_2", "SAMPLED_2","DERIV_3", "SAMPLED_3","DERIV_4", "SAMPLED_4")
finl<-as.data.frame(finl)
finl<-subset(finl, finl$SAMPLED_1 != 0)
finl<-subset(finl, finl$SAMPLED_2 != 0)
finl<-subset(finl, finl$SAMPLED_3 != 0)
finl<-subset(finl, finl$SAMPLED_4 != 0)
#output
write_tsv(as.data.frame(finl), file = args[5])




