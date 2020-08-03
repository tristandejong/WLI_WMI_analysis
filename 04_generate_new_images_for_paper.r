#install.packages("package:matrixStats")
#install.packages("circlize")

library(circlize)
library(UpSetR)
library("matrixStats")


rowMins = function(df){
  out = apply(df, 1, FUN=min)
  return(out)
}

rowMax = function(df){
  out = apply(df, 1, FUN=max)
  return(out)
}




################################
# Figure 1. Depth of coverage  #
################################

chr_order = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,"X","Y")


first = TRUE

headerNames = c("WLI\n IonProton ",
                "WLI\n Chromium",
                "WLI\n Illumina xTen",
                "WMI\n IonProton",
                "WMI\n Chromium",
                "WMI\n Illumina xTen")

fileNames = c("WLI_IonProton.bam.csv",
              "WLI_realn.bam.csv",
              "WLI_xTen.bam.csv",
              "WMI_IonProton.bam.csv",
              "WMI_realn.bam.csv",
              "WMI_xTen.bam.csv")
barplotCol = c("indianred1","red2","red4","steelblue2","blue","navyblue")

pdf("final_figure1.pdf", width = 10, height = 7)

par(mfrow=c(2,3))


first = TRUE
for( i in 1:6){
  
  
  print(paste("processing: ",i))
  depth = read.table(paste("table1/",fileNames[i],sep = ""), sep = ";", stringsAsFactors = TRUE, row.names = 1, header = T)
  
  depth = depth[chr_order,]
  barplot(depth[,3], names = rownames(depth), las = 2, 
          ylab = "Mean depth of coverage", xlab = "Chromosome",
          main = headerNames[i], 
          ylim = c(0,50), cex.lab=1.2, cex.axis = 1,
          col = barplotCol[i])  
  
  if (first == TRUE){
    total_depth = depth[,3]  
    total_coverage = sum(depth[,2])
    first = FALSE
  }else{
    total_depth = cbind(total_depth,depth[,3])
    total_coverage = c(total_coverage, sum(depth[,2]))
  }
}

dev.off()

################################
# LOAD IN THE DATA             #
################################

chr_list = c(1:20, "X", "Y")


first = TRUE
for (chr in chr_list){
    genotype_file = paste0("data/deepvariant_wmiwli_chr",chr,"_6_samples.gvcf.gt") 
    genotype = read.table(genotype_file, sep = "\t", stringsAsFactors = FALSE, row.names = 1, header = T)

    quality_file =  paste0("data/deepvariant_wmiwli_chr",chr,"_6_samples.gvcf.gq")
    quality = read.table(quality_file, sep = "\t", stringsAsFactors = FALSE, row.names = 1, header = T)

    depth_file =  paste0("data/deepvariant_wmiwli_chr",chr,"_6_samples.gvcf.ad")
    depth = read.table(depth_file, sep = "\t", stringsAsFactors = FALSE, row.names = 1, header = T)
    
    print(paste("chr",chr))
    print(paste("nrow: ", nrow(genotype)))
    
    ##################################################################################
    ## Check the number of different calls we have:
    table(c(genotype[,1], genotype[,2], genotype[,3], genotype[,4], genotype[,5], genotype[,6]))
    
    ##################################################################################
    
    if (first == TRUE){
      first = FALSE
      all_genotype = genotype
      all_quality = quality
      all_depth = depth
      
    }else{
      all_genotype = rbind(all_genotype, genotype)
      all_quality = rbind(all_quality, quality)
      all_depth = rbind(all_depth, depth)
    }
}

write.csv(all_genotype, "backup_genotype.csv")
write.csv(all_depth, "backup_depth.csv")
write.csv(all_quality, "backup_quality.csv")

#all_genotype = read.csv( "backup_genotype.csv", row.names = 1)
#all_depth = read.csv( "backup_depth.csv", row.names = 1)
#all_quality = read.csv( "backup_quality.csv", row.names = 1)

dim(all_genotype)
dim(all_depth)
dim(all_quality)



head(all_genotype)


#Investigation of sample: 1_106355362_A_C
all_genotype["1_106355362_A_C",]
all_quality["1_106355362_A_C",]

# Conclusion, there is an internal dispute, calling it ./. with a quality higher than 10.

# Investigation of sample: 8_58998305_C_CA
all_genotype["8_58998305_C_CA",]
all_quality["8_58998305_C_CA",]

# Conclusion, there is an internal dispute, calling it 0/1 with a quality higher than 10.


# Now comes the question: Do we want to keep a sample like this, as there clearly is a SNP of sorts...


############################################################################################
# 1. Group by different classes
############################################################################################


# We make our selection based on the highest quality call by deepvar. We can also st this to 
call_split = strsplit(all_genotype$call_qual,"_")
WLI = as.character(sapply(call_split, "[", 1))
WMI = as.character(sapply(call_split, "[", 2))


## The alt calls with proper opposition
WLI_ALT =       WLI %in% c("1/1","2/2","3/3", "4/4") &  WMI %in% c("0/0")
#all_genotype[WLI_ALT ,]
WMI_ALT =       WMI %in% c("1/1","2/2","3/3", "4/4") &  WLI %in% c("0/0")
#all_genotype[WMI_ALT ,]

# The Het calls with proper opposition
WLI_HET =       WLI %in% c("0/1","0/2","0/3","0/4") &  WMI %in% c("0/0")
#all_genotype[WLI_HET ,]
WMI_HET =       WMI %in% c("0/1","0/2","0/3","0/4") &  WLI %in% c("0/0")
#all_genotype[WMI_HET ,]


# . = n/a, M = Missing data, 
# P = Partial data, 
# I = gVCF input site is non-called, 
# D = insufficient Depth of coverage,
# - = unrepresentable overlapping deletion, 
# L = Lost/unrepresentable allele (other than deletion), 
# U = multiple Unphased variants present, 
# O = multiple Overlapping variants present, 
# 1 = site is Monoallelic, no assertion about presence of REF or ALT allele"

table(all_genotype$RNC_WLI_chromium)

RNC_WLI = all_genotype

# The Het calls with proper opposition
WLI_UNOPPOSED =       WMI %in% c("./.")
#all_genotype[WLI_UNOPPOSED ,]
WMI_UNOPPOSED =       WLI %in% c("./.")
#all_genotype[WMI_UNOPPOSED ,]

table(all_genotype[WMI_UNOPPOSED ,]$call_qual)


###########################################################################################################
## 1. Get tables of REF, HET, ALT
###########################################################################################################

# Whip them columns
depth_REF = ifelse(all_genotype[,1] == "0/0", all_depth[,1] , 0 )
for (i in 2:6){
  single_depth_REF = ifelse(all_genotype[,i] == "0/0", all_depth[,i] , 0 )
  depth_REF = cbind(depth_REF, single_depth_REF)
}
rownames(depth_REF) = rownames(all_depth)
colnames(depth_REF) = colnames(all_depth)

depth_HET = ifelse(all_genotype[,1]  %in% c("0/1","0/2","0/3","0/4"), all_depth[,1] , 0 )
for (i in 2:6){
  single_depth_HET = ifelse(all_genotype[,i]  %in% c("0/1","0/2","0/3","0/4"), all_depth[,i] , 0 )
  depth_HET = cbind(depth_HET, single_depth_HET)
}
rownames(depth_HET) = rownames(all_depth)
colnames(depth_HET) = colnames(all_depth)

depth_ALT = ifelse(all_genotype[,1] %in% c("0/0","1/1","2/2","3/3","4/4"), all_depth[,1] , 0 )
for (i in 2:6){
  single_depth_ALT = ifelse(all_genotype[,i] %in% c("0/0","1/1","2/2","3/3","4/4"), all_depth[,i] , 0 )
  depth_ALT = cbind(depth_ALT, single_depth_ALT)
}
rownames(depth_ALT) = rownames(all_depth)
colnames(depth_ALT) = colnames(all_depth)


head(depth_REF)
head(depth_HET)
head(depth_ALT)

########################################################

# Whip them columns # pt2
quality_REF = ifelse(all_genotype[,1] == "0/0", all_quality[,1] , 0 )
for (i in 2:6){
  single_quality_REF = ifelse(all_genotype[,i] == "0/0", all_quality[,i] , 0 )
  quality_REF = cbind(quality_REF, single_quality_REF)
}
rownames(quality_REF) = rownames(all_quality)
colnames(quality_REF) = colnames(all_quality)

quality_HET = ifelse(all_genotype[,1]  %in% c("0/1","0/2","0/3","0/4"), all_quality[,1] , 0 )
for (i in 2:6){
  single_quality_HET = ifelse(all_genotype[,i]  %in% c("0/1","0/2","0/3","0/4"), all_quality[,i] , 0 )
  quality_HET = cbind(quality_HET, single_quality_HET)
}
rownames(quality_HET) = rownames(all_quality)
colnames(quality_HET) = colnames(all_quality)

quality_ALT = ifelse(all_genotype[,1] %in% c("0/0","1/1","2/2","3/3","4/4"), all_quality[,1] , 0 )
for (i in 2:6){
  single_quality_ALT = ifelse(all_genotype[,i] %in% c("0/0","1/1","2/2","3/3","4/4"), all_quality[,i] , 0 )
  quality_ALT = cbind(quality_ALT, single_quality_ALT)
}
rownames(quality_ALT) = rownames(all_quality)
colnames(quality_ALT) = colnames(all_quality)


head(quality_REF)
head(quality_HET)
head(quality_ALT)



###########################################################################################################
## 1. Barplot of quality
###########################################################################################################

barPlot_color1 = c(rep("indianred1",3), rep("steelblue2",3)) #WMI_IonProton WMI_chromium WMI_xTen
barPlot_color2 = c(rep("red2",3),rep("blue",3)) #WMI_IonProton WMI_chromium WMI_xTen
barPlot_color3 = c(rep("red4",3),rep("navyblue",3)) #WMI_IonProton WMI_chromium WMI_xTen

barplotCol = c("indianred1","red2","red4","steelblue2","blue","navyblue")
barPlotNames = c("IonProton", "10x Chromium", "Illumina xTen", "IonProton",  "10x Chromium", "Illumina xTen")
par(mfrow=c(1,2), pty = "s")



###########################################
# Figure 0a. Number of high confidence SNPS
###########################################


nreads_10 = colSums(all_depth > 10)
nreads_20 = colSums(all_depth > 20)
nreads_30 = colSums(all_depth > 30)
nreads_40 = colSums(all_depth > 40)

nreads_10
nreads_20
nreads_30
nreads_40

barplot(nreads_10, las = 2, col = )
barplot(nreads_20, las = 2, add=TRUE)
barplot(nreads_30, las = 2, add=TRUE)


# Median depth of coverage.
plot(density(ifelse(depth_REF == 0, NA, depth_REF)[,1],na.rm = TRUE))

plot(density(ifelse(depth_HET == 0, NA, depth_HET)[,1],na.rm = TRUE))
plot(density(ifelse(depth_HET == 0, NA, depth_HET)[,2],na.rm = TRUE))
plot(density(ifelse(depth_HET == 0, NA, depth_HET)[,3],na.rm = TRUE))

plot(density(ifelse(depth_ALT == 0, NA, depth_ALT)[,1],na.rm = TRUE))


REF_d = colMedians(ifelse(depth_REF == 0, NA, depth_REF),na.rm = TRUE)
HET_d = colMedians(ifelse(depth_HET == 0, NA, depth_HET),na.rm = TRUE)
ALT_d = colMedians(ifelse(depth_ALT == 0, NA, depth_ALT),na.rm = TRUE)
ALL_d = colMedians(ifelse(as.matrix(all_depth) == 0, NA, as.matrix(all_depth)),na.rm = TRUE)

REF_d
HET_d
ALT_d
ALL_d

test = ifelse(as.matrix(all_depth) == 0, NA, as.matrix(all_depth))

plot(density(test[,1], na.rm = TRUE))
lines(density(test[,2], na.rm = TRUE))
lines(density(test[,3], na.rm = TRUE))



###########################################
# Figure 2a. Number of high confidence SNPS
###########################################
pdf("final_figure2.pdf", width = 10, height = 6)

par(mfrow=c(1,2), mar = c(10, 5, 4, 3)) # c(bottom, left, top, right))

ALT_10 = colSums(quality_ALT >= 10, na.rm =T)
ALT_20 = colSums(quality_ALT >= 20, na.rm =T)
ALT_30 = colSums(quality_ALT >= 30, na.rm =T)

barplot(ALT_10, add = F, col = barPlot_color1, las = 2, names = barPlotNames, ylim = c(0,450000),axes=FALSE, 
        main = "a                                                                ", las = 2)
axis(2,at=seq(0,500000, 50000),las = 2)
barplot(ALT_20, add = T, col = barPlot_color2, axes = FALSE, las = 2, names = barPlotNames,axisnames = FALSE)
barplot(ALT_30, add = T, col = barPlot_color3, axes = FALSE, las = 2, names = barPlotNames,axisnames = FALSE)

###########################################
# Figure 2b. Number of high confidence SNPS
###########################################

HET_10 = colSums(quality_HET >= 10, na.rm =T)
HET_20 = colSums(quality_HET >= 20, na.rm =T)
HET_30 = colSums(quality_HET >= 30, na.rm =T)

barplot(HET_10, add = F, col = barPlot_color1, las = 2, names = barPlotNames, ylim = c(0,120000),axes=FALSE, 
        main = "b                                                                ", las = 2)
axis(2,at=seq(0,500000, 10000),las = 2)
barplot(HET_20, add = T, col = barPlot_color2, axes = FALSE, las = 2, names = barPlotNames,axisnames = FALSE)
barplot(HET_30, add = T, col = barPlot_color3, axes = FALSE, las = 2, names = barPlotNames,axisnames = FALSE)

plot(x = 0)
 legend(x ="topleft", legend = c("WLI Quality > 10"," WLI Quality > 20", "WLI Quality > 30",
                                 "WMI Quality > 10"," WMI Quality > 20", "WMI Quality > 30"), col = barplotCol, pch = 15)

dev.off()

###########################################################################################################
## 2. Quality of the final calls
###########################################################################################################

max_quality_WLI = rowMax(all_quality[,1:3])
max_quality_WMI = rowMax(all_quality[,4:6])

## FOR THE WLI, The highest quality, in line with the call
WLI_ALT_QUALITY_10 = sum(max_quality_WLI[WLI_ALT] >= 10, na.rm =T)
WLI_ALT_QUALITY_20 = sum(max_quality_WLI[WLI_ALT] >= 20, na.rm =T)
WLI_ALT_QUALITY_30 = sum(max_quality_WLI[WLI_ALT] >= 30, na.rm =T)

WLI_HET_QUALITY_10 = sum(max_quality_WLI[WLI_HET] >= 10, na.rm =T)
WLI_HET_QUALITY_20 = sum(max_quality_WLI[WLI_HET] >= 20, na.rm =T)
WLI_HET_QUALITY_30 = sum(max_quality_WLI[WLI_HET] >= 30, na.rm =T)

WLI_UNOPPOSED_QUALITY_10 = sum(max_quality_WLI[WLI_UNOPPOSED] >= 10, na.rm =T)
WLI_UNOPPOSED_QUALITY_20 = sum(max_quality_WLI[WLI_UNOPPOSED] >= 20, na.rm =T)
WLI_UNOPPOSED_QUALITY_30 = sum(max_quality_WLI[WLI_UNOPPOSED] >= 30, na.rm =T)


## FOR THE WMI, The highest quality, in line with the call
WMI_ALT_QUALITY_10 = sum(max_quality_WMI[WMI_ALT] >= 10, na.rm =T)
WMI_ALT_QUALITY_20 = sum(max_quality_WMI[WMI_ALT] >= 20, na.rm =T)
WMI_ALT_QUALITY_30 = sum(max_quality_WMI[WMI_ALT] >= 30, na.rm =T)

WMI_HET_QUALITY_10 = sum(max_quality_WMI[WMI_HET] >= 10, na.rm =T)
WMI_HET_QUALITY_20 = sum(max_quality_WMI[WMI_HET] >= 20, na.rm =T)
WMI_HET_QUALITY_30 = sum(max_quality_WMI[WMI_HET] >= 30, na.rm =T)

WMI_UNOPPOSED_QUALITY_10 = sum(max_quality_WMI[WMI_UNOPPOSED] >= 10, na.rm =T)
WMI_UNOPPOSED_QUALITY_20 = sum(max_quality_WMI[WMI_UNOPPOSED] >= 20, na.rm =T)
WMI_UNOPPOSED_QUALITY_30 = sum(max_quality_WMI[WMI_UNOPPOSED] >= 30, na.rm =T)

par(mfrow=c(1,1), pty = "s",mar=c(8.1,4.1,4.1,2.1))

barplot(c(WLI_ALT_QUALITY_10,WLI_HET_QUALITY_10,WLI_UNOPPOSED_QUALITY_10, 
          WMI_ALT_QUALITY_10,WMI_HET_QUALITY_10,WMI_UNOPPOSED_QUALITY_10), 
        add = F, col = barPlot_color1, names = c("ALT","HET","UN-OPPOSED","ALT","HET","UN-OPPOSED"), 
        las = 2, ylim = c(0,60000), axes=FALSE, 
        main = "Number of calls per quality")
axis(2,at=seq(0,600000, 5000),las = 2)
barplot(c(WLI_ALT_QUALITY_20,WLI_HET_QUALITY_20,WLI_UNOPPOSED_QUALITY_20,
          WMI_ALT_QUALITY_20,WMI_HET_QUALITY_20,WMI_UNOPPOSED_QUALITY_20), add = T, col = barPlot_color2, las = 2, axes=FALSE)
barplot(c(WLI_ALT_QUALITY_30,WLI_HET_QUALITY_30,WLI_UNOPPOSED_QUALITY_30,
          WMI_ALT_QUALITY_30,WMI_HET_QUALITY_30,WMI_UNOPPOSED_QUALITY_30), add = T, col = barPlot_color3, las = 2, axes=FALSE)

# WLI_ALT
# WMI_ALT
# 
# WLI_HET
# WMI_HET
# 
# WLI_UNOPPOSED
# WMI_UNOPPOSED

###########################################################################################################
## 3. Venn diagram, do the methods agree?
###########################################################################################################

library(VennDiagram)
# library(tidyverse)
# library(hrbrthemes)
# library(tm)
# library(proustr)


WMI_color_2 = c("blue","steelblue2","navyblue") #WMI_IonProton WMI_chromium WMI_xTen
WLI_color_2 = c("red2","indianred1","red4") #WLI_IonProton WLI_chromium WLI_xTen 


#####################################################
## Venn diagram, REF
#####################################################
for (k in 1:4){
      WLI_Ion_REF = rownames(quality_REF[quality_REF[,1] >= k*10,])
      WLI_chrom_REF = rownames(quality_REF[quality_REF[,2] >= k*10,])
      WLI_xTen_REF = rownames(quality_REF[quality_REF[,3] >= k*10,])
      
      WMI_Ion_REF = rownames(quality_REF[quality_REF[,4] >= k*10,])
      WMI_chrom_REF = rownames(quality_REF[quality_REF[,5] >= k*10,])
      WMI_xTen_REF = rownames(quality_REF[quality_REF[,6] >= k*10,])
      
      n1 = paste("WLI ionProton", "(",length(WLI_Ion_REF),")")
      n2 = paste("WLI chromium", "(",length(WLI_chrom_REF),")")
      n3 = paste("WLI Xten", "(",length(WLI_xTen_REF),")")
      
      figure2Base = paste("figureX_Q",k, "_", sep = "")
      
      venn.diagram(
        x = list( "WLI ionProton" = WLI_Ion_REF,
                  "WLI chromium" = WLI_chrom_REF,
                  "WLI Xten" = WLI_xTen_REF),
        category.names = c(n1,n2,n3),
        filename = paste(figure2Base, '_REF.png',sep = ""),
        output = TRUE ,
        imagetype="png" ,
        height = 480 , 
        width = 480 , 
        resolution = 300,
        compression = "lzw",
        lwd = 1,
        scaled = FALSE,
        col=WLI_color_2,
        fill =   adjustcolor( WLI_color_2, alpha.f = 0.2),
        cex = 0.5,
        fontfamily = "sans",
        cat.cex = 0.3,
        cat.default.pos = "outer",
        cat.pos = c(-27, 27, 135),
        cat.dist = c(0.055, 0.055, 0.085),
        cat.fontfamily = "sans",
        cat.col = c(WLI_color_2),
        rotation = 1
      )
}

#####################################################
## Venn diagram, HET
#####################################################
for (k in 1:4){
  WLI_Ion_HET = rownames(quality_HET[quality_HET[,1] >= k*10,])
  WLI_chrom_HET = rownames(quality_HET[quality_HET[,2] >= k*10,])
  WLI_xTen_HET = rownames(quality_HET[quality_HET[,3] >= k*10,])
  
  WMI_Ion_HET = rownames(quality_HET[quality_HET[,4] >= k*10,])
  WMI_chrom_HET = rownames(quality_HET[quality_HET[,5] >= k*10,])
  WMI_xTen_HET = rownames(quality_HET[quality_HET[,6] >= k*10,])
  
  n1 = paste("WLI ionProton", "(",length(WLI_Ion_HET),")")
  n2 = paste("WLI chromium", "(",length(WLI_chrom_HET),")")
  n3 = paste("WLI Xten", "(",length(WLI_xTen_HET),")")
  
  figure2Base = paste("figureX_Q",k, "_", sep = "")
  
  venn.diagram(
    x = list( "WLI ionProton" = WLI_Ion_HET,
              "WLI chromium" = WLI_chrom_HET,
              "WLI Xten" = WLI_xTen_HET),
    category.names = c(n1,n2,n3),
    filename = paste(figure2Base, '_HET.png',sep = ""),
    output = TRUE ,
    imagetype="png" ,
    height = 480 , 
    width = 480 , 
    resolution = 300,
    compression = "lzw",
    lwd = 1,
    scaled = FALSE,
    col=WLI_color_2,
    fill =   adjustcolor( WLI_color_2, alpha.f = 0.2),
    cex = 0.5,
    fontfamily = "sans",
    cat.cex = 0.3,
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    cat.col = c(WLI_color_2),
    rotation = 1
  )
}

#####################################################
## Venn diagram, ALT
#####################################################


for (k in 1:4){
  WLI_Ion_ALT = rownames(quality_ALT[quality_ALT[,1] >= k*10,])
  WLI_chrom_ALT = rownames(quality_ALT[quality_ALT[,2] >= k*10,])
  WLI_xTen_ALT = rownames(quality_ALT[quality_ALT[,3] >= k*10,])
  
  WMI_Ion_ALT = rownames(quality_ALT[quality_ALT[,4] >= k*10,])
  WMI_chrom_ALT = rownames(quality_ALT[quality_ALT[,5] >= k*10,])
  WMI_xTen_ALT = rownames(quality_ALT[quality_ALT[,6] >= k*10,])
  
  n1 = paste("WLI ionProton", "(",length(WLI_Ion_ALT),")")
  n2 = paste("WLI chromium", "(",length(WLI_chrom_ALT),")")
  n3 = paste("WLI Xten", "(",length(WLI_xTen_ALT),")")
  
  figure2Base = paste("figureX_Q",k, "_", sep = "")
  
  venn.diagram(
    x = list( "WLI ionProton" = WLI_Ion_ALT,
              "WLI chromium" = WLI_chrom_ALT,
              "WLI Xten" = WLI_xTen_ALT),
    category.names = c(n1,n2,n3),
    filename = paste(figure2Base, '_ALT.png',sep = ""),
    output = TRUE ,
    imagetype="png" ,
    height = 480 , 
    width = 480 , 
    resolution = 300,
    compression = "lzw",
    lwd = 1,
    scaled = FALSE,
    col=WLI_color_2,
    fill =   adjustcolor( WLI_color_2, alpha.f = 0.2),
    cex = 0.5,
    fontfamily = "sans",
    cat.cex = 0.3,
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    cat.col = c(WLI_color_2),
    rotation = 1
  )
}

#####################################################################




###########################################################################################################
## 4. How often is the call called into doubt?
###########################################################################################################

WLI # This is the final call, based on quality

WMI # This is the final call, based on quality

# EDIT HERE. "HAO"

min_score = 10


# All we really want to filter, is that an ALT is not rejected by a REF, so we should account for that instead.

test1 = all_genotype[,1] == "0/0" & all_genotype[,1] != WLI & all_quality[,1] >= min_score
test2 = all_genotype[,2] == "0/0" &  all_genotype[,2] != WLI & all_quality[,2] >= min_score
test3 = all_genotype[,3] == "0/0" &  all_genotype[,3] != WLI & all_quality[,3] >= min_score

test4 = all_genotype[,4] == "0/0" & all_genotype[,4] != WMI & all_quality[,4] >= min_score
test5 = all_genotype[,5] == "0/0" & all_genotype[,5] != WMI & all_quality[,5] >= min_score
test6 = all_genotype[,6] == "0/0" & all_genotype[,6] != WMI & all_quality[,6] >= min_score

disqualifier = (test1 | test2 | test3 | test4 | test5 | test6)


#second_dispute
test11 = all_genotype[,1] != "0/0" & all_genotype[,1] != WLI & all_quality[,1] >= 30
test12 = all_genotype[,2] != "0/0" &  all_genotype[,2] != WLI & all_quality[,2] >= 30
test13 = all_genotype[,3] != "0/0" &  all_genotype[,3] != WLI & all_quality[,3] >= 30

test14 = all_genotype[,4] != "0/0" & all_genotype[,4] != WMI & all_quality[,4] >= 30
test15 = all_genotype[,5] != "0/0" & all_genotype[,5] != WMI & all_quality[,5] >= 30
test16 = all_genotype[,6] != "0/0" & all_genotype[,6] != WMI & all_quality[,6] >= 30

disqualifier2 = (test11 | test12 | test13 | test14 | test15 | test16)



# ALSO DISCUSS THIS WITH HAO.


# We now have our test that shows how many should be disqualified
## Should we trash calls that are disputed by other techniques?

###########################################################################################################
## 5a. Generate the final calls ALT
###########################################################################################################
final_WLI_ALT_10 = all_genotype[WLI_ALT & max_quality_WLI >=10 & max_quality_WLI < 20 & !disqualifier,]
final_WLI_ALT_20 = all_genotype[WLI_ALT & max_quality_WLI >=20 & max_quality_WLI < 30 & !disqualifier,]
final_WLI_ALT_30 = all_genotype[WLI_ALT & max_quality_WLI >=30 & !disqualifier,]

final_WMI_ALT_10 = all_genotype[WMI_ALT & max_quality_WMI >=10 & max_quality_WMI < 20 & !disqualifier,]
final_WMI_ALT_20 = all_genotype[WMI_ALT & max_quality_WMI >=20 & max_quality_WMI < 30 & !disqualifier,]
final_WMI_ALT_30 = all_genotype[WMI_ALT & max_quality_WMI >=30 & !disqualifier,]


final_WLI_ALT_30B = all_genotype[WLI_ALT & max_quality_WLI >=30 & !disqualifier & !disqualifier2,]
final_WMI_ALT_30B = all_genotype[WMI_ALT & max_quality_WMI >=30 & !disqualifier & !disqualifier2,]


dim(final_WLI_ALT_10)
dim(final_WLI_ALT_20)
dim(final_WLI_ALT_30)

dim(final_WMI_ALT_10)
dim(final_WMI_ALT_20)
dim(final_WMI_ALT_30)


final_QUALITY_WLI_ALT_30 = all_quality[WLI_ALT & max_quality_WLI >=30 & !disqualifier,]
final_QUALITY_WLI_ALT_30B = all_quality[WLI_ALT & max_quality_WLI >=30 & !disqualifier & !disqualifier2,]

final_WLI_HET_30B = all_genotype[WLI_HET & max_quality_WLI >=30 & !disqualifier & !disqualifier2,]



head(final_WLI_ALT_30)
head(final_QUALITY_WLI_ALT_30)


head(final_WLI_ALT_30B)
head(final_QUALITY_WLI_ALT_30B)

dim(final_WLI_ALT_30B)
dim(final_WLI_HET_30B)

###########################################################################################################
## 5b. Generate the final calls HET
###########################################################################################################
final_WLI_HET_10 = all_genotype[WLI_HET & max_quality_WLI >=10 & max_quality_WLI < 20 & !disqualifier,]
final_WLI_HET_20 = all_genotype[WLI_HET & max_quality_WLI >=20 & max_quality_WLI < 30 & !disqualifier,]
final_WLI_HET_30 = all_genotype[WLI_HET & max_quality_WLI >=30 & !disqualifier,]

final_WMI_HET_10 = all_genotype[WMI_HET & max_quality_WMI >=10 & max_quality_WMI < 20  & !disqualifier,]
final_WMI_HET_20 = all_genotype[WMI_HET & max_quality_WMI >=20 & max_quality_WMI < 30  & !disqualifier,]
final_WMI_HET_30 = all_genotype[WMI_HET & max_quality_WMI >=30 & !disqualifier,]

dim(final_WLI_HET_10)
dim(final_WLI_HET_20)
dim(final_WLI_HET_30)

dim(final_WMI_HET_10)
dim(final_WMI_HET_20)
dim(final_WMI_HET_30)

###########################################################################################################
## 5b. Generate the final calls UN-opposed
###########################################################################################################
# . = n/a, M = Missing data, 
# P = Partial data, 
# I = gVCF input site is non-called, 
# D = insufficient Depth of coverage,
# - = unrepresentable overlapping deletion, 
# L = Lost/unrepresentable allele (other than deletion), 
# U = multiple Unphased variants present, 
# O = multiple Overlapping variants present, 
# 1 = site is Monoallelic, no assertion about presence of REF or ALT allele"

allowed_unopposed = c( "--",  "-.", "-L", "..",  "II")#, "L-", "L.", "LL")
not_allowed = c("1.", "11","O.","OO","UU","1.", "11")


for(i in 10:15){
    print(table(all_genotype[,i]))
}

#Secondary disqualifier
testA = (all_quality[,1] < min_score)  | (all_quality[,1] >= min_score & all_genotype[,10] %in% allowed_unopposed) # maybe add II?
testB = (all_quality[,2] < min_score)  | (all_quality[,2] >= min_score & all_genotype[,11] %in% allowed_unopposed)
testC = (all_quality[,3] < min_score)  | (all_quality[,3] >= min_score & all_genotype[,12] %in% allowed_unopposed)

testD = (all_quality[,4] < min_score)  | (all_quality[,4] >= min_score & all_genotype[,13] %in% allowed_unopposed) # maybe add II?
testE = (all_quality[,5] < min_score)  | (all_quality[,5] >= min_score & all_genotype[,14] %in% allowed_unopposed)
testF = (all_quality[,6] < min_score)  | (all_quality[,6] >= min_score & all_genotype[,15] %in% allowed_unopposed)

disqualifierA = (testA | testB | testC)
disqualifierB = (testD | testE | testF)

final_WLI_UNOPPOSED_10 = all_genotype[WLI_UNOPPOSED & max_quality_WLI >=10 & max_quality_WLI < 20 & !disqualifier & disqualifierB,]
final_WLI_UNOPPOSED_20 = all_genotype[WLI_UNOPPOSED & max_quality_WLI >=20 & max_quality_WLI < 30 & !disqualifier & disqualifierB,]
final_WLI_UNOPPOSED_30 = all_genotype[WLI_UNOPPOSED & max_quality_WLI >=30 & !disqualifier & disqualifierB,]

final_WMI_UNOPPOSED_10 = all_genotype[WMI_UNOPPOSED & max_quality_WMI >=10 & max_quality_WMI < 20 & !disqualifier & disqualifierA,]
final_WMI_UNOPPOSED_20 = all_genotype[WMI_UNOPPOSED & max_quality_WMI >=20 & max_quality_WMI < 30 & !disqualifier & disqualifierA,]
final_WMI_UNOPPOSED_30 = all_genotype[WMI_UNOPPOSED & max_quality_WMI >=30 & !disqualifier & disqualifierA,]

dim(final_WLI_UNOPPOSED_10)
dim(final_WLI_UNOPPOSED_20)
dim(final_WLI_UNOPPOSED_30)

dim(final_WMI_UNOPPOSED_10)
dim(final_WMI_UNOPPOSED_20)
dim(final_WMI_UNOPPOSED_30)

# write a namelist
# generate other images

# write.table(rownames(final_WLI_ALT_30),"results/ALT_WLI_calls_names.list", row.names = FALSE, col.names = F, quote = FALSE)
# write.table(rownames(final_WMI_ALT_30),"results/ALT_WMI_calls_names.list", row.names = FALSE, col.names = F, quote = FALSE)
# 
# write.table(rownames(final_WLI_HET_30),"results/HET_WLI_calls_names.list", row.names = FALSE, col.names = F, quote = FALSE)
# write.table(rownames(final_WMI_HET_30),"results/HET_WMI_calls_names.list", row.names = FALSE, col.names = F, quote = FALSE)
# 
# write.table(rownames(final_WLI_UNOPPOSED_30),"results/UNOPPOSED_WLI_calls_names.list", row.names = FALSE, col.names = F, quote = FALSE)
# write.table(rownames(final_WMI_UNOPPOSED_30),"results/UNOPPOSED_WMI_calls_names.list", row.names = FALSE, col.names = F, quote = FALSE)

write.table(rownames(final_WLI_ALT_30),"results/ALT_WLI_calls_names.list", row.names = FALSE, col.names = F, quote = FALSE)
write.table(rownames(final_WMI_ALT_30),"results/ALT_WMI_calls_names.list", row.names = FALSE, col.names = F, quote = FALSE)

write.table(rownames(final_WLI_HET_30),"results/HET_WLI_calls_names.list", row.names = FALSE, col.names = F, quote = FALSE)
write.table(rownames(final_WMI_HET_30),"results/HET_WMI_calls_names.list", row.names = FALSE, col.names = F, quote = FALSE)

write.table(rownames(final_WLI_UNOPPOSED_30),"results/UNOPPOSED_WLI_calls_names.list", row.names = FALSE, col.names = F, quote = FALSE)
write.table(rownames(final_WMI_UNOPPOSED_30),"results/UNOPPOSED_WMI_calls_names.list", row.names = FALSE, col.names = F, quote = FALSE)


## Final ALT calls
dim(final_WLI_ALT_30)
dim(final_WMI_ALT_30)

dim(final_WLI_ALT_20)
dim(final_WMI_ALT_20)

## Final HET calls
dim(final_WLI_HET_30)
dim(final_WMI_HET_30)

dim(final_WLI_HET_20)
dim(final_WMI_HET_20)

## Final HET calls
dim(final_WLI_UNOPPOSED_30)
dim(final_WMI_UNOPPOSED_30)

dim(final_WLI_UNOPPOSED_20)
dim(final_WMI_UNOPPOSED_20)


###################################################################################
# 6. BARPLOT OF FINAL CUTOFF
###################################################################################


par(mfrow=c(1,1), pty = "s",mar=c(8.1,4.1,4.1,2.1))

barplot(c(nrow(final_WLI_ALT_10),nrow(final_WLI_HET_10),nrow(final_WLI_UNOPPOSED_10), 
          nrow(final_WMI_ALT_10),nrow(final_WMI_HET_10),nrow(final_WMI_UNOPPOSED_10)),
        add = F, col = barPlot_color1, names = c("ALT","HET","UN-OPPOSED","ALT","HET","UN-OPPOSED"), 
        las = 2, ylim = c(0,40000), axes=FALSE, 
        main = "Number of final calls per quality")
axis(2,at=seq(0,500000, 5000),las = 2)

barplot(c(nrow(final_WLI_ALT_20),nrow(final_WLI_HET_20),nrow(final_WLI_UNOPPOSED_20), 
          nrow(final_WMI_ALT_20),nrow(final_WMI_HET_20),nrow(final_WMI_UNOPPOSED_20)), add = T, col = barPlot_color2, las = 2, axes=FALSE)
barplot(c(nrow(final_WLI_ALT_30),nrow(final_WLI_HET_30),nrow(final_WLI_UNOPPOSED_30), 
          nrow(final_WMI_ALT_30),nrow(final_WMI_HET_30),nrow(final_WMI_UNOPPOSED_30)), add = T, col = barPlot_color3, las = 2, axes=FALSE)


head(final_WMI_HET_30)

# Investigation of HET calls!
table(final_WMI_HET_30[,1])
table(final_WMI_HET_30[,2])
table(final_WMI_HET_30[,3])

# Investigation of HET calls!
table(final_WMI_HET_30[,4])
table(final_WMI_HET_30[,5])
table(final_WMI_HET_30[,6])

quality[rownames(final_WMI_HET_30),]

################################################################################
## Circos, for depth
################################################################################
mySplit = strsplit(rownames(all_depth),"_")
chrom_number = as.character(sapply(mySplit, "[", 1))
chrom_position = as.numeric(sapply(mySplit, "[", 2))
chr = 1

subSelect = seq( 0, length(chrom_position), 100)
chrom_number_sub = chrom_number[subSelect]
chrom_position_sub = chrom_position[subSelect]

depth_WLI = rowSums(all_depth[subSelect,1:3])
depth_WMI = rowSums(all_depth[subSelect,4:6])

WLI_color = "red"
WMI_color = "blue"

plot(-100,-100, 
     xlim = c(0,max(chrom_position)),
     ylim = c(0,100),
     xlab = "chromosomal position", ylab = "smoothed coverage",
     pch = 16, cex = 0.5, main = paste("chr",chr))

# # WLI

# # WMI
s = (chrom_number_sub == chr & is.na(depth_WMI) == FALSE)
x <- chrom_position_sub[s]
y <- depth_WMI[s]
lo1 <- loess(y~x, span = 0.05)
lines(x = x , y = predict(lo1), col=WMI_color, lwd=2)

## Loop to create tables for circos





################################################################################
## Circos, for quality
################################################################################

subSelect = seq( 0, length(chrom_position), 100)
chrom_number_sub = chrom_number[subSelect]
chrom_position_sub = chrom_position[subSelect]

quality_WLI = all_quality[subSelect,1:3]
quality_WMI = all_quality[subSelect,4:6]

WLI_color_2 = c("pink","red","red4")

plot(-100,-100, 
     xlim = c(0,max(chrom_position)),
     ylim = c(0,50),
     xlab = "chromosomal position", ylab = "smoothed quality",
     pch = 16, cex = 0.5, main = paste("chr",chr))

# # WLI
for (i in 1:3){
  s = (chrom_number_sub == chr & is.na(quality_WLI[,i]) == FALSE)
  x <- chrom_position_sub[s]
  y <- quality_WLI[s,i]
  lo1 <- loess(y~x, span = 0.05)
  lines(x = x , y = predict(lo1), col=WLI_color_2[i], lwd=2)
}
# 


plot(-100,-100, 
     xlim = c(0,max(chrom_position)),
     ylim = c(0,50),
     xlab = "chromosomal position", ylab = "smoothed quality",
     pch = 16, cex = 0.5, main = paste("chr",chr))
for (i in 1:3){
  s = (chrom_number_sub == chr & is.na(quality_WMI[,i]) == FALSE)
  x <- chrom_position_sub[s]
  y <- quality_WMI[s,i]
  lo1 <- loess(y~x, span = 0.05)
  lines(x = x , y = predict(lo1), col=WMI_color_2[i], lwd=2)
}


################################################################################
## Circos, for FINAL ALT/HET CALLS
################################################################################
## Get the locations of ALT calls for WLI
WLI_call = as.character(sapply(strsplit(final_WLI_ALT_30$call_qual,"_"), "[", 1))
WLI_chr =  as.character(sapply(strsplit(rownames(final_WLI_ALT_30),"_"), "[", 1))
WLI_pos =  as.numeric(as.character(sapply(strsplit(rownames(final_WLI_ALT_30),"_"), "[", 2)))

## Get the locations of ALT calls for WMI
WMI_call = as.character(sapply(strsplit(final_WMI_ALT_30$call_qual,"_"), "[", 2))
WMI_chr =  as.character(sapply(strsplit(rownames(final_WMI_ALT_30),"_"), "[", 1))
WMI_pos =  as.numeric(as.character(sapply(strsplit(rownames(final_WMI_ALT_30),"_"), "[", 2)))

## Prep for drawing
draw_WLI_ALT = rep(1, length((WLI_call)))
draw_WMI_ALT = rep(-1, length((WMI_call)))
col_WLI_ALT = rep(adjustcolor( "red", alpha.f = 0.2), length((WLI_call)))
col_WMI_ALT = rep(adjustcolor( "blue", alpha.f = 0.2), length((WMI_call)))

## Create alt table for plotting
alt_table_1 = cbind(WLI_chr, WLI_pos, draw_WLI_ALT, col_WLI_ALT)
alt_table_2 = cbind(WMI_chr, WMI_pos, draw_WMI_ALT, col_WMI_ALT)
alt_table = rbind(alt_table_1,alt_table_2)

################################################################################

## Get the locations of ALT calls for WLI
WLI_call2 = as.character(sapply(strsplit(final_WLI_HET_30$call_qual,"_"), "[", 1))
WLI_chr2 =  as.character(sapply(strsplit(rownames(final_WLI_HET_30),"_"), "[", 1))
WLI_pos2 =  as.numeric(as.character(sapply(strsplit(rownames(final_WLI_HET_30),"_"), "[", 2)))

## Get the locations of HET calls for WMI
WMI_call2 = as.character(sapply(strsplit(final_WMI_HET_30$call_qual,"_"), "[", 2))
WMI_chr2 =  as.character(sapply(strsplit(rownames(final_WMI_HET_30),"_"), "[", 1))
WMI_pos2 =  as.numeric(as.character(sapply(strsplit(rownames(final_WMI_HET_30),"_"), "[", 2)))

## Prep for drawing
draw_WLI_HET = rep(1, length((WLI_call2)))
draw_WMI_HET = rep(-1, length((WMI_call2)))
col_WLI_HET = rep(adjustcolor( "red", alpha.f = 0.2), length((WLI_call2)))
col_WMI_HET = rep(adjustcolor( "blue", alpha.f = 0.2), length((WMI_call2)))

## Create HET table for plotting
het_table_1 = cbind(WLI_chr2, WLI_pos2, draw_WLI_HET, col_WLI_HET)
het_table_2 = cbind(WMI_chr2, WMI_pos2, draw_WMI_HET, col_WMI_HET)
het_table = rbind(het_table_1,het_table_2)



################################################################################
## Start building the Circos plot
################################################################################
mySplit = strsplit(rownames(all_depth),"_")
chrom_number = as.character(sapply(mySplit, "[", 1))
chrom_position = as.numeric(sapply(mySplit, "[", 2))
#chr = 1

subSelect = seq( 0, length(chrom_position), 10)
chrom_number_sub = chrom_number[subSelect]
chrom_position_sub = chrom_position[subSelect]

depth_WLI = rowSums(all_depth[subSelect,1:3])
depth_WMI = rowSums(all_depth[subSelect,4:6])

quality_WLI = rowMax(all_quality[subSelect,1:3])
quality_WMI =  rowMax(all_quality[subSelect,4:6])

# Get chr, 
table_circos = data.frame()

all_chromosomes = c(1:20, "X","Y")

for (chr in all_chromosomes){
  s = (chrom_number_sub == chr)
  t = which(s == TRUE)
  
  chr_num = chrom_number_sub[s]
  chr_pos = chrom_position_sub[s]
  smoothness = 200 / length(chr_pos)
  
  # depth for WMI
  chr_wmi_depth = depth_WMI[s]
  chr_wmi_loe = loess(chr_wmi_depth~chr_pos, span = smoothness)
  chr_wmi_depth_line = predict(chr_wmi_loe)
  
  # depth for WLI
  chr_wli_depth = depth_WLI[s]
  chr_wli_loe = loess(chr_wli_depth~chr_pos, span = smoothness)
  chr_wli_depth_line = predict(chr_wli_loe)
  
  # Quality for WMI
  #chr_wmi_quality = quality_WLI[s,]
  
  chr_wmi_quality = quality_WMI[s]
  chr_wmi_loe_0 = loess(chr_wmi_quality~chr_pos, span = smoothness)
  #chr_wmi_loe_1 = loess(chr_wmi_quality[,1]~chr_pos, span = smoothness)
  #chr_wmi_loe_2 = loess(chr_wmi_quality[,2]~chr_pos, span = smoothness)
  #chr_wmi_loe_3 = loess(chr_wmi_quality[,3]~chr_pos, span = smoothness)
  
  chr_wmi_quality_line_0 = predict(chr_wmi_loe_0)
  #chr_wmi_quality_line_1 = predict(chr_wmi_loe_1)
  #chr_wmi_quality_line_2 = predict(chr_wmi_loe_2)
  #chr_wmi_quality_line_3 = predict(chr_wmi_loe_3)
  
  # Quality for WLI
  chr_wli_quality = quality_WLI[s]
  chr_wli_loe_0 = loess(chr_wli_quality~chr_pos, span = smoothness)
  # chr_wli_loe_1 = loess(chr_wli_quality[,1]~chr_pos, span = smoothness)
  # chr_wli_loe_2 = loess(chr_wli_quality[,2]~chr_pos, span = smoothness)
  # chr_wli_loe_3 = loess(chr_wli_quality[,3]~chr_pos, span = smoothness)
  
  chr_wli_quality_line_0 = predict(chr_wli_loe_0)
  # chr_wli_quality_line_1 = predict(chr_wli_loe_1)
  # chr_wli_quality_line_2 = predict(chr_wli_loe_2)
  # chr_wli_quality_line_3 = predict(chr_wli_loe_3)
  
  test = cbind(chr_num,chr_pos,
               chr_wmi_depth_line,chr_wli_depth_line,
               chr_wmi_quality_line_0,
               chr_wli_quality_line_0)
  
  # test = cbind(chr_num,chr_pos,
  #              chr_wmi_depth_line,chr_wli_depth_line,
  #              chr_wmi_quality_line_1, chr_wmi_quality_line_2, chr_wmi_quality_line_3,
  #              chr_wli_quality_line_1, chr_wli_quality_line_2, chr_wli_quality_line_3, chr_wmi_quality_line_0, chr_wli_quality_line_0
  # )
  
  table_circos = rbind(table_circos, test)
  
}


################################################################################
## Initialize Circos plots
################################################################################


pdf("circos_plot_new_V6.pdf")


WLI_color = "red"
WMI_color = "blue"

head(table_circos)

table_circos_backup = table_circos

table_circos = table_circos[table_circos$chr_num != "Y",]

table_circos$chr_num = factor(table_circos$chr_num, levels = c('1','2','3','4','5','6','7','8','9','10',
                                                               '11','12','13','14','15','16','17','18',
                                                               '19','20',"X"))

cirle_factors_ALT = droplevels(as.factor(alt_table[,1]),"Y")
cirle_factors_ALT = factor(cirle_factors_ALT, levels = c('1','2','3','4','5','6','7','8','9','10',
                                                         '11','12','13','14','15','16','17','18',
                                                         '19','20',"X"))

cirle_factors_HET = droplevels(as.factor(het_table[,1]),"Y")
cirle_factors_HET = factor(cirle_factors_HET, levels = c('1','2','3','4','5','6','7','8','9','10',
                                                         '11','12','13','14','15','16','17','18',
                                                         '19','20',"X"))

depth_red = as.numeric(as.character(table_circos$chr_wli_depth_line))
depth_red = ifelse(depth_red > 60, 60, depth_red)
depth_blue = as.numeric(as.character(table_circos$chr_wmi_depth_line))
depth_blue = ifelse(depth_blue > 60, 60, depth_blue)


circos.clear()
circos.par("track.height" = 0.1, cell.padding = c(0, 0.4, 0, 0.4), start.degree = 0)
circos.initialize(factors = cirle_factors_ALT, x = alt_table[,2]) # cell.padding	c(0.02, 1.00, 0.02, 1.00))




######################

# #1. Draw the gene names
# circos.track(factors = as.numeric(get_genes_Chromosome), x = get_genes[,3], y = rep (0, length(get_genes_Chromosome) ),ylim = c(0,1),
#              track.height = 0.1, bg.border = NA,
#              panel.fun = function(x, y) {
# 
#              circos.text(x = x, y = y + uy(0, "mm"),
#                            labels = get_genes[,1],facing = "reverse.clockwise", niceFacing = TRUE,
#                            cex = 0.01, adj = c(1, 0.01))
#              })


# 1. Create the depth track: 
circos.track(as.factor(table_circos$chr_num),  ylim = c(19,61),
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(4, "mm"), 
                           CELL_META$sector.index,facing = "clockwise", niceFacing = TRUE, 
                           cex = 0.6, adj = c(1, 0.01))
             }
)


startOffset = 100000
endOffset = -100000

circos.trackLines(factors = as.factor(table_circos$chr_num), x = as.numeric(as.character(table_circos$chr_pos)), y =  depth_red, lwd = 0.4, col = adjustcolor( "red", alpha.f = 0.5))
circos.trackLines(factors = as.factor(table_circos$chr_num), x = as.numeric(as.character(table_circos$chr_pos)), y =  depth_blue, lwd = 0.4, col = adjustcolor( "blue", alpha.f = 0.5))



# 2. draw the snp positions
circos.track(factors = cirle_factors_ALT, x = as.numeric(alt_table[,2]), y = as.numeric(alt_table[,3]),
             track.height = 0.1, ylim = c(-1.05,1.05),
             panel.fun = function(x, y) {
               
               # circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(4, "mm"), 
               #             CELL_META$sector.index,facing = "clockwise", niceFacing = TRUE, 
               #             cex = 0.6, adj = c(1, 0.01))
               
               circos.rect(xleft = x,
                           ybottom = rep(0,length(y)),
                           xright = x+2500000, #5000000 is pretty good
                           ytop = y,
                           col = ifelse( y > 0, adjustcolor( "red", alpha.f = 0.1), adjustcolor( "blue", alpha.f = 0.1)),
                           border = NA)
             })


# 3. Draw the HET snps
circos.track(factors = cirle_factors_HET, x = as.numeric(het_table[,2]), y = as.numeric(het_table[,3]),
             track.height = 0.06, ylim = c(-1.05,1.05),
             panel.fun = function(x, y) {
               
               circos.rect(xleft = x,
                           ybottom = rep(0,length(y)),
                           xright = x+2500000, #5000000 is pretty good
                           ytop = y,
                           col = ifelse( y > 0, adjustcolor( "firebrick1", alpha.f = 0.1),adjustcolor("cornflowerblue", alpha.f = 0.1)),
                           border = NA)
             })


# 4. Quality track underneath
circos.track(as.factor(table_circos$chr_num), ylim = c(8,30))
circos.trackLines(factors = as.factor(table_circos$chr_num), x = as.numeric(as.character(table_circos$chr_pos)), y =  as.numeric(as.character(table_circos$chr_wmi_quality_line_0)), lwd = 0.1, col = adjustcolor( "blue", alpha.f = 0.5))
circos.trackLines(factors = as.factor(table_circos$chr_num), x = as.numeric(as.character(table_circos$chr_pos)), y =  as.numeric(as.character(table_circos$chr_wli_quality_line_0)), lwd = 0.1, col = adjustcolor( "red", alpha.f = 0.5))




dev.off()
circos.clear()




####################################
## Deletion overlap
####################################

# if (!require("BiocManager"))
#   install.packages("BiocManager")
# BiocManager::install("GenomicRanges")


library(GenomicRanges)

########### LOAD IN THE DELETION RANGES
wli_del_data = read.csv("data_deletions/wli_del_table.csv", sep = "\t")
wmi_del_data = read.csv("data_deletions/wmi_del_table.csv", sep = "\t")

wli_ranges = GRanges(seqnames = wli_del_data$X.CHROM, IRanges(wli_del_data$POS, end = wli_del_data$ALT ))
wmi_ranges = GRanges(seqnames = wmi_del_data$X.CHROM, IRanges(wmi_del_data$POS, end = wmi_del_data$ALT ))


########### LOAD IN THE UN-OPPOSED CALLS
WLI_unopposed_split = strsplit(rownames(final_WLI_UNOPPOSED_30),"_")
WLI_unopposed_chr = paste("chr",as.character(sapply(WLI_unopposed_split, "[", 1)),sep = "")
WLI_unopposed_pos = as.numeric(as.character(sapply(WLI_unopposed_split, "[", 2)))

WMI_unopposed_split = strsplit(rownames(final_WMI_UNOPPOSED_30),"_")
WMI_unopposed_chr = paste("chr",as.character(sapply(WMI_unopposed_split, "[", 1)),sep = "")
WMI_unopposed_pos = as.numeric(as.character(sapply(WMI_unopposed_split, "[", 2)))

wli_unopposed_ranges = GRanges(seqnames = WLI_unopposed_chr, IRanges(WLI_unopposed_pos, end = WLI_unopposed_pos+1 ))
wmi_unopposed_ranges = GRanges(seqnames = WMI_unopposed_chr, IRanges(WMI_unopposed_pos, end = WMI_unopposed_pos+1 ))


### Alright, lets see where the ranges overlap.
nrow(wli_del_data)            # Number of deletions in the wmi
nrow(wmi_del_data)            # Number of deletions in the wli
length(WLI_unopposed_pos)     # Number of un-opposed calls in WLI
length(WMI_unopposed_pos)     # Number of un-opposed calls in WMI


## Work with these
wli_ranges
wmi_ranges
wli_unopposed_ranges
wmi_unopposed_ranges

intersect(wmi_unopposed_ranges, wli_ranges)
intersect(wli_unopposed_ranges, wmi_ranges)





