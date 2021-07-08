

# load the SNPs WLI

#install.packages("devtools")
#devtools::install_github("bmansfeld/QTLseqr")
library("QTLseqr")

loadVCF = function(inFile){
    tmp_vcf<-readLines(inFile)
    tmp_vcf_data<-read.table(inFile, stringsAsFactors = FALSE)
    
    # filter for the columns names
    tmp_vcf<-tmp_vcf[-(grep("#CHROM",tmp_vcf)+1):-(length(tmp_vcf))]
    vcf_names<-unlist(strsplit(tmp_vcf[length(tmp_vcf)],"\t"))
    names(tmp_vcf_data)<-vcf_names
    head(tmp_vcf_data)
    ready_VCF = tmp_vcf_data
    return(ready_VCF)
}



# 1. Where are the hotspots?

runWindows = function(inData, windowSize){
    
    #windowSize = 10000
    outFrame = data.frame("chromosome" = "Test", "position" = 1, "start" = 1, "end" = 1, "numberOfSNPs" = 0)
    for (chromosome in c(1:20,"X","Y")){
        just1 = inData[inData$`#CHROM` == chromosome,]
        
        
        
        if(nrow(just1)==0){
          next
        }
        
        for (i in 1:nrow(just1)){
          #print(i)
          position = just1$POS[i]         # a SNP
          start = just1$POS[i] -1         # Search from here
          end = just1$POS[i] + windowSize # until here
          totalSNP = sum(just1$POS > start & just1$POS < end) # how many from all SNPs fall within the window?
          
          #print(just1$POS > start)
          #print(just1$POS < end)
          
          outline = data.frame("chromosome" = chromosome, 
                               "position" = position, 
                               "start" = start, 
                               "end" = end,
                               "numberOfSNPs" = totalSNP )
          outFrame = rbind(outFrame, outline)
        }
        plot(outFrame$start[outFrame$chromosome == chromosome], 
             outFrame$numberOfSNPs[outFrame$chromosome == chromosome], 
             pch = 16,
             main = paste("SNP enrichment Chr", chromosome, "\n", "window = ", windowSize,sep = ""))
        localMaximum = max(outFrame$numberOfSNPs[outFrame$chromosome == chromosome])
        if(localMaximum >= 5){
          print(paste("Notable amount: ",localMaximum," on chr",chromosome,"Window size:", windowSize))
        }
    }
    return(outFrame)
}


########################################################################
# Question: are there any hotspots?

WLI = loadVCF(inFile = "results_final_Snpeff/ALT_WLI_calls_30_ann.vcf")
#WLI = loadVCF(inFile = "results_final_vcf_merged/ALT_WLI_calls_30.vcf")


#hotSpots_WLI = runWindows(inData = WLI, windowSize = 100)
#hotSpots_WLI = runWindows(inData = WLI, windowSize = 1000)
#hotSpots_WLI = runWindows(inData = WLI, windowSize = 5000)
hotSpots_WLI = runWindows(inData = WLI, windowSize = 10000) # 10K
hotSpots_WLI = runWindows(inData = WLI, windowSize = 15000) # 15K
#hotSpots_WLI = runWindows(inData = WLI, windowSize = 25000) # 25K
#hotSpots_WLI = runWindows(inData = WLI, windowSize = 50000) # 50K
#hotSpots_WLI = runWindows(inData = WLI, windowSize = 60000) # 100K
hotSpots_WLI = runWindows(inData = WLI, windowSize = 500000) # 500K
#hotSpots_WLI = runWindows(inData = WLI, windowSize = 1000000) # 1M

hot_WLI = hotSpots_WLI[hotSpots_WLI$numberOfSNPs>=3,]
selection = which(WLI$POS %in% hot_WLI$position)
selection2 = unique(sort(c(selection, selection+1, selection+2)))
clustered_WLI = WLI[selection2,] 

clustered_WLI$location = sapply(strsplit(clustered_WLI$INFO,"|", fixed = TRUE), `[`, 2)
clustered_WLI$type = sapply(strsplit(clustered_WLI$INFO,"|", fixed = TRUE), `[`, 3)
clustered_WLI$nearestGene = sapply(strsplit(clustered_WLI$INFO,"|", fixed = TRUE), `[`, 4)
WLI_out = clustered_WLI[,c(1:6,16,17,18)]
write.csv(WLI_out, "WLI_hotspot_SNPS.csv")



########################################################################
WMI = loadVCF(inFile = "results_final_Snpeff/ALT_WMI_calls_30_ann.vcf")
#WMI = loadVCF(inFile = "results_final_vcf_merged/ALT_WMI_calls_30.vcf")


#hotSpots_WMI = runWindows(inData = WMI, windowSize = 100)
hotSpots_WMI = runWindows(inData = WMI, windowSize = 1000)
#hotSpots_WMI = runWindows(inData = WMI, windowSize = 3000)
#hotSpots_WMI = runWindows(inData = WMI, windowSize = 10000) # 10K
#hotSpots_WMI = runWindows(inData = WMI, windowSize = 25000) # 25K
#hotSpots_WMI = runWindows(inData = WMI, windowSize = 50000) # 50K
#hotSpots_WMI = runWindows(inData = WMI, windowSize = 100000) # 100K
#hotSpots_WMI = runWindows(inData = WMI, windowSize = 500000) # 500K
#hotSpots_WMI = runWindows(inData = WMI, windowSize = 1000000) # 1M

hot_WMI = hotSpots_WMI[hotSpots_WMI$numberOfSNPs>=5,]
selection = which(WMI$POS %in% hot_WMI$position)
selection2 = unique(sort(c(selection, selection+1, selection+2,selection+3,selection+4,selection+5)))
clustered_WMI = WMI[selection2,] 

clustered_WMI$location = sapply(strsplit(clustered_WMI$INFO,"|", fixed = TRUE), `[`, 2)
clustered_WMI$type = sapply(strsplit(clustered_WMI$INFO,"|", fixed = TRUE), `[`, 3)
clustered_WMI$nearestGene = sapply(strsplit(clustered_WMI$INFO,"|", fixed = TRUE), `[`, 4)
WMI_out = clustered_WMI[,c(1:6,16,17,18)]
write.csv(WMI_out, "WMI_hotspot_SNPS.csv")

#Retrieve hotspots


# answer: Few, we can look into more depth after meeting with Hao

########################################################################
# Question: Are there any sex specific SNPs?

WLI_X = WLI[WLI$`#CHROM` == "X",]
WLI_Y = WLI[WLI$`#CHROM` == "Y",]

WLI_X$location = sapply(strsplit(WLI_X$INFO,"|", fixed = TRUE), `[`, 2)
WLI_X$type = sapply(strsplit(WLI_X$INFO,"|", fixed = TRUE), `[`, 3)
WLI_X$nearestGene = sapply(strsplit(WLI_X$INFO,"|", fixed = TRUE), `[`, 4)
WLI_out = WLI_X[,c(1:6,16,17,18)]
write.csv(WLI_out, "WLI_sex_SNPS.csv")

WMI_X = WMI[WMI$`#CHROM` == "X",]
WMI_Y = WMI[WMI$`#CHROM` == "Y",]

WMI_X$location = sapply(strsplit(WMI_X$INFO,"|", fixed = TRUE), `[`, 2)
WMI_X$type = sapply(strsplit(WMI_X$INFO,"|", fixed = TRUE), `[`, 3)
WMI_X$nearestGene = sapply(strsplit(WMI_X$INFO,"|", fixed = TRUE), `[`, 4)
WMI_out = WMI_X[,c(1:6,16,17,18)]
write.csv(WMI_out, "WMI_sex_SNPS.csv")


dim(WLI_X)
table(WLI_X$location)

write.csv(WLI_out[WLI_out$location != "intergenic_region",], "WLI_X_impact.csv")


dim(WMI_X)
table(WMI_X$location)

write.csv(WMI_out[WMI_out$location != "intergenic_region",], "WMI_X_impact.csv")

########################################

#CHECK RNA OVERLAP:

WMI$location = sapply(strsplit(WMI$INFO,"|", fixed = TRUE), `[`, 2)
WMI$type = sapply(strsplit(WMI$INFO,"|", fixed = TRUE), `[`, 3)
WMI$nearestGene = sapply(strsplit(WMI$INFO,"|", fixed = TRUE), `[`, 4)

WMI$nearestGene1 = tolower(sapply(strsplit(WMI$nearestGene,"-", fixed = TRUE), `[`, 1))
WMI$nearestGene2 = tolower(sapply(strsplit(WMI$nearestGene,"-", fixed = TRUE), `[`, 2))


WLI$location = sapply(strsplit(WLI$INFO,"|", fixed = TRUE), `[`, 2)
WLI$type = sapply(strsplit(WLI$INFO,"|", fixed = TRUE), `[`, 3)
WLI$nearestGene = sapply(strsplit(WLI$INFO,"|", fixed = TRUE), `[`, 4)

WLI$nearestGene1 = tolower(sapply(strsplit(WLI$nearestGene,"-", fixed = TRUE), `[`, 1))
WLI$nearestGene2 = tolower(sapply(strsplit(WLI$nearestGene,"-", fixed = TRUE), `[`, 2))


########################################################################
## Microarray
########################################################################
micro_results = read.csv("data_RNA/23_DE_genes_WLIWMI.csv", sep = ";", header = FALSE)

headNames = micro_results[1,]
colnames(micro_results) = headNames
micro_results = micro_results[2:nrow(micro_results),]
micro_results$GeneSymbol = tolower(micro_results$GeneSymbol)

## Checking
table(WLI$nearestGene1 %in% micro_results$GeneSymbol)
table(WLI$nearestGene2 %in% micro_results$GeneSymbol)

table(WMI$nearestGene1 %in% micro_results$GeneSymbol)
table(WMI$nearestGene2 %in% micro_results$GeneSymbol)

WLI_simple = WLI[,c(1:7,16,17,19,20)]
WMI_simple = WMI[,c(1:7,16,17,19,20)]

WLI_simple$match1 =  (WLI_simple$nearestGene1 %in% micro_results$GeneSymbol)
WLI_simple$match2 =  (WLI_simple$nearestGene2 %in% micro_results$GeneSymbol)
WLI_simple$DE_gene = ifelse(WLI_simple$match1, WLI_simple$nearestGene1, WLI_simple$nearestGene2)

WMI_simple$match1 =  (WMI_simple$nearestGene1 %in% micro_results$GeneSymbol)
WMI_simple$match2 =  (WMI_simple$nearestGene2 %in% micro_results$GeneSymbol)
WMI_simple$DE_gene = ifelse(WMI_simple$match1, WMI_simple$nearestGene1, WMI_simple$nearestGene2)

WLI_micro_match = WLI_simple[ (WLI_simple$match1 == TRUE | WLI_simple$match2 == TRUE),]
WMI_micro_match = WMI_simple[ (WMI_simple$match1 == TRUE | WMI_simple$match2 == TRUE),]

##########################################################
## Find and bind
##########################################################
WLI_micro_match$sample = rep("WLI",nrow(WLI_micro_match))
WMI_micro_match$sample = rep("WMI",nrow(WMI_micro_match))

WLI_micro_match$source = rep("PMC3117129",nrow(WLI_micro_match))
WMI_micro_match$source = rep("PMC3117129",nrow(WMI_micro_match))

WLI_micro_match$tissue = rep("Amygdala",nrow(WLI_micro_match))
WMI_micro_match$tissue = rep("Amygdala",nrow(WMI_micro_match))

WLI_micro_match$method = rep("microarray",nrow(WLI_micro_match))
WMI_micro_match$method = rep("microarray",nrow(WMI_micro_match))

micro_match = rbind(WLI_micro_match,WMI_micro_match)

table(WLI_micro_match$location)
table(WMI_micro_match$location)

########################################################################
## RNA
########################################################################
dim(WLI_RNA_match)
dim(WMI_RNA_match)

head(WMI_RNA_match)

head(WMI_RNA_match)

## RNA SEQ META ANALYSIS:
RNA_sham_results = read.csv("data_RNA/RNA_seq_shamVsham.csv", sep = ";", header = FALSE)


WLI_simple$match1 =  (WLI_simple$nearestGene1 %in% tolower(RNA_sham_results$V1))
WLI_simple$match2 =  (WLI_simple$nearestGene2 %in% tolower(RNA_sham_results$V1))
WLI_simple$DE_gene = ifelse(WLI_simple$match1, WLI_simple$nearestGene1, WLI_simple$nearestGene2)

WMI_simple$match1 =  (WMI_simple$nearestGene1 %in% RNA_results$GeneSymbol)
WMI_simple$match2 =  (WMI_simple$nearestGene2 %in% RNA_results$GeneSymbol)
WMI_simple$DE_gene = ifelse(WMI_simple$match1, WMI_simple$nearestGene1, WMI_simple$nearestGene2)

WLI_RNA_match = WLI_simple[ (WLI_simple$match1 == TRUE | WLI_simple$match2 == TRUE),]
WMI_RNA_match = WMI_simple[ (WMI_simple$match1 == TRUE | WMI_simple$match2 == TRUE),]


##########################################################
## Find and bind
##########################################################
WLI_RNA_match$sample = rep("WLI",nrow(WLI_RNA_match))
WMI_RNA_match$sample = rep("WMI",nrow(WMI_RNA_match))

WLI_RNA_match$source = rep("PMC5786888",nrow(WLI_RNA_match))
WMI_RNA_match$source = rep("PMC5786888",nrow(WMI_RNA_match))

WLI_RNA_match$tissue = rep("Hippocampus",nrow(WLI_RNA_match))
WMI_RNA_match$tissue = rep("Hippocampus",nrow(WMI_RNA_match))

WLI_RNA_match$method = rep("RNA-seq",nrow(WLI_RNA_match))
WMI_RNA_match$method = rep("RNA-seq",nrow(WMI_RNA_match))

RNA_match = rbind(WLI_RNA_match,WMI_RNA_match)


outTable = rbind(micro_match, RNA_match)

dim(outTable) # There are 232 SNPs

length(unique(outTable$DE_gene)) # in proximity of 101 genes 
                                 # which were significantly differentially expressed in previous research.
                        
table(outTable$location) # 128 in intergenic regions, 95 intron variants, 
                         # 4 upstream gene variants, 3 downstream of genes, 2 in 3'prime UTR.

outTable$tissue

table(outTable$source)


write.csv(outTable, "Z_additional_S_table.csv")

table(WLI_micro_match$location)
#table(WMI_micro_match$location)


# Merge results tables.


