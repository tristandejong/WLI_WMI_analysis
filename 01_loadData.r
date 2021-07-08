
chr_list = c(1:20,"X", "Y")

#chr_list = c("X","Y")

rowMax = function(df){
  out = apply(df, 1, FUN=max)
  return(out)
}


first = TRUE

grandTotal = 0
rawTotal = 0

for (chr in chr_list){
    #pdf(paste("figures_chr",chr,".pdf",sep = ""),paper='A4r')#, width = 16, height = 5)
  
    print(paste("current chromosome",chr))
    ###############################################################################################
    ## genotype
    ###############################################################################################
    print("1. Reading file") 
    readfolder = ("filtered/")
    readFile = paste("deepvariant_wmiwli_chr",chr,"_6_samples.gvcf.gt",sep = "")
    
    gq = 0
    readFile3 = paste("deepvariant_wmiwli_chr",chr,"_6_samples.gvcf.gq",sep = "")
    gq = read.table(paste(readfolder,readFile3, sep = ""), sep = "\t", stringsAsFactors = TRUE, row.names = 1, header = T)
    gq = gq[,c("WLI_IonProton", "WLI_chromium", "WLI_xTen", "WMI_IonProton", "WMI_chromium", "WMI_xTen")]
    
    rawTotal = rawTotal + nrow(gq)
    
    gt = 0
    gt = read.table(paste(readfolder,readFile, sep = ""), sep = "\t", stringsAsFactors = TRUE, row.names = 1, header = T)
    gt = gt[,c("WLI_IonProton", "WLI_chromium", "WLI_xTen", "WMI_IonProton", "WMI_chromium", "WMI_xTen")]
    
    dp = 0
    readFile2 = paste("deepvariant_wmiwli_chr",chr,"_6_samples.gvcf.dp",sep = "")
    dp = read.table(paste(readfolder,readFile2, sep = ""), sep = "\t", stringsAsFactors = TRUE, row.names = 1, header = T)
    dp = dp[,c("WLI_IonProton", "WLI_chromium", "WLI_xTen", "WMI_IonProton", "WMI_chromium", "WMI_xTen")]
    
    
    gq_filter = rowMax(gq) >= 10
    gt_filter1 = rowSums(gt == "1/1") <= 5 
    gt_filter2 = rowSums(gt == "0/0") <= 5
    
    final_filter = gq_filter & gt_filter1 & gt_filter2
    
    gq = gq[final_filter,]
    gt = gt[final_filter,]
    dp = dp[final_filter,]

    #REF = ifelse(gt == "0/0" | gt == "0|0",1,0)
    #HET = ifelse(gt == "0/1" | gt == "0|1",1,0)
    #ALT = ifelse(gt == "1/1" | gt == "1|1",1,0)
    
    #ALL = ifelse(gt == "0/1" | gt == "0|1",1,0)
    #ALL = ifelse(gt == "1/1" | gt == "1|1",2, ALL)
    
    #REF = REF[rowSums(REF) != 0,] 
    #HET = HET[rowSums(HET) != 0,] 
    #ALT = ALT[rowSums(ALT) != 0,] 
    #ALL = ALL[rowSums(ALL) != 0,]
    
    grandTotal = grandTotal + nrow(gt)
    
    #write.csv(REF, file = paste("filtered/", readFile,".REF.filtered.csv", sep = ""))
    #write.csv(HET, file = paste("filtered/", readFile,".HET.filtered.csv", sep = ""))
    #write.csv(ALT, file = paste("filtered/", readFile,".ALT.filtered.csv", sep = ""))
    #write.csv(ALL, file = paste("filtered/", readFile,".ALL.filtered.csv", sep = ""))
    write.csv(gt, file = paste("filtered/", readFile,".filtered.csv", sep = ""))
    write.csv(dp, file = paste("filtered/", readFile2,".filtered.csv", sep = ""))
    write.csv(gq, file = paste("filtered/", readFile3,".filtered.csv", sep = ""))

}


print(paste("before filtering:", rawTotal))
print(paste("after filtering:", grandTotal))
print(paste("fraction of total:", grandTotal /rawTotal))   



