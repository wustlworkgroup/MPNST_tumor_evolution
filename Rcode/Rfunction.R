library(karyoploteR)
library(RColorBrewer)
library(GenomicRanges)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(dplyr)
load("/Volumes/lyu.yang/MAC_1/R/Public_database/Reference/Chrosmome_cytobank.rdata") # define the path for your reference


Gene_annotation_for_CNV=function(Data) {
  
  edb=EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
  
  library(ensembldb)
  #data_list=list()
  Data$Gene_list=""
  for (i in 1:nrow(Data)){
    nchr=Data[i,]$Chromosome
    
    chr=gsub("chr","" ,nchr)
    Start=Data[i,]$Start
    End=Data[i,]$End
    
    grf =GRangesFilter(GenomicRanges::GRanges(chr, ranges =IRanges:: IRanges(Start,End)), type = "any")
    ## Query genes:
    gn <- genes(edb, filter = grf)
    Genelist=gn%>%as.data.frame() %>%dplyr::filter(gene_biotype=="protein_coding") %>%.$gene_name
    
    #return(Genelist)
    # return(Data[i,"Gene_list"])
    if (length(Genelist>0)) {
      Data[i,"Gene_list"]=toString (Genelist)
      #data_list[[i]]=Genelist
    }
  }
  
  return(Data)
  
}


Gene_annotation_for_CNV_2=function(Data) {
  
  edb=EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
  
  library(ensembldb)
  #data_list=list()
  Data$Gene_list=""
  for (i in 1:nrow(Data)){
    nchr=Data[i,]$Chromosome
    
    chr=gsub("chr","" ,nchr)
    Start=Data[i,]$Start
    End=Data[i,]$End
    
    grf =GRangesFilter(GenomicRanges::GRanges(chr, ranges =IRanges:: IRanges(Start,End)), type = "any")
    ## Query genes:
    gn <- genes(edb, filter = grf)
    Genelist=gn%>%as.data.frame() %>%dplyr::filter(gene_biotype=="protein_coding") %>%.$gene_name
    
    #return(Genelist)
    # return(Data[i,"Gene_list"])
    if (length(Genelist>0)) {
      Data[i,"Gene_list"]=toString (Genelist)%>%gsub(" ","",.)
      #data_list[[i]]=Genelist
    }
  }
  
  return(Data)
  
}


Normal_gene_WEX_CN_2=function (WEX_germline_path) {
  # Internal_bin=1000
  # loading normal germline data (chr17)
  
  #Sample=list.files(WEX_germline_path,recursive = F) %>% .[grepl(Patient_ID,.)]
  #File_path=file.path(WEX_germline_path,Sample)
  seg_file=list.files(WEX_germline_path,pattern = "_segments.tsv",full.names = T) 
  SEG_data=readr::read_tsv(seg_file,show_col_types = FALSE)%>%tidyr::separate(ID,c("Type","Chromosome","Start","End"),sep="_")
  SEG_data$Start=as.numeric(SEG_data$Start)
  SEG_data$End=as.numeric(SEG_data$End)
  # ALT_seg=SEG_data%>%dplyr::filter(ALT!=".")
  
  Gene_name=Gene_annotation_for_CNV(SEG_data)
  colnames( Gene_name)[grepl("NP",colnames( Gene_name))]="NP"
  
  Gene_name2=Gene_name%>%dplyr::filter(NP>100)
  # fILTER BY NP
  
  return(Gene_name2)
  #Interval_file=file.path(File_path,paste0(Sample,"_intervals.tsv") )  
  #Interval_file=list.files(WEX_germline_path,pattern = "_intervals.tsv",full.names = T) 
  
}



 

BAF_densMode <- function(Data,nclust,minBAF=0,maxBAF=0.5,title){
  library(maftools)
  #td <- density(Data[,"BAF"])
  
  # maxDens <- which.max(td$y)
  # output=data.frame(x=td$x[maxDens], y=td$y[maxDens])
  #  p=ggplot(Data, aes_string(x)) + geom_density() + geom_vline(xintercept =output[[1]])+ggtitle(title) 
  
  Data = Data%>%dplyr::filter(BAF > minBAF&BAF < maxBAF)
  
  tsb.cluster = mclust::densityMclust(Data[,"BAF"], G = 1:nclust, verbose = FALSE,plot = T)
  
  Data$cluster = as.character(tsb.cluster$classification)
  abs.med.dev = abs(Data[, "BAF"] - median(Data[, "BAF"]))
  pat.mad = median(abs.med.dev) * 100
  pat.math = pat.mad * 1.4826/median((Data[, "BAF"]))
  Data$MATH = pat.math
  Data$MedianAbsoluteDeviation = pat.mad     
  # return(Data)
  # Data = maftools:::refineClusters(clusters = Data)
  
  #clust.dat = rbind(clust.dat, Data, fill = TRUE)
  
  # clust.dat.mean = clust.dat[, mean(x), by = .(cluster)]
  clust.dat.mean =Data%>%group_by(cluster) %>%summarise(meanBAF=mean(BAF))
  return(clust.dat.mean)
  
  
  # plot(density (DN_reads$RD))
  # DN_reads_mean=density(DN_reads$RD)$x[which.max(density(DN_reads$RD)$y)]
  # ggplot(DN_reads, aes(RD)) + geom_density() + geom_vline(xintercept = DN_reads_mean)
}

Tumor_purity_cal_by_BAF=function(BAF_value,Variant_type) {
  
  
  TP="Not determined"
  
  
  
  if (Variant_type=="LOH") {
    
    #TP=(2*BAF_value)-1
    TP=1-2*BAF_value%>%round(.,2)
    
  }
  
  if (Variant_type=="LOH+Trisomy") {
    
    #TP=(2*BAF_value)-1
    TP=(1-2*BAF_value)/(BAF_value+1)%>%round(.,2)
    
  }
  
  
  if(Variant_type=="Het_A_loss") {
    
    TP=2-(1/BAF_value)%>%round(2,.)
    
  }
  
  if(Variant_type=="Het_B_loss") {
    
    TP=(1-(2*BAF_value))/(1-BAF_value)%>%round(.,2)
    
  }
  
  # if(Variant_type=="Trisomy"){
  #  # TP=(2*RDR)-2
  #   TP=(1/BAF_value)-2
  # }
  
  # if(Variant_type=="Trisomy"){
  #   # B gain
  #   TP=(2*BAF_value-1)/(2-BAF_value)%>%round(.,2)
  # }  
  
  if(Variant_type=="Trisomy"){
    # B gain
    TP=(1/BAF_value)-2%>%round(.,2) # TP=0.5,BAF=0.4
  }
  
  # if((Variant_type=="Tetrasomy")&(BAF_value<0.45)){
  #   BAF=1-BAF_value
  #   # B gain
  #   TP1=(BAF-0.45)/(1-BAF)
  #   TP2=(0.5/BAF_value)-1
  #   TP=paste0("LOH_B:",round(TP1,2),",","LOH_A:",round(TP2,2)) 
  #   }
  
  
  
  return(TP)
  
}







group_numbers <- function(numbers, threshold = 0.1,start_with=1) {
  # Sort numbers to simplify grouping
  sorted_indices <- order(numbers)
  sorted_numbers <- numbers[sorted_indices]
  
  # Initialize the first group
  group_ids <- numeric(length(sorted_numbers))
  current_group <- start_with
  group_ids[1] <- current_group
  
  for (i in 2:length(sorted_numbers)) {
    if (sorted_numbers[i] - sorted_numbers[i - 1] < threshold) {
      # Add to the current group if the difference is less than the threshold
      group_ids[i] <- current_group
    } else {
      # Increment group number and start a new group
      current_group <- current_group + 1
      group_ids[i] <- current_group
    }
  }
  
  # Create a dataframe with original numbers and their group numbers
  result <- data.frame(
    Number = numbers,
    Group = group_ids[order(sorted_indices)]
  )
  
  return(result)
}

group_numbers_with_clonality <- function(data, clonality, threshold = 0.05) {
  # Sort the data
  #sorted_data <- sort(data)
  sorted_data=data
  # Initialize the group vector
  groups <- numeric(length(sorted_data))
  
  for (i in 1:length(sorted_data)) {
    # Find the closest group from clonality data based on the threshold
    if(!is.na(sorted_data[i])) {
      differences <- abs(clonality$Value - sorted_data[i])
      min_diff <- min(differences)
      if (min_diff < threshold) {
        closest_group <- clonality$Group[which.min(differences)]
        groups[i] <- closest_group
      }  else {
        # Assign to "other" group if no close group is found
        groups[i] <- "other"
      }
    } else { groups[i] <- "other"}
  }
  return(as.character(groups))
  # Create a dataframe with the original numbers and their group numbers
  result <- data.frame(Group = groups[order(match(data, sorted_data))],
                       Value = data
                       
  )
  
  return(result)
}


calculate_tumor_purity_by_BAF_LOH <- function(BAF_LOH) {
  # Ensure BAF_LOH is less than 0.5 and not equal to 0
  if (BAF_LOH >= 0.5 || BAF_LOH == 0) {
    print("BAF_LOH should be less than 0.5 and not equal to 0")
    tumor_purity=NA
  }
  
  # Calculate tumor purity based on BAF_LOH using the provided equation
  tumor_purity <- (1 - BAF_LOH * 2)
  
  return(tumor_purity)
}

calculate_tumor_purity_hetloss_by_RDR <- function(RDR) {
  # Function to calculate tumor purity for heterozygous loss based on RDR
  
  # Calculate tumor purity based on RDR using the equation
  tumor_purity <- (1-RDR)/0.5
  
  return(tumor_purity)
}

calculate_tumor_purity_tetrasomy_by_RDR <- function(rdr_tet) {
  tumor_purity <- rdr_tet - 1
  return(tumor_purity)
}

# Function to calculate tumor purity by RDR for trisomy
calculate_tumor_purity_trisomy_by_RDR <- function(rdr_tris) {
  tumor_purity <- 2 * (rdr_tris - 1)
  return(tumor_purity)
}

# Function to calculate tumor purity by BAF for trisomy
calculate_tumor_purity_trisomy_by_BAF <- function(BAF) {
  tumor_purity <- (2 * BAF - 1) / (BAF - 1)
  return(tumor_purity)
}

# Function to calculate tumor purity by BAF for tetrasomy
calculate_tumor_purity_tetrasomy_by_BAF <- function(BAF) {
  tumor_purity <- (1 - 2 * BAF) / (2 * BAF + 1)
  return(tumor_purity)
}

# Function to calculate tumor purity by BAF_LOH for tetrasomy with 3B
calculate_tumor_purity_tetrasomy_by_BAF_LOH_3B <- function(BAF_LOH) {
  if (BAF_LOH >= 0.5 || BAF_LOH == 0) {
    stop("BAF_LOH should be less than 0.5 and not equal to 0")
  }
  tumor_purity <- (1 / (2 * BAF_LOH)) - 1
  return(tumor_purity)
}

# Function to calculate tumor purity by BAF_LOH for tetrasomy with 4B
calculate_tumor_purity_tetrasomy_by_BAF_LOH_4B <- function(BAF_LOH) {
  if (BAF_LOH >= 0.5 || BAF_LOH == 0) {
    stop("BAF_LOH should be less than 0.5 and not equal to 0")
  }
  tumor_purity <- (1 - 2 * BAF_LOH) / (2 * BAF_LOH + 1)
  return(tumor_purity)
}

# Function to calculate tumor purity by BAF_LOH for trisomy
calculate_tumor_purity_trisomy_by_BAF_LOH <- function(BAF_LOH) {
  # Calculate tumor purity based on BAF_LOH using the provided equation
  tumor_purity <- (1 - BAF_LOH * 2) / BAF_LOH
  return(tumor_purity)
}


# Function to calculate tumor purity for different scenarios
calculate_tumor_purity <- function(RDR, BAF) {
  
  # Function to calculate tumor purity by RDR for tetrasomy
  
  
  
  # Calculate tumor purity for different scenarios
  results <- list()
  
  # Group by Trisomy
  results$Trisomy <- list(
    # Subgroup by Tumor_purity_by_RDR
    by_RDR = calculate_tumor_purity_trisomy_by_RDR(RDR),
    # Subgroup by Tumor_purity_by_BAF_no_LOH
    BY_BAF_no_LOH = calculate_tumor_purity_trisomy_by_BAF(BAF),
    # Subgroup by Tumor_purity_by_BAF_WITH_LOH
    BY_BAF_WITH_LOH = calculate_tumor_purity_trisomy_by_BAF_LOH(BAF)
  )
  
  # Group by Tetrasomy
  results$Tetrasomy <- list(
    # Subgroup by Tumor_purity_by_RDR
    by_RDR = calculate_tumor_purity_tetrasomy_by_RDR(RDR),
    # Subgroup by Tumor_purity_by_BAF_WITH_LOH
    BY_BAF_WITH_LOH_3B = tryCatch({
      calculate_tumor_purity_tetrasomy_by_BAF_LOH_3B(BAF)
    }, error = function(e) {
      NA
    }),
    BY_BAF_WITH_LOH_4B = tryCatch({
      calculate_tumor_purity_tetrasomy_by_BAF_LOH_4B(BAF)
    }, error = function(e) {
      NA
    })
  )
  
  # Convert results to a dataframe
  if (BAF<0.4) {
    results_df <- data.frame(
      Type = rep(c("Trisomy", "Tetrasomy"), each = 3),
      Subtype = c("by_RDR", "BY_BAF_no_LOH", "BY_BAF_WITH_LOH", "by_RDR", "BY_BAF_WITH_LOH_3B", "BY_BAF_WITH_LOH_4B"),
      Tumor_Purity = unlist(results)
    )} else {
      
      results_df <-data.frame(Type = rep(c( "Tetrasomy"), 3),Subtype = c( "by_RDR", "BY_BAF_WITH_LOH_3B", "BY_BAF_WITH_LOH_4B"),Tumor_Purity = unlist(results$Tetrasomy) )
    }
  
  if (RDR>1.5) {
    results_df=results_df%>%dplyr::filter(Type!="Trisomy")
    
  } 
  
  return(results_df)
}


# Load necessary libraries
# if (!requireNamespace("karyoploteR", quietly = TRUE)) {
#   stop("The 'karyoploteR' package is required but not installed.")
# }
# if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
#   stop("The 'GenomicRanges' package is required but not installed.")
# }
# if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
#   stop("The 'RColorBrewer' package is required but not installed.")
# }


Somatic_CN_plot=function(Somatic_SV,Germline_SV) {
  MAX_CN=max(Somatic_SV$CN)
  kp <- karyoploteR::plotKaryotype(genome="hg38",chromosomes =paste0("chr", 1:22),plot.type = 4 )
  kpDataBackground(kp, data.panel = 1, col="grey90")
  kpPlotRegions(kp,arm_q,col="white",r0 = 0, r1 =1)
  
  # kpDataBackground(kp)
  #Default axis
  
  #Axis on the right side of the data.panel
  kpAxis(kp, side = 1,labels = NA,numticks = 1)
  kpAxis(kp, r1=1, ymin=0, ymax =1, numticks = MAX_CN+1, cex=2, data.panel=1,tick.len = 1,labels=c(0:MAX_CN) ) 
  kpAbline(kp, h=c(1), col="gray10", ymin=-1, ymax=1, r0=0.5, r1=2/MAX_CN,lty=2)
  kpAbline(kp, h=c(1), col="gray10", ymin=-1, ymax=1, r0=0.5, r1=4/MAX_CN,lty=2)
  kpAbline(kp, h=c(1), col="gray10", ymin=-1, ymax=1, r0=0.5, r1=1/MAX_CN,lty=2)
  kpAbline(kp, h=c(1), col="gray10", ymin=-1, ymax=1, r0=0.5, r1=3/MAX_CN,lty=2)
  for ( i in 1:nrow(Somatic_SV)) {
    
    Somatic_SV_select=Somatic_SV[i,]
    Somatic_SV_torange=Somatic_SV_select[,c("chr","start","end")]%>%toGRanges()
    
    kpPlotRegions(kp,Somatic_SV_torange,col=Somatic_SV_select$color,r0 =Somatic_SV_select$r0, r1 = Somatic_SV_select$r1)
    kpText(kp, data.panel = 1,chr=Somatic_SV_select$chr, x=Somatic_SV_select$start, y=Somatic_SV_select$r1, labels=Somatic_SV_select$symbol,col=Somatic_SV_select$color, pos=1,cex=1)
    
  }
  
  for ( i in 1:nrow(Germline_SV)) {
    
    Germline_SV_select=Germline_SV[i,]
    Germline_SV_torange=Germline_SV_select[,c("chr","start","end")]%>%toGRanges()
    
    kpPlotRegions(kp,Germline_SV_torange,col=Germline_SV_select$color,r0 =Germline_SV_select$r0, r1 = Germline_SV_select$r1)
    kpText(kp, data.panel = 1,chr=Germline_SV_select$chr, x=Germline_SV_select$start, y=Germline_SV_select$r1, labels=Germline_SV_select$symbol,col=Germline_SV_select$color, pos=2)
    
  }
}



Somatic_CN_plot_2 <- function(Somatic_SV, Germline_SV = NULL,MAX_CN=4,Title) {
  
  # Create a color palette for different clonals
  unique_clonal <- as.integer(unique(Somatic_SV$Value)*100)
  
  num_colors <- length(unique_clonal)
  
  dark2_colors=colorRampPalette(c("white", "black"))(100)
  
  # Map CNV2 values to colors
  clonal_colors <- setNames(dark2_colors, 1:100)
  
  Somatic_SV$line_color <- clonal_colors[as.integer(Somatic_SV$Value*100)]
  Somatic_SV= Somatic_SV%>%dplyr::rename(chr=Chromosome, start=Start, end=End)
  
  # Determine the maximum copy number for axis scaling
  #MAX_CN <- max(Somatic_SV$CN, na.rm = TRUE)
  # Set up the karyotype plot
  kp <- plotKaryotype(genome = "hg38", chromosomes = paste0("chr", 1:22), plot.type = 4)
  kpDataBackground(kp, data.panel = 1, col = "grey90")
  kpPlotRegions(kp, arm_q, col = "white", r0 = 0, r1 = 1)
  
  # Add axis to the plot
  kpAxis(kp, side = 1,labels = NA,numticks = 1)
  kpAxis(kp, r1 = 1, ymin = 0, ymax = 1, numticks = MAX_CN + 1, cex = 2, data.panel = 1, tick.len = 1, labels = 0:MAX_CN)
  
  # Add horizontal lines
  for (i in 1:4) {
    kpAbline(kp, h = 1, col = "gray10", ymin = -1, ymax = 1, r0 = 0.5, r1 = i / MAX_CN, lty = 2)
  }
  
  # Plot somatic structural variants with different colors for each clonal
  for (i in seq_len(nrow(Somatic_SV))) {
    Somatic_SV_select <- Somatic_SV[i, ]
    
    Somatic_SV_torange <- Somatic_SV_select[, c("chr", "start", "end")] %>% toGRanges()
    
    if (Somatic_SV_select$CNV == "LOH") {
      kpPlotRegions(kp, Somatic_SV_torange, col= NA, border = Somatic_SV_select$line_color,r0 = (Somatic_SV_select$CN / MAX_CN) - 0.05, r1 = Somatic_SV_select$CN / MAX_CN, density = 20, angle = 45)
      
      
    } else if (grepl(",CDKN2A," ,Somatic_SV_select$Gene_list)&grepl("loss",Somatic_SV_select$CNV) ) {
      kpPlotRegions(kp, Somatic_SV_torange, col = Somatic_SV_select$line_color, border = Somatic_SV_select$line_color, r0 = (Somatic_SV_select$CN / MAX_CN) - 0.05, r1 = Somatic_SV_select$CN / MAX_CN)
      kpText(kp, data.panel = 1,col="red" ,chr =Somatic_SV_select$chr, x = Somatic_SV_select$start, y = Somatic_SV_select$CN/ MAX_CN+0.1, labels = paste0("CDKN2A_",Somatic_SV_select$CNV), pos = 2)
    } 
    else if (grepl("NF1" ,Somatic_SV_select$Gene_list)&grepl("loss",Somatic_SV_select$CNV)&is.null(Germline_SV) ) {
      kpPlotRegions(kp, Somatic_SV_torange, col = Somatic_SV_select$line_color, border = Somatic_SV_select$line_color, r0 = (Somatic_SV_select$CN / MAX_CN) - 0.05, r1 = Somatic_SV_select$CN / MAX_CN)
      kpText(kp, data.panel = 1,col="red" ,chr =Somatic_SV_select$chr, x = Somatic_SV_select$start, y = Somatic_SV_select$CN/ MAX_CN+0.1, labels = paste0("NF1_",Somatic_SV_select$CNV), pos = 2)
    } 
    
    else {
      kpPlotRegions(kp, Somatic_SV_torange, col = Somatic_SV_select$line_color, border = Somatic_SV_select$line_color, r0 = (Somatic_SV_select$CN / MAX_CN) - 0.05, r1 = Somatic_SV_select$CN / MAX_CN)
    }
    
    #    kpPlotRegions(kp, Somatic_SV_torange, col = Somatic_SV_select$line_color, border = Somatic_SV_select$line_color, r0 =  (Somatic_SV_select$CN/MAX_CN)-0.05, r1 = Somatic_SV_select$CN/MAX_CN)
    #  kpText(kp, data.panel = 1, chr = Somatic_SV_select$chr, x = Somatic_SV_select$start, y = Somatic_SV_select$r1, labels = Somatic_SV_select$symbol, col = Somatic_SV_select$line_color, pos = 1, cex = 1)
  }
  
  # Plot germline structural variants
  if (!is.null(Germline_SV)&nrow(Germline_SV)>0) {
    for (i in seq_len(nrow(Germline_SV))) {
      Germline_SV_select <- Germline_SV[i, ]
      Germline_SV_torange <- Germline_SV_select[, c("chr", "start", "end")] %>% toGRanges()
      kpPlotRegions(kp, Germline_SV_torange, col= "green", border = "green",r0 = (Germline_SV_select$cn / MAX_CN) - 0.05, r1 = Germline_SV_select$cn / MAX_CN, density = 60, angle = -45)
      kpText(kp, data.panel = 1, chr = Germline_SV_select$chr, x = Germline_SV_select$start, y = Germline_SV_select$cn/MAX_CN, labels = paste0("Germline_" ,Germline_SV_select$name), pos = 2)
      
    }  
  }
  
  Group=Somatic_SV[,c("Value","Group2")]%>%dplyr::filter(!grepl("Normal|other",Group2))%>%arrange(desc(Value))%>%unique()%>%mutate(Value=as.integer(Value*100))
  
  
  title(main = Title, cex.main = 1.5, line = 2)
  # Add y-axis label using base R graphics
  mtext("Copy number", side = 2, line = 2.5, cex = 1.2) # Adjust line and cex as needed
  
  # Set up layout for legends
  par(mfrow = c(2, 1), mar = c(2, 6, 3, 2), oma = c(0, 0, 2, 0)) # Adjust margins for layout
  
  # Plot empty frame to create space for legends
  plot.new()
  # legend("bottom", legend = legend_labels, fill = legend_colors, title = "Subclonal_proportion", cex = 0.8, border = "black", ncol = length(legend_labels))
  legend_labels <-paste0("Clonal_", unique(Group$Group2))
  legend_colors <- clonal_colors [as.character (Group$Value)]
  legend("topleft", legend = legend_labels, fill = legend_colors, title = "Subclonal: Proportion", cex = 0.4)  
  # LOH texture legend
  #legend("bottom", inset = c(0, 0), legend = "LOH", fill = NA, density = 20, angle = 45, title = "LOH", cex = 0.8, border = "black")
  legend("topleft", inset=c(0.1,0.1),  legend = c("LOH","Normal"), fill = c(NA,"black"), density = c(20,100), angle =c(45,0), title = "Normal and LOH", cex = 0.4, border = "black")
  
  
}

Somatic_seg_Raw_RDR_plot_2 <- function(Somatic_SV, Germline_SV = NULL) {
  load("/Volumes/lyu.yang/MAC_1/R/Public_database/Reference/Chrosmome_cytobank.rdata")
  
  Somatic_SV= Somatic_SV%>%dplyr::rename(chr=Chromosome, start=Start, end=End)
  MAX_RDR <- max(Somatic_SV$RDR, na.rm = TRUE)%>%ceiling()
  if (MAX_RDR<2) {
    MAX_RDR=2
  }
  # Set up the karyotype plot
  kp <- plotKaryotype(genome = "hg38", chromosomes = paste0("chr", 1:22), plot.type = 4)
  kpDataBackground(kp, data.panel = 1, col = "grey90")
  kpPlotRegions(kp, arm_q, col = "white", r0 = 0, r1 = 1)
  
  # Add axis to the plot
  
  # kpAxis(kp, side = 1,labels = NA,numticks = 1)
  # kpAxis(kp, r1 = 1, ymin = 0, ymax = 1, numticks = MAX_RDR  , cex = 2, data.panel = 1, tick.len = 1,labels = NA)
  
  y_ticks <- seq(0, MAX_RDR, by = 0.5)
  
  y_labels <- as.character(y_ticks)  # Display labels for 0, 0.5, 1, and MAX_RDR
  
  # Add y-axis to the plot
  kpAxis(kp, side = 1, labels = NA, numticks = 1)  # Side 1 may not need labels if it's redundant or for stylistic choices
  kpAxis(kp, r1 = 1, ymin = 0, ymax = 1, numticks = MAX_RDR*2+1,cex = 2, data.panel = 1, tick.len = 1,labels = y_ticks)
  
  
  # Add horizontal lines
  for (i in c(0.5,1.5,2)) {
    kpAbline(kp, h = 1, col = "gray10", ymin = -1, ymax = 1, r0 = 0.5, r1 = i / MAX_RDR, lty = 2)
    
  }
  kpAbline(kp, h = 1, col = "blue", ymin = -1, ymax = 1, r0 = 0.5, r1 = 1 / MAX_RDR, lty = 2,lwd=2)
  
  
  # Plot somatic structural variants with different colors for each clonal
  for (i in seq_len(nrow(Somatic_SV))) {
    Somatic_SV_select <- Somatic_SV[i, ]
    
    Somatic_SV_torange <- Somatic_SV_select[, c("chr", "start", "end")] %>% toGRanges()
    kpPlotRegions(kp, Somatic_SV_torange, col= "black", border =NA,r0 = (Somatic_SV_select$RDR / MAX_RDR) - 0.01, r1 = Somatic_SV_select$RDR/ MAX_RDR)
    
    
    
    #    kpPlotRegions(kp, Somatic_SV_torange, col = Somatic_SV_select$line_color, border = Somatic_SV_select$line_color, r0 =  (Somatic_SV_select$CN/MAX_CN)-0.05, r1 = Somatic_SV_select$CN/MAX_CN)
    #  kpText(kp, data.panel = 1, chr = Somatic_SV_select$chr, x = Somatic_SV_select$start, y = Somatic_SV_select$r1, labels = Somatic_SV_select$symbol, col = Somatic_SV_select$line_color, pos = 1, cex = 1)
  }
  
  # Plot germline structural variants
  if (!is.null(Germline_SV)) {
    for (i in seq_len(nrow(Germline_SV))) {
      Germline_SV_select <- Germline_SV[i, ]
      Germline_SV_torange <- Germline_SV_select[, c("chr", "start", "end")] %>% toGRanges()
      
      kpPlotRegions(kp, Germline_SV_torange, col = Germline_SV_select$color, r0 = Germline_SV_select$r0, r1 = Germline_SV_select$r1)
      kpText(kp, data.panel = 1, chr = Germline_SV_select$chr, x = Germline_SV_select$start, y = Germline_SV_select$r1, labels = Germline_SV_select$symbol, col = Germline_SV_select$color, pos = 2)
    }
  }
  
  
  title(main = "Somatic RDR Plot (before normalization)", cex.main = 1.5, line = 2)
  # Add y-axis label using base R graphics
  mtext("Read depth ratio", side = 2, line = 2.5, cex = 1.2) # Adjust line and cex as needed
}





Somatic_CN_plot=function(Somatic_SV,Germline_SV) {
  MAX_CN=max(Somatic_SV$CN)
  kp <- karyoploteR::plotKaryotype(genome="hg38",chromosomes =paste0("chr", 1:22),plot.type = 4 )
  kpDataBackground(kp, data.panel = 1, col="grey90")
  kpPlotRegions(kp,arm_q,col="white",r0 = 0, r1 =1)
  
  # kpDataBackground(kp)
  #Default axis
  
  #Axis on the right side of the data.panel
  kpAxis(kp, side = 1,labels = NA,numticks = 1)
  kpAxis(kp, r1=1, ymin=0, ymax =1, numticks = MAX_CN+1, cex=2, data.panel=1,tick.len = 1,labels=c(0:MAX_CN) ) 
  kpAbline(kp, h=c(1), col="gray10", ymin=-1, ymax=1, r0=0.5, r1=2/MAX_CN,lty=2)
  kpAbline(kp, h=c(1), col="gray10", ymin=-1, ymax=1, r0=0.5, r1=4/MAX_CN,lty=2)
  kpAbline(kp, h=c(1), col="gray10", ymin=-1, ymax=1, r0=0.5, r1=1/MAX_CN,lty=2)
  kpAbline(kp, h=c(1), col="gray10", ymin=-1, ymax=1, r0=0.5, r1=3/MAX_CN,lty=2)
  for ( i in 1:nrow(Somatic_SV)) {
    
    Somatic_SV_select=Somatic_SV[i,]
    Somatic_SV_torange=Somatic_SV_select[,c("chr","start","end")]%>%toGRanges()
    
    kpPlotRegions(kp,Somatic_SV_torange,col=Somatic_SV_select$color,r0 =Somatic_SV_select$r0, r1 = Somatic_SV_select$r1)
    kpText(kp, data.panel = 1,chr=Somatic_SV_select$chr, x=Somatic_SV_select$start, y=Somatic_SV_select$r1, labels=Somatic_SV_select$symbol,col=Somatic_SV_select$color, pos=1,cex=1)
    
  }
  
  for ( i in 1:nrow(Germline_SV)) {
    
    Germline_SV_select=Germline_SV[i,]
    Germline_SV_torange=Germline_SV_select[,c("chr","start","end")]%>%toGRanges()
    
    kpPlotRegions(kp,Germline_SV_torange,col=Germline_SV_select$color,r0 =Germline_SV_select$r0, r1 = Germline_SV_select$r1)
    kpText(kp, data.panel = 1,chr=Germline_SV_select$chr, x=Germline_SV_select$start, y=Germline_SV_select$r1, labels=Germline_SV_select$symbol,col=Germline_SV_select$color, pos=2)
    
  }
}

Somatic_seg_Norm_RDR_plot_2 <- function(Somatic_SV, Germline_SV = NULL) {
  load("/Volumes/lyu.yang/MAC_1/R/Public_database/Reference/Chrosmome_cytobank.rdata")
  
  Somatic_SV= Somatic_SV%>%dplyr::rename(chr=Chromosome, start=Start, end=End)
  MAX_RDR <- max(Somatic_SV$Norm_RDR, na.rm = TRUE)%>%ceiling()
  if (MAX_RDR<2) {
    MAX_RDR=2
  }
  # Set up the karyotype plot
  kp <- plotKaryotype(genome = "hg38", chromosomes = paste0("chr", 1:22), plot.type = 4)
  kpDataBackground(kp, data.panel = 1, col = "grey90")
  kpPlotRegions(kp, arm_q, col = "white", r0 = 0, r1 = 1)
  
  # Add axis to the plot
  
  # kpAxis(kp, side = 1,labels = NA,numticks = 1)
  # kpAxis(kp, r1 = 1, ymin = 0, ymax = 1, numticks = MAX_RDR  , cex = 2, data.panel = 1, tick.len = 1,labels = NA)
  
  y_ticks <- seq(0, MAX_RDR, by = 0.5)
  
  y_labels <- as.character(y_ticks)  # Display labels for 0, 0.5, 1, and MAX_RDR
  
  # Add y-axis to the plot
  kpAxis(kp, side = 1, labels = NA, numticks = 1)  # Side 1 may not need labels if it's redundant or for stylistic choices
  kpAxis(kp, r1 = 1, ymin = 0, ymax = 1, numticks = MAX_RDR*2+1,cex = 2, data.panel = 1, tick.len = 1,labels = y_ticks)
  
  
  # Add horizontal lines
  for (i in c(0.5,1.5,2)) {
    kpAbline(kp, h = 1, col = "gray10", ymin = -1, ymax = 1, r0 = 0.5, r1 = i / MAX_RDR, lty = 2)
    
  }
  kpAbline(kp, h = 1, col = "blue", ymin = -1, ymax = 1, r0 = 0.5, r1 = 1 / MAX_RDR, lty = 2,lwd=2)
  
  
  # Plot somatic structural variants with different colors for each clonal
  for (i in seq_len(nrow(Somatic_SV))) {
    Somatic_SV_select <- Somatic_SV[i, ]
    
    Somatic_SV_torange <- Somatic_SV_select[, c("chr", "start", "end")] %>% toGRanges()
    kpPlotRegions(kp, Somatic_SV_torange, col= "black", border =NA,r0 = (Somatic_SV_select$Norm_RDR / MAX_RDR) - 0.01, r1 = Somatic_SV_select$Norm_RDR/ MAX_RDR)
    
    
    
    #    kpPlotRegions(kp, Somatic_SV_torange, col = Somatic_SV_select$line_color, border = Somatic_SV_select$line_color, r0 =  (Somatic_SV_select$CN/MAX_CN)-0.05, r1 = Somatic_SV_select$CN/MAX_CN)
    #  kpText(kp, data.panel = 1, chr = Somatic_SV_select$chr, x = Somatic_SV_select$start, y = Somatic_SV_select$r1, labels = Somatic_SV_select$symbol, col = Somatic_SV_select$line_color, pos = 1, cex = 1)
  }
  
  # Plot germline structural variants
  if (!is.null(Germline_SV)) {
    for (i in seq_len(nrow(Germline_SV))) {
      Germline_SV_select <- Germline_SV[i, ]
      Germline_SV_torange <- Germline_SV_select[, c("chr", "start", "end")] %>% toGRanges()
      
      kpPlotRegions(kp, Germline_SV_torange, col = Germline_SV_select$color, r0 = Germline_SV_select$r0, r1 = Germline_SV_select$r1)
      kpText(kp, data.panel = 1, chr = Germline_SV_select$chr, x = Germline_SV_select$start, y = Germline_SV_select$r1, labels = Germline_SV_select$symbol, col = Germline_SV_select$color, pos = 2)
    }
  }
  
  
  title(main = "Somatic RDR Plot (after normalization)", cex.main = 1.5, line = 2)
  # Add y-axis label using base R graphics
  mtext("Read depth ratio", side = 2, line = 2.5, cex = 1.2) # Adjust line and cex as needed
}


Somatic_seg_BAF_plot_2 <- function(Somatic_SV, Germline_SV = NULL) {
  
  
  Somatic_SV= Somatic_SV%>%dplyr::rename(chr=Chromosome, start=Start, end=End)
  # Set up the karyotype plot
  kp <- plotKaryotype(genome = "hg38", chromosomes = paste0("chr", 1:22), plot.type = 4)
  kpDataBackground(kp, data.panel = 1, col = "grey90")
  kpPlotRegions(kp, arm_q, col = "white", r0 = 0, r1 = 1)
  
  # Add axis to the plot
  
  kpAxis(kp, side = 1,labels = NA,numticks = 1)
  kpAxis(kp, r1 = 1, ymin = 0, ymax = 0.5, numticks =3, cex = 2, data.panel = 1, tick.len = 1, labels = c(0,0.25,0.5))
  kpAbline(kp, h = 1, col = "gray10", ymin = -1, ymax = 1, r0 = 0.5, r1 = 0.25 / 0.5, lty = 2)
  
  # Plot somatic structural variants with different colors for each clonal
  for (i in seq_len(nrow(Somatic_SV))) {
    Somatic_SV_select <- Somatic_SV[i, ]
    
    Somatic_SV_torange <- Somatic_SV_select[, c("chr", "start", "end")] %>% toGRanges()
    kpPlotRegions(kp, Somatic_SV_torange, col= "black", border =NA,r0 = (Somatic_SV_select$BAF / 0.5) - 0.01, r1 = Somatic_SV_select$BAF/0.5)
    
    
    
    #    kpPlotRegions(kp, Somatic_SV_torange, col = Somatic_SV_select$line_color, border = Somatic_SV_select$line_color, r0 =  (Somatic_SV_select$CN/MAX_CN)-0.05, r1 = Somatic_SV_select$CN/MAX_CN)
    #  kpText(kp, data.panel = 1, chr = Somatic_SV_select$chr, x = Somatic_SV_select$start, y = Somatic_SV_select$r1, labels = Somatic_SV_select$symbol, col = Somatic_SV_select$line_color, pos = 1, cex = 1)
  }
  
  # Plot germline structural variants
  if (!is.null(Germline_SV)) {
    for (i in seq_len(nrow(Germline_SV))) {
      Germline_SV_select <- Germline_SV[i, ]
      Germline_SV_torange <- Germline_SV_select[, c("chr", "start", "end")] %>% toGRanges()
      
      kpPlotRegions(kp, Germline_SV_torange, col = Germline_SV_select$color, r0 = Germline_SV_select$r0, r1 = Germline_SV_select$r1)
      kpText(kp, data.panel = 1, chr = Germline_SV_select$chr, x = Germline_SV_select$start, y = Germline_SV_select$r1, labels = Germline_SV_select$symbol, col = Germline_SV_select$color, pos = 2)
    }
  }
  
  
  title(main = "Somatic BAF Plot", cex.main = 1.5, line = 2)
  mtext("BAF", side = 2, line = 2.5, cex = 1.2) # Adjust line and cex as needed
}


generate_normal_plots <- function(data, UMAP_data ,n_cluster=8, seed = 123456,Cluster_select,Bias="") {
  # Select relevant columns for UMAP
  # UMAP
  set.seed(seed)
  UMAP_out <- umap::umap(UMAP_data, n_neighbors = n_cluster, n_components = 2, metric = "euclidean")$layout %>%
    as.data.frame() %>%
    setNames(c("UMAP1", "UMAP2"))
  
  set.seed(seed)
  UMAP_out_kmeans <- kmeans(UMAP_out, centers = n_cluster, nstart = 20)
  
  # Combine original data with UMAP output and k-means clusters
  UMAP_out2 <- cbind(data, UMAP_out) %>%
    mutate(Cluster = as.character(UMAP_out_kmeans$cluster))
  
  # Create summary table
  df_summary <- UMAP_out2 %>%
    mutate(Cluster = factor(Cluster, levels = unique(Cluster))) %>%
    group_by(Cluster) %>%
    summarize(RDR = mean(RDR), BAF = mean(BAF)) %>%
    arrange(RDR)%>%mutate(CNV=ifelse(Cluster%in%Cluster_select ,"Normal","Other"))
  
  if (Bias=="") {
    Bias=df_summary%>%dplyr::filter(CNV=="Normal")%>%.$RDR
    
  } 
  
  
  
  # Annotate data with clusters
  UMAP_out2$Cluster <- factor(UMAP_out2$Cluster, levels = unique(UMAP_out2$Cluster))
  
  # Create label data frame
  label.df <- data.frame(Cluster = unique(UMAP_out2$Cluster), label = paste0("cluster: ", unique(UMAP_out2$Cluster)))
  
  # Normal_CNV
  UMAP_out3=UMAP_out2%>%left_join(.,df_summary[,c("Cluster","CNV")])
  UMAP_plot_Normal <- ggplot(UMAP_out3, aes(x = BAF, y = RDR, color = CNV)) +
    geom_point(size = 3)  +ylim(c(0,4))+geom_hline(yintercept=1, linetype="dashed", color = "red")+
    ggtitle(paste0("Normal cluster:", n_cluster,"in the raw data"))
  
  UMAP_out4=UMAP_out3%>%mutate(Norm_RDR=RDR/Bias)  
  UMAP_plot_Normal_bias <- ggplot(UMAP_out4, aes(x = BAF, y = Norm_RDR, color = CNV)) +
    geom_point(size = 3) +ylim(c(0,4))+geom_hline(yintercept=1, linetype="dashed", color = "red")+
    ggtitle(paste0("Normal cluster:", n_cluster,"in the data with bias normalization"))
  
  
  
  return(list(UMAP_plot_list = cowplot::plot_grid(UMAP_plot_Normal,UMAP_plot_Normal_bias), Norm_RDR = UMAP_out4,Bias_set=Bias))
}

generate_umap_plots <- function(data, UMAP_data ,n_cluster = 10, seed = 123456) {
  # Select relevant columns for UMAP
  # UMAP
  set.seed(seed)
  UMAP_out <- umap::umap(UMAP_data, n_neighbors = n_cluster, n_components = 2, metric = "euclidean")$layout %>%
    as.data.frame() %>%
    setNames(c("UMAP1", "UMAP2"))
  
  set.seed(seed)
  UMAP_out_kmeans <- kmeans(UMAP_out, centers = n_cluster, nstart = 20)
  
  # Combine original data with UMAP output and k-means clusters
  UMAP_out2 <- cbind(data, UMAP_out) %>%
    mutate(Cluster = as.character(UMAP_out_kmeans$cluster))
  
  # Create summary table
  df_summary <- UMAP_out2 %>%
    mutate(Cluster = factor(Cluster, levels = unique(Cluster))) %>%
    group_by(Cluster) %>%
    summarize(RDR = mean(RDR), BAF = mean(BAF)) %>%
    arrange(RDR)%>%mutate(CNV=case_when(BAF>0.45&RDR>0.85&RDR<1.05 ~"Normal",TRUE ~ "Other"))
  
  # Annotate data with clusters
  UMAP_out2$Cluster <- factor(UMAP_out2$Cluster, levels = unique(UMAP_out2$Cluster))
  
  # Print summary table
  summary_table <- knitr::kable(df_summary, caption = "UMAP cluster") %>%
    kableExtra::kable_styling(latex_options = c("scale_down"))
  print(summary_table)
  
  # Create label data frame
  label.df <- data.frame(Cluster = unique(UMAP_out2$Cluster), label = paste0("cluster: ", unique(UMAP_out2$Cluster)))
  
  label.df_2 <- UMAP_out2 %>%
    mutate(Cluster = factor(Cluster, levels = unique(Cluster))) %>%
    group_by(Cluster) %>%
    summarize(UMAP1 = min(UMAP1), UMAP2 = max(UMAP2)) %>%
    left_join(label.df, by = "Cluster")
  
  # UMAP plot
  UMAP_plot <- ggplot(UMAP_out2, aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(colour = as.factor(Cluster))) +
    theme(legend.position = "none") +
    ggrepel::geom_label_repel(data = label.df_2, aes(label = label), vjust = 1) +
    ggtitle("UMAP plot")
  
  
  
  #UMAP-based k-means plot
  UMAP_plot_RDR <- ggplot(UMAP_out2, aes(x = BAF, y = RDR, color = Cluster)) +
    geom_point(size = 3) +
    ggtitle("UMAP-based Kmeans")
  
  Plot_list=cowplot::plot_grid(UMAP_plot,UMAP_plot_RDR)
  
  return(list(UMAP_plot_list = Plot_list, summary_table = df_summary))
}

plot_density_with_vline <- function(data, x_col = "RDR") {
  ggplot(data, aes_string(x = x_col)) +
    geom_density(fill = "blue", alpha = 0.1) +
    geom_vline(xintercept = 1, color = "red", linetype = "dashed") +
    labs(x = x_col, y = "Density", title = paste("Density Plot of", x_col)) +
    theme_minimal()
}

find_peaks <- function(density_estimate) {
  y <- density_estimate$y
  peaks <- c(FALSE, diff(sign(diff(y))) == -2, FALSE)
  return(peaks)
}

density_with_hdi <- function(data, bandwidth = 0.03, credMass = 0.9, initial_group_number = 1) {
  library(HDInterval)
  
  # Compute the density estimate
  density_estimate <- density(data, bw = bandwidth)
  
  # Find peaks
  peaks <- find_peaks(density_estimate)
  peak_values <- density_estimate$x[peaks]
  
  # Compute HDI
  hdi_result <- hdi(density_estimate, credMass = credMass, allowSplit = TRUE)
  hdi_df <- as.data.frame(hdi_result)
  #%>%dplyr::filter(begin<0.85&begin>0.1)
  
  # Filter peaks to include only those within HDI intervals
  in_hdi <- sapply(peak_values, function(pv) any(hdi_df$begin <= pv & pv <= hdi_df$end))
  filtered_peak_values <- peak_values[in_hdi]
  
  # Filter out peak values greater than 0.9
  filtered_peak_values <- filtered_peak_values
  #[filtered_peak_values <= 0.9&filtered_peak_values>0.1 ]
  
  # Find the corresponding HDI interval for each peak
  min_values <- sapply(filtered_peak_values, function(pv) hdi_df$begin[which(hdi_df$begin <= pv & pv <= hdi_df$end)])
  max_values <- sapply(filtered_peak_values, function(pv) hdi_df$end[which(hdi_df$begin <= pv & pv <= hdi_df$end)])
  
  # Create a dataframe for filtered peak values with group numbers, min_value, and max_value
  peak_df <- data.frame(
    Value = round(filtered_peak_values,2),
    Group = seq(initial_group_number, length.out = length(filtered_peak_values)),
    min_value = min_values,
    max_value = max_values
  )
  
  # Plot the density with peaks and HDI highlighted
  plot(density_estimate, main = "Density with Peaks and HDI Highlighted")
  
  # Find y values for the filtered peaks
  peak_y_values_hi <- approx(density_estimate$x, density_estimate$y, xout = peak_df$Value)$y
  points(peak_df$Value, peak_y_values_hi, col = "red", pch = 19)
  
  peak_y_values <- approx(density_estimate$x, density_estimate$y, xout = hdi_df$begin)$y
  points(hdi_df$begin,  peak_y_values, col = "blue", pch = 19)
  
  peak_y_values <- approx(density_estimate$x, density_estimate$y, xout = hdi_df$end)$y
  points(hdi_df$end, peak_y_values, col = "green", pch = 19)
  
  peak_df2=peak_df%>%mutate(Y=peak_y_values_hi)
  # Return results as a list
  result_list <- list(
    density_estimate = density_estimate,
    peak_values_df = peak_df2,
    hdi_result = hdi_result
  )
  
  return(result_list)
}


calculate_subclone_purities <- function(data) {
  data <- data %>%
    mutate(Subclone_p = case_when(
      CNV == "LOH" ~ calculate_tumor_purity_by_BAF_LOH(BAF),
      grepl("Tetrasomy",CNV)   ~ calculate_tumor_purity_tetrasomy_by_RDR(Norm_RDR),
      grepl("Trisomy",CNV)  ~ calculate_tumor_purity_trisomy_by_RDR(Norm_RDR),
      CNV == "Het_loss" ~ calculate_tumor_purity_hetloss_by_RDR(Norm_RDR),
      CNV == "Homo_loss" ~ 1 - Norm_RDR,
      CNV == "Normal" ~ 1,
      TRUE ~ NA_real_  # In case CNV doesn't match any condition
    ))
  
  return(data)
}

Trisomy_model_predict <- function(model, RDR, BAF, threshold=0.5) {
  
  # Validate model input
  if (missing(model) || !inherits(model, "lm")) {
    stop("Please provide a valid linear model.")
  }
  
  # Create a new data frame for the specific point
  specific_point <- data.frame(RDR = RDR)
  
  # Predict the BAF for this specific RDR using the model
  predicted_BAF <- predict(model, newdata = specific_point)
  
  # Calculate the residual (difference between actual and predicted BAF)
  residual <- BAF - predicted_BAF
  
  # Determine if the point is close to the line based on the threshold
  is_close <- abs(residual) <= threshold
  
  # Output the result along with some informative data
  list(
    Predicted_BAF = predicted_BAF,
    Actual_BAF = BAF,
    Residual = residual,
    Is_Close_To_Line = is_close,
    Message = if (is_close) "The point is close to the regression line." else "The point is not close to the regression line."
  )
}


Gene_filter_from_segmentation_2=function( Gene_chorom, data) {
  nchr= Gene_chorom$chr
  colnames(data)=tolower(colnames(data))
  seg_data=data%>%dplyr::filter(chr==nchr)
  
  edb=EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
  
  library(ensembldb)
  chr=gsub("chr","" ,nchr)
  data2=data[,c("chr","start","end")]
  detail_region <-  regioneR::overlapRegions(Gene_chorom,data2)
  ## Query genes:
  Start_from=min(detail_region$startB)
  End_with=max(detail_region$endB)
  output=c(Start_from,End_with)
  Gene_filter=data %>%dplyr::filter(chr==nchr)%>% dplyr::filter(start>=output[1]) %>%dplyr:: filter(end<=output[2]) 
}



WEX_gene_somatic_GATK_TN_ratio_2=function (Gene_chorom,Patient_ID,PT_RD_file,Normal_GATK_path,Hatchat_path,Tumor_type){
  #PT_RD_file=Rd_files[grepl(Patient_ID,Rd_files )]
  RD_data=readr::read_tsv(PT_RD_file)%>%as.data.frame()%>%setNames(c("CHR","Other"))%>%dplyr::filter(CHR!="@SQ")%>%dplyr::filter(CHR!="@RG")%>%dplyr::filter(CHR!="CONTIG")%>%tidyr::separate(Other,c("START","END","RDR"),sep="\t")%>%mutate(RDR=2^(as.numeric(RDR)),START=as.integer(START),END=as.integer(END))
  
  Gene_filter=Gene_filter_from_segmentation_2(Gene_chorom,RD_data)  
  
  return(Gene_filter)
  
}


Chrom_GATK_Gene_check=function(Gene_chorom,PT_RD_file,Main){
  RD_data=readr::read_tsv(PT_RD_file)%>%as.data.frame()%>%setNames(c("CHR","Other"))%>%dplyr::filter(CHR!="@SQ")%>%dplyr::filter(CHR!="@RG")%>%dplyr::filter(CHR!="CONTIG")%>%tidyr::separate(Other,c("START","END","RDR"),sep="\t")%>%mutate(RDR=2^(as.numeric(RDR)),START=as.integer(START),END=as.integer(END))
  Tumor_bbc_reform=RD_data[,c("CHR","START","END","RDR")] 
  colnames(Tumor_bbc_reform)=tolower(colnames(Tumor_bbc_reform))
  nchr=Gene_chorom$chr
  
  s1 <- CopyNumberPlots:: loadSNPData(Tumor_bbc_reform)
  kp <- karyoploteR::plotKaryotype(genome="hg38",chromosomes=nchr,main=paste0(Main ,"Read depth ratio"))
  df=data.frame(chr=nchr, start=Gene_chorom$start,end=Gene_chorom$end)
  karyoploteR::kpAddLabels(kp, labels="Sample3", r0=2, r1=2.30, data.panel = 1, cex=2)
  karyoploteR::kpPlotRegions(kp, data=df, col="#F30943", border="#F30943", r0=0.3, r1=0.55)
  karyoploteR::kpAbline(kp, h=0.5, col="green")
  karyoploteR::kpAbline(kp, h=0.25, col="black")
  karyoploteR::kpText(kp, chr=nchr, x=Gene_chorom$start, y=0.75, labels=Gene_chorom$symbol, data.panel = 1)
  CopyNumberPlots::plotLRR(kp, s1,lrr.column ="rdr",r0=0, r1=1,labels ="Copy number",  ymin=0, ymax=2,add.axis=T, axis.cex=1,label.cex=1, points.cex = 0.5)
  
}


Plot_GATK_RD_in_chr_zoom=function (Gene_chorom,BBC,Main,ref_line ) {
  colnames(BBC)=toupper(colnames(BBC))
  Tumor_bbc_reform=BBC[,c("CHR","START","END","RDR")]
  colnames(Tumor_bbc_reform)=tolower(colnames(Tumor_bbc_reform))
  nchr=Gene_chorom$chr
  Gene_chorom=Gene_chorom[,c("chr","start","end")]
  detail.region <-  regioneR::toGRanges(Gene_chorom,BBC,Main,genome="hg38")
  s1 <- CopyNumberPlots:: loadSNPData(Tumor_bbc_reform)
  
  kp <- karyoploteR::plotKaryotype(genome="hg38",chromosomes=nchr,main=paste(Main,Gene_chorom$symbol, "Read depth ratio"),zoom=detail.region)
  karyoploteR::kpAddLabels(kp, labels="Sample3", r0=2, r1=2.30, data.panel = 1, cex=2)
  #karyoploteR::kpPlotRegions(kp, data=df, col="#F30943", border="#F30943", r0=0.3, r1=0.55)
  # karyoploteR::kpAbline(kp, h=0.5, col="green")
  karyoploteR::kpAbline(kp, h=ref_line, col="black")
  karyoploteR::kpText(kp, chr=nchr, x=Gene_chorom$start+10000, y=0.75, labels=Gene_chorom$symbol, data.panel = 1,cex=0.8)
  karyoploteR::kpText(kp, chr=nchr, x=Gene_chorom$start+15000, y=0.25, labels=as.character(Gene_chorom$start), data.panel = 1,cex=0.8)
  karyoploteR::kpText(kp, chr=nchr, x=Gene_chorom$end-15000, y=0.25, labels=as.character(Gene_chorom$end), data.panel = 1,cex=0.8)
  CopyNumberPlots::plotLRR(kp,s1, lrr.column ="rdr",r0=0, r1=1,labels ="Read depth ratio", ymin=0, ymax=2,add.axis=T, axis.cex=1,label.cex=1, points.cex =1.5)
  
  
}

add_gene_names_to_BAF <- function(BAF_data) {
  
  # Load EnsDb annotation database
  edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
  
  # Ensure column names are in lower case
  colnames(BAF_data) <- tolower(colnames(BAF_data))
  
  for (i in seq_len(nrow(BAF_data))) {
    # Extract chromosome, start, and end positions
    nchr <- gsub("chr", "", BAF_data[i, "chr"])  # Remove "chr" prefix if present
    start <- BAF_data[i, "position"]%>%as.integer()
    end <- BAF_data[i, "position"]%>%as.integer()
    
    
    # Create a GRanges object for this region
    region_gr <- GRanges(seqnames = nchr, ranges = IRanges(start = start, end = end))
    
    # Create a GRangesFilter for the region
    grf <- GRangesFilter(region_gr, type = "any")
    
    # Query genes that overlap with the region using the EnsDb database
    genes_in_region <- genes(edb, filter = grf)
    gene_info=genes_in_region%>%as.data.frame()%>%.[,c(6:8)]%>%dplyr::filter(grepl("protein_coding" ,gene_biotype))%>%dplyr::filter(!grepl("RP11" ,gene_name))
    # Filter for protein-coding genes and extract gene names
    # gene_names <- genes_in_region %>%
    #   as.data.frame() %>%
    #   filter(gene_biotype == "protein_coding") %>%
    #   pull(gene_name)
    
    # Add gene names to the BAF_data (if any)
    if (nrow(gene_info) ==1) {
      
      BAF_data[i, "gene_id"] =  gene_info[,1]
      BAF_data[i, "gene_name"] =  gene_info[,2]
      BAF_data[i,"gene_biotype"]=gene_info[,3]
    } else if(nrow(gene_info) >1){
      BAF_data[i, "gene_id"] =  gene_info[,1]%>%toString()
      BAF_data[i, "gene_name"] =  gene_info[,2]%>%toString()
      BAF_data[i,"gene_biotype"]=gene_info[,3]%>%toString()
    }
    else{
      BAF_data[i, "gene_id"] =  ""
      BAF_data[i, "gene_name"] =  ""
      BAF_data[i,"gene_biotype"]=""
    }
  }
  
  return(BAF_data)
}

label_bed_with_genes <- function(bed_data) {
  colnames(bed_data)=tolower(colnames(bed_data))
  
  # Ensure the bed_data has the necessary columns
  if (!all(c("chr", "start", "end") %in% colnames(bed_data))) {
    stop("Input data must have columns named 'chr', 'start', and 'end'")
  }
  
  # Load the human genome database
  edb <- EnsDb.Hsapiens.v86
  
  # Create an empty column for gene names in the input bed_data
  # Loop through each row in the bed data
  for (i in seq_len(nrow(bed_data))) {
    # Extract chromosome, start, and end positions
    nchr <- gsub("chr", "", bed_data[i, "chr"])  # Remove "chr" prefix if present
    start <- bed_data[i, "start"]
    end <- bed_data[i, "end"]
    
    # Create a GRanges object for this region
    region_gr <- GRanges(seqnames = nchr, ranges = IRanges(start = start, end = end))
    
    # Create a GRangesFilter for the region
    grf <- GRangesFilter(region_gr, type = "any")
    
    # Query genes that overlap with the region using the EnsDb database
    genes_in_region <- genes(edb, filter = grf)
    gene_info=genes_in_region%>%as.data.frame()%>%.[,c(6:8)]%>%dplyr::filter(grepl("protein_coding" ,gene_biotype))%>%dplyr::filter(!grepl("RP11" ,gene_name))
    # Filter for protein-coding genes and extract gene names
    # gene_names <- genes_in_region %>%
    #   as.data.frame() %>%
    #   filter(gene_biotype == "protein_coding") %>%
    #   pull(gene_name)
    
    # Add gene names to the bed_data (if any)
    if (nrow(gene_info) ==1) {
      
      bed_data[i, "gene_id"] =  gene_info[,1]
      bed_data[i, "gene_name"] =  gene_info[,2]
      bed_data[i,"gene_biotype"]=gene_info[,3]
    } else if(nrow(gene_info) >1){
      bed_data[i, "gene_id"] =  gene_info[,1]%>%toString()
      bed_data[i, "gene_name"] =  gene_info[,2]%>%toString()
      bed_data[i,"gene_biotype"]=gene_info[,3]%>%toString()
    }
    else{
      bed_data[i, "gene_id"] =  ""
      bed_data[i, "gene_name"] =  ""
      bed_data[i,"gene_biotype"]=""
    }
  }
  
  return(bed_data)
}

detect_segments_in_rdr <- function(RDR_data, lmax = 10) {
  # Initialize necessary variables
  chri <- as.character(RDR_data$chr[1])  # Assuming seqnames contains chromosome info
  num_positions <- nrow(RDR_data)
  rdr_values <- RDR_data$rdr  # Extract RDR values
  gene_names <- RDR_data$gene_name  # Extract gene names if available
  
  # Initialize output matrix
  final_call <- matrix(ncol = 7)
  
  # Pad with zeros for segment lengths
  y <- c(rdr_values, rep(0, lmax))
  
  i <- rep(1:num_positions, rep(lmax, num_positions))
  j <- rep(1:lmax, num_positions) + i
  yact <- rep(0, length(i))
  
  # Compute cumulative sums for yact
  for (k in 1:num_positions) {
    yact[(lmax * k - (lmax - 1)):(lmax * k)] <- cumsum(y[k:(k + lmax)])[-1]
  }
  
  # Filter valid indices
  i <- i[j <= num_positions]
  j <- j[j <= num_positions]
  yact <- yact[j <= num_positions]
  
  # Prevent yact values from dropping too low
  yact[yact < 20] <- 20
  
  # Calculate the log-likelihood ratio (lratio)
  lratio <- (yact - mean(yact))^2 / mean(yact)  # Squared difference from the mean
  
  # Filter positive log-likelihood ratios and detect segments
  if (sum(lratio > 0) > 0) {
    final_mat <- cbind(i, j, yact, lratio)
    final_mat <- final_mat[lratio > 0, ]
    final_mat <- final_mat[order(-final_mat[, 4]), ]
    
    # Resolve overlapping segments by keeping the highest scoring
    s <- 1
    while (s <= nrow(final_mat)) {
      row_start <- final_mat[s, 1]
      row_end <- final_mat[s, 2]
      row_sel <- (final_mat[, 1] <= row_end & final_mat[, 2] >= row_start)
      row_sel[s] <- FALSE
      final_mat <- final_mat[!row_sel, ]
      if (is.vector(final_mat)) final_mat <- t(as.matrix(final_mat))
      s <- s + 1
    }
    
    final_mat <- round(final_mat, digits = 3)
    
    # Model selection with mBIC
    loglikeij <- cumsum(final_mat[, 4])
    mBIC <- rep(NA, length(loglikeij))
    for (s in 1:nrow(final_mat)) {
      tau <- sort(unique(c(as.vector(final_mat[1:s, 1:2]), 1, num_positions)))
      P <- length(tau) - 2
      mbic <- loglikeij[s]
      mbic <- mbic - 0.5 * sum(log(tau[2:length(tau)] - tau[1:(length(tau) - 1)]))
      mbic <- mbic + (0.5 - P) * log(num_positions)
      mBIC[s] <- mbic
    }
    
    mBIC <- round(mBIC, digits = 3)
    if (mBIC[1] > 0) {
      final_call <- cbind(
        rep(chri, nrow(final_mat)),
        rep(gene_names[1], nrow(final_mat)),  # Assuming all positions are from the same gene
        final_mat[1:which.max(mBIC), ],
        mBIC[1:which.max(mBIC)]
      )
    }
  }
  
  # Return the segments
  final_call <- as.data.frame(final_call)
  colnames(final_call) <- c("Chromosome", "Gene", "Start_Index", "End_Index", "Observed_Count", "Log_Likelihood", "mBIC")
  return(final_call)
}


summarize_rdr_by_gene <- function(RDR_data) {
  chr_anno=RDR_data[,c("chr","gene_name")]%>%unique()
  # Summarize by grouping on 'gene_name' and calculating mean RDR, min start, and max end
  summarized_data <- RDR_data %>%
    dplyr::group_by(gene_name) %>%
    dplyr::summarise(
      rdr = mean(rdr, na.rm = TRUE),
      start = min(start),
      end = max(end)
    )%>%dplyr::filter(gene_name!="")%>%left_join(.,chr_anno)
  return(summarized_data)
}




summarize_baf_by_gene <- function(BAF_data) {
  # Summarize by grouping on 'gene_name' and calculating mean RDR, min start, and max end
  summarized_data <- BAF_data %>%
    dplyr::group_by(gene_name) %>%
    dplyr::summarise(
      baf = mean(baf, na.rm = TRUE),
      N_baf=n()
      
    )%>%dplyr::filter(gene_name!="")
  return(summarized_data)
}

detect_segments_in_rdr_by_gene <- function(rdr_data, threshold = 0.1, min_segment_size = 1) {
  # Ensure the input data is arranged by chromosome and start position
  rdr_data <- rdr_data %>%
    dplyr::arrange(chr, start)
  
  # Initialize variables for the first segment
  segment_start <- rdr_data$start[1]
  segment_end <- rdr_data$end[1]
  segment_genes <- rdr_data$gene_name[1]
  segment_rdr_values <- rdr_data$rdr[1]
  
  segments <- list()
  segment_index <- 1
  
  # Iterate through the rdr_data rows
  for (i in 2:nrow(rdr_data)) {
    current_gene <- rdr_data$gene_name[i]
    current_rdr <- rdr_data$rdr[i]
    
    # Check if the current gene RDR value is within the threshold compared to the last RDR value
    if (abs(current_rdr - segment_rdr_values) <= threshold) {
      # Continue the segment
      segment_end <- rdr_data$end[i]
      segment_genes <- c(segment_genes, current_gene)
      segment_rdr_values <- c(segment_rdr_values, current_rdr)
    } else {
      # If the segment is large enough, save it
      if (length(segment_genes) >= min_segment_size) {
        segments[[segment_index]] <- data.frame(
          chr = rdr_data$chr[i - 1],
          start = segment_start,
          end = segment_end,
          genes = paste(segment_genes, collapse = ","),
          mean_rdr = mean(segment_rdr_values),
          segment_size = length(segment_genes)
        )
        segment_index <- segment_index + 1
      }
      
      # Start a new segment with the current gene
      segment_start <- rdr_data$start[i]
      segment_end <- rdr_data$end[i]
      segment_genes <- current_gene
      segment_rdr_values <- current_rdr
    }
  }
  
  # Add the last segment if it meets the size threshold
  if (length(segment_genes) >= min_segment_size) {
    segments[[segment_index]] <- data.frame(
      chr = rdr_data$chr[nrow(rdr_data)],
      start = segment_start,
      end = segment_end,
      genes = paste(segment_genes, collapse = ","),
      mean_rdr = mean(segment_rdr_values),
      segment_size = length(segment_genes)
    )
  }
  
  # Combine all segments into a single data frame
  final_segments <- do.call(rbind, segments)
  
  # Define a function to compute log-likelihood
  log_likelihood <- function(segment_rdr_values) {
    segment_mean <- mean(segment_rdr_values)
    ll <- sum(dnorm(segment_rdr_values, mean = segment_mean, sd = sd(segment_rdr_values), log = TRUE))
    return(ll)
  }
  
  # Calculate mBIC for each segment using log-likelihood
  n <- nrow(rdr_data)
  final_segments$mBIC <- sapply(1:nrow(final_segments), function(s) {
    segment_size <- final_segments$segment_size[s]
    #segment_rdr_values <- as.numeric(unlist(strsplit(final_segments$mean_rdr[s], ",")))
    segment_rdr_values <- as.numeric(final_segments$mean_rdr)
    log_lik <- log_likelihood(segment_rdr_values)
    
    P <- segment_size  # Number of parameters
    mBIC <- log_lik - 0.5 * P * log(n)
    return(mBIC)
  })
  
  # Select the segments with the maximum mBIC
  best_segments <- final_segments[final_segments$mBIC == max(final_segments$mBIC), ]
  
  return(list(final_segments,best_segments))
}



detect_segments_in_rdr <- function(rdr_data, threshold = 0.1, min_segment_size = 3) {
  # Sort data by chromosome and start position
  rdr_data <- rdr_data %>%
    dplyr::arrange(chr, start)
  
  # Initialize variables
  segment_start <- rdr_data$start[1]
  segment_end <- rdr_data$end[1]
  segment_genes <- rdr_data$gene_name[1]
  segment_rdr_values <- rdr_data$rdr[1]
  
  segments <- list()
  segment_index <- 1
  
  for (i in 2:nrow(rdr_data)) {
    current_gene <- rdr_data$gene_name[i]
    current_rdr <- rdr_data$rdr[i]
    
    # Get the last RDR value in the current segment for comparison
    last_segment_rdr <- tail(segment_rdr_values, 1)
    
    # Check if the RDR value is within the threshold compared to the previous gene
    if (abs(current_rdr - last_segment_rdr) <= threshold) {
      # Continue the segment
      segment_end <- rdr_data$end[i]
      segment_genes <- c(segment_genes, current_gene)
      segment_rdr_values <- c(segment_rdr_values, current_rdr)
    } else {
      # If the current segment is smaller than the minimum size, merge with the previous one
      if (length(segment_genes) < min_segment_size && segment_index > 1) {
        # Merge with the previous segment
        segments[[segment_index - 1]]$end <- segment_end
        segments[[segment_index - 1]]$genes <- paste(segments[[segment_index - 1]]$genes, paste(segment_genes, collapse = ","))
        segments[[segment_index - 1]]$mean_rdr <- mean(c(segments[[segment_index - 1]]$mean_rdr, segment_rdr_values))
        segments[[segment_index - 1]]$segment_size <- segments[[segment_index - 1]]$segment_size + length(segment_genes)
      } else {
        # Otherwise, end the current segment and start a new one
        if (length(segment_genes) >= min_segment_size) {
          segments[[segment_index]] <- data.frame(
            chr = rdr_data$chr[i-1],
            start = segment_start,
            end = segment_end,
            genes = paste(segment_genes, collapse = ","),
            mean_rdr = mean(segment_rdr_values),
            segment_size = length(segment_genes)
          )
          segment_index <- segment_index + 1
        }
      }
      
      # Start a new segment
      segment_start <- rdr_data$start[i]
      segment_end <- rdr_data$end[i]
      segment_genes <- current_gene
      segment_rdr_values <- current_rdr
    }
  }
  
  # Add the last segment, and if too small, merge it with the previous one
  if (length(segment_genes) >= min_segment_size) {
    segments[[segment_index]] <- data.frame(
      chr = rdr_data$chr[nrow(rdr_data)],
      start = segment_start,
      end = segment_end,
      genes = paste(segment_genes, collapse = ","),
      mean_rdr = mean(segment_rdr_values),
      segment_size = length(segment_genes)
    )
  } else if (segment_index > 1) {
    # If the last segment is small, merge with the previous segment
    segments[[segment_index - 1]]$end <- segment_end
    segments[[segment_index - 1]]$genes <- paste(segments[[segment_index - 1]]$genes, paste(segment_genes, collapse = ","))
    segments[[segment_index - 1]]$mean_rdr <- mean(c(segments[[segment_index - 1]]$mean_rdr, segment_rdr_values))
    segments[[segment_index - 1]]$segment_size <- segments[[segment_index - 1]]$segment_size + length(segment_genes)
  }
  
  # Combine all segments into a single data frame
  final_segments <- do.call(rbind, segments)
  final_segments$genes <- sapply(final_segments$genes, function(gene_str) {
    paste(unique(unlist(strsplit(gene_str, ","))), collapse = ",")    
  })
  
  return(final_segments)
  # # Model Selection using mBIC
  # log_likelihood <- function(segment_rdr_values) {
  #   # Calculate the log-likelihood for the segment
  #   segment_mean <- mean(segment_rdr_values)
  #   ll <- sum(dnorm(segment_rdr_values, mean = segment_mean, sd = sd(segment_rdr_values), log = TRUE))
  #   return(ll)
  # }
  
  # Calculate mBIC for each segment
  # n <- nrow(rdr_data)
  # final_segments$mBIC <- sapply(1:nrow(final_segments), function(s) {
  #   segment_size <- final_segments$segment_size[s]
  #   segment_rdr_values <- unlist(strsplit(final_segments$mean_rdr[s], ","))
  #   log_lik <- log_likelihood(segment_rdr_values)
  #   
  #   P <- segment_size  # Number of parameters (segment size)
  #   mBIC <- log_lik - 0.5 * P * log(n)
  #   return(mBIC)
  # })
  # 
  # # Select the best segment based on the maximum mBIC
  # best_segments <- final_segments[final_segments$mBIC == max(final_segments$mBIC), ]
  # 
  # return(best_segments)
}



# Function to plot the segments detected by detect_segments_in_rdr
plot_detected_segments <- function(RD_raw_data, segments,BAF_data,Gene_chorom=NULL,plot_zoom=F) {
  # Load necessary libraries
  library(karyoploteR)
  library(ggplotify)
  library(GenomicRanges)
  library(CopyNumberPlots)
  library(ensembldb)
  library(EnsDb.Hsapiens.v86)
  library(cowplot)
  
  # Call the GATK_Chrom_Gene_RD_2 function to get the necessary data
  colnames(RD_raw_data)=tolower(colnames(RD_raw_data))
  Tumor_bbc_reform=RD_raw_data[,c("chr","start","end","rdr")]
  colnames(Tumor_bbc_reform)=tolower(colnames(Tumor_bbc_reform))
  # Set plot parameters
  pp <- karyoploteR::getDefaultPlotParams(plot.type =5 )
  pp$data1height <- 100
  pp$data2height <- 100
  pp$ideogramheight <- 0
  pp$topmargin <- 10
  
  # Extract data from the test object
  nchr=RD_raw_data$chr%>%unique()
  df <- Tumor_bbc_reform
  
  s1 <- CopyNumberPlots:: loadSNPData(Tumor_bbc_reform)
  
  if (!is.null(Gene_chorom) & (isTRUE(plot_zoom)  ) ) {
    
   detail.region <-  regioneR::toGRanges(Gene_chorom,Tumor_bbc_reform,Main="",genome="hg38")
  
   kp <- karyoploteR::plotKaryotype(genome = "hg38", chromosomes = nchr, plot.params = pp,plot.type = 3,zoom = detail.region)
   point_cex=2
  } else if (!is.null(Gene_chorom)) {
    detail.region <-  regioneR::toGRanges(Gene_chorom,Tumor_bbc_reform,Main="",genome="hg38")
    kp <- karyoploteR::plotKaryotype(genome = "hg38", chromosomes = nchr, plot.params = pp,plot.type = 3) 
    point_cex=0.3
    kpPlotRegions(kp, data = detail.region, col = "#F30943", r0 = 0, r1 = 0.05, border = "blue")
    
  } else {
    kp <- karyoploteR::plotKaryotype(genome = "hg38", chromosomes = nchr, plot.params = pp,plot.type = 3) 
    point_cex=0.3
  }
  karyoploteR::kpAbline(kp, h = 0.5, col = "green")
  karyoploteR::kpAbline(kp, h = 0.25, col = "black")
  # Iterate through segments and plot them using kpPlotRegions
  for (i in 1:nrow(segments)) {
    Data=segments[i,]
    segment_data <- GRanges(seqnames = segments$chr[i],
                            ranges = IRanges(start = segments$start[i], end = segments$end[i]),
                            names = segments$genes[i])
    
    kpPlotRegions(kp, data = segment_data, col = "#F30943", r0 = (Data$rdr)/2-0.02, r1 = (Data$rdr)/2+0.02, border = "#F30943")
    
    # Add text for each segment showing the mean RDR
    if (!is.null(Gene_chorom)&(isTRUE(plot_zoom)  )) {
     kpText(kp, chr = segments$chr[i], x = (segments$start[i] + segments$end[i]) / 2,
            y = 0.5, labels =as.character(Data$genes), cex = 0.8)
 
      
    }
  }
 
  
  # Plot using the chromosome data
  
  CopyNumberPlots::plotLRR(kp, s1, lrr.column = "rdr", r0 = 0, r1 = 1, labels = "Read depth ratio", ymin = 0, ymax = 2, add.axis = TRUE, axis.cex = 1, label.cex = 1, points.cex = point_cex, data.panel = 1) 
  #p1 <- ggplotify::as.ggplot(plot_karyotype_grob_chr(nchr, Sample_name, df))
  #  p1 <- ggplotify::as.ggplot(function() {plot_karyotype_chr (nchr, Sample_name, df, s1)} )
  
  Tumor_baf=BAF_data[,c("chr","start","end","baf")]%>%as.data.frame()
  s2 <- CopyNumberPlots:: loadSNPData(Tumor_baf)
  karyoploteR::kpAbline(kp, h = 0.5, col = "green",data.panel = 2)
  CopyNumberPlots::plotBAF(kp,s2, r0=0, r1=1, labels = "BAF", points.col = "orange", points.cex = point_cex, add.axis = TRUE, points.pch = 16,label.cex = 1, axis.cex = 1, data.panel = 2)
  # return(p1)
  
  
}



