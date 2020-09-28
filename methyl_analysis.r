#script: Differential methylation analysis
#author: Shruti Srivastava


library(methylKit)
library(genomation)

file.list=list('ALS3_793_hg19/ALS3_CG_XYremoved_report.txt',
               'ALS4_794_hg19/ALS4_CG_XYremoved_report.txt',
               'ALS6_798_hg19/ALS6_CG_XYremoved_report.txt',
               'ALS7_797_hg19/ALS7_CG_XYremoved_report.txt',
               'ALS8_799_hg19/ALS8_CG_XYremoved_report.txt',
               'ACRI22_795_hg19/ACRI22_CG_XYremoved_report.txt',
               'ACRI247_796_hg19/ACRI247_CG_XYremoved_report.txt',
               'Prost018_791_hg19/Prost018_CG_XYremoved_report.txt'
)

myobj8 = methRead(location = file.list, sample.id = list("ALS3","ALS4","ALS6", "ALS7", "ALS8", "ACRI22", "ACRI247", "Prost018"), 
                  assembly = "hg19", dbtype = "tabix", pipeline = "bismarkCytosineReport",
                  header = FALSE, skip = 0, sep = "\t", resolution = "region", 
                  treatment = c(1,1,1,1,1,0,0,0), dbdir = "ALS_methyl_call_tile_500", mincov = 10, )

tile500_noXY = tileMethylCounts(myobj8, win.size=500, step.size=500, cov.bases = 0, mc.cores = 1)
tile500_noXY_filtered = filterByCoverage(tile500_noXY, lo.count = 10, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)
tile500_noXY_normalized = normalizeCoverage(tile500_noXY_filtered, method = "median")
tileMeth_2 = unite(tile500_noXY_normalized)

pm8=percMethylation(tileMeth_2)
mds8=matrixStats::rowSds(pm8)

names(mds8)<- 1:length(mds8)

mds8

hist(mds8,col="green",xlab="Std.dev.per CpG")

pdf(file="MethDiff_Analysis_v7_tiled_region/Methylkit.02.ALS_Meth_CorrelationPlot_tile500_v8.pdf",
    width=12, height=12); getCorrelation(tileMeth_2, plot=TRUE); dev.off()

pdf(file="MethDiff_Analysis_v7_tiled_region/Methylkit.03.ALS_PCASamples.ctreeplot_tile500_v8.pdf", 
    width=12, height=12); clusterSamples(tileMeth_2, dist="correlation", method="ward", plot=TRUE); dev.off()

pdf(file="MethDiff_Analysis_v7_tiled_region/Methylkit.03.ALS_PCASamples.screeplot_tile500_v8.pdf",
    width=12,height=12); PCASamples(tileMeth_2, screeplot=TRUE); dev.off()

pdf(file="MethDiff_Analysis_v7_tiled_region/Methylkit.03.ALS_PCASamples.pcaxyplot_tile500_v8.pdf",
    width=12,height=12); PCASamples(tileMeth_2); dev.off()

myDiff_tile500_8=calculateDiffMeth(tileMeth_2)

myDiff25p_hyper_tile500_8 <- getMethylDiff(myDiff_tile500_8,difference=25,qvalue=0.01, type= "hyper",)
myDiff25p_hypo_tile500_8 <- getMethylDiff(myDiff_tile500_8,difference=25,qvalue=0.01, type="hypo",)

pdf("MethDiff_Analysis_v7_tiled_region/Methylkit.04.DiffMethPerChr_tile500_v8.pdf",
    width=12,height=12); diffMethPerChr(myDiff_8,plot=TRUE,qvalue.cutoff=0.01,meth.cutoff=25,
                                        exclude = c("chr6_ssto_hap7","chr6_mcf_hap5","chr6_cox_hap2","chr6_mann_hap4","chr6_apd_hap1","chr6_qbl_hap6","chr6_dbb_hap3","chr17_ctg5_hap1","chr4_ctg9_hap1","chr1_gl000192_random",
                                                    "chrUn_gl000225","chr4_gl000194_random","chr4_gl000193_random","chr9_gl000200_random",
                                                    "chrUn_gl000222","chrUn_gl000212","chr7_gl000195_random","chrUn_gl000223",
                                                    "chrUn_gl000224","chrUn_gl000219","chr17_gl000205_random","chrUn_gl000215",
                                                    "chrUn_gl000216","chrUn_gl000217","chr9_gl000199_random","chrUn_gl000211",
                                                    "chrUn_gl000213","chrUn_gl000220","chrUn_gl000218","chr19_gl000209_random",
                                                    "chrUn_gl000221","chrUn_gl000214","chrUn_gl000228","chrUn_gl000227",
                                                    "chr1_gl000191_random","chr19_gl000208_random","chr9_gl000198_random",
                                                    "chr17_gl000204_random","chrUn_gl000233","chrUn_gl000237","chrUn_gl000230",
                                                    "chrUn_gl000242","chrUn_gl000243","chrUn_gl000241","chrUn_gl000236",
                                                    "chrUn_gl000240","chr17_gl000206_random","chrUn_gl000232","chrUn_gl000234",
                                                    "chr11_gl000202_random","chrUn_gl000238","chrUn_gl000244","chrUn_gl000248",
                                                    "chr8_gl000196_random","chrUn_gl000249","chrUn_gl000246","chr17_gl000203_random",
                                                    "chr8_gl000197_random","chrUn_gl000245","chrUn_gl000247","chr9_gl000201_random",
                                                    "chrUn_gl000235","chrUn_gl000239","chr21_gl000210_random","chrUn_gl000231",
                                                    "chrUn_gl000229","chrM","chrUn_gl000226","chr18_gl000207_random")); dev.off()


bedgraph(myDiff25p_hyper_tile500_8,file.name="MethDiff_Analysis_v7_tiled_region/diff.bedgraph_hyper_tile500_v8",col.name="meth.diff")
bedgraph(myDiff25p_hypo_tile500_8,file.name="MethDiff_Analysis_v7_tiled_region/diff.bedgraph_hypo_tile500_v8",col.name="meth.diff")


#**************************************************************************Annotation******************************************************************************

gene.obj=genomation::readTranscriptFeatures("UCSC_hg19_gene_genepred_ncbirefseq.txt")

diffCpG_regions_ann_hypo_v8 = annotateWithGeneParts(as(myDiff25p_hypo_tile500_8,"GRanges"), gene.obj)
diffCpG_regions_ann_hyper_v8 = annotateWithGeneParts(as(myDiff25p_hyper_tile500_8,"GRanges"), gene.obj)

pdf("MethDiff_Analysis_v7_tiled_region/Methylkit.05.hypoCpG_annotate_tile500_v8.pdf",
    width=12,height=12); plotTargetAnnotation(diffCpG_regions_ann_hypo_v8); dev.off()
pdf("MethDiff_Analysis_v7_tiled_region/Methylkit.05.hyperCpG_annotate_tile500_v8.pdf",
    width=12,height=12); plotTargetAnnotation(diffCpG_regions_ann_hyper_v8); dev.off()

write.table(getAssociationWithTSS(diffCpG_regions_ann_hypo_v8), file="MethDiff_Analysis_v7_tiled_region/TSS_hypo_tile500.v8.txt")
write.table(getAssociationWithTSS(diffCpG_regions_ann_hyper_v8), file="MethDiff_Analysis_v7_tiled_region/TSS_hyper_tile500.v8.txt")

tss.assoc_v8 <- getAssociationWithTSS(diffCpG_regions_ann_hypo_v8)
tss.assoc_hyper_v8 <- getAssociationWithTSS(diffCpG_regions_ann_hyper_v8)

pdf("MethDiff_Analysis_v7_tiled_region/Methylkit.06.TSS_dist_hypoCpG_tile500_v8.pdf",
    width=12,height=12); hist(tss.assoc_v8$dist.to.feature[abs(tss.assoc_v8$dist.to.feature)<=100000],
                              main="distance to nearest TSS",xlab="distance in bp",breaks=50,col="brown4"); dev.off()

pdf("MethDiff_Analysis_v7_tiled_region/Methylkit.06.TSS_dist_hyperCpG_tile500_v8.pdf",
    width=12,height=12); hist(tss.assoc_hyper_v8$dist.to.feature[abs(tss.assoc_hyper_v8$dist.to.feature)<=100000],
                              main="distance to nearest TSS",xlab="distance in bp",breaks=50,col="brown4"); dev.off()

write.table(getMembers(diffCpG_regions_ann_hyper_v8), file="MethDiff_Analysis_v7_tiled_region/CpG_hypermeth_tile500_v8.txt")
write.table(getMembers(diffCpG_regions_ann_hypo_v8), file="MethDiff_Analysis_v7_tiled_region/CpG_hypometh_tile500_v8.txt")

CpG_region_hypo_v8 <- tss.assoc_v8[,1]
CpG_region_hyper_v8 <- tss.assoc_hyper_v8[,1]

write.table(myDiff25p_hypo_tile500_8[CpG_region_hypo_v8,], file="MethDiff_Analysis_v7_tiled_region/Hypo_cpg_region_tile500_v8.tsv", sep='\t', quote=FALSE)
write.table(myDiff25p_hyper_tile500_8[CpG_region_hyper_v8,], file="MethDiff_Analysis_v7_tiled_region/Hyper_cpg_region_tile500_v8.tsv", sep='\t', quote=FALSE)

#annotation of differentially methylated hyper and hypo CpG
cpg.obj <- readFeatureFlank("/Raw_data/shruti/methylation/temp/dragen_ALS/regulation_CpG_anno.bed.txt", feature.flank.name=c("CpGi","shores"))

diffCpGann_hyper_tile500_v8 = annotateWithFeatureFlank(as(myDiff25p_hyper_tile500_8,"GRanges"),
                                    cpg.obj$CpGi,cpg.obj$shores,
                                    feature.name="CpGi",flank.name="shores")

write.table(getMembers(diffCpGann_hyper_tile500_v8), file="MethDiff_Analysis_v7_tiled_region/diffhyper_CpGi_shore_v8.txt")
plotTargetAnnotation(diffCpGann_hyper_tile500_v8, main = "Hypermethylated CpG annotation")

diffCpGann_hypo_tile500_v8 = annotateWithFeatureFlank(as(myDiff25p_hypo_tile500_8,"GRanges"),
                                                       cpg.obj$CpGi,cpg.obj$shores,
                                                       feature.name="CpGi",flank.name="shores")
write.table(getMembers(diffCpGann_hypo_tile500_v8), file="MethDiff_Analysis_v7_tiled_region/diffhypo_CpGi_shore_v8.txt")
plotTargetAnnotation(diffCpGann_hypo_tile500_v8, main = "Hypomethylated CpG annotation")

CpG_annotation_hyper_v8 <- read.table("MethDiff_Analysis_v7_tiled_region/diffhyper_CpGi_shore_v8.txt", header=T,
                                      stringsAsFactors = F, row.names = 1)
CpG_annotation_hypo_v8 <- read.table("MethDiff_Analysis_v7_tiled_region/diffhypo_CpGi_shore_v8.txt", header=T, 
                                     stringsAsFactors = F, row.names = 1)

#adding exon promoter and intron info in the annotation for hypermethylated CpG
hypermeth_CpG <- getData(myDiff25p_hyper_tile500_8)

myDiff25p_hyper_tile500_8_CpGann <- merge(hypermeth_CpG, CpG_annotation_hyper_v8, by=0, all=TRUE)

hyperdiff_exon_prom_intron <- read.table("MethDiff_Analysis_v7_tiled_region/CpG_hypermeth_tile500_v8.txt", header=T,
                                    stringsAsFactors = F, row.names = 1)
hyper_annotate <- cbind(myDiff25p_hyper_tile500_8_CpGann, hyperdiff_exon_prom_intron)
write.table(hyper_annotate, file="hyper_annotate_v8.txt",sep='\t', quote=FALSE)


#adding exon promoter and intron info in the annotation for hypomethylated CpG
hypometh_CpG <- getData(myDiff25p_hypo_tile500_8)

myDiff25p_hypo_tile500_8_CpGann <- merge(hypometh_CpG, CpG_annotation_hypo_v8, by=0, all=TRUE)

hypodiff_exon_prom_intron <- read.table("MethDiff_Analysis_v7_tiled_region/CpG_hypometh_tile500_v8.txt", header=T,
                                         stringsAsFactors = F, row.names = 1)
hypo_annotate <- cbind(myDiff25p_hypo_tile500_8_CpGann, hypodiff_exon_prom_intron)
write.table(hypo_annotate, file="hypo_annotate_v8.txt",sep='\t', quote=FALSE)

library(UpSetR)

hyper_ann <- read.csv("upset_plot_file/hyper_cpg_upsetplot.v2.bed", header=T, sep = "\t")
upset(hypo_ann, sets = c("cpgi", "shores", "prom", "exon", "intron"), sets.bar.color = "#56B4E9",
      +       order.by = "freq", empty.intersections = "on")
hyper_ann <- read.csv("upset_plot_file/hyper_cpg_upsetplot.v2.bed", header=T, sep = "\t")
upset(hyper_ann, sets = c("cpgi", "shores", "prom", "exon", "intron"), sets.bar.color = "#56B4E9",
        +       order.by = "freq", empty.intersections = "on")
        
##############################################################################################################################################
############################################################ Heatmaps ########################################################################
##############################################################################################################################################

#Heatmap for hypomethylated regions
query_hypo_rows_v8<-paste(myDiff25p_hypo_tile500_8[1:295]$chr,":",myDiff25p_hypo_tile500_8[1:295]$start, "-", myDiff25p_hypo_tile500_8[1:295]$end, sep = "")
#dataframe for all DMLoci
all_loci_samples_v8 <- getData(tileMeth_2)
# placing row name as chr:start
rownames(all_loci_samples_v8)<- paste(all_loci_samples_v8$chr, ":", all_loci_samples_v8$start, "-", all_loci_samples_v8$end, sep = "")
# getting only those lines from all_loci_samples that are present in query_rows
query_hypo_res_v8<- all_loci_samples_v8[rownames(all_loci_samples_v8) %in% query_hypo_rows_v8,]
# made a copy
query_hypo_res_natozero_v8 <- query_hypo_res_v8
# replacing na with 0
# query_hypo_res_natozero_v8[is.na(query_hypo_res_natozero_v8)] <- 0
# dividing column having Cs to Coverage and storing in a columns with meth1,2,3...so on as their column names

query_hypo_res_natozero_v8$meth1 <- query_hypo_res_v8[,6]/query_hypo_res_v8[,5]
query_hypo_res_natozero_v8$meth2 <- query_hypo_res_v8[,9]/query_hypo_res_v8[,8]
query_hypo_res_natozero_v8$meth3 <- query_hypo_res_v8[,12]/query_hypo_res_v8[,11]
query_hypo_res_natozero_v8$meth4 <- query_hypo_res_v8[,15]/query_hypo_res_v8[,14]
query_hypo_res_natozero_v8$meth5 <- query_hypo_res_v8[,18]/query_hypo_res_v8[,17]
query_hypo_res_natozero_v8$meth6 <- query_hypo_res_v8[,21]/query_hypo_res_v8[,20]
query_hypo_res_natozero_v8$meth7 <- query_hypo_res_v8[,24]/query_hypo_res_v8[,23]
query_hypo_res_natozero_v8$meth8 <- query_hypo_res_v8[,27]/query_hypo_res_v8[,26]

# again converted na to zero because above computation generated NA values
#query_hypo_res_natozero_v8[is.na(query_hypo_res_natozero_v8)] <- 0
#round the ratio to 2 decimal place; it just rewritten col 29 to 36 in the dataframe while removing other columns
query_hypo_res_natozero_v8<-round(query_hypo_res_natozero_v8[,(29:36)], digits = 2)
#renaming the columns in the dataframe just created above
colnames(query_hypo_res_natozero_v8)[1] = "ALS3"
colnames(query_hypo_res_natozero_v8)[2] = "ALS4"
colnames(query_hypo_res_natozero_v8)[3] = "ALS6"
colnames(query_hypo_res_natozero_v8)[4] = "ALS7"
colnames(query_hypo_res_natozero_v8)[5] = "ALS8"
colnames(query_hypo_res_natozero_v8)[6] = "ACRI22"
colnames(query_hypo_res_natozero_v8)[7] = "ACRI247"
colnames(query_hypo_res_natozero_v8)[8] = "Prost018"
#saves the dataframe to a file

write.table(query_hypo_res_natozero_v8, file="query_hypo_res_natozero_v8.tsv", sep='\t', quote=FALSE)

#Annotation_columns creation
condition<- c("f-ALS","s-ALS","s-ALS","s-ALS","s-ALS","control","control","control")
sex <- c("female","male","female","female","female","male","female","male")

#created a dataframe patients
data.frame(patients, stringsAsFactors = True)
patients <- data.frame(condition,sex)
rownames(patients)<- c("ALS3", "ALS4", "ALS6", "ALS7", "ALS8", "ACRI22", "ACRI247", "Prost018")
write.table(patients, file="ALSpatients.csv", sep=',', quote=FALSE)

library(pheatmap)

CpGloci_hypo_v8 <- read.csv("query_hypo_res_natozero_v8.tsv", header=T, sep='\t',
                             row.names = 1)
#annotation file for column
annotation_col <- read.csv("ALSpatients.csv", header = T, row.names = 1)

#creating annotation_row file
dt<- read.table("CpG_hypometh_tile500_v8.txt", header = T, stringsAsFactors = F, row.names = 1)
head(dt)
v <- ifelse(dt$prom == 1, "Promoter", "Other")
v[dt$exon == 1]<- "Exon"
v[dt$intron == 1]<- "Intron"
dt$new <- v
head(dt)
#write.table(dt, file="diffCpGregion_hypo_v8.csv", sep=',', quote=FALSE)
#"diffCpGregion_hyper_v6.csv" has more column and I only wanted "new" column 
#with Promoter/exon/Intron/Others, so I reformed

copy_dt <- dt

#just adding the correcsponding position chr:start as row name in my copy_dt.
#did not shuffle the dataframe, just pasted the sequence of chr:start in the
#same order as it is present in the query_res_natozero
rownames(copy_dt)<-paste(rownames(query_hypo_res_natozero_v8))

#similarly copied the steps in dt
rownames(dt)<-paste(rownames(query_hypo_res_natozero_v8))

#changed the rowname from "new" to "region"
colnames(copy_dt)[4] = "region"
undesired <- c("prom","exon","intron")
copy_dt_regions <- copy_dt %>% select(-one_of(undesired))
write.table(copy_dt_regions, file="diffCpGregion_hypo_v8.csv", sep=',', quote=FALSE)

#Annoation file for row
annotation_row <- read.csv("diffCpGregion_hypo_v8.csv", header = T,
                           row.names = 1)

str(annotation_row)
str(CpGloci_hypo_v8)
str(annotation_col)

#color coding for the legends
annotation_colors <- list(condition = c('control'="cornflowerblue", 'f-ALS'="#FCBBA1",
                                        's-ALS'= "#FC9272"), sex = c('female' = "orange", 'male' = "brown"), 
                          region = c('Promoter'="blue", 'Exon'="purple", 'Intron'="pink",
                                     'Other' = "light yellow"))

pheatmap(CpGloci_hypo_v8, fontsize = 10, fontsize_row = 0.1, fontsize_col = 8, cluster_rows = T, 
         cluster_cols = T, annotation_row = annotation_row, annotation_col = annotation_col, cellwidth = 25, cellheight = 1.5, angle_row = 45,
         main = "Differential Methylation across ALS and Control samples", angle_col = 45,
         treeheight_row = 50, treeheight_col = 5, annotation_colors = annotation_colors, na_col = "grey",
         color = colorRampPalette(c("white", "yellow","green", "blue"))(100))




#Heatmap for hypermethylated regions

query_rows_v8<-paste(myDiff25p_hypo_tile500_8[1:792]$chr,":",myDiff25p_hypo_tile500_8[1:792]$start, "-", myDiff25p_hypo_tile500_8[1:792]$end, sep = "")
#dataframe for all DMLoci
all_loci_samples_v8 <- getData(tileMeth_2)
# placing row name as chr:start
rownames(all_loci_samples_v8)<- paste(all_loci_samples_v8$chr, ":", all_loci_samples_v8$start, "-", all_loci_samples_v8$end, sep = "")
# getting only those lines from all_loci_samples that are present in query_rows
query_rows_v8<- all_loci_samples_v8[rownames(all_loci_samples_v8) %in% query_rows_v8,]
#made a copy
query_hypo_res_natozero_v8 <- query_rows_v8
# replacing na with 0
#query_hypo_res_natozero_v8[is.na(query_hypo_res_natozero_v8)] <- 0
#dividing col having Cs to Coverage and storing in a columns meth1,2,3...so on
query_res_natozero_v8$meth1 <- query_rows_v8[,6]/query_rows_v8[,5]
query_res_natozero_v8$meth2 <- query_rows_v8[,9]/query_rows_v8[,8]
query_res_natozero_v8$meth3 <- query_rows_v8[,12]/query_rows_v8[,11]
query_res_natozero_v8$meth4 <- query_rows_v8[,15]/query_rows_v8[,14]
query_res_natozero_v8$meth5 <- query_rows_v8[,18]/query_rows_v8[,17]
query_res_natozero_v8$meth6 <- query_rows_v8[,21]/query_rows_v8[,20]
query_res_natozero_v8$meth7 <- query_rows_v8[,24]/query_rows_v8[,23]
query_res_natozero_v8$meth8 <- query_rows_v8[,27]/query_rows_v8[,26]
# again converted na to zero because above computation generated NA values
#query_res_natozero_v8[is.na(query_res_natozero_v8)] <- 0
#round the ratio to 2 decimal place; it just rewritten col 29 to 36 in the dataframe while removing other columns
query_res_natozero_v8<-round(query_res_natozero_v8[,(29:36)], digits = 2)
#renaming the columns in the dataframe just created above
colnames(query_res_natozero_v8)[1] = "ALS3"
colnames(query_res_natozero_v8)[2] = "ALS4"
colnames(query_res_natozero_v8)[3] = "ALS6"
colnames(query_res_natozero_v8)[4] = "ALS7"
colnames(query_res_natozero_v8)[5] = "ALS8"
colnames(query_res_natozero_v8)[6] = "ACRI22"
colnames(query_res_natozero_v8)[7] = "ACRI247"
colnames(query_res_natozero_v8)[8] = "Prost018"
#saves the dataframe to a file

write.table(query_res_natozero_v8, file="query_res_natozero_v8.tsv", sep='\t', quote=FALSE)
#Annotation_columns creation
condition<- c("f-ALS","s-ALS","s-ALS","s-ALS","s-ALS","control","control","control")
sex <- c("female","male","female","female","female","male","female","male")
#created a dataframe patients
data.frame(patients, stringsAsFactors = True)
patients <- data.frame(condition,sex)
rownames(patients)<- c("ALS3", "ALS4", "ALS6", "ALS7", "ALS8", "ACRI22", "ACRI247", "Prost018")
write.table(patients, file="ALSpatients.csv", sep=',', quote=FALSE)

library(pheatmap)
CpGloci_hyper_v8 <- read.csv("query_res_natozero_v8.tsv", header=T, sep='\t',
                             row.names = 1)

#annotation file for column
annotation_col <- read.csv("ALSpatients.csv", header = T, row.names = 1)

#creating annotation_row file
dt<- read.table("CpG_hypometh_tile500_v8.txt", header = T, stringsAsFactors = F, row.names = 1)
head(dt)
v <- ifelse(dt$prom == 1, "Promoter", "Other")
v[dt$exon == 1]<- "Exon"
v[dt$intron == 1]<- "Intron"
dt$new <- v
head(dt)
#write.table(dt, file="diffCpGregion_hypo_v8.csv", sep=',', quote=FALSE)
#"diffCpGregion_hyper_v6.csv" has more column and i only wanted "new" column 
#with Promoter/exon/Intron/Others, so i did

copy_dt <- dt

#just adding the correcsponding position chr:start as row name in my copy_dt.
#did not shuffle the dataframe, just pasted the sequence of chr:start in the
#same order as it is present in the query_res_natozero
rownames(copy_dt)<-paste(rownames(query_res_natozero_v8))

#similarly copied the action in dt
rownames(dt)<-paste(rownames(query_res_natozero_v8))

#changed the rowname from "new" to "region"
colnames(copy_dt)[4] = "region"
undesired <- c("prom","exon","intron")
copy_dt_regions <- copy_dt %>% select(-one_of(undesired))
write.table(copy_dt_regions, file="diffCpGregion_hypo_v8.csv", sep=',', quote=FALSE)

#Annoation file for row
annotation_row <- read.csv("diffCpGregion_hypo_v8.csv", header = T,
                           row.names = 1)

str(annotation_row)
str(CpGloci_hyper_v8)
str(annotation_col)

#color coding for the legends
annotation_colors <- list(condition = c('control'="cornflowerblue", 'f-ALS'="#FCBBA1",
                                        's-ALS'= "#FC9272"), sex = c('female' = "orange", 'male' = "brown"), 
                          region = c('Promoter'="blue", 'Exon'="purple", 'Intron'="pink",
                                     'Other' = "light yellow"))

pheatmap(CpGloci_hyper_v8, fontsize = 10, fontsize_row = 0.1, fontsize_col = 8, cluster_rows = T, 
         cluster_cols = T, annotation_row = annotation_row, annotation_col = annotation_col, cellwidth = 25, cellheight = 1.5, angle_row = 45,
         main = "Differential Methylation across ALS and Control samples", angle_col = 45,
         treeheight_row = 50, treeheight_col = 5, annotation_colors = annotation_colors, na_col = "grey",
         color = colorRampPalette(c("white", "yellow","green", "blue"))(100))(base)
