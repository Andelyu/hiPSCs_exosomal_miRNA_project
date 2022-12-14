
options(stringsAsFactors = F)

#########################################################################################
#------------------------------prepare stage----------------------------#
#########################################################################################
if(! require(ggplotify)){install.packages('ggplotify')}
if(! require(BiocManager)){install.packages('BiocManager')}
if(! require(DESeq2)){BiocManager::install("DESeq2")}
if(! require(edgeR)){BiocManager::install("edgeR")}
if(! require(ggplot2)){install.packages("ggplot2")}
if(! require(genefilter)){install.packages('genefilter')}
if(! require(pheatmap)){install.packages('pheatmap')}

install.packages("https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.2.tar.gz", repos=NULL, type="source")
install.packages("https://www.bioconductor.org/packages/release/bioc/src/contrib/DOSE_3.22.1.tar.gz", repos=NULL, type="source")
if(! require(clusterProfiler)){ install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/clusterProfiler_4.4.4.tar.gz", repos=NULL, type="source")}
if(! require(org.Hs.eg.db)){  BiocManager::install("org.Hs.eg.db",ask = F,update = F)}
if(! require(Rgraphviz)){ BiocManager::install('Rgraphviz')}
if(! require(VennDiagram)){ install.packages('VennDiagram')}
if(! require(ggplotify)){ install.packages('ggplotify')}
if(! require(multiMiR)){ install.packages('multiMiR')}

get_deinf_fun=function(data,coldata,outfile){
  # standard analysis
  rownames(coldata)=colnames(data)
  dds <- DESeqDataSetFromMatrix(data, coldata, ~ type)
  dds_filter <- dds[ rowSums(counts(dds))>6, ] 
  dds <- DESeq(dds_filter)
  res <- results(dds,contrast = c('type','hiPSC','Soma'))
  #table(res$padj<0.05)
  res_deseq <- res[order(res$padj),]
  diff_gene_deseq2 <- subset(res_deseq, padj<0.05 & abs(log2FoldChange)>1)
  summary(diff_gene_deseq2)
  res_diff_data <- merge(as.data.frame(diff_gene_deseq2),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
  rownames(res_diff_data)=res_diff_data$Row.names
  write.csv(res_diff_data,file =outfile,row.names = F,quote = F)
  return(DataFrame(res_diff_data))
}


setwd('./miRNA_data')
#----------------------get the target file----------------------#
miRNA_exp_count_df=read.csv('./raw_data/miRNAs_expressed_all_samples.csv',sep = "\t")

# calculate the mean value for the miRNA with multiple pre-miRNA
df_list=sapply(unique(miRNA_exp_count_df$X.miRNA), function(x){
  tmp=miRNA_exp_count_df[miRNA_exp_count_df$X.miRNA==x,]
  tmp=round(colMeans(tmp[,-1]))
  data.frame(tmp)
})
miRNA_exp_count_df=do.call(rbind,df_list)
colnames(miRNA_exp_count_df)=colnames(miRNA_exp_df)
rownames(miRNA_exp_count_df)=gsub(pattern = '.tmp',replacement = "",rownames(miRNA_exp_count_df))
head(miRNA_exp_count_df)


dir.create('reseults/res_pic')
dir.create('reseults/res_data')

#-------PCA analysis-----------#
#row: sampleids,col:variates
data_df <- data.frame(t(log2(cpm(miRNA_exp_count_df)+1)))
data_df <- data_df[,colMeans(data_df)>1]
sdata_df <- scale(data_df)
pdata_df<-prcomp(sdata_df)
screeplot(pdata_df,npcs = 20,type = 'line')
pca_pdata_df <- data.frame(pdata_df$x)
pca_pdata_df$sample<- rownames(pca_pdata_df);
pca_pdata_df$celltype <- 'Soma'; pca_pdata_df$celltype[grep('i',pca_pdata_df$sample)]='hiPSC'
pca_pdata_df$source='F'
pca_pdata_df$source[grep(pattern = 'U',x = pca_pdata_df$sample)]='U'
pca_pdata_df$source[grep(pattern = 'B',x = pca_pdata_df$sample)]='B'

p1=ggplot(data = pca_pdata_df,aes(x=PC1,y=PC2,color=source,shape=celltype))+geom_point(size=4)+
  geom_text(aes(label=sample,x = PC1,y = PC2+2) )+theme_classic()+theme(text = element_text(family = 'serif'))+scale_colour_brewer(palette = 'Set2')
p1
ggsave(plot = p1,filename = 'results/res_pic/Exsample_miRNA_pca_res.pdf',width = 6,height = 6,dpi = 300)
dev.off()

#----------------Sample Correlation ---------------#
filt_inf=miRNA_exp_count_df[rowMeans(miRNA_exp_count_df)>1,]
filt_inf<-log2(cpm(filt_inf)+1)
anno_df=pca_pdata_df[,c("celltype",'source')]

p=as.ggplot(pheatmap(cor(filt_inf),main = 'Correlation heatmap',annotation_col= anno_df,color =colorRampPalette(colors =c('gray','firebrick3'))(100)))+theme(text = element_text(family = 'serif'))         
ggsave(p,filename = 'results/res_pic/Exsample_miRNA_correlation_res.pdf',width = 6,height = 4)


#########################################################################################
#------------------------------DE miRNA analysis ----------------------#
#########################################################################################
#------------extract B/F/U group data ------#
exp_B=miRNA_exp_count_df[,1:6]; head(exp_B);
exp_F=miRNA_exp_count_df[,7:12];head(exp_F)
exp_U=miRNA_exp_count_df[,13:18];head(exp_U)
coldata=data.frame(type=factor(rep(c('hiPSC','Soma'),each=3),levels=c('hiPSC','Soma'))) # hiPSC vs Soma

outfile="results/res_data/Bips_B_exp_de_miRNA.csv"
de_exp_B=get_deinf_fun(exp_B,coldata,outfile)
dim(de_exp_B)

outfile='results/res_data/Fip_F_gene_exp_degs.csv'
de_exp_F=get_deinf_fun(exp_F,coldata,outfile)

outfile='results/res_data/exp_de_Uips_U.csv'
de_exp_U=get_deinf_fun(exp_U,coldata,outfile)

{
  #-----------------------------venn.diagram--------------------------#
  U_B_F_id=list(Uips_U=rownames(de_exp_U),Fips_F=rownames(de_exp_F),Bips_B=rownames(de_exp_B))
  all_dem_ids=c(rownames(de_exp_U),rownames(de_exp_F),rownames(de_exp_B))
  
  venn.plot=venn.diagram(U_B_F_id,
                         resolution = 500, imagetype = "tiff", alpha=c(0.5,0.5,0.5),
                         fill=c("#fbb4ae", "#b3cde3", "#ccebc5"), cat.fontface=10,fontfamily=3,
                         main="Intersection of de_genes identified among three cell lines",
                         main.cex = 1, main.fontface = 2, main.fontfamily = 3,width = 6,height = 6,
                         filename = NULL)
  
  grid.draw(venn.plot)
    
  #------------shared DE miRNA -----------#
  U_B_F=intersect(intersect(rownames(de_exp_U),rownames(de_exp_B)),rownames(de_exp_F))
  
  
  comm_de_exp=data.frame(Bips_B=de_exp_B[U_B_F,]$log2FoldChange,Uips_U=de_exp_U[U_B_F,]$log2FoldChange,Fips_F=de_exp_F[U_B_F,]$log2FoldChange)
  rownames(comm_de_exp)=U_B_F
  comm_de_exp=comm_de_exp[order(comm_de_exp$Bips_B,decreasing = T),] 
  pheatmap(comm_de_exp,main = 'shared DEM log2FoldChange heatmap',cluster_rows =T,show_rownames =T,color =colorRampPalette(colors = c('navy','white','firebrick3'))(100))
  write.csv(file='results/res_data/shared_dem_BFU.csv',comm_de_exp)
  annotation_color_list=list('source'=c('B'=palette('Set3')[1],'F'=palette('Set3')[2],'U'=palette('Set3')[3]),"celltype"=c('hiPSC'=palette('Set3')[4],'Soma'=palette('Set3')[5]))
  p= as.ggplot(pheatmap(log2(cpm(miRNA_exp_ncount_df)[U_B_F,]+1),cluster_cols = F,color =colorRampPalette(colors = c('navy','white','firebrick3'))(100),
                        annotation_col = data.frame(row.names =rownames(anno_df) ,celltype=factor(anno_df$celltype),source=factor(anno_df$source)),annotation_colors  =annotation_color_list ))
  p=p+theme(text = element_text(family = 'serif') )
  p
  ggsave(p,filename = 'results/res_pic/shared_DE_miRNA_exp_heatmap.pdf',width = 8,height = 12,dpi = 300)
  
}

#########################################################################################
#--------------------------DE analysis of  RNA-Seq ------------------------#
#########################################################################################
gene_inf=read.csv('raw_data/allsample_gene_count_df.csv',header = T)
rownames(gene_inf)=gene_inf$X;gene_inf=gene_inf[,-1]
gene_inf=gene_inf[,sort(colnames(gene_inf))]
head(gene_inf)
nr_gene_inf_df=data.frame(log2(cpm(gene_inf)+1))
nr_gene_inf_df[,'gene']=as.character(mapIds(org.Hs.eg.db,keys = rownames(nr_gene_inf_df),column = 'SYMBOL',keytype = 'ENSEMBL'))

# hiPSC/soma marker genes expression
hiPSC_markers=c('DPPA2','GDF3','TERT','NANOG','ZFP42','FGF4','LIN28A','DPPA4','POU5F1','SOX2')
differentiate_markers=c('AFP','NR2F2','ANPEP','SOX17')
tpm_gene_df=read.csv('raw_data/allsample_gene_TPM_df.csv',header = T)
rownames(tpm_gene_df)=tpm_gene_df$X;tpm_gene_df=tpm_gene_df[,-1]
tpm_gene_df$gene=nr_gene_inf_df[rownames(tpm_gene_df),'gene']
  
marker_df=tpm_gene_df[match(c(hiPSC_markers,differentiate_markers),tpm_gene_df$gene),]
rownames(marker_df)=marker_df$gene;marker_df=marker_df[,-18]
row_anno=data.frame(row.names =rownames(marker_df),type=c(rep('pluripotency',length(hiPSC_markers)),rep('differentiation',length(differentiate_markers))))
p=as.ggplot(pheatmap(marker_df,cluster_rows = F,scale = 'row',annotation_row = row_anno,color =colorRampPalette(colors = c('navy','white','firebrick3'))(100)))+theme(text = element_text(family = 'serif'))
ggsave(p,filename = 'results/res_pic/marker_gene_exp_heatmap.pdf',width = 6,height = 4)


# pca anlysis
data_df <- t(nr_gene_inf_df[,-18])
data_df <- data_df[,colMeans(data_df)>1]
sdata_df <- scale(data_df)
pdata_df<-prcomp(sdata_df)
screeplot(pdata_df,npcs = 20,type = 'line')
pca_pdata_df <- data.frame(pdata_df$x)
pca_pdata_df$sample<- rownames(pca_pdata_df);
pca_pdata_df$celltype <- 'Soma'; pca_pdata_df$celltype[grep('i',pca_pdata_df$sample)]='hiPSC'
pca_pdata_df$source='F'
pca_pdata_df$source[grep(pattern = 'U',x = pca_pdata_df$sample)]='U'
pca_pdata_df$source[grep(pattern = 'B',x = pca_pdata_df$sample)]='B'

p1=ggplot(data = pca_pdata_df,aes(x=PC1,y=PC2,color=source,shape=celltype))+geom_point(size=4)+
   theme_classic()+theme(text = element_text(family = 'serif'))+scale_colour_brewer(palette = 'Set2')
p1
ggsave(plot = p1,filename = 'results/res_pic/Exsample_gene_pca_res.pdf',width = 6,height = 6,dpi = 300)

#----------------Sample Correlation ---------------#
filt_inf=data.frame(nr_gene_inf_df[rowMeans(nr_gene_inf_df[,-18])>1,-18])
anno_df=pca_pdata_df[,c("celltype",'source')]
p=as.ggplot(pheatmap(cor(filt_inf),main = 'Correlation heatmap',annotation_col= anno_df,color =colorRampPalette(colors =c('gray','firebrick3'))(100)))+theme(text = element_text(family = 'serif'))
p         
ggsave(p,filename = 'results/res_pic/sample_gene_correlation_res.pdf',width = 6,height = 4)


coldata=data.frame(type=as.factor(c(rep('hiPSC',3),rep('Soma',3),rep('hiPSC',2),rep('Soma',3),rep('hiPS',3),rep('Soma',3))))
#-----------------------------group DE ---------------------------------#
outfile='results/res_data/gene_exp_deg_basedon_ips_somatic.csv'
deg_basedon_celltype=get_deinf_fun(gene_inf,coldata =coldata,outfile = outfile )
pheatmap(log2(data.frame(deg_basedon_celltype[1:50,8:24])+1),show_rownames = F)

#-----------------------------group type DE ---------------------------------#
coldata=data.frame(type=factor(rep(c('hiPSC','Soma'),each=3),levels=c('hiPSC','Soma')),row.names =colnames(gene_inf[,1:6]) )
outfile='results/res_data/Bip_B_gene_exp_degs.csv'                                     
deg_exp_B=get_deinf_fun(gene_inf[,1:6],coldata =coldata, outfile = outfile)

coldata=data.frame(type=factor(c(rep(c('hiPSC','Soma'),each=2),'Soma'),levels=c('hiPSC','Soma')),row.names =colnames(gene_inf[,7:11]) )
outfile='results/res_data/Fips_F_gene_exp_degs.csv'                                      
deg_exp_F=get_deinf_fun(gene_inf[,7:11],coldata =coldata, outfile = outfile)

coldata=data.frame(type=factor(rep(c('hiPSC','Soma'),each=3),levels=c('hiPSC','Soma')),row.names =colnames(gene_inf[,12:17]) )
outfile='results/res_data/Uips_U_gene_exp_degs.csv'
deg_exp_U=get_deinf_fun(gene_inf[,c(12:17)],coldata =coldata, outfile = outfile)

{
  #-----------------------------venn.diagram--------------------------#
  
  U_B_F_id=list(Uips_U=rownames(deg_exp_U),Fips_F=rownames(deg_exp_F),Bips_B=rownames(deg_exp_B))
  all_dem_ids=c(rownames(deg_exp_U),rownames(deg_exp_F),rownames(deg_exp_B))
  length(uniq_all_dem_ids)

  p=venn.diagram(U_B_F_id,
               resolution = 300, imagetype = "tiff", alpha=c(0.5,0.5,0.5),
               fill=c("#fbb4ae", "#b3cde3", "#ccebc5"), cat.fontface=10,fontfamily=3,
               main.cex = 1, main.fontface = 2, main.fontfamily = 3,
               filename = NULL)
  grid.draw(p)
  # width=6,heigth=6,file=shared_DEGs_venn.pdf
  
  shared_Up_ids=list(Uips_U=rownames(deg_exp_U[deg_exp_U$log2FoldChange>0,]),Fips_F=rownames(deg_exp_F[deg_exp_F$log2FoldChange>0,]),Bips_B=rownames(deg_exp_B[deg_exp_B$log2FoldChange>0,]))
  p1=venn.diagram(shared_Up_ids,
                 resolution = 300, imagetype = "tiff", alpha=c(0.5,0.5,0.5),main = 'Shared Up DEGs',
                 fill=c("#fbb4ae", "#b3cde3", "#ccebc5"), cat.fontface=10,fontfamily=3,
                 main.cex = 1, main.fontface = 2, main.fontfamily = 3,
                 filename = NULL)
  dev.off()
  grid.draw(p1) # width=6,heigth=6,file=shared_Up_DEGs_venn.pdf
  
  shared_Dn_ids=list(Uips_U=rownames(deg_exp_U[deg_exp_U$log2FoldChange<0,]),Fips_F=rownames(deg_exp_F[deg_exp_F$log2FoldChange<0,]),Bips_B=rownames(deg_exp_B[deg_exp_B$log2FoldChange<0,]))
  p2=venn.diagram(shared_Dn_ids,
                  resolution = 300, imagetype = "tiff", alpha=c(0.5,0.5,0.5),main = 'Shared Dn DEGs',
                  fill=c("#fbb4ae", "#b3cde3", "#ccebc5"), cat.fontface=10,fontfamily=3,
                  main.cex = 1, main.fontface = 2, main.fontfamily = 3,
                  filename = NULL)
  dev.off()
  grid.draw(p2) # width=6,heigth=6,shared_Dn_DEGs_venn.pdf
  
  #  enrichment analysis :GO/KEGG
  group_deg_df=rbind(rbind(deg_exp_B[,c('log2FoldChange','Row.names')],deg_exp_F[,c('log2FoldChange','Row.names')]),deg_exp_U[,c('log2FoldChange','Row.names')])
  rownames(group_deg_df)=NULL;group_deg_df=data.frame(group_deg_df)
  group_deg_df$group=c(rep('Bips_B',dim(deg_exp_B)[1]),rep('Fips_F',dim(deg_exp_F)[1]),rep('Uips_U',dim(deg_exp_U)[1]))
  group_deg_df$entrzid= mapIds(x = org.Hs.eg.db,keys =group_deg_df$Row.names,keytype = 'ENSEMBL',column = 'ENTREZID' )
  group_deg_df=group_deg_df[order(group_deg_df$log2FoldChange,decreasing = T),]
  all_deg_gsego_BP=compareCluster(geneClusters = entrzid|log2FoldChange~group,data = group_deg_df,fun ='gseGO' ,ont = "BP",OrgDb='org.Hs.eg.db')

  all_deg_gsego_BP_top= data.frame(rbind(all_deg_gsego_BP@compareClusterResult[1:20,],all_deg_gsego_BP@compareClusterResult[all_deg_gsego_BP@compareClusterResult$Cluster=='Fips_F',][1:20,],
                                         all_deg_gsego_BP@compareClusterResult[all_deg_gsego_BP@compareClusterResult$Cluster=='Uips_U',][1:20,]))
  all_deg_gsego_BP_top$Gene.Ratio= round(as.numeric(sapply(all_deg_gsego_BP_top$core_enrichment,function(x){length(strsplit(x = x,split = '/')[[1]])}))/all_deg_gsego_BP_top$setSize,2)
  all_deg_gsego_BP_top$Description=factor(all_deg_gsego_BP_top$Description,levels = unique(c(all_deg_gsego_BP_top[all_deg_gsego_BP_top$Cluster=='Bips_B','Description'],
                                         all_deg_gsego_BP_top[all_deg_gsego_BP_top$Cluster=='Fips_F' & (! all_deg_gsego_BP_top$Description %in% all_deg_gsego_BP_top[all_deg_gsego_BP_top$group=='Uips_U','Description']),'Description'],
                                         all_deg_gsego_BP_top$Description )))                                                                                       
  p=ggplot(all_deg_gsego_BP_top,aes(x=group,y=Description,color=NES,size=Gene.Ratio))+geom_point()+scale_color_gradient(low = 'blue',high = 'red')+theme_classic()+theme(text = element_text(family = 'serif'));p
  ggsave(p,filename = 'results/res_pic/de_gene_gseGOBP.pdf',width = 10,height = 10,dpi = 300)
  write.csv(data.frame(all_deg_gsego_BP@compareClusterResult),file = 'results/res_data/all_deg_gsego_BP.csv')
  
  # gseKEGG  analysis 
  all_deg_gseKEGG=compareCluster(geneClusters = entrzid|log2FoldChange~group,data = group_deg_df,fun ='gseKEGG' ,organism="hsa",eps=0)
  all_deg_gseKEGG_top= data.frame(rbind(all_deg_gseKEGG@compareClusterResult[1:20,],all_deg_gseKEGG@compareClusterResult[all_deg_gseKEGG@compareClusterResult$Cluster=='Fips_F',][1:20,],
                                        all_deg_gseKEGG@compareClusterResult[all_deg_gseKEGG@compareClusterResult$Cluster=='Uips_U',][1:20,]))
  all_deg_gseKEGG_top$Gene.Ratio= round(as.numeric(sapply(all_deg_gseKEGG_top$core_enrichment,function(x){length(strsplit(x = x,split = '/')[[1]])}))/all_deg_gseKEGG_top$setSize,2)
  all_deg_gseKEGG_top$Description=factor(all_deg_gseKEGG_top$Description,levels = unique(c(all_deg_gseKEGG_top[all_deg_gseKEGG_top$Cluster=='Bips_B','Description'],
                                                                                           all_deg_gseKEGG_top[all_deg_gseKEGG_top$Cluster=='Fips_F' & (! all_deg_gseKEGG_top$Description %in% all_deg_gseKEGG_top[all_deg_gseKEGG_top$group=='Uips_U','Description']),'Description'],
                                                                                           all_deg_gseKEGG_top$Description )))                                                                                       
  p=ggplot(all_deg_gseKEGG_top,aes(x=group,y=Description,color=NES,size=Gene.Ratio))+geom_point()+scale_color_gradient(low = 'blue',high = 'red')+theme_classic()+theme(text = element_text(family = 'serif'));p
  ggsave(p,filename = 'results/res_pic/de_gene_gseKEGG.pdf',width = 6,height = 8,dpi = 300)
  write.csv(data.frame(all_deg_gseKEGG@compareClusterResult),file = 'results/res_data/all_deg_gseKEGG.csv')
    
 }


#########################################################################################
#----------------------------------miRNA targeted genes-------------------------#
#########################################################################################
all_mirids=unique(c(rownames(de_exp_U),rownames(de_exp_F),rownames(de_exp_B)))
all_geneids=c(rownames(deg_exp_U),rownames(deg_exp_F),rownames(deg_exp_B))
all_geneids_sym=as.character(mapIds(org.Hs.eg.db,keys = unique(all_geneids), keytype = 'ENSEMBL',column = 'SYMBOL'))

all_up_mirids_list=list('Bips_B'=rownames(de_exp_B[de_exp_B$log2FoldChange > 0,]), 'Fips_F'=rownames(de_exp_F[de_exp_F$log2FoldChange > 0,]), 'Uips_U'=rownames(de_exp_U[de_exp_U$log2FoldChange > 0,]))
all_dn_mirids_list=list('Bips_B'=rownames(de_exp_B[de_exp_B$log2FoldChange < 0,]), 'Fips_F'=rownames(de_exp_F[de_exp_F$log2FoldChange < 0,]), 'Uips_U'=rownames(de_exp_U[de_exp_U$log2FoldChange < 0,]))

all_up_deg_list=list('Bips_B'=rownames(deg_exp_B[deg_exp_U$log2FoldChange > 0,]),'Fips_F'=rownames(deg_exp_F[deg_exp_U$log2FoldChange > 0,]),'Uips_U'=rownames(deg_exp_U[deg_exp_U$log2FoldChange > 0,]))
all_dn_deg_list=list('Bips_B'=rownames(deg_exp_B[deg_exp_U$log2FoldChange < 0,]),'Fips_F'=rownames(deg_exp_F[deg_exp_U$log2FoldChange < 0,]),'Uips_U'=rownames(deg_exp_U[deg_exp_U$log2FoldChange < 0,]))

Up_dem2deg_validated_results=list();Up_dem2deg_predicted_results=list()
Dn_dem2deg_validated_results=list();Dn_dem2deg_predicted_results=list()

# 比较消耗时间，注意及时保存
library(multiMiR )
for(i in 1:3){
  Up_dem2deg_validated_results[[i]]=get_multimir(org = 'hsa',mirna = all_up_mirids_list[[i]],target =all_dn_deg_list[[i]],table   = 'validated',summary = TRUE)
  Up_dem2deg_predicted_results[[i]]=get_multimir(org = 'hsa',mirna = all_up_mirids_list[[i]],target =all_dn_deg_list[[i]],table   = 'predicted',summary = TRUE)
  Dn_dem2deg_validated_results[[i]]=get_multimir(org = 'hsa',mirna = all_dn_mirids_list[[i]],target =all_up_deg_list[[i]],table   = 'validated',summary = TRUE)
  Dn_dem2deg_predicted_results[[i]]=get_multimir(org = 'hsa',mirna = all_dn_mirids_list[[i]],target =all_up_deg_list[[i]],table   = 'predicted',summary = TRUE)
}

Up_dem2deg_validated_results_df=rbind(Up_dem2deg_validated_results[[1]]@summary,Up_dem2deg_validated_results[[2]]@summary,Up_dem2deg_validated_results[[3]]@summary)
Up_dem2deg_validated_results_df$group=c(rep('B_up',dim(Up_dem2deg_validated_results[[1]]@summary)[1]),rep('F_up',dim(Up_dem2deg_validated_results[[2]]@summary)[1]),rep('U_up',dim(Up_dem2deg_validated_results[[3]]@summary)[1]))
Dn_dem2deg_validated_results_df=rbind(Dn_dem2deg_validated_results[[1]]@summary,Dn_dem2deg_validated_results[[2]]@summary,Dn_dem2deg_validated_results[[3]]@summary)
Dn_dem2deg_validated_results_df$group=c(rep('B_dn',dim(Dn_dem2deg_validated_results[[1]]@summary)[1]),rep('F_dn',dim(Dn_dem2deg_validated_results[[2]]@summary)[1]),rep('U_dn',dim(Dn_dem2deg_validated_results[[3]]@summary)[1]))

dem2deg_validated_results_df=rbind(Up_dem2deg_validated_results_df,Dn_dem2deg_validated_results_df)

Up_dem2deg_predicted_results_df=rbind(Up_dem2deg_predicted_results[[1]]@summary,Up_dem2deg_predicted_results[[2]]@summary,Up_dem2deg_predicted_results[[3]]@summary)
Up_dem2deg_predicted_results_df$group=c(rep('B_up',dim(Up_dem2deg_predicted_results[[1]]@summary)[1]),rep('F_up',dim(Up_dem2deg_predicted_results[[2]]@summary)[1]),rep('U_up',dim(Up_dem2deg_predicted_results[[3]]@summary)[1]))
Dn_dem2deg_predicted_results_df=rbind(Dn_dem2deg_predicted_results[[1]]@summary,Dn_dem2deg_predicted_results[[2]]@summary,Dn_dem2deg_predicted_results[[3]]@summary)
Dn_dem2deg_predicted_results_df$group=c(rep('B_dn',dim(Dn_dem2deg_predicted_results[[1]]@summary)[1]),rep('F_dn',dim(Dn_dem2deg_predicted_results[[2]]@summary)[1]),rep('U_dn',dim(Dn_dem2deg_predicted_results[[3]]@summary)[1]))

dem2deg_predicted_results_df=rbind(Up_dem2deg_predicted_results_df,Dn_dem2deg_predicted_results_df)
## no less 4: >=4 databases
dem2deg_predicted_results_df=dem2deg_predicted_results_df[dem2deg_predicted_results_df$all.sum>3,]

write.csv(dem2deg_validated_results_df,file = 'results/res_data/dem2deg_validated_results_df.csv')
write.csv(dem2deg_predicted_results_df,file = 'results/res_data/dem2deg_predicted_results_df.csv')

table(dem2deg_validated_results_df$group)
table(dem2deg_predicted_results_df$group)

dem2deg_results_df=rbind(dem2deg_validated_results_df[,c('mature_mirna_id','target_ensembl','group')],dem2deg_predicted_results_df[,c('mature_mirna_id','target_ensembl','group')])
dem2deg_results_df=unique(dem2deg_results_df)
dem2deg_results_df$group=factor(dem2deg_results_df$group,levels = unique(dem2deg_results_df$group))
dem2deg_results_df$type='Bips_B';dem2deg_results_df$type[grep(dem2deg_results_df$group,pattern = 'F')]='Fips_F';dem2deg_results_df$type[grep(dem2deg_results_df$group,pattern = 'U')]='Uips_U'
dem2deg_exp_df=rbind(group_deg_df[group_deg_df$group=='Bips_B' & group_deg_df$Row.names %in% dem2deg_results_df[dem2deg_results_df$type=='Bips_B','target_ensembl'],],
                     group_deg_df[group_deg_df$group=='Fips_F' & group_deg_df$Row.names %in% dem2deg_results_df[dem2deg_results_df$type=='Fips_F','target_ensembl'],],
                     group_deg_df[group_deg_df$group=='Uips_U' & group_deg_df$Row.names %in% dem2deg_results_df[dem2deg_results_df$type=='Uips_U','target_ensembl'],]
                     )

dem2deg_gsego_BP=compareCluster(geneClusters = entrzid|log2FoldChange~group,data = dem2deg_exp_df,fun ='gseGO' ,ont = "BP",OrgDb='org.Hs.eg.db',eps=0,nPermSimple = 10000)
dem2deg_gsego_BP_top= data.frame(rbind(dem2deg_gsego_BP@compareClusterResult[1:20,],dem2deg_gsego_BP@compareClusterResult[dem2deg_gsego_BP@compareClusterResult$Cluster=='Fips_F',][1:20,],
                                       dem2deg_gsego_BP@compareClusterResult[dem2deg_gsego_BP@compareClusterResult$Cluster=='Uips_U',][1:20,]))

dem2deg_gsego_BP_top$Gene.Ratio= round(as.numeric(sapply(dem2deg_gsego_BP_top$core_enrichment,function(x){length(strsplit(x = x,split = '/')[[1]])}))/dem2deg_gsego_BP_top$setSize,2)
dem2deg_gsego_BP_top$Description=factor(dem2deg_gsego_BP_top$Description,levels = unique(c(dem2deg_gsego_BP_top[dem2deg_gsego_BP_top$Cluster=='Bips_B','Description'],
                                                                                           dem2deg_gsego_BP_top[dem2deg_gsego_BP_top$Cluster=='Fips_F' & (! dem2deg_gsego_BP_top$Description %in% dem2deg_gsego_BP_top[dem2deg_gsego_BP_top$group=='Uips_U','Description']),'Description'],
                                                                                           dem2deg_gsego_BP_top$Description )))                                                                                       
p=ggplot(dem2deg_gsego_BP_top,aes(x=group,y=Description,color=NES,size=Gene.Ratio))+geom_point()+scale_color_gradient(low = 'blue',high = 'red')+theme_classic()+theme(text = element_text(family = 'serif'));p
ggsave(p,filename = 'results/res_pic/dem2deg_gseGOBP.pdf',width = 6,height = 8,dpi = 300)
write.csv(data.frame(dem2deg_gsego_BP@compareClusterResult),file = 'results/res_data/dem2deg_gsego_BP.csv')


dem2deg_gseKEGG=compareCluster(geneClusters = entrzid|log2FoldChange~group,data = dem2deg_exp_df,fun ='gseKEGG' ,organism='hsa',eps=0)
dem2deg_gseKEGG_top= data.frame(rbind(dem2deg_gseKEGG@compareClusterResult[1:20,],dem2deg_gseKEGG@compareClusterResult[dem2deg_gseKEGG@compareClusterResult$Cluster=='Fips_F',][1:20,],
                                      dem2deg_gseKEGG@compareClusterResult[dem2deg_gseKEGG@compareClusterResult$Cluster=='Uips_U',][1:20,]))

dem2deg_gseKEGG_top$Gene.Ratio= round(as.numeric(sapply(dem2deg_gseKEGG_top$core_enrichment,function(x){length(strsplit(x = x,split = '/')[[1]])}))/dem2deg_gseKEGG_top$setSize,2)
dem2deg_gseKEGG_top$Description=factor(dem2deg_gseKEGG_top$Description,levels = unique(c(dem2deg_gseKEGG_top[dem2deg_gseKEGG_top$Cluster=='Bips_B','Description'],
                                                                                          dem2deg_gseKEGG_top[dem2deg_gseKEGG_top$Cluster=='Fips_F' & (! dem2deg_gseKEGG_top$Description %in% dem2deg_gseKEGG_top[dem2deg_gseKEGG_top$group=='Uips_U','Description']),'Description'],
                                                                                          dem2deg_gseKEGG_top$Description )))                                                                                       
p=ggplot(dem2deg_gseKEGG_top,aes(x=group,y=Description,color=NES,size=Gene.Ratio))+geom_point()+scale_color_gradient(low = 'blue',high = 'red')+theme_classic()+theme(text = element_text(family = 'serif'));p
ggsave(p,filename = 'results/res_pic/dem2deg_gseKEGG.pdf',width = 6,height = 8,dpi = 300)
write.csv(data.frame(dem2deg_gseKEGG@compareClusterResult),file = 'results/res_data/dem2deg_gseKEGG.csv')


#  calculate the shared pathways of DEGs and DEMs2DEGs 
shared_go_list=list();go_number_list=list()
shared_kegg_list=list();kegg_number_list=list()
for (type in unique(dem2deg_gseKEGG@compareClusterResult$group)){
  DEGs_go_pathwaysid=all_deg_gsego_BP@compareClusterResult[all_deg_gsego_BP@compareClusterResult$group==type,'ID']
  targeted_DEGs_go_pathwaysid=dem2deg_gsego_BP@compareClusterResult[dem2deg_gsego_BP@compareClusterResult$group==type,'ID']
  shared_go_list[[type]]=length(intersect(DEGs_go_pathwaysid,targeted_DEGs_go_pathwaysid))
  go_number_list[[type]]=c(length(DEGs_go_pathwaysid),length(targeted_DEGs_go_pathwaysid))
  
  DEGs_kegg_pathwaysid=all_deg_gseKEGG@compareClusterResult[all_deg_gseKEGG@compareClusterResult$group==type,'ID']
  targeted_DEGs_kegg_pathwaysid=dem2deg_gseKEGG@compareClusterResult[dem2deg_gseKEGG@compareClusterResult$group==type,'ID']
  shared_kegg_list[[type]]=length(intersect(DEGs_kegg_pathwaysid,targeted_DEGs_kegg_pathwaysid))
  kegg_number_list[[type]]=c(length(DEGs_kegg_pathwaysid),length(targeted_DEGs_kegg_pathwaysid))
}

go_summary_df=rbind(data.frame(shared_go_list),data.frame(go_number_list))
rownames(go_summary_df)=c('shared','DEGs','DEMs2DEGs')
go_summary_df
kegg_summary_df=rbind(data.frame(shared_kegg_list),data.frame(kegg_number_list))
rownames(kegg_summary_df)=c('shared','DEGs','DEMs2DEGs')
kegg_summary_df
write.csv(kegg_summary_df,file='results/res_data/kegg_summary_df.csv',quote=F)


