
library(tidyverse)
library(dplyr)


setwd('~/Desktop/projects/phase_separation/')


#10 gene clusters abtained from enrichmentmap results:
#chromatin:
#	chromatin organization
#	chromatin organization regulation
#	histone modification
#	histone modification regulation
#	DNA binding.
#transcription:
#	mRNA regulation
#	RNAP
#	RNA splicing
#cytoskeleton:
#	microtubule/actin
#	spindle


#1. AA preferences of GSEA enriched clusters: ----------

#step 1: clean up the gene list whose ratio of idr is over 0.3:
#
#awk '{print $1 "\t"$2}' DNA_binding.txt | grep -v Description | awk '{print $0 "\tDNA_binding"}' > cluster_cleanup.txt
#awk '{print $1 "\t"$2}' RNAP.txt| grep -v Description | awk '{print $0 "\tRNAP"}' >> cluster_cleanup.txt
#awk '{print $1 "\t"$2}' RNA_splice.txt | grep -v Description | awk '{print $0 "\tRNA_splice"}' >> cluster_cleanup.txt
#awk '{print $1 "\t"$2}'  chromatin_organization.txt | grep -v Description | awk '{print $0 "\tchromatin_organization"}' >> cluster_cleanup.txt
#awk '{print $1 "\t"$2}'  chromatin_organization_regulation.txt | grep -v Description | awk '{print $0 "\tchromatin_organization_regulation"}' >> cluster_cleanup.txt
#awk '{print $1 "\t"$2}'  histone_modifications.txt | grep -v Description | awk '{print $0 "\thistone_modification"}' >> cluster_cleanup.txt
#awk '{print $1 "\t"$2}'  histone_modifications_regulation.txt | grep -v Description | awk '{print $0 "\thistone_modification_regulation"}' >> cluster_cleanup.txt
#awk '{print $1 "\t"$2}'  mRNA_regulation.txt | grep -v Description | awk '{print $0 "\tmRNA_regulation"}' >> cluster_cleanup.txt
#awk '{print $1 "\t"$2}'  microtubule_actin.txt | grep -v Description | awk '{print $0 "\tmicrotubule_actin"}' >> cluster_cleanup.txt
#awk '{print $1 "\t"$2}'  spindle.txt | grep -v Description | awk '{print $0 "\tspindle"}' >> cluster_cleanup.txt


# step2. calculate the AA distribution of all IDR regions of given gene list:

cluster_gs <- read_tsv('./results/gsea/enrichment_map/cluster_gene_list/cluster_cleanup.txt')


cluster_names <- unique(cluster_gs$enrichmap_cluster)

#subtract the AAs of IDR regions:
IDR_list <- list() #store 10 clusters


for (i in cluster_names) {
	print(i)
	tmp <- list() #store AAs of each proteins in the cluster
	gene_ls <- cluster_gs[cluster_gs$enrichmap_cluster == i & 
						  	cluster_gs$norm_mean_value >= 0.3 & 
						  	!is.na(cluster_gs$norm_mean_value), ]$gene
	for (J in gene_ls) {
		file_path <- sprintf('./results/iupred2a_results/human/%s.iupred2a.txt', J)
		gene_score <- read.table(file_path)
	
#		AAs_in_IDR <- gene_score[mean(gene_score$iupred2_score + gene_score$achore2_score) >= 0.5,]$amino_acid
		AAs_in_IDR <- gene_score[mean(gene_score$V3 + gene_score$V4) >= 0.5,]$V2
		tmp[[J]] <- AAs_in_IDR
	}
	IDR_list[[i]] <- tmp
}


IDR_list$DNA_binding %>% length()

#check the gene numbers:
sapply(IDR_list, length)
filter(cluster_gs, norm_mean_value >= 0.3 & !is.na(norm_mean_value)) %>%
	group_by(enrichmap_cluster) %>% count()




AA_tibbles <- IDR_list$DNA_binding %>% unlist() %>% table() %>% as_tibble()
names(AA_tibbles) <- c('AAs', 'DNA_binding')

for (i in cluster_names[2:10]) {
	print(i)
	counts <- IDR_list[[i]] %>% unlist() %>% table() %>% as_tibble()
	names(counts) <- c('AAs', i)
	AA_tibbles <- left_join(AA_tibbles, counts, by = 'AAs')
}

AA_tibbles

AA_tibbles_ratio <- AA_tibbles[2:11] %>% apply(2, function(x){return(x/sum(x))})
rownames(AA_tibbles_ratio) <- AA_tibbles$AAs



pheatmap::pheatmap(AA_tibbles_ratio, cluster_rows = F, cluster_cols = F)

pheatmap::pheatmap(AA_tibbles_ratio)
pheatmap::pheatmap(AA_tibbles_ratio %>% t(), las.anno =2)

pheatmap::pheatmap(AA_tibbles_ratio %>% t(), 
				   color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

#AA charges:	
#positive: R(Arg)	K(Lys)	H(His)
#negative: D(Asp)	E(Glu)
#polar: S(Ser)	T(Thr)	N(Asn)	Q(Gln)
#special: C(Cys)	G(Gly)	P(Pro)
#Hydrophobic: A(Ala)	V(Val)	I(Ile)	L(Leu)	M(Met)	F(Phe)	Y(Tyr)	W(Trp)



#2. IDR ratio and protein length --------------------------------------------
library(Biostrings)
cluster_gs <- read_tsv('./results/gsea/enrichment_map/cluster_gene_list/cluster_cleanup.txt')
protein_fa <- readDNAStringSet('./results/fa/human/all_protein_coding_genes_without_asterix_aa_sequence.fa')

protein_length <- tibble(gene = names(protein_fa), length = width(protein_fa))


#add gene length info:
cluster_gs <- left_join(cluster_gs, protein_length)

layout(mat = matrix(1))

plot(cluster_gs$length, cluster_gs$norm_mean_value, pch = 19, cex =2,
	 col = as.factor(cluster_gs$enrichmap_cluster),
	 xlim = c(0, 2000), ylim = c(0, 1))

legend(1400, 1, unique(cluster_gs$enrichmap_cluster), col = 1:10)





layout(mat = matrix(1:10, nrow = 5))

for (i in unique(cluster_gs$enrichmap_cluster)) {
	tmp <- cluster_gs[cluster_gs$enrichmap_cluster ==i, ]
	plot(tmp$length, tmp$norm_mean_value, pch = 19,
		 col = as.factor(tmp$enrichmap_cluster),
		 xlim = c(0, 2000), ylim = c(0, 1),
		 main = i)
	
}



#3. motif enrichment analysis --------------------------

cluster_gs <- read_tsv('./results/gsea/enrichment_map/cluster_gene_list/cluster_cleanup.txt')

chromatin_gene_list <- filter(cluster_gs, 
		norm_mean_value >=0.3 & 
		!(is.na(norm_mean_value)) &
		(enrichmap_cluster %in% c('chromatin_organization', 
								  'chromatin_organization_regulation', 
								  'histone_modification', 
								  'histone_modification_regulation', 
								  'DNA_binding')))

chromatin_gene_list$enrichmap_cluster %>% table()

#write_tsv(chromatin_gene_list[, 'gene'], col_names = F,
#		  './results/gsea/enrichment_map/cluster_gene_list/chromatin_gene_list.txt')
#


cluster_gs$enrichmap_cluster %>% table()

transcription_gene_list <- filter(cluster_gs, 
							  norm_mean_value >=0.3 & 
							  	!(is.na(norm_mean_value)) &
							  	(enrichmap_cluster %in% c('mRNA_regulation', 
							  							  'RNA_splice', 
							  							  'RNAP')))
transcription_gene_list$enrichmap_cluster %>% table()
write_tsv(transcription_gene_list[, 'gene'], col_names = F,
		  './results/gsea/enrichment_map/cluster_gene_list/transcription_gene_list.txt')



cytoskeleton_gene_list <- filter(cluster_gs, 
								  norm_mean_value >=0.3 & 
								  	!(is.na(norm_mean_value)) &
								  	(enrichmap_cluster %in% c('microtubule_actin', 
								  							  'spindle')))
cytoskeleton_gene_list$enrichmap_cluster %>% table()
write_tsv(cytoskeleton_gene_list[, 'gene'], col_names = F,
		  './results/gsea/enrichment_map/cluster_gene_list/cytoskeleton_gene_list.txt')


unique(chromatin_gene_list$gene) %in% unique(transcription_gene_list$gene) %>% table()


#calculate the motifs occurance time:
