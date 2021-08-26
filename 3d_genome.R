
library(tidyverse)
library(dplyr)
library(magrittr)
setwd('~/Desktop/projects/phase_separation/')
library(RColorBrewer)
display.brewer.all()
palette(c(brewer.pal(12, 'Set3'), brewer.pal(8, 'Pastel1')))



# 1. gsea analysis with enrichmentmap gene sets --------------------------------
#result path:
#	/Users/xuejia/gsea_home/output/phase_separation/phase_sepration_10w_permutation.GseaPreranked.1577773895962
#parameters:
#	10000 permutations


# 2. select the 3d genome related gene sets enriched in gsea ----------------------
ps_enrichmentmap <- read_csv('./results/gsea/gsea_report_for_na_pos_1577773895962.csv')

#only keep GOBP(biological process) terms:
ps_enrichmentmap %<>% filter(grepl('GOBP', NAME)) %>%
					mutate(NAME = tolower(NAME) %>% 
							 	strsplit(., split = '%') %>% 
							 	sapply(., function(x) return(x[[1]]))) %>%
					select(c('NAME', 'SIZE', 'FDR q-val'))

#get 3d structure related terms:
#	key words: chromatin, rna, histone, dna
struct_related <- 
	filter(ps_enrichmentmap, grepl('chromatin', ps_enrichmentmap$NAME)) %>%
	rbind(filter(ps_enrichmentmap, grepl('rna', ps_enrichmentmap$NAME))) %>%
	rbind(filter(ps_enrichmentmap, grepl('histone', ps_enrichmentmap$NAME))) %>%
	rbind(filter(ps_enrichmentmap, grepl('dna', ps_enrichmentmap$NAME))) %>%
	distinct(NAME, .keep_all = T) %>%
	arrange(`FDR q-val`) %>%
	filter(!grepl('negative', NAME))  %>%
	filter(!grepl('positive', NAME))  %>%
	filter(!grepl('spliceosome', NAME))  %>%
	filter(!grepl('stress', NAME)) %>%
	filter(!grepl('dna-templated', NAME)) %>%
	filter(!grepl('transesterification', NAME)) %>%
	filter(!grepl('senescence', NAME)) %>%
	filter(!grepl('mrna splic', NAME)) %>%
	filter(!grepl('mrna processing', NAME))%>%
	filter(!grepl('metabolic', NAME))%>%
	filter(!grepl('ercc6', NAME))%>%
	filter(!grepl('promoter', NAME))%>%
	filter(!grepl('capped', NAME))%>%
	filter(!grepl('events', NAME))%>%
	filter(!grepl('internal', NAME))%>%
	filter(!grepl('pri-mirna', NAME))%>%
	filter(!grepl("3'-end", NAME))%>%
	filter(!grepl("rdna", NAME))%>%
	filter(!grepl("stability", NAME))%>%
	filter(!grepl("atp-dependent", NAME))%>%
	filter(!grepl("destabilization", NAME))%>%
	filter(!grepl("stabilization", NAME)) %>%
	filter(!grepl("mrna transcription by", NAME))

struct_related

write_tsv(struct_related,
		'./results/3d_genome/enrichment_map_gene_set/3d_struc__enriched_term.txt')


# 3. generate a clean table of gene and its 3d-related cluster(from step2) ------
# 3.1 get a clean table: one gene and its cluster per line ----
#3d-related clusters:
#	wc -l 3d_struc__enriched_term.txt     #162
#for i in `seq 2 1 162` 
#do
#cluster=$(sed -n ${i}p 3d_struc__enriched_term.txt| cut -f1)
#grep -i "^${cluster}" /Users/xuejia/data/public/gsea/enrichment_map_genesets/human/nov_1_2018/Human_GOBP_AllPathways_no_GO_iea_November_01_2018_symbol.gmt.txt | 
#	awk 'BEGIN{FS="\t"}{for(i=3; i <= NF; i++){print $2"\t"$i}}' >> gene_cluster_clean.txt
#done
#
# 3.2 add the iupred score for each gene --------
clean_cluster <- read_tsv('./results/3d_genome/enrichment_map_gene_set/gene_cluster_clean.txt',
						  col_names = c('cluster', 'gene'))
clean_cluster <- mutate(clean_cluster, cluster = tolower(cluster))

#remove na:
is.na(clean_cluster$gene) %>% table()
clean_cluster <- clean_cluster[!is.na(clean_cluster$gene),]

is_scores <- read_tsv('./results/phase_percentag_all_protein_coding.txt')
is_scores <- mutate(is_scores,mean_score = (percent_of_iupred2a + percent_of_anchor2)/2)


clean_cluster <- left_join(clean_cluster, is_scores, by = c('gene' = 'gene_symbol'))
filter(clean_cluster, gene == 'CTCF')
write_tsv(clean_cluster, './results/3d_genome/enrichment_map_gene_set/gene_cluster_clean_iupred.txt')

clean_cluster_rank <- arrange(clean_cluster, desc(mean_score))
write_tsv(clean_cluster_rank, './results/3d_genome/enrichment_map_gene_set/gene_cluster_clean_iupred_rank.txt')




# 4. calcuate the IDR-rich gene numbers for each 3d-related cluster -------
clean_cluster <- read_tsv('./results/3d_genome/enrichment_map_gene_set/gene_cluster_clean_iupred.txt')
struct_related <- read_tsv('./results/3d_genome/enrichment_map_gene_set/3d_struc__enriched_term.txt')

# 4.1 calculate number of genes with IDR region >=0.3 in each cluster ------
gene_counts_0.3 <- sapply(struct_related$NAME, function(x){
	gene_list <- clean_cluster$gene[clean_cluster$cluster == x]
	gene_score <- filter(is_scores, gene_symbol %in% gene_list)
	return(dim(filter(gene_score, mean_score >= 0.3))[1])
}) 

gene_counts_0.3_dt <- tibble(NAME = names(gene_counts_0.3),
							 gene_counts_0.3 = gene_counts_0.3)

struct_related_counts <- left_join(struct_related, gene_counts_0.3_dt)

#write_tsv(struct_related_counts, 
#	'./results/3d_genome/enrichment_map_gene_set/3d_struc__enriched_term_counts.txt')

# 4.2 barplot of 3d-related gene clusters ------------
library(RColorBrewer)
display.brewer.all()
palette(c(brewer.pal(12, 'Set3'), brewer.pal(8, 'Pastel1')))

struct_related_counts <- read_tsv(
	'./results/3d_genome/enrichment_map_gene_set/3d_struc__enriched_term_counts.txt'
)

struct_related_counts %<>% filter(NAME != 'chromatin assembly') %>%
	filter(NAME != 'chromatin disassembly') %>%
	filter(!grepl('transcription regulatory', NAME))

struct_related_counts$`FDR q-val`[1] <- 10e-6
par(mar = c(4,17,2,1))
#heatmap of top30 gene sets
barplot(struct_related_counts$`FDR q-val`[1:20] %>%
			log10 %>% multiply_by(-1) %>% rev,
		horiz = T,col = 1:20,  xlab = '-log10(FDR q-value)', las =2,
		names.arg = struct_related_counts$NAME[1:20]  %>%  rev)

# 4.3 calculate the AA distribution of all IDR regions of top gene clusters: ------

cluster_gs <- read_tsv('./results/gsea/enrichment_map/cluster_gene_list/cluster_cleanup.txt')
cluster_names <- unique(cluster_gs$enrichmap_cluster)


clean_cluster


cluster_names  <- struct_related_counts$NAME[1:20] 
#subtract the AAs of IDR regions:
IDR_list <- list() #store  clusters

for (i in cluster_names) {
	print(i)
	tmp <- list() #store AAs of each proteins in the cluster
	gene_ls <- filter(clean_cluster,
					  cluster == i & mean_score >= 0.3)$gene
	
	for (J in gene_ls) {
		file_path <- sprintf('./results/iupred2a_results/human/%s.iupred2a.txt', J)
		gene_score <- read.table(file_path)
		#		
		AAs_in_IDR <- gene_score[mean(gene_score$V3 + gene_score$V4) >= 0.5,]$V2
		tmp[[J]] <- AAs_in_IDR
	}
	IDR_list[[i]] <- tmp
}

#check gene numbers for each cluster:
IDR_list$DNA_binding %>% length()
sapply(IDR_list, length)


AA_tibbles <- IDR_list$`regulation of rna splicing` %>% unlist() %>% table() %>% as_tibble()
names(AA_tibbles) <- c('AAs', 'regulation of rna splicing')

for (i in cluster_names[2:20]) {
	print(i)
	counts <- IDR_list[[i]] %>% unlist() %>% table() %>% as_tibble()
	names(counts) <- c('AAs', i)
	AA_tibbles <- left_join(AA_tibbles, counts, by = 'AAs')
}

AA_tibbles

AA_tibbles_ratio <- AA_tibbles[2:21] %>% apply(2, function(x){return(x/sum(x))})
rownames(AA_tibbles_ratio) <- AA_tibbles$AAs


pheatmap::pheatmap(AA_tibbles_ratio %>% t(), 
				   color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

#AA charges:	
#positive: R(Arg)	K(Lys)	H(His)
#negative: D(Asp)	E(Glu)
#polar: S(Ser)	T(Thr)	N(Asn)	Q(Gln)
#special: C(Cys)	G(Gly)	P(Pro)
#Hydrophobic: A(Ala)	V(Val)	I(Ile)	L(Leu)	M(Met)	F(Phe)	Y(Tyr)	W(Trp)


annotation_col = data.frame(
	AAproperty = factor(c('hydrophobic', 'special', 'negative', 'negative', 'hydrophobic',
						  'special', 'positive','hydrophobic', 'positive', 'hydrophobic',
						  'hydrophobic', 'polar', 'special', 'polar', 'positive',
						  'polar', 'polar', 'hydrophobic', 'hydrophobic', 'hydrophobic'))
)

rownames(annotation_col) <- colnames(t(AA_tibbles_ratio))
pheatmap::pheatmap(AA_tibbles_ratio %>% t(), 
				   annotation_col = annotation_col,
				   color = colorRampPalette(c("navy", "white", "firebrick3"))(100))


#normalize with total AA freq:
AA_tibbles.norm <- AA_tibbles[2:21] %>% as.matrix()
rownames(AA_tibbles.norm) <- AA_tibbles$AAs
AA_tibbles.norm <- t(AA_tibbles.norm)

for (i in names(aasumratio)) {
	AA_tibbles.norm[, colnames(AA_tibbles.norm) ==i] = 
		AA_tibbles.norm[, colnames(AA_tibbles.norm) ==i]/aasumratio[names(aasumratio) == i]
}

pheatmap::pheatmap(apply(AA_tibbles.norm, 1,function(x){return(x/sum(x))}) %>% t , 
				   annotation_col = anno.col.df,
				   color = colorRampPalette(c("navy", "white", "firebrick3"))(100))












#randomly selected proteins:
library(Biostrings)
samplenum <- sapply(IDR_list, length) %>% as.numeric()

ensemble_peps <- readAAStringSet('./results/fa/human/all_protein_coding_genes_without_asterix_aa_sequence.fa')

#the total aa frequency among whole proteins:

aasumcount <- alphabetFrequency(ensemble_peps) %>% apply(2, sum) %>% head(20)
aasumratio <- aasumcount /mean(aasumcount)


t1 <- ensemble_peps[sample( 19807, samplenum[1])] %>% alphabetFrequency() %>% colSums() 
t2 <- ensemble_peps[sample( 19807, samplenum[2])] %>% alphabetFrequency() %>% colSums() 

rbind(t1, t2) %>% as_tibble()


aafreq_perm <- ensemble_peps[sample( 19807, samplenum[1])] %>% alphabetFrequency() %>% colSums()

for (i in samplenum[2:length(samplenum)]) {
	tmp <- ensemble_peps[sample( 19807, i)] %>% alphabetFrequency() %>% colSums()
	aafreq_perm <- rbind(aafreq_perm, tmp)
}

aafreq_perm <- aafreq_perm[, 1:20] %>% apply(1, function(x){return(x/sum(x))})
aafreq_perm <- t(aafreq_perm)
apply(aafreq_perm, 1, sum)


anno.col <- tibble(colname = colnames(aafreq_perm))
cr.col <- tibble(colname = rownames(AA_tibbles_ratio), 
				 AAproperty = factor(c('hydrophobic', 'special', 'negative', 'negative', 'hydrophobic',
				 					  'special', 'positive','hydrophobic', 'positive', 'hydrophobic',
				 					  'hydrophobic', 'polar', 'special', 'polar', 'positive',
				 					  'polar', 'polar', 'hydrophobic', 'hydrophobic', 'hydrophobic')))

anno.col <- left_join(anno.col, cr.col)
anno.col.df <- as.data.frame(anno.col[,2])
rownames(anno.col.df) <- anno.col$colname

pheatmap::pheatmap(aafreq_perm , 
				   annotation_col = anno.col.df,
				   color = colorRampPalette(c("navy", "white", "firebrick3"))(100))


# normalize the aa with total frequency:
aasumratio[names(aasumratio) == 'A']
aafreq_perm.norm <- aafreq_perm[, 1:20]

for (i in names(aasumratio)) {
	aafreq_perm.norm[, colnames(aafreq_perm.norm) ==i] = 
		aafreq_perm.norm[, colnames(aafreq_perm.norm) ==i]/aasumratio[names(aasumratio) == i]
}

pheatmap::pheatmap(apply(aafreq_perm.norm, 1,function(x){return(x/sum(x))}) %>% t , 
				   annotation_col = anno.col.df,
				   color = colorRampPalette(c("navy", "white", "firebrick3"))(100))





# 5. motif enrichment analysis --------------------------------------------
clean_cluster <- read_tsv('./results/3d_genome/enrichment_map_gene_set/gene_cluster_clean_iupred.txt')
struct_related <- read_tsv('./results/3d_genome/enrichment_map_gene_set/3d_struc__enriched_term.txt')

struct_related %<>% filter(`FDR q-val` <= 0.05)
tail(struct_related)

clean_gene_filter <- filter(clean_cluster,
							cluster %in% struct_related$NAME)
# 5.1 collect the clean gene list occurred in the top gene clusters ------
#	all genes with q-val <= 0.05
#	subgroups:
#		chromatin/dna/organization
#		RNA/polymerase
#		histone

select(clean_gene_filter, gene) %>% distinct()

# write_tsv(select(clean_gene_filter, gene) %>% distinct(),
# 		  './results/3d_genome/enrichment_map_gene_set/motif_david/all_3d_related_genes.txt')

#chromatin:
chromatin_gene_list <- rbind(
	filter(clean_gene_filter, grepl('chromatin', cluster, ignore.case = T)),
	filter(clean_gene_filter, grepl('dna', cluster, ignore.case = T)),
	filter(clean_gene_filter, grepl('organization', cluster))
)
tail(chromatin_gene_list)
# write_tsv(select(chromatin_gene_list, gene) %>% distinct(),
# 		  './results/3d_genome/enrichment_map_gene_set/motif_david/chromatin_3d_related_genes.txt')

#RNA:
RNA_gene_list <- rbind(
	filter(clean_gene_filter, grepl('RNA', cluster, ignore.case = T)),
	filter(clean_gene_filter, grepl('polymerase', cluster, ignore.case = T)))

tail(RNA_gene_list)
# write_tsv(select(RNA_gene_list, gene) %>% distinct(),
# 		  './results/3d_genome/enrichment_map_gene_set/motif_david/RNA_3d_related_genes.txt')

#histone:
histone_gene_list <- filter(clean_gene_filter, grepl('histone', cluster, ignore.case = T))

tail(histone_gene_list)
# write_tsv(select(histone_gene_list, gene) %>% distinct(),
# 		  './results/3d_genome/enrichment_map_gene_set/motif_david/histone_3d_related_genes.txt')


#overlaps between 3 groups:
table(unique(chromatin_gene_list$gene) %in% unique(RNA_gene_list$gene)) #F705 T 150
table( unique(RNA_gene_list$gene) %in% unique(chromatin_gene_list$gene) ) #F875 T150

#all histone related genes occur in chromatin related gene list:
table(unique(chromatin_gene_list$gene) %in% unique(histone_gene_list$gene)) #F460 T395
table(unique(histone_gene_list$gene) %in% unique(chromatin_gene_list$gene)) #F0 T395

table(unique(histone_gene_list$gene) %in% unique(RNA_gene_list$gene)) #F331 T64
table(unique(RNA_gene_list$gene) %in% unique(histone_gene_list$gene) ) #961 T64


# 5.2 subtract IDR-rich genes (mean score >= 0.2) ---------------

clean_gene_filter_0.2 <- filter(clean_gene_filter, mean_score >= 0.2)
fivenum(clean_gene_filter_0.2$mean_score)

chromatin_gene_list_score_0.2 <- filter(chromatin_gene_list, mean_score >= 0.2)
chromatin_gene_list_score_0.2$mean_score %>% fivenum()

RNA_gene_list_score_0.2 <- filter(RNA_gene_list, mean_score >= 0.2)
fivenum(RNA_gene_list_score_0.2$mean_score)

histone_gene_list_0.2 <- filter(histone_gene_list, mean_score >= 0.2)
fivenum(histone_gene_list_0.2$mean_score)


#write_tsv(select(clean_gene_filter_0.2, gene) %>% distinct(),
# 		  './results/3d_genome/enrichment_map_gene_set/motif_david/all_3d_related_genes_score_0.2.txt')
#
#write_tsv(select(chromatin_gene_list_score_0.2, gene) %>% distinct(),
#		  './results/3d_genome/enrichment_map_gene_set/motif_david/chromatin_3d_related_genes_score_0.2.txt')
#
#write_tsv(select(RNA_gene_list_score_0.2, gene) %>% distinct(),
#		  './results/3d_genome/enrichment_map_gene_set/motif_david/RNA_3d_related_genes_score_0.2.txt')
#
#write_tsv(select(histone_gene_list_0.2, gene) %>% distinct(),
#		  './results/3d_genome/enrichment_map_gene_set/motif_david/histone_3d_related_genes_score_0.2.txt')
#


# 5.3 david enrichment analysis -------------------------------------------
# use david to analysis the protein domains of all/chromatin/rna/histone proteins 
#	and the IDR-riched genes



# 5.3 motif number stats of 3d-related genes ------------------------------
# statatic the motif numbers of proteins grouped as in 5.2.

#Step1: calculate the motif numbers occured in each gene(SMART):
#cut -f1,3,4 smart_clean.txt| sort | uniq -c |sed 's/^[ ]*//; s/ /\t/' > smart_counts_per_gene.txt

counts_smart <- read_tsv('./results/interproscan/clean_data/smart_counts_per_gene.txt',
						 col_names = c('count', 'gene', 'method', 'motif'))

#step2: intersect with gene list 

clean_gene_filter_0.2 %<>% select(c('gene', "mean_score")) %>% filter(!duplicated(gene))

motif_smart_count2 <- left_join(clean_gene_filter_0.2, counts_smart) %>% 
	filter(count >=2) %>%
	group_by(motif) %>%
	count() %>% arrange(desc(n))


par(mar = c(4,18, 2,1))
barplot(motif_smart_count2$n[1:20] %>% rev(), horiz = T, col = 1:20,
		names.arg = motif_smart_count2$motif[1:20] %>% rev(), las =2,
		xlab = 'protein numbers (IDR region >= 20%) involved in chromatin organization',
		main = 'repeated motifs (>=2 in one protein)')



#get gene lists that contain at least 2 same domains:

clean_gene_filter_0.2_smart_count2 <- 
	left_join(clean_gene_filter_0.2, counts_smart) %>% filter(count >=2)

write_tsv(select(clean_gene_filter_0.2_smart_count2, gene) %>% distinct(),
'./results/3d_genome/enrichment_map_gene_set/motif_david/all_3d_related_genes_score_0.2_smart_count2.txt')

rna_motif_genes <- filter(clean_gene_filter_0.2_smart_count2, grepl('WD40', motif)) %>% distinct(gene)


#barplot of all motifs enriched by david--------
motifs_david_iupred2_0.2_all <- read_tsv('./results/3d_genome/enrichment_map_gene_set/motif_david/results/iupred_score_0.2/all_3d_related_genes_iupred_score_0.2_smart.txt')

par(mar = c(4,8,2,1))

barplot(motifs_david_iupred2_0.2_all$PValue[1:20] %>%
		 	log10() %>% multiply_by(-1) %>% sort(decreasing = F), 
		 horiz = T, col = 1:20,
		names.arg = gsub('^.*:', '', x =motifs_david_iupred2_0.2_all$Term[1:20]) %>% rev(), 
		las =2,
		xlab = '-log10(q-value)',
		main = 'all motif enrichment analysis of proteins related to 3d (IDR-region >=20%)')


# calculate the gene numbers containing more than 2 same motifs -----------

#1: calculate the motif numbers occured in each gene(SMART):
#	cut -f1,3,4 smart_clean.txt| sort | uniq -c |sed 's/^[ ]*//; s/ /\t/' > smart_counts_per_gene.txt
#2: calculate the protein numbers which contain repeated motifs:
#	awk '$1 >=2' smart_counts_per_gene.txt| cut -f2 | sort | uniq | wc -l, 4045
#3: calculate the numbers of well known repeated motifs:
#awk '$1 >=2' smart_counts_per_gene.txt | cut -f4 | sort | uniq -c | sort -k 1nr | 
#   sed 's/^[ ]*//; s/ /\t/' > repeated_motif_numbers_smart.txt
# 
#4: calculate the numbers of well known repeated motifs:
#awk '$1 >=2' smart_counts_per_gene.txt | grep 'RNA recognition motif domain' | cut -f2 | sort | uniq -c | wc -l
#awk '$1 >=2' smart_counts_per_gene.txt | grep 'SH3' | cut -f2 | sort | uniq -c | wc -l
#awk '$1 >=2' smart_counts_per_gene.txt | grep 'SH2' | cut -f2 | sort | uniq -c | wc -l
#	RNA recognition motif domain	110
#	SH3 domain	42
#	SH2 domain	13


repeated_motifs <- read_tsv('./results/interproscan/clean_data/repeated_motif_numbers_smart.txt',
                            col_names = c('repeated_numbers', 'motif'))
repeated_motifs

par(mar = c(6, 16, 2,1))

barplot(repeated_motifs$repeated_numbers[1:20] %>% rev(), 
        horiz = T, col = 1:20,
        names.arg = repeated_motifs$motif[1:20]%>% rev(), 
        las =2,
        xlab = 'protein numbers',
        main = 'numbers of proteins containing repeated motifs')


#Barplot (fig 1.B)
selected_motifs <- rbind(repeated_motifs[1:3,],
                         filter(repeated_motifs, 
                    motif %in% c('RNA recognition motif domain', 'SH3 domain', 'SH2 domain')))

par(mar = c(8, 6, 2,1))
library(RSkittleBrewer)
palette(RSkittleBrewer("smarties"))
barplot(selected_motifs$repeated_numbers,
        col = 1:6, las =2, ylim = c(0, 800), 
        cex.axis  =1.2,  cex.lab = 2,
        names.arg = selected_motifs$motif,
        ylab = 'protein numbers')










