### Installing Packages

install.packages("tidyverse")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
install.packages("readr")
install.packages("seqinr")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("decontam")
install.packages("ape")
install.packages("vegan")
install.packages("RColorBrewer")
install.packages("remotes")
remotes::install_github("microbiome/microbiome") `force=TRUE`
if (!requireNamespace("BiocManager", quietly = TRUE))
 install.packages("BiocManager")
 BiocManager::install("DESeq2")

 ### Loading Library 
 
 library(tidyverse)
library ("phyloseq") 
library("readr") 
library("seqinr") 
library("decontam") 
library("ape") 
library("vegan") 
library("RColorBrewer") 
library("microbiome") 
library("DESeq2") 
 
 ### Importing SRA Run Table

 SraRunTable <- read_delim("SraRunTable.txt", delim = ",")
 
 ###Importing Results from QIIME2
 
 count_table <- read_tsv(file="Exports/table/table.tsv", skip = 1)
 count_table <- column_to_rownames(count_table, var = colnames(count_table)[1])
 tree = read_tree("Exports/exported-tree/tree.nwk")
 fasta <- read.fasta(file = "Exports/rep-seqs.fasta/dna-sequences.fasta")
 taxonomy <- read_tsv(file="/Users/leneemason/Desktop/export2/taxonomy.tsv")
 
 ### Checking Sequencing Depth with Rarefraction curves
 
 count_table_df <- as.data.frame(count_table)
 rarecurve(t(count_table_df), step=100, cex=0.5, ylab="ASVs", label=T)

 ### Removing Singletons
 
 count_table_no_singletons <- filter(count_table,rowSums(count_table)>1)

 ### Modifying the Taxonomy Table
 
 head(taxonomy)
 
 taxonomy_mod <-  taxonomy %>%
   mutate(taxonomy=str_replace_all(string=Taxon, pattern="D_\\d*\\__", replacement="")) %>%
   mutate(taxonomy=str_replace_all(string=taxonomy, pattern=";$", replacement="")) %>%
   separate(taxonomy, into=c("Domain", "Phylum", "Class", "Order", "Family", "Genus","Species"), sep=";") %>%
   select (-Taxon, -Confidence) %>%
   column_to_rownames(var = 'Feature ID')
 
head(taxonomy_mod)

### Pull all into Phyloseq Object

ASV =   otu_table(data.frame(count_table_no_singletons), taxa_are_rows =  TRUE)
TAX =   tax_table(as.matrix(taxonomy_mod))
META    =   sample_data(data.frame(Bioproject_mod, row.names = Bioproject_mod$`Library Name`))

head(taxa_names(TAX))
head(taxa_names(ASV))
head(taxa_names(tree))
head(sample_names(ASV))
head(sample_names(META))

ps <- phyloseq(ASV,TAX,META,tree)

rank_names(ps)

unique(tax_table(ps)[, "Domain"])
table(tax_table(ps)[, "Domain"], exclude = NULL)

ps <- subset_taxa(ps, !is.na(Domain) & !Domain %in% c("Unassigned", "Eukaryota"))
table(tax_table(ps)[, "Domain"], exclude = NULL)

table(tax_table(ps)[, "Phylum"], exclude = NULL)

ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c(""))

table(tax_table(ps)[, "Phylum"], exclude = NULL)

pick_new_outgroup <- function(tree.unrooted){
  require("magrittr")
  require("data.table")
  require("ape")
  treeDT <- 
    cbind(
      data.table(tree.unrooted$edge),
      data.table(length = tree.unrooted$edge.length)
    )[1:Ntip(tree.unrooted)] %>% 
    cbind(data.table(id = tree.unrooted$tip.label))
  new.outgroup <- treeDT[which.max(length)]$id
  return(new.outgroup) }
my.tree <- phy_tree(ps)
out.group <- pick_new_outgroup(my.tree)


out.group
new.tree1 <- ape::root(my.tree, outgroup=out.group, resolve.root=TRUE)
new.tree2 <- ape::multi2di(new.tree1)
phy_tree(ps) <- new.tree2
phy_tree(ps)

phylumGlommed = tax_glom(ps, "Phylum")

colourCount = length(table(tax_table(ps)[, "Phylum"], exclude = NULL))
getPalette = colorRampPalette(brewer.pal(9, "Spectral"))
PhylaPalette = getPalette(colourCount)

plot_bar(phylumGlommed, x = "Sample", fill = "Phylum") + 
  scale_fill_manual(values = PhylaPalette)

ps_ra <- microbiome::transform(ps, transform = "compositional")
head(otu_table(ps_ra))

phylumGlommed_RA = tax_glom(ps_ra, "Phylum")

plot_bar(phylumGlommed_RA, x = "Sample", fill = "Phylum") + 
  scale_fill_manual(values = PhylaPalette)

#### Proteobacteria 

ps_proteo_ra <- subset_taxa(ps_ra, Phylum == "Proteobacteria")

colourCount = length(table(tax_table(ps_proteo_ra)[, "Order"], exclude = NULL))
getPalette = colorRampPalette(brewer.pal(9, "Spectral"))
OrderPalette = getPalette(colourCount)

orderGlommed_RA = tax_glom(ps_proteo_ra, "Order")

plot_bar(orderGlommed_RA,, fill = "Order") + 
  scale_fill_manual(values = OrderPalette)

#### Actinobacteria

ps_actino_ra <- subset_taxa(ps_ra, Phylum == "Actinobacteria")

colourCount = length(table(tax_table(ps_actino_ra)[, "Order"], exclude = NULL))
getPalette = colorRampPalette(brewer.pal(9, "Spectral"))
OrderPalette = getPalette(colourCount)

orderGlommed_RA = tax_glom(ps_actino_ra, "Order")

plot_bar(orderGlommed_RA, fill = "Order") + 
  scale_fill_manual(values = OrderPalette)

###The Chlamydiae

ps_chlamy_ra <- subset_taxa(ps_ra, Phylum == "Chlamydiae")

colourCount = length(table(tax_table(ps_chlamy_ra)[, "Family"], exclude = NULL))
getPalette = colorRampPalette(brewer.pal(9, "Spectral"))
FamilyPalette = getPalette(colourCount)

familyGlommed_RA = tax_glom(ps_chlamy_ra, "Family")

plot_bar(familyGlommed_RA, fill = "Family") + 
  scale_fill_manual(values = FamilyPalette)

### Ordinations

ps_hellinger <- microbiome::transform(ps, transform = "hellinger")
head(otu_table(ps_hellinger))

out.pcoa <- ordinate(ps_hellinger, method = "PCoA", distance = "bray")

pcoa_plot = plot_ordination(ps_hellinger, out.pcoa, color ="Temperature", shape = "TissueType") +
  geom_point(size = 3) 
pcoa_plot

evals <- out.pcoa$values$Eigenvalues

pcoa_plot.scaled = plot_ordination(ps_hellinger, out.pcoa, color ="Temperature", shape = "TissueType") +
  geom_point(size = 3) +
  coord_fixed(sqrt(evals[2] / evals[1]))

pcoa_plot.scaled

out.pcoa <- ordinate(ps_hellinger, method = "PCoA", distance = "wunifrac")

wuf_pcoa_plot = plot_ordination(ps_hellinger, out.pcoa, color ="Temperature", shape = "TissueType") +
  geom_point(size = 3) +
  coord_fixed(sqrt(evals[2] / evals[1]))


wuf_pcoa_plot

out.nmds <- ordinate(ps_hellinger, method = "NMDS", distance = "bray")

### Differential abundance with DeSeq2

otu_table(ps)+1
otu_table(ps) <- otu_table(ps)+1
ps_deseq <- phyloseq_to_deseq2(ps, ~Temperature)
ps_deseq <- DESeq(ps_deseq)

### Comparing Different Factors

### Temperature, Time-Zero to 26 degrees

deseq_res_temp_T0_26 <- results(ps_deseq, alpha=0.01, contrast=c("Temperature", "Time-zero", "26"))
summary(deseq_res_temp_T0_26)
sigtab_deseq_res_temp_T0_26 <- deseq_res_temp_T0_26[which(deseq_res_temp_T0_26$padj < 0.01), ]
summary(deseq_res_temp_T0_26)

sigtab_deseq_res_temp_T0_26_with_tax <- cbind(as(sigtab_deseq_res_temp_T0_26, "data.frame"), as(tax_table(ps)[row.names(sigtab_deseq_res_temp_T0_26), ], "matrix"))
sigtab_deseq_res_temp_T0_26_with_tax[order(sigtab_deseq_res_temp_T0_26_with_tax$baseMean, decreasing=T), ]

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
x = tapply(sigtab_deseq_res_temp_T0_26_with_tax$log2FoldChange, sigtab_deseq_res_temp_T0_26_with_tax$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab_deseq_res_temp_T0_26_with_tax$Family = factor(as.character(sigtab_deseq_res_temp_T0_26_with_tax$Family), levels=names(x))
ggplot(sigtab_deseq_res_temp_T0_26_with_tax, aes(x=Family, y=log2FoldChange, color=Class)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))


T0_26_results <-sigtab_deseq_res_temp_T0_26_with_tax[order(sigtab_deseq_res_temp_T0_26_with_tax$log2FoldChange, decreasing=T), ]
write_csv(x = T0_26_results,"T0_26_results.csv")


### Temperature, 26 degrees to 30 degrees

deseq_res_temp_T26_30 <- results(ps_deseq, alpha=0.01, contrast=c("Temperature", "26", "30"))
summary(deseq_res_temp_T26_30)
sigtab_deseq_res_temp_T0_26 <- deseq_res_temp_T26_30[which(deseq_res_temp_T26_30$padj < 0.01), ]
summary(deseq_res_temp_T26_30)

sigtab_deseq_res_temp_T26_30_with_tax <- cbind(as(sigtab_deseq_res_temp_T26_30, "data.frame"), as(tax_table(ps)[row.names(sigtab_deseq_res_temp_T26_30), ], "matrix"))
sigtab_deseq_res_temp_T26_30_with_tax[order(sigtab_deseq_res_temp_T26_30_with_tax$baseMean, decreasing=T), ]

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
x = tapply(sigtab_deseq_res_temp_T26_30_with_tax$log2FoldChange, sigtab_deseq_res_temp_T26_30_with_tax$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab_deseq_res_temp_T26_30_with_tax$Family = factor(as.character(sigtab_deseq_res_temp_T26_30_with_tax$Family), levels=names(x))
ggplot(sigtab_deseq_res_temp_T26_30_with_tax, aes(x=Family, y=log2FoldChange, color=Class)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))


T26_30_results <-sigtab_deseq_res_temp_T26_30_with_tax[order(sigtab_deseq_res_temp_T26_30_with_tax$log2FoldChange, decreasing=T), ]
write_csv(x = T26_30_results,"T26_30_results.csv")

### TissueType, Gut to Feces

ps_deseq <- phyloseq_to_deseq2(ps, ~TissueType)
ps_deseq <- DESeq(ps_deseq)

deseq_res_tissue_Gut_Feces <- results(ps_deseq, alpha=0.01, contrast=c("TissueType", "Gut", "Feces"))
summary(deseq_res_tissue_Gut_Feces)
sigtab_deseq_res_tissue_Gut_Feces <- deseq_res_tissue_Gut_Feces[which(deseq_res_tissue_Gut_Feces$padj < 0.01), ]
summary(deseq_res_tissue_Gut_Feces)

sigtab_deseq_res_tissue_Gut_Feces_with_tax <- cbind(as(sigtab_deseq_res_tissue_Gut_Feces, "data.frame"), as(tax_table(ps)[row.names(sigtab_deseq_res_tissue_Gut_Feces), ], "matrix"))
sigtab_deseq_res_tissue_Gut_Feces[order(sigtab_deseq_res_tissue_Gut_Feces_with_tax$baseMean, decreasing=T), ]

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
x = tapply(sigtab_deseq_res_tissue_Gut_Feces_with_tax$log2FoldChange, sigtab_deseq_res_tissue_Gut_Feces_with_tax$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab_deseq_res_tissue_Gut_Feces_with_tax$Family = factor(as.character(sigtab_deseq_res_tissue_Gut_Feces_with_tax$Family), levels=names(x))
ggplot(sigtab_deseq_res_tissue_Gut_Feces_with_tax, aes(x=Family, y=log2FoldChange, color=Class)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))


Gut_Feces_results <-sigtab_deseq_res_tissue_Gut_Feces_with_tax[order(sigtab_deseq_res_tissue_Gut_Feces_with_tax$log2FoldChange, decreasing=T), ]
write_csv(x = Gut_Feces_results,"Gut_Feces_results.csv")

### TissueType, Feces to Seawater

ps_deseq <- phyloseq_to_deseq2(ps, ~TissueType)
ps_deseq <- DESeq(ps_deseq)

deseq_res_tissue_Feces_SWTR <- results(ps_deseq, alpha=0.01, contrast=c("TissueType", "Feces", "SWTR"))
summary(deseq_res_tissue_Feces_SWTR)
sigtab_deseq_res_tissue_Feces_SWTR <- deseq_res_tissue_Feces_SWTR[which(deseq_res_tissue_Feces_SWTR$padj < 0.01), ]
summary(deseq_res_tissue_Feces_SWTR)

sigtab_deseq_res_tissue_Feces_SWTR_with_tax <- cbind(as(sigtab_deseq_res_tissue_Feces_SWTR, "data.frame"), as(tax_table(ps)[row.names(sigtab_deseq_res_tissue_Feces_SWTR), ], "matrix"))
sigtab_deseq_res_tissue_Feces_SWTR[order(sigtab_deseq_res_tissue_Feces_SWTR_with_tax$baseMean, decreasing=T), ]

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
x = tapply(sigtab_deseq_res_tissue_Feces_SWTR_with_tax$log2FoldChange, sigtab_deseq_res_tissue_Feces_SWTR_with_tax$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab_deseq_res_tissue_Feces_SWTR_with_tax$Family = factor(as.character(sigtab_deseq_res_tissue_Feces_SWTR_with_tax$Family), levels=names(x))
ggplot(sigtab_deseq_res_tissue_Feces_SWTR_with_tax, aes(x=Family, y=log2FoldChange, color=Class)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))


Feces_SWTR_results <-sigtab_deseq_res_tissue_Feces_SWTR_with_tax[order(sigtab_deseq_res_tissue_Feces_SWTR_with_tax$log2FoldChange, decreasing=T), ]
write_csv(x = Feces_SWTR_results,"Feces_SWTR_results.csv")

### Exporting 

descript_file<-writeLines (c ("The phyloseq object in this file are data processed from Brothers et al. 2018 (doi: 10.1098/rspb.2018.0340), NCBI BioProject # PRJNA376395.
The ASV count table was produced using Qiime2 (v. 2020.02), calling DADA2 for denoising, merging, and ASV inference, and a Silva v132 Naive 
Bayes classifier for taxonomy calling.
Samples were collected from Common Sea Urchin (Lytechinus variegatus) from Eagle Harbor, Michigan in May 2015. Authors performed a simulated heat wave in mesocosms to test
impacts on microbiome.
Samples were extracted from the gut, feces, seawater, coelimic fluid, feed, and digested pellet at both temperatures. 
The samples correspond to a time-zero (field) collection and from two different heat stress conditions (26˚C and 30˚C).
The metadata in this file come from the SraRunTable (Bioprojectmod) and have been modified to clearly indicate temperature and tissue type data.
The otu_table object in the phyloseq object contains absolute abundaces of ASVs, without any tranformation or agglomeration at higher taxonomic levels.

qiime2 analysis file (done in Jupyter notebook) =  Analysis Notebook.ipynb" ))

save(ps, descript_file, file = "PRJNA376395_qiime2_processed_data.RData")
