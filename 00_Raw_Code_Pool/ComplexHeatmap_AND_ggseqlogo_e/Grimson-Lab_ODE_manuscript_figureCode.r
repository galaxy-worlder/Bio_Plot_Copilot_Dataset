### This is the R code needed to generate all necessary figures for the manuscript "Pervasive cis-regulatory co-option of a transposable 
### element family reinforces cell identity across the mouse immune system" by Jason Chobirko. Enjoy!
#
### This code requires that you've downloaded and analyzed all of the data for each of the different datasets. It requires various 
### intermediate files, such as a 450+ million row object reporting correlations between ATAC-seq and RNA-seq across 80+ datasets! 
#
### Otherwise, I hope this code is well commented and provides the necessary insights for your own analyses. Good luck! 

### NOTE: All the below R code was ran with R version 4.1.3 on a server operating on Rocky Linux 9.0, with any particular package requirement
### being loaded in at the beginning of each figure section (or via adding lbzip2 to the path beforehand for large objects) 

################
### Figure 1 ###
################

# Load all the relevant libraries you need for this section
suppressPackageStartupMessages({ 
	library(DESeq2); library(Biostrings); library(dplyr); library(ggtree); library(RColorBrewer)
	library(data.table); setDTthreads(8); library(ggplot2); library(ggbeeswarm); library(cowplot)
})


## A 

### NOTE: This figure was created entirely in illustrator using the total samples for each cell type in ATAC-seq as defined in 
### the paper Yoshida H et al 2019 (https://pmc.ncbi.nlm.nih.gov/articles/PMC6785993/)


## B & C

# Load in the list of TE-derived peak ids 
tes <- unique(read.table("atac/peak.te.int")$V4)

# Load in all ImmGen Peaks via ATAC-seq to calculate CPM
tmp <- read.table("atac/immgen.iter")$V1
f.m <- as.data.frame(fread("atac/summit_unif_peaks_counts.txt")); colnames(f.m) <- c("chr", "start", "end", "id", tmp )
rownames(f.m) <- with(f.m, paste0(chr, ":", start, "-", end, "|", id)); l <- f.m$end - f.m$start

# Load in the lut so you know all the different samples to iterate over
lut.pca <- read.table("atac/immgen_heatmap_lut.txt", sep = "\t", header = T, comment.char = "")

# Get each of the separate objects that contain either the (likely) TE-derived peaks or non-TE-derived peaks
f.te <- f.m[ f.m$id %in% tes, rownames(lut.pca) ]; f.nte <- f.m[ !f.m$id %in% tes, rownames(lut.pca) ]

# Now to generate the PCA plots!
coldata <- data.frame("tmp" = rownames(lut.pca), "condition" = gsub("_[1-9]", "", tmp), "type" = gsub(".*_", "rep", tmp), row.names = 1)
ddste <- DESeqDataSetFromMatrix(countData = f.te, colData = lut.pca, design = ~ name)
ddsnte <- DESeqDataSetFromMatrix(countData = f.nte, colData = lut.pca, design = ~ name)

val <- unique(lut.pca$col); names(val) <- unique(lut.pca$name)
val <- val[c("Stem", "B", "Tgd", "Tab", "Act T", "ILC", "Mo", "MF/GN", "DC", "Stroma")]

d <- plotPCA(vst(ddste), intgroup = c("name"), ntop = 1000, returnData = T); pV <- round(100 * attr(d, "percentVar"), 1)
d$cell <- factor(lut.pca[ d$name.1, "name"], levels = names(val) )
d2 <- plotPCA(vst(ddsnte), intgroup = c("name"), ntop = 1000, returnData = T); pV2 <- round(100 * attr(d2, "percentVar"), 1)
d2$cell <- factor(lut.pca[ d2$name.1, "name"], levels = names(val) )

g <- list()
g[[1]] <- ggplot(d, aes(PC1, PC2, color = cell, alpha = 0.88)) + geom_point(size = 3.03, shape = 16) + xlab(paste0("PC1: ", pV[1], "% variance")) + 
	ylab(paste0("PC2: ", pV[2], "% variance")) + ggtitle("TE peaks") + theme_classic() + scale_color_manual(values = val, name = "Cell type") +
	theme(legend.key.height = unit(0.14, "in"), plot.title = element_text(size = 8, color = "black", hjust = 0.5), axis.title = element_text(size = 8, color = "black"), axis.text = element_text(size = 8, color = "black"), legend.text = element_text(size = 8, color = "black"), legend.title = element_text(size = 8, color = "black"), axis.ticks.length = unit(.04, "in"), axis.ticks = element_line(linewidth = 1/.75/.pt, color = "black"), axis.line = element_line(linewidth = 1/.75/.pt, color = "black")) +
	guides(color = guide_legend(override.aes = list(alpha = 1, shape = 16, size = 3.03)), alpha = "none")

### NOTE: It came out "upside-down" relative to the other plot, so make PC2 negative to compensate 
g[[2]] <- ggplot(d2, aes(PC1, -PC2, color = cell, alpha = 0.88)) + geom_point(size = 3.03, shape = 16, show.legend = F) + xlab(paste0("PC1: ", pV2[1], "% variance")) + 
	ylab(paste0("PC2: ", pV2[2], "% variance")) + ggtitle("nonTE peaks") + theme_classic() + scale_color_manual(values = val) +
	theme(plot.title = element_text(size = 8, color = "black", hjust = 0.5), axis.title = element_text(size = 8, color = "black"), axis.text = element_text(size = 8, color = "black"), axis.ticks.length = unit(.04, "in"), axis.ticks = element_line(linewidth = 1/.75/.pt, color = "black"), axis.line = element_line(linewidth = 1/.75/.pt, color = "black"))

# Create the first row of panels! Leave space for the illustrator figure
g1 <- plot_grid( NULL, plot_grid( plotlist = g, labels = c("B", "C"), label_size = 12, rel_widths = c(1,0.7) ), labels = c("A", ""), rel_widths = c(1,3.2), label_size = 12 )

# Remove all the objects you no longer need
rm( coldata, d, d2, g, ddsnte, ddste, f.m, f.nte, f.te, l, lut.pca, pV, pV2, tes, tmp, val )


## D

# Start by reading in the summarized shuffle enrichment results for ImmGen datasets
res <- read.table("atac/TEshuffle_results.txt", header = T)

# Overwrite the column for the sample name
res$samp <- gsub("_onlyEnhancer", "", res$samp)

# Overwrite the fisher's test p-values with their "bonferroni-adjusted" values and filter based on those
res$fep <- p.adjust(res$fisher_enriched_padj, method = "bonferroni")
p.bonf <- res %>% filter(perm_enriched_pval <= 0.05 & fep <= 0.05 & repClass %in% c("DNA", "LINE", "SINE", "LTR"))

# Read in the look-up-table that has all the colors and file names
lut.ig <- read.table("atac/immgen_heatmap_lut.txt", header = T, sep = "\t", comment.char = "")

# Read in all of the consensus sequences for Dfam v3.6 (taken from RepeatMasker when loading in Dfam annotations from https://www.dfam.org/releases/Dfam_3.6/families/Dfam_curatedonly.h5.gz)
cons <- readDNAStringSet("rm/dfam3.6_cons.fa")

# Iterate through the families that are significant, and find the right consensus sequence to use for an alignment to order the x-axis
fams <- unique(p.bonf$repName)
for (f in 1:length(fams)) {
	if ( f == 1 ) {
		tmp <- cons[which(grepl(paste0("^", fams[f], "#"), names(cons)))]
	} else {
		# If no hits (ie it's a LINE with consensus sequence chunks for orf2/5'/3') then check first for orf2, then 5', then 3'
		if ( length( which(grepl(paste0("^", fams[f], "#"), names(cons))) ) == 0 ) {
			if ( length ( which(grepl(paste0("^", fams[f], "_orf2#"), names(cons))) ) == 1 ) {
				tmp <- c(tmp, cons[which(grepl(paste0("^", fams[f], "_orf2#"), names(cons)))])
			} else if ( length ( which(grepl(paste0("^", fams[f], "_5end#"), names(cons))) ) == 1 ) {
				tmp <- c(tmp, cons[which(grepl(paste0("^", fams[f], "_5end#"), names(cons)))])
			} else {
				tmp <- c(tmp, cons[which(grepl(paste0("^", fams[f], "_3end#"), names(cons)))])
			}
		} else {
			tmp <- c(tmp, cons[which(grepl(paste0("^", fams[f], "#"), names(cons)))])
		}
	}
}

names(tmp) <- gsub("_3end", "", gsub("_5end", "", gsub("_orf2", "", names(tmp))))

# Output the consensus sequences that match your current list of elements
writeXStringSet(tmp, "tmp.fa", append=FALSE, compress=FALSE, format="fasta")

# Generate sequence alignments between all the above consensus sequences to get the order for the current figure
system("export PATH=/programs/mafft/bin:$PATH; export OMP_NUM_THREADS=16")
system("mafft --ep 0.123 --quiet --thread 12 tmp.fa > tmp.afa")
system("clipkit tmp.afa -o tmp.clipkit -l")
system("/programs/FastTree-2.1.11/FastTreeMP -nt -gtr -gamma -seed 8 -log tmp.clipkit.ftmp.log < tmp.clipkit > tmp.clipkit.ftmp.nwk")
system("rm tmp.fa tmp.afa tmp.clipkit")

# Now load in the tree output from above so that you can use it to order your dotplot. 
ggtree_plot <- ggtree(read.tree("tmp.clipkit.ftmp.nwk"), branch.length = "none") + geom_tiplab()

# Reorder the dotplot by the clusters generated above
p.bonf$repName <- factor(p.bonf$repName, levels = gsub("#.*", "", get_taxa_name(ggtree_plot)))
lab <- get_taxa_name(ggtree_plot); classlab <- gsub("-consensus", "", gsub("/.*", "", gsub(".*#", "", lab)))
lut <- brewer.pal(4, "Set1"); names(lut) <- c("DNA", "LINE", "SINE", "LTR")

### NOTE: You'll need to manually add in break blocks so that you can get the requested plot visualization. What a pain :'(
val <- c( rownames(subset(lut.ig, name == "Stem")), "b1", "b2", rownames(subset(lut.ig, name == "B")), "b3", "b4", rownames(subset(lut.ig, name == "Tgd")), "b5", "b6", rownames(subset(lut.ig, name == "Tab")), "b7", "b8", rownames(subset(lut.ig, name == "Act T")), "b9", "b10", rownames(subset(lut.ig, name == "ILC")), "b11", "b12", rownames(subset(lut.ig, name == "Mo")), "b13", "b14", rownames(subset(lut.ig, name == "MF/GN")), "b15", "b16", rownames(subset(lut.ig, name == "DC")), "b17", "b18", rownames(subset(lut.ig, name == "Stroma")) )
p.bonf$lab <- factor(p.bonf$samp, levels = val); val <- unique(lut.ig[,c("name","col")]); tmp <- val[,"col"]; names(tmp) <- val[,"name"]; val <- tmp; rm(tmp)
p.bonf$fill <- lut.ig[p.bonf$samp, "name"]

### NOTE: To manually create a label for the x-axis coloring, you need to create some fake data to plot below and use the fill parameter to 
### add in the legend. Yeah
f1 <- data.frame(x = "ORR1E", y = "proB.CLP.BM", f = factor( unique(p.bonf$repClass), levels = c("DNA", "LINE", "SINE", "LTR", "Retroposon"))) 

# You created an object that contains the coordinates with which to generate rectangular boxes (around the whole chunk?) in the below plot. 
# Yay
p2 <- data.frame()
for (h in names(val)) {
	a <- which(levels(p.bonf$lab) %in% rownames(subset(lut.ig, name == h)))
	p2 <- rbind(p2, data.frame(ymin = min(a)-0.5, ymax = max(a)+0.5, xmin = 0, xmax = length(levels(p.bonf$repName))+1, fill = h))
}

# Now make that plot! 
g2 <- ggplot( arrange(p.bonf, desc(fep)), aes(x=lab, y=repName, color = -log10(fep), size = bound_uniq), alpha = 0.75) + geom_point(data = f1, aes(x = y, y = x, fill = f), alpha = 0, size = 5, shape = 22, stroke = 1/.75/.pt, inherit.aes = F) + scale_fill_manual(values = lut) + 
	geom_segment(data = p2, aes(y = length(unique(p.bonf$repName))+1, x = ymin - 0.5, xend = ymax + 0.5), color = val[p2$fill], inherit.aes = F, show.legend = F, linewidth = 2/.75/.pt) + geom_point(aes(size = ifelse( bound_uniq > 625, 625, ifelse( bound_uniq < 25, 25, bound_uniq )))) + scale_size(breaks = c(25, 325, 625), limits = c(25,625), range = c(1,3), labels = c("25", "325", "625")) + cowplot::theme_minimal_grid(line_size = 0.5/.75/.pt, color = "gray69") +
        ylab("") + xlab("") + scale_color_gradient(low = "gray69", high = "dark orange", limits = c(-log10(0.05),-log10(1e-45)), oob = scales::squish, name = "-log10\nenriched\np-value", breaks = c(-log10(0.05), 15, 30, 45), labels = c("1.3", "15", "30", "45+")) +
        scale_x_discrete(drop = F, breaks = rownames(lut.ig)) + labs(size = "Accessible\nelements") + coord_cartesian(ylim = c(1,length(levels(p.bonf$repName))), clip = "off") + guides(color = guide_colorbar(order = 1, barheight = unit(0.5, "in"), barwidth = unit(0.1, "in")), size = guide_legend(order = 2), fill = guide_legend(title = "TE class", override.aes = list(alpha = 1, size = 3), order = 3), alpha = "none") +
	theme(legend.text = element_text(size = 8, color = "black"), legend.key.height = unit(0.15, "in"), legend.title = element_text(size = 8, color = "black"), plot.margin = margin(r = 5, t = 35, unit = "pt"), axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(color = lut[classlab], size = 8)) +
	geom_text(data = p2, aes(y = length(unique(p.bonf$repName))+2, x = ((ymax + ymin) / 2), label = gsub("embryonic_facial_prominence", "efp", fill)), inherit.aes = F, show.legend = F, size = 8/.pt, color = "black")

# Remove all the objects you no longer need
rm( a, classlab, cons, f, f1, fams, ggtree_plot, h, lab, lut, lut.ig, p2, res, val )
system("rm tmp.clipkit.ftmp.nwk")


## E

a <- as.data.frame(data.table::fread("dnase/TEshuffle_results.txt"))
a$tissue <- factor(gsub("_[1-9+].*", "", a$samp), levels = unique(gsub("_[1-9+].*", "", a$samp))); a$lab <- gsub("_", "\n", a$tissue)

# Get a sense of how many total samples there are for each tissue type
totis <- as.data.frame(table(unique(a[,c("samp", "lab")])$lab), row.names = 1)

# Remove all of the tissues where you only have a single replicate
a <- a[ which(!gsub("_[1-9+].*", "", a$lab) %in% rownames(totis[totis$Freq == 1, , drop = F]) ), ]

# Rename some of the tissues for ease of visualization
a$lab <- gsub("gastrocnemius", "gastr.", a$lab); a$lab <- gsub("adrenal\ngland", "adr.", a$lab)
a$lab <- gsub("embryonic\nfacial\nprominence", "efp", a$lab); a$lab <- gsub("urinary\nbladder", "bladder", a$lab)
a$lab <- gsub("\ncerebral\ncortex", " cc", a$lab); a$lab <- gsub("cerebellum", "cereb.", a$lab)
a$lab <- gsub("layer\nof\nhippocampus", "hipp.", a$lab); a$lab <- gsub("frontal\ncortex", "fr. cor.", a$lab)
a$lab <- gsub("\n", " ", a$lab); a$lab <- factor(a$lab, levels = unique(sort(a$lab)))

# Instead of it being alphabetical, make it alphabetical within each kind of tissue
b <- unique(a$lab); a$lab <- factor(a$lab, levels = c("cereb.", "forebrain", "fr. cor.", "gastr.", "hindbrain", "hipp.", "left cc", "midbrain", "neural tube", "right cc", "b1", "adr.", "bladder", "embryo", "efp", "heart", "kidney", "limb", "ovary", "testis", "b2", "liver", "lung", "spleen", "thymus"))

# Make some dotplots showing the distribution of significant instances where ORR1E/ORR1D2 are enriched in DNase-seq datasets
a$sig <- factor(ifelse(a$perm_enriched_pval <= 0.05, "Enriched", ifelse(a$perm_depleted_pval <= 0.05, "Depleted", "Not sig")), levels = c("Enriched", "Not sig", "Depleted"))

# You created an object that contains the coordinates with which to generate rectangular boxes in the below plot. Yay!
tmp <- data.frame("name" = c("Brain", "Non-Brain", "Immune")); tmp$samp <- list( c("cereb.", "forebrain", "fr. cor.", "gastr.", "hindbrain", "hipp.", "left cc", "midbrain", "neural tube", "right cc"), c("adr.", "bladder", "embryo", "efp", "heart", "kidney", "limb", "ovary", "testis"), c("liver", "lung", "spleen", "thymus"))
p2 <- data.frame(); tmp$x <- 0
for ( f in 1:nrow(tmp) ) {
	d <- which(levels(a$lab) %in% tmp[f,"samp"][[1]])
	p2 <- rbind(p2, data.frame(x = min(d), xend = max(d), y = -4.5))
	tmp[f,"x"] <- mean(d)
}

# Generate the figure
p <- subset(a, repName == "ORR1E") 
g3 <- ggplot(p, aes(x = lab, y = log2(bound_uniq/exp_uniq))) + geom_hline(yintercept = 0, lty = "31", color = "gray69", linewidth = 1/.75/.pt) + geom_beeswarm(data = p, aes(x = lab, y = log2(bound_uniq/exp_uniq), color = sig), size = 1, shape = 16, inherit.aes = F, cex = 0.5) +
	stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="crossbar", width=0.5, color = "gray69", fatten = 2, linewidth = 1/.75/.pt) + scale_color_manual(values = c("Enriched" = "#ff7518", "Not sig" = "black"), name = "") + 
	theme_classic() + scale_x_discrete(drop = F, breaks = b) + xlab("") + ylab(expr("log"[2]*"(obs/exp)")) + guides(color = guide_legend(override.aes=list(size = 2))) +
	coord_cartesian(ylim = c(-1.1,2.5), clip = "off") + geom_segment(data = p2, aes(x = x, xend = xend, y = y), color = "black", inherit.aes = F, show.legend = F, linewidth = 1/.75/.pt) +
	geom_text(data = tmp, aes(x = x, y = -4.9, label = name), size = 8/.pt, color = "black", inherit.aes = F, show.legend = F) +
	theme(plot.margin=unit(c(5.5,5.5,40.5,20.5), 'pt'), text = element_text(size = 8, color = "black"), axis.text.x = element_text(size = 8, color = "black", angle = 45, vjust = 1, hjust = 1), axis.title.y = element_text(size = 8, color = "black"), axis.text.y = element_text(size = 8, color = "black"), legend.text = element_text(size = 8, color = "black"), axis.ticks.length = unit(.04, "in"), axis.ticks = element_line(linewidth = 1/.75/.pt, color = "black"), axis.line = element_line(linewidth = 1/.75/.pt, color = "black"))

# Save the figures! 
pdf("figures/fig1be.pdf", height = 6.8, width = 6.9)
plot_grid( g1, g2, g3, ncol = 1, labels = c("", "D", "E"), rel_heights = c(1, 1.4, 1), label_size = 12 )
dev.off()

# Remove all the objects you no longer need
rm( a, b, d, f, p, p2, tmp, totis, g1, g2, g3, p.bonf )



################
### Figure 2 ###
################

# Add lbzip2 to your path and make sure that the R version is correct
# export PATH=/home/jdc397/lbzip2-2.5/bin:$PATH; module load R/4.1.3-r9; R --vanilla

# Load all the relevant libraries you need for this section
suppressPackageStartupMessages({ 
	library(data.table); setDTthreads(8); library(ggpubr); library(fastSave); library(cluster); library(RColorBrewer)
	library(ggplot2); library(cowplot); library(dplyr); library(circlize); library(ComplexHeatmap); library(ggthemes)
	library(scales); library(tidyr); library(doParallel); library(stringr);library(ggseqlogo)
})

# First, load in all ImmGen Peaks via ATAC-seq to calculate CPM
tmp <- read.table("atac/immgen.iter")$V1
f.m <- as.data.frame(fread("atac/summit_unif_peaks_counts.txt")); colnames(f.m) <- c("chr", "start", "end", "id", tmp )
rownames(f.m) <- with(f.m, paste0(chr, ":", start, "-", end, "|", id)); l <- f.m$end - f.m$start

# Load in the lut so you know all the different samples to iterate over
lut.a <- read.table("atac/immgen_heatmap_lut.txt", sep = "\t", header = T, comment.char = "")

# Generate log2(CPM+1) for *all* the peaks
t.c <- as.data.frame(matrix(0,nrow=nrow(f.m),ncol=1), row.names = rownames(f.m) )
l2.atac.cpm <- log2( bind_cols(lapply( rownames(lut.a), function(h) {
	s <- f.m[ , h, drop = F ]
	t.c[ , 1 ] <- 1e6 * s[ , 1, drop=F ] / l / sum( s[ , 1, drop=F ] / l )  
	colnames(t.c) <- h; return(t.c)
}) ) + 1 )

# Save the above object so you don't have to regenerate it for further analyses
write.table(l2.atac.cpm, file = "atac/rObject_l2.atac.cpm.txt", sep = "\t", quote = F)

# Run "qnorm" to generate the quantile normalized values on the above object (https://pypi.org/project/qnorm/)
system("qnorm atac/rObject_l2.atac.cpm.txt > atac/rObject_l2.atac.cpm.qnorm.txt")
system("rm atac/rObject_l2.atac.cpm.txt")

# Load in the qnorm output
l2.atac.cpm.qnorm <- data.frame(fread("atac/rObject_l2.atac.cpm.qnorm.txt", header = T), row.names = 1)

# Save the above object as a compressed file so that it saves space and you don't have to regenerate it! Woot
save.lbzip2(l2.atac.cpm.qnorm, file = "atac/rObject_l2.atac.cpm.qnorm.RDataFS", n.cores = 16)

# Load in the list providing the peaks which are overlapping ORR1E/ORR1D2 with at least one summit across one dataset
te.int <- subset(as.data.frame(fread("atac/peak.te.int")), V8 %in% c("ORR1E", "ORR1D2")); te.int$id <- with(te.int, paste0(V5, ":", V6, "-", V7, "|", V8))
te.int$pid <- with(te.int, paste0(V1, ":", V2, "-", V3, "|", V4))

# Subset only the peaks that are assigned to ORR1E
r1e.atac <- l2.atac.cpm.qnorm[unique(subset(te.int, V8 == "ORR1E")$pid), rownames(lut.a)]

# Generate the color scale that recapitulates Pheatmap's default since it's cool
col.fun <- colorRamp2(seq(-3, 3, length.out = 7), rev(brewer.pal(n = 7, name = "RdYlBu")))

# Get the cell type order and colors ready for the heatmap ordering
col.a <- unique(lut.a$col); names(col.a) <- unique(lut.a$name)

# Time for the same figure but using the CPM of the peaks these ORR1E elements underlie! 
set.seed(8)
tmp <- draw( Heatmap(lut.a$name, row_title_rot = 0, col = structure(col.a, names = names(col.a)), name = "Cell Type", show_row_names = F, width = unit(3, "mm")) +
        Heatmap( scale(t(r1e.atac[ , rownames(lut.a) ])), column_km = 8, use_raster = T, column_km_repeats = 200, show_column_dend = F, cluster_column_slices = FALSE, cluster_rows = F, show_row_names = F, column_title = "ORR1E", show_column_names = F, name = "CS.CPM", col = col.fun, border = T) )
r1e.col.cpm <- column_order(tmp)$CS.CPM

tmp <- as.data.frame(matrix("", ncol = 1, nrow = nrow(r1e.atac) ))
for (i in 1:length(r1e.col.cpm)) {
	tmp[r1e.col.cpm[[i]],"V1"] <- i 
}

# NEW ORDER: 1 - Stem & Prog, 2 - B, 3 - T, 4 - NK, 5 - ILC3, 6 - MF/GN & Mo, 7 - DC, 8 - Random
val <- data.frame("i" = 1:8, "o" = c(1,7,2,6,5,8,3,4), row.names = 1)
r1e.col.cpm <- factor(val[tmp$V1, "o"], levels = c(1:8))

# To figure out which column you'll add the annotation to, you need to calculate the floored average for each block
a <- unlist(lapply(unique(lut.a$name), function(i){ floor(mean(which(lut.a$name == i))) }))
b <- rep("", nrow(lut.a)); b[a] <- unique(lut.a$name)

g1 <- grid.grabExpr( draw( Heatmap(t(scale(t(r1e.atac[ , rownames(lut.a) ]))), use_raster = T, 
	top_annotation = HeatmapAnnotation("Name" = anno_text(b, rot = 45, just = "left", location = unit(1, "mm"), gp = gpar(fontsize = 8, color = "black")), "Cell Type" = lut.a$name, col = list("Cell Type" = col.a), annotation_name_side = "left", annotation_name_gp = gpar(col = "black", fontsize = 8), show_legend = c(F,F), simple_anno_size = unit(3, "mm")), 
	heatmap_legend_param = list(direction = "horizontal", title_position = "topcenter", legend_width = unit(1, "in"), legend_height = unit(0.1, "in"), title = "log2(CPM+1)", labels_gp = gpar(col = "black", fontsize = 8), title_gp = gpar(col = "black", fontsize = 8)), 
	row_split = r1e.col.cpm, row_title = c("Stem\n& Prog", "B", "T", "NK", "ILC3", "MF", "DC", "Bkgd"), row_title_rot = 0, show_row_dend = F, 
	cluster_row_slices = FALSE, cluster_columns = F, show_row_names = F, show_column_names = F, row_title_gp = gpar(col = "black", fontsize = 8),
	col = col.fun, border = T, row_gap = unit(.04, "in")),
	heatmap_legend_side = "bottom", padding = unit(c(0.5,0.5,0.5,1), "cm") ) )


### NOTE: While Figure 2A is done, do the same analysis for ORR1D2 to get a complete list of both clusters for downstream analyses here

# Load the subset for ORR1D2
r1d2.atac <- l2.atac.cpm.qnorm[unique(subset(te.int, V8 == "ORR1D2")$pid), rownames(lut.a)]

set.seed(8)
tmp <- draw( Heatmap(lut.a$name, row_title_rot = 0, col = structure(col.a, names = names(col.a)), name = "Cell Type", show_row_names = F, width = unit(3, "mm")) +
        Heatmap(scale(t(r1d2.atac[ , rownames(lut.a) ])), column_km = 6, use_raster = T, column_km_repeats = 200, show_column_dend = F, cluster_column_slices = FALSE, cluster_rows = F, show_row_names = F, column_title = "ORR1D2", show_column_names = F, name = "CS.CPM", col = col.fun, border = T) )
r1d2.col.cpm <- column_order(tmp)$CS.CPM

tmp <- as.data.frame(matrix("", ncol = 1, nrow = nrow(r1d2.atac) ))
for (i in 1:length(r1d2.col.cpm)) {
	tmp[r1d2.col.cpm[[i]],"V1"] <- i 
}

# NEW ORDER: 1 - Stem & Prog, 2 - B, 3 - T, 4 - NK, 5 - MF/DC, 6 - Random
val <- data.frame("i" = 1:6, "o" = c(2,6,5,4,1,3), row.names = 1)
r1d2.col.cpm <- factor(val[tmp$V1, "o"], levels = c(1:6))

# Now that you have both ORR1E and ORR1D2 clusters, save it so you can use it in downstream analyses. Yeah! 
names(r1e.col.cpm) <- rownames(r1e.atac); names(r1d2.col.cpm) <- rownames(r1d2.atac)
ode <- te.int[ , c(5:8,13:14) ]; colnames(ode)[1:4] <- c("chr", "start", "end", "repName"); rownames(ode) <- ode$id
ode$cpm.c <- ifelse( ode$pid %in% names(r1e.col.cpm) & ode$repName == "ORR1E", r1e.col.cpm[ode$pid], r1d2.col.cpm[ode$pid] )
write.table(ode, "atac/rObject_ode.txt", quote = F, sep = "\t", row.names = F)

# Remove all the objects you no longer need
rm(a, b, col.a, col.fun, f.m, i, l, l2.atac.cpm.qnorm, lut.a, r1d2.atac, r1d2.col.cpm, r1e.atac, r1e.col.cpm, t.c, te.int, tmp, val, ode)



## B

# Load the qnorm output and annotations for ode elements
load.lbzip2("atac/rObject_l2.atac.cpm.qnorm.RDataFS", n.cores = 16)
ode <- read.table("atac/rObject_ode.txt", header = T)

# Get a list of all the samples that each cluster are most associated with so you have this information handy. Also save it so you can just
# load it in later. You wanted to have more than a majority of elements have this threshold, so go with an even 50%. Yeah
samp.names.1e <- samp.names.1d2 <- list(); a.thresh <- 1; b.thresh <- 0.50
for (f in 1:7 ) {
	s <- unique(subset(ode, repName == "ORR1E" & cpm.c == f)$pid); tmp <- colSums( t(scale(t(l2.atac.cpm.qnorm[ s, ]))) > a.thresh ) / length(s)
	samp.names.1e[[f]] <- c( names(tmp)[which(tmp >= b.thresh)] )

	if ( f < 6 ) {
		s <- unique(subset(ode, repName == "ORR1D2" & cpm.c == f)$pid); tmp <- colSums( t(scale(t(l2.atac.cpm.qnorm[ s, ]))) > a.thresh ) / length(s)
		samp.names.1d2[[f]] <- c( names(tmp)[which(tmp >= b.thresh)] )
	} else if ( f == 6 ) {
		samp.names.1d2[[f]] <- c( colnames(l2.atac.cpm.qnorm) )
	}
}
samp.names.1e[[8]] <- colnames(l2.atac.cpm.qnorm); rm(a.thresh, b.thresh, f, s, tmp)

# Save the objects so that you don't need to redo this again
save.lbzip2(samp.names.1e, file = "atac/rObject_samp.names.1e.RDataFS", n.cores = 16)
save.lbzip2(samp.names.1d2, file = "atac/rObject_samp.names.1d2.RDataFS", n.cores = 16)

# Load in the homer motif enrichment analyses
df.mot.enrich <- data.frame()
for (f in 1:7) { 
	if ( f < 6 ) {
		tmp <- read.table(paste0("homer/mouse/orr1d2_km_", f, ".vs.6.txt"), sep = "\t", skip = 1); tmp$clust <- f; tmp$samp1 <- "orr1d2"; tmp$samp2 <- "bkgd" 
		df.mot.enrich <- rbind(df.mot.enrich, tmp)
		tmp <- read.table(paste0("homer/mouse/orr1d2_km_6.vs.", f, ".txt"), sep = "\t", skip = 1); tmp$clust <- f; tmp$samp1 <- "bkgd"; tmp$samp2 <- "orr1d2"
		df.mot.enrich <- rbind(df.mot.enrich, tmp) 		
	}

	tmp <- read.table(paste0("homer/mouse/orr1e_km_", f, ".vs.8.txt"), sep = "\t", skip = 1); tmp$clust <- f; tmp$samp1 <- "orr1e"; tmp$samp2 <- "bkgd" 
	df.mot.enrich <- rbind(df.mot.enrich, tmp)
	tmp <- read.table(paste0("homer/mouse/orr1e_km_8.vs.", f, ".txt"), sep = "\t", skip = 1); tmp$clust <- f; tmp$samp1 <- "bkgd"; tmp$samp2 <- "orr1e"
	df.mot.enrich <- rbind(df.mot.enrich, tmp) 		
}

# With the individual comparisons loaded in, wrangle the object for plotting purposes
df.mot.enrich <- unique(df.mot.enrich); rm(tmp); df.mot.enrich$V7 <- as.numeric(gsub("%", "", df.mot.enrich$V7)); df.mot.enrich$V9 <- as.numeric(gsub("%", "", df.mot.enrich$V9))
df.mot.enrich <- separate(df.mot.enrich, col = "V1", into = c("tmp", "data", "origin"), sep = "/")
df.mot.enrich <- separate(df.mot.enrich, col = "tmp", into = c("motif", "motif.fam"), sep = "\\(")
df.mot.enrich$motif.fam.short <- gsub("\\?", "", gsub("\\?,", "", gsub(").*", "", df.mot.enrich$motif.fam)))
df.mot.enrich$motif.fam.short <- ifelse( is.na(df.mot.enrich$motif.fam.short) | df.mot.enrich$motif.fam.short == "", "Other", df.mot.enrich$motif.fam.short )
df.mot.enrich[df.mot.enrich$motif.fam.short == "forkhead","motif.fam.short"] <- "Forkhead"
df.mot.enrich[df.mot.enrich$motif.fam.short == "ETS:IRF","motif.fam.short"] <- "ETS"
df.mot.enrich$motif.fam <- paste0("(", df.mot.enrich$motif.fam); df.mot.enrich$motif.uniq <- paste0(df.mot.enrich$motif, df.mot.enrich$motif.fam, "/", df.mot.enrich$data)

# Adjust the p-values separately for ORR1D2 and ORR1E!
df.mot.enrich$padj <- 0
a <- which(grepl("1e", df.mot.enrich$samp1) | grepl("1e", df.mot.enrich$samp2))
df.mot.enrich[ a, "padj" ] <- p.adjust(exp(df.mot.enrich[ a, "V4" ]), method = "BH")
a <- which(grepl("1d2", df.mot.enrich$samp1) | grepl("1d2", df.mot.enrich$samp2))
df.mot.enrich[ a, "padj" ] <- p.adjust(exp(df.mot.enrich[ a, "V4" ]), method = "BH")

# Load in the python-generated HOMER motif clusters
mot.clust <- read.table("homer/homer_mot_cluster_python.txt", header = T, sep = "\t", row.names = 1)

# Add the python motif cluster information 
df.mot.enrich$motif.clust <- mot.clust[ df.mot.enrich$motif.uniq, "cluster" ]

# Load in the gene expression information and the lut for the TF names to gene names
gex <- read.csv("rna/refseqCurated_tpm.csv", header = T, row.names = 1)
tf.lut <- read.table("homer/motif_lut.txt", sep = "\t", header = T)

### NOTE: The above file "motif_lut.txt" was the effort of yours truly manually going through *every* homer motif and assigning each to 
### gene name(s) for the purpose of the below expression analysis

# Iterate through "df.mot.enrich" and determine whether the factor with the given motif is actually expressed across all samples that correspond
# to a given cell type. Yeah!  
df.mot.enrich$expr.all <- unlist( mclapply(1:nrow(df.mot.enrich), function(f){
	s <- subset(tf.lut, data == df.mot.enrich[f,"data"] & motif == df.mot.enrich[f,"motif"] & V2 == df.mot.enrich[f,"V2"])
	if ( grepl("1e", df.mot.enrich[f,"samp1"]) | grepl("1e", df.mot.enrich[f,"samp2"]) ) { a <- samp.names.1e } else { a <- samp.names.1d2 }
	if ( grepl(",", s$gene.name) ) {
		u <- unlist( a[[df.mot.enrich[f,"clust"]]] ); u <- u[u %in% colnames(gex)]
		splt <- sapply(unlist(str_split(s$gene.name, pattern = ",")), function(g){ return( sum(gex[which(gsub(".*\\|", "", rownames(gex)) %in% g), u ] >= 1) == length(u) ) })
		return( ifelse( sum( splt ) == length( splt ), T, F ) )
	} else if ( s$gene.name %in% gsub(".*\\|", "", rownames(gex)) ) {
		u <- a[[df.mot.enrich[f,"clust"]]]; u <- u[u %in% colnames(gex)]
		return( sum(gex[which(gsub(".*\\|", "", rownames(gex)) %in% s$gene.name), u ] >= 1) == length(u) )
	} else {
		return( NA )
	}
}, mc.cores = 32) )

# Save the above object so that you don't need to regenerate it in the future
write.table(df.mot.enrich, file = "homer/rObject_df.mot.enrich.txt", sep = "\t", quote = F, row.names = F)

# Now get a list of all the significant TF motifs that are enriched in accessible ORR1s # & V7 >= 10
tf.list.1e <- unique(subset(df.mot.enrich, padj <= 5e-2 & V7 >= 10 & expr.all == T & samp1 == "orr1e")$motif.uniq)

# Get the list of the top 3 TFs (each given motif cluster - as per your version of Vierstra code - is limited to top factor)
a <- sapply(1:7, function(i) {
	b <- subset(df.mot.enrich, samp1 == "orr1e" & motif.uniq %in% tf.list.1e & clust == i & padj <= 0.05) %>% group_by(motif.clust) %>% arrange(padj) %>% slice_head(n = 1) %>% ungroup() %>% arrange(desc(V5)) %>% head(3) %>% as.data.frame()
	return( as.character(b$motif.uniq) )
})
b <- unique(as.character(unlist(a))) 

# Subset the main file to get a filtered list along with necessary columns for the plot 
p <- subset(df.mot.enrich, padj <= 0.05 & motif.uniq %in% b & (samp1 == "orr1e" | samp2 == "orr1e")) %>% arrange(padj)
p$lrat <- ifelse(p$samp1 == "orr1e", p$V5, -p$V5); p.lut <- data.frame(row.names = 1:7, "name" = c("Stem\n& Prog", "B", "T", "NK", "ILC3", "MF", "DC"))

# Before you plot, turn it into a matrix to perform hierarchical clustering on and get the columns rows. Yeah
p$lab <- factor(p.lut[p$clust,"name"], levels = p.lut$name)
a <- reshape2::dcast(p, motif ~ lab, value.var = "padj"); a[is.na(a)] <- 1; row.names(a) <- a$motif; a$motif <- NULL 
dist_mat <- dist(scale(a), method = 'euclidean'); hclust_avg <- hclust(dist_mat, method = 'average')
# plot(hclust_avg)

# Set the factor order for the columns based on the hclust output. Remove the motifs "Ets1-distal" and "PU.1-IRF" as they are duplicates and
# also switch the places of "CTCF" and "EBF1" 
p$motif <- factor(p$motif, levels = hclust_avg$labels[hclust_avg$order][c(1:5,8,7,6,10:14,16)])
p <- subset(p, motif %in% hclust_avg$labels[hclust_avg$order][c(-9,-15)])

# Now reset the label levels for the plot
p$lab <- factor(p.lut[p$clust,"name"], levels = rev(p.lut$name))
f1 <- data.frame(x = factor("DC", levels = rev(p.lut$name)), y = factor("PU.1", levels = levels(p$motif)), f = factor(c("Yes", "No"), levels = c("Yes", "No")))

g2 <- ggplot(p, aes(x = motif, y = lab, size = -log10(padj), fill = lrat)) + geom_point(data = f1, aes(x = y, y = x, color = f), pch = 21, alpha = 0, inherit.aes = F) + geom_point(data = subset(p, expr.all == T), pch = 21, stroke = 1/.pt/.5, color = "black") + geom_point(data = subset(p, expr.all == F), pch = 21, color = "transparent") + theme_bw() + 
	scale_fill_gradient2(low = "blue", mid = "white", high = "dark orange", limits = c(-1.5,2.5), oob = squish, breaks = c(-1,0,1,2), labels = c("-1", "0", "1", "2"), name = expr("log"[2]*"(obs/bkgd)"~"ratio")) + 
	scale_size_area(limits = c(-log10(0.05),15), oob = squish, max_size = 5, breaks = c(-log10(0.05), 5, 10, 15), labels = c("5e-2", "1e-5", "1e-10", "<1e-15"), name = "Adj p-val") +
	scale_color_manual(name = "Expressed", labels = c("Yes", ""), values = c("black", "transparent")) +
	labs(x = "", y = "") + guides(fill = guide_colorbar(order = 1, title.position="top", title.hjust = 0.5, barwidth = unit(0.5, "in"), barheight = unit(.1, "in")), size = guide_legend(override.aes = list(pch = 21, fill = "black", color = "transparent"), order = 2, title.position="top", title.hjust = 0.5), color = guide_legend(override.aes = list(pch = 21, size = 5, fill = "white", color = c("black", "transparent"), stroke = 1/.pt/.5, alpha = 1), order = 3, title.position="top", title.hjust = 0.5)) +
	coord_fixed(clip = "off") + scale_x_discrete(position = "top", guide = guide_axis(angle = 45)) + scale_y_discrete(position = "right") +
	theme(panel.border = element_blank(), text = element_text(size = 8, color = "black"), axis.text = element_text(size = 8, color = "black"), axis.ticks = element_blank(), axis.text.y = element_text(size = 8, color = "black"), legend.position="bottom" )

# Load in the HOMER motif matrices to make their sequence logos
tmp <- monaLisa::homerToPFMatrixList("atac/homerFilt.motifs")
pfms <- sapply(tmp, function(i){return(i@profileMatrix)})
names(pfms) <- sapply(tmp, function(i){return( gsub("/Homer", "", i@ID) )}); rm(tmp)

d <- lapply(levels(p$motif), function(i) {return( pfms[[ which( names(pfms) == subset(p, motif == i)[1,"motif.uniq"] ) ]] )} )

# Make the reverse compliment of Fli1 to better match the other ETS motifs
d[[3]] <- universalmotif::motif_rc(d[[3]])@motif


# You couldn't find a way to append the sequence logos on top of the above figure, so make a separate one to combine them manually in Illustator
pdf("figures/fig2_logos.pdf")
for (f in 1:length(d)) {
	grid.newpage() 
	print(ggseqlogo(d[[f]], method = "probability", ncol = 1) + theme_void(), vp = viewport(angle = 90))
}
dev.off()

# Plot all the R-based panels. The ChIP-seq heatmaps are generated using DeepTools, and both these figures need lots of cleanup using illustrator 
pdf("figures/fig2ab.pdf", height = 4, width = 7)
plot_grid( g1, g2, labels = c("A", "B"), rel_widths = c(.85, 1), ncol = 2, label_size = 12 )
dev.off()

# Remove all the objects you no longer need
rm(df.mot.enrich, dist_mat, f1, g1, g2, hclust_avg, p, p.lut, tf.list.1e, samp.names.1e, f, gex, mot.clust, tf.lut, samp.names.1d2)


## C 

### NOTE: The code for this section is instead in the file "bashCode" since it uses DeepTools for the heatmap generation and row ordering



################
### Figure 3 ###
################

# export PATH=/home/jdc397/lbzip2-2.5/bin:$PATH; module load R/4.1.3-r9; R --vanilla

# Load all the relevant libraries you need for this section
suppressPackageStartupMessages({ 
	library(ComplexHeatmap); library(ggthemes); library(dplyr); library(ggplot2); library(cowplot); library(TFBSTools)
	library(PWMEnrich); registerCoresPWMEnrich(16); library(monaLisa); library(Biostrings); library(fastSave); library(doParallel)
	library(scatterpie)
}) 


## A

# Figure out how many ODEs are bound in at least 2 of the 6 ChIP-seq experiments, and of those how many are bound by at least PU.1 in progenitors
tmp <- c("PU.1_prog", "Pax5_B", "Runx1_CD8+", "RORg_Th17", "PU.1_MF", "IRF8_CD8+DC")
a <- read.table("chip/figure3_ode_int_chip_vals.txt", sep = "\t"); colnames(a) <- c("chr", "start", "end", "repName", "id", "cpm.c", tmp)
rownames(a) <- a$id; b <- lapply(tmp, function(i){ return( rownames(a[a[,i] >= 1 & a$repName == "ORR1E", ]) ) }); names(b) <- tmp

# Use ComplexHeatmap to generate the underlying object for the Upset plot
m <- make_comb_mat(b); m <- m[comb_size(m) >= 10]; ss <- set_size(m); cs <- comb_size(m)

# Create the color theme you'll be using, since it has plenty of good orange
col.c <- c(tableau_color_pal(type = "ordered-diverging", palette = "Sunset-Sunrise Diverging")(7), "gray69")

# Load in the ODE annotations
ode <- read.table("atac/rObject_ode.txt", sep = "\t", header = T)[,c(1:5,7)]; rownames(ode) <- ode$id

bt <- bind_rows(
	lapply(names(comb_degree(m)), function(i){
		data.frame( "id" = extract_comb(m, comb_name = i), "code" = i ) 
	})
)
bt2 <- bt

bt$cpm.c <- ode[bt$id,"cpm.c"]; bt <- table(bt$code, bt$cpm.c); bt <- matrix(bt, ncol=ncol(bt), dimnames=dimnames(bt))
bt <- bt[names(comb_degree(m)),]

us <- UpSet(m, set_order = order(ss), comb_order = order(comb_degree(m), -cs), pt_size = unit(8, "pt"), lwd = 2.86/.75/.pt,
	top_annotation = HeatmapAnnotation("ORR1E overlaps" = anno_barplot(bt, ylim = c(0, max(rowSums(bt))*1.15), 
		border = FALSE, gp = gpar(fill = col.c), height = unit(1.3, "in"), axis_param = list(gp = gpar(fontsize = 8, lwd = 2.86/.75/.pt, color = "black"))),
        	annotation_name_side = "left", annotation_name_rot = 90, annotation_name_gp= gpar(fontsize = 8, color = "black")),
	right_annotation = NULL,
	row_names_gp = gpar(fontsize = 8, color = "black"),
)
us <- draw(us); od <- column_order(us)

g1 <- grid.grabExpr( {
	draw(us);
	decorate_annotation("ORR1E overlaps", {
	    grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(3, "pt"), 
        	default.units = "native", just = c("left", "centre"), 
	        gp = gpar(fontsize = 8, col = "black"), rot = 90)
	})
})

rm( tmp, m, a, b, col.c, ode, bt, bt2, us, od, cs, ss )



## B & C

# Load in the binding information across ODEs
tmp <- c("PU.1_prog", "Pax5_B", "Runx1_CD8+", "RORg_Th17", "PU.1_MF", "IRF8_CD8+DC")
a <- read.table("chip/figure3_ode_int_chip_vals.txt", sep = "\t"); colnames(a) <- c("chr", "start", "end", "repName", "id", "cpm.c", tmp)
rownames(a) <- a$id

# Create the color theme you'll be using
col.c <- c(tableau_color_pal(type = "ordered-diverging", palette = "Sunset-Sunrise Diverging")(7), "gray69")

g <- list()

# What's the fraction of ORR1D2s bound by at least one factor?
p <- bind_rows( lapply(1:8, function(i) {
	return( data.frame( "cpm.c" = i, "bound" = sum( rowSums( subset(a, repName == "ORR1E" & cpm.c == i)[,tmp] ) > 0 ), "total" = nrow(subset(a, repName == "ORR1E" & cpm.c == i)) ) )
}) )
p$frac <- round(p$bound / p$total * 100, 2); p$cpm.c <- factor(p$cpm.c, levels = 8:1)

# Of the above, how many of them are bound by at least PU.1? 
q <- bind_rows( lapply(1:8, function(i) {
	return( data.frame( "cpm.c" = i, "bound" = nrow(subset(a, repName == "ORR1E" & cpm.c == i & PU.1_prog >= 1)), "total" = nrow(subset(a, repName == "ORR1E" & cpm.c == i)) ) )
}) )
q$frac <- round(q$bound / q$total * 100, 2); q$cpm.c <- factor(q$cpm.c, levels = 8:1)

# Run stats to see if the proportion of cell type-specific ODEs bound is significantly higher than background
sig <- bind_rows( lapply(1:7, function(i) {
	return( data.frame( "cpm.c" = i, "sig" = fisher.test( matrix(c( p[i,"bound"], p[i,"total"] - p[i,"bound"], p[8,"bound"], p[8,"total"] - p[8,"bound"] ), nrow = 2, byrow = T), alternative = "greater" )$p.value, "frac" = p[i,"frac"] ) )
}) )
sig$sig <- p.adjust(sig$sig, method = "BH"); sig$lab <- ifelse( sig$sig <= 0.05, "*", "" ); sig$cpm.c <- factor(sig$cpm.c, levels = 7:1)

r <- 86
g[[1]] <- ggplot(p, aes(x = cpm.c, y = frac, fill = cpm.c)) + geom_col(show.legend = F) + geom_col(data = q, aes(x = cpm.c, y = frac), fill = "transparent", color = "black", linewidth = 1/.75/.pt, show.legend = F) +
	scale_y_continuous(breaks = c(0,25,50,75), limits = c(0,r), labels = c("0%", "25%", "50%", "75%"), expand = c(0,0)) + scale_fill_manual(values = rev(col.c)) + theme_classic() + 
	labs(x = "", y = "% bound by at least one TF") + scale_x_discrete(labels = rev(c("Stem\n& Prog", "B", "T", "NK", "ILC3", "MF", "DC", "Bkgd"))) + geom_text(aes(x = cpm.c, y = frac + (r * 0.07), label = bound), size = 8/.pt, color = "black", angle = 270) + 
	geom_text(data = sig, aes(x = cpm.c, y = frac + (r * 0.13), label = lab), angle = 270, color = "dark orange", size = 8/.pt, fontface = "bold", inherit.aes = F) + geom_text(data = q, aes(x = cpm.c, y = frac - (r * 0.07), label = bound), size = 8/.pt, color = "black", angle = 270) + 
	theme(axis.title = element_text(size = 8, color = "black"), axis.text = element_text(size = 8, color = "black"), axis.ticks.length.x = unit(.04, "in"), axis.ticks.y = element_blank(), axis.ticks.x = element_line(linewidth = 1/.75/.pt, color = "black"), axis.line = element_line(linewidth = 1/.75/.pt, color = "black")) + 
	guides(fill = guide_legend(reverse=TRUE)) + coord_flip()


# What's the fraction of ORR1Es bound by at least two factors? 
p <- bind_rows( lapply(1:8, function(i) {
	return( data.frame( "cpm.c" = i, "bound" = length(which(rowSums(subset(a, repName == "ORR1E" & cpm.c == i)[,7:12] >= 1) >= 2)), "total" = nrow(subset(a, repName == "ORR1E" & cpm.c == i)) ) )
}) )
p$frac <- round(p$bound / p$total * 100, 2); p$cpm.c <- factor(p$cpm.c, levels = 8:1)

# Of the above, how many of them are bound by at least PU.1? 
q <- bind_rows( lapply(1:8, function(i) {
	return( data.frame( "cpm.c" = i, "bound" = length(which(rowSums(subset(a, repName == "ORR1E" & cpm.c == i & PU.1_prog >= 1)[,8:12] >= 1) >= 1)), "total" = nrow(subset(a, repName == "ORR1E" & cpm.c == i)) ) )
}) )
q$frac <- round(q$bound / q$total * 100, 2); q$cpm.c <- factor(q$cpm.c, levels = 8:1)

# Run stats to see if the proportion of cell type-specific ODEs bound is significantly higher than background
sig <- bind_rows( lapply(1:7, function(i) {
	return( data.frame( "cpm.c" = i, "sig" = fisher.test( matrix(c( p[i,"bound"], p[i,"total"] - p[i,"bound"], p[8,"bound"], p[8,"total"] - p[8,"bound"] ), nrow = 2, byrow = T), alternative = "greater" )$p.value, "frac" = p[i,"frac"] ) )
}) )
sig$sig <- p.adjust(sig$sig, method = "BH"); sig$lab <- ifelse( sig$sig <= 0.05, "*", "" ); sig$cpm.c <- factor(sig$cpm.c, levels = 7:1)

r <- 62
g[[2]] <- ggplot(p, aes(x = cpm.c, y = frac, fill = cpm.c)) + geom_col(show.legend = F) + geom_col(data = q, aes(x = cpm.c, y = frac), fill = "transparent", color = "black", linewidth = 1/.75/.pt, show.legend = F) +
	scale_y_continuous(breaks = c(0,20,40,60), limits = c(0,r), labels = c("0%", "20%", "40%", "60%"), expand = c(0,0)) + scale_fill_manual(values = rev(col.c)) + theme_classic() + 
	labs(x = "", y = "% bound by at least two TFs") + scale_x_discrete(labels = rep("", 8) ) + geom_text(aes(x = cpm.c, y = frac + (r * 0.045), label = bound), size = 8/.pt, color = "black", angle = 270) + 
	geom_text(data = sig, aes(x = cpm.c, y = frac + (r * 0.09), label = lab), angle = 270, color = "dark orange", size = 8/.pt, fontface = "bold", inherit.aes = F) + geom_text(data = q, aes(x = cpm.c, y = frac - (r * 0.045), label = bound), size = 8/.pt, color = "black", angle = 270) + 
	theme(axis.title = element_text(size = 8, color = "black"), axis.text = element_text(size = 8, color = "black"), axis.ticks.length.x = unit(.04, "in"), axis.ticks.y = element_blank(), axis.ticks.x = element_line(linewidth = 1/.75/.pt, color = "black"), axis.line = element_line(linewidth = 1/.75/.pt, color = "black")) + 
	guides(fill = guide_legend(reverse=TRUE)) + coord_flip()

g1.5 <- plot_grid( g1, plot_grid( plotlist = g, labels = c("B", "C"), ncol = 2, rel_widths = c(1,.85), label_size = 12 ), nrow = 1, labels = c("A", ""), label_size = 12 )

rm(col.c, a, g, p, q, r, sig, tmp, g1)



## D

### NOTE: Below code is for generating piecharts to incorporate into the output of the phylogenetic tree using Iroki (https://www.iroki.net/),
### 	which is a very nice tree viewer. The piecharts and tree were combined in Illustrator otherwise

# Load in the list of ODEs and ODE-derived peaks
ode <- read.table("atac/rObject_ode.txt", sep = "\t", header = T); ode$cpm.c <- factor(ode$cpm.c, levels = 1:8); ode$peak <- gsub(".*\\|", "", ode$pid)
rownames(ode) <- ode$id

# For each clade, you want to figure out how many of each ATAC-seq cluster there are to make the piecharts. Get started!
a <- read.table("figures/r1e_iroki_metadata.txt", header = T); a$cpm.c <- ode[gsub("_", ":", a$name), "cpm.c"]
p <- data.frame(clade = 1:max(as.numeric(gsub("k_", "", a$branch_color))))
for (f in 1:8) {
	p[,LETTERS[f]] <- sapply(1:nrow(p), function(i){nrow(subset(a, branch_color == paste0("k_", i) & cpm.c == f))})
}
p$tot <- sapply(1:nrow(p), function(i){sum(p[i,2:ncol(p)])})

# Create the color theme you're using
col.c <- c(tableau_color_pal(type = "ordered-diverging", palette = "Sunset-Sunrise Diverging")(7), "gray69") 

pdf("figures/fig3d_phyloClades_piecharts.pdf", height = 2.5, width = 5)
ggplot() + geom_scatterpie(aes(x = clade*2, y = 4, r = 0.75), data = subset(p, clade < 7), cols = LETTERS[1:8], show.legend = F) + 
	geom_scatterpie(aes(x = (clade-6)*2, y = 1, r = 0.75), data = subset(p, clade > 6), cols = LETTERS[1:8], show.legend = F) + coord_fixed() + 
	scale_fill_manual(values = col.c, labels = c("A" = "Stem\n& Prog", "B" = "B", "C" = "T", "D" = "NK", "E" = "ILC3", "F" = "MF", "G" = "DC", "H" = "Bkgd"), name = "ATAC-seq\nCluster") + theme_void() +
	geom_text(data = subset(p, clade < 7), aes(x = clade*2, y = 5, label = tot), fontface = "bold", size = 8/.pt) + geom_text(data = subset(p, clade > 6), aes(x = (clade-6)*2, y = 2, label = tot), fontface = "bold", size = 8/.pt) + 
	geom_text(data = subset(p, clade < 7), aes(x = clade*2, y = 3, label = clade), fontface = "bold", size = 8/.pt) + geom_text(data = subset(p, clade > 6), aes(x = (clade-6)*2, y = 0, label = clade), fontface = "bold", size = 8/.pt) + 
	ylim(0,7)
dev.off()

# Include the proportion of all elements assigned to clusters to get a sense of the overall distribution
p1 <- data.frame(clade = "all")
for (f in 1:8) {
	p1[,LETTERS[f]] <- nrow(subset(a, cpm.c == f))
}
p1$tot <- sum(p1[,2:ncol(p1)])

pdf("figures/fig3d_all_piechart.pdf", height = 3, width = 5)
ggplot() + geom_scatterpie(aes(x = 1, y = 1, r = 0.75), data = p1, cols = LETTERS[1:8]) + coord_fixed() + 
	scale_fill_manual(values = col.c, labels = c("A" = "Stem\n& Prog", "B" = "B", "C" = "T", "D" = "NK", "E" = "ILC3", "F" = "MF", "G" = "DC", "H" = "Bkgd"), name = "ATAC-seq\ncluster") + theme_void() +
	geom_text(data = p1, aes(x = 1, y = 1.9, label = tot), size = 8/.pt) + ylim(0,3) + theme(legend.text = element_text(size = 8, color = "black"))
dev.off()

# Remove objects you no longer need 
rm(a, f, p, p1)



## E 

# export PATH=/home/jdc397/lbzip2-2.5/bin:$PATH; module load R/4.1.3-r9; R --vanilla

# Load in the packages you'll need throughout
suppressPackageStartupMessages({ 
	library(PWMEnrich); registerCoresPWMEnrich(16); library(Biostrings); library(fastSave); library(ggplot2)
}) 


tmp <- monaLisa::homerToPFMatrixList("homer/homerFilt.motifs")
mot.thresh <- sapply(tmp, function(i){return(log(2^i@tags$log2cut))}); pfms <- sapply(tmp, function(i){return(i@profileMatrix)})
names(pfms) <- names(mot.thresh) <- sapply(tmp, function(i){return( gsub("/Homer", "", i@ID) )})
meta <- data.frame(row.names = sapply(tmp, function(i){return(i@ID)}), "log.thresh" = mot.thresh, "proto.thresh" = mot.thresh * .5, "mot.length" = sapply(tmp, function(i){return(ncol(i@profileMatrix))}))

# Load in the consensus sequences for all of the ORR1 and related subfamily consensus sequences
cons <- readDNAStringSet("rm/dfam3.6_cons.fa"); cons <- cons[grepl("ORR1E|ORR1D2|ORR1F|ORR1B1#|ORR1B2#|MTD", names(cons)),]
names(cons) <- gsub("#.*", "", names(cons))

# Use the object "pfms" to calculate the scores for each motif over each consensus sequence chosen. Woot woot
scores.con <- motifScores(cons, pfms, raw.scores=T)
save.lbzip2(scores.con, file = "homer/rObject_scores.con.RDataFS", n.cores = 16)

# Load in the python-generated HOMER motif clusters
mot.clust <- read.table("homer/homer_mot_cluster_python.txt", header = T, sep = "\t", row.names = 1)

# Use the "meta" object to create a "base.mot" object that contains all the above consensus sequences
base.con.mot <- data.frame()
for (f in 1:length(scores.con)) {
	for (g in 1:ncol(scores.con[[f]])) {
		
		a <- which( log(scores.con[[f]][,g]) >= meta[ colnames(scores.con[[f]])[g], "proto.thresh" ] )
		if ( length(a) != 0 ) {
			tmp <- data.frame( "seq" = gsub("#.*", "", names(scores.con)[f]), "mot" = colnames(scores.con[[f]])[g], "motif" = gsub( "\\(.*", "", colnames(scores.con[[f]])[g] ),
				"start" = ifelse( a > nrow(scores.con[[f]]) / 2, a - (nrow(scores.con[[f]]) / 2), a ),
				"end" = ifelse( a > nrow(scores.con[[f]]) / 2, a - (nrow(scores.con[[f]]) / 2) + (meta[colnames(scores.con[[f]])[g],"mot.length"]) - 1, a + (meta[colnames(scores.con[[f]])[g],"mot.length"]) - 1 ),
				"strand" = ifelse( a > nrow(scores.con[[f]]) / 2, "-", "+" ), "score" = log(scores.con[[f]][a,g]), 
				"qual" = ifelse( log(scores.con[[f]][a,g]) >= meta[colnames(scores.con[[f]])[g],"log.thresh"], "mot", ifelse( log(scores.con[[f]][a,g]) >= meta[colnames(scores.con[[f]])[g],"proto.thresh"], "proto", "no" ) ) )
			base.con.mot <- rbind(base.con.mot, tmp)
		}
	}
}
rm(f,g); rownames(base.con.mot) <- 1:nrow(base.con.mot); base.con.mot$mot.clust <- mot.clust[ base.con.mot$mot, "cluster" ]

# Save the above object to not have to generate again
write.table(base.con.mot, file = "homer/rObject_base.con.mot.txt", row.names = F, quote = F, sep = "\t")

aa <- c("ORR1B1" = 6, "ORR1B2" = 5, "ORR1D2" = 4, "ORR1E" = 3, "ORR1F" = 2, "MTD" = 1)
bb <- data.frame(row.names = rev(c("ORR1B1", "ORR1B2", "ORR1D2", "ORR1E", "ORR1F", "MTD")), xmin = 1, xmax = width(cons[rev(c("ORR1B1", "ORR1B2", "ORR1D2", "ORR1E", "ORR1F", "MTD"))]))
dd <- subset(base.con.mot, motif == "PU.1")

g3 <- ggplot(bb) + geom_rect(aes(xmin = xmin, xmax = xmax, ymin = (1:length(aa))-0.05, ymax = (1:length(aa))+0.05 )) +
	geom_segment(data = subset(dd, qual == "proto"), aes(x = (start+end)/2, y = ifelse(strand == "+", aa[seq]+0.05, aa[seq]-0.05), yend = ifelse(strand == "+", aa[seq]+0.4, aa[seq]-0.4), color = score), linewidth = 1.4, inherit.aes = F, show.legend = F) +
	geom_point(data = subset(dd, qual == "mot" & strand == "+"), aes(x = (start+end)/2, y = aa[seq]+0.35, fill = score), shape = 25, size = 2, color = "black", stroke = 1/.5/.pt, inherit.aes = F) +
	geom_point(data = subset(dd, qual == "mot" & strand == "-"), aes(x = (start+end)/2, y = aa[seq]-0.35, fill = score), shape = 24, size = 2, color = "black", stroke = 1/.5/.pt, inherit.aes = F) +
	theme_classic() + scale_y_continuous(breaks = 1:length(aa), labels = rev(names(aa))) + labs(y = "Consensuses", x = "Consensus length")  +
	scale_fill_gradient(name = "PU.1 motif\nstrength", low = "black", high = "dark orange", limits = c(0,max(dd$score))) + scale_color_gradient(low = "black", high = "dark orange", limits = c(0,max(dd$score))) +
	guides(color = element_blank(), fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) + 
	theme(legend.key.height = unit(0.15, "in"), legend.key.width = unit(0.1, "in"), text = element_text(size = 8, color = "black"), axis.ticks.length=unit(.04, "in"), axis.ticks = element_line(linewidth = 1/.75/.pt), axis.line = element_line(linewidth = 1/.75/.pt, color = "black"), axis.text = element_text(size = 8, color = "black"), legend.text = element_text(size = 8, color = "black") )


## F

# export PATH=/home/jdc397/lbzip2-2.5/bin:$PATH; module load R/4.1.3-r9; R --vanilla

# Load in the packages you'll need throughout
suppressPackageStartupMessages({ 
	library(Biostrings); library(doParallel); library(fastSave); library(PWMEnrich); registerCoresPWMEnrich(16); library(ggplot2)
}) 

# Load in the consensus sequences and combined ODE .fa files
cons <- readDNAStringSet("rm/dfam3.6_cons.fa"); cons <- cons[grepl("ORR1E|ORR1D2", names(cons)),]
comb <- readDNAStringSet("rm/ode.fa")

# Read in the DNA substitution matrix (that you found from multiple alignment softwares, such as MEME)
dnafull <- as.matrix(read.table("rm/DNAfull.mat", header = T))

res <- mclapply(1:length(comb), function(i) {
	test <- pairwiseAlignment(comb[i], cons[which(grepl(gsub(".*\\|", "", names(comb)[i]), names(cons)))], substitutionMatrix = dnafull, gapOpening = 10, gapExtension = -0.5, type = "global")
	val <- as.data.frame(matrix(0,nrow=2,ncol=width(alignedPattern(test))), row.names = c(names(comb[i]), "cons"))
	val[1,] <- as.matrix(alignedPattern(test)); val[2,] <- as.matrix(alignedSubject(test))
	return(val)
}, mc.cores = 32)
names(res) <- names(comb)

lut <- as.data.frame(matrix(-1,nrow=length(comb),ncol=max(width(comb))), row.names = names(comb))
lut <- bind_rows( mclapply(1:nrow(lut), function(f) {
	cc <- 1; cb <- 1; val <- lut[f,]; clen <- width(cons[grepl(gsub("_.*", "", gsub(".*\\|", "", rownames(lut)[f])), names(cons))])
	for ( g in 1:ncol(res[[f]]) ) { 
		if ( res[[f]][1,g] == "-" ) { # Sequence is gapped, so increment the consensus counter
			cc <- cc + 1
		} else { # Sequence is not gapped. Append the matricies			
			if (cc > clen) { cc <- clen } # Make sure consensus doesn't go above its maximum value (when sequence extends past consensus?)

			if ( res[[f]][2,g] == "-" ) { # Consensus is gapped, add NA instead
				val[1,cb] <- NA;
			} else { # No gaps, so increment consensus counter
				val[1,cb] <- cc; cc <- cc + 1
			}
			# With another base dealt with, append the current base counter
			cb <- cb + 1
		}
	}
	return(val)
}, mc.cores = 16) )

# Write the above objects to disk so you don't have to regenerate them in the future
save.lbzip2(res, file = "rm/rObject_res.RDataFS", n.cores = 16)
save.lbzip2(lut, file = "rm/rObject_lut.RDataFS", n.cores = 16)

# Load in the python-generated HOMER motif clusters
mot.clust <- read.table("homer/homer_mot_cluster_python.txt", header = T, sep = "\t", row.names = 1)

# Load in the HOMER motif file and get the tresholds you need to call a given motif
tmp <- monaLisa::homerToPFMatrixList("homer/homerFilt.motifs")
mot.thresh <- sapply(tmp, function(i){return(log(2^i@tags$log2cut))}); pfms <- sapply(tmp, function(i){return(i@profileMatrix)})
names(pfms) <- names(mot.thresh) <- sapply(tmp, function(i){return( gsub("/Homer", "", i@ID) )})
meta <- data.frame(row.names = sapply(tmp, function(i){return(i@ID)}), "log.thresh" = mot.thresh, "proto.thresh" = mot.thresh * .5, "mot.length" = sapply(tmp, function(i){return(ncol(i@profileMatrix))}))

# Figure out the location of each motif for each ODE element
scores.ode <- motifScores(comb, pfms, raw.scores=T); rm(tmp)

# Save the above object to disk so you don't have to regenerate it! 
# WARNING: Will be an ~14Gb sized object (that's R lists for you) so take caution
save.lbzip2(scores.ode, file = "homer/rObject_scores.ode.RDataFS", n.cores = 16)

# Generate an object to contain all of the motifs and protomotif sequences
b <- data.frame()
base.mot <- bind_rows(mclapply(1:length(scores.ode), function(f) {
	for (g in 1:ncol(scores.ode[[f]])) {
		
		a <- which( log(scores.ode[[f]][,g]) >= meta[ colnames(scores.ode[[f]])[g], "proto.thresh" ] )
		if ( length(a) != 0 ) {
			tmp <- data.frame( "seq" = gsub("#.*", "", names(scores.ode)[f]), "mot" = colnames(scores.ode[[f]])[g], "start" = ifelse( a > nrow(scores.ode[[f]]) / 2, a - (nrow(scores.ode[[f]]) / 2), a ),
				"end" = ifelse( a > nrow(scores.ode[[f]]) / 2, a - (nrow(scores.ode[[f]]) / 2) + (meta[colnames(scores.ode[[f]])[g],"mot.length"]) - 1, a + (meta[colnames(scores.ode[[f]])[g],"mot.length"]) - 1 ),
				"strand" = ifelse( a > nrow(scores.ode[[f]]) / 2, "-", "+" ), "score" = log(scores.ode[[f]][a,g]), 
				"qual" = ifelse( log(scores.ode[[f]][a,g]) >= meta[colnames(scores.ode[[f]])[g],"log.thresh"], "mot", "proto" ) )
			b <- rbind(b, tmp)
		}
	}
	b
}, mc.cores = 32))
rownames(base.mot) <- 1:nrow(base.mot); base.mot$mot.clust <- mot.clust[ base.mot$mot, "clust" ]

# Save the above object
save.lbzip2(base.mot, file = "homer/rObject_base.mot.RDataFS", n.cores = 16)

# Load in the list of ODEs and ODE-derived peaks
ode <- read.table("atac/rObject_ode.txt", sep = "\t", header = T); ode$cpm.c <- factor(ode$cpm.c, levels = 1:8); ode$peak <- gsub(".*\\|", "", ode$pid)
rownames(ode) <- ode$id

# Create an object that contains the sequences at each position relative to their respective consensus 
tmp <- as.data.frame(matrix("-", nrow = 1, ncol = max(width(cons))))
sa.base <- bind_rows(mclapply(names(comb), function(i) {
	a <- tmp
	
	# For the length of the current element, assign the elements to the matrix using "lut"
	for (f in 1:sum(res[[i]][1,] != "-")) {
		if( !is.na(lut[i,f]) ) { a[1, lut[i,f]] <- res[[i]][1,lut[i,f]] }
	}
	a
}, mc.cores = 32) )
rownames(sa.base) <- names(comb)

# Save the above object
write.table(sa.base, "rm/rObject_sa.base.txt", quote = F, col.names = F, sep = "\t")

# Create the object whose data will underly the plot, but also for every single factor
ode.con.mot.per <- mclapply( names(pfms), function(i) {

	a <- subset(base.mot, grepl("1E", seq) & qual == "mot" & mot == i); a$cpm.c <- ode[a$seq, "cpm.c"]

	tmp <- matrix(0, nrow = 8, ncol = 359)
	for (f in 1:8) {
		b <- subset(a, cpm.c == f)
		if ( nrow(b) != 0 ) { 
			for (g in 1:nrow(b)) {
				d <- b[g,"seq"]; e <- lut[ d, b[ g, "start" ]:b[ g, "end" ] ]; e <- e[!is.na(e)] 
				for (h in e) { tmp[f, h] <- tmp[f, h] + 1 } 
			}
		}

		# Now, divide the total number of instances in each position by the total number of elements with sequence
		tmp[f,] <- tmp[f,] / colSums(sa.base[ subset(ode, cpm.c == f & repName == "ORR1E")$id, 1:359 ] != "-")
	}

	return( tmp )
}, mc.cores = 16)
names(ode.con.mot.per) <- names(pfms)

# Save the above object to not have to make it again
save.lbzip2(ode.con.mot.per, file = "homer/rObject_ode.con.mot.per.RDataFS", n.cores = 16)

# Load in the object containing all the motif sites for each relevant consensus sequence
base.con.mot <- read.table("homer/rObject_base.con.mot.txt", header = T, sep = "\t")

# Create a named vector to use for renaming the cluster numbers into something meaningful! 
lut.vec <- c("1" = "Stem\n& Prog", "2" = "B", "3" = "T", "4" = "NK", "5" = "ILC3", "6" = "MF", "7" = "DC", "8" = "Bkgd")

# Set up the objects that will underly the ggplot, specifically the PU.1 motif
a <- which(grepl("PU.1\\(", names(pfms))) 
base.p <- reshape2::melt(ode.con.mot.per[[a]]); colnames(base.p) <- c("cluster", "x", "frac")
base.p$fc <- base.p$frac - rep(subset(base.p, cluster == 8)$frac, each = 8)
base.p$lab <- factor(lut.vec[base.p$clust], levels = lut.vec)

val <- as.data.frame(matrix("", nrow = 1, ncol = 359))
a <- subset(base.con.mot, mot.clust == 33 & seq == "ORR1E" & qual == "proto")
for (f in 1:nrow(a)) { 	val[ 1, a[f,"start"]:a[f,"end"] ] <- "p" }

a <- subset(base.con.mot, mot.clust == 33 & seq == "ORR1E" & qual == "mot")
for (f in 1:nrow(a)) { 	val[ 1, a[f,"start"]:a[f,"end"] ] <- "m" }

base.rect <- data.frame(); cur <- ""; st <- 1.5
for (f in 1:ncol(val)) {
	if ( f == ncol(val) ) {	base.rect <- rbind(base.rect, data.frame( "mot.clust" = 33, "xmin" = st - 0.5, "xmax" = f, "qual" = cur )); next }
	if ( val[1,f] != cur ) { 
		base.rect <- rbind(base.rect, data.frame( "mot.clust" = 33, "xmin" = st - 0.5, "xmax" = f - 0.5, "qual" = cur ))
		st <- f; cur <- val[1,f]
	}
}
base.rect <- subset(base.rect, qual != "" & xmax - xmin > 1)

g4 <- ggplot(subset(base.p, cluster %in% c(3,5:6)), aes(x = x, y = fc * 100)) + geom_hline(yintercept = 0, linewidth = 0.5/.75/.pt, color = "black") + geom_col() + 
	geom_rect(data = base.rect, aes(xmin = xmin, xmax = xmax, ymin = -6, ymax = 21, fill = qual), color = "transparent", alpha = 0.25, inherit.aes = F, show.legend = F) + 
	facet_wrap(~lab, ncol = 1) + theme_minimal() + scale_fill_manual(values = c("dark orange", "dodgerblue")) +
	scale_x_continuous(limits = c(0.5,359.5), breaks = c(100,200,300), expand = c(0,0)) + ylab("PU.1 motif presence\ncluster % - Bkgd %") + xlab("Position in ORR1E consensus") + 
	theme(text = element_text(size = 8, color = "black"), strip.text.x = element_text(size = 8, color = "black"), axis.title = element_text(size = 8, color = "black"), axis.text = element_text(size = 8, color = "black"), panel.grid.major = element_line(linewidth = 0.5/.75/.pt), panel.grid.minor = element_line(linewidth = 0.25/.75/.pt))

g3.5 <- plot_grid( NULL, plot_grid( g3, g4, ncol = 1, labels = c("E", "F"), label_size = 12, rel_heights = c(1,1.25) ), labels = c("D", ""), label_size = 12, ncol = 2 )

# Time to finally combine everything together
pdf("figures/fig3acef.pdf", height = 6.5, width = 7)
plot_grid( g1.5, g3.5, ncol = 1, labels = c("", ""), rel_heights = c(1, 1.5) )
dev.off()

# You couldn't find a good way to put a small plot showing RORg inside of figure 3, so do it manually in Illustrator
a <- which(names(ode.con.mot.per) == "RORgt(NR)/EL4-RORgt.Flag-ChIP-Seq(GSE56019)")
b <- data.frame("frac" = apply(ode.con.mot.per[[ a ]][,121:130], 1, max) * 100, "clus" = lut.vec); b$clus <- factor(b$clus, levels = rev(b$clus))

# Create the color theme you're using
col.c <- c(tableau_color_pal(type = "ordered-diverging", palette = "Sunset-Sunrise Diverging")(7), "gray69") 
names(col.c) <- lut.vec

pdf("figures/fig3f_subPanel.pdf", height = 5, width = 3.5)
ggplot(subset(b, clus %in% c("T", "ILC3", "MF")), aes(x = frac, y = clus, fill = clus, label = round(frac,1))) + geom_col(show.legend = F) + 
	scale_fill_manual(values = col.c) + theme_classic() + ylab("") + xlab("% with ROR motif") + 
	scale_x_continuous(breaks = c(0,25,50), limits = c(0,70), labels = c(0,25,50), expand = c(0,0)) + 
	geom_text(nudge_x = .5, hjust = 0, size = 8/.pt, color = "black") + 
	theme(axis.title = element_text(size = 8, color = "black"), axis.text = element_text(size = 8, color = "black"), axis.ticks.length.x = unit(.04, "in"), axis.ticks.y = element_blank(), axis.ticks.x = element_line(linewidth = 1/.75/.pt, color = "black"), axis.line = element_line(linewidth = 1/.75/.pt, color = "black"))
dev.off()



################
### Figure 4 ###
################

# export PATH=/home/jdc397/lbzip2-2.5/bin:$PATH; module load R/4.1.3-r9; R --vanilla

suppressPackageStartupMessages({
        library(data.table); setDTthreads(8); library(ggplot2); library(dplyr); library(fastSave); library(doParallel); library(ggthemes)
	library(fgsea); library(GenomicRanges)
})

## B

### NOTE: The code for generating the "cor.mat" object can be found in the file "Rcode_cor.r", since it takes a while with lots of cores 
### and memory. In either case, it creates and saves a file containing ~460 million rows... Oy vey. In either case, it is required for 
### this section and beyond, and also generates underlying objects that aren't created below

# Load in the big correlation object 
load.lbzip2(file = "atac/rObject_cor.mat.RDataFS", n.cores = 16)

# Now split up the above object for downstream analyses
cor.mat.filt <- subset(cor.mat, filt.abc == T | filt.s == T)
save.lbzip2(cor.mat.filt, file = "atac/rObject_cor.mat.filt.RDataFS", n.cores = 16)

cor.mat.mil <- subset(cor.mat, abs(dist) <= 1e6)
save.lbzip2(cor.mat.mil, file = "atac/rObject_cor.mat.mil.RDataFS", n.cores = 16)

# Identify proportions of peaks that passed the thresholding across ORR1E and everything else
p <- as.data.frame(table(subset( cor.mat.mil, g1 != "orr1d2" )$g2))$Freq; names(p) <- c(1:8, "nonte", "te")
meta <- rbind( data.frame("g" = factor(c(1:8, "nonte", "te"), levels = c(1:8, "nonte", "te")), 
	"frac" = as.data.frame(table(subset(cor.mat.mil, filt.s == T & g1 != "orr1d2" & cor.s >= 0)$g2))$Freq / p, "num" = as.data.frame(table(subset(cor.mat.mil, filt.s == T & g1 != "orr1d2" & cor.s >= 0)$g2))$Freq, "dir" = "pos"), 
	data.frame("g" = factor(c(1:8, "nonte", "te"), levels = c(1:8, "nonte", "te")), "frac" = as.data.frame(table(subset(cor.mat.mil, filt.s == T & g1 != "orr1d2" & cor.s < 0)$g2))$Freq / p, "num" = as.data.frame(table(subset(cor.mat.mil, filt.s == T & g1 != "orr1d2" & cor.s < 0)$g2))$Freq, "dir" = "neg"))
meta$dir <- factor(meta$dir, levels = c("neg", "pos"))

b <- meta["8","num"]; b2 <- meta["81","num"]
sig <- rbind( data.frame("g" = factor(c(1:8, "nonte", "te"), levels = c(1:8, "nonte", "te")), "dir" = "pos", "rel" = sapply(c(1:8, "te", "nonte"), function(i){ ifelse( subset(meta, dir == "pos" & g == i)$frac >= subset(meta, dir == "pos" & g == 8)$frac, "up", "down" )}), "sig" = c( sapply( 1:7, function(i){ a <- nrow(subset(cor.mat.mil, g1 == "orr1e" & g2 == i & filt.s == T & cor.s >= 0)); fisher.test( matrix(c( a, p[i] - a, b, p[8] - b ), nrow = 2), alternative = "two.sided" )$p.value } ), NA, sapply( c("te", "nonte"), function(i){ a <- nrow(subset(cor.mat.mil, g2 == i & filt.s == T & cor.s >= 0)); fisher.test( matrix(c( a, p[i] - a, b, p[8] - b ), nrow = 2), alternative = "two.sided" )$p.value } )) ),
	data.frame("g" = factor(c(1:8, "nonte", "te"), levels = c(1:8, "nonte", "te")), "dir" = "neg", "rel" = sapply(c(1:8, "te", "nonte"), function(i){ ifelse( subset(meta, dir == "neg" & g == i)$frac >= subset(meta, dir == "neg" & g == 8)$frac, "up", "down" )}), "sig" = c( sapply( 1:7, function(i){ a <- nrow(subset(cor.mat.mil, g1 == "orr1e" & g2 == i & filt.s == T & cor.s < 0)); fisher.test( matrix(c( a, p[i] - a, b2, p[8] - b2 ), nrow = 2), alternative = "two.sided" )$p.value } ), NA, sapply( c("te", "nonte"), function(i){ a <- nrow(subset(cor.mat.mil, g2 == i & filt.s == T & cor.s >= 0)); fisher.test( matrix(c( a, p[i] - a, b, p[8] - b ), nrow = 2), alternative = "two.sided" )$p.value } )) ))
sig$sig <- p.adjust(sig$sig, method = "BH")
sig$lab <- ifelse( sig$sig <= 0.05, "*", "" ); sig$color <- ifelse( sig$rel == "up", "do", "b" ); sig$dir <- factor(sig$dir, levels = c("neg", "pos"))

g1 <- ggplot(meta, aes(x = frac * 100, y = rev(g), fill = dir)) + geom_col(position = position_dodge(), width = .6) +
	scale_fill_manual(values = c("pos" = "dark orange", "neg" = "black"), labels = c("pos" = "Pos", "neg" = "Neg"), name = "Parity") + 
	scale_x_continuous(labels = seq(0,2,0.5), breaks = seq(0,2,0.5), expand = c(0,0)) +
	geom_text(aes(x = (frac * 100) + 0.02, y = rev(g), label = num), hjust = "left", position = position_dodge(.65), size = 8/.pt) + theme_classic() +
	coord_cartesian(xlim = c(0,2.4), clip = "off") + geom_text(data = sig, aes(x = -0.1, y = rev(g), label = lab, color = color), hjust = "left", fontface = "bold", position = position_dodge(.65), size = 8/.pt, vjust = 0.8, show.legend = F) + 
	annotate(geom = "text", x = -.15, y = rev(sig$g[1:10]), label = c("Stem\n& Prog", "B", "T", "NK", "ILC3", "MF", "DC", "Bkgd", "nonTE", "TE"), size = 8/.pt, hjust = "right") + 
	scale_color_manual(values = c("do" = "dark orange", "b" = "blue")) + labs(x = "Filtered peak-gene pairs (%)", y = "ORR1E ATAC-seq clusters") +
	guides( fill = guide_legend(reverse = T) ) + scale_y_discrete(labels = rep("", 10)) +
	theme( text = element_text(size = 8, color = "black"), axis.title.y = element_text(size = 8, color = "black", margin = margin(r = 35)), axis.text = element_text(size = 8, color = "black"), legend.key.size = unit(.08, 'in'), axis.ticks.length = unit(.04, "in"), axis.ticks.y = element_blank(), axis.ticks.x = element_line(linewidth = 1/.75/.pt, color = "black"), legend.text = element_text(size = 8, color = "black") )

# Remove all the objects you no longer need
rm(b, b2, cor.mat.mil, meta, p, sig)


## C

# Load in the filtered R correlation object using the "fastSave" library
load.lbzip2(file = "atac/rObject_cor.mat.filt.RDataFS", n.cores = 16)

# Load in the list of ODEs and ODE-derived peaks
ode <- read.table("atac/rObject_ode.txt", sep = "\t", header = T); ode$cpm.c <- factor(ode$cpm.c, levels = 1:8); ode$peak <- gsub(".*\\|", "", ode$pid)
rownames(ode) <- ode$id

# Load in the list of TE-derived peak ids 
tes <- unique(read.table("atac/peak.te.int")$V4)

# Load in all ImmGen ATAC-seq peaks
f.m <- as.data.frame(fread("atac/summit_unif_peaks_counts.txt"))[,1:4]; colnames(f.m) <- c("chr", "start", "end", "id")
rownames(f.m) <- with(f.m, paste0(chr, ":", start, "-", end, "|", id)); f.m$class <- ifelse( f.m$id %in% gsub(".*\\|", "", ode$pid), "ode", ifelse( f.m$id %in% tes, "te", "nonte" ) )

# Create the color theme you'll be using
col.c <- c(tableau_color_pal(type = "ordered-diverging", palette = "Sunset-Sunrise Diverging")(7), "gray69", "gray47", "gray25")

# Calculate the enrichment ratio for each cluster for each ODE. Formula is: (# of DNA-DNA contacts in group / total DNA-DNA contacts) / (total bp in group / genome size)
gs <- 2725537669
tc <- sum(cor.mat.filt$filt.abc == T)

b <- bind_rows( lapply(1:8, function(i) {
	return( data.frame( "g" = i, "total" = nrow(subset(cor.mat.filt, filt.abc == T & g1 == "orr1e" & g2 == i)), "bp" = sum(with(subset(ode, cpm.c == i & repName == "ORR1E"), end - start + 1)) ) )
}) )
b <- rbind( b, data.frame( "g" = "nonTE", "total" = nrow(subset(cor.mat.filt, filt.abc == T & g1 == "nonte")), "bp" = sum(with(subset(f.m, class == "nonte"), end - start + 1)) ) )
b <- rbind( b, data.frame( "g" = "TE", "total" = nrow(subset(cor.mat.filt, filt.abc == T & g1 == "te")), "bp" = sum(with(subset(f.m, class == "te"), end - start + 1)) ) )
b$score <- with( b, (total / tc) / (bp / gs) )

b$g <- factor(b$g, levels = c("TE", "nonTE", 8:1))

g2 <- ggplot(b, aes(x = score, y = g, fill = g)) + geom_col(show.legend = F) + scale_fill_manual(values = rev(col.c)) + 
	geom_text(aes(x = score + (90 * .02), y = g, label = total), size = 8/.pt, color = "black", hjust = 0) + theme_classic() +
	scale_x_continuous(name = "T cell DNA-DNA contact\nenrichment score", limits = c(0,90), breaks = seq(0,80,20), expand = c(0,0)) +
	scale_y_discrete(name = "ORR1E ATAC-seq clusters", labels = rev(c("Stem\n& Prog", "B", "T", "NK", "ILC3", "MF", "DC", "Bkgd", "nonTE", "TE"))) + 
	theme( text = element_text(size = 8, color = "black"), axis.text = element_text(size = 8, color = "black"), axis.ticks.length = unit(.04, "in"), axis.ticks.y = element_blank(), axis.ticks = element_line(linewidth = 1/.75/.pt, color = "black"), axis.line = element_line(linewidth = 1/.75/.pt, color = "black") )

# Remove all the objects you no longer need
rm(b, col.c, cor.mat.filt, f.m, gs, ode, tc, tes)



## D

# Load in the filtered R correlation object using the "fastSave" library
load.lbzip2(file = "atac/rObject_cor.mat.filt.RDataFS", n.cores = 16)

# Set the random seed to apply to mclapply as well! 
RNGkind("L'Ecuyer-CMRG")

skip.streams <- function(n) {
  x <- .Random.seed
  for (i in seq_len(n))
    x <- nextRNGStream(x)
  assign('.Random.seed', x, pos=.GlobalEnv)
}

# Get a subset of all the relevant ode-associated peak-gene linkages
filt <- arrange(subset(cor.mat.filt, (filt.abc == T | filt.s == T) & abs(dist) > 1000 & g1 %in% c("orr1e", "orr1d2")), desc(abs(cor.s)))

# Load in the cell type pathways you got from a cool CITE-seq paper (https://www.sciencedirect.com/science/article/pii/S2211124723003157?via%3Dihub) 
# but remove the neutrophil-derived groups since those aren't present in your cells
newGO <- gmtPathways("atac/Konturek-Ciesla - 2023 - HSPC marker genes.txt")
newGO <- newGO[!grepl("Neu", names(newGO))]

set.seed(18)
skip.streams(8)
r1e.filt.gsea <- mclapply(1:8, function(f) {
	a <- subset(filt, g2 == f & g1 == "orr1e"); a <- group_by(a, gene) %>% arrange(desc(abs(cor.s))) %>% slice_head(n=1) %>% as.data.frame() %>% arrange(desc(cor.s))
	b <- a$cor.s; names(b) <- a$gene
	return( fgseaMultilevel(pathways = newGO, stats = b, minSize=5) %>% as.data.frame() %>% arrange(desc(NES)) )
}, mc.cores = 8)

skip.streams(6)
r1d2.filt.gsea <- mclapply(1:6, function(f) {
	a <- subset(filt, g2 == f & g1 == "orr1d2"); a <- group_by(a, gene) %>% arrange(desc(abs(cor.s))) %>% slice_head(n=1) %>% as.data.frame() %>% arrange(desc(cor.s))
	a[is.na(a$cor.s),"cor.s"] <- 0; b <- a$cor.s; names(b) <- a$gene
	return( fgseaMultilevel(pathways = newGO, stats = b, minSize=5) %>% as.data.frame() %>% arrange(desc(NES)) )
}, mc.cores = 6)


filt.abc <- subset(filt, filt.abc == T)

# Generate the new version of the GSEA plots
tmp <- data.frame()
for (f in 1:length(r1e.filt.gsea)) {
	if ( nrow( r1e.filt.gsea[[f]] ) != 0 ) {
		temp <- r1e.filt.gsea[[f]]; temp$g <- f
		tmp <- rbind(tmp, temp)
	}
}
tmp$g <- factor(tmp$g, levels = (8:1)); rm(temp)

# Have size be tmp-value and color be NES
g3 <- ggplot(subset(tmp, padj <= 0.05), aes(x = pathway, y = g, size = -log10(padj), color = NES)) + geom_point() + theme_minimal_grid(line_size = 1/.75/.pt, color = "#EBEBEB") + 
	scale_y_discrete(name = "ORR1E ATAC-seq clusters", drop = F, labels = rev(c("Stem\n& Prog", "B", "T", "NK", "ILC3", "MF", "DC", "Bkgd"))) +
	scale_color_gradient2(low = "blue", mid = "white", high = "dark orange") + scale_x_discrete(name = "", guide = guide_axis(angle = 90)) + scale_size(name = "-log10\nadj pval", range = c(1,5), breaks = seq(2,6,2)) +
	guides(fill = "none", color = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(0.1, "in"), barheight = unit(.65, "in"))) + theme( text = element_text(size = 8, color = "black"), text = element_text(size = 8, color = "black"), axis.text = element_text(size = 8, color = "black"), axis.ticks.length = unit(.04, "in"), axis.ticks = element_line(linewidth = 1/.75/.pt, color = "#EBEBEB") )

# Remove the objects you no longer need
rm(filt, tmp, f, newGO, r1e.filt.gsea, r1d2.filt.gsea, skip.streams, cor.mat.filt)


## E & F

# Load in the filtered R correlation object using the "fastSave" library
load.lbzip2(file = "atac/rObject_cor.mat.filt.RDataFS", n.cores = 16)

# Load in the list of ODEs and ODE-derived peaks
ode <- read.table("atac/rObject_ode.txt", sep = "\t", header = T); ode$cpm.c <- factor(ode$cpm.c, levels = 1:8); ode$peak <- gsub(".*\\|", "", ode$pid)
rownames(ode) <- ode$id

# Subset for the top ODE candidate interactions
top <- arrange(subset(cor.mat.filt, filt.abc == T & filt.s == T & abs(dist) > 1000 & g1 %in% c("orr1e", "orr1d2")), desc(abs(cor.s)))

# Load in the values to make all the pretty colors
lut.a <- read.table("atac/immgen_heatmap_lut.txt", sep = "\t", header = T, comment.char = "")
p.col <- lut.a$col; names(p.col) <- rownames(lut.a)

### NOTE: Reminder that objects such as "rObject_l2.gex.qnorm.RDataFS" are created using the correlation.r script

# Load in the quantile normalized ATAC-seq values and RNA-seq values
load.lbzip2(file = "atac/rObject_l2.atac.cpm.qnorm.RDataFS", n.cores = 16)
load.lbzip2(file = "rna/rObject_l2.gex.qnorm.RDataFS", n.cores = 16)

# Make sure all the rows in the RNA-seq object are present in the ATAC-seq object
l2.atac.cpm.qnorm <- l2.atac.cpm.qnorm[ , colnames(l2.gex.qnorm) ]; rownames(l2.atac.cpm.qnorm) <- gsub(".*\\|", "", rownames(l2.atac.cpm.qnorm))

# Create a lookup table to convert cluster number to name
r1e <- c("Stem & Prog", "B", "T", "NK", "ILC3", "MF", "DC", "Bkgd")

s <- subset(top, gene == "Zc3h12a")
a <- data.frame("ATAC" = as.numeric(l2.atac.cpm.qnorm[ s[1,peak], ]), "RNA" = as.numeric(l2.gex.qnorm[ s[1,gene], ]), "cell" = colnames(l2.atac.cpm.qnorm))
g4 <- ggplot(a, aes(x = ATAC, y = RNA)) + geom_point(aes(color = cell), stroke = NA, size = 6/.pt, alpha = 0.75) + ggtitle(paste0(r1e[ as.numeric(s[1,g2]) ], " Cluster | rho = ", round(s[,cor.s], 2), " | ", abs(ceiling(s[,dist] / 1000)), "kb ", ifelse( s[1,dist] <= 0, "down", "up"), "stream")) +
 	scale_color_manual(values = p.col) + guides(color = "none", alpha = "none") + labs(x = paste0("log2(CPM)\n", subset(ode, peak == s[1,peak] & cpm.c == s[1,g2])$id[1]), y = paste0("log2(TPM) ", gsub(".*\\|", "", s[1,gene]))) + theme_classic() +
	coord_cartesian(ylim = c( max(min(a$RNA) - 0.1, 0), max(a$RNA) + 0.1)) + theme(plot.title = element_text(size = 8, color = "black", hjust = 0.5), axis.text = element_text(size = 8, color = "black"), axis.title = element_text(size = 8, color = "black"), axis.ticks.length = unit(.04, "in"), axis.ticks = element_line(linewidth = 1/.75/.pt, color = "black"))

s <- subset(top, gene == "Prkcq")
a <- data.frame("ATAC" = as.numeric(l2.atac.cpm.qnorm[ s[1,peak], ]), "RNA" = as.numeric(l2.gex.qnorm[ s[1,gene], ]), "cell" = colnames(l2.atac.cpm.qnorm))
g5 <- ggplot(a, aes(x = ATAC, y = RNA)) + geom_point(aes(color = cell), stroke = NA, size = 6/.pt, alpha = 0.75) + ggtitle(paste0(r1e[ as.numeric(s[1,g2]) ], " Cluster | rho = ", round(s[,cor.s], 2), " | ", abs(ceiling(s[,dist] / 1000)), "kb ", ifelse( s[1,dist] <= 0, "down", "up"), "stream")) +
 	scale_color_manual(values = p.col) + guides(color = "none", alpha = "none") + labs(x = paste0("log2(CPM)\n", subset(ode, peak == s[1,peak] & cpm.c == s[1,g2])$id[1]), y = paste0("log2(TPM) ", gsub(".*\\|", "", s[1,gene]))) + theme_classic() +
	coord_cartesian(ylim = c( max(min(a$RNA) - 0.1, 0), max(a$RNA) + 0.1)) + theme(plot.title = element_text(size = 8, color = "black", hjust = 0.5), axis.text = element_text(size = 8, color = "black"), axis.title = element_text(size = 8, color = "black"), axis.ticks.length = unit(.04, "in"), axis.ticks = element_line(linewidth = 1/.75/.pt, color = "black"))


pdf("figures/fig4_cowplot.pdf", height = 5.5, width = 7)
plot_grid( plot_grid( NULL, g1, g2, labels = c("A", "B", "C"), align = "h", nrow = 1, rel_widths = c(.25, .394, .34 ), label_size = 12 ), plot_grid( g3, g4, g5, labels = c("D", "E", "F"), nrow = 1, rel_widths = c(.369, .308, .323 ), label_size = 12 ), rel_heights = c(.57, .43), nrow = 2)
dev.off()



################
### Figure 5 ###
################

# Make sure to add "lbzip2" to your environment prior to starting R
export PATH=/home/jdc397/lbzip2-2.5/bin:$PATH

# Load in the packages you'll need throughout this section
suppressPackageStartupMessages({ 
	library(fastSave); library(orthogene); library(data.table); setDTthreads(8); library(doParallel); library(dplyr); library(fgsea)
	library(cowplot); library(ggplot2); library(ggthemes); library(scatterpie)
}) 

## A & B

# Load in the R correlation object using the "fastSave" library
load.lbzip2(file = "atac/rObject_cor.mat.filt.RDataFS", n.cores = 16)

# Generate a look-up-table that contains the samples you'll use for each different cell type in human and mouse
lut.hm <- data.frame("mouse" = list(), "human" = list())
lut.hm["nCD8","mouse"][[1]] <- list(c("T.8.Nve.Sp", "T.8.Th")); lut.hm["nCD8","human"][[1]] <- list(c("CD8_monaco", "CD8_schmiedel", "CD8_hpa"))
lut.hm["B","mouse"][[1]] <- list(c("B.Sp")); lut.hm["B","human"][[1]] <- list(c("B_schmiedel", "B_monaco", "B_hpa"))
lut.hm["nCD4","mouse"][[1]] <- list(c("T.4.Nve.Sp", "T.4.Th")); lut.hm["nCD4","human"][[1]] <- list(c("CD4_monaco", "CD4_schmiedel", "CD4_hpa"))
lut.hm["NK","mouse"][[1]] <- list(c("NK.27mi11bpl.Sp", "NK.27pl11bmi.Sp", "NK.27pl11bpl.Sp")); lut.hm["NK","human"][[1]] <- list(c("NK_monaco", "NK_schmiedel", "NK_hpa"))
lut.hm["aCD8","mouse"][[1]] <- list(c("T8.TN.P14.Sp", "T8.MP.LCMV.d7.Sp", "T8.TE.LCMV.d7.Sp")); lut.hm["aCD8","human"][[1]] <- list(c("Act.CD8_schmiedel"))
lut.hm["aCD4","mouse"][[1]] <- list(c("T.4.aCD3plCD40.18h.Sp")); lut.hm["aCD4","human"][[1]] <- list(c("Act.CD4_schmiedel"))
lut.hm["gdT","mouse"][[1]] <- list(c("Tgd.Sp", "Tgd.g1.1pld1.LN", "Tgd.g2pld17.LN", "Tgd.g2pld1.LN")); lut.hm["gdT","human"][[1]] <- list(c("gdT_monaco", "gdT_hpa"))
lut.hm["mDC","mouse"][[1]] <- list(c("DC.4pl.Sp")); lut.hm["mDC","human"][[1]] <- list(c("mDC_monaco", "mDC_hpa"))
lut.hm["GN","mouse"][[1]] <- list(c("GN.Sp", "GN.Thio.PC")); lut.hm["GN","human"][[1]] <- list(c("basophil_monaco", "basophil_hpa", "eosinophil_hpa"))
lut.hm["pDC","mouse"][[1]] <- list(c("DC.pDC.Sp")); lut.hm["pDC","human"][[1]] <- list(c("pDC_monaco", "pDC_hpa"))
lut.hm["C.Mo","mouse"][[1]] <- list(c("Mo.6CplIImi.Bl")); lut.hm["C.Mo","human"][[1]] <- list(c("C.Mo_monaco", "C.Mo_hpa", "C.Mo_schmiedel"))
lut.hm["NC.Mo","mouse"][[1]] <- list(c("Mo.6CmiIImi.Bl")); lut.hm["NC.Mo","human"][[1]] <- list(c("NC.Mo_monaco", "NC.Mo_hpa", "NC.Mo_schmiedel"))

# Now save the above lut so you can load it later as needed
save.lbzip2(lut.hm, file = "orthology/rObject_lut.hm.RDataFS", n.cores = 2)

# Load in the RNA-seq data from ImmGen that is quantile normalized log2 CPM
load.lbzip2("rna/rObject_l2.gex.qnorm.RDataFS", n.cores = 16)

# Scale the above object so you can more reasonably compare the values between mouse and human
l2.gex.qnorm <- scale(l2.gex.qnorm)

# Load in the RNA-seq data from the Human Protein Atlas that is quantile normalized log2 CPM
load.lbzip2("human/rObject_l2.hpa.qnorm.RDataFS", n.cores = 16)

# Scale the above object so you can more reasonably compare the values between mouse and human
l2.hpa.qnorm <- scale(l2.hpa.qnorm)

# Use the power of orthogene to figure out all of the orthologous/syntenic genes between mouse and human
orth.mh <- convert_orthologs(rownames(l2.gex.qnorm), input_species = "mouse", output_species = "human", non121_strategy = "drop_both_species", gene_output = "columns", method = "homologene", drop_nonorths = T)
orth.mg <- convert_orthologs(rownames(l2.gex.qnorm), input_species = "mouse", output_species = "human", non121_strategy = "drop_both_species", gene_output = "columns", method = "gprofiler", drop_nonorths = T)
orth.mb <- convert_orthologs(rownames(l2.gex.qnorm), input_species = "mouse", output_species = "human", non121_strategy = "drop_both_species", gene_output = "columns", method = "babelgene", drop_nonorths = T)

# Merge them all together to get a list of the orthogonal genes!
orth.m <- unique( rbind(orth.mh, orth.mg, orth.mb) ) 

# Use the reciprocal approach to make sure you get everything
orth.hh <- convert_orthologs( unique(gsub(".*\\|", "", rownames(l2.hpa.qnorm))), input_species = "human", output_species = "mouse", non121_strategy = "drop_both_species", gene_output = "columns", method = "homologene", drop_nonorths = T)
orth.hg <- convert_orthologs( unique(gsub(".*\\|", "", rownames(l2.hpa.qnorm))), input_species = "human", output_species = "mouse", non121_strategy = "drop_both_species", gene_output = "columns", method = "gprofiler", drop_nonorths = T)
orth.hb <- convert_orthologs( unique(gsub(".*\\|", "", rownames(l2.hpa.qnorm))), input_species = "human", output_species = "mouse", non121_strategy = "drop_both_species", gene_output = "columns", method = "babelgene", drop_nonorths = T)

# Merge them all together to get a list of the orthogonal genes!
orth.h <- unique( rbind(orth.hh, orth.hg, orth.hb) ) 

# Combine them all together so that you have the full list. No unique rownames unfortunately...
orth <- unique(data.frame( "input_gene" = c(orth.m$input_gene, orth.h$ortholog_gene), "ortholog_gene" = c(orth.m$ortholog_gene, orth.h$input_gene) ))
rm(orth.m, orth.h)

# Use the power of orthogene to figure out all of the complicated orthologous/syntenic gene relationships between mouse and human
orth.mh <- convert_orthologs(rownames(l2.gex.qnorm), input_species = "mouse", output_species = "human", non121_strategy = "keep_both_species", gene_output = "columns", method = "homologene", drop_nonorths = T)
orth.mg <- convert_orthologs(rownames(l2.gex.qnorm), input_species = "mouse", output_species = "human", non121_strategy = "keep_both_species", gene_output = "columns", method = "gprofiler", drop_nonorths = T)
orth.mb <- convert_orthologs(rownames(l2.gex.qnorm), input_species = "mouse", output_species = "human", non121_strategy = "keep_both_species", gene_output = "columns", method = "babelgene", drop_nonorths = T)

# Merge them all together to get a list of the orthogonal (not necessarily 1:1 here) genes!
all.m <- unique( rbind(orth.mh, orth.mg, orth.mb) ) 

# Use the reciprocal approach to make sure you get everything
orth.hh <- convert_orthologs( unique(gsub(".*\\|", "", rownames(l2.hpa.qnorm))), input_species = "human", output_species = "mouse", non121_strategy = "keep_both_species", gene_output = "columns", method = "homologene", drop_nonorths = T)
orth.hg <- convert_orthologs( unique(gsub(".*\\|", "", rownames(l2.hpa.qnorm))), input_species = "human", output_species = "mouse", non121_strategy = "keep_both_species", gene_output = "columns", method = "gprofiler", drop_nonorths = T)
orth.hb <- convert_orthologs( unique(gsub(".*\\|", "", rownames(l2.hpa.qnorm))), input_species = "human", output_species = "mouse", non121_strategy = "keep_both_species", gene_output = "columns", method = "babelgene", drop_nonorths = T)

# Merge them all together to get a full list of the orthogonal genes (not necessarily 1:1 here, can be quite complicated relationships!)
all.h <- unique( rbind(orth.hh, orth.hg, orth.hb) ) 

# Combine them all together so that you have the full list. No unique rownames unfortunately...
all <- unique(data.frame( "input_gene" = c(all.m$input_gene, all.h$ortholog_gene), "ortholog_gene" = c(all.m$ortholog_gene, all.h$input_gene) ))
rm(orth.mh, orth.mg, orth.mb, all.m, orth.hh, orth.hg, orth.hb, all.h)

# Now make a list of all the possible combos of orthologous, mouse-specific, or human specific
gene.scores <- data.frame("mouse.name" = orth$input_gene, "mouse.score" = 0, "human.name" = orth$ortholog_gene, "human.score" = 0, "type" = "orth") 

# "Mouse-specific"
tmp <- rownames(l2.gex.qnorm)[ which( !rownames(l2.gex.qnorm) %in% gene.scores$mouse.name ) ] 
tmp2 <- data.frame("mouse.name" = tmp, "mouse.score" = 0, "human.name" = NA, "human.score" = 0)
tmp2$type <- ifelse( tmp2$mouse.name %in% all$input_gene, "comp.mouse", "mouse" )
gene.scores <- rbind(gene.scores, tmp2 )

# "Human-specific"
tmp <- unique(gsub(".*\\|", "", rownames(l2.hpa.qnorm))[ which( !gsub(".*\\|", "", rownames(l2.hpa.qnorm)) %in% gene.scores$human.name ) ])
tmp2 <- data.frame("mouse.name" = NA, "mouse.score" = 0, "human.name" = tmp, "human.score" = 0)
tmp2$type <- ifelse( tmp2$human.name %in% all$ortholog_gene, "comp.human", "human" )
gene.scores <- rbind(gene.scores, tmp2)

# Check to see which of the provided genes are actually in the gene expression matricies
hs <- unique( gsub(".*\\|", "", rownames(l2.hpa.qnorm)) )
for (f in 1:nrow(gene.scores)) {
	if( gene.scores[f,"type"] == "orth" & !gene.scores[f,"mouse.name"] %in% rownames(l2.gex.qnorm) ) { gene.scores[f,"type"] <- "orth.human" }
	if( gene.scores[f,"type"] == "orth" & !gene.scores[f,"human.name"] %in% hs ) { gene.scores[f,"type"] <- "orth.mouse" }
}

# Make a small version of "cor.mat" so that it can be subsetted quicker and without much memory usage
smol.dt <- as.data.table(subset(cor.mat.filt, g1 %in% c("orr1e", "orr1d2"))); setkey(smol.dt, gene)

# Before you generate the GSEA data, get the underlying information for the mouse vs human delta TPM going! 
tmp <- data.frame(matrix(0, nrow = 0, ncol = ncol(gene.scores) + 6)); colnames(tmp) <- c(colnames(gene.scores), c("r1e.clust", "r1d2.clust", "r1e.con.clust", "r1e.spear", "r1d2.con.clust", "r1d2.spear"))
p <- bind_rows( mclapply(1:nrow(gene.scores), function(i) {
	tmp[1,1:5] <- gene.scores[i,1:5]

	# Find all the peaks linked with the current gene and classify the peak based on whether there's an ODE there or not
	s <- smol.dt[J(gene.scores[i,"mouse.name"])][filt.abc == T | filt.s == T] 

	# If there are multiple 1es and there are different clusters within them
	if ( sum(s$g1 == "orr1e") > 1 & length(unique(subset(s, g1 == "orr1e")$g2)) > 1 ) {
		tmp[1,"r1e.clust"] <- "mult"
		tmp[1,"r1e.con.clust"] <- arrange(subset(s, g1 == "orr1e"), desc(abs(cor.s)))$g2[1]
		tmp[1,"r1e.spear"] <- arrange(subset(s, g1 == "orr1e"), desc(abs(cor.s)))$cor.s[1]

	# If there are multiple 1es that are in the same cluster or just filt single 1e
	} else if ( ( sum(s$g1 == "orr1e") > 1 & length(unique(subset(s, g1 == "orr1e")$g2)) == 1 ) | sum(s$g1 == "orr1e") == 1 ) {
		tmp[1,"r1e.clust"] <- tmp[1,"r1e.con.clust"] <- subset(s, g1 == "orr1e")$g2[1]	
		tmp[1,"r1e.spear"] <- arrange(subset(s, g1 == "orr1e"), desc(abs(cor.s)))$cor.s[1]	

	# Otherwise, just have it be other to incorporate both TE and non-TE-derived interactions
	} else {
		tmp[1,"r1e.clust"] <- tmp[1,"r1e.con.clust"] <- "other"; tmp[1,"r1e.spear"] <- NA
	}

	# Now repeat the above with no comments, but for ORR1D2!
	if ( sum(s$g1 == "orr1d2") > 1 & length(unique(subset(s, g1 == "orr1d2")$g2)) > 1 ) {
		tmp[1,"r1d2.clust"] <- "mult"
		tmp[1,"r1d2.con.clust"] <- arrange(subset(s, g1 == "orr1d2"), desc(abs(cor.s)))$g2[1]; tmp[1,"r1d2.spear"] <- arrange(subset(s, g1 == "orr1d2"), desc(abs(cor.s)))$cor.s[1]
	} else if ( ( sum(s$g1 == "orr1d2") > 1 & length(unique(subset(s, g1 == "orr1d2")$g2)) == 1 ) | sum(s$g1 == "orr1d2") == 1 ) {
		tmp[1,"r1d2.clust"] <- tmp[1,"r1d2.con.clust"] <- subset(s, g1 == "orr1d2")$g2[1]; tmp[1,"r1d2.spear"] <- arrange(subset(s, g1 == "orr1d2"), desc(abs(cor.s)))$cor.s[1]	
	} else {
		tmp[1,"r1d2.clust"] <- tmp[1,"r1d2.con.clust"] <- "other"; tmp[1,"r1d2.spear"] <- NA
	}

	# Finally, output the current sample
	return(tmp)

}, mc.cores = 16) )

# Make sure the GSEA are reproducible 
RNGkind("L'Ecuyer-CMRG")
skip.streams <- function(n) {
  x <- .Random.seed
  for (i in seq_len(n))
	x <- nextRNGStream(x)
  assign('.Random.seed', x, pos=.GlobalEnv)
}

# It's time to iterate over each of the various cell types in "lut.hm" to get the GSEA results! Yeah
res <- bind_rows( mclapply(rownames(lut.hm), function(g) {
	val <- p
	val$mouse.score <- unlist( mclapply(1:nrow(val), function(i) {
		# Only operate on gene pairs that are present in mouse, ie *don't* have "human" in the type column
		if ( !grepl("human$", val[i,"type"]) ) {

			# Find which row in the gene expression object contain the current gene (since there's only one here for sure)
			v <- which( rownames(l2.gex.qnorm) == val[i,"mouse.name"] ); w <- unlist(lut.hm[g,"mouse"])
			
			# If there are multiple samples, average them. Otherwise take the single value
			if ( length( w ) > 1 ) {
				return( mean( l2.gex.qnorm[ v, w ] ) )
			} else {
				return( l2.gex.qnorm[ v, w ] )
			}

		} else {
			return( 0 )
		}
	}, mc.cores = 16) )

	# Since there are some entries with multiple matching names, you have to potentially account for them...
	val$human.score <- unlist( mclapply(1:nrow(val), function(i) {
		# Only operate on gene pairs that are present in human, ie *don't* have "mouse" in the type column
		if ( !grepl("mouse$", val[i,"type"]) ) {

			# Find which row(s) in the gene expression object contain the current gene
			v <- which( gsub(".*\\|", "", rownames(l2.hpa.qnorm)) == val[i,"human.name"] ); w <- unlist(lut.hm[g,"human"])

			# If there are multiple samples, average them. Otherwise take the single value as needed
			if ( length(v) > 1 ) {
				if ( length( w ) > 1 ) {
					return( mean( rowMeans( l2.hpa.qnorm[ v, w ] ) ) )
				} else {
					return( mean( l2.hpa.qnorm[ v, w ] ) )
				}
			} else if ( length(v) == 1 ) {
				if ( length( w ) > 1 ) {
					return( mean( l2.hpa.qnorm[ v, w ] ) )
				} else {
					return( l2.hpa.qnorm[ v, w ] )
				}
			}
		} else {
			return( 0 )
		}
	}, mc.cores = 16) )

	# Generate the gene set-containing object 
	pws <- list()

	# You decided to also include all possible "mouse-specific" genes, including those with not a 1:1 ortholog or those without human
	# expression data
	s <- unique(arrange(subset(val, grepl("mouse$", type)), desc(mouse.score))$mouse.name)
	pws[["combMouseUp"]] <- unique(arrange(subset(val, grepl("mouse$", type) & mouse.score > 0), desc(mouse.score))$mouse.name)
	pws[["combMouseDown"]] <- unique(arrange(subset(val, grepl("mouse$", type) & mouse.score < 0), desc(mouse.score))$mouse.name)

	# Include everything that is orthologous and up in mouse compared to human, while including the opposite for depletion
	val$delta.mvh <- val$mouse.score - val$human.score
	pws[["orthMouseUp_pos"]] <- unique(arrange(subset(val, type == "orth" & delta.mvh > 0 & mouse.score > 0), desc(delta.mvh))$mouse.name)

	# Now look at those that are down in mouse and up in human, potentially due to a repressive effect in mouse (less clear but eh)
	tmp <- unique(arrange(subset(val, type == "orth" & delta.mvh < 0 & mouse.score < 0), desc(delta.mvh))$mouse.name); s <- length(tmp)
	pws[["orthMouseDown_neg"]] <- tmp

	# Finally, save the above objects so you don't have to generate them again
	write.table(val, file = paste0("orthology/mouse_vs_human_orthology_", g, ".txt"), quote = F, sep = "\t", row.names = F)
	save.lbzip2(pws, file = paste0("orthology/object_orth_gsea_pathways_", g, ".RDataFS"), n.cores = 16)

	# Now it's time for the GSEA analyses! Yeah
	set.seed(18)
	skip.streams(8)
	r1e <- mclapply(1:8, function(f) {
		a <- subset(smol.dt, g2 == f & g1 == "orr1e" & !is.na(cor.s) & (filt.abc == T | filt.s == T)); a <- group_by(a, gene) %>% arrange(desc(abs(cor.s))) %>% slice_head(n=1) %>% as.data.frame() %>% arrange(desc(cor.s))
		b <- a$cor.s; names(b) <- a$gene
		return( fgseaMultilevel(pathways = pws, stats = b, minSize=3) %>% as.data.frame() %>% mutate("cell.type" = g, "ode" = "orr1e", "cluster" = f) %>% arrange(desc(NES)) )
	}, mc.cores = 8)

	skip.streams(6)
	r1d2 <- mclapply(1:6, function(f) {
		a <- subset(smol.dt, g2 == f & g1 == "orr1d2" & !is.na(cor.s) & (filt.abc == T | filt.s == T)); a <- group_by(a, gene) %>% arrange(desc(abs(cor.s))) %>% slice_head(n=1) %>% as.data.frame() %>% arrange(desc(cor.s))
		b <- a$cor.s; names(b) <- a$gene
		return( fgseaMultilevel(pathways = pws, stats = b, minSize=3) %>% as.data.frame() %>% mutate("cell.type" = g, "ode" = "orr1d2", "cluster" = f) %>% arrange(desc(NES)) )
	}, mc.cores = 6)

	return( rbind( bind_rows(r1e), bind_rows(r1d2) ) )

}, mc.cores = 5) )

# With everything together in a single object, save it for the GSEA plots
save.lbzip2(res, file = "orthology/rObject_mouseVsHumanOrthologyGsea_res.RDataFS", n.cores = 16)

# Remove all the objects you no longer need
rm(gene.scores,l2.gex.qnorm,l2.hpa.qnorm,smol.dt,p,orth,all)

# Load in all of the assigned mouse vs human delta TPMs across just NK and GN for ORR1E, with the rest plotted in the supplement
l <- list.files(pattern = "orthology/mouse_vs_human_orthology_*"); l <- l[grepl("NK|GN", l)]; tmp <- data.frame()
for (f in l) {
	a <- read.table(f, sep = "\t", header = T); a$samp <- gsub(".txt", "", gsub(".*_", "", f))
	tmp <- rbind(tmp, a)
}
mvh <- tmp; rm(l, tmp, a)

# Load the look-up-table that contains the samples you'll use for each different cell type in mouse and human
# load.lbzip2(file = "orthology/rObject_lut.hm.RDataFS", n.cores = 2)

# Get that color palette ready (with an extra value for all other genes!) 
col.c <- c(tableau_color_pal(type = "ordered-diverging", palette = "Sunset-Sunrise Diverging")(7), "gray69", "gray51") 

# Generate the distribution of positively correlated ODEs for each ATAC-seq cluster for each cell type to identify whether the target genes are significantly 
# more expressed in mouse than in human. Include all putatively orthologous genes! Also, don't include ORR1E GN and NK in the supplemental figure
# and instead have them by themselves with the numbers
h <- list()
for ( f in sort(unique(mvh$samp), decreasing = T) ) {

	# Subset the object containing all the gene expression score information for the current cell type
	tmp <- subset(mvh, samp == f & type == "orth" & ((r1e.con.clust != "other" & r1e.spear >= 0) | r1e.con.clust == "other")); tmp$lab <- factor(tmp$r1e.con.clust, levels = c("other",8:1))

	# Generate the object and calculate stats for plotting 
	meta <- data.frame("lab" = factor(c(1:8,"other"), levels = c("other",8:1)), 
		"p.val" = c(p.adjust(sapply(1:8, function(i){wilcox.test(subset(tmp, r1e.con.clust == i)$delta.mvh, subset(tmp, r1e.con.clust == "other")$delta.mvh, alternative = "two.sided")$p.value}), method = "BH"), 1), 
		"delta.mvh" = rep(-1.9,9), 
		"num" = as.numeric(table(tmp$r1e.con.clust)),
		"dir" = c(ifelse( sapply(1:8, function(i){ median(subset(tmp, r1e.con.clust == i)$delta.mvh) }) <= median(subset(tmp, r1e.con.clust == "other")$delta.mvh), "less", "greater" ), NA) )
	meta$sig <- with(meta, ifelse(p.val < 0.0005, "***", ifelse(p.val < 0.005, "**", ifelse(p.val < 0.05, "*", "") ) ) )
	meta[meta$lab == "other","num"] <- paste0("~", round(meta[meta$lab == "other","num"]/1000,0), "k")

	# Create the base plot and determine whether to add y-axis labels (yes for "NK", no for "GN")
	base <- ggplot(tmp, aes(x = lab, y = delta.mvh, color = lab, fill = lab)) + theme_classic() + geom_hline(yintercept = median(subset(tmp, lab == "other")$delta.mvh), linetype = "31", linewidth = 1/.75/.pt) +
		geom_boxplot(width = 0.5, outlier.color = NA, fill = "white", show.legend = F, linewidth = 1/.75/.pt, fatten = 1) + scale_color_manual(values = col.c[9:1]) + 
		new_scale_color() + geom_text(data = meta, aes(x = lab, y = delta.mvh, label = sig, color = dir), size = 8/.pt, fontface = "bold", angle = 270, show.legend = F) +
		scale_color_manual(values = c("less" = "black", "greater" = "dark orange")) + geom_text(data = meta, aes(x = lab, y = delta.mvh - 0.3, label = num), color = "black", size = 8/.pt, show.legend = F) + 
		coord_flip(ylim = c(-2.5,2.5)) + labs( x = "", y = expression(Delta~"TPM (Mouse - Human)") ) +
		ggtitle(paste0("ORR1E: ", f, " expression")) + theme(plot.title = element_text(size = 8, color = "black", hjust = 0.5), axis.text = element_text(size = 8, color = "black"), axis.title = element_text(size = 8), axis.ticks.length = unit(.04, "in"), axis.ticks = element_line(linewidth = 1/.75/.pt, color = "black"), axis.line = element_line(linewidth = 1/.75/.pt, color = "black"))

	# Add the cell type plot to the list and selectively include y-axis labels for "NK" only
	if ( f == "NK") {
		h[[ length(h) + 1 ]] <- base + scale_x_discrete(labels = rev(c("Stem\n& Prog", "B", "T", "NK", "ILC3", "MF", "DC", "Bkgd", "All Other\nGenes")))
	} else {
		h[[ length(h) + 1 ]] <- base + scale_x_discrete( labels = rep("", length(levels(meta$lab))) )
	}
}


## C 

### NOTE: This figure was generated entirely in Illustrator. Yay?


## D & F

load.lbzip2(file = "orthology/rObject_mouseVsHumanOrthologyGsea_res.RDataFS", n.cores = 16)
load.lbzip2(file = "orthology/rObject_lut.hm.RDataFS", n.cores = 2)

# Generate the new version of the GSEA plots for mouse-specific up/down and mouse orthologous up/down
tmp <- subset(res, pathway == "combMouseUp" & ode == "orr1e")
tmp$g <- factor(tmp$cluster, levels = (8:1)); tmp$cell.type <- factor(tmp$cell.type, levels = c("nCD4", "aCD4", "nCD8", "aCD8", "gdT", "NK", "B", "", "C.Mo", "NC.Mo", "GN", "mDC", "pDC"))

# Generate the plot for mouse-specific genes for ORR1E
g1.m <- ggplot(subset(tmp, padj > 0.05 & padj <= 0.25), aes(x = cell.type, y = g, size = -log10(padj), color = NES)) + geom_point(pch = 1, stroke = 1/.pt/.5) + geom_point(data = subset(tmp, padj <= 0.05), aes(x = cell.type, y = g, size = -log10(padj), color = NES), stroke = 1/.pt/.5) + 
	theme_minimal_grid(line_size = 0.5/.75/.pt, color = "gray69") + scale_y_discrete(drop = F, labels = rev(c("Stem\n& Prog", "B", "T", "NK", "ILC3", "MF", "DC", "Bkgd"))) +
	scale_color_gradient2(low = "blue", mid = "white", high = "dark orange") + scale_x_discrete(guide = guide_axis(angle = 90), breaks = rownames(lut.hm), drop = F) + scale_size(range = c(1,3), breaks = c(1,2), name = "-log10\nadj pval") +
	ylab("") + xlab("Cell type expression data") + ggtitle("Mouse-specific") + theme( legend.text = element_text(size = 8, color = "black"), legend.title = element_text(size = 8, color = "black"), plot.title = element_text(size = 8, color = "black", hjust = 0.5, face = "plain"), axis.text = element_text(size = 8, color = "black"), axis.title = element_text(size = 8, color = "black"), axis.ticks.length = unit(.04, "in"), axis.ticks = element_line(linewidth = .5/.75/.pt, color = "gray69") ) +
	guides(fill = "none", size = guide_legend(override.aes = list(stroke = 1/.pt/.5), order = 1), color = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(0.1, "in"), barheight = unit(.65, "in"), order = 2))

# Now for the orthologous genes for ORR1E
tmp <- subset(res, pathway == "orthMouseUp_pos" & ode == "orr1e")
tmp$g <- factor(tmp$cluster, levels = (8:1)); tmp$cell.type <- factor(tmp$cell.type, levels = c("nCD4", "aCD4", "nCD8", "aCD8", "gdT", "NK", "B", "", "C.Mo", "NC.Mo", "GN", "mDC", "pDC"))

g1.o <- ggplot(subset(tmp, padj > 0.05 & padj <= 0.25), aes(x = cell.type, y = g, size = -log10(padj), color = NES)) + geom_point(pch = 1, stroke = 1/.pt/.5) + geom_point(data = subset(tmp, padj <= 0.05), aes(x = cell.type, y = g, size = -log10(padj), color = NES), stroke = 1/.pt/.5) + 
	theme_minimal_grid(line_size = 0.5/.75/.pt, color = "gray69") + scale_y_discrete(drop = F, labels = rev(c("Stem\n& Prog", "B", "T", "NK", "ILC3", "MF", "DC", "Bkgd"))) +
        scale_color_gradient2(low = "blue", mid = "white", high = "dark orange") + scale_x_discrete(guide = guide_axis(angle = 90), breaks = rownames(lut.hm), drop = F) + scale_size(range = c(1,3), name = "-log10\nadj pval") +
        ylab("") + xlab("Cell type expression data") + ggtitle("Orthologous") + theme( legend.text = element_text(size = 8, color = "black"), legend.title = element_text(size = 8, color = "black"), plot.title = element_text(size = 8, color = "black", hjust = 0.5, face = "plain"), axis.text = element_text(size = 8, color = "black"), axis.title = element_text(size = 8, color = "black"), axis.ticks.length = unit(.04, "in"), axis.ticks = element_line(linewidth = .5/.75/.pt, color = "gray69") ) +
	guides(fill = "none", size = guide_legend(override.aes = list(stroke = 1/.pt/.5), order = 1), color = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(0.1, "in"), barheight = unit(.65, "in"), order = 2))


## E

# Load in the filtered R correlation object using the "fastSave" library
load.lbzip2(file = "atac/rObject_cor.mat.filt.RDataFS", n.cores = 16)

# Focus on ORR1E here, so assign all ORR1D2 links to the "te" category
tmp <- cor.mat.filt; tmp[ tmp$g1 == "orr1d2", "g1"] <- tmp[ tmp$g1 == "orr1d2", "g2"] <- "te"

# Arbitrarilly load in one of the previous objects that contains orthology information for gene targets (any cell type would do)
b <- fread(file = "orthology/mouse_vs_human_orthology_nCD8.txt"); setkey(b, mouse.name)

# Use the orthology information saved in the above file to add that to the filtered correlation file per gene
tmp$type <- unlist(mclapply(1:nrow(tmp), function(i){ arrange(b[J(tmp[i,"gene"])], desc(abs(delta.mvh)))[,type][1] }, mc.cores = 64))
tmp$type <- ifelse( grepl("orth", tmp$type) | tmp$type == "comp.mouse", "orth", tmp$type )

# Generate the table of *unique* ODE cluster / peak class and synteny status (not including parity of the correlation)
tmp.filt <- unique(subset(tmp, !is.na(type))[,c("gene","g2","type")])
a <- as.data.frame(with(tmp.filt, table(g2, type)))

# Add the genomic distribution of mouse-specific to orthologous genes to the above object for comparison (calculated from above "b" file)
sm <- sum(subset(b, !is.na(mouse.name))$type == "mouse"); so <- sum(subset(b, !is.na(mouse.name))$type != "mouse")
a <- rbind( a, data.frame("g2" = rep("genome", 2), "type" = c("mouse", "orth"), "Freq" = c(sm,so) ) )

# Create the "frac" column by either dividing by the total size of each "g2", or split by "g2" and "Var3" which is the parity of the correlation
a$frac.all <- a$Freq / c( rep( table( tmp.filt$g2), 2 ), rep(sum(sm,so), 2) ) * 100

# Figure out (via fisher's test) if there are fewer mouse-specific ODE-gene links (out of protein-coding genes!) than the genome proportion
b <- table( tmp.filt$g2 ); b[ length(b) + 1 ] <- sum(sm,so); names(b)[length(b)] <- "genome"
sig <- data.frame( "g2" = factor(c(1:8,"nonte","te"), levels = c("genome", "te", "nonte", 8:1)), "sig" = p.adjust( sapply( c(1:8,"nonte","te"), function(i){ a1 <- subset(a, g2 == i & type == "mouse")$Freq; fisher.test( matrix(c( a1, b[i] - a1, sm, so ), nrow = 2), alternative = "less" )$p.value }), method = "BH") )
sig$lab <- ifelse( sig$sig <= 0.05, "*", "" )

# Assign factor levels and generate the plot!
a$g2 <- factor(a$g2, levels = c("genome", "te", "nonte", 8:1)); a$type <- factor(a$type, levels = c("orth", "mouse"))
g4 <- ggplot(a, aes(x = frac.all, y = g2, fill = type)) + geom_vline(xintercept = subset(a, g2 == "genome" & type == "mouse")$frac.all, linetype = "31", linewidth = 1/0.75/.pt) + geom_col(width = 0.5) + scale_fill_manual(values = c("#000000", "#FF7518"), breaks = c("orth", "mouse"), labels = c("orth" = "Orthologous", "mouse" = "Mouse-Specific")) +
	scale_y_discrete(name = "", labels = NULL) + scale_x_continuous(name = "Filtered interactions (%)", labels = c("0","25","50","75","100"), expand = c(0,0)) + theme_classic() +
	coord_cartesian(clip = "off", xlim = c(0,100)) + annotate(geom = "text", x = -9, y = levels(a$g2), label = rev(c("Stem\n& Prog", "B", "T", "NK", "ILC3", "MF", "DC", "Bkgd", "nonTE", "TE", "Genome")), size = 8/.pt, hjust = "right") +
	annotate(geom = "text", x = -2, y = rev(levels(a$g2)), label = c(sig$lab, ""), fontface = "bold", size = 8/.pt, hjust = "right", vjust = .75) +
	theme( legend.position = "top", legend.key.size = unit(.5, 'cm'), legend.text = element_text(size = 8, color = "#000000"), legend.title = element_text(size = 8, color = "#000000"), axis.text = element_text(size = 8, color = "black"), axis.title.y = element_text(size = 8, color = "black", margin = margin(r = 35)), axis.title.x = element_text(size = 8, color = "#000000"), axis.ticks.length = unit(.04, "in"), axis.ticks.y = element_blank(), axis.ticks.x = element_line(linewidth = 1/.75/.pt, color = "black"), axis.line = element_line(linewidth = 1/.75/.pt, color = "black") ) +
	labs(fill = "")

# Now combine all the panels together into the main figure. It will require some Illustrator work to shmush A & B together and clean up D, E, & F
pdf("figures/fig5abdf.pdf", width = 7.5, height = 6)
plot_grid( plot_grid( plot_grid( plotlist = h, ncol = 2, nrow = 1, labels = c("A", "B"), align = "hv", label_size = 12 ), NULL, ncol = 2, rel_widths = c(1, 0.33), labels = c("", "C"), label_size = 12), plot_grid( g1.o, g4, g1.m, ncol = 3, nrow = 1, labels = c("D", "E", "F"), rel_widths = c(1, 0.65, 1) ), nrow = 2 )
dev.off()



################
### Figure 6 ###
################

# export PATH=/home/jdc397/lbzip2-2.5/bin:$PATH; module load R/4.1.3-r9; R --vanilla

# Load all the relevant libraries you need for this section
suppressPackageStartupMessages({ 
	library(Biostrings); library(data.table); setDTthreads(8); library(tidyr); library(doParallel); library(dplyr); library(fastSave); 
	library(ggthemes); library(ggplot2); library(PWMEnrich); library(TFBSTools); library(monaLisa); library(patchwork); library(scales)
	library(cowplot)
}) 

## A

# Load in the list of ODEs and ODE-derived peaks
ode <- read.table("atac/rObject_ode.txt", sep = "\t", header = T); ode$cpm.c <- factor(ode$cpm.c, levels = 1:8); ode$peak <- gsub(".*\\|", "", ode$pid)
rownames(ode) <- ode$id

# Load in the large file containing the intersection of 60-way phastCons scores and TEs
pc <- as.data.table(fread("orthology/ode.phastCons60way.int.txt.gz", sep = "\t")); setkey(pc, V4)

# Iterate over each element in the "pc" object 
it <- unique(pc$V4)
tmp <- as.data.frame(matrix(-1,nrow=1,ncol=max(pc$V3 - pc$V2 + 1) ))
pc.out <- bind_rows( mclapply(it, function(i) {
	val <- tmp; rownames(val) <- i; s <- pc[J(i)]; s$V8 <- ifelse( s$V8 < s$V2, s$V2, s$V8 ); s$V9 <- ifelse( s$V9 > s$V3, s$V3, s$V9 )

	# Have a case for no rows
	if ( nrow(s) == 0 ) { return(val); next }

	# For each row in the subset object, append the value (phastCons score) at "V10" to the correct location in "val"
	for ( f in 1:nrow(s) ) {

		# Get the start and end coordinates for the current interval
		ist <- s[f,"V8"] - s[f,"V2"] + 1; ied <- s[f,"V9"] - s[f,"V2"]

		# Use the start and end to get the right coordinates for the matrix
		val[ i, ist:ied ] <- as.numeric(s[ f, "V10" ])
	}

	# Finally, return the matrix so that it can be combined into one big happy matrix
	return( val )
}, mc.cores = 64) )
rm(tmp)

# Save the above so you don't have to regenerate it
save.lbzip2(pc.out, file = "orthology/rObject_pc.out.RDataFS", n.cores = 16)

# Turn all of the "-1" values into NAs
pc.out[pc.out == -1] <- NA

# Calculate the average conservation across each element and append it since that's generally helpful to have in one place
pc.avg <- data.frame("repName" = gsub(".*\\|", "", rownames(pc.out)), "avg" = rowMeans(pc.out, na.rm = T), "num" = rowSums(!is.na(pc.out)))

# Get the clusters for each ODE
pc.avg$cpm.c <- NA
for (f in 1:nrow(pc.avg)) {
	if ( rownames(pc.avg)[f] %in% rownames(ode) ) {
		pc.avg[f,"cpm.c"] <- ode[ rownames(pc.avg)[f], "cpm.c" ] 
	}
}

# Read in the full-length non-accessible ODEs to use as a "neutral" background
fl <- read.table("rm/fl_ode.txt")

s <- subset(pc.avg, grepl("1E", rownames(pc.avg)) & !is.na(cpm.c))
a <- s; a$repName <- paste0(a$repName, "_", a$cpm.c); a <- rbind(a, subset(pc.avg, repName == "ORR1E" & rownames(pc.avg) %in% fl$V4)) 
a$repName <- factor(a$repName, levels = c("ORR1E",paste0("ORR1E_", 8:1)))
meta <- data.frame("repName" = factor( levels(a$repName), levels = levels(a$repName) ), 
	"p.val" = c(NA, p.adjust(sapply( levels(a$repName)[-1], function(i){wilcox.test(subset(a, repName == i)$avg, subset(a, repName == "ORR1E")$avg, alternative = "two.sided")$p.value}), method = "BH")), 
	"num" = as.data.frame(table(a$repName))$Freq, 
	"y" = c(NA, sapply( levels(a$repName)[-1], function(i){max(subset(a, repName == i)$avg)+0.01})) )
meta$sig <- with(meta, ifelse(p.val < 0.0005, "***", ifelse(p.val < 0.005, "**", ifelse(p.val < 0.05, "*", "") ) ) )
meta$dir <- c(NA, sapply( levels(a$repName)[-1], function(i){ if (median(subset(a, repName == i)$avg) - median(subset(a, repName == "ORR1E")$avg) > 0) {return("up")} else {return("down")} }))
meta[1,"num"] <- paste0("~", round(meta[1,"num"]/1000,0), "k")

# Create the color theme for each comparison
col.c <- c(tableau_color_pal(type = "ordered-diverging", palette = "Sunset-Sunrise Diverging")(7), "gray69", "gray58") 

g1 <- ggplot(a, aes(x = repName, y = avg, fill = repName)) + geom_violin(linewidth = 1/.75/.pt, show.legend = F) + geom_hline(yintercept = median(subset(a, repName == "ORR1E")$avg), lty = "31", lwd = 1/.75/.pt) + 
	geom_boxplot(width = 0.2, outlier.shape = NA, fill = "white", fatten = 1/.75/.pt) + scale_fill_manual(values = col.c[9:1]) + theme_classic() + scale_x_discrete(name = "", labels = rev(c("Stem\n& Prog", "B", "T", "NK", "ILC3", "MF", "DC", "Bkgd", "Full Length\nORR1E"))) +
	labs(y = "Average phastCons score") + geom_text(data = meta, aes(x = repName, y = -0.04, label = num), size = 8/.pt) + geom_text(data = meta, aes(x = repName, y = y, label = sig, color = dir), size = 8/.pt, angle = 270, hjust = 0.5, fontface = "bold", show.legend = F) + 
	theme(axis.text = element_text(size = 8, color = "black"), axis.title = element_text(size = 8, color = "black"), axis.ticks.length = unit(.04, "in"), axis.ticks = element_line(linewidth = 1/.75/.pt, color = "black"), axis.line = element_line(linewidth = 1/.75/.pt, color = "black")) +
	ylim(c(-0.05,1.02)) + coord_flip() + scale_color_manual(values = c("up" = "dark orange", "down" = "black"))

rm(meta, a, s, f, fl, pc.avg, pc.out, col.c, ode)


## B

# Load in the packages you'll need for dealing with orthologous ODE sequence between mouse and rat the correlations
suppressPackageStartupMessages({
	library(fastSave); library(tidyr); library(dplyr); library(doParallel); library(stringr); library(ggthemes); library(scales); library(ggplot2)
	library(Biostrings); library(ComplexHeatmap)
})

# Generate some objects that will contain the look-up table that you need to go from genomic to consensus coordinates for ORR1E and ORR1D2!
comb <- readDNAStringSet("orthology/ode_s1000.fa"); rn7 <- readDNAStringSet("ode_s1000_match_rn7.fa")

# Read in the created substitution matrix to recapitulate the NEEDLE alignments
dnafull <- as.matrix(read.table("rm/DNAfull.mat", header = T))

res <- mclapply(1:length(rn7), function(i) {
	test <- pairwiseAlignment(rn7[i], comb[which(names(comb) == names(rn7)[i])], substitutionMatrix = dnafull, gapOpening = 10, gapExtension = -0.5, type = "global")
	val <- as.data.frame(matrix(0,nrow=2,ncol=width(alignedPattern(test))), row.names = c(names(rn7[i]), "mm10"))
	val[1,] <- as.matrix(alignedPattern(test)); val[2,] <- as.matrix(alignedSubject(test))
	return(val)
}, mc.cores = 12)

# The difference here is that you only want to know where the rat equivalent matches up with the mouse version. liftOver is quite lenient when
# it comes to providing an equivalent sized region, so you want to trim it down to what actually matters
bool <- as.data.frame(matrix(0,nrow=length(rn7),ncol=max(width( comb[ which(names(comb) %in% names(rn7)) ] ))), row.names = names(rn7))
lut <- as.data.frame(matrix(-1,nrow=length(rn7),ncol=max(width( comb[ which(names(comb) %in% names(rn7)) ] ))), row.names = names(rn7))
tmp <- mclapply(1:nrow(lut), function(f) {
	cc <- 0; cb <- 1; cr <- 0; val <- list(lut[f,],bool[f,]); clen <- width(comb[ which(names(comb) == rownames(lut)[f]) ])
	while ( cc < clen ) {
		if ( res[[f]][1,cb] != "-" ) { # rat sequence is *not* gapped, increment the counter
			cr <- cr + 1
		}

		if ( res[[f]][2,cb] != "-" ) { # mm10 sequence is *not* gapped, append!
			cc <- cc + 1
			if ( res[[f]][1,cb] == "-" ) { # rat sequence is gapped
				m <- 0
			} else { # No gaps between mouse and rat
				if ( res[[f]][1,cb] == res[[f]][2,cb] ) { # check to see if the current rat sequence is identical to the mouse sequence. 
					m <- 1
				} else { # it's not the same, so false
					m <- -1
				}
			}

			# Mark the current mouse sequence with whatever it currently is in rat
			val[[1]][1,cc] <- cr

			# Mark whether the current mouse base matches the sequence in rat, or has no matching sequence at all...
			val[[2]][1,cc] <- m
		} 

		# Append the current position in the alignment
		cb <- cb + 1
	}
	return(val)
}, mc.cores = 12)
lut <- bind_rows(lapply(tmp, `[[`, 1))
bool <- bind_rows(lapply(tmp, `[[`, 2))
rm(tmp)

# Start by reading in the list of ORR1 intervals 
ode <- read.table("atac/rObject_ode.txt", sep = "\t", header = T); ode$cpm.c <- factor(ode$cpm.c, levels = 1:8); ode$peak <- gsub(".*\\|", "", ode$pid)

# Load in the file that has all of the start and stop (relative to mm10) to filter each syntenic match for motif enrichment
a <- read.table("orthologyode_blast_s1000.bed"); a$stop <- unlist(mclapply(1:nrow(lut), function(i){sum(lut[i,] != "-1")-1000}, mc.cores = 12))
a$flank <- unlist(mclapply(1:nrow(a), function(i){ 
	b <- 0

	# Left flank math
	if ( a[i,"V7"] < 1000 ) { 
		if ( a[i,"V8"] >= 1000 ) { b <- b + 1001 - a[i,"V7"] }
		if ( a[i,"V8"] < 1000 ) { b <- b + ( a[i,"V8"] - a[i,"V7"] + 1 ) }	 
	} 
	# Right flank math
	if ( a[i,"V8"] > a[i,"stop"] ) { 
		if ( a[i,"V7"] <= a[i,"stop"] ) { b <- b + ( a[i,"V8"] - a[i,"stop"] + 1 ) } 
		if ( a[i,"V7"] > a[i,"stop"] ) { b <- b + ( a[i,"V8"] - a[i,"V7"] + 1 ) } 
	}

	# Return the amount of sequence in either flank
	return( b )
}, mc.cores = 16) )

a$ele <- unlist(mclapply(1:nrow(a), function(i){ 
	# Determine what value to use as the beginning 
	if ( a[i,"V7"] <= 1001 ) { # The start includes the full beginning of the element. Use 1001 as the beginning
		beg <- 1000 
	} else { # The start does not include the full beginning of the element. Use start as the beginning
		beg <- a[i,"V7"] - 1
	}

	# Determine what value to use as the end
	if ( a[i,"V8"] >= a[i,"stop"] ) { # The end extends beyond the stop boundary. Use stop as the end of the element
		e <- a[i,"stop"]
	} else { # The end does not include the stop boundary. Use end as the end of the element
		e <- a[i,"V8"]
	}

	# Return the amount of sequence in the syntenic ODE
	if ( e - beg > 0 ) { return( e - beg ) } else { return( 0 ) }

}, mc.cores = 16) )

# Finally, figure out which elements to filter and which to keep
a$filt <- with(a, ifelse( ele > (stop - 1000) / 2 & flank >= 500, "pass", "fail" ) )

# Generate some objects that will contain the look-up table that you need to go from genomic to consensus coordinates for ORR1E and ORR1D2!
rn7 <- readDNAStringSet("orthology/ode_blast_s1000.fa")

# For each cluster, iterate through all of the elements, subseq them, and then write them to .fa file
for (f in 1:8) {

	s <- subset(ode, repName == "ORR1E" & cpm.c == f & id %in% rownames(lut) & id %in% subset(a, filt == "pass")$V4); it <- unique(s$id)
	s1 <- DNAStringSet()
	for (g in it) {
		h <- subset(s, id == g)[1,]; if ( lut[ g, 1000 ] == 0 ) { i <- 1 } else { i <- lut[ g, 1000 ] }
		s1 <- c(s1, subseq( rn7[g], start = i, end = lut[ g, subset(a, V4 == g)[1,"stop"] ] ) )
	}
	writeXStringSet(s1, paste0("orthology/ORR1E_blast_cpm.c_", f, ".fa"))

	if (f < 7) {
		s <- subset(ode, repName == "ORR1D2" & cpm.c == f & id %in% rownames(lut) & id %in% subset(a, filt == "pass")$V4); it <- unique(s$id)
		s1 <- DNAStringSet()
		for (g in it) {
			h <- subset(s, id == g)[1,]; if ( lut[ g, 1000 ] == 0 ) { i <- 1 } else { i <- lut[ g, 1000 ] }
			s1 <- c(s1, subseq( rn7[g], start = i, end = lut[ g, subset(a, V4 == g)[1,"stop"] ] ) )
		}
		writeXStringSet(s1, paste0("orthology/ORR1D2_blast_cpm.c_", f, ".fa"))
	}
}


### NOTE: The below output is created from running blast between mm10 and rn7 on the ODEs with accessibility across the ImmGen data. Woot woot 

df.tmp <- data.frame()
for (f in 1:7) { 
	if ( f < 6 ) {
		tmp <- read.table(paste0("homer/rat/orr1d2_", f, ".vs.6.txt"), sep = "\t", skip = 1); tmp$clust <- f; tmp$samp1 <- "orr1d2"; tmp$samp2 <- "bkgd" 
		df.tmp <- rbind(df.tmp, tmp)
		tmp <- read.table(paste0("homer/rat/orr1d2_6.vs.", f, ".txt"), sep = "\t", skip = 1); tmp$clust <- f; tmp$samp1 <- "bkgd"; tmp$samp2 <- "orr1d2"
		df.tmp <- rbind(df.tmp, tmp) 		
	}

	tmp <- read.table(paste0("homer/rat/orr1e_", f, ".vs.8.txt"), sep = "\t", skip = 1); tmp$clust <- f; tmp$samp1 <- "orr1e"; tmp$samp2 <- "bkgd" 
	df.tmp <- rbind(df.tmp, tmp)
	tmp <- read.table(paste0("homer/rat/orr1e_8.vs.", f, ".txt"), sep = "\t", skip = 1); tmp$clust <- f; tmp$samp1 <- "bkgd"; tmp$samp2 <- "orr1e"
	df.tmp <- rbind(df.tmp, tmp) 		
}

# With the individual comparisons loaded in, get rid of any duplicate motifs (ie the RORgt one that is rude) and wrangle the object for 
# plotting purposes
df.tmp <- unique(df.tmp); rm(tmp); df.tmp$V7 <- as.numeric(gsub("%", "", df.tmp$V7)); df.tmp$V9 <- as.numeric(gsub("%", "", df.tmp$V9))
df.tmp <- separate(df.tmp, col = "V1", into = c("tmp", "data", "origin"), sep = "/")
df.tmp <- separate(df.tmp, col = "tmp", into = c("motif", "motif.fam"), sep = "\\(")
df.tmp$motif.fam.short <- gsub("\\?", "", gsub("\\?,", "", gsub(").*", "", df.tmp$motif.fam)))
df.tmp$motif.fam.short <- ifelse( is.na(df.tmp$motif.fam.short) | df.tmp$motif.fam.short == "", "Other", df.tmp$motif.fam.short )
df.tmp[df.tmp$motif.fam.short == "forkhead","motif.fam.short"] <- "Forkhead"
df.tmp[df.tmp$motif.fam.short == "ETS:IRF","motif.fam.short"] <- "ETS"
df.tmp$motif.fam <- paste0("(", df.tmp$motif.fam); df.tmp$motif.uniq <- paste0(df.tmp$motif, df.tmp$motif.fam, "/", df.tmp$data)

# Adjust the p-values!
df.tmp$padj <- p.adjust(exp(df.tmp$V4), method = "BH")

# Finally, create a copy of the above "df.tmp" object so you can overwrite it below!
df.mot.enrich <- df.tmp; rm(df.tmp)

# Load in the version of the above done in mouse
df.mot.enrich.m <- read.table(file = "homer/rObject_df.mot.enrich.txt", sep = "\t", header = T)

# In the above rat-based object "df.mot.enrich", add a column that determines whether that enrichment exists in mouse or not
df.mot.enrich$mouse <- F
for (f in 1:nrow(df.mot.enrich)) {
	sub.m <- subset(df.mot.enrich.m, motif.uniq == df.mot.enrich[f,"motif.uniq"] & clust == df.mot.enrich[f,"clust"] & samp1 == df.mot.enrich[f,"samp1"] & samp2 == df.mot.enrich[f,"samp2"])
	if ( nrow(sub.m) > 0 ) {
		if ( sub.m[1, "padj"] <= 0.05 ) {
			df.mot.enrich[f,"mouse"] <- T
		} 
	} 
}

# Save the above object so that you don't need to regenerate it in the future
write.table(df.mot.enrich, file = "homer/rObject_syntenicRatODE_df.mot.enrich.txt", sep = "\t", quote = F, row.names = F)

# Just gather the same motifs as are in the mouse motif list. Yeah
mtp <- c("PU.1", "SpiB", "Fli1", "KLF6", "RUNX2", "EBF1", "E2A", "CTCF", "Maz", "ZFX", "Reverb", "RORgt", "KLF3", "PU.1:IRF8")
p <- subset(df.mot.enrich, padj <= 0.05 & motif %in% mtp & (samp1 == "orr1e" | samp2 == "orr1e")) %>% group_by(motif,clust) %>% arrange(padj) %>% slice_head(n = 1) %>% as.data.frame()
p$prat <- -log10(p$padj); p$lrat <- ifelse(p$samp1 == "orr1e", p$V5, -p$V5)
p.lut <- data.frame(row.names = 1:7, "name" = c("Stem\n& Prog", "B", "T", "NK", "ILC3", "MF", "DC")); p$lab <- factor(p.lut[p$clust,"name"], levels = rev(p.lut$name))
p$motif <- factor(p$motif, levels = mtp)

# Create fake data to make sure all the legend values are present in the final plot
f1 <- data.frame(x = factor("DC", levels = rev(p.lut$name)), y = factor("PU.1", levels = levels(p$motif)), f = factor(c("Yes", "No"), levels = c("Yes", "No")))

g2 <- ggplot(p, aes(x = motif, y = lab, size = -log10(padj), fill = lrat, color = mouse)) + geom_point(data = f1, aes(x = y, y = x, color = f), pch = 21, alpha = 0, inherit.aes = F) + geom_point(data = subset(p, mouse == T), pch = 21, stroke = 1/.pt/.5, color = "black") + geom_point(data = subset(p, mouse == F), pch = 21, color = "transparent") + theme_bw() + 
	scale_fill_gradient2(low = "blue", mid = "white", high = "dark orange", limits = c(-1.5,1.5), oob = squish, breaks = c(-1,0,1), labels = c("-1", "0", "1"), name = "log2 ratio") + 
	scale_size_area(limits = c(-log10(0.05),9), oob = squish, max_size = 5, breaks = c(-log10(0.05), 3, 6, 9), labels = c("5x10-2", "1x10-3", "1x10-6", "<1x10-9"), name = "Adjusted p-value") +
	scale_color_manual(name = "Enrichment\nin mouse?", labels = c("Yes", ""), values = c("black", "transparent")) + scale_y_discrete(limits = rev(p.lut$name)) +
	theme(panel.border = element_blank(), text = element_text(size = 8, color = "black"), axis.text = element_text(size = 8, color = "black"), axis.title = element_text(size = 8, color = "black"), axis.ticks = element_blank(), legend.title = element_text(size = 8, color = "black"), legend.text = element_text(size = 8, color = "black"), legend.position = "bottom" ) + labs(x = "", y = "") +
	labs(color = "Expressed") + guides(fill = guide_colorbar(order = 1, title.position="top", title.hjust = 0.5, barwidth = unit(0.5, "in"), barheight = unit(.1, "in")), size = guide_legend(override.aes = list(pch = 21, fill = "black", color = "transparent"), order = 2, title.position="top", title.hjust = 0.5), color = guide_legend(override.aes = list(pch = 21, fill = "white", size = 5, color = c("black", "transparent"), stroke = 1/.pt/.5, alpha = 1), order = 3, title.position="top", title.hjust = 0.5)) +
	coord_fixed(clip = "off") + scale_x_discrete(position = "top", guide = guide_axis(angle = 45), drop = F)

# With the above generated, remove all the objects you no longer need
rm(p, mtp, p.lut, df.mot.enrich.m, df.mot.enrich, f1)



## C & D

# Set the number of cores to use for speed
setDTthreads(8)

# Load in the consensus sequences and combined ODE .fa files
cons <- readDNAStringSet("rm/dfam3.6_cons.fa"); cons <- cons[grepl("ORR1E|ORR1D2", names(cons)),]
comb <- c(readDNAStringSet("rm/ode.fa"), readDNAStringSet("rm/fl_ode.fa"))

# Read in the created substitution matrix
dnafull <- as.matrix(read.table("rm/DNAfull.mat", header = T))

res <- mclapply(1:length(comb), function(i) {
	test <- pairwiseAlignment(comb[i], cons[which(grepl(gsub(".*\\|", "", names(comb)[i]), names(cons)))], substitutionMatrix = dnafull, gapOpening = 10, gapExtension = -0.5, type = "global")
	val <- as.data.frame(matrix(0,nrow=2,ncol=width(alignedPattern(test))), row.names = c(names(comb[i]), "cons"))
	val[1,] <- as.matrix(alignedPattern(test)); val[2,] <- as.matrix(alignedSubject(test))
	return(val)
}, mc.cores = 16)

lut.all <- as.data.frame(matrix(-1,nrow=length(comb),ncol=max(width(comb))), row.names = names(comb))
lut.all <- bind_rows( mclapply(1:nrow(lut.all), function(f) {
	cc <- 1; cb <- 1; val <- lut.all[f,]; clen <- width(cons[grepl(gsub("_.*", "", gsub(".*\\|", "", rownames(lut.all)[f])), names(cons))])
	for ( g in 1:ncol(res[[f]]) ) { 
		if ( res[[f]][1,g] == "-" ) { # Sequence is gapped, so increment the consensus counter
			cc <- cc + 1
		} else { # Sequence is not gapped. Append the matricies			
			if (cc > clen) { cc <- clen } # Make sure consensus doesn't go above its maximum value (when sequence extends past consensus?)

			if ( res[[f]][2,g] == "-" ) { # Consensus is gapped, add NA instead
				val[1,cb] <- NA;
			} else { # No gaps, so increment consensus counter
				val[1,cb] <- cc; cc <- cc + 1
			}
			# With another base dealt with, append the current base counter
			cb <- cb + 1
		}
	}
	return(val)
}, mc.cores = 16) )

# Write the above object to disk so you don't have to regenerate it again in the future
save.lbzip2(lut.all, file = "rm/rObject_lut.all.RDataFS", n.cores = 16)

# Load in the HOMER motifs 
tmp <- homerToPFMatrixList("homer/homerFilt.motifs")
mot.thresh <- sapply(tmp, function(i){return(log(2^i@tags$log2cut))}); pfms <- sapply(tmp, function(i){return(i@profileMatrix)})
names(pfms) <- names(mot.thresh) <- sapply(tmp, function(i){return( gsub("/Homer", "", i@ID) )})

### NOTE: This works for everything >0, except for OCT:OCT-short which has a long motif and a negative (!) score threshold. Gross.
meta <- data.frame(row.names = sapply(tmp, function(i){return(i@ID)}), "log.thresh" = mot.thresh, "proto.thresh" = mot.thresh * .5, "mot.length" = sapply(tmp, function(i){return(ncol(i@profileMatrix))}))

# Load in the object containing all of the motif locations along each consensus sequence
base.con.mot <- read.table("homer/rObject_base.con.mot.txt", header = T, sep = "\t")

# Load in the python-generated HOMER motif clusters
mot.clust <- read.table("homer/homer_mot_cluster_python.txt", header = T, sep = "\t", row.names = 1)

# Regenerate the object that *only* includes the ORR1E motifs that you care about so that you can easily make the figure. 
bmc <- data.frame()
for (f in c(33,131,79,114,8,75,129,85,12,11,36)) {
	s <- base.con.mot[base.con.mot$mot %in% rownames(mot.clust[mot.clust$cluster == f,,drop=F]),]
	s <- subset(s, qual != "no" & seq == "ORR1E")
	tmp <- as.data.frame(matrix("",nrow=1,ncol=359))
	for (g in which(s$qual == "proto")) { tmp[1,s[g,"start"]:s[g,"end"]] <- s[g,"qual"] }
	for (g in which(s$qual == "mot")) { tmp[1,s[g,"start"]:s[g,"end"]] <- s[g,"qual"] }

	# Now iterate through the columns and output a row that matches the stretch of a given character. Yeah!
	cur <- ""; st <- 1; w <- rev(c(33,131,79,114,8,75,129,85,12,11,36))
	for (g in 1:ncol(tmp)) {
		if ( tmp[1,g] != cur ) { # The current stretch has ended. Append to the main object and start another one
			bmc <- rbind(bmc, data.frame("seq" = "ORR1E", "mot.clust" = f, "qual" = cur, "xmin" = st - 0.5, "xmax" = g - 0.5, ymin = which(w == f)*.5, ymax = which(w == f)*.5 ) )
			st <- g; cur <- tmp[1,g]
		} 
	}
}
bmc <- subset(bmc, qual != ""); rm(s,tmp,f,g,cur,st,w)

# Generate the base ggplot for the motif ancestral and proto-motif locations 
bmc$xmin <- ifelse( bmc$xmin <= 1, 1, bmc$xmin )
bmc.p <- ggplot(bmc, aes(x=xmin, xend=xmax, y=ymin, yend=ymax, color = qual)) + geom_segment(linewidth = 3/.75/.pt) + scale_y_continuous(breaks = c(1:11)*.5, labels = rev(c("ETS", "KLF", "RUNX", "EBF1", "E2A", "CTCF", "Maz", "ZFX", "Reverb", "ROR", "PU.1:IRF"))) + theme_void() +
	coord_cartesian(xlim = c(1,359)) + scale_x_continuous(labels = rep("", 5), breaks = c(1,100,200,300,359)) + scale_color_manual(name = "Motifs in\nConsensus", labels = c("mot" = "Ancestral", "proto" = "Proto"), values = c("firebrick2", "dodgerblue")) +
	theme(axis.text = element_text(size = 8, color = "black"), legend.title = element_text(size = 8, color = "black"), legend.text = element_text(size = 8, color = "black"), aspect.ratio = 1/2)

# Load in the list of ODEs alongside the full-length elements
ode <- read.table("atac/rObject_ode.txt", sep = "\t", header = T)[,c(1:5,7)]; ode$cpm.c <- factor(ode$cpm.c, levels = 1:8); rownames(ode) <- ode$id
fl <- read.table("rm/fl_ode.txt"); colnames(fl)[1:4] <- c("chr", "start", "end", "id"); fl$repName <- gsub(".*\\|", "", fl$id); fl$cpm.c <- "fl"
ode <- rbind(ode, fl[,c(1:3,7,4,8)]); ode$cpm.c <- factor(ode$cpm.c, levels = c(1:8,"fl")) 

# Load in the large file containing the intersection of 60-way phastCons scores and TEs and only take ones within the above object
pc <- subset(as.data.table(fread("orthology/ode.phastCons60way.int.txt.gz", sep = "\t")), V4 %in% ode$id); setkey(pc, V4)

# Iterate through "pc" to create the underlying matrix for plotting things! As a reminder for phastCons scores, not all 
# elements will be present. The majority will have values, but some are NA cause they're just in mouse. Woot? 
it <- unique(pc$V4)
tmp <- as.data.frame(matrix(-1,nrow=1,ncol=max( width(cons) ) ))
pc.con.out <- bind_rows( mclapply(it, function(i) {
	val <- tmp; rownames(val) <- i; s <- pc[J(i),]; s$V8 <- ifelse( s$V8 < s$V2, s$V2, s$V8 ); s$V9 <- ifelse( s$V9 > s$V3, s$V3, s$V9 )

	# Have a case for no rows
	if ( nrow(s) == 0 ) { return(val); next }

	# For each row in the subset object, append the value (phastCons score) at "V10" to the correct location in "val"
	for ( f in 1:nrow(s) ) {

		# Get the start and end coordinates for the current interval
		ist <- s[f,V8] - s[f,V2] + 1; ied <- s[f,V9] - s[f,V2]

		# Check to see if any of the returned values in "lut" are NA. If yes, remove them. 
		r <- unique(unlist(lut.all[ i, ist:ied ])); r <- r[ !is.na(r) ]

		# If the length of "r" is not 0, set the values specified in "r" to the current value
		if ( length(r) != 0 ) { val[ i, r ] <- as.numeric(s[ f, "V10" ]) }
	}

	# Finally, return the matrix so that it can be combined into one big happy matrix
	return( val )
}, mc.cores = 64) )
rm(tmp, it, pc)

# Time to save the above object so you don't have to regenerate it
save.lbzip2(pc.con.out, file = "orthology/rObject_pc.con.out.RDataFS", n.cores = 16)

# Turn all of the "-1" values into NAs
pc.con.out[pc.con.out == -1] <- NA

p <- data.frame(); p2 <- data.frame()
for ( f in c(1:8,"fl") ) {
	a <- unique(subset(ode, grepl("ORR1E", id) & cpm.c == f)$id)
	
	for (g in 1:359) {
		ci <- confint(lm(pc.con.out[a,g] ~ 1), level = 0.95)
		p <- rbind(p, data.frame( "x" = g, "ymin" = ci[1], "ymax" = ci[2], "repName" = "orr1e", "group" = f ))
	}

	ks <- ksmooth(1:359, colMeans(pc.con.out[a,], na.rm = T), "normal", bandwidth = 10)
	p2 <- rbind(p2, data.frame("x" = ks$x, "y" = ks$y, "repName" = "orr1e", "group" = f))

	if ( f < 7 | f == "fl" ) {
		a <- unique(subset(ode, grepl("ORR1D2", id) & cpm.c == f)$id)
		for (g in 1:370) {
			ci <- confint(lm(pc.con.out[a,g] ~ 1), level = 0.95)
			p <- rbind(p, data.frame( "x" = g, "ymin" = ci[1], "ymax" = ci[2], "repName" = "orr1d2", "group" = f ))
		}

		ks <- ksmooth(1:370, colMeans(pc.con.out[a,], na.rm = T), "normal", bandwidth = 10)
		p2 <- rbind(p2, data.frame("x" = ks$x, "y" = ks$y, "repName" = "orr1d2", "group" = f))
	}
}
p$group <- factor(p$group, levels = c(1:8,"fl")); p2$group <- factor(p2$group, levels = c(1:8,"fl"))

# Get the colors you want for the below plot
col.c <- c(tableau_color_pal(type = "ordered-diverging", palette = "Sunset-Sunrise Diverging")(7), "gray69", "gray58")

# Create a named vector to use for renaming the cluster numbers into something meaningful! 
lut.vec <- c("1" = "Stem\n& Prog", "2" = "B", "3" = "T", "4" = "NK", "5" = "ILC3", "6" = "MF", "7" = "DC", "8" = "Bkgd", "fl" = "Full Length\nBkgd")

# Get the plots for the Stem & Prog and MF ODEs
g <- list()
peaks <- zoo::rollapply( zoo::as.zoo(subset(p2, repName == "orr1e" & group == 1)$y), 5, function(x) which.max(x)==2)
a <- bmc.p + geom_vline(xintercept = which(peaks == T)+1, linetype = "31", color = "black", linewidth = 1/.75/.pt)
p1 <- ggplot(subset(p, group %in% c(1,"fl") & repName == "orr1e")) + geom_ribbon(aes(x = x, group = group, ymin = ymin, ymax = ymax, fill = group), alpha = 0.5) + 
	geom_line(data = subset(p2, group %in% c(1,"fl") & repName == "orr1e"), aes(x = x, y = y, group = group, color = group), linewidth = 1/.75/.pt, alpha = 0.85, inherit.aes = F) + 
	scale_y_continuous(breaks = c(0.05,0.165,0.28), limits = c(min(subset(p, repName == "orr1e")$ymin), max(subset(p, repName == "orr1e")$ymax))) + scale_x_continuous(breaks = c(1,100,200,300,359)) + 
	geom_vline(xintercept = which(peaks == T)+1, linetype = "31", color = "black", linewidth = 0.75/.75/.pt) + scale_fill_manual(values = col.c[c(1,9)], name = "ORR1E Cluster", labels = lut.vec[c(1,9)]) + 
	scale_color_manual(values = col.c[c(1,9)], name = "ORR1E Cluster", labels = lut.vec[c(1,9)]) + theme_classic() + labs(y = "Average phastCons score", x = "Position in ORR1E Consensus") + 
	theme(aspect.ratio = 1/2, axis.title.y = element_text(size = 8, angle = 90, vjust = 10, color = "black"), axis.title.x = element_text(size = 8, color = "black"), axis.text = element_text(size = 8, color = "black"), legend.title = element_text(size = 8, color = "black"), legend.text = element_text(size = 8, color = "black"), axis.ticks.length = unit(.04, "in"), axis.ticks = element_line(linewidth = 1/.75/.pt, color = "black"), axis.line = element_line(linewidth = 1/.75/.pt, color = "black"))
g[[ length(g) + 1 ]] <- a / plot_spacer() / p1 + plot_layout(heights = c(2,-.5,2), guides = "collect")

peaks <- zoo::rollapply( zoo::as.zoo(subset(p2, repName == "orr1e" & group == 6)$y), 5, function(x) which.max(x)==2)
a <- bmc.p + geom_vline(xintercept = which(peaks == T)+1, linetype = "31", color = "black", linewidth = 1/.75/.pt)
p1 <- ggplot(subset(p, group %in% c(6,"fl") & repName == "orr1e")) + geom_ribbon(aes(x = x, group = group, ymin = ymin, ymax = ymax, fill = group), alpha = 0.5) + 
	geom_line(data = subset(p2, group %in% c(6,"fl") & repName == "orr1e"), aes(x = x, y = y, group = group, color = group), linewidth = 1/.75/.pt, alpha = 0.85, inherit.aes = 1) + 
	scale_y_continuous(breaks = c(0.05,0.165,0.28), limits = c(min(subset(p, repName == "orr1e")$ymin), max(subset(p, repName == "orr1e")$ymax))) + scale_x_continuous(breaks = c(1,100,200,300,359)) + 
	geom_vline(xintercept = which(peaks == T)+1, linetype = "31", color = "black", linewidth = 0.75/.75/.pt) + scale_fill_manual(values = col.c[c(6,9)], name = "ORR1E Cluster", labels = lut.vec[c(6,9)]) + 
	scale_color_manual(values = col.c[c(6,9)], name = "ORR1E Cluster", labels = lut.vec[c(6,9)]) + theme_classic() + labs(y = "Average phastCons score", x = "Position in ORR1E Consensus") + 
	theme(aspect.ratio = 1/2, axis.title.y = element_text(size = 8, angle = 90, vjust = 10, color = "black"), axis.title.x = element_text(size = 8, color = "black"), axis.text = element_text(size = 8, color = "black"), legend.title = element_text(size = 8, color = "black"), legend.text = element_text(size = 8, color = "black"), axis.ticks.length = unit(.04, "in"), axis.ticks = element_line(linewidth = 1/.75/.pt, color = "black"), axis.line = element_line(linewidth = 1/.75/.pt, color = "black"))
g[[ length(g) + 1 ]] <- a / plot_spacer() / p1 + plot_layout(heights = c(2,-.5,2), guides = "collect")


### NOTE: There were some errors involving being unable to convert object of class numeric into a grob, so you just saved each plot as a .pdf
### and combined them in Illustrator. RIP for not having everything *nearly* done in just R...
pdf("figures/fig6a_violinplot.pdf", height = 2.9, width = 3.4)
plot_grid( g1, labels = c("A"), label_size = 12 )
dev.off()

pdf("figures/fig6b_dotplot.pdf", height = 3.8, width = 4.4)
plot_grid( g2, labels = c("B"), label_size = 12 )
dev.off()

### NOTE: For the below figure, it has dotted lines at all local maxima as sfigs 9-10 do. Those were manually deleted and rectangles highlighting
### notable maxima were manually added in Illustrator, since it would be very difficult to add a rectangle between both figures in R as is
pdf("figures/fig6cd_combplot.pdf", height = 3.5, width = 8.6)
plot_grid( plotlist = g, ncol = 2, nrow = 1, labels = c("C", "D"), label_size = 12 )
dev.off()


### HERE


##############################
### Supplementary Figure 1 ###
##############################

# export PATH=/home/jdc397/lbzip2-2.5/bin:$PATH; module load R/4.1.3-r9; R --vanilla

# Load all the relevant libraries you need for this section
suppressPackageStartupMessages({ 
	library(ggbeeswarm); library(tidyr); library(ggplot2); library(RColorBrewer); library(data.table); setDTthreads(8)
	library(dplyr); library(Biostrings); library(ggtree); library(cowplot)
})

## A & B & Supplemental Table 1

g <- list()

# Load in the file containing all accessible TEs, all TE annotations, and perform some clean up
ate <- read.table("accessible.tes.atac.txt")
all.te <- data.frame(fread("../../../../022124_reAnno_mm10_dfam.txt")); all.te$id <- with(all.te, paste0(V1, ":", V2, "-", V3, "|", V4))
all.te$qual <- ifelse( all.te$id %in% ate$V4, "open", "closed" ); all.te$a <- all.te$V5
all.te <- separate(all.te, col = "a", into = c("repClass", "repFamily"), sep = "/")
all.te$repClass <- gsub("\\?", "", all.te$repClass)

# Manually change the class of "Retroposon" families (RMER1A/B/C) into LINE for simplicity, since they're L1-dependent and don't have LTRs
all.te$repClass <- ifelse( all.te$repClass == "Retroposon", "LINE", all.te$repClass )

# Now create a factor level for everything of note for the downstream plots 
all.te$filt.rc <- ifelse( all.te$repClass %in% c("DNA", "LTR", "LINE", "SINE"), all.te$repClass, "Other" )

# Generate the underlying plot object of the fraction of all accessible TEs separated by class
tab <- with(subset(all.te, filt.rc != "Other"), as.data.frame(table(repClass, qual))); tab$Frac <- tab$Freq / as.numeric(table(subset(all.te, filt.rc != "Other")$repClass)) * 100
tab$repClass <- factor(tab$repClass, levels = rev(c("DNA", "LINE", "SINE", "LTR")))
g[[1]] <- ggplot( subset(tab, qual == "open"), aes(x = Frac, y = repClass, fill = repClass) ) + geom_col(show.legend = F) + scale_fill_manual(values = rev(brewer.pal(4, "Set1"))) +
	scale_x_continuous(name = "Total accessible TE class annotations (%)", limits = c(0,5.9), expand = c(0,0), breaks = seq(0,5,1), labels = seq(0,5,1)) + labs(y = "TE class") +
	theme_classic() + geom_text(aes(x = Frac + .06, y = repClass, label = Freq), size = 8/.pt, hjust = "left") +
	theme(text = element_text(size = 8, color = "black"), axis.text = element_text(size = 8, color = "black"), axis.title = element_text(size = 8, color = "black"), axis.ticks.length = unit(.04, "in"), axis.ticks = element_line(linewidth = 1/.75/.pt, color = "black"), axis.line = element_line(linewidth = 1/.75/.pt, color = "black") )

rm(all.te, ate, tab)



# Read in the original "te.bound" object that includes summit intersections against TEs
tb <- fread("te.bound"); tb$peak.id <- with(tb, paste0(V1, ":", V2, "-", V3))
tb <- separate(tb, col = "V8", into = c("repClass", "repFamily"), sep = "/")
tb$repClass <- ifelse( tb$repClass == "Retroposon", "LINE", tb$repClass )
tb <- subset(tb, repClass %in% c("DNA", "LINE", "SINE", "LTR"))

# Iterate over each sample to gather the per-class and total proportions of TE-derived peaks
a <- unique(tb$V12); p <- data.frame()
for (f in a) {
	tmp <- read.table(paste0("/workdir/jdc397/1_currentWork/8_teCactus/3_summits/", f, "_summits.bed")); tmp$peak.id <- with(tmp, paste0(V1, ":", V2, "-", V3))
	val <- nrow(tmp); tab <- as.data.frame(table(subset(tb, V12 == f)$repClass)); tab$Var1 <- as.character(tab$Var1)
	tab[ nrow(tab) + 1, "Var1" ] <- "Total"; tab[ nrow(tab), "Freq" ] <- sum(tab$Freq[1:(nrow(tab) - 1)])
	tab$Frac <- tab$Freq / val * 100; tab$Sample <- f; p <- rbind(p, tab)
}

# Prepare the input to plot the mouse proportion of TEs and genes (bps) and also use it to incorporate its information into the above object
te <- subset(read.table("/workdir/jdc397/1_currentWork/8_teCactus/00_clean/mm10_reMask_repClass_bps.txt"), V1 != "Other")
te$frac <- round((te$V2 / 2725537669)*100, 1); te$V1 <- factor(te$V1, levels = c("DNA", "LINE", "SINE", "LTR"))
te <- arrange(te, V1); te$color <- RColorBrewer::brewer.pal(4,"Set1"); rownames(te) <- te$V1

# Add the genome fraction, calculate the difference in summit overlap to genome proportion and output the information for supplemental data 1
p$Genome.Frac <- ifelse( p$Var1 %in% te$V1, te[ as.character(p$Var1), "V2" ] / 2725537669 * 100, sum(te$V2) / 2725537669 * 100 )
p$Frac.Diff <- p$Frac - p$Genome.Frac
colnames(p)[1:2] <- c("repClass", "Summit.Overlap")
write.table(p, file = "tableS1_teSummitFrac.txt", row.names = F, col.names = T, quote = F, sep = "\t")

# Calculate the range of proportion of sample summits that are putatively derived from TEs for the paper
# a <- group_by(subset(p, repClass != "Total"), Sample) %>% summarize(s = sum(Frac)); range(a$s)
# 6.614276 25.594045

# Calculate the number of samples where at least one TE class (or total) is overrepresented
# length(unique(subset(p, Frac.Diff > 0)$Sample))
# 2

# Now add the plot to the list
g[[2]] <- ggplot(te, aes(x = V1, y = frac, fill = V1)) + geom_col(show.legend = F) + scale_x_discrete() + scale_y_continuous(limits = c(0,24.5), labels = seq(0,20,5), breaks = c(seq(0,20,5)), expand = c(0,0)) +
	geom_text(data = te, aes(x = V1, y = frac + 0.9, label = frac), color = "black", size = 8/.pt) + theme_classic() + scale_fill_manual(values = RColorBrewer::brewer.pal(4,"Set1")) +
	ylab("Genome fraction (%)") + xlab("") + theme(axis.text = element_text(size = 8, color = "black"), axis.title = element_text(size = 8, color = "black"), axis.ticks.length = unit(.04, "in"), axis.ticks = element_line(linewidth = 1/.75/.pt, color = "black"), axis.line = element_line(linewidth = 1/.75/.pt, color = "black"))

# Plot the first row of panels for this figure 
g1 <- plot_grid( plotlist = g, ncol = 2, labels = c("A","B"), label_size = 12 )

# Remove all the objects you no longer need
rm(a, f, g, p, tab, tb, te, tmp, val)


## C

# Read in all the sample outputs for the only "enhancer" peaks, and adjust the fisher output for all of the tests you ran as part of this
res <- read.table("/workdir/jdc397/1_currentWork/8_teCactus/00_clean/1_data/1_atac/encode/shuffleResults/TEshuffle_results.txt", header = T)

# Overwrite the column for the sample name
res$samp <- gsub("_onlyEnhancer", "", res$samp)

# Overwrite the fisher's test p-values with their "bonferroni-adjusted" values and filter based on those
res$fep <- p.adjust(res$fisher_enriched_padj, method = "bonferroni")
res$repClass <- ifelse( res$repClass == "Retroposon", "LINE", res$repClass )
p.bonf <- res %>% filter(perm_enriched_pval <= 0.05 & fep <= 0.05 & repClass %in% c("DNA", "LINE", "SINE", "LTR"))

# Read in all of the consensus sequences for dFam 
cons <- readDNAStringSet("/workdir/jdc397/1_currentWork/8_teCactus/dfam3.6_cons.fa")

# Iterate through the families that are significant, and find the right consensus sequence to use for an alignment to order the x-axis
fams <- unique(p.bonf$repName)
for (f in 1:length(fams)) {
	if ( f == 1 ) {
		tmp <- cons[which(grepl(paste0("^", fams[f], "#"), names(cons)))]
	} else {
		# If no hits (ie it's a LINE with consensus sequence chunks for orf2/5'/3') then check first for orf2, then 5', then 3'
		if ( length( which(grepl(paste0("^", fams[f], "#"), names(cons))) ) == 0 ) {
			if ( length ( which(grepl(paste0("^", fams[f], "_orf2#"), names(cons))) ) == 1 ) {
				tmp <- c(tmp, cons[which(grepl(paste0("^", fams[f], "_orf2#"), names(cons)))])
			} else if ( length ( which(grepl(paste0("^", fams[f], "_5end#"), names(cons))) ) == 1 ) {
				tmp <- c(tmp, cons[which(grepl(paste0("^", fams[f], "_5end#"), names(cons)))])
			} else {
				tmp <- c(tmp, cons[which(grepl(paste0("^", fams[f], "_3end#"), names(cons)))])
			}
		} else {
			tmp <- c(tmp, cons[which(grepl(paste0("^", fams[f], "#"), names(cons)))])
		}
	}
}

names(tmp) <- gsub("_3end", "", gsub("_5end", "", gsub("_orf2", "", names(tmp))))

# Output the consensus sequences that match your current list of elements
writeXStringSet(tmp, "tmp.fa", append=FALSE, compress=FALSE, format="fasta")

# Generate sequence alignments between all the above consensus sequences to get the order for the current figure
system("export OMP_NUM_THREADS=16")
system("/programs/mafft/bin/mafft --ep 0.123 --quiet --thread 12 tmp.fa > tmp.afa")
system("clipkit tmp.afa -o tmp.clipkit -l")
system("/programs/FastTree-2.1.11/FastTreeMP -nt -gtr -gamma -seed 8 -log tmp.clipkit.ftmp.log < tmp.clipkit > tmp.clipkit.ftmp.nwk")
system("rm tmp.fa tmp.afa tmp.clipkit")

# Now load in the tree output from above so that you can use it to order your dotplot. Also, rotate the branches so that ORR1E & ORR1D2 are
# next to each other, since when you look at the branch length (and annotations) ORR1B1 is younger than either
ggtree_plot <- ggtree(read.tree("tmp.clipkit.ftmp.nwk"), branch.length = "none") 

# Reorder the dotplot by the clusters generated above
p.bonf$repName <- factor(p.bonf$repName, levels = gsub("#.*", "", get_taxa_name(ggtree_plot)))
lab <- get_taxa_name(ggtree_plot); classlab <- gsub("-consensus", "", gsub("/.*", "", gsub(".*#", "", lab))); classlab[classlab == "Retroposon"] <- "LINE"
lut <- brewer.pal(4, "Set1"); names(lut) <- c("DNA", "LINE", "SINE", "LTR")

# Create an object to organize the different tissues and their order
s <- sort(unique(res$samp)); v <- unique(gsub("_[e,p][0-9].*", "", s)) 
lut.ig <- data.frame(row.names = s, "name" = gsub("_[e,p][0-9].*", "", s), "lab" = gsub(".*_", "", s), "col" = "#000000", "source" = "ENCODE")

# Manually add in break blocks so that you can get the desired plot visualization. What a pain :'(
val <- c( rownames(subset(lut.ig, name == "embryonic_facial_prominence")), "b1", "b2", rownames(subset(lut.ig, name == "forebrain")), "b3", "b4", rownames(subset(lut.ig, name == "heart")), "b5", "b6", rownames(subset(lut.ig, name == "hindbrain")), "b7", "b8", rownames(subset(lut.ig, name == "intestine")), "b9", "b10", rownames(subset(lut.ig, name == "kidney")), "b11", "b12", rownames(subset(lut.ig, name == "limb")), "b13", "b14", rownames(subset(lut.ig, name == "liver")), "b15", "b16", rownames(subset(lut.ig, name == "lung")), "b17", "b18", rownames(subset(lut.ig, name == "midbrain")), "b19", "b20", rownames(subset(lut.ig, name == "neural_tube")), "b21", "b22", rownames(subset(lut.ig, name == "stomach")) )
p.bonf$samp <- factor(p.bonf$samp, levels = val); val <- unique(lut.ig[,c("name","col")]); tmp <- val[,"col"]; names(tmp) <- val[,"name"]; val <- tmp; rm(tmp)
p.bonf$lab <- lut.ig[p.bonf$samp, "lab"]; p.bonf$fill <- lut.ig[p.bonf$samp, "col"]

### NOTE: To create the desired label for the x-axis coloring, you need to create some fake data to plot below and use the fill parameter to 
### add in the legend. Yeah!
f1 <- data.frame(x = "RMER15", y = "embryonic_facial_prominence_e11p5", f = factor( unique(p.bonf$repClass), levels = c("DNA", "LINE", "SINE", "LTR"))) 

# You created an object that contains the coordinates with which to generate rectangular boxes in the below plot. Yay!
p2 <- data.frame()
for (h in names(val)) {
	a <- which(levels(p.bonf$samp) %in% rownames(subset(lut.ig, name == h)))
	p2 <- rbind(p2, data.frame(ymin = min(a)-0.5, ymax = max(a)+0.5, xmin = 0, xmax = length(levels(p.bonf$repName))+1, fill = h))
}

# Now make that plot! 
g2 <- ggplot( arrange(p.bonf, desc(fep)), aes(x=samp, y=repName, color = -log10(fep), size = bound_uniq), alpha = 0.75) + 
	geom_point(data = f1, aes(x = y, y = x, fill = f), alpha = 0, size = 5, shape = 22, stroke = 1/.75/.pt, inherit.aes = F) + 
	scale_fill_manual(values = lut) + geom_segment(data = p2, aes(y = length(unique(p.bonf$repName))+1, x = ymin - 0.5, xend = ymax + 0.5), color = val[p2$fill], inherit.aes = F, show.legend = F, linewidth = 2/.75/.pt) + 
	geom_point(aes(size = ifelse( bound_uniq > 650, 650, ifelse( bound_uniq < 50, 50, bound_uniq )))) + scale_size(breaks = c(50, 350, 650), limits = c(50,650), range = c(1,3), labels = c("50", "350", "650")) + 
	cowplot::theme_minimal_grid(line_size = 0.5/.75/.pt, color = "gray69") + ylab("") + xlab("") + scale_color_gradient(low = "gray69", high = "dark orange", limits = c(-log10(0.05),-log10(1e-45)), oob = scales::squish, name = "-log10\nenriched\np-value", breaks = c(-log10(0.05), 15, 30, 45), labels = c("1.3", "15", "30", "45+")) +
        scale_x_discrete(drop = F, breaks = rownames(lut.ig)) + coord_cartesian(ylim = c(1,length(levels(p.bonf$repName))), clip = "off") +
	labs(size = "Accessible\nelements") + guides(color = guide_colorbar(order = 1, barheight = unit(0.5, "in"), barwidth = unit(0.1, "in")), size = guide_legend(order = 2), fill = guide_legend(title = "TE class", override.aes = list(alpha = 1, size = 3), order = 3), alpha = "none") +
	theme(legend.text = element_text(size = 8, color = "black"), legend.key.height = unit(0.15, "in"), legend.title = element_text(size = 8, color = "black"), plot.margin = margin(r = 5, t = 35, unit = "pt"), axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(color = lut[classlab], size = 8)) +
	geom_text(data = p2, aes(y = length(unique(p.bonf$repName))+1.6, x = ((ymax + ymin) / 2), label = gsub("_", " ", gsub("embryonic_facial_prominence", "efp", fill))), angle = 30, hjust = 0, inherit.aes = F, show.legend = F, size = 8/.pt, color = "black")

# Remove all the objects you no longer need
rm(f1, p2, a, h, fams, val, v, lab, lut.ig, classlab, cons, ggtree_plot, lut, s, res, f)



## D

a <- as.data.frame(data.table::fread("../../2_dnase/012925_TEshuffle_results.txt"))
a$tissue <- factor(gsub("_[1-9+].*", "", a$samp), levels = unique(gsub("_[1-9+].*", "", a$samp))); a$lab <- gsub("_", "\n", a$tissue)

# Get a sense of how many total samples there are for each tissue type
totis <- as.data.frame(table(unique(a[,c("samp", "lab")])$lab), row.names = 1)

# Remove all of the tissues where you only have a single replicate
a <- a[ which(!gsub("_[1-9+].*", "", a$lab) %in% rownames(totis[totis$Freq == 1, , drop = F]) ), ]

# Rename some of the tissues for ease of visualization
a$lab <- gsub("gastrocnemius", "gastr.", a$lab); a$lab <- gsub("adrenal\ngland", "adr.", a$lab)
a$lab <- gsub("embryonic\nfacial\nprominence", "efp", a$lab); a$lab <- gsub("urinary\nbladder", "bladder", a$lab)
a$lab <- gsub("\ncerebral\ncortex", " cc", a$lab); a$lab <- gsub("cerebellum", "cereb.", a$lab)
a$lab <- gsub("layer\nof\nhippocampus", "hipp.", a$lab); a$lab <- gsub("frontal\ncortex", "fr. cor.", a$lab)
a$lab <- gsub("\n", " ", a$lab); a$lab <- factor(a$lab, levels = unique(sort(a$lab)))

# Make it alphabetical within each kind of tissue
b <- unique(a$lab); a$lab <- factor(a$lab, levels = c("cereb.", "forebrain", "fr. cor.", "gastr.", "hindbrain", "hipp.", "left cc", "midbrain", "neural tube", "right cc", "b1", "adr.", "bladder", "embryo", "efp", "heart", "kidney", "limb", "ovary", "testis", "b2", "liver", "lung", "spleen", "thymus"))

# Identify which samples are enriched, depleted or not significant
a$sig <- factor(ifelse(a$perm_enriched_pval <= 0.05, "Enriched", ifelse(a$perm_depleted_pval <= 0.05, "Depleted", "Not sig")), levels = c("Enriched", "Not sig", "Depleted"))

# You created an object that contains the coordinates with which to generate rectangular boxes in the below plot. Yay!
tmp <- data.frame("name" = c("Brain", "Non-Brain", "Immune")); tmp$samp <- list( c("cereb.", "forebrain", "fr. cor.", "gastr.", "hindbrain", "hipp.", "left cc", "midbrain", "neural tube", "right cc"), c("adr.", "bladder", "embryo", "efp", "heart", "kidney", "limb", "ovary", "testis"), c("liver", "lung", "spleen", "thymus"))
p2 <- data.frame(); tmp$x <- 0
for ( f in 1:nrow(tmp) ) {
	d <- which(levels(a$lab) %in% tmp[f,"samp"][[1]])
	p2 <- rbind(p2, data.frame(x = min(d), xend = max(d), y = -6.65))
	tmp[f,"x"] <- mean(d)
}

# Generate the figure
p <- subset(a, repName == "ORR1D2") 
g3 <- ggplot(p, aes(x = lab, y = log2(bound_uniq/exp_uniq))) + geom_hline(yintercept = 0, lty = "31", color = "gray69", linewidth = 1/.75/.pt) + geom_beeswarm(data = p, aes(x = lab, y = log2(bound_uniq/exp_uniq), color = sig), size = 1, shape = 16, inherit.aes = F, cex = 0.5) +
	stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="crossbar", width=0.5, color = "gray69", fatten = 2, linewidth = 1/.75/.pt) + scale_color_manual(values = c("Enriched" = "#ff7518", "Not sig" = "black", "Depleted" = "dodgerblue"), name = "ORR1D2") + 
	theme_classic() + scale_x_discrete(drop = F, breaks = b) + xlab("") + ylab(expr("log"[2]*"(obs/exp)")) + guides(color = guide_legend(override.aes=list(size = 2))) +
	coord_cartesian(ylim = c(-2.5,2.5), clip = "off") + geom_segment(data = p2, aes(x = x, xend = xend, y = y), color = "black", inherit.aes = F, show.legend = F, linewidth = 1/.75/.pt) +
	geom_text(data = tmp, aes(x = x, y = -7.2, label = name), size = 8/.pt, color = "black", inherit.aes = F, show.legend = F) +
	theme(plot.margin=unit(c(5.5,5.5,40.5,20.5), 'pt'), text = element_text(size = 8, color = "black"), axis.text.x = element_text(size = 8, color = "black", angle = 45, vjust = 1, hjust = 1), axis.title.y = element_text(size = 8, color = "black"), axis.text.y = element_text(size = 8, color = "black"), legend.text = element_text(size = 8, color = "black"), axis.ticks.length = unit(.04, "in"), axis.ticks = element_line(linewidth = 1/.75/.pt, color = "black"), axis.line = element_line(linewidth = 1/.75/.pt, color = "black"))

pdf("sfig1.pdf", height = 9.75, width = 6.5)
plot_grid( g1, g2, g3, labels = c("", "C", "D"), ncol = 1, rel_heights = c(1.75, 5.5, 2), label_size = 12 )
dev.off()

rm(a, b, p, totis, d, f, g1, g2, g3, p.bonf, p2, tmp)




##############################
### Supplementary Figure 2 ###
##############################

# export PATH=/home/jdc397/lbzip2-2.5/bin:$PATH; module load R/4.1.3-r9; R --vanilla

# Load all the relevant libraries you need for this section
suppressPackageStartupMessages({ 
	library(data.table); setDTthreads(8); library(ggpubr); library(fastSave); library(cluster); library(RColorBrewer)
	library(ggplot2); library(cowplot); library(dplyr); library(circlize); library(ComplexHeatmap); library(ggthemes)
	library(scales)
})


## A

# Load in the lut so you know all the different samples to iterate over
lut.a <- read.table("../../immgen_heatmap_lut.txt", sep = "\t", header = T, comment.char = "")

# Load the qnorm output 
load.lbzip2(file = "rObject_l2.atac.cpm.qnorm.RDataFS", n.cores = 16)

# Load in the list providing the peaks which are overlapping ORR1E/ORR1D2 with at least one summit across one dataset
te.int <- subset(as.data.frame(fread("peak.te.int")), V8 %in% c("ORR1E", "ORR1D2")); te.int$id <- with(te.int, paste0(V5, ":", V6, "-", V7, "|", V8))
te.int$pid <- with(te.int, paste0(V1, ":", V2, "-", V3, "|", V4))

# Generate the color scale that recapitulates Pheatmap's default since it's cool
col.fun <- colorRamp2(seq(-3, 3, length.out = 7), rev(brewer.pal(n = 7, name = "RdYlBu")))

# Get the cell type order and colors ready for the heatmap ordering
col.a <- unique(lut.a$col); names(col.a) <- unique(lut.a$name)

# Load the subset for ORR1D2
r1d2.atac <- l2.atac.cpm.qnorm[unique(subset(te.int, V8 == "ORR1D2")$pid), rownames(lut.a)]

set.seed(8)
tmp <- draw( Heatmap(lut.a$name, row_title_rot = 0, col = structure(col.a, names = names(col.a)), name = "Cell Type", show_row_names = F, width = unit(3, "mm")) +
        Heatmap(scale(t(r1d2.atac[ , rownames(lut.a) ])), column_km = 6, use_raster = T, column_km_repeats = 200, show_column_dend = F, cluster_column_slices = FALSE, cluster_rows = F, show_row_names = F, column_title = "ORR1D2", show_column_names = F, name = "CS.CPM", col = col.fun, border = T) )
r1d2.col.cpm <- column_order(tmp)$CS.CPM

tmp <- as.data.frame(matrix("", ncol = 1, nrow = nrow(r1d2.atac) ))
for (i in 1:length(r1d2.col.cpm)) {
	tmp[r1d2.col.cpm[[i]],"V1"] <- i 
}

# NEW ORDER: 1 - Stem & Prog, 2 - B, 3 - T, 4 - NK, 5 - MF/DC, 6 - Random
val <- data.frame("i" = 1:6, "o" = c(2,6,5,4,1,3), row.names = 1)
r1d2.col.cpm <- factor(val[tmp$V1, "o"], levels = c(1:6))

# To figure out which column you'll add the annotation to, you need to calculate the floored average for each block
a <- unlist(lapply(unique(lut.a$name), function(i){ floor(mean(which(lut.a$name == i))) }))
b <- rep("", nrow(lut.a)); b[a] <- unique(lut.a$name)

### NOTE: The above names are not figure-final. They need greek symbols and you found it easier to just add them in illustrator

g1 <- grid.grabExpr( draw( Heatmap(t(scale(t(r1d2.atac[ , rownames(lut.a) ]))), use_raster = T, 
	top_annotation = HeatmapAnnotation("Name" = anno_text(b, rot = 45, just = "left", location = unit(1, "mm"), gp = gpar(fontsize = 8, color = "black")), "Cell Type" = lut.a$name, col = list("Cell Type" = col.a), annotation_name_side = "left", annotation_name_gp = gpar(col = "black", fontsize = 8), show_legend = c(F,F), simple_anno_size = unit(3, "mm")), 
	heatmap_legend_param = list(direction = "horizontal", title_position = "topcenter", legend_width = unit(1, "in"), title = "log2(CPM+1)", labels_gp = gpar(col = "black", fontsize = 8), title_gp = gpar(col = "black", fontsize = 8)), 
	row_split = r1d2.col.cpm, row_title = c("Stem\n& Prog", "B", "T", "NK", "MF\n& DC", "Bkgd"), row_title_rot = 0, show_row_dend = F, 
	cluster_row_slices = FALSE, cluster_columns = F, show_row_names = F, show_column_names = F, row_title_gp = gpar(col = "black", fontsize = 8),
	col = col.fun, border = T, row_gap = unit(3, "mm")),
	heatmap_legend_side = "bottom", padding = unit(c(0.5,0.5,0.5,1), "cm") ) )

rm(a, b, col.a, col.fun, i, l2.atac.cpm.qnorm, lut.a, r1d2.atac, r1d2.col.cpm, te.int, tmp, val)


## B & C

# Load in the lut so you know all the different samples to iterate over
lut.a <- read.table("../../immgen_heatmap_lut.txt", sep = "\t", header = T, comment.char = "")

# Load the qnorm output 
load.lbzip2(file = "rObject_l2.atac.cpm.qnorm.RDataFS", n.cores = 16)

# Load in the list providing the peaks which are overlapping ORR1E/ORR1D2 with at least one summit across one dataset
te.int <- subset(as.data.frame(fread("peak.te.int")), V8 %in% c("ORR1E", "ORR1D2")); te.int$id <- with(te.int, paste0(V5, ":", V6, "-", V7, "|", V8))
te.int$pid <- with(te.int, paste0(V1, ":", V2, "-", V3, "|", V4))

# Load the subset for ORR1E and ORR1D2
r1e.atac <- l2.atac.cpm.qnorm[unique(subset(te.int, V8 == "ORR1E")$pid), rownames(lut.a)]
r1d2.atac <- l2.atac.cpm.qnorm[unique(subset(te.int, V8 == "ORR1D2")$pid), rownames(lut.a)]

# Get the list for cowplot ready
g <- list()

# Start with ORR1E
set.seed(8)
gap_stat <- clusGap( scale(t(r1e.atac[ , rownames(lut.a) ])), kmeans, K.max = 10, B = 200)
gap <- gap_stat$Tab[, "gap"]; se <- gap_stat$Tab[, "SE.sim"]; decr <- diff(gap) <= 0
df <- as.data.frame(gap_stat$Tab, stringsAsFactors = TRUE); df$clusters <- as.factor(1:nrow(df)); df$ymin <- gap - se; df$ymax <- gap + se
g[[1]] <- ggline(df, x = "clusters", y = "gap", group = 1, color = "dark orange", point.size = 1) + geom_vline(xintercept = 8, linetype = "31", color = "black") + 
	geom_errorbar(aes_string(ymin = "ymin", ymax = "ymax"), width = 0.5, color = "black") + labs(y = "Gap statistic", x = "Number of clusters", title = "ORR1E") +
	scale_y_continuous(breaks = c(0.47, 0.52, 0.57)) + theme(plot.title = element_text(size = 8, color = "black", hjust = 0.5), axis.title = element_text(size = 8, color = "black"), axis.text = element_text(size = 8, color = "black"), 
	aspect.ratio = 1, axis.ticks.length = unit(.04, "in"), axis.ticks = element_line(linewidth = 1/.75/.pt, color = "black"), axis.line = element_line(linewidth = 1/.75/.pt, color = "black"))

# Now for ORR1D2!
set.seed(8)
gap_stat <- clusGap( scale(t(r1d2.atac[ , rownames(lut.a) ])), kmeans, K.max = 10, B = 200)
gap <- gap_stat$Tab[, "gap"]; se <- gap_stat$Tab[, "SE.sim"]; decr <- diff(gap) <= 0
df <- as.data.frame(gap_stat$Tab, stringsAsFactors = TRUE); df$clusters <- as.factor(1:nrow(df)); df$ymin <- gap - se; df$ymax <- gap + se
g[[2]] <- ggline(df, x = "clusters", y = "gap", group = 1, color = "dark orange", point.size = 1) + geom_vline(xintercept = 6, linetype = "31", color = "black") + 
	geom_errorbar(aes_string(ymin = "ymin", ymax = "ymax"), width = 0.5, color = "black") + labs(y = "Gap statistic", x = "Number of clusters", title = "ORR1D2") +
	scale_y_continuous(breaks = c(0.44, 0.49, 0.54)) + theme(plot.title = element_text(size = 8, color = "black", hjust = 0.5), axis.title = element_text(size = 8, color = "black"), axis.text = element_text(size = 8, color = "black"), 
	aspect.ratio = 1, axis.ticks.length = unit(.04, "in"), axis.ticks = element_line(linewidth = 1/.75/.pt, color = "black"), axis.line = element_line(linewidth = 1/.75/.pt, color = "black"))

# Plot the top right side of panels for this figure 
g2 <- plot_grid( plotlist = g, ncol = 2, labels = c("B","C"), label_size = 12 )

# Remove all the objects you no longer need
rm(decr, df, gap, gap_stat, lut.a, r1d2.atac, r1e.atac, se, l2.atac.cpm.qnorm, g, te.int)


## D

# Load the HOMER motif enrichment object
df.mot.enrich <- read.table(file = "homer/rObject_df.mot.enrich.txt", sep = "\t", header = T)

# Find the motifs that have a decent amount of occurances and are expressed in the accessible cell types
tf.list.1d2 <- unique(subset(df.mot.enrich, padj <= 5e-2 & V7 >= 10 & expr.all == T & samp1 == "orr1d2")$motif.uniq)

# Get the list of the top 3 TFs (each given motif cluster - as per your version of Vierstra code - is limited to top factor)
a <- sapply(1:5, function(i) {
	b <- subset(df.mot.enrich, samp1 == "orr1d2" & motif.uniq %in% tf.list.1d2 & clust == i & padj <= 0.05) %>% group_by(motif.clust) %>% arrange(padj) %>% slice_head(n = 1) %>% ungroup() %>% arrange(desc(V5)) %>% head(3) %>% as.data.frame()
	return( as.character(b$motif.uniq) )
})
b <- unique(as.character(unlist(a))) 


# Subset the main file to get a filtered list along with necessary columns for the plot 
p <- subset(df.mot.enrich, padj <= 0.05 & motif.uniq %in% b & (samp1 == "orr1d2" | samp2 == "orr1d2")) %>% arrange(padj)
p$lrat <- ifelse(p$samp1 == "orr1d2", p$V5, -p$V5); p.lut <- data.frame(row.names = 1:5, "name" = c("Stem\n& Prog", "B", "T", "NK", "MF\n& DC"))

# Before you plot, turn it into a matrix to perform hierarchical clustering on and get the columns rows. Yeah
p$lab <- factor(p.lut[p$clust,"name"], levels = p.lut$name)
a <- reshape2::dcast(p, motif ~ lab, value.var = "padj"); a[is.na(a)] <- 1; row.names(a) <- a$motif; a$motif <- NULL 
dist_mat <- dist(scale(a), method = 'euclidean'); hclust_avg <- hclust(dist_mat, method = 'average')
# plot(hclust_avg)

# Set the factor order for the columns based on the hclust output. Remove the motifs "RUNX2" as its duplicative of RUNX1 and reorder the columns
# to approximate a sliding view from top to bottom 
p$motif <- factor(p$motif, levels = hclust_avg$labels[hclust_avg$order][c(5:6,8,4,9:13,1:3)])
p <- subset(p, motif %in% hclust_avg$labels[hclust_avg$order][c(-7)])

# Now reset the label levels for the plot
p$lab <- factor(p.lut[p$clust,"name"], levels = rev(p.lut$name))
f1 <- data.frame(x = factor("MF\n& DC", levels = rev(p.lut$name)), y = factor("ERG", levels = levels(p$motif)), f = factor(c("Yes", "No"), levels = c("Yes", "No")))

# Create the various themes you're thinking of using. Although you really don't want to use the rainbow colors, it is pretty distinguishable
col.c <- c(tableau_color_pal(type = "ordered-diverging", palette = "Sunset-Sunrise Diverging")(7)[1:5], "gray69"); names(col.c) = c(p.lut$name, "bkgd")

g3 <- ggplot(p, aes(x = motif, y = lab, size = -log10(padj), fill = lrat)) + geom_point(data = f1, aes(x = y, y = x, color = f), pch = 21, alpha = 0, inherit.aes = F) + geom_point(data = subset(p, expr.all == T), pch = 21, stroke = 1/.pt/.5, color = "black") + geom_point(data = subset(p, expr.all == F), pch = 21, color = "transparent") + theme_bw() + 
	scale_fill_gradient2(low = "blue", mid = "white", high = "dark orange", limits = c(-1.5,1.5), oob = squish, breaks = c(-1,0,1), labels = c("-1", "0", "1"), name = "log2 ratio") + 
	scale_size_area(limits = c(-log10(0.05),15), oob = squish, max_size = 5, breaks = c(-log10(0.05), 5, 10), labels = c("5e-2", "1e-5", "<1e-10"), name = "Adj p-val") +
	scale_color_manual(name = "Expressed", labels = c("Yes", ""), values = c("black", "transparent")) +
	theme(panel.border = element_blank(), text = element_text(size = 8, color = "black"), axis.text = element_text(size = 8, color = "black"), axis.ticks = element_blank(), axis.text.y = element_text(size = 8, color = "black"), legend.position="bottom" ) + labs(x = "", y = "") +
	guides(fill = guide_colorbar(order = 1, title.position="top", title.hjust = 0.5, barwidth = unit(0.5, "in"), barheight = unit(.1, "in")), size = guide_legend(override.aes = list(pch = 21, fill = "black", color = "transparent"), order = 2, title.position="top", title.hjust = 0.5), color = guide_legend(override.aes = list(pch = 21, size = 5, fill = "white", color = c("black", "transparent"), stroke = 1/.pt/.5, alpha = 1), order = 3, title.position="top", title.hjust = 0.5)) +
	coord_fixed(clip = "off") + scale_x_discrete(position = "top", guide = guide_axis(angle = 45)) + scale_y_discrete(position = "right")

# Combine panels for easier arranging via cowplot
g2.5 <- plot_grid( g2, g3, labels = c("", "D"), rel_heights = c(.6,1), ncol = 1, label_size = 12 )

# Finally, plot all the R-based panels. The heatmaps are generated using DeepTools and needs additional cleanup using illustrator 
pdf("sfig2ad_cowplot.pdf", height = 4.5, width = 7)
plot_grid( g1, g2.5, labels = c("A", ""), rel_widths = c(.85, 1), ncol = 2, label_size = 12 )
dev.off()

# Remove all the objects you no longer need
rm(col.c, df.mot.enrich, dist_mat, f1, g1, g2, g2.5, g3, hclust_avg, p, p.lut, tf.list.1d2)


## E 

### NOTE: The code for this section is instead in the file "code_final" since it uses DeepTools for the heatmap generation and row ordering




##############################
### Supplementary Figure 3 ###
##############################

# export PATH=/home/jdc397/lbzip2-2.5/bin:$PATH; module load R/4.1.3-r9; R --vanilla

# Load all the relevant libraries you need for this section
suppressPackageStartupMessages({ 
	library(ComplexHeatmap); library(ggthemes); library(dplyr); library(ggplot2); library(cowplot); library(PWMEnrich); library(TFBSTools)
	library(monaLisa)
}) 


## A

# Figure out how many ODEs are bound in at least 2 of the 6 ChIP-seq experiments, and of those how many are bound by at least PU.1 in progenitors
tmp <- c("PU.1_prog", "Pax5_B", "Runx1_CD8+", "RORg_Th17", "PU.1_MF", "IRF8_CD8+DC")
a <- read.table("figure3_ode_int_chip_vals.txt", sep = "\t"); colnames(a) <- c("chr", "start", "end", "repName", "id", "cpm.c", tmp)
rownames(a) <- a$id; b <- lapply(tmp, function(i){ return( rownames(a[a[,i] >= 1 & a$repName == "ORR1D2", ]) ) }); names(b) <- tmp

m <- make_comb_mat(b); m <- m[comb_size(m) >= 10]; ss <- set_size(m); cs <- comb_size(m)

# Create the various themes you're thinking of using. Although you really don't want to use the rainbow colors, it is pretty distinguishable
col.c <- c(tableau_color_pal(type = "ordered-diverging", palette = "Sunset-Sunrise Diverging")(7)[1:5], "gray69")

# Load in the ODE annotations
ode <- read.table("cpm_heatmap_ode_anno.txt", sep = "\t", header = T)[,c(1:5,7)]; rownames(ode) <- ode$id

bt <- bind_rows(
	lapply(names(comb_degree(m)), function(i){
		data.frame( "id" = extract_comb(m, comb_name = i), "code" = i ) 
	})
)
bt2 <- bt

bt$cpm.c <- ode[bt$id,"cpm.c"]; bt <- table(bt$code, bt$cpm.c); bt <- matrix(bt, ncol=ncol(bt), dimnames=dimnames(bt))
bt <- bt[names(comb_degree(m)),]

us <- UpSet(m, set_order = order(ss), comb_order = order(comb_degree(m), -cs), pt_size = unit(8, "pt"), lwd = 1.34,
	top_annotation = HeatmapAnnotation("ORR1D2 Overlaps" = anno_barplot(bt, ylim = c(0, max(rowSums(bt))*1.15), 
		border = FALSE, gp = gpar(fill = col.c), height = unit(1, "in"), axis_param = list(gp = gpar(fontsize = 8, color = "black"))),
        	annotation_name_side = "left", annotation_name_rot = 90, annotation_name_gp= gpar(fontsize = 8, color = "black")),
	right_annotation = NULL,
	row_names_gp = gpar(fontsize = 8, color = "black"),
)
us <- draw(us); od <- column_order(us)

g1 <- grid.grabExpr( {
	draw(us);
	decorate_annotation("ORR1D2 Overlaps", {
	    grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(3, "pt"), 
        	default.units = "native", just = c("left", "centre"), 
	        gp = gpar(fontsize = 8, col = "black"), rot = 90)
	})
})

rm(a, b, bt, bt2, col.c, cs, m, od, ode, ss, tmp, us)


## B & C

# Load in the binding information across ODEs
tmp <- c("PU.1_prog", "Pax5_B", "Runx1_CD8+", "RORg_Th17", "PU.1_MF", "IRF8_CD8+DC")
a <- read.table("figure3_ode_int_chip_vals.txt", sep = "\t"); colnames(a) <- c("chr", "start", "end", "repName", "id", "cpm.c", tmp)
rownames(a) <- a$id

# Create the color theme you'll be using
col.c <- c(tableau_color_pal(type = "ordered-diverging", palette = "Sunset-Sunrise Diverging")(7)[1:5], "gray69")

g <- list()

# What's the fraction of ORR1D2s bound by at least one factor?
p <- bind_rows( lapply(1:6, function(i) {
	return( data.frame( "cpm.c" = i, "bound" = sum( rowSums( subset(a, repName == "ORR1D2" & cpm.c == i)[,tmp] ) > 0 ), "total" = nrow(subset(a, repName == "ORR1D2" & cpm.c == i)) ) )
}) )
p$frac <- round(p$bound / p$total * 100, 2); p$cpm.c <- factor(p$cpm.c, levels = 6:1)

# Of the above, how many of them are bound by at least PU.1? 
q <- bind_rows( lapply(1:6, function(i) {
	return( data.frame( "cpm.c" = i, "bound" = sum( rowSums( subset(a, repName == "ORR1D2" & cpm.c == i & PU.1_prog >= 1)[,tmp] ) > 0 ), "total" = nrow(subset(a, repName == "ORR1D2" & cpm.c == i)) ) )
}) )
q$frac <- round(q$bound / q$total * 100, 2); q$cpm.c <- factor(q$cpm.c, levels = 6:1)

# Run stats to see if the proportion of cell type-specific ODEs bound is significantly higher than background
sig <- bind_rows( lapply(1:5, function(i) {
	return( data.frame( "cpm.c" = i, "sig" = fisher.test( matrix(c( p[i,"bound"], p[i,"total"] - p[i,"bound"], p[6,"bound"], p[6,"total"] - p[6,"bound"] ), nrow = 2, byrow = T), alternative = "greater" )$p.value, "frac" = p[i,"frac"] ) )
}) )
sig$sig <- p.adjust(sig$sig, method = "BH"); sig$lab <- ifelse( sig$sig <= 0.05, "*", "" ); sig$cpm.c <- factor(sig$cpm.c, levels = 5:1)

r <- 82
g[[1]] <- ggplot(p, aes(x = cpm.c, y = frac, fill = cpm.c)) + geom_col(show.legend = F) + geom_col(data = q, aes(x = cpm.c, y = frac), fill = "transparent", color = "black", linewidth = 1/.75/.pt, show.legend = F) +
	scale_y_continuous(breaks = seq(0, 75, 25), limits = c(0,r), labels = seq(0, 75, 25), expand = c(0,0)) + scale_fill_manual(values = rev(col.c)) + theme_classic() + 
	labs(x = "", y = "% bound by at least one TF") + scale_x_discrete(labels = rev(c("Stem\n& Prog", "B", "T", "NK", "MF\n& DC", "Bkgd"))) + geom_text(aes(x = cpm.c, y = frac + (r * 0.07), label = bound), size = 8/.pt, color = "black", angle = 270) + 
	geom_text(data = sig, aes(x = cpm.c, y = frac + (r * 0.13), label = lab), angle = 270, color = "dark orange", size = 8/.pt, fontface = "bold", inherit.aes = F) + geom_text(data = q, aes(x = cpm.c, y = frac - (r * 0.07), label = bound), size = 8/.pt, color = "black", angle = 270) + 
	theme(axis.title = element_text(size = 8, color = "black"), axis.text = element_text(size = 8, color = "black"), axis.ticks.length.x = unit(.04, "in"), axis.ticks.y = element_blank(), axis.ticks.x = element_line(linewidth = 1/.75/.pt, color = "black"), axis.line = element_line(linewidth = 1/.75/.pt, color = "black")) + 
	guides(fill = guide_legend(reverse=TRUE)) + coord_flip()

# What's the fraction of ORR1D2s bound by at least two factors? 
p <- bind_rows( lapply(1:6, function(i) {
	return( data.frame( "cpm.c" = i, "bound" = length(which(rowSums(subset(a, repName == "ORR1D2" & cpm.c == i)[,7:12] >= 1) >= 2)), "total" = nrow(subset(a, repName == "ORR1D2" & cpm.c == i)) ) )
}) )
p$frac <- round(p$bound / p$total * 100, 2); p$cpm.c <- factor(p$cpm.c, levels = 6:1)

# Of the above, how many of them are bound by at least PU.1? 
q <- bind_rows( lapply(1:6, function(i) {
	return( data.frame( "cpm.c" = i, "bound" = length(which(rowSums(subset(a, repName == "ORR1D2" & cpm.c == i & PU.1_prog >= 1)[,8:12] >= 1) >= 1)), "total" = length(which(rowSums(subset(a, repName == "ORR1D2" & cpm.c == i)[,7:12] >= 1) >= 2)) ) )
}) )
q$frac <- round(q$bound / p$total * 100, 2); q$cpm.c <- factor(q$cpm.c, levels = 6:1)

# Run stats to see if the proportion of cell type-specific ODEs bound is significantly higher than background
sig <- bind_rows( lapply(1:5, function(i) {
	return( data.frame( "cpm.c" = i, "sig" = fisher.test( matrix(c( p[i,"bound"], p[i,"total"] - p[i,"bound"], p[6,"bound"], p[6,"total"] - p[6,"bound"] ), nrow = 2, byrow = T), alternative = "greater" )$p.value, "frac" = p[i,"frac"] ) )
}) )
sig$sig <- p.adjust(sig$sig, method = "BH"); sig$lab <- ifelse( sig$sig <= 0.05, "*", "" ); sig$cpm.c <- factor(sig$cpm.c, levels = 5:1)

r <- 58
g[[2]] <- ggplot(p, aes(x = cpm.c, y = frac, fill = cpm.c)) + geom_col(show.legend = F) + geom_col(data = q, aes(x = cpm.c, y = frac), fill = "transparent", color = "black", linewidth = 1/.75/.pt, show.legend = F) +
	scale_y_continuous(breaks = c(0,20,40), limits = c(0,r), labels = seq(0, 40, 20), expand = c(0,0)) + scale_fill_manual(values = rev(col.c)) + theme_classic() + 
	labs(x = "", y = "% bound by at least two TFs") + scale_x_discrete(labels = rep("", 6) ) + geom_text(aes(x = cpm.c, y = frac + (r * 0.045), label = bound), size = 8/.pt, color = "black", angle = 270) + 
	geom_text(data = sig, aes(x = cpm.c, y = frac + (r * 0.09), label = lab), angle = 270, color = "dark orange", size = 8/.pt, fontface = "bold", inherit.aes = F) + geom_text(data = q, aes(x = cpm.c, y = frac - (r * 0.045), label = bound), size = 8/.pt, color = "black", angle = 270) + 
	theme(axis.title = element_text(size = 8, color = "black"), axis.text = element_text(size = 8, color = "black"), axis.ticks.length.x = unit(.04, "in"), axis.ticks.y = element_blank(), axis.ticks.x = element_line(linewidth = 1/.75/.pt, color = "black"), axis.line = element_line(linewidth = 1/.75/.pt, color = "black")) + 
	guides(fill = guide_legend(reverse=TRUE)) + coord_flip()

# Combine panels for easier arranging via cowplot
g2 <- plot_grid( plotlist = g, labels = c("B", "C"), ncol = 2, rel_widths = c(1,.85), label_size = 12 )

# Remove the objects you no longer need
rm(col.c, g, p, q, r, sig, tmp)


## For the paper, calculate how many of the ODEs are bound by at least 2 TFs
# length(which(rowSums(subset(a, repName == "ORR1E")[,7:12] >= 1) >= 2)) / nrow(subset(a, repName == "ORR1E")) * 100
# 30.3681
# 
# length(which(rowSums(subset(a, repName == "ORR1D2")[,7:12] >= 1) >= 2)) / nrow(subset(a, repName == "ORR1D2")) * 100
# 27.93462


## For the paper, calculate how many of the ODEs bound by at least 2 TFs are bound by PU.1 in progenitor cells
# sum(a[names(which(rowSums(subset(a, repName == "ORR1E")[,7:12] >= 1) >= 2)),7] > 0) / length(which(rowSums(subset(a, repName == "ORR1E")[,7:12] >= 1) >= 2)) * 100
# 89.77273
# 
# sum(a[names(which(rowSums(subset(a, repName == "ORR1D2")[,7:12] >= 1) >= 2)),7] > 0) / length(which(rowSums(subset(a, repName == "ORR1D2")[,7:12] >= 1) >= 2)) * 100
# 83.24468



## D 

### NOTE: The majority of panel D visualization happened in illustrator (to visualize clades) and unix (FastTree to generate phylogeny) and
### all relevant code will be found in the file "code_final"

# Load in the list of ODEs and ODE-derived peaks
ode <- read.table("rObject_ode.txt", sep = "\t", header = T); ode$cpm.c <- factor(ode$cpm.c, levels = 1:8); ode$peak <- gsub(".*\\|", "", ode$pid)
rownames(ode) <- ode$id

# For each clade, you want to figure out how many of each ATAC-seq cluster there are to make the piecharts. Get started!
a <- read.table("r1d2_iroki_metadata.txt", header = T); a$cpm.c <- ode[gsub("_", ":", a$name), "cpm.c"]
p <- data.frame(clade = 1:max(as.numeric(gsub("k_", "", a$branch_color))))
for (f in 1:6) {
	p[,LETTERS[f]] <- sapply(1:nrow(p), function(i){nrow(subset(a, branch_color == paste0("k_", i) & cpm.c == f))})
}
p$tot <- sapply(1:nrow(p), function(i){sum(p[i,2:ncol(p)])})

# Create the various themes you're thinking of using
col.c <- c(tableau_color_pal(type = "ordered-diverging", palette = "Sunset-Sunrise Diverging")(7)[1:5], "gray69") 

pdf("sfig3d_phyloClades_piecharts.pdf", height = 2.5, width = 5)
ggplot() + geom_scatterpie(aes(x = clade*2, y = 3, r = 0.75), data = p, cols = LETTERS[1:6], show.legend = F) + coord_fixed() + 
	scale_fill_manual(values = col.c) + theme_void() + geom_text(data = p, aes(x = clade*2, y = 4, label = tot), color = "black", size = 8/.pt) + 
	geom_text(data = p, aes(x = clade*2, y = 2, label = clade), color = "black", size = 8/.pt) + ylim(0,5)
dev.off()

# 
# geom_text(data = subset(p, clade > 6), aes(x = (clade-6)*2, y = 2, label = tot), fontface = "bold", size = 8/.pt) + 
# geom_text(data = subset(p, clade < 7), aes(x = clade*2, y = 3, label = clade), fontface = "bold", size = 8/.pt) + geom_text(data = subset(p, clade > 6), aes(x = (clade-6)*2, y = 0, label = clade), fontface = "bold", size = 8/.pt) + 

# Include the proportion of all elements assigned to clusters to get a sense of the overall distribution
p1 <- data.frame(clade = "all")
for (f in 1:6) {
	p1[,LETTERS[f]] <- nrow(subset(a, cpm.c == f))
}
p1$tot <- sum(p1[,2:ncol(p1)])

pdf("sfig3d_all_piechart.pdf", height = 2, width = 4)
ggplot() + geom_scatterpie(aes(x = 1, y = 1, r = 0.75), data = p1, cols = LETTERS[1:6]) + coord_fixed() + 
	scale_fill_manual(values = col.c, labels = c("A" = "Stem\n& Prog", "B" = "B", "C" = "T", "D" = "NK", "E" = "MF\n& DC", "F" = "Bkgd"), name = "ATAC-seq\ncluster") + theme_void() +
	geom_text(data = p1, aes(x = 1, y = 2, label = tot), size = 8/.pt) + ylim(0,3) + theme(legend.text = element_text(size = 8, color = "black"))
dev.off()

# Remove unneccesary objects
rm(a, f, p, p1)




## E-I

# Load in the packages you'll need throughout
suppressPackageStartupMessages({ 
	library(Biostrings); library(data.table); library(tidyr); library(doParallel); library(dplyr); ; library(ggthemes); library(ggplot2)
	library(fastSave); library(patchwork)
}) 

tmp <- homerToPFMatrixList("homerFilt.motifs")
mot.thresh <- sapply(tmp, function(i){return(log(2^i@tags$log2cut))}); pfms <- sapply(tmp, function(i){return(i@profileMatrix)})
names(pfms) <- names(mot.thresh) <- sapply(tmp, function(i){return( gsub("/Homer", "", i@ID) )})
meta <- data.frame(row.names = sapply(tmp, function(i){return(i@ID)}), "log.thresh" = mot.thresh, "proto.thresh" = mot.thresh * .5, "mot.length" = sapply(tmp, function(i){return(ncol(i@profileMatrix))}))

# Load in the consensus sequences for all of the ORR1 and related subfamily consensus sequences
cons <- readDNAStringSet("/workdir/jdc397/1_currentWork/8_teCactus/dfam3.6_cons.fa"); cons <- cons[grepl("ORR1E|ORR1D2|ORR1F|ORR1B1#|ORR1B2#|MTD", names(cons)),]
names(cons) <- gsub("#.*", "", names(cons))

# Load in the object containing all of the motif locations along each consensus sequence
base.con.mot <- read.table("rObject_base.con.mot.txt", header = T, sep = "\t")

a <- c("ORR1B1" = 6, "ORR1B2" = 5, "ORR1D2" = 4, "ORR1E" = 3, "ORR1F" = 2, "MTD" = 1)
b <- data.frame(row.names = rev(c("ORR1B1", "ORR1B2", "ORR1D2", "ORR1E", "ORR1F", "MTD")), xmin = 1, xmax = width(cons[rev(c("ORR1B1", "ORR1B2", "ORR1D2", "ORR1E", "ORR1F", "MTD"))]))

# To deal with a separate E2A motif having the same shorthand name, remove it from the above object first
base.con.mot <- subset(base.con.mot, mot != "E2A(bHLH),near_PU.1/Bcell-PU.1-ChIP-Seq(GSE21512)")

# Generate plots for the following motifs mentioned in the paper
n <- c("RORgt", "E2A", "ETS:RUNX", "RUNX", "CTCF")

g.1 <- g.2 <- list() 
for ( f in n ) {
	d <- subset(base.con.mot, motif == f)
	tmp <- ggplot(b) + geom_rect(aes(xmin = xmin, xmax = xmax, ymin = (1:length(a))-0.05, ymax = (1:length(a))+0.05 )) +
		geom_segment(data = subset(d, qual == "proto"), aes(x = (start+end)/2, y = ifelse(strand == "+", a[seq]+0.05, a[seq]-0.05), yend = ifelse(strand == "+", a[seq]+0.4, a[seq]-0.4), color = score), linewidth = 1.4, inherit.aes = F, show.legend = F) +
		geom_point(data = subset(d, qual == "mot" & strand == "+"), aes(x = (start+end)/2, y = a[seq]+0.35, fill = score), shape = 25, size = 2, color = "black", stroke = 1/.5/.pt,, inherit.aes = F) +
		geom_point(data = subset(d, qual == "mot" & strand == "-"), aes(x = (start+end)/2, y = a[seq]-0.35, fill = score), shape = 24, size = 2, color = "black", stroke = 1/.5/.pt, inherit.aes = F) +
		theme_classic() + scale_y_continuous(breaks = 1:length(a), labels = rev(names(a))) + labs(y = "Consensuses", x = "Consensus length")  +
		scale_fill_gradient(name = paste0(f, "\nmotif score"), low = "black", high = "dark orange", limits = c(0,max(d$score))) + scale_color_gradient(low = "black", high = "dark orange", limits = c(0,max(d$score))) +
		guides(color = element_blank(), fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) + 
		theme(legend.key.height = unit(0.15, "in"), legend.key.width = unit(0.1, "in"), aspect.ratio = 1/2, text = element_text(size = 8, color = "black"), axis.ticks.length=unit(.04, "in"), axis.ticks = element_line(linewidth = 1/.75/.pt), axis.line = element_line(linewidth = 1/.75/.pt, color = "black"), axis.text = element_text(size = 8, color = "black"), legend.text = element_text(size = 8, color = "black") )

	if ( f %in% c("RORgt", "E2A", "ETS:RUNX") ) {
		g.1[[ length(g.1) + 1 ]] <- tmp
	} else {
		g.2[[ length(g.2) + 1 ]] <- tmp
	}
}
rm(tmp)

# Based on current paneling, you'll have a group of three vertical plots and the remaining two horizontal.  
g4 <- plot_grid( plotlist = g.1, labels = c("E", "F", "G"), ncol = 1, label_size = 12 )
g5 <- plot_grid( plotlist = g.2, labels = c("H", "I"), ncol = 2, label_size = 12 )

# Finally, plot all the R-based panels. The heatmaps are generated using DeepTools and needs additional cleanup using illustrator 
pdf("sfig3acei_cowplot.pdf", height = 8.5, width = 7)
plot_grid( plot_grid(g1, g2, labels = c("A", ""), label_size = 12), plot_grid( NULL, g4 ), g5, rel_heights = c(1, 1.75, 0.75), ncol = 1, label_size = 12 )
dev.off()

# Remove all the objects you don't need
rm(a, b, base.con.mot, cons, d, f, g.1, g.2, meta, mot.thresh, n, pfms)


# Code for *all* the different motif plots
# 
# pdf("base.con.mot.plots.pdf", height = 8, width = 7)
# for ( f in unique(base.con.mot$motif) ) {
# 	d <- subset(base.con.mot, motif == f)
# 	print( ggplot(b) + geom_rect(aes(xmin = xmin, xmax = xmax, ymin = (1:length(a))-0.05, ymax = (1:length(a))+0.05 )) +
# 		geom_segment(data = subset(d, qual == "proto"), aes(x = (start+end)/2, y = ifelse(strand == "+", a[seq]+0.05, a[seq]-0.05), yend = ifelse(strand == "+", a[seq]+0.4, a[seq]-0.4), color = score), linewidth = 1.4, inherit.aes = F, show.legend = F) +
# 		geom_point(data = subset(d, qual == "mot" & strand == "+"), aes(x = (start+end)/2, y = a[seq]+0.35, fill = score), shape = 25, size = 2, color = "black", stroke = 1/.5/.pt,, inherit.aes = F) +
# 		geom_point(data = subset(d, qual == "mot" & strand == "-"), aes(x = (start+end)/2, y = a[seq]-0.35, fill = score), shape = 24, size = 2, color = "black", stroke = 1/.5/.pt, inherit.aes = F) +
# 		theme_classic() + scale_y_continuous(breaks = 1:length(a), labels = rev(names(a))) + labs(y = "Consensus sequences", x = "Consensus length")  +
# 		scale_fill_gradient(name = paste0(f, "\nmotif score"), low = "black", high = "dark orange", limits = c(0,max(d$score))) + scale_color_gradient(low = "black", high = "dark orange", limits = c(0,max(d$score))) +
# 		guides(color = element_blank(), fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(0.1, "in"), barheight = unit(.5, "in"))) + 
# 		theme(aspect.ratio = 1/2, text = element_text(size = 8, color = "black"), axis.ticks.length=unit(.04, "in"), axis.ticks = element_line(linewidth = 1/.75/.pt), axis.line = element_line(linewidth = 1/.75/.pt, color = "black"), axis.text = element_text(size = 8, color = "black"), legend.text = element_text(size = 8, color = "black") ) )
# }
# dev.off()



##############################
### Supplementary Figure 4 ###
##############################

# export PATH=/home/jdc397/lbzip2-2.5/bin:$PATH; module load R/4.1.3-r9; R --vanilla

# Load all the relevant libraries you need for this section
suppressPackageStartupMessages({ 
	library(data.table); library(tidyr); library(dplyr); library(ggplot2); library(doParallel); library(scatterpie); library(fastSave)
	library(Biostrings); library(patchwork)
}) 

## A & B

# Load the look-up table for ODEs and their respective consensus sequences
load.lbzip2(file = "rObject_lut.RDataFS", n.cores = 16)

# Load in the python-generated HOMER motif clusters
mot.clust <- read.table("../../../homer_mot_cluster_python.txt", header = T, sep = "\t", row.names = 1)

# Load in the list of ODEs and ODE-derived peaks
ode <- read.table("cpm_heatmap_ode_anno.txt", sep = "\t", header = T); ode$cpm.c <- factor(ode$cpm.c, levels = 1:8); ode$peak <- gsub(".*\\|", "", ode$pid)
rownames(ode) <- ode$id

# Load in the results of the HOMER motif enrichment
df.mot.enrich <- read.table(file = "homer/rObject_df.mot.enrich.txt", sep = "\t", header = T)

# Load in the sequence of the consensus sequences you (sometimes) care about
cons <- readDNAStringSet("/workdir/jdc397/1_currentWork/8_teCactus/dfam3.6_cons.fa") # ; cons <- cons[grepl("ORR1E|ORR1D2|ORR1F|ORR1B1#|ORR1B2#|MTD", names(cons)),]
cons <- cons[grepl("ORR1E|ORR1D2", names(cons)),]; names(cons) <- gsub("#.*", "", names(cons))

# Load in *all* motifs found across ODE loci and filter for the ones above
load.lbzip2(file = "rObject_ode.mot.RDataFS", n.cores = 16); ode.mot <- subset(ode.mot, repName %in% names(cons))

# Load in the list of motifs and proto-motifs for ODE consensuses and filter for the ones above
base.con.mot <- subset(read.table(file = "rObject_base.con.mot.txt", sep = "\t", header = T), seq %in% names(cons))

# Generate the table of ODEs per CPM cluster
tmp.tab <- c( table( unique(subset(ode, repName == "ORR1E")[,c("id","cpm.c")])$cpm.c), table(subset(ode, repName == "ORR1D2")$cpm.c)[1:6] )

# Generate a temporary file that you use as the base for the below loop
con.tmp <- as.data.frame(matrix( -1, nrow = length(cons), ncol = max(width(cons)) ), row.names = names(cons))
for (f in 1:nrow(con.tmp)) {
	con.tmp[ rownames(con.tmp)[f], 1:(width(cons)[f]) ] <- 0
}
rm(f)

# Now for the big boy method! Iterate through each mot along ORR1E first and figure out the proportions as mentioned. Woot woot
it <- sort(unique(subset(ode.mot, !is.na(mot.clust))$motif.uniq))
ode.anc.base <- bind_rows( mclapply(it, function(f) {
	tmp <- con.tmp; mc <- mot.clust[f,"cluster"]

	# Start with the proto-motifs (if they exist!
	s <- base.con.mot[base.con.mot$mot.clust == mc & base.con.mot$qual == "proto",]
	if ( nrow(s) != 0 ) {
		for (g in 1:nrow(s)) {
			tmp[ s[g,"seq"], s[g,"start"]:s[g,"end"] ] <- 1
		}
	}

	# Now for the consensus motif matches (if they exist!)
	s <- base.con.mot[base.con.mot$mot.clust == mc & base.con.mot$qual == "mot",]
	if ( nrow(s) != 0 ) {
		for (g in 1:nrow(s)) {
			tmp[ s[g,"seq"], s[g,"start"]:s[g,"end"] ] <- 2
		}
	}

	### NOTE: To make things easier for ggscatterpie, you followed their naming scheme for columns: "ancestral" = "A", "gained" = "B", "scratch" = "C" 

	# With the temporary consensus object created, it's time to calculate the proportion for each cluster of elements! Good luck
	tmp.ode <- data.frame("repName" = rep(c("ORR1E", "ORR1D2"), c(8,6)), "motif" = f, "cpm.c" = c(1:8, 1:6), "total" = 0, "A" = 0, "B" = 0, "C" = 0)
	rownames(tmp.ode) <- with(tmp.ode, paste0(repName, "_", cpm.c)) 

	for (h in 1:8) {
		# iterate through each cluster, assuming there are hits!
		s <- subset(ode.mot, motif.uniq == f & cpm.c == h) 
		if ( nrow(s) != 0 ) {
			for (g in 1:nrow(s)) {
				# First, get the motif location in the element in the respective cluster and find the consensus value calculated above at those
				# positions
				a <- as.numeric(lut[ s[g,"V1"], s[g,"start"]:s[g,"end"] ]); l <- length(a); a <- a[!is.na(a)] 
				b <- paste0(s[g,"repName"], "_", s[g,"cpm.c"]); d <- as.numeric(tmp[ s[g,"repName"], a ])

				# Now that those values are acquired within "d", check to see if at least half of those base pairs contain "2" (motif), "1" 
				# (proto), or either "0/-1" (novel)
				if ( sum(d == 2) >= floor(l/2) ) {
					tmp.ode[b,"A"] <- tmp.ode[b,"A"] + 1
				} else if ( sum(d == 1) >= floor(l/2) ) {
					tmp.ode[b,"B"] <- tmp.ode[b,"B"] + 1
				} else {
					tmp.ode[b,"C"] <- tmp.ode[b,"C"] + 1
				}
				tmp.ode[b,"total"] <- tmp.ode[b,"total"] + 1
			}

		}
	}

	# Now move over the results you've gathered to the combined object to save the results and also plot them!
	# tmp.ode$num.total <- c( table(subset(ode, repName == "ORR1E")$cpm.c), table(subset(ode, repName == "ORR1D2")$cpm.c)[1:6] ); rownames(tmp.ode) <- NULL
	tmp.ode$num.total <- tmp.tab; rownames(tmp.ode) <- NULL

	return( tmp.ode )
}, mc.cores = 32) )

# With all that fancy analysis done, save the above object to disk so it doesn't need to be regenerated and make those plots!
write.table(ode.anc.base, file = "rObject_ode.anc.base.txt", sep = "\t", quote = F, row.names = F)


# Otherwise, load in the above table
# ode.anc.base <- read.table("rObject_ode.anc.base.txt", sep = "\t", header = T)

ode.anc.base$cpm.c <- factor(ode.anc.base$cpm.c, levels = (1:8))

# Add in the fraction of elements that have the motif from the HOMER enrichment analyses
ode.anc.base$uniq.num <- 0; ode.anc.base$sig <- "not"
for (f in 1:nrow(ode.anc.base)) {
	d <- ode.anc.base[f, "cpm.c"]
	if ( ( ode.anc.base[f, "cpm.c"] == 8 & ode.anc.base[f, "repName"] == "ORR1E" ) | ( ode.anc.base[f, "cpm.c"] == 6 & ode.anc.base[f, "repName"] == "ORR1D2" ) ) {
		ode.anc.base[f, "uniq.num"] <- subset(df.mot.enrich, motif.uniq == ode.anc.base[f, "motif"] & samp1 == tolower(ode.anc.base[f, "repName"]) & clust == 1)$V9
		ode.anc.base[f, "sig"] <- NA
	} else {
		ode.anc.base[f, "uniq.num"] <- subset(df.mot.enrich, motif.uniq == ode.anc.base[f, "motif"] & samp1 == tolower(ode.anc.base[f, "repName"]) & clust == d)$V7

		# Check whether the motif is enriched, depleted, or not significant relative to background 
		if ( nrow(subset(df.mot.enrich, motif.uniq == ode.anc.base[f, "motif"] & samp1 == tolower(ode.anc.base[f, "repName"]) & clust == d & V5 > 0 & padj <= 0.05) != 0) ) {
			ode.anc.base[f, "sig"] <- "enr"
		} else if ( nrow(subset(df.mot.enrich, motif.uniq == ode.anc.base[f, "motif"] & samp2 == tolower(ode.anc.base[f, "repName"]) & clust == d & V5 > 0 & padj <= 0.05) != 0) ) {
			ode.anc.base[f, "sig"] <- "dep"
		} 
	}
}
rm(f,d)

# Now get a list of all the significant TF motifs that are enriched in accessible ORR1s # & V7 >= 10
tf.list.1e <- unique(subset(df.mot.enrich, padj <= 5e-2 & V7 >= 10 & expr.all == T & samp1 == "orr1e")$motif.uniq)

# Get the list of the top 3 TFs (each given motif cluster - as per your version of Vierstra code - is limited to top factor)
a <- sapply(1:7, function(i) {
	b <- subset(df.mot.enrich, samp1 == "orr1e" & motif.uniq %in% tf.list.1e & clust == i & padj <= 0.05) %>% group_by(motif.clust) %>% arrange(padj) %>% slice_head(n = 1) %>% ungroup() %>% arrange(desc(V5)) %>% head(3) %>% as.data.frame()
	return( as.character(b$motif.uniq) )
})
b <- unique(as.character(unlist(a)))[c(-6,-15)][rev(c(13,3,2,8,1,5,6,4,7,11,10,9,12,14))]
d <- 14:1; names(d) <- b

p <- subset(ode.anc.base, motif %in% names(d) & repName == "ORR1E"); p$lab <- d[p$motif]
p$c.lab <- 9 - as.numeric(p$cpm.c)

g1 <- ggplot() + geom_scatterpie(aes( x = lab, y = c.lab, r = 0.45 ), data = p, cols = LETTERS[1:3]) + coord_fixed() + geom_point(data = p, aes(x = lab, y = c.lab, color = sig), size = 5.5) + geom_point(data = p, aes(x = lab, y = c.lab), size = 4, color = "white") + 
	geom_text(data = p, aes(x = lab, y = c.lab, label = round(uniq.num)), color = "black", fontface = "bold", size = 8/.pt) + scale_fill_manual(values = c("#33A02C", "#B2df8A", "#CAB2D6"), name = "Consensus\nmotif", labels = c("A" = "Ancestral", "B" = "Gained", "C" = "de novo")) + 
	scale_color_manual(name = "Significance", breaks = c("enr", "not", "dep"), labels = c("not" = "No", "enr" = "Enriched", "dep" = "Depleted"), values = c("not" = "black", "enr" = "orange", "dep" = "blue"), na.value = "white") + 
	scale_x_continuous(labels = rev(gsub("\\(.*", "", names(d))), breaks = 1:length(d), name = "", position = "top", guide = guide_axis(angle = 45)) + 
	scale_y_continuous(labels = rev(c("Stem\n& Prog", "B", "T", "NK", "ILC3", "MF", "DC", "Bkgd")), breaks = 1:8, name = "", position = "right") +
	theme_classic() + guides(color = guide_legend(override.aes = list(size = 2))) + theme(text = element_text(size = 8, color = "black"), axis.text = element_text(size = 8, color = "black"), axis.ticks.length=unit(.04, "in"), axis.ticks = element_line(linewidth = 1/.75/.pt), legend.text = element_text(size = 8, color = "black") )


# Now get a list of all the significant TF motifs that are enriched in accessible ORR1D2s
tf.list.1d2 <- unique(subset(df.mot.enrich, padj <= 5e-2 & V7 >= 10 & expr.all == T & samp1 == "orr1d2")$motif.uniq)

# Get the list of the top 3 TFs (each given motif cluster - as per your version of Vierstra code - is limited to top factor)
a <- sapply(1:5, function(i) {
	b <- subset(df.mot.enrich, samp1 == "orr1d2" & motif.uniq %in% tf.list.1d2 & clust == i & padj <= 0.05) %>% group_by(motif.clust) %>% arrange(padj) %>% slice_head(n = 1) %>% ungroup() %>% arrange(desc(V5)) %>% head(3) %>% as.data.frame()
	return( as.character(b$motif.uniq) )
})
b <- unique(as.character(unlist(a)))[c(11,13,12,8,9,5,6,7,4,3,10,2)]
d <- 12:1; names(d) <- b

p <- subset(ode.anc.base, motif %in% names(d) & repName == "ORR1D2"); p$lab <- d[p$motif]
p$c.lab <- 7 - as.numeric(p$cpm.c)

g2 <- ggplot() + geom_scatterpie(aes( x = lab, y = c.lab, r = 0.45 ), data = p, cols = LETTERS[1:3]) + coord_fixed() + geom_point(data = p, aes(x = lab, y = c.lab, color = sig), size = 5.5) + geom_point(data = p, aes(x = lab, y = c.lab), size = 4, color = "white") + 
	geom_text(data = p, aes(x = lab, y = c.lab, label = round(uniq.num)), color = "black", fontface = "bold", size = 8/.pt) + scale_fill_manual(values = c("#33A02C", "#B2df8A", "#CAB2D6"), name = "Consensus\nmotif", labels = c("A" = "Ancestral", "B" = "Gained", "C" = "de novo")) + 
	scale_color_manual(name = "Significance", breaks = c("enr", "not", "dep"), labels = c("not" = "No", "enr" = "Enriched", "dep" = "Depleted"), values = c("not" = "black", "enr" = "orange", "dep" = "blue"), na.value = "white") + 
	scale_x_continuous(labels = rev(gsub("\\(.*", "", names(d))), breaks = 1:length(d), name = "", position = "top", guide = guide_axis(angle = 45)) + 
	scale_y_continuous(labels = rev(c("Stem\n& Prog", "B", "T", "NK", "MF\n& DC", "Bkgd")), breaks = 1:6, name = "", position = "right") +
	theme_classic() + guides(color = guide_legend(override.aes = list(size = 2)), fill = guide_legend(override.aes = list(size = 2))) + theme(text = element_text(size = 8, color = "black"), axis.text = element_text(size = 8, color = "black"), axis.ticks.length=unit(.04, "in"), axis.ticks = element_line(linewidth = 1/.75/.pt), legend.text = element_text(size = 8, color = "black") )

# Finally, plot the figures! Patchwork works better in only showing one set of legends, so use that here 
pdf("sfig4ab_patchwork.pdf", height = 7, width = 6.5)
g1 / g2 + plot_layout(guides = "collect") + plot_annotation(tag_levels = c("A", "B")) & theme(plot.tag = element_text(size = 12, face = "bold"))
dev.off()

# pdf("sfig4ab_cowplot.pdf", height = 7, width = 6.5)
# plot_grid(g1, g2, labels = c("A", "B"), ncol = 1, label_size = 12)
# dev.off()

rm(a, b, base.con.mot, con.tmp, cons, d, df.mot.enrich, it, lut, mot.clust, ode, ode.anc.base, ode.mot, p, tf.list.1d2, tf.list.1e, tmp.tab)



##############################
### Supplementary Figure 5 ###
##############################

# export PATH=/home/jdc397/lbzip2-2.5/bin:$PATH; module load R/4.1.3-r9; R --vanilla

# Load all the relevant libraries you need for this section
suppressPackageStartupMessages({ 
	library(data.table); setDTthreads(8); library(dplyr); library(fastSave); library(doParallel); library(ggplot2); library(cowplot)
	library(grid); library(fgsea); library(ggthemes)
})


## A

# Load in the full R correlation object "cor.mat" using the "fastSave" library
load.lbzip2(file = "rObject_cor.mat.RDataFS", n.cores = 32)

g1 <- ggplot() + geom_violin(data = cor.mat, aes(x = cor.s, y = 3)) + geom_violin(data = subset(cor.mat, abs(dist) <= 1e6), aes(x = cor.s, y = 2)) + 
	geom_violin(data = subset(cor.mat, filt.s == T), aes(x = cor.s, y = 1)) + scale_y_continuous(breaks = c(1,2,3), labels = c("filtered", "within 1Mb", "all"), name = "CRE-gene correlation categories") + 
	scale_x_continuous(limits = c(-1,1), name = "Spearman correlation coefficient") + theme_classic() + 
	theme(text = element_text(size = 8, color = "black"), axis.title = element_text(size = 8, color = "black"), axis.text = element_text(size = 8, color = "black"), axis.ticks.y = element_blank(), axis.ticks.x = element_line(linewidth = unit(0.4, "in"), color = "black"))


## B

# Load in the subseted R correlation object "cor.mat.mil" using the "fastSave" library
load.lbzip2(file = "rObject_cor.mat.mil.RDataFS", n.cores = 16)

p <- as.data.frame(table(subset( cor.mat.mil, g1 != "orr1e" )$g2))$Freq; names(p) <- c(1:6, "nonte", "te")
meta <- rbind( data.frame("g" = factor(c(1:6, "nonte", "te"), levels = c(1:6, "nonte", "te")), 
	"frac" = as.data.frame(table(subset(cor.mat.mil, filt.s == T & g1 != "orr1e" & cor.s >= 0)$g2))$Freq / p, "num" = as.data.frame(table(subset(cor.mat.mil, filt.s == T & g1 != "orr1e" & cor.s >= 0)$g2))$Freq, "dir" = "pos"), 
	data.frame("g" = factor(c(1:6, "nonte", "te"), levels = c(1:6, "nonte", "te")), "frac" = as.data.frame(table(subset(cor.mat.mil, filt.s == T & g1 != "orr1e" & cor.s < 0)$g2))$Freq / p, "num" = as.data.frame(table(subset(cor.mat.mil, filt.s == T & g1 != "orr1e" & cor.s < 0)$g2))$Freq, "dir" = "neg"))
meta$dir <- factor(meta$dir, levels = c("neg", "pos"))

b <- nrow(subset(cor.mat.mil, g1 == "orr1d2" & g2 == 6 & filt.s == T & cor.s >= 0)); b2 <- nrow(subset(cor.mat.mil, g1 == "orr1d2" & g2 == 6 & filt.s == T & cor.s < 0))
sig <- rbind( data.frame("g" = factor(c(1:6, "nonte", "te"), levels = c(1:6, "nonte", "te")), "dir" = "pos", "rel" = sapply(c(1:6, "te", "nonte"), function(i){ ifelse( subset(meta, dir == "pos" & g == i)$frac >= subset(meta, dir == "pos" & g == 6)$frac, "up", "down" )}), "sig" = c( sapply( 1:5, function(i){ a <- nrow(subset(cor.mat.mil, g1 == "orr1d2" & g2 == i & filt.s == T & cor.s >= 0)); fisher.test( matrix(c( a, p[i] - a, b, p[6] - b ), nrow = 2, byrow = T), alternative = "two.sided" )$p.value } ), NA, sapply( c("te", "nonte"), function(i){ a <- nrow(subset(cor.mat.mil, g2 == i & filt.s == T & cor.s >= 0)); fisher.test( matrix(c( a, p[i] - a, b, p[6] - b ), nrow = 2, byrow = T), alternative = "two.sided" )$p.value } )) ),
	data.frame("g" = factor(c(1:6, "nonte", "te"), levels = c(1:6, "nonte", "te")), "dir" = "neg", "rel" = sapply(c(1:6, "te", "nonte"), function(i){ ifelse( subset(meta, dir == "neg" & g == i)$frac >= subset(meta, dir == "neg" & g == 6)$frac, "up", "down" )}), "sig" = c( sapply( 1:5, function(i){ a <- nrow(subset(cor.mat.mil, g1 == "orr1d2" & g2 == i & filt.s == T & cor.s < 0)); fisher.test( matrix(c( a, p[i] - a, b2, p[6] - b2 ), nrow = 2, byrow = T), alternative = "two.sided" )$p.value } ), NA, sapply( c("te", "nonte"), function(i){ a <- nrow(subset(cor.mat.mil, g2 == i & filt.s == T & cor.s >= 0)); fisher.test( matrix(c( a, p[i] - a, b, p[6] - b ), nrow = 2, byrow = T), alternative = "two.sided" )$p.value } )) ))
sig$sig <- p.adjust(sig$sig, method = "BH")
sig$lab <- ifelse( sig$sig <= 0.05, "*", "" ); sig$color <- ifelse( sig$rel == "up", "do", "b" ); sig$dir <- factor(sig$dir, levels = c("neg", "pos"))

tmp <- ggplot(meta, aes(x = frac * 100, y = rev(g), fill = dir)) + geom_col(position = position_dodge(), width = .6) +
	scale_fill_manual(values = c("pos" = "dark orange", "neg" = "black"), labels = c("pos" = "Pos", "neg" = "Neg"), name = "Parity") + 
	scale_x_continuous(labels = seq(0, 2.5, 0.5), breaks = seq(0, 2.5, 0.5), expand = c(0,0)) +
	geom_text(aes(x = (frac * 100) + 0.02, y = rev(g), label = num), hjust = "left", position = position_dodge(.65), size = 8/.pt) + theme_classic() +
	coord_cartesian(xlim = c(0,2.55), clip = "off") + geom_text(data = sig, aes(x = -0.1, y = rev(g), label = lab, color = color), hjust = "left", fontface = "bold", position = position_dodge(.65), size = 8/.pt, vjust = 0.8, show.legend = F) +
	annotate(geom = "text", x = -.15, y = rev(sig$g[1:8]), label = c("Stem\n& Prog", "B", "T", "NK", "MF\n& DC", "Bkgd", "nonTE", "TE"), size = 8/.pt, hjust = "right") + 
	scale_color_manual(values = c("do" = "dark orange", "b" = "blue")) + labs(x = "Filtered CRE-gene pairs (%)", y = "ORR1D2 ATAC-seq clusters") +
	guides( fill = guide_legend(reverse = T) ) + scale_y_discrete(labels = rep("", 8)) +
	theme( text = element_text(size = 8, color = "black"), axis.title.y = element_text(size = 8, color = "black", margin = margin(r = 35)), axis.text = element_text(size = 8, color = "black"), legend.key.size = unit(.08, 'in'), axis.ticks.length = unit(.04, "in"), axis.ticks.y = element_blank(), axis.ticks.x = element_line(linewidth = 1/.75/.pt, color = "black"), legend.text = element_text(size = 8, color = "black") )

g2 <- plot_grid( tmp, labels = "B", label_size = 12 )

rm( b, b2, meta, p, sig, tmp, cor.mat.mil )


## C & D


# Load in the subseted R correlation object "cor.mat.mil" using the "fastSave" library
load.lbzip2(file = "rObject_cor.mat.mil.RDataFS", n.cores = 16)

# Now do total regulatory interactions (pos + neg) and compare against non-TE and TE together
cor.mat.mil$g2 <- ifelse( cor.mat.mil$g2 %in% c("nonte", "te"), "comb", cor.mat.mil$g2 )

p <- as.data.frame(table(subset(cor.mat.mil, g1 != "orr1d2" )$g2))$Freq; names(p) <- c(1:8, "comb")
meta <- rbind( data.frame("g" = factor(c(1:8, "comb"), levels = c(1:8, "comb")), 
	"frac" = as.data.frame(table(subset(cor.mat.mil, filt.s == T & g1 != "orr1d2")$g2))$Freq / p, 
	"num" = as.data.frame(table(subset(cor.mat.mil, filt.s == T & g1 != "orr1d2")$g2))$Freq) )

b <- meta["comb","num"]
sig <- rbind( data.frame("g" = factor(c(1:8, "comb"), levels = c(1:8, "comb")), 
	"rel" = sapply(c(1:8, "comb"), function(i){ ifelse( subset(meta, g == i)$frac >= subset(meta, g == "comb")$frac, "up", "down" )}), 
	"sig" = c(sapply( 1:8, function(i){ a <- meta[i,"num"]; fisher.test( matrix(c( a, p[i] - a, b, p["comb"] - b ), nrow = 2, byrow = T), alternative = "two.sided" )$p.value } ), NA)) )
sig$sig <- p.adjust(sig$sig, method = "BH")
sig$lab <- ifelse( sig$sig <= 0.05, "*", "" ); sig$lab[9] <- ""; sig$color <- ifelse( sig$rel == "up", "do", "b" )

# Create the color theme
col.c <- rev(c(tableau_color_pal(type = "ordered-diverging", palette = "Sunset-Sunrise Diverging")(7), "gray69", "gray57")); names(col.c) <- c(1:8, "comb")

# Make the plot!
g3 <- ggplot(meta, aes(x = frac * 100, y = rev(g), fill = rev(g))) + geom_col(width = .6, show.legend = F) +
	scale_fill_manual(values = rev(col.c)) + 
	scale_x_continuous(labels = seq(0, 3, 1), breaks = c(0,1,2,3), expand = c(0,0)) +
	geom_text(aes(x = (frac * 100) + 0.04, y = rev(g), label = num), hjust = "left", size = 8/.pt) + theme_classic() +
	coord_cartesian(xlim = c(0,4), clip = "off") + geom_text(data = sig, aes(x = -0.1, y = rev(g), label = lab, color = color), hjust = "left", fontface = "bold", size = 8/.pt, vjust = 0.8, show.legend = F) +
	annotate(geom = "text", x = -.15, y = rev(sig$g), label = c("Stem\n& Prog", "B", "T", "NK", "ILC3", "MF", "DC", "Bkgd", "All other"), size = 8/.pt, hjust = "right") + 
	scale_color_manual(values = c("do" = "dark orange", "b" = "blue")) + labs(x = "Filtered CRE-gene pairs (%)", y = "ORR1E ATAC-seq clusters") +
	guides( fill = guide_legend(reverse = T) ) + scale_y_discrete(labels = rep("", 9)) +
	theme( text = element_text(size = 8, color = "black"), axis.title.y = element_text(size = 8, color = "black", margin = margin(r = 35)), axis.text = element_text(size = 8, color = "black"), legend.key.size = unit(.08, 'in'), axis.ticks.length = unit(.04, "in"), axis.ticks.y = element_blank(), axis.ticks.x = element_line(linewidth = 1/.75/.pt, color = "black"), legend.text = element_text(size = 8, color = "black") )



p <- as.data.frame(table(subset(cor.mat.mil, g1 != "orr1e" )$g2))$Freq; names(p) <- c(1:6, "comb")
meta <- rbind( data.frame("g" = factor(c(1:6, "comb"), levels = c(1:6, "comb")), 
	"frac" = as.data.frame(table(subset(cor.mat.mil, filt.s == T & g1 != "orr1e")$g2))$Freq / p, 
	"num" = as.data.frame(table(subset(cor.mat.mil, filt.s == T & g1 != "orr1e")$g2))$Freq) )

b <- meta["comb","num"]
sig <- rbind( data.frame("g" = factor(c(1:6, "comb"), levels = c(1:6, "comb")), 
	"rel" = sapply(c(1:6, "comb"), function(i){ ifelse( subset(meta, g == i)$frac >= subset(meta, g == "comb")$frac, "up", "down" )}), 
	"sig" = c(sapply( 1:6, function(i){ a <- meta[i,"num"]; fisher.test( matrix(c( a, p[i] - a, b, p["comb"] - b ), nrow = 2, byrow = T), alternative = "two.sided" )$p.value } ), NA)) )
sig$sig <- p.adjust(sig$sig, method = "BH")
sig$lab <- ifelse( sig$sig <= 0.05, "*", "" ); sig$lab[7] <- ""; sig$color <- ifelse( sig$rel == "up", "do", "b" )

# Create the color theme
col.c <- rev(c(tableau_color_pal(type = "ordered-diverging", palette = "Sunset-Sunrise Diverging")(7)[1:5], "gray69", "gray57")); names(col.c) <- c(1:6, "comb")

# Make the plot for ORR1D2
g4 <- ggplot(meta, aes(x = frac * 100, y = rev(g), fill = rev(g))) + geom_col(width = .6, show.legend = F) +
	scale_fill_manual(values = rev(col.c)) + 
	scale_x_continuous(labels = seq(0, 4, 1), breaks = c(0,1,2,3,4), expand = c(0,0)) +
	geom_text(aes(x = (frac * 100) + 0.04, y = rev(g), label = num), hjust = "left", size = 8/.pt) + theme_classic() +
	coord_cartesian(xlim = c(0,4.1), clip = "off") + geom_text(data = sig, aes(x = -0.1, y = rev(g), label = lab, color = color), hjust = "left", fontface = "bold", size = 8/.pt, vjust = 0.8, show.legend = F) +
	annotate(geom = "text", x = -.15, y = rev(sig$g), label = c("Stem\n& Prog", "B", "T", "NK", "MF\n& DC", "Bkgd", "All other"), size = 8/.pt, hjust = "right") + 
	scale_color_manual(values = c("do" = "dark orange", "b" = "blue")) + labs(x = "Filtered CRE-gene pairs (%)", y = "ORR1D2 ATAC-seq clusters") +
	guides( fill = guide_legend(reverse = T) ) + scale_y_discrete(labels = rep("", 7)) +
	theme( text = element_text(size = 8, color = "black"), axis.title.y = element_text(size = 8, color = "black", margin = margin(r = 35)), axis.text = element_text(size = 8, color = "black"), legend.key.size = unit(.08, 'in'), axis.ticks.length = unit(.04, "in"), axis.ticks.y = element_blank(), axis.ticks.x = element_line(linewidth = 1/.75/.pt, color = "black"), legend.text = element_text(size = 8, color = "black") )


## E

# Load in the list of ODEs and ODE-derived peaks
ode <- read.table("rObject_ode.txt", sep = "\t", header = T); ode$cpm.c <- factor(ode$cpm.c, levels = 1:8); ode$peak <- gsub(".*\\|", "", ode$pid)
rownames(ode) <- ode$id

# Load in the subseted R correlation object "cor.mat.filt" using the "fastSave" library
load.lbzip2(file = "rObject_cor.mat.filt.RDataFS", n.cores = 16)

# Load in the list of TE-derived peak ids 
tes <- unique(read.table("peak.te.int")$V4)

# Load in all ImmGen ATAC-seq peaks
f.m <- as.data.frame(fread("summit_unif_peaks_counts.txt"))[,1:4]; colnames(f.m) <- c("chr", "start", "end", "id")
rownames(f.m) <- with(f.m, paste0(chr, ":", start, "-", end, "|", id)); f.m$class <- ifelse( f.m$id %in% gsub(".*\\|", "", ode$pid), "ode", ifelse( f.m$id %in% tes, "te", "nonte" ) )

# Create the color theme
col.c <- c(tableau_color_pal(type = "ordered-diverging", palette = "Sunset-Sunrise Diverging")(7), "gray69", "gray47", "gray25")

# Calculate the enrichment ratio for each cluster for each ODE. Formula is: (# of DNA-DNA contacts in group / total DNA-DNA contacts) / (total bp in group / genome size)
gs <- 2725537669
tc <- sum(cor.mat.filt$filt.abc == T)

b <- bind_rows( lapply(1:6, function(i) {
	return( data.frame( "g" = i, "total" = nrow(subset(cor.mat.filt, filt.abc == T & g1 == "orr1d2" & g2 == i)), "bp" = sum(with(subset(ode, cpm.c == i & repName == "ORR1D2"), end - start + 1)) ) )
}) )
b <- rbind( b, data.frame( "g" = "nonTE", "total" = nrow(subset(cor.mat.filt, filt.abc == T & g1 == "nonte")), "bp" = sum(with(subset(f.m, class == "nonte"), end - start + 1)) ) )
b <- rbind( b, data.frame( "g" = "TE", "total" = nrow(subset(cor.mat.filt, filt.abc == T & g1 == "te")), "bp" = sum(with(subset(f.m, class == "te"), end - start + 1)) ) )
b$frac <- with( b, (total / tc) / (bp / gs) )

b$g <- factor(b$g, levels = c("TE", "nonTE", 6:1))

r <- 82
g5 <- ggplot(b, aes(x = frac, y = g, fill = g)) + geom_col(show.legend = F) + scale_fill_manual(values = rev(col.c[c(1:5,8:10)])) + 
	geom_text(aes(x = frac + (r * .02), y = g, label = total), size = 8/.pt, color = "black", hjust = 0) + theme_classic() +
	scale_x_continuous(name = "T cell DNA-DNA contact\nenrichment score", limits = c(0,r), breaks = c(0,20,40,60,80), expand = c(0,0)) +
	scale_y_discrete(name = "ORR1D2 ATAC-seq clusters", labels = rev(c("Stem\n& Prog", "B", "T", "NK", "MF\n& DC", "Bkgd", "nonTE", "TE"))) + 
	theme( text = element_text(size = 8, color = "black"), axis.text = element_text(size = 8, color = "black"), axis.ticks.length = unit(.04, "in"), axis.ticks.y = element_blank(), axis.ticks = element_line(linewidth = 1/.75/.pt, color = "black"), axis.line = element_line(linewidth = 1/.75/.pt, color = "black") )

# g3 <- plot_grid( tmp, labels = "C", label_size = 12 )


## For the discussion, calculate the fraction of the genome taken up by ODEs with cell type-specific accessibility patterns overall

# Gather the genome sizes for accessible ORR1E loci
b.2 <- bind_rows( lapply(1:8, function(i) {
	return( data.frame( "g" = i, "total" = nrow(subset(cor.mat.filt, filt.abc == T & g1 == "orr1e" & g2 == i)), "bp" = sum(with(subset(ode, cpm.c == i & repName == "ORR1E"), end - start + 1)) ) )
}) )

# Add the lengths of non-background ORR1E and ORR1D2 loci and divide them by the mouse genome size
sum(b.2[1:7,"bp"], b[1:5,"bp"]) / gs * 100
# 0.03251443

rm( b, b.2, col.c, cor.mat.mil, f.m, gs, ode, tc, tes )



## F

# Load in the packages you'll need for dealing with *all* the correlations
suppressPackageStartupMessages({
	library(ggplot2); library(data.table); setDTthreads(8); library(fastSave); library(fgsea); library(doParallel); library(dplyr)
	library(cowplot); 
})

# Load in the R correlation object using the "fastSave" library
load.lbzip2(file = "rObject_cor.mat.filt.RDataFS", n.cores = 16)

# Set the random seed to apply to mclapply as well! 
RNGkind("L'Ecuyer-CMRG")

skip.streams <- function(n) {
  x <- .Random.seed
  for (i in seq_len(n))
    x <- nextRNGStream(x)
  assign('.Random.seed', x, pos=.GlobalEnv)
}

# Get a subset of all the relevant ode-associated peak-gene linkages
filt <- arrange(subset(cor.mat.filt, (filt.abc == T | filt.s == T) & abs(dist) > 1000 & g1 %in% c("orr1e", "orr1d2")), desc(abs(cor.s)))

# Load in the HSC pathways, but remove the neutrophil-derived groups since those aren't present in your cells
newGO <- gmtPathways("/workdir/jdc397/1_currentWork/8_teCactus/Konturek-Ciesla - 2023 - HSPC marker genes.txt")
newGO <- newGO[!grepl("Neu", names(newGO))]

set.seed(18)
skip.streams(8)
r1e.filt.gsea <- mclapply(1:8, function(f) {
	a <- subset(filt, g2 == f & g1 == "orr1e"); a <- group_by(a, gene) %>% arrange(desc(abs(cor.s))) %>% slice_head(n=1) %>% as.data.frame() %>% arrange(desc(cor.s))
	b <- a$cor.s; names(b) <- a$gene
	return( fgseaMultilevel(pathways = newGO, stats = b, minSize=5) %>% as.data.frame() %>% arrange(desc(NES)) )
}, mc.cores = 8)

skip.streams(6)
r1d2.filt.gsea <- mclapply(1:6, function(f) {
	a <- subset(filt, g2 == f & g1 == "orr1d2"); a <- group_by(a, gene) %>% arrange(desc(abs(cor.s))) %>% slice_head(n=1) %>% as.data.frame() %>% arrange(desc(cor.s))
	a[is.na(a$cor.s),"cor.s"] <- 0; b <- a$cor.s; names(b) <- a$gene
	return( fgseaMultilevel(pathways = newGO, stats = b, minSize=5) %>% as.data.frame() %>% arrange(desc(NES)) )
}, mc.cores = 6)


# Now for ORR1D2
tmp <- data.frame()
for (f in 1:length(r1d2.filt.gsea)) {
	if ( nrow( r1d2.filt.gsea[[f]] ) != 0 ) {
		temp <- r1d2.filt.gsea[[f]]; temp$g <- f
		tmp <- rbind(tmp, temp)
	}
}
tmp$g <- factor(tmp$g, levels = (6:1))

# Have size be tmp-value and color be NES
g6 <- ggplot(subset(tmp, padj <= 0.05), aes(x = pathway, y = g, size = -log10(padj), color = NES)) + geom_point() + theme_minimal_grid(line_size = 1/.75/.pt, color = "#EBEBEB") + 
	scale_y_discrete(name = "ORR1D2 ATAC-seq clusters", drop = F, labels = rev(c("Stem\n& Prog", "B", "T", "NK", "MF\n& DC", "Bkgd"))) +
	scale_color_gradient2(low = "blue", mid = "white", high = "dark orange") + scale_x_discrete(name = "", guide = guide_axis(angle = 90)) + scale_size(range = c(1,5)) +
	guides(fill = "none", color = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(0.15, "in"), barheight = unit(.5, "in"))) + theme( text = element_text(size = 8, color = "black"), axis.text = element_text(size = 8, color = "black"), axis.ticks.length = unit(.04, "in"), axis.ticks = element_line(linewidth = 1/.75/.pt, color = "#EBEBEB") )

pdf("sfig5af.pdf", height = 9, width = 7)
plot_grid( plot_grid( g1, g2, ncol = 2, labels = c("A", ""), label_size = 12 ), plot_grid( g3, g4, ncol = 2, labels = c("C", "D"), label_size = 12 ), plot_grid( g5, g6, ncol = 2 , labels = c("E", "F"), label_size = 12 ), nrow = 3 )
dev.off()

rm( f, filt, newGO, r1d2.filt.gsea, r1e.filt.gsea, skip.streams, temp )

### NOTE: You manually added Fig. S5G in Illustrator since it came from IGV and was *much* easier to just add on later


###################################
### Supplementary Figures 6 & 7 ###
###################################

# export PATH=/home/jdc397/lbzip2-2.5/bin:$PATH; module load R/4.1.3-r9; R --vanilla

# Load all the relevant libraries you need for this section
suppressPackageStartupMessages({ 
	library(fastSave); library(ggplot2); library(ggthemes); library(cowplot); library(ggnewscale)
})

# Load in all of the assigned mouse vs human delta TPMs across each cell type
l <- list.files(pattern = "mouse_vs_human_orthology_*"); tmp <- data.frame()
for (f in l) {
	a <- read.table(f, sep = "\t", header = T); a$samp <- gsub(".txt", "", gsub(".*_", "", f))
	tmp <- rbind(tmp, a)
}
mvh <- tmp; rm(l, tmp, a)

# Load the look-up-table that contains the samples you'll use for each different cell type in mouse and human
load.lbzip2(file = "rObject_lut.hm.RDataFS", n.cores = 2)

# Get that color palette ready (with an extra value for all other genes!) 
col.c <- c(tableau_color_pal(type = "ordered-diverging", palette = "Sunset-Sunrise Diverging")(7), "gray69", "gray51") 

# Generate the distribution of positively correlated ODEs for each ATAC-seq cluster for each cell type to identify whether the target genes are significantly 
# more expressed in mouse than in human. Include all putatively orthologous genes! Also, don't include ORR1E GN and NK in the supplemental figure
# and instead have them by themselves with the numbers
g <- h <- list()
for ( f in rownames(lut.hm) ) {

	# Subset the object containing all the gene expression score information for the current cell type
	tmp <- subset(mvh, samp == f & type == "orth" & ((r1e.con.clust != "other" & r1e.spear >= 0) | r1e.con.clust == "other")); tmp$lab <- factor(tmp$r1e.con.clust, levels = c("other",8:1))

	# Generate the object and calculate stats for plotting 
	meta <- data.frame("lab" = factor(c(1:8,"other"), levels = c("other",8:1)), 
		"p.val" = c(p.adjust(sapply(1:8, function(i){wilcox.test(subset(tmp, r1e.con.clust == i)$delta.mvh, subset(tmp, r1e.con.clust == "other")$delta.mvh, alternative = "two.sided")$p.value}), method = "BH"), 1), 
		"delta.mvh" = rep(-1.9,9), 
		"num" = as.numeric(table(tmp$r1e.con.clust)),
		"dir" = c(ifelse( sapply(1:8, function(i){ median(subset(tmp, r1e.con.clust == i)$delta.mvh) }) <= median(subset(tmp, r1e.con.clust == "other")$delta.mvh), "less", "greater" ), NA) )
	meta$sig <- with(meta, ifelse(p.val < 0.0005, "***", ifelse(p.val < 0.005, "**", ifelse(p.val < 0.05, "*", "") ) ) )
	meta[meta$lab == "other","num"] <- paste0("~", round(meta[meta$lab == "other","num"]/1000,0), "k")

	# Check to see if the current value is "NK" or "GN" to add the plot or not
	if ( !(f %in% c("NK", "GN")) ) {

		g[[ length(g) + 1 ]] <- ggplot(tmp, aes(x = lab, y = delta.mvh, color = lab, fill = lab)) + theme_classic() + geom_hline(yintercept = median(subset(tmp, lab == "other")$delta.mvh), linetype = "31", linewidth = 1/.75/.pt) +
			geom_boxplot(width = 0.5, outlier.color = NA, fill = "white", show.legend = F, linewidth = 1/.75/.pt, fatten = 1) + scale_color_manual(values = col.c[9:1]) + 
			new_scale_color() + geom_text(data = meta, aes(x = lab, y = delta.mvh, label = sig, color = dir), size = 8/.pt, fontface = "bold", angle = 270, show.legend = F) +
			scale_color_manual(values = c("less" = "black", "greater" = "dark orange")) + 
			scale_x_discrete(labels = rev(c("Stem\n& Prog", "B", "T", "NK", "ILC3", "MF", "DC", "Bkgd", "All Other\nGenes"))) + coord_flip(ylim = c(-2.5,2.5)) + labs( x = "", y = expression(Delta~"TPM (M - H)") ) +
			ggtitle(paste0("ORR1E: ", f)) + theme(plot.title = element_text(size = 8, color = "black", hjust = 0.5), axis.text = element_text(size = 8, color = "black"), axis.title = element_text(size = 8), axis.ticks.length = unit(.04, "in"), axis.ticks = element_line(linewidth = 1/.75/.pt, color = "black"), axis.line = element_line(linewidth = 1/.75/.pt, color = "black"))
	}
}


# Now for ORR1D2 (same as above but with no comments and no removal of main plots, woot woot)
for ( f in rownames(lut.hm) ) {
	tmp <- subset(mvh, samp == f & type == "orth" & ((r1d2.con.clust != "other" & r1d2.spear >= 0) | r1d2.con.clust == "other"))
	tmp$lab <- factor(tmp$r1d2.con.clust, levels = c("other",6:1))
	meta <- data.frame("lab" = factor(c(1:6,"other"), levels = c("other",6:1)), 
		"p.val" = c(p.adjust(sapply(1:6, function(i){wilcox.test(subset(tmp, r1d2.con.clust == i)$delta.mvh, subset(tmp, r1d2.con.clust == "other")$delta.mvh, alternative = "two.sided")$p.value}), method = "BH"), 1), 
		"delta.mvh" = rep(-1.9,7), 
		"num" = as.numeric(table(tmp$r1d2.con.clust)),
		"dir" = c(ifelse( sapply(1:6, function(i){ median(subset(tmp, r1d2.con.clust == i)$delta.mvh) }) <= median(subset(tmp, r1d2.con.clust == "other")$delta.mvh), "less", "greater" ), NA) )
	meta$sig <- with(meta, ifelse(p.val < 0.0005, "***", ifelse(p.val < 0.005, "**", ifelse(p.val < 0.05, "*", "") ) ) )
	meta[meta$lab == "other","num"] <- paste0("~", round(meta[meta$lab == "other","num"]/1000,0), "k")

	g[[ length(g) + 1 ]] <- ggplot(tmp, aes(x = lab, y = delta.mvh, color = lab, fill = lab)) + theme_classic() + geom_hline(yintercept = median(subset(tmp, lab == "other")$delta.mvh), linetype = "31", linewidth = 1/.75/.pt) +
		geom_boxplot(width = 0.5, outlier.color = NA, fill = "white", show.legend = F, linewidth = 1/.75/.pt, fatten = 1) + scale_color_manual(values = col.c[c(9:8,5:1)]) + 
		new_scale_color() + geom_text(data = meta, aes(x = lab, y = delta.mvh, label = sig, color = dir), size = 8/.pt, fontface = "bold", angle = 270, show.legend = F) +
		scale_color_manual(values = c("less" = "black", "greater" = "dark orange")) + 
		scale_x_discrete(labels = rev(c("Stem\n& Prog", "B", "T", "NK", "MF\n& DC", "Bkgd", "All Other\nGenes"))) + coord_flip(ylim = c(-2.5,2.5)) + labs( x = "", y = expression(Delta~"TPM (M - H)") ) +
		ggtitle(paste0("ORR1D2: ", f)) + theme(plot.title = element_text(size = 8, color = "black"), axis.text = element_text(size = 8, color = "black"), axis.title = element_text(size = 8), axis.ticks.length = unit(.04, "in"), axis.ticks = element_line(linewidth = 1/.75/.pt, color = "black"), axis.line = element_line(linewidth = 1/.75/.pt, color = "black"))
}

### NOTE: For the sfig 6 legend, ORR1E numbers are: 
# Stem & Prog: 83, B: 81, T: 174, NK: 146, ILC3: 22, MF: 147, DC: 79, Bkgd: 38, Other: 15609
# 
### ORR1D2:
# Stem & Prog: 54, B: 61, T: 83, NK: 94, MF & DC: 115, Bkgd: 22, Other: 16155

# Plot the supplemental panels for sfig 6
pdf("sfig6av_deltaTPM_orthPosCorMvH_boxplots.pdf", height = 12, width = 8) 
plot_grid( plotlist = g, ncol = 5, nrow = 5, labels = "AUTO", label_size = 12 )
dev.off()


# Now, do the above but for sfig 7, which includes all of the negatively correlated linkages!  

g <- h <- list()
for ( f in rownames(lut.hm) ) {
	tmp <- subset(mvh, samp == f & type == "orth" & ((r1e.con.clust != "other" & r1e.spear < 0) | r1e.con.clust == "other")); tmp$lab <- factor(tmp$r1e.con.clust, levels = c("other",8:1))
	meta <- data.frame("lab" = factor(c(1:8,"other"), levels = c("other",8:1)), 
		"p.val" = c(p.adjust(sapply(1:8, function(i){wilcox.test(subset(tmp, r1e.con.clust == i)$delta.mvh, subset(tmp, r1e.con.clust == "other")$delta.mvh, alternative = "two.sided")$p.value}), method = "BH"), 1), 
		"delta.mvh" = rep(-1.9,9), 
		"num" = as.numeric(table(tmp$r1e.con.clust)),
		"dir" = c(ifelse( sapply(1:8, function(i){ median(subset(tmp, r1e.con.clust == i)$delta.mvh) }) <= median(subset(tmp, r1e.con.clust == "other")$delta.mvh), "less", "greater" ), NA) )
	meta$sig <- with(meta, ifelse(p.val < 0.0005, "***", ifelse(p.val < 0.005, "**", ifelse(p.val < 0.05, "*", "") ) ) )
	meta[meta$lab == "other","num"] <- paste0("~", round(meta[meta$lab == "other","num"]/1000,0), "k")

	g[[ length(g) + 1 ]] <- ggplot(tmp, aes(x = lab, y = delta.mvh, color = lab, fill = lab)) + theme_classic() + geom_hline(yintercept = median(subset(tmp, lab == "other")$delta.mvh), linetype = "31", linewidth = 1/.75/.pt) +
		geom_boxplot(width = 0.5, outlier.color = NA, fill = "white", show.legend = F, linewidth = 1/.75/.pt, fatten = 1) + scale_color_manual(values = col.c[9:1]) + 
		new_scale_color() + geom_text(data = meta, aes(x = lab, y = delta.mvh, label = sig, color = dir), size = 8/.pt, fontface = "bold", angle = 270, show.legend = F) +
		scale_color_manual(values = c("less" = "black", "greater" = "dark orange")) + 
		scale_x_discrete(labels = rev(c("Stem\n& Prog", "B", "T", "NK", "ILC3", "MF", "DC", "Bkgd", "All Other\nGenes"))) + coord_flip(ylim = c(-2.5,2.5)) + labs( x = "", y = expression(Delta~"TPM (M - H)") ) +
		ggtitle(paste0("ORR1E: ", f)) + theme(plot.title = element_text(size = 8, color = "black", hjust = 0.5), axis.text = element_text(size = 8, color = "black"), axis.title = element_text(size = 8), axis.ticks.length = unit(.04, "in"), axis.ticks = element_line(linewidth = 1/.75/.pt, color = "black"), axis.line = element_line(linewidth = 1/.75/.pt, color = "black"))

}

# Now for ORR1D2
for ( f in rownames(lut.hm) ) {
	tmp <- subset(mvh, samp == f & type == "orth" & ((r1d2.con.clust != "other" & r1d2.spear < 0) | r1d2.con.clust == "other"))
	tmp$lab <- factor(tmp$r1d2.con.clust, levels = c("other",6:1))
	meta <- data.frame("lab" = factor(c(1:6,"other"), levels = c("other",6:1)), 
		"p.val" = c(p.adjust(sapply(1:6, function(i){wilcox.test(subset(tmp, r1d2.con.clust == i)$delta.mvh, subset(tmp, r1d2.con.clust == "other")$delta.mvh, alternative = "two.sided")$p.value}), method = "BH"), 1), 
		"delta.mvh" = rep(-1.9,7), 
		"num" = as.numeric(table(tmp$r1d2.con.clust)),
		"dir" = c(ifelse( sapply(1:6, function(i){ median(subset(tmp, r1d2.con.clust == i)$delta.mvh) }) <= median(subset(tmp, r1d2.con.clust == "other")$delta.mvh), "less", "greater" ), NA) )
	meta$sig <- with(meta, ifelse(p.val < 0.0005, "***", ifelse(p.val < 0.005, "**", ifelse(p.val < 0.05, "*", "") ) ) )
	meta[meta$lab == "other","num"] <- paste0("~", round(meta[meta$lab == "other","num"]/1000,0), "k")

	g[[ length(g) + 1 ]] <- ggplot(tmp, aes(x = lab, y = delta.mvh, color = lab, fill = lab)) + theme_classic() + geom_hline(yintercept = median(subset(tmp, lab == "other")$delta.mvh), linetype = "31", linewidth = 1/.75/.pt) +
		geom_boxplot(width = 0.5, outlier.color = NA, fill = "white", show.legend = F, linewidth = 1/.75/.pt, fatten = 1) + scale_color_manual(values = col.c[c(9:8,5:1)]) + 
		new_scale_color() + geom_text(data = meta, aes(x = lab, y = delta.mvh, label = sig, color = dir), size = 8/.pt, fontface = "bold", angle = 270, show.legend = F) +
		scale_color_manual(values = c("less" = "black", "greater" = "dark orange")) + 
		scale_x_discrete(labels = rev(c("Stem\n& Prog", "B", "T", "NK", "MF\n& DC", "Bkgd", "All Other\nGenes"))) + coord_flip(ylim = c(-2.5,2.5)) + labs( x = "", y = expression(Delta~"TPM (M - H)") ) +
		ggtitle(paste0("ORR1D2: ", f)) + theme(plot.title = element_text(size = 8, color = "black"), axis.text = element_text(size = 8, color = "black"), axis.title = element_text(size = 8), axis.ticks.length = unit(.04, "in"), axis.ticks = element_line(linewidth = 1/.75/.pt, color = "black"), axis.line = element_line(linewidth = 1/.75/.pt, color = "black"))
}

### NOTE: For the sfig 7 legend, ORR1E numbers are: 
# Stem & Prog: 32, B: 94, T: 149, NK: 96, ILC3: 16, MF: 98, DC: 40, Bkgd: 43, Other: 15609
# 
### ORR1D2:
# Stem & Prog: 17, B: 46, T: 103, NK: 82, MF & DC: 84, Bkgd: 31, Other: 16155

# Plot the supplemental panels for sfig 7
pdf("sfig7ax_deltaTPM_orthNegCorMvH_boxplots.pdf", height = 12, width = 8) 
plot_grid( plotlist = g, ncol = 5, nrow = 5, labels = "AUTO", label_size = 12 )
dev.off()



##############################
### Supplementary Figure 8 ###
##############################

# export PATH=/home/jdc397/lbzip2-2.5/bin:$PATH; module load R/4.1.3-r9; R --vanilla

# Load in the packages you'll need throughout this section
suppressPackageStartupMessages({ 
	library(fastSave); library(orthogene); library(data.table); setDTthreads(8); library(doParallel); library(dplyr); library(fgsea)
	library(cowplot); library(ggplot2); library(ggthemes); library(scatterpie)
}) 


## A & C

# Load the look-up-table that contains the samples you'll use for each different cell type in mouse and human
load.lbzip2(file = "rObject_lut.hm.RDataFS", n.cores = 2)

# Load the GSEA results between mouse and human
load.lbzip2(file = "rObject_mouseVsHumanOrthologyGsea_res.RDataFS", n.cores = 16)

# Subset the orthologous/mouse-specific results for ORR1D2
tmp <- subset(res, pathway == "combMouseUp" & ode == "orr1d2")
tmp$g <- factor(tmp$cluster, levels = (6:1)); tmp$cell.type <- factor(tmp$cell.type, levels = c("nCD4", "aCD4", "nCD8", "aCD8", "gdT", "NK", "B", "", "C.Mo", "NC.Mo", "GN", "mDC", "pDC"))

# Generate the GSEA for ORR1D2
g1.m <- ggplot(subset(tmp, padj > 0.05 & padj <= 0.25), aes(x = cell.type, y = g, size = -log10(padj), color = NES)) + geom_point(pch = 1, stroke = 1/.5/.pt) + 
	geom_point(data = subset(tmp, padj <= 0.05), aes(x = cell.type, y = g, size = -log10(padj), color = NES)) + theme_minimal_grid(line_size = 1/.75/.pt) + 
	scale_y_discrete(drop = F, labels = rev(c("Stem\n& Prog", "B", "T", "NK", "MF\n& DC", "Bkgd"))) + scale_color_gradient2(low = "blue", mid = "white", high = "dark orange") + 
	scale_x_discrete(guide = guide_axis(angle = 90), breaks = rownames(lut.hm), drop = F) + scale_size(range = c(1,4)) +
        ylab("") + xlab("Cell Type Expression Data") + ggtitle("Mouse-specific genes") + guides(fill = "none", size = guide_legend(override.aes = list(stroke = 1/.75/.pt)), color = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(0.09, "in"), barheight = unit(.66, "in"))) + 
	theme(legend.title = element_text(size = 8, color = "black"), legend.text = element_text(size = 8, color = "black"), plot.title = element_text(size = 8, color = "black", hjust = 0.5), axis.text = element_text(size = 8, color = "black"), axis.title = element_text(size = 8, color = "black"))

tmp <- subset(res, pathway == "orthMouseUp_pos" & ode == "orr1d2")
tmp$g <- factor(tmp$cluster, levels = (6:1)); tmp$cell.type <- factor(tmp$cell.type, levels = c("nCD4", "aCD4", "nCD8", "aCD8", "gdT", "NK", "B", "", "C.Mo", "NC.Mo", "GN", "mDC", "pDC"))

g1.o <- ggplot(subset(tmp, padj > 0.05 & padj <= 0.25), aes(x = cell.type, y = g, size = -log10(padj), color = NES)) + geom_point(pch = 1, stroke = 1/.5/.pt) + 
	geom_point(data = subset(tmp, padj <= 0.05), aes(x = cell.type, y = g, size = -log10(padj), color = NES)) + theme_minimal_grid(line_size = 1/.75/.pt) + 
	scale_y_discrete(drop = F, labels = rev(c("Stem\n& Prog", "B", "T", "NK", "MF\n& DC", "Bkgd"))) + scale_color_gradient2(low = "blue", mid = "white", high = "dark orange") + 
	scale_x_discrete(guide = guide_axis(angle = 90), breaks = rownames(lut.hm), drop = F) + scale_size(range = c(1,4)) +
	ylab("") + xlab("Cell Type Expression Data") + ggtitle("Orthologous genes") + guides(fill = "none", size = guide_legend(override.aes = list(stroke = 1/.75/.pt)), color = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(0.09, "in"), barheight = unit(.66, "in"))) +
	theme(legend.title = element_text(size = 8, color = "black"), legend.text = element_text(size = 8, color = "black"), plot.title = element_text(size = 8, color = "black", hjust = 0.5), axis.text = element_text(size = 8, color = "black"), axis.title = element_text(size = 8, color = "black"))

# Save the plot of the above to manually add them to the supplement
pdf("sfig8ac_mouseVsHumanMouseUpGseaSigEnrich_dotplot.pdf", height = 2.9, width = 3.25)
g1.o; g1.m
dev.off()



## B

# Load in the filtered R correlation object using the "fastSave" library
load.lbzip2(file = "rObject_cor.mat.filt.RDataFS", n.cores = 16)

# Focus on ORR1D2 here, so assign all ORR1E links to the "te" category
tmp <- cor.mat.filt; tmp[ tmp$g1 == "orr1e", "g1"] <- tmp[ tmp$g1 == "orr1e", "g2"] <- "te"

# Arbitrarilly load in one of the previous objects that contains orthology information for gene targets (any cell type would do)
b <- fread(file = "prior_orth_noProteinCodingFilter/mouse_vs_human_orthology_nCD8.txt"); setkey(b, mouse.name)

# Use the orthology information saved in the above file to add that to the filtered correlation file per gene
tmp$type <- unlist(mclapply(1:nrow(tmp), function(i){ arrange(b[J(tmp[i,"gene"])], desc(abs(delta.mvh)))[,type][1] }, mc.cores = 64))
tmp$type <- ifelse( grepl("orth", tmp$type) | tmp$type == "comp.mouse", "orth", tmp$type )

# Generate the table of *unique* ODE cluster / peak class and synteny status (not including parity of the correlation)
tmp.filt <- unique(subset(tmp, !is.na(type))[,c("gene","g2","type")])
a <- as.data.frame(with(tmp.filt, table(g2, type)))

# Add the genomic distribution of mouse-specific to orthologous genes to the above object for comparison (calculated from above "b" file)
sm <- sum(subset(b, !is.na(mouse.name))$type == "mouse"); so <- sum(subset(b, !is.na(mouse.name))$type != "mouse")
a <- rbind( a, data.frame("g2" = rep("genome", 2), "type" = c("mouse", "orth"), "Freq" = c(sm,so) ) )

# Create the "frac" column by either dividing by the total size of each "g2", or split by "g2" and "Var3" which is the parity of the correlation
a$frac.all <- a$Freq / c( rep( table( tmp.filt$g2), 2 ), rep(sum(sm,so), 2) ) * 100

# Figure out (via fisher's test) if there are fewer mouse-specific ODE-gene links (out of protein-coding genes!) than the genome proportion
b <- table( tmp.filt$g2 ); b[ length(b) + 1 ] <- sum(sm,so); names(b)[length(b)] <- "genome"
sig <- data.frame( "g2" = factor(c(1:6,"nonte","te"), levels = c("genome", "te", "nonte", 6:1)), "sig" = p.adjust( sapply( c(1:6,"nonte","te"), function(i){ a1 <- subset(a, g2 == i & type == "mouse")$Freq; fisher.test( matrix(c( a1, b[i] - a1, sm, so ), nrow = 2), alternative = "less" )$p.value }), method = "BH") )
sig$lab <- ifelse( sig$sig <= 0.05, "*", "" )

# Assign factor levels and generate the plot!
a$g2 <- factor(a$g2, levels = c("genome", "te", "nonte", 6:1)); a$type <- factor(a$type, levels = c("orth", "mouse"))
pdf("sfig8b_orr1d2GeneTargetType_barplot.pdf", width = 1.75, height = 3)
ggplot(a, aes(x = frac.all, y = g2, fill = type)) + geom_vline(xintercept = subset(a, g2 == "genome" & type == "mouse")$frac.all, linetype = "31", linewidth = 1/0.75/.pt) + geom_col(width = 0.5) + scale_fill_manual(values = c("#000000", "#FF7518"), breaks = c("orth", "mouse"), labels = c("orth" = "Orthologous", "mouse" = "Mouse-Specific")) +
	scale_y_discrete(name = "", labels = NULL) + scale_x_continuous(name = "Filtered interactions (%)", labels = c("0","25","50","75","100"), expand = c(0,0)) + theme_classic() +
	coord_cartesian(clip = "off", xlim = c(0,100)) + annotate(geom = "text", x = -9, y = levels(a$g2), label = rev(c("Stem\n& Prog", "B", "T", "NK", "MF\n& DC", "Bkgd", "nonTE", "TE", "Genome")), size = 8/.pt, hjust = "right") +
	annotate(geom = "text", x = -2, y = rev(levels(a$g2)), label = c(sig$lab, ""), fontface = "bold", size = 8/.pt, hjust = "right", vjust = .75) +
	theme( legend.position = "top", legend.text = element_text(size = 8, color = "#000000"), legend.title = element_text(size = 8, color = "#000000"), axis.text = element_text(size = 8, color = "black"), axis.title.y = element_text(size = 8, color = "black", margin = margin(r = 35)), axis.title.x = element_text(size = 8, color = "#000000"), axis.ticks.length = unit(.04, "in"), axis.ticks.y = element_blank(), axis.ticks.x = element_line(linewidth = 1/.75/.pt, color = "black"), axis.line = element_line(linewidth = 1/.75/.pt, color = "black") ) +
	labs(fill = "")
dev.off()




## D

# Load all the relevant libraries you need for this section
suppressPackageStartupMessages({ 
	library(Biostrings); library(data.table); library(tidyr); library(doParallel); library(dplyr); library(fastSave); library(ggthemes); library(ggplot2)
	library(PWMEnrich); library(TFBSTools); library(monaLisa); library(patchwork)
}) 

# Load in the list of ODEs and ODE-derived peaks
ode <- read.table("cpm_heatmap_ode_anno.txt", sep = "\t", header = T); ode$cpm.c <- factor(ode$cpm.c, levels = 1:8); ode$peak <- gsub(".*\\|", "", ode$pid)
rownames(ode) <- ode$id

# Load the R object containing the phastCons scores across ODEs
load.lbzip2(file = "rObject_pc.out.RDataFS", n.cores = 32)

# Turn all of the "-1" values into NAs
pc.out[pc.out == -1] <- NA

# Calculate the average conservation across each element and append it since that's generally helpful to have in one place
pc.avg <- data.frame("repName" = gsub(".*\\|", "", rownames(pc.out)), "avg" = rowMeans(pc.out, na.rm = T), "num" = rowSums(!is.na(pc.out)))

# Get the clusters for each ODE
pc.avg$cpm.c <- NA
for (f in 1:nrow(pc.avg)) {
	if ( rownames(pc.avg)[f] %in% rownames(ode) ) {
		pc.avg[f,"cpm.c"] <- ode[ rownames(pc.avg)[f], "cpm.c" ] 
	}
}

# Read in the full-length non-accessible ODEs to use as a "neutral" background
fl <- read.table("fl_ode.txt")

# Get the colors you want for the below plot
col.c <- c(tableau_color_pal(type = "ordered-diverging", palette = "Sunset-Sunrise Diverging")(7), "gray69", "gray58")

## Now for ORR1D2!
s <- subset(pc.avg, grepl("1D2", rownames(pc.avg)) & !is.na(cpm.c))
a <- s; a$repName <- paste0(a$repName, "_", a$cpm.c); a <- rbind(a, subset(pc.avg, repName == "ORR1D2" & rownames(pc.avg) %in% fl$V4)) 
a$repName <- factor(a$repName, levels = c("ORR1D2",paste0("ORR1D2_", 6:1)))
meta <- data.frame("repName" = factor( levels(a$repName), levels = levels(a$repName) ), 
	"p.val" = c(NA, p.adjust(sapply( levels(a$repName)[-1], function(i){wilcox.test(subset(a, repName == i)$avg, subset(a, repName == "ORR1D2")$avg, alternative = "two.sided")$p.value}), method = "BH")), 
	"num" = as.data.frame(table(a$repName))$Freq, 
	"y" = c(NA, sapply( levels(a$repName)[-1], function(i){max(subset(a, repName == i)$avg)+0.01})) )
meta$sig <- with(meta, ifelse(p.val < 0.0005, "***", ifelse(p.val < 0.005, "**", ifelse(p.val < 0.05, "*", "") ) ) )
meta$dir <- c(NA, sapply( levels(a$repName)[-1], function(i){ if (median(subset(a, repName == i)$avg) - median(subset(a, repName == "ORR1D2")$avg) > 0) {return("up")} else {return("down")} }))
meta[1,"num"] <- paste0("~", round(meta[1,"num"]/1000,0), "k")

g4 <- ggplot(a, aes(x = repName, y = avg, fill = repName)) + geom_violin(linewidth = 1/.75/.pt, show.legend = F) + geom_hline(yintercept = median(subset(a, repName == "ORR1D2")$avg), lty = "31", lwd = 1/.75/.pt) + 
	geom_boxplot(width = 0.2, outlier.shape = NA, fill = "white", fatten = 1/.75/.pt) + scale_fill_manual(values = col.c[c(9:8,5:1)]) + theme_classic() + scale_x_discrete(name = "", labels = rev(c("Stem\n& Prog", "B", "T", "NK", "MF\n& DC", "Bkgd", "Full Length\nORR1D2"))) +
	labs(y = "Average phastCons score") + geom_text(data = meta, aes(x = repName, y = -0.04, label = num), size = 8/.pt) + geom_text(data = meta, aes(x = repName, y = y, label = sig, color = dir), size = 8/.pt, angle = 270, hjust = 0.25, fontface = "bold", show.legend = F) + 
	theme(axis.text = element_text(size = 8, color = "black"), axis.title = element_text(size = 8, color = "black"), axis.ticks.length = unit(.04, "in"), axis.ticks = element_line(linewidth = 1/.75/.pt, color = "black"), axis.line = element_line(linewidth = 1/.75/.pt, color = "black")) +
	ylim(c(-0.05,1.02)) + coord_flip() + scale_color_manual(values = c("up" = "dark orange", "down" = "black"))


## E

# Load the object containing all the syntenic ODEs and their motif enrichments in rat
df.mot.enrich <- read.table(file = "rObject_syntenicRatODE_df.mot.enrich.txt", sep = "\t", header = T)

## Now for ORR1D2
mtp <- c("ERG", "RUNX1", "EBF1", "ZEB1", "Elk1", "GATA", "Maz", "KLF3", "Tbet", "KLF6", "Smad3", "PU.1")
p <- subset(df.mot.enrich, padj <= 0.05 & motif %in% mtp & (samp1 == "orr1d2" | samp2 == "orr1d2")) %>% group_by(motif,clust) %>% arrange(padj) %>% slice_head(n = 1) %>% as.data.frame()
p$prat <- -log10(p$padj); p$lrat <- ifelse(p$samp1 == "orr1d2", p$V5, -p$V5)
p.lut <- data.frame(row.names = 1:5, "name" = c("Stem\n& Prog", "B", "T", "NK", "MF\n& DC")); p$lab <- factor(p.lut[p$clust,"name"], levels = rev(p.lut$name))
p$motif <- factor(p$motif, levels = mtp)

# Create fake data to make sure all the legend values are present in the final plot
f1 <- data.frame(x = factor("MF\n& DC", levels = rev(p.lut$name)), y = factor("PU.1", levels = levels(p$motif)), f = factor(c("Yes", "No"), levels = c("Yes", "No")))

g5 <- ggplot(p, aes(x = motif, y = lab, size = -log10(padj), fill = lrat)) + geom_point(data = f1, aes(x = y, y = x, color = f), pch = 21, alpha = 0, inherit.aes = F) + geom_point(data = subset(p, mouse == T), pch = 21, stroke = 1/.pt/.5, color = "black") + geom_point(data = subset(p, mouse == F), pch = 21, color = "transparent") + theme_bw() + 
	scale_fill_gradient2(low = "blue", mid = "white", high = "dark orange", limits = c(-1.5,1.5), oob = squish, breaks = c(-1,0,1), labels = c("-1", "0", "1"), name = "log2 ratio") + 
	scale_size_area(limits = c(-log10(0.05),3), oob = squish, max_size = 5, breaks = c(-log10(0.05), 2, 3), labels = c("5x10-2", "1x10-2", "<1x10-3"), name = "Adjusted p-value") +
	scale_color_manual(name = "Enrichment\nin mouse?", labels = c("Yes", ""), values = c("black", "transparent")) + scale_y_discrete(limits = rev(p.lut$name)) +
	theme(panel.border = element_blank(), text = element_text(size = 8, color = "black"), axis.text = element_text(size = 8, color = "black"), axis.title = element_text(size = 8, color = "black"), axis.ticks = element_blank(), legend.title = element_text(size = 8, color = "black"), legend.text = element_text(size = 8, color = "black"), legend.position = "bottom" ) + labs(x = "", y = "") +
	labs(color = "Expressed") + guides(fill = guide_colorbar(order = 1, title.position="top", title.hjust = 0.5, barwidth = unit(0.5, "in"), barheight = unit(.1, "in")), size = guide_legend(override.aes = list(pch = 21, fill = "black", color = "transparent"), order = 2, title.position="top", title.hjust = 0.5), color = guide_legend(override.aes = list(pch = 21, fill = "white", size = 5, color = c("black", "transparent"), stroke = 1/.pt/.5, alpha = 1), order = 3, title.position="top", title.hjust = 0.5)) +
	coord_fixed(clip = "off") + scale_x_discrete(position = "top", guide = guide_axis(angle = 45), drop = F)



####################################
### Supplementary Figures 9 & 10 ###
####################################

# export PATH=/home/jdc397/lbzip2-2.5/bin:$PATH; module load R/4.1.3-r9; R --vanilla

# Load all the relevant libraries you need for this section
suppressPackageStartupMessages({ 
	library(Biostrings); library(data.table); library(tidyr); library(doParallel); library(dplyr); library(fastSave); library(ggthemes); library(ggplot2)
	library(PWMEnrich); library(TFBSTools); library(monaLisa); library(patchwork)
}) 

# Set the number of cores to use for speed
setDTthreads(8); registerCoresPWMEnrich(16)

# Load in the consensus sequences and combined ODE .fa files
cons <- readDNAStringSet("/workdir/jdc397/1_currentWork/8_teCactus/dfam3.6_cons.fa"); cons <- cons[grepl("ORR1E|ORR1D2", names(cons)),]
comb <- c(readDNAStringSet("ode.fa"), readDNAStringSet("fl_ode.fa"))

# Load in the lookup-table with all the ODEs and the full-length background ones as well
load.lbzip2(file = "rObject_lut.all.RDataFS", n.cores = 16)

tmp <- homerToPFMatrixList("homerFilt.motifs")
mot.thresh <- sapply(tmp, function(i){return(log(2^i@tags$log2cut))}); pfms <- sapply(tmp, function(i){return(i@profileMatrix)})
names(pfms) <- names(mot.thresh) <- sapply(tmp, function(i){return( gsub("/Homer", "", i@ID) )})

### NOTE: This works for everything >0, except for OCT:OCT-short which has a long motif and a negative (!) score threshold. Gross.
meta <- data.frame(row.names = sapply(tmp, function(i){return(i@ID)}), "log.thresh" = mot.thresh, "proto.thresh" = mot.thresh * .5, "mot.length" = sapply(tmp, function(i){return(ncol(i@profileMatrix))}))

load.lbzip2(file = "rObject_scores.RDataFS", n.cores = 16)

# load.lbzip2(file = "rObject_base.mot.RDataFS", n.cores = 16)
base.con.mot <- read.table("rObject_base.con.mot.txt", header = T, sep = "\t")

# Load in the python-generated HOMER motif clusters
mot.clust <- read.table("../../../homer_mot_cluster_python.txt", header = T, sep = "\t", row.names = 1)

# Regenerate the object that *only* includes the ORR1E motifs that you care about so that you can easily make the figure. 
bmc <- data.frame()
for (f in c(33,131,79,114,8,75,129,85,12,11,36)) {
	s <- base.con.mot[base.con.mot$mot %in% rownames(mot.clust[mot.clust$cluster == f,,drop=F]),]
	s <- subset(s, qual != "no" & seq == "ORR1E")
	tmp <- as.data.frame(matrix("",nrow=1,ncol=359))
	for (g in which(s$qual == "proto")) { tmp[1,s[g,"start"]:s[g,"end"]] <- s[g,"qual"] }
	for (g in which(s$qual == "mot")) { tmp[1,s[g,"start"]:s[g,"end"]] <- s[g,"qual"] }

	# Now iterate through the columns and output a row that matches the stretch of a given character. Yeah!
	cur <- ""; st <- 1; w <- rev(c(33,131,79,114,8,75,129,85,12,11,36))
	for (g in 1:ncol(tmp)) {
		if ( tmp[1,g] != cur ) { # The current stretch has ended. Append to the main object and start another one
			bmc <- rbind(bmc, data.frame("seq" = "ORR1E", "mot.clust" = f, "qual" = cur, "xmin" = st - 0.5, "xmax" = g - 0.5, ymin = which(w == f)*.5, ymax = which(w == f)*.5 ) )
			st <- g; cur <- tmp[1,g]
		} 
	}
}
bmc <- subset(bmc, qual != ""); rm(s,tmp,f,g,cur,st,w)

# Generate the base ggplot for the motif ancestral and proto-motif locations 
bmc$xmin <- ifelse( bmc$xmin <= 1, 1, bmc$xmin )
bmc.p <- ggplot(bmc, aes(x=xmin, xend=xmax, y=ymin, yend=ymax, color = qual)) + geom_segment(linewidth = 3/.75/.pt) + scale_y_continuous(breaks = c(1:11)*.5, labels = rev(c("ETS", "KLF", "RUNX", "EBF1", "E2A", "CTCF", "Maz", "ZFX", "Reverb", "ROR", "PU.1:IRF"))) + theme_void() +
	coord_cartesian(xlim = c(1,359)) + scale_x_continuous(labels = rep("", 5), breaks = c(1,100,200,300,359)) + scale_color_manual(name = "Motifs in\nConsensus", labels = c("mot" = "Ancestral", "proto" = "Proto"), values = c("firebrick2", "dodgerblue")) +
	theme(axis.text = element_text(size = 8, color = "black"), legend.title = element_text(size = 8, color = "black"), legend.text = element_text(size = 8, color = "black"), aspect.ratio = 1/2)

# Load in the list of ODEs and ODE-derived peaks
ode <- read.table("rObject_ode.txt", sep = "\t", header = T)[,c(1:5,7)]; ode$cpm.c <- factor(ode$cpm.c, levels = 1:8); rownames(ode) <- ode$id
fl <- read.table("fl_ode.txt"); colnames(fl)[1:4] <- c("chr", "start", "end", "id"); fl$repName <- gsub(".*\\|", "", fl$id); fl$cpm.c <- "fl"
ode <- rbind(ode, fl[,c(1:3,7,4,8)]); ode$cpm.c <- factor(ode$cpm.c, levels = c(1:8,"fl")) 

# Load in the object containing all the phastCons values for each ODE cluster (including full-length background)
load.lbzip2(file = "rObject_pc.con.out.RDataFS", n.cores = 16)

# Turn all of the "-1" values into NAs
pc.con.out[pc.con.out == -1] <- NA

p <- data.frame(); p2 <- data.frame()
for ( f in c(1:8,"fl") ) {
	a <- unique(subset(ode, grepl("ORR1E", id) & cpm.c == f)$id)
	
	for (g in 1:359) {
		ci <- confint(lm(pc.con.out[a,g] ~ 1), level = 0.95)
		p <- rbind(p, data.frame( "x" = g, "ymin" = ci[1], "ymax" = ci[2], "repName" = "orr1e", "group" = f ))
	}

	ks <- ksmooth(1:359, colMeans(pc.con.out[a,], na.rm = T), "normal", bandwidth = 10)
	p2 <- rbind(p2, data.frame("x" = ks$x, "y" = ks$y, "repName" = "orr1e", "group" = f))

	if ( f < 7 | f == "fl" ) {
		a <- unique(subset(ode, grepl("ORR1D2", id) & cpm.c == f)$id)
		for (g in 1:370) {
			ci <- confint(lm(pc.con.out[a,g] ~ 1), level = 0.95)
			p <- rbind(p, data.frame( "x" = g, "ymin" = ci[1], "ymax" = ci[2], "repName" = "orr1d2", "group" = f ))
		}

		ks <- ksmooth(1:370, colMeans(pc.con.out[a,], na.rm = T), "normal", bandwidth = 10)
		p2 <- rbind(p2, data.frame("x" = ks$x, "y" = ks$y, "repName" = "orr1d2", "group" = f))
	}
}
p$group <- factor(p$group, levels = c(1:8,"fl")); p2$group <- factor(p2$group, levels = c(1:8,"fl"))

# Get the colors you want for the below plot
col.c <- c(tableau_color_pal(type = "ordered-diverging", palette = "Sunset-Sunrise Diverging")(7), "gray69", "gray58")

# Create a named vector to use for renaming the cluster numbers into something meaningful! 
lut.vec <- c("1" = "Stem\n& Prog", "2" = "B", "3" = "T", "4" = "NK", "5" = "ILC3", "6" = "MF", "7" = "DC", "8" = "Bkgd", "fl" = "Full Length\nBkgd")

# Generate all of the plots for ORR1E you don't already have in main Figure 6
g <- list()
for (f in c(2:5,7:8)) {

	peaks <- zoo::rollapply( zoo::as.zoo(subset(p2, repName == "orr1e" & group == f)$y), 5, function(x) which.max(x)==2)

	a <- bmc.p + geom_vline(xintercept = which(peaks == T)+1, linetype = "31", color = "black", linewidth = 0.75/.75/.pt)

	p1 <- ggplot(subset(p, group %in% c(f,"fl") & repName == "orr1e")) + geom_ribbon(aes(x = x, group = group, ymin = ymin, ymax = ymax, fill = group), alpha = 0.5) + 
		geom_line(data = subset(p2, group %in% c(f,"fl") & repName == "orr1e"), aes(x = x, y = y, group = group, color = group), linewidth = 1/.75/.pt, alpha = 0.85, inherit.aes = F) + 
		scale_y_continuous(breaks = c(0.05,0.165,0.28), limits = c(min(subset(p, repName == "orr1e")$ymin), max(subset(p, repName == "orr1e")$ymax))) + scale_x_continuous(breaks = c(1,100,200,300,359)) + 
		geom_vline(xintercept = which(peaks == T)+1, linetype = "31", color = "black", linewidth = 0.75/.75/.pt) + scale_fill_manual(values = col.c[c(f,9)], name = "ORR1E Cluster", labels = lut.vec[c(f,9)]) + 
		scale_color_manual(values = col.c[c(f,9)], name = "ORR1E Cluster", labels = lut.vec[c(f,9)]) + theme_classic() + labs(y = "Average phastCons score", x = "Position in ORR1E Consensus") + 
		theme(aspect.ratio = 1/2, axis.title.y = element_text(size = 8, angle = 90, vjust = 10, color = "black"), axis.title.x = element_text(size = 8, color = "black"), axis.text = element_text(size = 8, color = "black"), legend.title = element_text(size = 8, color = "black"), legend.text = element_text(size = 8, color = "black"), axis.ticks.length = unit(.04, "in"), axis.ticks = element_line(linewidth = 1/.75/.pt, color = "black"), axis.line = element_line(linewidth = 1/.75/.pt, color = "black"))

	g[[ length(g) + 1 ]] <- a / plot_spacer() / p1 + plot_layout(heights = c(2,-.5,2), guides = "collect")
}

pdf("sfig9af_phastConsEnrichedMotif_combplot.pdf", height = 9.75, width = 8.6)
plot_grid( plotlist = g, ncol = 2, nrow = 3, labels = "AUTO", label_size = 12 )
dev.off()


# Now do the same but for ORR1D2

# Regenerate the object that *only* includes the ORR1D2 motifs that you care about so that you can easily make the figure. 
bmc <- data.frame()
for (f in c(33,79,114,6,54,129,131,127,97)) {
	s <- base.con.mot[base.con.mot$mot %in% rownames(mot.clust[mot.clust$cluster == f,,drop=F]),]
	s <- subset(s, qual != "no" & seq == "ORR1D2")
	tmp <- as.data.frame(matrix("",nrow=1,ncol=359))
	for (g in which(s$qual == "proto")) { tmp[1,s[g,"start"]:s[g,"end"]] <- s[g,"qual"] }
	for (g in which(s$qual == "mot")) { tmp[1,s[g,"start"]:s[g,"end"]] <- s[g,"qual"] }

	# Now iterate through the columns and output a row that matches the stretch of a given character. Yeah!
	cur <- ""; st <- 1; w <- rev(c(33,79,114,6,54,129,131,127,97))
	for (g in 1:ncol(tmp)) {
		if ( tmp[1,g] != cur ) { # The current stretch has ended. Append to the main object and start another one
			bmc <- rbind(bmc, data.frame("seq" = "ORR1D2", "mot.clust" = f, "qual" = cur, "xmin" = st - 0.5, "xmax" = g - 0.5, ymin = which(w == f)*.5, ymax = which(w == f)*.5 ) )
			st <- g; cur <- tmp[1,g]
		} 
	}
}
bmc <- subset(bmc, qual != ""); rm(s,tmp,f,g,cur,st,w)

# Generate the base ggplot for the motif ancestral and proto-motif locations 
bmc$xmin <- ifelse( bmc$xmin <= 1, 1, bmc$xmin )
bmc.p <- ggplot(bmc, aes(x=xmin, xend=xmax, y=ymin, yend=ymax, color = qual)) + geom_segment(linewidth = 3/.75/.pt) + scale_y_continuous(breaks = c(1:9)*.5, labels = rev(c("ETS", "RUNX", "EBF1", "ZEB1", "GATA", "Maz", "KLF", "Tbet", "Smad3"))) + theme_void() +
	coord_cartesian(xlim = c(1,370)) + scale_x_continuous(labels = rep("", 5), breaks = c(1,100,200,300,370)) + scale_color_manual(name = "Motifs in\nConsensus", labels = c("mot" = "Ancestral", "proto" = "Proto"), values = c("firebrick2", "dodgerblue")) +
	theme(axis.text = element_text(size = 8, color = "black"), legend.title = element_text(size = 8, color = "black"), legend.text = element_text(size = 8, color = "black"), aspect.ratio = 1/2)

col.c <- c(tableau_color_pal(type = "ordered-diverging", palette = "Sunset-Sunrise Diverging")(7)[1:5], "gray69", "gray58")
lut.vec <- c("1" = "Stem\n& Prog", "2" = "B", "3" = "T", "4" = "NK", "5" = "MF\n& DC", "6" = "Bkgd", "fl" = "Full Length\nBkgd")

g <- list()
for (f in c(1:6)) {

	peaks <- zoo::rollapply( zoo::as.zoo(subset(p2, repName == "orr1d2" & group == f)$y), 5, function(x) which.max(x)==2)

	a <- bmc.p + geom_vline(xintercept = which(peaks == T)+1, linetype = "31", color = "black", linewidth = 0.75/.75/.pt)

	p1 <- ggplot(subset(p, group %in% c(f,"fl") & repName == "orr1d2")) + geom_ribbon(aes(x = x, group = group, ymin = ymin, ymax = ymax, fill = group), alpha = 0.5) + 
		geom_line(data = subset(p2, group %in% c(f,"fl") & repName == "orr1d2"), aes(x = x, y = y, group = group, color = group), linewidth = 1/.75/.pt, alpha = 0.85, inherit.aes = F) + 
		scale_y_continuous(breaks = c(0.07,0.19,0.31), limits = c(min(subset(p, repName == "orr1d2")$ymin), max(subset(p, repName == "orr1d2")$ymax))) + scale_x_continuous(breaks = c(1,100,200,300,370)) + 
		geom_vline(xintercept = which(peaks == T)+1, linetype = "31", color = "black", linewidth = 0.75/.75/.pt) + scale_fill_manual(values = col.c[c(f,7)], name = "ORR1D2 Cluster", labels = lut.vec[c(f,7)]) + 
		scale_color_manual(values = col.c[c(f,7)], name = "ORR1D2 Cluster", labels = lut.vec[c(f,7)]) + theme_classic() + labs(y = "Average phastCons score", x = "Position in ORR1D2 Consensus") + 
		theme(aspect.ratio = 1/2, axis.title.y = element_text(size = 8, angle = 90, vjust = 10, color = "black"), axis.title.x = element_text(size = 8, color = "black"), axis.text = element_text(size = 8, color = "black"), legend.title = element_text(size = 8, color = "black"), legend.text = element_text(size = 8, color = "black"), axis.ticks.length = unit(.04, "in"), axis.ticks = element_line(linewidth = 1/.75/.pt, color = "black"), axis.line = element_line(linewidth = 1/.75/.pt, color = "black"))

	g[[ length(g) + 1 ]] <- a / plot_spacer() / p1 + plot_layout(heights = c(2,-.5,2), guides = "collect")
}

pdf("sfig10af_phastConsEnrichedMotif_combplot.pdf", height = 9.75, width = 8.6)
plot_grid( plotlist = g, ncol = 2, nrow = 3, labels = "AUTO", label_size = 12 )
dev.off()












#################################
### Supplementary Table 2 & 3 ###
#################################

# export PATH=/home/jdc397/lbzip2-2.5/bin:$PATH; module load R/4.1.3-r9; R --vanilla

# Load all the relevant libraries you need for this section
suppressPackageStartupMessages({ 
	library(ComplexHeatmap); library(ggthemes); library(dplyr); library(ggplot2); library(cowplot)
}) 


tmp <- c("PU.1_prog", "Pax5_B", "Runx1_CD8+", "RORg_Th17", "PU.1_MF", "IRF8_CD8+DC")
a <- read.table("figure3_ode_int_chip_vals.txt", sep = "\t"); colnames(a) <- c("chr", "start", "end", "repName", "id", "cpm.c", tmp)
rownames(a) <- a$id

# Rows are clusters and columns are binding datasets
p.1e <- as.data.frame(matrix(0, nrow = 8, ncol = length(tmp)), row.names = 1:8); colnames(p.1e) <- tmp; p.1d2 <- p.1e[1:6,]
q.1e <- p.1e; q.1d2 <- p.1d2

for ( i in 1:8 ) {

	k <- nrow(subset(a, repName == "ORR1E" & cpm.c == i))
	for ( j in tmp ) {
		p.1e[ i,j ] <- sum( subset(a, repName == "ORR1E" & cpm.c == i)[ , j ] > 0 )
		q.1e[ i,j ] <- round( sum( subset(a, repName == "ORR1E" & cpm.c == i)[ , j ] > 0 ) / k * 100, 2 )
	}

	if ( i <= 6 ) {
		k <- nrow(subset(a, repName == "ORR1D2" & cpm.c == i))
		for ( j in tmp ) {
			p.1d2[ i,j ] <- sum( subset(a, repName == "ORR1D2" & cpm.c == i)[ , j ] > 0 )
			q.1d2[ i,j ] <- round( sum( subset(a, repName == "ORR1D2" & cpm.c == i)[ , j ] > 0 ) / k * 100, 2 )
		}
	}
}

rownames(p.1e) <- rownames(q.1e) <- c("Stem & Prog", "B", "T", "NK", "ILC3", "MF", "DC", "Bkgd")
rownames(p.1d2) <- rownames(q.1d2) <- c("Stem & Prog", "B", "T", "NK", "MF & DC", "Bkgd")

# write.table(p.1e, "sTable1a_orr1e_bound_num.txt", quote = F, sep = "\t"); write.table(q.1e, "sTable1b_orr1e_bound_frac.txt", quote = F, sep = "\t")
# write.table(p.1d2, "sTable2a_orr1d2_bound_num.txt", quote = F, sep = "\t"); write.table(q.1d2, "sTable2b_orr1d2_bound_frac.txt", quote = F, sep = "\t")

write.table(q.1e, "tableS2_orr1eChipFracBound.txt", quote = F, sep = "\t")
write.table(q.1d2, "tableS3_orr1d2ChipFracBound.txt", quote = F, sep = "\t")



##############################
### Supplementary Analyses ###
##############################

## Measuring enrichment of ABC-linkages in the ATAC-RNA correlation linkages

# export PATH=/home/jdc397/lbzip2-2.5/bin:$PATH; module load R/4.1.3-r9; R --vanilla

library(fastSave); library(regioneR); library(data.table); setDTthreads(8)

# First, take the filtered correlations and output all of the ATAC-RNA linkages as a .bed file
load.lbzip2(file = "rObject_cor.mat.filt.RDataFS", n.cores = 16)

# Output a file that contains each of the linkages with TSS coordinates first and then peak coordinates second
write.table( tidyr::separate( tidyr::separate( with(subset(cor.mat.filt, filt.s == T), data.frame("g2" = gene.coords, "g" = gene, "p" = peak.coords, "id" = peak) ), col = "p", into = c("chr", "start", "end"), sep = "[:-]" ), col = "g2", into = c("gChr", "gStart", "gEnd"), sep = "[:-]" ), file = "atac_rna_linkages.txt", quote = F, sep = "\t", row.names = F, col.names = F )

# Load in the current genome sizes and TSS sites. 
txdb <- GenomicFeatures::makeTxDbFromGFF("/workdir/jdc397/1_currentWork/8_teCactus/00_clean/1_data/1_atac/mm10.refseqCurated.gtf")
tss <- promoters(txdb, upstream = 0, downstream = 3)

# Get a list of all the tss regions and write them to file!
df <- data.frame(seqnames = seqnames(tss), starts = start(tss)-1, ends = end(tss)); df <- subset(df, !grepl("Un|random|alt", seqnames))
write.table(unique(df), file = "tmp", row.names = F, col.names = F, sep = "\t", quote = F)
system(paste0("sort -k1,1 -k2,2n tmp | awk '{OFS=\"\\t\"; print $0, \"tss_\" ++count}' > refSC_mm10_tss.bed"))
system("rm tmp")

# Figure out which of the TSSs overlap an ATAC-RNA linkage and only use those for the shuffle
system("intersectBed -wa -a refSC_mm10_tss.bed -b atac_rna_linkages.txt | uniq > tmp")
a <- read.table("tmp"); colnames(a) <- c("chr", "start", "end", "id")
# write.table(a, file = "atac_rna_cor_tss.txt", sep = "\t", quote = F, row.names = F, col.names = F)
tss <- makeGRangesFromDataFrame(a, starts.in.df.are.0based=T, keep.extra.columns=T)

# Save the "tss" object so you can use it for the shuffle 
saveRDS(tss, file = "atac_rna_cor_tss.rds")

# Now get a data.table object where each row is a tss (from above) and the entries are a list of peaks that connect with that region
system("intersectBed -wa -wb -a atac_rna_linkages.txt -b tmp | cut -f8,12 | sort | uniq > int")
lut <- as.data.table(a); setkey(lut, id)
lut$val <- list()

# The below takes a few minutes, but is still *much* faster than normal. You can look into "set()" to try to make it even faster if possible...
int <- read.table("int"); colnames(int) <- c("peak_id", "tss_id")
for (f in 1:nrow(int)) {
	lut[ J(int[f,"tss_id"]), val := c(lut[J(int[f,"tss_id"]),val][[1]], int[f,"peak_id"]) ] 
}

# Save the "lut" object so you don't have to regenerate it again
saveRDS(lut, file = "atac_rna_cor_lut.rds")

# Remove all the objects/files you don't need anymore
system("rm tmp int"); rm(a, int, f, df)




