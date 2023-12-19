if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("enrichplot")
BiocManager::install("ReactomePA")

install.packages("gridExtra")
library(ggplot2)
library(gridExtra)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ReactomePA)

version
package_version <- packageVersion("clusterProfiler")
print(package_version)
# ----------------------------------------------------------------------------
# GSE152075 top 40 DEGs
#(ENTREZ ID of "DDX58" is "23586")
# DDX58 change to RIGI 
# gene list
gene_list <- c("IFI44L", "XAF1", "IFIT1", "OAS3", "OAS2", "IFIT3", "IFIT2", "RSAD2", "IGFBP2", "RIGI",
               "GBP1", "TRIM22", "EPSTI1", "MX2", "CD163", "CMPK2", "HERC6", "SAMD9", "CXCL10", "GBP4",
               "CRIP1", "PARP9", "RPLP1", "DDX60", "IFI44", "IFIT5", "RPS21", "RPS8", "FPR3", "PCSK5",
               "SAMD9L", "DDX60L", "OASL", "RPL13A", "CD300E", "PLA2G7", "ZEB2", "SBK1", "PRDX5", "RRAD")
fold_change <- c(4.5056, 3.3942, 4.3856, 3.7143, 3.5389, 4.0541, 4.3022, 3.626, -3.5954, 3.4351, 3.0214, 2.5574, 
                 3.2102, 2.7623, 3.3016, 3.376, 2.9625, 3.0783, 4.9611, 2.8892, -3.6812, 2.7013, -3.5888, 2.8777, 
                 2.9774, 2.8347, -3.2165, -2.9648, 3.6188, 3.0945, 3.1236, 3.0557, 3.593, -2.8992, 3.2175, 4.4151, 
                 2.835, 2.7974, -2.9724, -3.2282)
names(fold_change) <- gene_list

# Go analysis
bp_enrich <- enrichGO(gene         = gene_list,
                OrgDb        = org.Hs.eg.db,  
                keyType      = "SYMBOL",
                ont          = "BP",  
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05)

mf_enrich <- enrichGO(gene         = gene_list,
                 OrgDb        = org.Hs.eg.db,
                 keyType      = "SYMBOL",
                 ont          = "MF",
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05)

cc_enrich <- enrichGO(gene         = gene_list,
                 OrgDb        = org.Hs.eg.db,
                 keyType      = "SYMBOL",
                 ont          = "CC",
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05)

all_enrich <- enrichGO(gene         = gene_list,
                 OrgDb        = org.Hs.eg.db,  
                 keyType      = "SYMBOL",
                 ont          = "ALL",  
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05)

# Plot
pbp <- cnetplot(bp_enrich, circular = TRUE, foldChange = fold_change, colorEdge = TRUE, 
         node_label="gene", cex_label_category = 0.2, 
         cex_label_gene = 0.6, 
         node_size = 1) +
         theme(text = element_text(size = 8), 
         legend.text = element_text(size = 10)) +
         guides(colour = guide_legend(override.aes = list(size = 6)))
#         ggtitle("GO Analysis BP")
pbp <- pbp + guides(colour = guide_legend(ncol = 2))
print(pbp)

cnetplot(mf_enrich, circular = TRUE, foldChange = fold_change, colorEdge = TRUE, 
         node_label="gene", cex_label_category = 0.2, 
         cex_label_gene = 0.6, 
         node_size = 1) +
         theme(text = element_text(size = 8), 
         legend.text = element_text(size = 10)) +
         guides(colour = guide_legend(override.aes = list(size = 6)))

cnetplot(cc_enrich, circular = TRUE, foldChange = fold_change, colorEdge = TRUE, 
         node_label="gene", cex_label_category = 0.2, 
         cex_label_gene = 0.6, 
         node_size = 1) +
         theme(text = element_text(size = 8), 
         legend.text = element_text(size = 10)) +
         guides(colour = guide_legend(override.aes = list(size = 6)))


pall <- cnetplot(all_enrich, circular = TRUE, foldChange = fold_change, colorEdge = TRUE, 
         node_label="gene", cex_label_category = 0.2, 
         cex_label_gene = 0.6, 
         node_size = 1) +
         theme(text = element_text(size = 8), 
         legend.text = element_text(size = 10)) +
         guides(colour = guide_legend(override.aes = list(size = 6)))
pall <- pall + guides(colour = guide_legend(ncol = 2))
print(pall)

# KEGG
# Transfer gene symbol to ENTREZ ID
gene_list_KEGG <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# chech which gene can't be mapped 
unmapped_genes <- setdiff(gene_list, gene_list_KEGG$SYMBOL)

# Output the gene names which can't be mapped
if(length(unmapped_genes) > 0) {
  message("the gene names which can't be mapped：", paste(unmapped_genes, collapse = ", "))
}

# KEGG analysis
kegg_enrich <- enrichKEGG(gene  = gene_list_KEGG$ENTREZID,
                  organism   = 'hsa',
                  keyType    = 'kegg',
                  pAdjustMethod = 'BH',
                  qvalueCutoff = 0.05)

# Transfer ENTREZ ID to gene symbol
kegg_enrich <- setReadable(kegg_enrich, 'org.Hs.eg.db', 'ENTREZID')

# Plot
pkegg <- cnetplot(kegg_enrich, circular = TRUE, foldChange = fold_change, colorEdge = TRUE, 
         node_label="gene", cex_label_category = 0.2, 
         cex_label_gene = 0.6, 
         node_size = 1) +
         theme(text = element_text(size = 8), 
         legend.text = element_text(size = 10)) +
         guides(colour = guide_legend(override.aes = list(size = 6)))
pkegg <- pkegg + guides(colour = guide_legend(ncol = 2))
print(pkegg)

# Reactome
gene_list_entrez <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Reactome analysis
reactome_enrich <- enrichPathway(gene=gene_list_entrez$ENTREZID, organism = "human")

# Transfer ENTREZ ID to gene symbol
reactome_enrich <- setReadable(reactome_enrich, 'org.Hs.eg.db', 'ENTREZID')

# Plot
cnetplot(reactome_enrich, circular = TRUE, foldChange = fold_change, colorEdge = TRUE, 
         node_label="gene", cex_label_category = 0.2, 
         cex_label_gene = 0.6, 
         node_size = 1) +
         theme(text = element_text(size = 8), 
         legend.text = element_text(size = 10)) +
         guides(colour = guide_legend(override.aes = list(size = 6)))
# ----------------------------------------------------------------------------
# GSE152075 top Mito and organ DEGs
# gene list
gene_listMAO <- c("CXCL10", "RRAD", "USP18", "ATP5F1E", "CIB1", "C3AR1", "CYBB", "PROS1", "ACE2", "STUB1", "UQCRQ", "ND6", 
               "SLC8A1", "NDUFV1", "COX5A", "FLT1", "NDUFA13", "NDUFB7", "BAD", "ATP5ME", "NDUFAB1", "LAMB2", "SOCS3", 
               "PHB2", "TFF3", "KLF15")
fold_changeMAO <- c(4.9611, -3.2282, 2.0591, -2.2805, -3.1719, 2.4597, 2.9142, -1.9321, 1.8788, -2.2597, -2.152, 1.5331, 
                    2.7206, -2.4526, -1.5757, 3.0567, -2.1662, -2.1203, -1.9653, -1.9723, -2.0451, -2.0764, -2.0597, -1.7793, 
                    -2.6523, -1.9999)
names(fold_changeMAO) <- gene_listMAO  
  
# Go analysis
bp_enrichMAO <- enrichGO(gene         = gene_listMAO,
                      OrgDb        = org.Hs.eg.db,  
                      keyType      = "SYMBOL",
                      ont          = "BP",         
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05)

mf_enrichMAO <- enrichGO(gene         = gene_listMAO,
                      OrgDb        = org.Hs.eg.db,
                      keyType      = "SYMBOL",
                      ont          = "MF",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05)

cc_enrichMAO <- enrichGO(gene         = gene_listMAO,
                      OrgDb        = org.Hs.eg.db,
                      keyType      = "SYMBOL",
                      ont          = "CC",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05)

all_enrichMAO <- enrichGO(gene         = gene_listMAO,
                         OrgDb        = org.Hs.eg.db,  
                         keyType      = "SYMBOL",
                         ont          = "ALL",         
                         pAdjustMethod = "BH",
                         qvalueCutoff = 0.05)

# Plot
pbpMAO <- cnetplot(bp_enrichMAO, circular = TRUE, foldChange = fold_changeMAO, colorEdge = TRUE, 
         node_label="gene", cex_label_category = 0.2, 
         cex_label_gene = 0.6, 
         node_size = 1) +
         theme(text = element_text(size = 8), 
         legend.text = element_text(size = 10)) +
         guides(colour = guide_legend(override.aes = list(size = 6)))
pbpMAO <- pbpMAO + guides(colour = guide_legend(ncol = 2)) 
print(pbpMAO)

cnetplot(mf_enrichMAO, circular = TRUE, foldChange = fold_changeMAO, colorEdge = TRUE, 
         node_label="gene", cex_label_category = 0.2, 
         cex_label_gene = 0.6, 
         node_size = 1) +
         theme(text = element_text(size = 8), 
         legend.text = element_text(size = 10)) +
         guides(colour = guide_legend(override.aes = list(size = 6)))

pcc <- cnetplot(cc_enrichMAO, circular = TRUE, foldChange = fold_changeMAO, colorEdge = TRUE, 
         node_label="gene", cex_label_category = 0.2, 
         cex_label_gene = 0.6, 
         node_size = 1) +
         theme(text = element_text(size = 8), 
         legend.text = element_text(size = 10)) +
         guides(colour = guide_legend(override.aes = list(size = 6)))
pcc <- pcc + guides(colour = guide_legend(ncol = 2))
print(pcc)

pallMAO <- cnetplot(all_enrichMAO, circular = TRUE, foldChange = fold_changeMAO, colorEdge = TRUE, 
            node_label="gene", cex_label_category = 0.2, 
            cex_label_gene = 0.6, 
            node_size = 1) +
            theme(text = element_text(size = 8), 
            legend.text = element_text(size = 10)) +
            guides(colour = guide_legend(override.aes = list(size = 6)))
pallMAO <- pallMAO + guides(colour = guide_legend(ncol = 2)) 
print(pallMAO)

# KEGG
# Transfer gene symbol to ENTREZ ID
gene_list_KEGGMAO <- bitr(gene_listMAO, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# chech which gene can't be mapped 
unmapped_genesMAO <- setdiff(gene_listMAO, gene_list_KEGGMAO$SYMBOL)

# Output the gene names which can't be mapped
if(length(unmapped_genesMAO) > 0) {
  message("the gene names which can't be mapped：", paste(unmapped_genesMAO, collapse = ", "))
}

# KEGG analysis
kegg_enrichMAO <- enrichKEGG(gene  = gene_list_KEGGMAO$ENTREZID,
                          organism   = 'hsa',
                          keyType    = 'kegg',
                          pAdjustMethod = 'BH',
                          qvalueCutoff = 0.05)

# Transfer ENTREZ ID to gene symbol
kegg_enrichMAO <- setReadable(kegg_enrichMAO, 'org.Hs.eg.db', 'ENTREZID')

# Plot
cnetplot(kegg_enrichMAO, circular = TRUE, foldChange = fold_changeMAO, colorEdge = TRUE, 
              node_label="gene", cex_label_category = 0.2, 
              cex_label_gene = 0.6, 
              node_size = 1) +
              theme(text = element_text(size = 8), 
              legend.text = element_text(size = 10)) +
              guides(colour = guide_legend(override.aes = list(size = 6)))

# Reactome
gene_list_entrezMAO <- bitr(gene_listMAO, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Reactome analysis
reactome_enrichMAO <- enrichPathway(gene=gene_list_entrezMAO$ENTREZID, organism = "human")

# Transfer ENTREZ ID to gene symbol
reactome_enrichMAO <- setReadable(reactome_enrichMAO, 'org.Hs.eg.db', 'ENTREZID')

# Plot
# change line between "chemiosmotic" and "coupling"
reactome_enrichMAO@result$Description <- gsub("chemiosmotic coupling", "chemiosmotic\ncoupling", reactome_enrichMAO@result$Description, fixed = TRUE)

# cnetplot
prect <- cnetplot(reactome_enrichMAO, circular = TRUE, foldChange = fold_changeMAO, colorEdge = TRUE, 
              node_label="gene", cex_label_category = 0.2, 
              cex_label_gene = 0.6, 
              node_size = 1) +
              theme(text = element_text(size = 8), 
              legend.text = element_text(size = 10)) +
              guides(colour = guide_legend(override.aes = list(size = 6)))

# 显示图表
print(prect)
