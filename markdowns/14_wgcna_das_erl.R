# title: "JDM WGCNA"
# author: "Elisha Roberson"
# date: "2020-10-17"

# Setup
library( here )
library( tidyverse )
library( ggrepel )
library( reshape2 )
library( cowplot )
library( WGCNA )
library( dynamicTreeCut )

# source
source( here( 'src', 'shared_r_functions.R' ) )

# variable setup & multiple processors
wgcna_figure_directory <- file.path( figure_dir, 'wgcna' )
wgcna_output_dir <- here( "results", "wgcna" )

options( stringsAsFactors = FALSE )
enableWGCNAThreads()

# important reference variables
image_width <- 1280
image_height <- 720
powers <- c( 1:30 )
minModuleSize <- 10
MEDissThresh <- 0.25 # 75% correlation
pvalue_cutoff_modules <- 0.05
this_set <- paste0( 'PBMC' )
variance_quantile <- 0.10
scale_free_topology_cutoff <- 0.85
cytoscape_threshold <- 0.15 # --> after some hacking
#cytoscape_threshold <- 0.02 # --> from the tutorial

# Gene annotations
ensembl <- read_tsv( file = here( 'info', 'biomart_gene_transcript_map.txt.gz' ), col_names = TRUE ) %>%
  dplyr::rename( .data = ., gene_id = `Gene stable ID` ) %>%
  dplyr::rename( .data = ., symbol = `Gene name` ) %>%
  dplyr::rename( .data = ., gene_biotype = `Gene type` ) %>%
  select( .data = ., gene_id, symbol, gene_biotype ) %>%
  as.data.frame( . ) %>%
  unique( . )

ensembl_ids <- pull( ensembl, gene_id )

# 01 - Sample info and expression import
info <- read_csv( file = here( "info", "JDM_RNASeq_ReadGroup_Info.csv" ) ) %>%
  select( -RGID ) %>%
  filter( Status %in% c( 'JDM_Untreated_PBMC', 'JDM_Inactive_PBMC' ) ) %>%
  as.data.frame( . ) %>%
  unique( . )
rownames( info ) = info$RGSM

# trait info
trait_data <- read_csv( file = here( 'data', 'JDM_RNASeq_DAS_Info.csv' ) ) %>%
  as.data.frame( . ) %>%
  column_to_rownames( 'RGSM' )

# load regularized logarithm of the data
temp_expression <- read_tsv(
  file = here( 'results', 'unspliced_vst_expression.tsv.gz' ),
  col_names = TRUE,
  col_types = cols( .default = col_double(), gene_id = col_character() ) )

transposed_expression <- temp_expression[ , c( 'gene_id', rownames( info ) ) ] %>%
  as.data.frame( . ) %>%
  filter( gene_id %in% ensembl_ids ) %>%
  column_to_rownames( "gene_id" ) %>%
  as.matrix( x = . ) %>%
  t( x = . )

rm( temp_expression )

keep_names <- intersect( x = rownames( transposed_expression ), y = rownames( info ) ) %>%
  intersect( x = ., y = rownames( trait_data ) ) %>%
  sort( . )

info <- info[ keep_names, ]
trait_data <- trait_data[ keep_names, ]
transposed_expression <- transposed_expression[ keep_names, ]

# We've taken this down to just the samples we want to test
# but we can remove uninformative genes and perhaps have better results
variances <- apply( X = transposed_expression, MARGIN = 2, var, na.rm = TRUE )
variance_cutoff <- quantile( variances, variance_quantile ) %>% unname( . )
remove_column_index <- which( variances < variance_cutoff )

transposed_expression <- transposed_expression[ , -c( remove_column_index ) ]

# WGCNA tree
sampleTree <- transposed_expression %>%
  dist( . ) %>%
  hclust( ., method = 'average' )

jpeg(
  width = image_width,
  height = image_height ,
  filename = file.path( wgcna_figure_directory, paste0( this_set, "_initialSampleTree.jpeg" ) ),
  type = 'cairo' )
plot( sampleTree, sub="", xlab="" )
dev.off()

# recluster
sample_tree_2 <- transposed_expression %>%
  dist( . ) %>%
  hclust( ., method = 'average' )

trait_colors <- numbers2colors( trait_data )

jpeg(
  width = image_width,
  height = image_height ,
  filename = file.path( wgcna_figure_directory, paste0( this_set, "_dendroAndTraitColors.jpeg" ) ),
  type = 'cairo' )
plotDendroAndColors( sample_tree_2, trait_colors, groupLabels = names( trait_data ) )
dev.off()

# 02 - Setting thresholds and processing networks
sft <- pickSoftThreshold( transposed_expression,
                          dataIsExpr = TRUE,
                          powerVector = powers,
                          verbose = 5 )

# scale-free topology fit index as a function of the soft-thresholding power
independence_df <- data.frame(
  fitIndices = sft$fitIndices[ ,1 ],
  scaleFreeTopology = -sign( sft$fitIndices[ , 3] ) * sft$fitIndices[ ,2 ],
  connectivity = sft$fitIndices[ ,5 ],
  text = as.character( powers )
) %>%
  mutate( .data = ., rounded = round( scaleFreeTopology, 2 ) )

autocalculated_soft_power <- independence_df %>%
  filter( .data = ., rounded >= scale_free_topology_cutoff ) %>%
  arrange( fitIndices ) %>%
  .[ 1, "fitIndices" ]

stopifnot( !is.na( autocalculated_soft_power ) )

independence_df <- independence_df %>%
  mutate( isThreshold = case_when(
    text == as.character( autocalculated_soft_power ) ~ "Soft-threshold",
    TRUE ~ "Not Threshold" ) ) %>%
  mutate( isThreshold = factor( isThreshold, levels=c( "Not Threshold", "Soft-threshold" ) ) )

# thresholding plot
independence_ggplot <- ggplot( data = independence_df, aes( x = fitIndices, y = scaleFreeTopology, colour = isThreshold ) ) +
  geom_point() +
  geom_text_repel( aes( label = text ), size = 5, direction = 'y' ) +
  #geom_text( aes( label = text ), size = 5 ) +
  geom_hline( mapping = aes( yintercept = scale_free_topology_cutoff ), linetype = 2 ) +
  theme_bw() +
  theme( legend.position = 'none' ) +
  xlab( "Soft Threshold (power)" ) +
  ylab( parse( text = paste0( "Scale~Free~Toplogy~Model~Fit~(signed~R^2)" ) ) ) +
  ggtitle( "Scale independence" ) +
  gg_no_grid +
  gg_center_title

# mean connectivity as fxn of soft-thresholding power
connectivity_ggplot <- ggplot( data = independence_df, aes( x = fitIndices, y = connectivity, colour = isThreshold ) ) +
  geom_point() +
  geom_text_repel( aes( label = text ), size = 5, direction = 'y' ) +
  #geom_text( aes( label = text ), size = 5 ) +
  theme_bw() +
  theme( legend.position = 'none' ) +
  xlab( "Soft Threshold (power)" ) +
  ylab( "Mean connectivity\n" ) +
  ggtitle( "Mean connectivity" ) +
  gg_no_grid +
  gg_center_title

# ATTN-PLOT independence
jpeg(
  width = image_width,
  height = image_height ,
  filename = file.path( wgcna_figure_directory, paste0( this_set, "_independence.jpeg" ) ),
  type = 'cairo' )
independence_ggplot
dev.off()

# ATTN-PLOT connectivity
jpeg(
  width = image_width,
  height = image_height ,
  filename = file.path( wgcna_figure_directory, paste0( this_set, "_connectivity.jpeg" ) ),
  type = 'cairo' )
connectivity_ggplot
dev.off()

# ATTN-PLOT both
jpeg(
  width = image_width,
  height = image_height ,
  filename = file.path( wgcna_figure_directory, paste0( this_set, "_indepedence_connectivity_cowplot.jpeg" ) ),
  type = 'cairo' )
plot_grid( independence_ggplot, connectivity_ggplot, labels = c( 'A', 'B' ) )
dev.off()

cat( paste0( "The scale-independence will be chosen at about the ", scale_free_topology_cutoff, " threshold.\n" ) )
cat( paste0( "This is at approximately power ", autocalculated_soft_power, " for this dataset.\n" ) )

# run adjacencies
collectGarbage()

softPower <- autocalculated_soft_power

adjacency_out = adjacency( datExpr = transposed_expression, power = softPower, corFnc = bicor )

# convert adjacency to TOM
collectGarbage()

dissTOM <- 1.0 - TOMsimilarity( adjacency_out )

# wgcna plot TOM
geneTree <- dissTOM %>%
  as.dist( . ) %>%
  hclust( ., method = 'average' )

jpeg(
  width = image_width,
  height = image_height ,
  filename = file.path( wgcna_figure_directory, paste0( this_set, "_genetree_dissimilarity_TOM.jpeg" ) ),
  type = 'cairo' )

sizeGrWindow( 6, 4 )
plot( geneTree )

dev.off()

# select modules
dynamicMods = cutreeDynamic(
  dendro = geneTree,
  distM = dissTOM,
  deepSplit = 2,
  pamRespectsDendro = FALSE,
  minClusterSize = minModuleSize )

# wgcna gene dynamic colors
dynamicColors <- labels2colors( dynamicMods )
table( dynamicColors )

jpeg(
  width = image_width,
  height = image_height ,
  filename = file.path( wgcna_figure_directory, paste0( this_set, "_gene_dynamic_colors.jpeg" ) ),
  type = 'cairo' )

sizeGrWindow( 8,6 )
plotDendroAndColors(
  geneTree,
  dynamicColors,
  "Dynamic Tree Cut",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene dendrogram and module colors" )

dev.off()

# wgcna merge similar
MEList <- moduleEigengenes( transposed_expression, colors = dynamicColors )
MEs <- MEList$eigengenes

# dissimilarity
MEDiss <- 1.0 - bicor( MEs )
METree <- MEDiss %>%
  as.dist( . ) %>%
  hclust( ., method = 'average' )

jpeg(
  width = image_width,
  height = image_height ,
  filename = file.path( wgcna_figure_directory, paste0( this_set, "_eigengene_clustering.jpeg" ) ),
  type = 'cairo' )

sizeGrWindow( 7, 6 )
plot( METree, main = "Clustering of module eigengenes", xlab = "", sub = "" )

dev.off()

# merge
merge <- mergeCloseModules( transposed_expression, dynamicColors, cutHeight = MEDissThresh, verbose = 3 )

mergedColors = merge$colors
mergedMEs <- merge$newMEs

# wgcna plot merged clusters

jpeg(
  width = image_width,
  height = image_height ,
  filename = file.path( wgcna_figure_directory, paste0( this_set, "_MergedClusters.jpeg" ) ),
  type = 'cairo' )

sizeGrWindow( 12,9 )
plotDendroAndColors(
  geneTree,
  cbind( dynamicColors, mergedColors ),
  c( "Dynamic Tree Cut", "Merged dynamic" ),
  dendroLabels =  FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

dev.off()

# prep and save info
moduleColors <- mergedColors
colorOrder = c( "grey", standardColors( 50 ) )
moduleLabels <- match( moduleColors, colorOrder ) - 1
MEs <- mergedMEs

# save some of the these lengthy computations
save( MEs, moduleLabels, moduleColors, geneTree,
      file = file.path( wgcna_output_dir,  paste0( this_set, "_step_by_step_network.RData" ) ) )

# 03 - relate modules to external information
# define ngenes samples
nGenes <- ncol( transposed_expression )
nSamples <- nrow( transposed_expression )

# recalc MEs with color labels
MEs0 <- moduleEigengenes( transposed_expression, moduleColors )$eigengenes
MEs <- orderMEs( MEs0 )

moduleTraitCor <- bicor( x = MEs, y = trait_data, use = 'p' )
moduleTraitPvalue <- corPvalueStudent( cor = moduleTraitCor, nSamples = nSamples )

# wgcna module trait graphic
textMatrix <- paste0(
  signif( moduleTraitCor, 2 ),
  "\n(",
  signif( moduleTraitPvalue, 1 ),
  ")" )

dim( textMatrix ) <- dim( moduleTraitCor )

# display correlation values within heatmap plot
jpeg(
  width = image_width,
  height = image_height ,
  filename = file.path( wgcna_figure_directory, paste0( this_set, "_correlationHeatmap.jpeg" ) ),
  type = 'cairo' )

sizeGrWindow( 10, 6 )
par( mar = c( 6, 8.5, 3, 3 ) )
labeledHeatmap( Matrix = moduleTraitCor,
                xLabels = names( trait_data ),
                yLabels = names( MEs ),
                colorLabels = FALSE,
                colors = blueWhiteRed( 50 ),
                textMatrix = textMatrix,
                setStdMargins = FALSE,
                cex.text = 0.5,
                zlim = c( -1, 1 ),
                main = "Module-trait relationships" )
dev.off()

## Module membership
# gene signif module member
modNames <- names( MEs )

geneModuleMembership <- bicor( x = transposed_expression, y = MEs, use = 'p' ) %>%
  as.data.frame( . )

MMPvalue <- geneModuleMembership %>%
  as.matrix( x = . ) %>%
  corPvalueStudent( cor = ., nSamples = nSamples ) %>%
  as.data.frame( x = . )

names( geneModuleMembership ) <- paste0( "MM_", modNames )
names( MMPvalue ) <- paste0( "p.MM_", modNames )

## Trait associations
# trait associations

trait_list <- names( trait_data )
pvalue_colors <- rownames( moduleTraitPvalue )

signif_colors_for_later <- c()

for( trait_idx in 1:length( trait_list ) )
{
  curr_trait <- trait_list[ trait_idx ]

  trait_df <- trait_data[ , curr_trait, drop = FALSE ]

  geneTraitSignificance <- bicor( x = transposed_expression, y = trait_df, use = 'p' ) %>%
    as.data.frame( . )

  GSPvalue <- as.matrix( geneTraitSignificance ) %>%
    corPvalueStudent( cor = ., nSamples = nSamples ) %>%
    as.data.frame( . )

  names( geneTraitSignificance ) <- paste0( "GS.corr.Status" )
  names( GSPvalue ) <- paste0( "p.GS.Status" )

  # now process each signif module for the trait
  signif_colors <- pvalue_colors[ which( moduleTraitPvalue[ , curr_trait ] < pvalue_cutoff_modules ) ]
  signif_colors_for_later = c( signif_colors_for_later, signif_colors )

  tmp <- data.frame(
    gene_id = colnames( transposed_expression ),
    moduleColor = moduleColors,
    geneTraitSignificance,
    GSPvalue ) %>%
    merge( x = ., y = ensembl, by = "gene_id", all.x = TRUE ) %>%
    select( ., gene_id, symbol, gene_biotype, moduleColor, GS.corr.Status, p.GS.Status )

  modOrder <- MEs %>%
    bicor( x = ., y = trait_df, use = 'p' ) %>%
    -abs( x = . ) %>%
    order( . )

  for ( mod_idx in 1:ncol( geneModuleMembership ) )
  {
    oldNames <- names( tmp )
    curr_col <- modOrder[ mod_idx ]
    tmp <- data.frame( tmp,
                       geneModuleMembership[ , curr_col ],
                       MMPvalue[ , curr_col ] )
    names( tmp ) <- c(
      oldNames,
      paste0( "MM.", modNames[ curr_col ] ),
      paste0( "p.MM.", modNames[ curr_col ] ) )
  }

  fname <- file.path(
    wgcna_output_dir,
    paste0( this_set, "_WGCNA_table_", curr_trait, ".tsv" ) )

  fname_down <- file.path(
    wgcna_output_dir,
    paste0( this_set, "_WGCNA_table_", curr_trait, "_down_genes.tsv" ) )

  fname_up <- file.path(
    wgcna_output_dir,
    paste0( this_set, "_WGCNA_table_", curr_trait, "_up_genes.tsv" ) )

  tmp_signif_cor_all <- tmp %>%
    filter( moduleColor %in% str_replace( string = signif_colors, pattern = "^ME", replacement = "" ) ) %>%
    arrange( p.GS.Status, GS.corr.Status )

  if ( nrow( tmp_signif_cor_all ) > 0 ) {
    tmp_signif_cor_all %>%
      write_tsv( x = ., path = fname )

    tmp_signif_cor_down <- tmp_signif_cor_all %>%
      filter( p.GS.Status < 0.05 & GS.corr.Status < 0.0 ) %>%
      arrange( p.GS.Status, GS.corr.Status )

    if ( nrow( tmp_signif_cor_down ) > 0 ) {
      stopifnot( tmp_signif_cor_down$GS.corr.Status < 0.0 )

      tmp_signif_cor_down %>%
        write_tsv( x = ., path = fname_down )
    }

    tmp_signif_cor_up <- tmp_signif_cor_all %>%
      filter( p.GS.Status < 0.05 & GS.corr.Status > 0.0 ) %>%
      arrange( p.GS.Status, GS.corr.Status )

    if ( nrow( tmp_signif_cor_up ) > 0 ) {
      stopifnot( tmp_signif_cor_up$GS.corr.Status > 0.0 )

      tmp_signif_cor_up %>%
        write_tsv( x = ., path = fname_up )
    }
  }
}

# 04 - interfacing network analysis of other data such as functional annotation and gene ontology
# Can be manually run through gProfileR

# 05 - network visualization using WGCNA functions
# Skip

# 06 - exporting a gene network to external visualization software

# wgcna export for gephi
signif_colors_for_later <- signif_colors_for_later %>%
  unique( . ) %>%
  str_replace( string = ., pattern = "^ME", replacement = "" )

TOM <- TOMsimilarityFromExpr(
  datExpr = transposed_expression,
  power = softPower,
  corType = 'bicor',
  networkType = 'unsigned' )

modules <- signif_colors_for_later
gene_list <- colnames( transposed_expression )

inModule <- match( moduleColors, modules ) %>%
  is.finite( . )

modProbes <- gene_list[ inModule ]
modGenes <- data.frame( gene_id = modProbes ) %>%
  merge( x = ., y = ensembl, by = "gene_id", all.x = TRUE ) %>%
  mutate( ., symbol = case_when(
    is.na( symbol ) ~ gene_id,
    TRUE ~ symbol ) ) %>%
  pull( symbol )

modTOM <- TOM[ inModule, inModule ]

dimnames( modTOM ) <- list( modProbes, modProbes )

edgeFilePath <- file.path(
  wgcna_output_dir,
  paste0( this_set, "_Cytoscape-edges.txt" ) )

gephiEdgeFilePath <- file.path(
  wgcna_output_dir,
  paste0( this_set, "_Gephi-edges.csv.gz" ) )

nodeFilePath <- file.path(
  wgcna_output_dir,
  paste0( this_set, "_Cytoscape-nodes.txt" ) )

gephiNodeFilePath <- file.path(
  wgcna_output_dir,
  paste0( this_set, "_Gephi-nodes.csv.gz" ) )

cyt <- exportNetworkToCytoscape(
  adjMat = modTOM,
  edgeFile = edgeFilePath,
  nodeFile = nodeFilePath,
  weighted = TRUE,
  threshold = cytoscape_threshold,
  nodeNames = modGenes,
  altNodeNames = modProbes,
  nodeAttr = moduleColors[ inModule ] )

# setup for GEPHI import
read_tsv( file = edgeFilePath ) %>%
  dplyr::rename( ., Source = fromNode ) %>%
  dplyr::rename( ., Target = toNode ) %>%
  select( ., Source, Target, weight ) %>%
  write_csv( x = ., path = gephiEdgeFilePath )

# warning - there may be a comma-separated list embedded in the membership column.
# it wasn't in this case, but fair warning.
read_tsv( file = nodeFilePath ) %>%
  dplyr::rename( Id = nodeName ) %>%
  dplyr::rename( network = `nodeAttr[nodesPresent, ]` ) %>%
  select( Id, network ) %>%
  write_csv( x = ., path = gephiNodeFilePath )

# Session info
Sys.time()

getwd()

sessionInfo()

