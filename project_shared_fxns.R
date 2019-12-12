####################
# libraries needed #
####################
library( here )
library( ggplot2 ) # for ggplot modifiers
library( plyr )
library( dplyr )
library( ggrepel )
library( purrr )

#############
# constants #
#############
present_count_cutoff <- 5
min_detected_samples <- 5
fold_change_for_graphs_y_gprofiler <- 1.5

###############
# directories #
###############
path_to_info <- here::here( 'info' )
path_to_data <- here::here( 'data' )
path_to_results <- here::here( 'results' )
wgcna_output_path <- here::here( 'results', "wgcna_output" )
path_to_figures <- here::here( 'results', 'figures' )
path_to_gprofiler_figures <- here::here( 'results', 'gprofile_figures' )
path_to_gprofiles <- here::here( 'results', 'gprofiler_output_maxSize500' )
path_to_gene_figures <- here::here( 'results', 'individual_genes_all' )
path_to_gene_figures_untreated <- here::here( 'results', 'individual_genes_untreated' )

dir.create( here::here( 'info' ), showWarnings = FALSE )
dir.create( path_to_results, showWarnings = FALSE )
dir.create( wgcna_output_path, showWarnings = FALSE )
dir.create( path_to_data, showWarnings = FALSE )
dir.create( path_to_figures, showWarnings = FALSE )
dir.create( path_to_gprofiler_figures, showWarnings = FALSE )
dir.create( path_to_gprofiles, showWarnings = FALSE )
dir.create( path_to_gene_figures, showWarnings = FALSE )
dir.create( path_to_gene_figures_untreated, showWarnings = FALSE )

##############
# file paths #
##############
project_rgid_filename <- "JDM_RNASeq_ReadGroup_Info.csv"
clinical_info_filename <- "JDM_RNASeq_DAS_Info.csv"

project_rgid_path <- here::here( 'info', project_rgid_filename )
clinical_info_path <- here::here( 'data', clinical_info_filename )
merged_count_rgid_path <- here::here( 'data', "counts_by_rgid.csv.gz" )
merged_count_path <- here::here( 'results', "rgsm_merged_counts.tsv" )
das_category_path <- here::here( 'results', "das_categories.csv" )
sex_from_expression_path <- here::here( 'results', "sample_sex_clustering.csv" )

ensembl_annotation_path <- here::here( 'info', 'gene_name_maps.tsv.gz' )

sample_intensities_path <- here::here( 'results', "all_sample_deseq2_logrratio.csv" )
paired_treated_samples_path <- here::here( 'results', "paired_treated_sample_deseq2_logrratio.csv" )
paired_pre_post_samples_path <- here::here( 'results', "paired_pre_post_psoriasis_sample_deseq2_logrratio.csv" )

genes_with_sex_linked_effect <- here::here( 'results', "sex_linked_genes.csv" )

rtqpcr_raw_data_path <- here::here( 'data', 'OtOf_RTqPCR_tidy.csv' )

##############################
# Download data for analysis #
##############################
file_link_map <- list(
'read_group_info' = list( 'https://ndownloader.figshare.com/files/15188312', project_rgid_path ),
'gene_counts' = list( 'https://ndownloader.figshare.com/files/15188324', merged_count_rgid_path ),
'das_scores' = list( 'https://ndownloader.figshare.com/files/15188321', clinical_info_path ),
'rtqpcr' = list ( 'https://ndownloader.figshare.com/files/15188327', rtqpcr_raw_data_path )
)

##################################
# Differentially expressed genes #
##################################
deg_file_untreated_control <- here::here( 'results', "deg_untreated_vs_control.csv" )

deg_file_inactive_control <- here::here( 'results', "deg_inactive_vs_control.csv" )

deg_file_inactive_untreated_paired <- here::here( 'results', "deg_untreated_vs_inactive_paired.csv" )

deg_file_prepsor_control <- here::here( 'results', "deg_prepsor_vs_control.csv" )

deg_file_withpsor_control <- here::here( 'results', "deg_psor_vs_control.csv" )

deg_file_prepsor_psor_paired <- here::here( 'results', "deg_prepsor_vs_psor_paired.csv" )

####################
# gprofiler inputs #
####################
gprofiler_untreated_control_up <- here::here( 'results', "gprofiler_input_untreated_control_up.csv" )
gprofiler_untreated_control_down <- here::here( 'results', "gprofiler_input_untreated_control_down.csv" )

gprofiler_inactive_control_up <- here::here( 'results', "gprofiler_input_inactive_control_up.csv" )
gprofiler_inactive_control_down <- here::here( 'results', "gprofiler_input_inactive_control_down.csv" )

gprofiler_inactive_untreated_paired_up <- here::here( 'results', "gprofiler_input_untreated_v_inactive_paired_up.csv" )
gprofiler_inactive_untreated_paired_down <- here::here( 'results', "gprofiler_input_untreated_v_inactive_paired_down.csv" )

gprofiler_prepsor_control_up <- here::here( 'results', "gprofiler_input_prepsor_control_up.csv" )
gprofiler_prepsor_control_down <-here::here( 'results', "gprofiler_input_prepsor_control_down.csv" )

gprofiler_withpsor_control_up <- here::here( 'results', "gprofiler_input_withpsor_control_up.csv" )
gprofiler_withpsor_control_down <- here::here( 'results', "gprofiler_input_withpsor_control_down.csv" )

gprofiler_prepsor_psor_paired_up <- here::here( 'results', "gprofiler_input_prepsor_psor_paired_up.csv" )
gprofiler_prepsor_psor_paired_down <- here::here( 'results', "gprofiler_input_prepsor_psor_paired_down.csv" )

###############################
# color blind friendly colors #
###############################
colorBlindPalette <- c( "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7" )
severityPalette <- c( "#808080", "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c" )
blues5Palette <- c( '#ffffcc', '#a1dab4', '#41b6c4', '#2c7fb8', '#253494' )
greens5Palette <- c( '#ffffcc', '#c2e699', '#78c679', '#31a354', '#006837' )
purples5Palette <- c( '#feebe2', '#fbb4b9', '#f768a1', '#c51b8a', '#7a0177' )
reds5Palette <- c( '#ffffb2', '#f3cc5c', '#fd8d3c', '#f03b20', '#bd0026' )
reds3Palette <- c( '#fee0d2', '#fc9272', '#de2d26' )

####################
# ggplot constants #
####################
gg_bigger_texts <- theme(
    axis.title = element_text( size=22 ),
    axis.text = element_text( size=20 ),
    legend.text = element_text( size=14 ),
    legend.title = element_text( size=15 ),
    plot.title = element_text( size=22 ),
	strip.text = element_text( size=15 )
)

gg_merged_texts <- theme(
    axis.title = element_text( size=18 ),
    axis.text = element_text( size=17 ),
    legend.text = element_text( size=13 ),
    legend.title = element_text( size=14 ),
    plot.title = element_text( size=18 ),
	strip.text = element_text( size=13 )
)

gg_no_legend <- theme(
	legend.position='none'
)

gg_no_grid <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)

gg_no_x_grid <- theme(
  panel.grid.major.x = element_blank() )

gg_no_y_grid <- theme(
  panel.grid.major.y = element_blank() )

gg_center_title <- theme(
	plot.title = element_text( hjust = 0.5 )
)

gg_gprofile_quadplot = theme(
    axis.title = element_text( size=14 ),
    axis.text = element_text( size=9 ),
    plot.title = element_text( size=15 )
)

###########################
# Detect values in matrix #
###########################
detect <- function( input_vector, cutoff, upper=TRUE )
{
	if ( upper == TRUE )
	{
		return( length( which( input_vector >= cutoff ) ) )
	} else
	{
		return( length( which( input_vector <= cutoff ) ) )
	}
}

################################
# pretty format for DE results #
################################
format_deseq_results <- function( deseq_results, annotation_df ) {
	as.data.frame( deseq_results ) %>%
	rownames_to_column( "gene_id" ) %>%
	filter( !is.na( padj ) ) %>%
	filter( !is.na( pvalue ) ) %>%
	mutate( log2FC = log2FoldChange ) %>%
	mutate( FC = 2^log2FC ) %>%
	mutate( FC = case_when(
		FC < 1.0 ~ -1.0 / FC,
		TRUE ~ FC )
	) %>%
	mutate( pval = pvalue ) %>%
	mutate( qval = padj ) %>%
	select( gene_id, FC, log2FC, pval, qval ) %>%
	merge( annotation_df, ., by="gene_id" ) %>%
	arrange( qval, pval ) %>%
	return( . )
}

##################################
# Report differential expression #
##################################
summarize_pretty_de_results <- function( pretty_de, fc_cutoff ) {
  fc_cutoff <- abs( fc_cutoff )

	down <- pretty_de %>%
		filter( qval < 0.05 & FC < 0.0 ) %>%
		pull( gene_id ) %>%
		length( . )

	down_fc <- pretty_de %>%
		filter( qval < 0.05 & FC < ( -1.0 * fc_cutoff ) ) %>%
		pull( gene_id ) %>%
		length( . )

	up <- pretty_de %>%
		filter( qval < 0.05 & FC > 0.0 ) %>%
		pull( gene_id ) %>%
		length( . )

	up_fc <- pretty_de %>%
		filter( qval < 0.05 & FC > fc_cutoff ) %>%
		pull( gene_id ) %>%
		length( . )

	data.frame( Change = c( 'Down', 'Down - FC', 'Up', 'Up - FC' ), n = c( down, down_fc, up, up_fc ) ) %>%
		return( . )
}

#################
# volcano plots #
#################
make_ggplot_volcano <- function( deg_dataframe, case_name, control_name, axis_steps = 2, fold_change_cutoff = 1.5, qvalue_cutoff = 0.05, max_label = 30 )
{
	##############################
	# set significance threshold #
	##############################
	deg_dataframe <- deg_dataframe %>%
	mutate( Significant = case_when(
		qval < qvalue_cutoff & abs( FC ) >= fold_change_cutoff ~ "Large",
		qval < qvalue_cutoff ~ "Modest",
		TRUE ~ "Not" ) ) %>%
	mutate( Significant = factor( Significant, levels=c( "Not", "Modest", "Large" ) ) )

	################################
	# set values for square x axis #
	################################
	x_volcano_value <- ( abs( deg_dataframe$log2FC[ is.finite( deg_dataframe$log2FC ) ] ) + 0.051 ) %>%
		max( . ) %>%
		round( ., 1 )

	if ( x_volcano_value < 1.0 ) {
		x_volcano_value = 1.0
	}

	x_num_for_limits <- round( x_volcano_value, 0 )

	x_volcano_low <- x_volcano_value * -1
	x_volcano_high <- x_volcano_value

	x_break_list <- seq( -1 * x_num_for_limits, x_num_for_limits, by = axis_steps )

	##############
	# plot lines #
	##############
	horizontal_line <- log10( qvalue_cutoff ) * -1
	vertical_line_1 <- log2( fold_change_cutoff )
	vertical_line_2 <- vertical_line_1 * -1

	###################################
	# actually make the volcano plots #
	###################################
	plot_volcano <- ggplot( deg_dataframe, aes( x=log2FC, y=-log10( qval ), colour=Significant ) ) +
		scale_colour_manual( values = c( "darkgray", "blue", "red" ) ) +
		scale_x_continuous( limits = c( x_volcano_low, x_volcano_high ), breaks = x_break_list ) +
		theme_bw() +
		gg_bigger_texts +
		gg_no_legend +
		gg_no_grid +
		gg_center_title +
		geom_point( size=1.2 ) +
		geom_hline( yintercept = horizontal_line, linetype=2 ) +
		geom_vline( xintercept=c( vertical_line_1, vertical_line_2 ), linetype=2 ) +
		geom_text_repel( data=subset( deg_dataframe, Significant == "Large" )[c(1:max_label),], colour="black", aes( label=symbol ), size=3 ) +
		xlab( parse( text=paste0( "log[2]~(", case_name, "/", control_name, ")" ) ) ) +
		ylab( parse( text = paste0( "-log[10]~(Adj.~p-value)" ) ) )

	return( plot_volcano )
}

###################################
# make lists of up and down genes #
###################################
direction_gene_ids_from_de <- function( deg_df, qval_cutoff = 0.05 ) {
  # gene_id lists
  tested_genes_gid <- pull( deg_df, gene_id )

  deg_gid <- filter( deg_df, qval < qval_cutoff ) %>%
    pull( gene_id )

  down_gid <- filter( deg_df, qval < qval_cutoff & FC < 0.0 ) %>%
    pull( gene_id )

  up_gid <- filter( deg_df, qval < qval_cutoff & FC > 0.0 ) %>%
    pull( gene_id )

  unchanged_gid <- filter( deg_df, qval >= qval_cutoff ) %>%
    pull( gene_id )

  # gene class counts
  total_size <- length( tested_genes_gid )

  de_size <- length( deg_gid )

  up_size <- length( up_gid )

  down_size <- length( down_gid )

  # some validity checking
  stopifnot( de_size == up_size + down_size )
  stopifnot( all( tested_genes_gid %in% c( down_gid, up_gid, unchanged_gid ) ) )
  stopifnot( !any( deg_gid %in% unchanged_gid ) )

  # build the output list
  list( 'down_gene_ids' = down_gid,
        'up_gene_ids' = up_gid,
        'unchanged_gene_ids' = unchanged_gid,
        'deg_gid' = deg_gid,
        'all_tested_gid' = tested_genes_gid,
        'num_tested' = total_size,
        'num_de' = de_size,
        'num_up' = up_size,
        'num_down' = down_size ) %>%
    return( . )
}

#############################
# Set before / after status #
#############################
setup_alluvial_table_paired_data <- function( left_direction_list, right_direction_list ) {
  # Left should be case / control.
  # Right should be case timepoint 2 / case timepoint 1
  up2improved <- intersect( left_direction_list[[ 'up_gene_ids' ]], right_direction_list[[ 'down_gene_ids' ]] ) %>% length( . )

  up2worsened <- intersect( left_direction_list[[ 'up_gene_ids' ]], right_direction_list[[ 'up_gene_ids' ]] ) %>% length( . )

  up2unchanged <- setdiff( left_direction_list[[ 'up_gene_ids' ]], right_direction_list[[ 'deg_gid' ]] ) %>% length( . )

  down2improved <- intersect( left_direction_list[[ 'down_gene_ids' ]], right_direction_list[[ 'up_gene_ids' ]] ) %>% length( . )

  down2worsened <- intersect( left_direction_list[[ 'down_gene_ids' ]], right_direction_list[[ 'down_gene_ids' ]] ) %>% length( . )

  down2unchanged <- setdiff( left_direction_list[[ 'down_gene_ids' ]], right_direction_list[[ 'deg_gid' ]] ) %>% length( . )

  unchanged2up <- setdiff( right_direction_list[[ 'up_gene_ids' ]], left_direction_list[[ 'deg_gid' ]] ) %>% length( . )

  unchanged2down <- setdiff( right_direction_list[[ 'down_gene_ids' ]], left_direction_list[[ 'deg_gid' ]] ) %>% length( . )

  unchanged2unchanged <- union( left_direction_list[[ 'all_tested_gid' ]], right_direction_list[[ 'all_tested_gid' ]] ) %>%
    setdiff( ., union( left_direction_list[[ 'deg_gid' ]], right_direction_list[[ 'deg_gid' ]] ) ) %>% length( . )

  data.frame( Left = "Up", Right = "Improved", n = up2improved ) %>%
    rbind( ., data.frame( Left = "Up", Right = "Worsened", n = up2worsened ) ) %>%
    rbind( ., data.frame( Left = "Up", Right = "Unchanged", n = up2unchanged ) ) %>%
    rbind( ., data.frame( Left = "Down", Right = "Improved", n = down2improved ) ) %>%
    rbind( ., data.frame( Left = "Down", Right = "Worsened", n = down2worsened ) ) %>%
    rbind( ., data.frame( Left = "Down", Right = "Unchanged", n = down2unchanged ) ) %>%
    rbind( ., data.frame( Left = "Unchanged", Right = "Up", n = unchanged2up ) ) %>%
    rbind( ., data.frame( Left = "Unchanged", Right = "Down", n = unchanged2down ) ) %>%
    rbind( ., data.frame( Left = "Unchanged", Right = "Unchanged", n = unchanged2unchanged ) ) %>%
    return( . )
}

old_setup_alluvial_table_paired_data <- function( left_direction_list, right_direction_list ) {
  # Left should be case / control.
  # Right should be case timepoint 2 / case timepoint 1
  up2improved <- which( left_direction_list[[ 'up_gene_ids' ]] %in% right_direction_list[[ 'down_gene_ids' ]] ) %>% length( . )

  up2worsened <- which( left_direction_list[[ 'up_gene_ids' ]] %in% right_direction_list[[ 'up_gene_ids' ]] ) %>% length( . )

  ## ATTN
  up2unchanged <- left_direction_list[[ 'num_up' ]] - ( length( up2improved ) + length( up2worsened ) )

  down2improved <- which( left_direction_list[[ 'down_gene_ids' ]] %in% right_direction_list[[ 'up_gene_ids' ]] ) %>% length( . )

  down2worsened <- which( left_direction_list[[ 'down_gene_ids' ]] %in% right_direction_list[[ 'down_gene_ids' ]] ) %>% length( . )

  ## ATTN
  down2unchanged <- left_direction_list[[ 'num_down' ]] - ( length( down2improved ) + length( down2worsened ) )

  unchanged2up <- which( !( right_direction_list[[ 'up' ]] %in% left_direction_list[[ 'full_de_list']] ) ) %>%
    length( . )

  unchanged2down <- which( !( right_direction_list[[ 'down' ]] %in% left_direction_list[[ 'full_de_list' ]] ) ) %>%
    length( . )

  data.frame( Left = "Up", Right = "Improved", n = up2improved ) %>%
    rbind( ., data.frame( Left = "Up", Right = "Worsened", n = up2worsened ) ) %>%
    rbind( ., data.frame( Left = "Up", Right = "Unchanged", n = up2unchanged ) ) %>%
    rbind( ., data.frame( Left = "Down", Right = "Improved", n = down2improved ) ) %>%
    rbind( ., data.frame( Left = "Down", Right = "Worsened", n = down2worsened ) ) %>%
    rbind( ., data.frame( Left = "Down", Right = "Unchanged", n = down2unchanged ) ) %>%
    rbind( ., data.frame( Left = "Unchanged", Right = "Up", n = unchanged2up ) ) %>%
    rbind( ., data.frame( Left = "Unchanged", Right = "Down", n= unchanged2down ) ) %>%
    return( . )
}

#######################################################
# string list of gene_ids to a string list of symbols #
#######################################################
convertRefseqList2SymbolList <- function( input_string, name_frame, string_sep = ",", out_sep=";" ) {
  str_split( string = input_string, pattern = string_sep ) %>%
    unlist( "," ) %>%
    map_chr( .x = ., ~ name_frame[ .x, 'symbol' ] ) %>%
    sort( . ) %>%
    paste( ., collapse = out_sep ) %>%
    return( . )
}

#################################################
# process raw gprofile into an updated gprofile #
#################################################
parse_raw_gprofile <- function( gprof_file_path, symbol_df, title.truncate = 65, title.wrap = 35 ) {
  gprof_colnames = c( "ignore", "signf", "pvalue", "T", "Q", "Q&T", "Q&T/Q", "Q&T/T", "term_ID", "t_type", "t group", "t_name", "t depth", "Q&T list" )

  new_filename_string <- str_replace( gprof_file_path, pattern = "v00", "v01" )

  read_tsv( file = gprof_file_path, comment="#", col_names = gprof_colnames, col_types = c("ccdiiiddccicic") ) %>%
    mutate( name = str_replace( t_name, "^ +", "" ) ) %>%
    select( -ignore, -signf ) %>%
    mutate( symbols = map( `Q&T list`, convertRefseqList2SymbolList, name_frame = symbol_df ) ) %>%
    unnest( . ) %>%
    mutate( title = paste0( term_ID, ' - ', t_name ) ) %>%
    mutate( title = str_replace_all( title, pattern = "\\; +match class\\: +[0-9]+", replacement = "" ) ) %>%
    mutate( title = str_trunc( title, width = title.truncate ) ) %>%
    mutate( title = str_wrap( title, width = title.wrap ) ) %>%
    select( pvalue, term_ID, t_type, t_name, symbols, `Q&T list`, title ) %>%
    mutate( `Q&T list` = str_replace( string = `Q&T list`, pattern = ",", replacement = ";" ) ) %>%
    arrange( pvalue ) %>%
    write_tsv( path = new_filename_string ) %>%
    return( . )
}

###########################
# Make multitude of plots #
###########################
gprofiler_ggplot_bonanza <- function( full.down.tbl, full.up.tbl, contrast_name, gprofile_figure_path ) {
  if ( !is.null( dim( full.down.tbl ) ) ) {
    plot.down.tbl <- full.down.tbl %>%
      ddply( .variables = c( "t_type" ), .fun = top_n, wt = pvalue, n = -6 )

    ################################
    # individual decreased pathway #
    ################################
    bp.plt.down <- filter( plot.down.tbl, t_type == "BP" ) %>%
      arrange( desc( pvalue ) ) %>%
      mutate( title = factor( title, levels=title ) ) %>%
      ggplot( aes( y = -log10( pvalue ), x = title ) ) +
      theme_classic() +
      gg_bigger_texts +
      gg_center_title +
      coord_flip() +
      geom_col() +
      geom_hline( yintercept = -log10( 0.05 ), linetype=2, colour="white"  ) +
      ylab( parse( text="-log[10](P-value)" ) ) +
      xlab( "Pathway\n" ) +
      ggtitle( "GO Biological Process" )

    cc.plt.down <- filter( plot.down.tbl, t_type == "CC" ) %>%
      arrange( desc( pvalue ) ) %>%
      mutate( title = factor( title, levels=title ) ) %>%
      ggplot( aes( y = -log10( pvalue ), x = title ) ) +
      theme_classic() +
      gg_bigger_texts +
      gg_center_title +
      coord_flip() +
      geom_col() +
      geom_hline( yintercept = -log10( 0.05 ), linetype=2, colour="white"  ) +
      ylab( parse( text="-log[10](P-value)" ) ) +
      xlab( "Pathway\n" ) +
      ggtitle( "GO Cellular Compartment" )

    kegg.plt.down  <- filter( plot.down.tbl, t_type == "keg" ) %>%
      arrange( desc( pvalue ) ) %>%
      mutate( title = factor( title, levels=title ) ) %>%
      ggplot( aes( y = -log10( pvalue ), x = title ) ) +
      theme_classic() +
      gg_bigger_texts +
      gg_center_title +
      coord_flip() +
      geom_col() +
      geom_hline( yintercept = -log10( 0.05 ), linetype=2, colour="white"  ) +
      ylab( parse( text="-log[10](P-value)" ) ) +
      xlab( "Pathway\n" ) +
      ggtitle( "KEGG" )

    react.plt.down <- filter( plot.down.tbl, t_type == "rea" ) %>%
      arrange( desc( pvalue ) ) %>%
      mutate( title = factor( title, levels=title ) ) %>%
      ggplot( aes( y = -log10( pvalue ), x = title ) ) +
      theme_classic() +
      gg_bigger_texts +
      gg_center_title +
      coord_flip() +
      geom_col() +
      geom_hline( yintercept = -log10( 0.05 ), linetype=2, colour="white"  ) +
      ylab( parse( text="-log[10](P-value)" ) ) +
      xlab( "Pathway\n" ) +
      ggtitle( "Reactome" )

    jpeg( file = file.path( gprofile_figure_path, paste0( contrast, "_BP_down.jpeg" ) ), width=1280, height=1024, res=75 )
    print( bp.plt.down )
    dev.off()

    jpeg( file = file.path( gprofile_figure_path, paste0( contrast, "_CC_down.jpeg" ) ), width=1280, height=1024, res=75 )
    print( cc.plt.down )
    dev.off()

    jpeg( file = file.path( gprofile_figure_path, paste0( contrast, "_KEGG_down.jpeg" ) ), width=1280, height=1024, res=75 )
    print( kegg.plt.down )
    dev.off()

    jpeg( file = file.path( gprofile_figure_path, paste0( contrast, "_REAC_down.jpeg" ) ), width=1280, height=1024, res=75 )
    print( react.plt.down )
    dev.off()

    bp.plt.down <- bp.plt.down + gg_gprofile_quadplot

    cc.plt.down <- cc.plt.down + gg_gprofile_quadplot

    kegg.plt.down <- kegg.plt.down + gg_gprofile_quadplot

    react.plt.down <- react.plt.down + gg_gprofile_quadplot

    multi.plt.down <- plot_grid( bp.plt.down, cc.plt.down, kegg.plt.down, react.plt.down, labels=c( 'A', 'B', 'C', 'D' ), ncol=2, nrow=2 )

    jpeg( file = file.path( gprofile_figure_path, paste0( contrast, "_multi_down.jpeg" ) ), width=1280, height=1024, res=75 )
    print( multi.plt.down )
    dev.off()

    ##############
    # top n down #
    ##############
    topn.down.tbl <- full.down.tbl %>%
      filter( t_type != "tf" ) %>%
      top_n( ., wt=pvalue, n = -10 ) %>%
      arrange( desc( pvalue ) ) %>%
      mutate( title = factor( title, levels=title ) )

    num.down.topn <- dim( topn.down.tbl )[1]

    topn.down.plt <- ggplot( topn.down.tbl, aes( y = -log10( pvalue ), x = title ) ) +
      theme_classic() +
      gg_bigger_texts +
      gg_center_title +
      coord_flip() +
      geom_col() +
      geom_hline( yintercept = -log10( 0.05 ), linetype=2, colour="white" ) +
      ylab( parse( text="-log[10](P-value)" ) ) +
      xlab( "Pathway\n" ) +
      ggtitle( paste0( "Top ", num.down.topn, " Decreased Enrichments" ) ) +
      theme( axis.text = element_text( size=16 ) )

    jpeg( file = file.path( gprofile_figure_path, paste0( contrast, "_topn_down.jpeg" ) ), width=1280, height=1024, res=75 )
    print( topn.down.plt )
    dev.off()
  }

  if ( !is.null( dim( full.up.tbl ) ) ) {
    plot.up.tbl <- full.up.tbl %>%
      ddply( .variables = c( "t_type" ), .fun = top_n, wt = pvalue, n = -6 )

    ################################
    # individual increased pathway #
    ################################
    bp.plt.up <- filter( plot.up.tbl, t_type == "BP" ) %>%
      arrange( desc( pvalue ) ) %>%
      mutate( title = factor( title, levels=title ) ) %>%
      ggplot( aes( y = -log10( pvalue ), x = title ) ) +
      theme_classic() +
      gg_bigger_texts +
      gg_center_title +
      coord_flip() +
      geom_col() +
      geom_hline( yintercept = -log10( 0.05 ), linetype=2, colour="white"  ) +
      ylab( parse( text="-log[10](P-value)" ) ) +
      xlab( "Pathway\n" ) +
      ggtitle( "GO Biological Process" )

    cc.plt.up <- filter( plot.up.tbl, t_type == "CC" ) %>%
      arrange( desc( pvalue ) ) %>%
      mutate( title = factor( title, levels=title ) ) %>%
      ggplot( aes( y = -log10( pvalue ), x = title ) ) +
      theme_classic() +
      gg_bigger_texts +
      gg_center_title +
      coord_flip() +
      geom_col() +
      geom_hline( yintercept = -log10( 0.05 ), linetype=2, colour="white"  ) +
      ylab( parse( text="-log[10](P-value)" ) ) +
      xlab( "Pathway\n" ) +
      ggtitle( "GO Cellular Compartment" )

    kegg.plt.up  <- filter( plot.up.tbl, t_type == "keg" ) %>%
      arrange( desc( pvalue ) ) %>%
      mutate( title = factor( title, levels=title ) ) %>%
      ggplot( aes( y = -log10( pvalue ), x = title ) ) +
      theme_classic() +
      gg_bigger_texts +
      gg_center_title +
      coord_flip() +
      geom_col() +
      geom_hline( yintercept = -log10( 0.05 ), linetype=2, colour="white"  ) +
      ylab( parse( text="-log[10](P-value)" ) ) +
      xlab( "Pathway\n" ) +
      ggtitle( "KEGG" )

    react.plt.up <- filter( plot.up.tbl, t_type == "rea" ) %>%
      arrange( desc( pvalue ) ) %>%
      mutate( title = factor( title, levels=title ) ) %>%
      ggplot( aes( y = -log10( pvalue ), x = title ) ) +
      theme_classic() +
      gg_bigger_texts +
      gg_center_title +
      coord_flip() +
      geom_col() +
      geom_hline( yintercept = -log10( 0.05 ), linetype=2, colour="white"  ) +
      ylab( parse( text="-log[10](P-value)" ) ) +
      xlab( "Pathway\n" ) +
      ggtitle( "Reactome" )

    jpeg( file = file.path( gprofile_figure_path, paste0( contrast, "_BP_up.jpeg" ) ), width=1280, height=1024, res=75 )
    print( bp.plt.up )
    dev.off()

    jpeg( file = file.path( gprofile_figure_path, paste0( contrast, "_CC_up.jpeg" ) ), width=1280, height=1024, res=75 )
    print( cc.plt.up )
    dev.off()

    jpeg( file = file.path( gprofile_figure_path, paste0( contrast, "_KEGG_up.jpeg" ) ), width=1280, height=1024, res=75 )
    print( kegg.plt.up )
    dev.off()

    jpeg( file = file.path( gprofile_figure_path, paste0( contrast, "_REAC_up.jpeg" ) ), width=1280, height=1024, res=75 )
    print( react.plt.up )
    dev.off()

    bp.plt.up <- bp.plt.up + gg_gprofile_quadplot

    cc.plt.up <- cc.plt.up + gg_gprofile_quadplot

    kegg.plt.up <- kegg.plt.up + gg_gprofile_quadplot

    react.plt.up <- react.plt.up + gg_gprofile_quadplot

    multi.plt.up <- plot_grid( bp.plt.up, cc.plt.up, kegg.plt.up, react.plt.up, labels=c( 'A', 'B', 'C', 'D' ), ncol=2, nrow=2 )

    jpeg( file = file.path( gprofile_figure_path, paste0( contrast, "_multi_up.jpeg" ) ), width=1280, height=1024, res=75 )
    print( multi.plt.up )
    dev.off()

    ##############
    # top n up #
    ##############
    topn.up.tbl <- full.up.tbl %>%
      filter( t_type != "tf" ) %>%
      top_n( ., wt=pvalue, n = -10 ) %>%
      arrange( desc( pvalue ) ) %>%
      mutate( title = factor( title, levels=title ) )

    num.up.topn <- dim( topn.up.tbl )[1]

    topn.up.plt <- ggplot( topn.up.tbl, aes( y = -log10( pvalue ), x = title ) ) +
      theme_classic() +
      gg_bigger_texts +
      gg_center_title +
      coord_flip() +
      geom_col() +
      geom_hline( yintercept = -log10( 0.05 ), linetype=2, colour="white" ) +
      ylab( parse( text="-log[10](P-value)" ) ) +
      xlab( "Pathway\n" ) +
      ggtitle( paste0( "Top ", num.up.topn, " Increased Enrichments" ) ) +
      theme( axis.text = element_text( size=16 ) )

    jpeg( file = file.path( gprofile_figure_path, paste0( contrast, "_topn_up.jpeg" ) ), width=1280, height=1024, res=75 )
    print( topn.up.plt )
    dev.off()
  }
}
