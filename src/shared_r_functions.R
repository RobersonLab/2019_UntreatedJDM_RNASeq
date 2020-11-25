####################
# libraries needed #
####################
library( here )
library( tidyverse )
library( ggrepel )

#############
# constants #
#############
min_reads_detection_threshold <- 20
min_detected_samples_for_de <- 3

###############
# directories #
###############
dir.create( here( "results" ), showWarnings = FALSE )
dir.create( here( "results", "figures" ), showWarnings = FALSE )
dir.create( here( "results", "r_objects" ), showWarnings = FALSE )
dir.create( here( "results", "wgcna" ), showWarnings = FALSE )
dir.create( here( "results", "figures", "wgcna" ), showWarnings = FALSE )
dir.create( here( "results", "gprofiler" ), showWarnings = FALSE )
dir.create( here( "results", "figures", "gprofiler" ), showWarnings = FALSE )
dir.create( here( "results", "figures", "individual_gene_plots" ), showWarnings = FALSE )

figure_dir <- here( "results", "figures" )

##############################
# Download data for analysis #
##############################
file_link_map <- list(
'read_group_info' = list( 'https://ndownloader.figshare.com/files/24856892', here( "info", "JDM_RNASeq_ReadGroup_Info.csv" ) ),
'gene_counts' = list( 'https://ndownloader.figshare.com/files/24857000', here( "data", "rgid_unspliced_saf.csv.gz" ) ),
'das_scores' = list( 'https://ndownloader.figshare.com/articles/8150702', here( 'data', "JDM_RNASeq_DAS_Info.csv" ) ),
'rtqpcr' = list ( 'https://ndownloader.figshare.com/files/15188327', here( 'data', 'OtOf_RTqPCR_tidy.csv' ) )
)

###############################
# color blind friendly colors #
###############################
# from colorbrewer2.org       #
###############################
colorBlindPalette <- c( "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7" )

blues5Palette <- c( '#ffffcc', '#a1dab4', '#41b6c4', '#2c7fb8', '#253494' )
greens5Palette <- c( '#ffffcc', '#c2e699', '#78c679', '#31a354', '#006837' )
purples5Palette <- c( '#feebe2', '#fbb4b9', '#f768a1', '#c51b8a', '#7a0177' )

reds5Palette <- c( '#ffffb2', '#f3cc5c', '#fd8d3c', '#f03b20', '#bd0026' )
reds4Palette <- c( '#fef0d9', '#fdcc8a', '#fc8d59', '#d7301f' )
reds3Palette <- c( '#fee0d2', '#fc9272', '#de2d26' )

colorBlindPrintSafeDivergent1 <- c( "#a6611a", "#dfc27d", "#f5f5f5", "#80cdc1", "#018571" )
colorBlindPrintSafeDivergent2 <- c( "#7b3294", "#c2a5cf", "#f7f7f7", "#a6dba0", "#008837" )
colorBlindPrintSafeDivergent3 <- c( "#ca0020", "#f4a582", "#f7f7f7", "#92c5de", "#0571b0" )
colorBlindPrintSafeDivergent4 <- c( "#d7191c", "#fdae61", "#ffffbf", "#abd9e9", "#2c7bb6" )

colorBlindCopySafeSequential1 <- c( "#f0f0f0", "#bdbdbd", "#636363" )
colorBlindCopySafeSequential2 <- c( "#fef0d9", "#fdcc8a", "#fc8d59", "#d7301f" )

####################
# ggplot modifiers #
####################
gg_bigger_texts = theme(
  axis.title = element_text( size = 22 ),
  axis.text = element_text( size = 20 ),
  legend.text = element_text( size = 14 ),
  legend.title = element_text( size = 15 ),
  plot.title = element_text( size = 22 ),
  strip.text.x = element_text( size = 17, margin = margin( b = 5, t = 5 ) ),
  strip.text.y = element_text( size = 15 )
)

gg_multiplot_texts = theme(
  axis.title = element_text( size = 20 ),
  axis.text = element_text( size = 18 ),
  legend.text = element_text( size = 12 ),
  legend.title = element_text( size = 13 ),
  plot.title = element_text( size = 20 ),
  strip.text.x = element_text( size = 16, margin = margin( b = 5, t = 5 ) ),
  strip.text.y = element_text( size = 15 )
)

gg_reduce_pathway_text = theme(
  axis.title = element_text( size=14 ),
  axis.text.y = element_text( size=8 ),
  axis.text.x = element_text( size=10 ),
  plot.title = element_text( size=15 )
)

gg_no_legend = theme(
  legend.position='none'
)

gg_no_grid = theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)

gg_no_x_grid = theme(
  panel.grid.major.x = element_blank() )

gg_no_y_grid = theme(
  panel.grid.major.y = element_blank() )

gg_center_title = theme(
  plot.title = element_text( hjust = 0.5 )
)

gg_no_x_label = theme(
  axis.title.x = element_blank()
)

gg_no_y_label = theme(
  axis.title.y = element_blank()
)

gg_angled_x_text = theme (
  axis.text.x = element_text( angle = 45, vjust = 1, hjust = 1, color = 'black' )
)

###########################
# Detect values in matrix #
###########################
detect <- function( input_vector, cutoff, upper = TRUE )
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
    mutate( FoldChange = 2^log2FoldChange ) %>%
    mutate( FoldChange = case_when(
      FoldChange < 1.0 ~ -1 / FoldChange,
      TRUE ~ FoldChange )
    ) %>%
    dplyr::rename( pval = pvalue ) %>%
    dplyr::rename( qval = padj ) %>%
    select( gene_id, FoldChange, pval, qval, baseMean, log2FoldChange, lfcSE, stat ) %>%
    merge( annotation_df, ., by="gene_id" ) %>%
    mutate_at( ., .vars = c( 'FoldChange', 'pval', 'qval', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat' ), .funs = as.double ) %>%
    arrange( qval, pval, gene_id ) %>%
    return( . )
}

##################################
# Report differential expression #
##################################
summarize_pretty_de_results <- function( pretty_de, fc_cutoff ) {
  down <- pretty_de %>%
    filter( qval < 0.05 & FoldChange < 1.0 ) %>%
    pull( gene_id ) %>%
    length( . )

  down_fc <- pretty_de %>%
    filter( qval < 0.05 & FoldChange < ( -1 * fc_cutoff ) ) %>%
    pull( gene_id ) %>%
    length( . )

  up <- pretty_de %>%
    filter( qval < 0.05 & FoldChange > 1.0 ) %>%
    pull( gene_id ) %>%
    length( . )

  up_fc <- pretty_de %>%
    filter( qval < 0.05 & FoldChange > fc_cutoff ) %>%
    pull( gene_id ) %>%
    length( . )

  data.frame( Change = c( 'Down', paste0( 'Down - min ', sprintf( fmt="%.2f", fc_cutoff ) ), 'Up', paste0( 'Up - min ', sprintf( fmt="%.2f", fc_cutoff ) ) ), n = c( down, down_fc, up, up_fc ) ) %>%
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
      qval < qvalue_cutoff & abs( FoldChange ) >= fold_change_cutoff ~ "Large",
      qval < qvalue_cutoff ~ "Modest",
      TRUE ~ "Not" ) ) %>%
    mutate( Significant = factor( Significant, levels=c( "Not", "Modest", "Large" ) ) )

  ################################
  # set values for square x axis #
  ################################
  x_volcano_value <- ( abs( deg_dataframe$log2FoldChange[ is.finite( deg_dataframe$log2FoldChange ) ] ) + 0.051 ) %>%
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
  plot_volcano <- ggplot( deg_dataframe, aes( x=log2FoldChange, y=-log10( qval ), colour=Significant ) ) +
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

#################################################
# process raw gprofile into an updated gprofile #
#################################################
convertRefseqList2SymbolList <- function( input_string, name_frame, string_sep = ",", out_sep=";" ) {
  str_split( string = input_string, pattern = string_sep ) %>%
    unlist( "," ) %>%
    map_chr( .x = ., ~ name_frame[ .x, 'symbol' ] ) %>%
    sort( . ) %>%
    paste( ., collapse = out_sep ) %>%
    return( . )
}

parse_raw_gprofile <- function( gprof_file_path, symbol_df, title.truncate = 65, title.wrap = 35 ) {
  new_filename_string <- str_replace( gprof_file_path, pattern = ".csv", "_processed.tsv" )

  this_tibble <- read_csv( file = gprof_file_path, comment="#" ) %>%
    mutate( name = str_replace( term_name, "^ +", "" ) ) %>%
    mutate( symbols = map( intersections, convertRefseqList2SymbolList, name_frame = symbol_df ) ) %>%
    unnest( cols = c( symbols ) ) %>%
    mutate( title = paste0( term_id, ' - ', term_name ) ) %>%
    mutate( title = str_replace_all( string = title, pattern = "\\; +match class\\: +[0-9]+", replacement = "" ) ) %>%
    mutate( title = str_trunc( title, width = title.truncate ) ) %>%
    mutate( title = str_wrap( title, width = title.wrap ) ) %>%
    select( adjusted_p_value, source, term_id, term_name, symbols, intersections, title ) %>%
    mutate( intersections = str_replace_all( string = intersections, pattern = ",", replacement = ";" ) ) %>%
    arrange( adjusted_p_value )

  select( this_tibble, -title ) %>%
    write_tsv( path = new_filename_string )

  return( this_tibble )
}

###################################
# make lists of up and down genes #
###################################
direction_gene_ids_from_de <- function( deg_df, qval_cutoff = 0.05 ) {
  # gene_id lists
  tested_genes_gid <- pull( deg_df, gene_id )

  deg_gid <- filter( deg_df, qval < qval_cutoff ) %>%
    pull( gene_id )

  down_gid <- filter( deg_df, qval < qval_cutoff & FoldChange < 0.0 ) %>%
    pull( gene_id )

  up_gid <- filter( deg_df, qval < qval_cutoff & FoldChange > 0.0 ) %>%
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
  up2improved <- intersect( left_direction_list[[ 'up_gene_ids' ]],
                            right_direction_list[[ 'down_gene_ids' ]] ) %>%
    length( . )

  up2worsened <- intersect( left_direction_list[[ 'up_gene_ids' ]],
                            right_direction_list[[ 'up_gene_ids' ]] ) %>%
    length( . )

  up2unchanged <- setdiff( left_direction_list[[ 'up_gene_ids' ]],
                           right_direction_list[[ 'deg_gid' ]] ) %>%
    length( . )

  down2improved <- intersect( left_direction_list[[ 'down_gene_ids' ]],
                              right_direction_list[[ 'up_gene_ids' ]] ) %>%
    length( . )

  down2worsened <- intersect( left_direction_list[[ 'down_gene_ids' ]],
                              right_direction_list[[ 'down_gene_ids' ]] ) %>%
    length( . )

  down2unchanged <- setdiff( left_direction_list[[ 'down_gene_ids' ]],
                             right_direction_list[[ 'deg_gid' ]] ) %>%
    length( . )

  unchanged2up <- setdiff( right_direction_list[[ 'up_gene_ids' ]],
                           left_direction_list[[ 'deg_gid' ]] ) %>%
    length( . )

  unchanged2down <- setdiff( right_direction_list[[ 'down_gene_ids' ]],
                             left_direction_list[[ 'deg_gid' ]] ) %>%
    length( . )

  unchanged2unchanged <- union( left_direction_list[[ 'all_tested_gid' ]],
                                right_direction_list[[ 'all_tested_gid' ]] ) %>%
    setdiff( ., union( left_direction_list[[ 'deg_gid' ]],
                       right_direction_list[[ 'deg_gid' ]] ) ) %>%
    length( . )

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
