require(data.table)
require(ggplot2)
require(ggseqlogo)
require(ggforce)
require(gridExtra)
require(ggpubr)
require(optparse)
suppressPackageStartupMessages(require(ComplexHeatmap))
suppressPackageStartupMessages(require(circlize))
suppressPackageStartupMessages(require(MASS))
suppressPackageStartupMessages(require(patchwork))

options( error = traceback, nwarnings = 10000 )
# source( "/home/ahcorcha/repos/tools/vanilla_JAMS/src/Methyl_ChIP_ftns.R" ) # ahcorcha


############################################################################## #
################################ Read I/O. ################################### #
################################################################################
option_list = list(
  make_option(c("-e", "--experiment"), type="character",
              default="CTCF_HEK293_GSM2026781_smallest_in_R",
              help="Experiment ID", metavar="character"),

  make_option(c("-f", "--flanking"), type="integer", default=20,
              help="length of flanking sequence around the motif", metavar="character"),

  make_option(c("-i", "--input_dir"), type="character", metavar="character",
              default="./data/CTCF_demo/02_formatted_data/smallest_demo",
              help="Input directory with PFM, methylation counts etc ..."),

  make_option(c("-o", "--output_dir"), type="character",
              default="./data/CTCF_demo/03_output",
              help="", metavar="character"),
  
  make_option(c("-s", "--script_path"), type="character",
              default="/home/ahcorcha/repos/tools/vanilla_JAMS/src",
              help="", metavar="character")
  );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser); rm(option_list, opt_parser)

source( paste0( opt$script_path, "/Methyl_ChIP_ftns.R" ) )

### Create output directory
outdir <- paste0( opt$output_dir, "/JAMS_", opt$experiment, 
                  "_flanking_", opt$flanking, "bps" )

dir.create( outdir, recursive = T, showWarnings = FALSE )

################################################################################
### Redirect all stdout to a file
log_file <- file( paste0( outdir,"/", opt$experiment, "_log.txt" ), open = "wt" )

sink( log_file, append=TRUE, type = "output" )
sink( log_file, append=TRUE, type = "message" )


cat(paste0("\nExperiment: ", opt$experiment, "\n", "Flanking: ", opt$flanking, "\n",
           "Input dir: ", opt$input_dir, "\n", "Output dir: ", opt$output_dir, "\n") )


### Check if input directory exists.
if ( !dir.exists(opt$input_dir) ){
  cat("Input directory doesn't exist, check --input_dir.\n")
  } else{
  
  ############################################################################################################ #
  ########################################       Load and format data      ################################### #
  ##############################################################################################################
  cat(paste0("Loading data ... \n") )
  ## Load data from input directory.
  flanking <- opt$flanking
  input_root <- opt$input_dir
  
  methyl_chip_data <- load_data(opt$input_dir)

  cat( paste0( "motif length: ", methyl_chip_data$pfm_length, "\n" ) )
  
  ###################################### #
  ###    Check methylation coverage.   # #
  methyl <- methyl_chip_data[["methyl"]][,c(-1)]
  nonmethyl <- methyl_chip_data[["nonmethyl"]][,c(-1)]
  meth_coverage <- Reduce("+", list(methyl, nonmethyl) )
  reads_per_row <- rowSums(meth_coverage, na.rm = T)
  covered_peaks <- names(reads_per_row[reads_per_row >= 1])
  
  methyl <- methyl_chip_data[["methyl"]][,c(-1)]
  nonmethyl <- methyl_chip_data[["nonmethyl"]][,c(-1)]
  meth_coverage <- Reduce("+", list(methyl, nonmethyl) )
  reads_per_row <- rowSums(meth_coverage, na.rm = T)
  covered_peaks <- names(reads_per_row[reads_per_row >= 0])
  
  methyl_chip_data$acc <- methyl_chip_data$acc[ methyl_chip_data$acc$Name %in% covered_peaks , ]
  methyl_chip_data$affimx <- methyl_chip_data$affimx[ methyl_chip_data$affimx$Name %in% covered_peaks , ]
  methyl_chip_data$target <- methyl_chip_data$target[ methyl_chip_data$target$Name %in% covered_peaks , ]
  methyl_chip_data$seq <- methyl_chip_data$seq[ methyl_chip_data$seq$Name %in% covered_peaks , ]
  methyl_chip_data$CpGs <- methyl_chip_data$CpGs[ methyl_chip_data$CpGs$Name %in% covered_peaks , ]
  methyl_chip_data$methyl <- methyl_chip_data$methyl[ methyl_chip_data$methyl$Name %in% covered_peaks , ]
  methyl_chip_data$nonmethyl <- methyl_chip_data$nonmethyl[ methyl_chip_data$nonmethyl$Name %in% covered_peaks , ]
  
  rm(methyl, nonmethyl, meth_coverage)
  
  ##############################################################################################################
  
  
  cat(paste0("Num. of Peaks: ", dim(methyl_chip_data$target)[1], "\n" ) ) 
  
  ############################################################################################################ #
  ########################################          Length vs tags         ################################### #
  ##############################################################################################################
  ### Plot Tags vs length
  cat(paste0("Ploting: tags vs length ... \n"))
  plot_target <- ggplot( data = methyl_chip_data$target, aes(x = length, y = tags) ) + 
    stat_density2d(aes(fill = ..density..^0.25), 
                   geom = "tile", contour = FALSE, n = 200) +
    scale_fill_continuous(low = "white", high = "dodgerblue4") + 
    xlab("length") + ylab("tags") + 
    ggtitle("length vs tags of ChIP-seq peaks") +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(filename = paste0(outdir, "/", opt$experiment, "_tags_vs_lenght.pdf"), 
         plot = plot_target, height = 7, width = 10)
  rm(plot_target)
  ##############################################################################################################
  
  ############################################################################################################ #
  ######################################## Plot: DNA acc. by bins vs tags  ################################### #
  ##############################################################################################################
  ### Plot log of dna accessibility by bins
  dna_acc_by_bins <- calculate_mean_by_bins( acc = methyl_chip_data$acc, 
                                             pfm_length = methyl_chip_data$pfm_length )
  
  acc_motif_plot_name <- paste0( outdir, "/", opt$experiment, "_DNA_accessibility_around_motif.pdf" )
  plot_dna_acc_around_motif( dna_acc_by_bins = dna_acc_by_bins, methyl_chip_data = methyl_chip_data,
                             acc_motif_plot_name = acc_motif_plot_name )
  
  ############################################################################################################ #
  #####################################    dna_acc_vs_tags   ################################################# #
  ht_acc_plot <- paste0( outdir, "/", opt$experiment, "_dna_acc_vs_tags_heatmap.pdf" )
  plot_heatmap_dna_acc_vs_target( methyl_chip_data = methyl_chip_data, 
                                  dna_acc_by_bins = dna_acc_by_bins, 
                                  ht_acc_plot = ht_acc_plot )
  
  ##############################################################################################################
  
  
  ############################################################################################################ #
  ########################################       Sequence vs tags          ################################### #
  ##############################################################################################################
  ht_seq_plot <- paste0( outdir, "/", opt$experiment, "_sequence_vs_tags_heatmap.pdf" )
  plot_heatmap_sequence_vs_target( methyl_chip_data = methyl_chip_data, ht_seq_plot = ht_seq_plot )
  
  ##############################################################################################################
  
  
  ############################################################################################################ #
  ########################################       Methylation vs tags          ################################ #
  ##############################################################################################################
  ht_methyl_plot <- paste0( outdir, "/", opt$experiment, "_meth_vs_tags_heatmap_max_coverage.pdf" )
  plot_heatmap_methyl_vs_target( methyl_chip_data, ht_methyl_plot = ht_methyl_plot, show_max_depth = T )
  
  ht_methyl_plot <- paste0( outdir, "/", opt$experiment, "_meth_vs_tags_heatmap_sensitive_coverage.pdf" )
  plot_heatmap_methyl_vs_target( methyl_chip_data, ht_methyl_plot = ht_methyl_plot, show_max_depth = F )
  
  
  ######################## Methylation Stats
  cat(paste0("Ploting: methylation stats ...  \n"))
  plot_methylation_stats( methyl_chip_data = methyl_chip_data, 
                          pseudocount = 0.0001, cutoff = 10,
                          outdir = outdir, experiment = opt$experiment )
  
  ##############################################################################################################
  
  
  ############################################################################################################ #
  ########################################           CpG-only              ################################### #
  ##############################################################################################################
  methyl_chip_data$target$tags <- as.factor(methyl_chip_data$target$tags)
  methyl_chip_data$target$pulldown.tag.old <- as.factor(methyl_chip_data$target$pulldown.tag.old)
  methyl_chip_data$target$ctrl.tag.old <- as.factor(methyl_chip_data$target$ctrl.tag.old)
  
  ### Fit model
  cat( paste0( "CpG model ...  \nAIC of null model: ", 
               glm( cbind( methyl_chip_data$target$pulldown.tag.old,
                           methyl_chip_data$target$ctrl.tag.old) ~ 1, 
                    family = binomial( link = "logit" ) )$aic, 
               "\nBuilding CpG-only model ...\n" ) )
  
  fit_CpG_only <- sequence_fit( methyl_chip_data$seq, methyl_chip_data$CpGs, 
                                methyl_chip_data$methyl, methyl_chip_data$nonmethyl, 
                                methyl_chip_data$acc, methyl_chip_data$target, 
                                methyl_chip_data$pfm_length, flanking = opt$flanking,
                                CpG_only = T )
  
  num_rank_warnings <- sum( grepl( pattern = "prediction from a rank-deficient fit may be misleading", 
                            x = names( warnings( ) ) ) )
  
  if( num_rank_warnings >= 1 ){
    
    cat("One or more warnings of type: prediction from a rank-deficient fit may be misleading\nDONE\n")
    sink( type = "output" )
    sink( type = "message" )
    close( log_file )
    
  } else {
    
    #################################### ##################################### #
    fit_CpG_only$t_train <- as.data.frame( fit_CpG_only$t_train )
    fit_CpG_only$c_train <- as.data.frame( fit_CpG_only$c_train )
    
    rownames( fit_CpG_only$t_train ) <- rownames( fit_CpG_only$XX_train )
    rownames( fit_CpG_only$c_train ) <- rownames( fit_CpG_only$XX_train )
    
    write.table( x = fit_CpG_only$XX_train, sep = "\t", quote = FALSE, 
                 row.names = TRUE, col.names = TRUE,
                 file = paste0( outdir, "/", opt$experiment, "_XX_train.tab" ) )
    
    write.table( x = fit_CpG_only$t_train, sep = "\t", quote = FALSE, 
                 row.names = TRUE, col.names = TRUE,
                 file = paste0( outdir, "/", opt$experiment, "_t_train.tab" ) )
    
    write.table( x = fit_CpG_only$c_train, sep = "\t", quote = FALSE, 
                 row.names = TRUE, col.names = TRUE,
                 file = paste0( outdir, "/", opt$experiment, "_c_train.tab" ) )
    ##################################### ###################################### #
    
    
    cat(paste0("Mean of correlations: ", round( mean(fit_CpG_only$correlation), digits = 6),
               ", sd: ", round( sd(fit_CpG_only$correlation), digits = 4 ),
               ", se: ", round( se(fit_CpG_only$correlation), digits = 4 ),
               "\nAIC: " , round( fit_CpG_only$fit$aic, digits = 4 ), "\n") )
    
    ## create a jpeg file with the scatter plot of predicted vs. observed binding strength 
    pdf(file = paste0(outdir, "/", opt$experiment, "_mCpG_only_model_total_cv_scatterplot.pdf"), 
        width = 7,height = 7)
    
    smoothScatter(fit_CpG_only$all_c_real, fit_CpG_only$all_c_predicted,
                  xlab = "log( c_pdwn_observed / c_ctrl_observed )", 
                  ylab = "log( c_pdwn_predicted / c_ctrl_predicted )", 
                  main = paste0(opt$experiment, "\n CpG model 10-fold crossvalidation in ", 
                                round(cor(fit_CpG_only$all_c_predicted, 
                                          fit_CpG_only$all_c_real), 4), 
                                ", n = ", length( fit_CpG_only$all_c_predicted ) ) )
    null_message <- dev.off()
    
    cat ( paste0( "Correlation between model predictions and binding strength cv: ",
                  round( cor(fit_CpG_only$all_c_predicted, 
                             fit_CpG_only$all_c_real), digits = 5 ), "\n" ) )
    
    ## Write coefficients.
    write.table( summary(fit_CpG_only$fit)[["coefficients"]], 
                 file = paste0(outdir, "/", opt$experiment, "_mCpG_only_model_coefficients.txt"),
                 sep = "\t", row.names = T, col.names = NA, quote = F )
    
    
    label <- paste0(opt$experiment, "_mCpG_only_model")

    # ## ahcorcha ############################################################################################  ##
    model.mCpG_only <- write.sequence.model.av.met( fit_CpG_only, outdir, label )
    
    cat( paste0( "Correlation between control model and initial motif: ",
                 round( cor( methyl_chip_data$rcade_pwm,
                             as.vector( t( model.mCpG_only$motif_control[1:4,]) ) ), digits = 6 ), "\n\n" ) )
    
    
    cat( paste0( "Correlation between pulldown model and initial motif: ",
                 round( cor( methyl_chip_data$rcade_pwm,
                             as.vector( t( model.mCpG_only$motif_pulldown[1:4,]) ) ), digits = 6 ), "\n\n" ) )
    
    dna_acc_plot_name <- paste0( outdir, "/", opt$experiment, "_dna_coefficients.pdf" )
    plot_dna_acc_coefficients( fit_CpG_only$fit, dna_acc_plot_name )
    
    
    ## Save model
    save(fit_CpG_only, file = paste0(outdir, "/", opt$experiment, "_fit_CpG_only_model.rda"))
    
    
    
    
    
    # source( "./src/Methyl_ChIP_ftns.R" )
    
    paste0( outdir, "/", label, "_coefficients_with_FDR.txt")
    
    write.TF.background.coeffs.sep( fit_CpG_only, outdir, label )
     
    write.table( x =  coefs, file = paste0( outdir, "/", label, "_coefficients_with_FDR.txt"),
                 quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE )
    
    
    
    
    #######################################################################################################  ###
    methyl_chip_data_name <- paste0( outdir, "/", opt$experiment, "_methyl_chip_data.rda" )
    saveRDS( methyl_chip_data, file = methyl_chip_data_name )
    
    dna_acc_by_bins_name <- paste0( outdir, "/", opt$experiment, "_dna_acc_by_bins.rda" )
    saveRDS( dna_acc_by_bins, file = dna_acc_by_bins_name )
    
    fit_CpG_only_x.Met_name <- paste0( outdir, "/", opt$experiment, "_fit_CpG_only_x.Met.rda" )
    saveRDS( fit_CpG_only$x.Met, file = fit_CpG_only_x.Met_name  )
    
    
    motif_pulldown_name <- paste0( outdir, "/", opt$experiment, "_motif_pulldown.rda" )
    saveRDS( model.mCpG_only$motif_pulldown, file = motif_pulldown_name  )
    
    
    
    
    #######################################################################  ###    

    all_dat <- prepare_dat( methyl_chip_data$seq, methyl_chip_data$CpGs, 
                            methyl_chip_data$methyl, methyl_chip_data$nonmethyl, 
                            methyl_chip_data$acc, methyl_chip_data$target, 
                            methyl_chip_data$pfm_length, flanking = opt$flanking,
                            CpG_only = T )
    
    
    X_data_name <- paste0( outdir, "/", opt$experiment, "_predictors_response.rda" )
    saveRDS( all_dat, file = X_data_name )
    
    
    # heatmap_meth_affinity( all_dat = all_dat,
    #                        fit_object = fit_CpG_only, 
    #                        acc_data = dna_acc_by_bins,
    #                        outdir = outdir, 
    #                        experiment_name = opt$experiment )
    
    
    
    # top_acc <- 0.25
    # out_ht_path <- paste0( outdir, "/", opt$experiment, "_heatmap_all_data_top_",  gsub("\\.","_", top_acc), ".pdf" )
    # 
    # create_hm_all_data( methyl_chip_data = methyl_chip_data,
    #                     meth_data = fit_CpG_only$x.Met,
    #                     acc_data = dna_acc_by_bins,
    #                     pulldown_motif =, model.mCpG_only$motif_pulldown
    #                     out_ht_path = out_ht_path,
    #                     experiment =  opt$experiment,
    #                     top_acc = top_acc )
  
    cat("DONE\n")
    sink( type = "output" )
    sink( type = "message" )
    close( log_file )
  
  }
}

