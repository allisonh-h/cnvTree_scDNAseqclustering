#' Perform pqArm clustering step in scDNA-seq data clustering Workflow
#'
#' This function executes the pqArm clustering step in the single-cell DNA  
#' sequencing (scDNA-seq) data clustering workflow, grouping cells based on  
#' chromosomal arm-level copy number variations.
#'
#' @param input A named list where each element is a `GRanges` object representing 
#'  a single cell.
#' @param pqArm_file In-build cytoband template for selection: `hg38`, `hg19`, 
#'  `mm10`, `mm39`. Or a filepath of a table for cytoband information seen on  
#'  Giemsa-stained chromosomes. It should include the following columns:
#'  
#'   - `chrom`: Reference sequence chromosome or scaffold.
#'   - `chromStart`: Start position in genoSeq.
#'   - `chromEnd`: End position in genoSeq.
#'   - `name`: Name of cytogenetic band.
#'   - `gieStain`: Giemsa stain results.
#'   
#' @param cluster The assigned cluster label (k = 2)
#'
#' @return A data frame recording the pqArm clustering results for each cell, 
#'  tracking the clustering history at this stage.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'  file_path <- system.file("extdata", "example_data.rds", package = "cnvTree")
#'  config <- read_yaml(config_path)
#'  input_df <- changeFormat(input_dir_DNA = config$input_dir_DNA, cores = 2)
#'  pqArm_result <- NEW_pqArmClustering(input = input_df, 
#'                                      pqArm_file = config$pqArm_file) 
#' }
#'
NEW_pqArmClustering <- function(input, pqArm_file, cluster, sexchromosome)
{
  # Locked variable
  config_path_hid <- system.file("extdata", "cnvTree_config_hid.yaml", package = "cnvTree")
    # mark to remove

  config_hid <- read_yaml(config_path_hid)
  cluster = config_hid$cluster
  sexchromosome = config_hid$sexchromosome
  
  fucStep <- paste0(" 7.1_cnvTree_", config_hid$v_num)
  DebugMsg(fucStep, "start", msg = config_hid$msg)
  
  message("=== Step 02: pqArm Clustering ===")

  Clustering_output <- data.frame(cluster = cluster, cell = names(input))

  pqArm_result <- NULL
  for(Label in unique(Clustering_output$cluster)) {
    ptm <- startTimed("pqArm Clustering ...")
    Smooth_pqCN <- NEW_pqArm_CN(input = input, 
                                Cluster_label = Label, 
                                Clustering_output = Clustering_output, 
                                pqArm_file = pqArm_file,
                                sexchromosome = sexchromosome)
    pqArm_cluster <- pqArm_clustering(matrix = Smooth_pqCN, Label = Label)
    pqArm_cluster_summary <- NEW_pqArm_clustering_summary(matrix = pqArm_cluster, 
                                                          Label = Label)
    #LH_02102025: changed to NEW_pqArm_clustering_summary()
    pqArm_cluster <- merge(pqArm_cluster, pqArm_cluster_summary)
    pqArm_result <- rbind(pqArm_result, pqArm_cluster)
    endTimed(ptm)
  }
  
  DebugMsg(fucStep, "end", msg = config_hid$msg)
  
  return(pqArm_result)
}


#' Reclustering method. Perform Cluster Consolidation in scDNA-seq data clustering 
#' workflow
#'
#' This function executes the Cluster Consolidation step in the single-cell DNA 
#' sequencing (scDNA-seq) data clustering workflow, consolidating clusters based 
#' on differences in chromosomal arm-level copy number variations.
#'
#' @param input A named list where each element is a `GRanges` object representing 
#'  a single cell.
#' @param pqArm_output A data frame recording the pqArm clustering results for 
#'  each cell. This table tracks the clustering history up to this stage.
#' @param pqArm_file In-build cytoband template for selection: `hg38`, `hg19`, 
#'  `mm10`, `mm39`. Or a filepath of a table for cytoband information seen on 
#'  Giemsa-stained chromosomes. It should include the following columns:
#'  
#'   - `chrom`: Reference sequence chromosome or scaffold.
#'   - `chromStart`: Start position in genoSeq.
#'   - `chromEnd`: End position in genoSeq.
#'   - `name`: Name of cytogenetic band.
#'   - `gieStain`: Giemsa stain results.
#'   
#' @param difratio_chr A numeric value defining the threshold for acceptable 
#'  difference ratios across different chromosomes.
#'
#' @return A data frame recording both pqArm clustering and consolidating results 
#'  for each cell, maintaining the clustering history across steps.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'  file_path <- system.file("extdata", "example_data.rds", package = "cnvTree")
#'  config <- read_yaml(config_path)
#'  input_df <- changeFormat(input_dir_DNA = config$input_dir_DNA, cores = 2)
#'  pqArm_result <- NEW_pqArmClustering(input = input_df, 
#'                                     pqArm_file = config$pqArm_file) 
#'  Consolidation_result <- clusterConsolidation(input = input_df, 
#'                                               pqArm_output = pqArm_result, 
#'                                               pqArm_file = config$pqArm_file)
#' }
#'
clusterConsolidation <- function(input, pqArm_output, pqArm_file, difratio_chr)
{
  config_path_hid <- system.file("extdata", "cnvTree_config_hid.yaml", package = "cnvTree")
  config_hid <- read_yaml(config_path_hid)
  
  fucStep <- paste0(" 7.2_cnvTree_", config_hid$v_num)
  DebugMsg(fucStep, "start", msg = config_hid$msg)
  
  # Locked variable
  difratio_chr <- config_hid$difratio_chr
  
  message("=== Step 03: Cluster consolidation ===")
  
  Recluster_Output <- NULL
  for(Label in unique(pqArm_output$cluster)) {
    ptm <- startTimed("Consolidating clusters to make cluster bigger ...")
    pqArm_merge_target <- NULL
    New_pqArm_clustering <- pqArm_recluster(pqArm_cluster = pqArm_output, 
                                            Cluster = Label)

    # potential merge target must more than one (with less10 clusters( >=2 & <10), 
    # or more than one more10 clusters(>2))
    if (length(New_pqArm_clustering) > 1) {
      pqArm_dif <- pqArm_reclustering_dif(
                      input = input, 
                      pqArm_recluster_sim = New_pqArm_clustering,
                      pqArm_cluster = pqArm_output, 
                      Cluster = Label, 
                      pqArm_file = pqArm_file)
      
      pqArm_merge_target <- pqArm_reclusterBy_ratio_target(
                                pqArm_cluster = pqArm_output,
                                Cluster = Label, 
                                pqReclsut_sim = pqArm_dif, 
                                difratio_chr = difratio_chr)
    }

    # is.null(pqArm_merge_target)
    if (is.null(pqArm_merge_target) == FALSE) {
      Recluster_CellID <- pqArm_recluster_result(pqArm_cluster = pqArm_output, 
                                                 pqReclsut_target = pqArm_merge_target)
      Recluster_Output <- rbind(Recluster_Output, Recluster_CellID)
    } else {
      Recluster_CellID <- pqArm_output %>%
        dplyr::filter(.data$cluster %in% Label) %>%
        dplyr::mutate(Recluster_pattern = .data$pqArm_pattern,
                      Recluster_cellnum = .data$pqArm_cellnum,
                      Recluster_cluster = .data$pqArm_cluster)
      Recluster_Output <- rbind(Recluster_Output, Recluster_CellID)
    }
    endTimed(ptm)
  }

  DebugMsg(fucStep, "end", msg = config_hid$msg)
  
  return(Recluster_Output)
}


#' Perform subclustering in scDNA-seq data clustering workflow
#'
#' This function executes the subclustering step in the single-cell DNA sequencing 
#' (scDNA-seq) data clustering workflow, refining clusters based on copy number 
#' variations (CNVs) within subpopulations of cells.
#'
#' @param input A named list where each element is a `GRanges` object representing 
#'  a single cell.
#' @param Consolidating_output A data frame recording pqArm clustering and consolidating 
#'  results for each cell. This table tracks the clustering history up to this stage.
#' @param min_cell An integer specifying the minimum number of cells required for 
#'  a cluster to be retained.
#' @param overlap_region An integer representing the genomic region where copy 
#'  number frequently changes.
#' @param dif_ratio A numeric value defining the tolerance threshold for copy 
#'  number differences between cells within a cluster.
#'
#' @return A list containing two data frames:
#'   - `final_cluster_output`: Records the clustering history from the pqArm, 
#'  consolidating, and subclustering steps.
#'   - `Subclone_CN`: Records each subclone's unique chromosome segment template 
#'   and its copy number. It includes the following columns:
#'     - `chr`: Chromosome name (chr1, chr2, ...).
#'     - `start`: Start position of the segment.
#'     - `end`: End position of the segment.
#'     - `region`: The defined region index
#'     - `Subclone`: Subclone identifier.
#'     - `CN`: Copy number value of the segment.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'  file_path <- system.file("extdata", "example_data.rds", package = "cnvTree")
#'  config <- read_yaml(config_path)
#'  input_df <- changeFormat(input_dir_DNA = config$input_dir_DNA, cores = 2)
#'  pqArm_result <- NEW_pqArmClustering(input = input_df, 
#'                                     pqArm_file = config$pqArm_file) 
#'  Consolidation_result <- clusterConsolidation(input = input_df, 
#'                                               pqArm_output = pqArm_result, 
#'                                               pqArm_file = config$pqArm_file)
#'  Subclone_output <- SubClustering(input = input_df,
#'                                   Consolidating_output = Consolidation_result)
#' }
#'
SubClustering <- function(input, Consolidating_output, min_cell, overlap_region, 
                          dif_ratio)
{
  config_path_hid <- system.file("extdata", "cnvTree_config_hid.yaml", package = "cnvTree")
  config_hid <- read_yaml(config_path_hid)
  fucStep <- paste0(" 7.3_cnvTree_", config_hid$v_num)
  DebugMsg(fucStep, "start", msg = config_hid$msg)
  
  # Locked variable
  min_cell = config_hid$min_cell
  overlap_region = config_hid$overlap_region
  dif_ratio = config_hid$dif_ratio
  
  message("=== Step 04: Subclustering ===")
  
  Final_output <- list()
  Subclone_cluster <- NULL
  Subclone_CNr <- NULL
  cell_clustering <- NULL

  # from overlap_region to number of bins in coverage
  binsize <- as.data.frame(input[[1]]$bins)$width[1]
  overlap_bp <- round((overlap_region/binsize), digits = 0)

  # Select enough cell number Reclusters, avoiding keep filtering
  R <- Consolidating_output %>% dplyr::filter(.data$Recluster_cellnum >= min_cell)
  
  for (Recluster_label in unique(R$Recluster_cluster)) {
    ptm <- startTimed("Subclustering for Recluster ", Recluster_label, " ... \n")
    breakpoints <- collect_cluster_bp(input = input, 
                                      Clustering_output = R, 
                                      Recluster_label = Recluster_label)
    Cell_num <- R %>% 
                dplyr::filter(.data$Recluster_cluster == Recluster_label) %>% 
                dplyr::pull(.data$Recluster_cellnum)
    bp <- output_bp_covers(Template = breakpoints, 
                           binsize = binsize, 
                           overlap = overlap_bp, 
                           overlap_times = Cell_num[1] * 0.5)#overlap_times = 0.5*cell
    if (is.null(bp)) { #LH: added 12282025
      next  #LH: added 12282025
    }  #LH: added 12282025
    
    event_region <- bp_events(input = input, Template = bp, binsize = binsize)
    consensus_bp_template <- bp_region(event = event_region, binsize = binsize)
    cell_CNregion <- Region_CN(input = input, 
                               Reclustering_output = Consolidating_output,
                               Recluster_label = Recluster_label, 
                               events = consensus_bp_template)

    if (is.null(cell_clustering) == TRUE) {
      cell_clustering <- Subclone_clustering(CN_incells_input = cell_CNregion, 
                                             event_region = consensus_bp_template,
                                             dif_ratio = dif_ratio, 
                                             Subclone_num = 0)
    } else {
      Subclone_num = max(Subclone_cluster$Subclone)
      cell_clustering <- Subclone_clustering(CN_incells_input= cell_CNregion, 
                                             event_region= consensus_bp_template,
                                             dif_ratio = dif_ratio, 
                                             Subclone_num = Subclone_num)
    }
    #endTimed(ptm)

    Subclone_cluster <- rbind(Subclone_cluster, cell_clustering)

    # Output: 1. Copy number in each region in each subclone
    s_CN <- Subclone_CNregion(sep_region = consensus_bp_template, 
                              CN_region = cell_CNregion, 
                              each_subclone = cell_clustering,
                              min_cell = min_cell, 
                              output = "SubcloneRegionCN")
    Subclone_CNr <- rbind(Subclone_CNr, s_CN)
  }

  colnames(Subclone_cluster) <- c("Subclone_cluster", "cellID", "Subclone_cellnum")
  Final_output <- list(final_cluster_output = dplyr::left_join(Consolidating_output, 
                                                               Subclone_cluster, 
                                                               by = "cellID"),
                                                               Subclone_CN = Subclone_CNr)

  DebugMsg(fucStep, "end", msg = config_hid$msg)
  
  return(Final_output)
}


#' Outputting scDNA clustering final outputs. Generate final output in scDNA-seq 
#' data clustering workflow
#'
#' This function produces the final output for the single-cell DNA sequencing (scDNA-seq)
#' clustering workflow, integrating clustering results and generating visualizations.
#'
#' @param input A named list where each element is a `GRanges` object representing 
#'  a single cell.
#' @param Summary A list containing two data frames:
#' 
#'   - `final_cluster_output`: Records the clustering history from the pqArm, 
#'   consolidating, and subclustering steps.
#'   - `Subclone_CN`: Records each subclone's unique chromosome segment template 
#'   and its copy number. It includes the following columns:
#'     - `chr`: Chromosome name (chr1, chr2, ...).
#'     - `start`: Start position of the segment.
#'     - `end`: End position of the segment.
#'     - `region`: The defined region index
#'     - `Subclone`: Subclone identifier.
#'     - `CN`: Copy number value of the segment.
#'     
#' @param pqArm_file In-build cytoband template for selection: `hg38`, `hg19`, 
#'  `mm10`, `mm39`. Or a filepath of a table for cytoband information seen on 
#'  Giemsa-stained chromosomes. It should include the following columns:
#'  
#'   - `chrom`: Reference sequence chromosome or scaffold.
#'   - `chromStart`: Start position in genoSeq.
#'   - `chromEnd`: End position in genoSeq.
#'   - `name`: Name of cytogenetic band.
#'   - `gieStain`: Giemsa stain results.
#'   
#' @param output_dir A character string specifying the directory where the output 
#'  files will be saved.
#' @param consecutive_region A numeric value defining the minimum length threshold 
#'  for filtering copy number variation (CNV) regions.
#' @param cellcutoff A numeric value specifying the minimum number of cells required 
#'  for a cluster to be retained.
#' @param sexchromosome A logical value. If `TRUE`, the output includes plots 
#'  with sex chromosome copy number information.
#' @param smoothheatmap A logical value. If `TRUE`, the output applies smoothing 
#'  over a 10⁶ bp range in chromosome copy number visualizations.
#'
#' @return The function generates final clustering results and visualizations, 
#'  saving them to the specified output directory.
#' 
#' @export
#'
#' @examples
#' \dontrun{
#'  file_path <- system.file("extdata", "example_data.rds", package = "cnvTree")
#'  config <- read_yaml(config_path)
#'  input_df <- changeFormat(input_dir_DNA = config$input_dir_DNA, cores = 2)
#'  pqArm_result <- NEW_pqArmClustering(input = input_df, 
#'                                     pqArm_file = config$pqArm_file) 
#'  Consolidation_result <- clusterConsolidation(input = input_df, 
#'                                               pqArm_output = pqArm_result, 
#'                                               pqArm_file = config$pqArm_file)
#'  Subclone_output <- SubClustering(input = input_df,
#'                                   Consolidating_output = Consolidation_result)
#'  scDNA_Output(input = input_df,
#'               Summary = Subclone_output,
#'               output_dir = config$output_dir,
#'               pqArm_file = config$pqArm_file
#'               cellcutoff = config$cellcutoff,
#'               smoothheatmap = config_hid$smoothheatmap)
#'  }
#'
scDNA_Output <- function(input, Summary, pqArm_file, output_dir, 
                         consecutive_region, cellcutoff,
                         sexchromosome, smoothheatmap)
{
  config_path_hid <- system.file("extdata", "cnvTree_config_hid.yaml", package = "cnvTree")
  config_hid <- read_yaml(config_path_hid)
  fucStep <- paste0(" 7.4_cnvTree_", config_hid$v_num)
  DebugMsg(fucStep, "start", msg = config_hid$msg)
  
  # Locked variable
  consecutive_region = config_hid$consecutive_region 
  sexchromosome = config_hid$sexchromosome
  
  message("=== Step 05: Output cnvTree results ===")
  timestamp <- format(Sys.time(), "%m%d_%H")
  
  # 1. cellID summary
  filename = paste0("/", timestamp, "_", "cnvTree.scDNAseq_grouping")
  writeOutput(data = Summary$final_cluster_output, 
              filename = filename, path = output_dir)
  message("Output /cnvTree.scDNAseq_grouping.txt is done.")

  # 2. Each region copy number to subclone
  filename = paste0("/", timestamp, "_", "cnvTree.scDNAseq_SubcloneRegionCN")
  writeOutput(data = Summary$Subclone_CN, 
              filename = filename, 
              path = output_dir)
  message("Output /cnvTree.scDNAseq_SubcloneRegionCN.txt is done.")

  # 3. All CNV region in the sample
  CNV_Data <- Total_cnvRegion(input = input, 
                              Template = Summary$Subclone_CN, 
                              pqArm_file = pqArm_file, 
                              consecutive_region = consecutive_region)
  
  if (nrow(CNV_Data) != 0){
    CNV_Data <- cnvRegion.toPQarm(FILE = CNV_Data, pqArm_file = pqArm_file)
    filename = paste0("/", timestamp, "_", "cnvTree.scDNAseq_DefinedCNVregion")
    writeOutput(data = CNV_Data, 
                filename = filename, 
                path = output_dir)
    message("Output /cnvTree.scDNAseq_DefinedCNVregion.txt is done.")
  } else {
    message("Warning: Based on the length of consecutive_region, 
             no CNV regions remain in the sample. \n
             Consider lowering the consecutive_region parameter. \n
             The file cnvTree.scDNAseq_DefinedCNVregion.txt and 
             cnvTree.scDNAseq_DefinedCNVregion.txt was not generated.")
  }
  # 4. DNA superimpose
  if (nrow(CNV_Data) != 0){
    DNA_superimpose <- scDNA.superimpose(Template = Summary, DefinedCNVs = CNV_Data)
    DNA_superimpose <- scDNA.clustering(Template = DNA_superimpose)
    filename = paste0("/", timestamp, "_", "cnvTree.scDNAseq_DNAcluster")
    writeOutput(data = DNA_superimpose$DNA_cluster, 
                filename = filename, 
                path = output_dir)
    message("Output /cnvTree.scDNAseq_DNAcluster.txt is done.")
  } else {
    message("Warning: Based on the length of consecutive_region, 
             no CNV regions remain in the sample. \n
             Consider lowering the consecutive_region parameter. \n
             The file cnvTree.scDNAseq_DNAcluster.txt and 
             cnvTree.scDNAseq_DefinedCNVregion.txt were not generated.")
  }

  # 5. Final cluster heatmap
  FILEname = paste0("/", timestamp, "_", "cnvTree.scDNAseq_Grouping_fig.pdf")
  Totalcluster_pdf(Input = input, 
                   Template = Summary$final_cluster_output,
                   pqArm_file = pqArm_file, 
                   cellcutoff = cellcutoff, 
                   step = "Subclone", 
                   sexchromosome = sexchromosome,
                   FILEpath = output_dir, 
                   FILEname = FILEname)
  message("Output /cnvTree.scDNAseq_Grouping_fig.pdf is done.")

  # 6. cluster heatmap with dendrogram
  FILEname = paste0("/", timestamp, "_", "cnvTree.scDNAseq_heatmap.png")
  NEW_scDNA_CNVpattern(input = input, 
                       final_cluster = Summary$final_cluster_output, 
                       cellcutoff = cellcutoff,
                       sexchromosome = sexchromosome, 
                       smoothheatmap = smoothheatmap,
                       pqArm_file = pqArm_file, 
                       FILEpath = output_dir, 
                       FILEname = FILEname)
  message("Output /cnvTree.scDNAseq_heatmap.png is done.")

  print(paste("LH: Result saved to", output_dir))
  
  return(DNA_superimpose)
  
  DebugMsg(fucStep, "end", msg = config_hid$msg)
}


#' One-Step scDNA-seq Cell Clustering Pipeline
#'
#' This function executes the entire single-cell DNA sequencing (scDNA-seq) clustering 
#' workflow in one step, including NEW_pqArm_clustering, consolidating, subclustering, 
#' and final CNV-based output generation.
#'
#' @param input_dir_DNA A named list where each element is a `GRanges` object 
#'  representing a single cell.
#' @param pqArm_file In-build cytoband template for selection: `hg38`, `hg19`, 
#'  `mm10`, `mm39`. Or a filepath of a table for cytoband information seen on 
#'  Giemsa-stained chromosomes. It should include the following columns:
#'
#'   - `chrom`: Reference sequence chromosome or scaffold.
#'   - `chromStart`: Start position in genoSeq.
#'   - `chromEnd`: End position in genoSeq.
#'   - `name`: Name of cytogenetic band.
#'   - `gieStain`: Giemsa stain results.
#'   
#' @param difratio_chr A numeric value defining the threshold for acceptable 
#'  difference ratios across different chromosomes.
#' @param output_dir A character string specifying the directory where output 
#'  files will be saved.
#' @param cellcutoff A numeric value specifying the minimum number of cells required 
#'  for a cluster to be retained in the final output.   
#' @param smoothheatmap A logical value. If `TRUE`, the final output applies 
#'  smoothing over a 10⁶ bp range in chromosome copy number visualizations.
#' @param min_cell An integer specifying the minimum cell count required for a 
#'  cluster to be included in the output.
#' @param overlap_region An integer representing genomic regions where copy number 
#'  frequently changes in the subclone clustering step.
#' @param dif_ratio A numeric value defining the tolerance threshold for copy 
#'     number differences between cells within a cluster.
#' @param consecutive_region A numeric value defining the minimum CNV region length 
#'.   threshold for filtering CNV events in the final output.
#' @param sexchromosome A logical value. If `TRUE`, the final output includes 
#'    plots with sex chromosome copy number information.
#'
#' @return The function performs complete clustering and CNV analysis, saving final 
#'  results and visualizations in the specified output directory.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   file_path <- system.file("extdata", "example_data.rds", package = "cnvTree")
#'   config <- read_yaml(config_path)
#'   selected_groups <- run_cnvTree_Pipeline(config = config, 
#'                                           output_dir = config$output_dir,
#'                                           input_dir_RNA = config$input_dir_RNA,
#'                                           RNAdataSource = config$RNAdataSource)
#'                            
#'   cnvTree_scDNAclustering(input_dir_DNA = config$input_dir_DNA,
#'                           output_dir = config$output_dir,
#'                           pqArm_file = config$pqArm_file, 
#'                           cellcutoff = config$cellcutoff,
#'                           smoothheatmap = config_hid$smoothheatmap)
#' }
#'
cnvTree_scDNAclustering <- function(input_dir_DNA, 
                                    pqArm_file, 
                                    output_dir,
                                    cellcutoff,
                                    smoothheatmap,
                                    difratio_chr,
                                    min_cell, 
                                    overlap_region, 
                                    dif_ratio,
                                    consecutive_region, 
                                    sexchromosome) 
{
  config_path_hid <- system.file("extdata", "cnvTree_config_hid.yaml", package = "cnvTree")
  config_hid <- read_yaml(config_path_hid)
  # Locked variable
  difratio_chr = config_hid$difratio_chr
  min_cell = config_hid$min_cell 
  overlap_region = config_hid$overlap_region
  dif_ratio = config_hid$dif_ratio
  consecutive_region = config_hid$consecutive_region
  smoothheatmap = config_hid$smoothheatmap
  sexchromosome = config_hid$sexchromosome
  
  # step 00: 
  input <- ProcessHmmList(input_dir_DNA = input_dir_DNA)
  
  ptm <- startTimed("Start cnvTree_scDNAclustering...")
  
  # step 01: pqArm clustering
  pqArm_result <- NEW_pqArmClustering(input = input, pqArm_file = pqArm_file)
  
  # step 02: Re-clustering
  Consolidation_result <- clusterConsolidation(input = input, 
                                               pqArm_output = pqArm_result, 
                                               pqArm_file = pqArm_file, 
                                               difratio_chr = difratio_chr)
  # step 03: Subclone clustering
  Subclone_output <- SubClustering(input = input,
                                   Consolidating_output = Consolidation_result,
                                   min_cell = min_cell, 
                                   overlap_region = overlap_region, 
                                   dif_ratio = dif_ratio)
  # step 04: Output
  scDNA_Output(input = input,
               Summary = Subclone_output,
               output_dir = output_dir,
               pqArm_file = pqArm_file,
               consecutive_region = consecutive_region,
               cellcutoff = cellcutoff, ## LH: inconsistent
               sexchromosome = sexchromosome,
               smoothheatmap = smoothheatmap)
  
  endTimed(ptm)
}


#' One-Step scDNA-seq Cell Clustering Pipeline
#'
#' This function executes the entire single-cell DNA sequencing (scDNA-seq) clustering 
#' workflow in one step, including pqArm clustering, consolidating, subclustering, 
#' and final CNV-based output generation.
#'
#' @param input_dir_DNA A named list where each element is a `GRanges` object 
#'   representing a single cell.
#' @param pqArm_file In-build cytoband template for selection: `hg38`, `hg19`, 
#'   `mm10`, `mm39`. Or a filepath of a table for cytoband information seen on 
#'   Giemsa-stained chromosomes. It should include the following columns:
#'
#'   - `chrom`: Reference sequence chromosome or scaffold.
#'   - `chromStart`: Start position in genoSeq.
#'   - `chromEnd`: End position in genoSeq.
#'   - `name`: Name of cytogenetic band.
#'   - `gieStain`: Giemsa stain results.
#'   
#' @param difratio_chr A numeric value defining the threshold for acceptable 
#'   difference ratios across different chromosomes.
#' @param min_cell An integer specifying the minimum cell count required for a 
#'   cluster to be included in the output.
#' @param overlap_region An integer representing genomic regions where copy number 
#'   frequently changes in the subclone clustering step.
#' @param dif_ratio A numeric value defining the tolerance threshold for copy 
#'   number differences between cells within a cluster.
#' @param output_dir A character string specifying the directory where output 
#'   files will be saved.
#' @param consecutive_region A numeric value defining the minimum CNV region length 
#'.  threshold for filtering CNV events in the final output.
#' @param cellcutoff A numeric value specifying the minimum number of cells required 
#'   for a cluster to be retained in the final output.
#' @param sexchromosome A logical value. If `TRUE`, the final output includes 
#'   plots with sex chromosome copy number information.
#' @param smoothheatmap A logical value. If `TRUE`, the final output applies 
#'   smoothing over a 10⁶ bp range in chromosome copy number visualizations.
#'
#' @return The function performs complete clustering and CNV analysis, saving final 
#'   results and visualizations in the specified output directory.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'#' file_path <- system.file("extdata", "example_data.rds", package = "cnvTree")
#' config <- read_yaml(config_path)
#' selected_groups <- run_cnvTree_Pipeline(config = config, 
#'                                         output_dir = config$output_dir,
#'                                         input_dir_RNA = config$input_dir_RNA,
#'                                         RNAdataSource = config$RNAdataSource)
#' cnvTree_scDNAclustering_df(input_dir_DNA = config$input_dir_DNA, 
#'                            output_dir = config$output_dir,
#'                            pqArm_file = config$pqArm_file, 
#'                            cellcutoff = config$cellcutoff,
#'                            smoothheatmap = config$smoothheatmap)
#'  }                          
#'
cnvTree_scDNAclustering_df <- function(input_dir_DNA, 
                                       pqArm_file, 
                                       output_dir,
                                       cellcutoff,
                                       smoothheatmap) 
{
  config_path_hid <- system.file("extdata", "cnvTree_config_hid.yaml", package = "cnvTree")
  config_hid <- read_yaml(config_path_hid)
  # Locked variable
  difratio_chr = config_hid$difratio_chr
  min_cell = config_hid$min_cell 
  overlap_region = config_hid$overlap_region
  dif_ratio = config_hid$dif_ratio
  consecutive_region = config_hid$consecutive_region 
  sexchromosome = config_hid$sexchromosome
  cores = config_hid$cores
  smoothheatmap = config_hid$smoothheatmap
  
  message("Start data transform...")
  # step 00: data transform
  input_df <- changeFormat(input_dir_DNA = input_dir_DNA)

  ptm <- startTimed("Start cnvTree_scDNAclustering...")
  
  # step 01: pqArm clustering
  pqArm_result <- NEW_pqArmClustering(input = input_df, pqArm_file = pqArm_file)
  
  # step 02: Re-clustering
  Consolidation_result <- clusterConsolidation(input = input_df, 
                                               pqArm_output = pqArm_result, 
                                               pqArm_file = pqArm_file, 
                                               difratio_chr = difratio_chr)
  
  # step 03: Subclone clustering
  Subclone_output <- SubClustering(input = input_df,
                                   Consolidating_output = Consolidation_result,
                                   min_cell = min_cell, 
                                   overlap_region = overlap_region, 
                                   dif_ratio = dif_ratio)
  
  # step 04: Output
  scDNA_Output(input = input_df,
               Summary = Subclone_output,
               output_dir = output_dir,
               pqArm_file = pqArm_file,
               consecutive_region = consecutive_region,
               cellcutoff = cellcutoff,
               sexchromosome = sexchromosome,
               smoothheatmap = smoothheatmap)
  
  endTimed(ptm)
}
