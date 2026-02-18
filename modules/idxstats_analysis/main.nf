process ANALYZE_IDXSTATS {
    tag "Analyzing idxstats and loci distribution"
    label 'process_small'
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project}/idxstats_analysis", mode: 'copy'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-70a811298b15e02937b3477ebadce5062ab47a53:8e6675e8830568a8642779e4ba23787ea188a067-0':
        'quay.io/biocontainers/mulled-v2-70a811298b15e02937b3477ebadce5062ab47a53:8e6675e8830568a8642779e4ba23787ea188a067-0' }"

    input:
    path idxstats_files

    output:
    path "*.{pdf,jpg,png}", emit: plots
    path "reads_matrix.txt", emit: matrix
    path "*.{csv,txt}", emit: loci_qc
    path "*sex_id_results.tsv", emit: sexid
    
    script:
    """
    #!/usr/bin/env Rscript

    # Load required libraries
    library(tidyverse)
    library(viridis)
    library(vegan)
    library(cowplot)

    # Constants
    TOTAL_PANEL_LOCI <- 204
    MIN_READS <- 10
    POOR_PERFORMANCE_THRESHOLD <- 0.5
    sex_id_male_threshold <- $params.male_sexid_threshold # must have more than this many sdy_I183 reads mapped
    sex_id_min_total_reads <- $params.sexid_min_reads # must have at least this many total reads mapped or else sample is unknown sex
    sex_id_female_threshold <- $params.female_sexid_threshold # must have fewer than this many sdy_I183 reads mapped

    ## Functions ---------------
    # Function to process a idxstats files
    process_idxstats <- function(file, ind, n_loci) {
      df <- read.table(file, nrows = n_loci, stringsAsFactors = FALSE)
      df <- df[df\$V1 != "*", ] # remove rows with "*" in the first column
      df\$ind <- ind
      return(df)
    }

    # Sex ID function
    determine_sex <- function(df, sexid_locus = "sdy_I183", min_total_reads, min_male_prop, max_female_prop) {
      df <- df %>%
        filter(str_detect(ind, "_Chinook_FullPanel_VGLL3Six6LFARWRAP", negate = TRUE)) # just in case any fgm data is missed
      total_read_counts <- df %>%
        group_by(ind) %>%
        summarize(total_reads = sum(reads))
      sex_marker_reads <- df %>%
        filter(loc == sexid_locus) %>%
        select(ind, sex_marker_reads = reads)
      sex_id_df <- total_read_counts %>%
        left_join(sex_marker_reads, by = 'ind') %>%
        mutate(
          sex_marker_read_prop = sex_marker_reads / total_reads,
          inferred_sex = case_when(
            sex_marker_read_prop > min_male_prop & total_reads > min_total_reads ~ "Male",
            sex_marker_read_prop < max_female_prop & total_reads > min_total_reads ~ "Female",
            TRUE ~ "Unknown"
            )
          )
      return(sex_id_df)
    }
    
    # List all idxstats files in the input directory
    target_mapped_files <- tibble(
        file = list.files(pattern = "_full_idxstats.txt\$", full.names = TRUE)
        ) %>%
      mutate(
        ind = str_remove(basename(file), "_((full|transition)_)?idxstats.txt")
      )

    fgm_files <- tibble(
        file = list.files(pattern = "_Chinook_FullPanel_VGLL3Six6LFARWRAP_idxstats.txt\$", full.names = TRUE)
        ) %>%
      mutate(
        ind = str_remove(basename(file), "_Chinook_FullPanel_VGLL3Six6LFARWRAP_idxstats.txt")
      )

    file_list <- bind_rows(target_mapped_files, fgm_files)
    
    # Process all files
    idx_list <- map2(file_list\$file, file_list\$ind, ~ process_idxstats(.x, .y, TOTAL_PANEL_LOCI))
    
    # Combine into a single data frame
    idx_df <- bind_rows(idx_list)
    names(idx_df) <- c("loc", "len", "reads", "unmapd", "ind")
    
    # Summarize and order data
    locs <- idx_df %>%
      group_by(loc) %>%
      summarise(total = sum(reads)) %>%
      arrange(total)
    
    inds <- idx_df %>%
      group_by(ind) %>%
      summarise(total = sum(reads)) %>%
      arrange(total)
    
    idx_df <- idx_df %>%
      mutate(
        loc = factor(loc, levels = locs\$loc),
        ind = factor(ind, levels = inds\$ind)
      )
    # Generate Sex ID info
    sex_id_info <- determine_sex(
      idx_df,
      min_total_reads = sex_id_min_total_reads,
      min_male_prop = sex_id_male_threshold,
      max_female_prop = sex_id_female_threshold
      )
    write_tsv(sex_id_info, "${params.project}_sex_id_results.tsv")

    # Generate idxstats plots
    plot_reads_by_ind <- ggplot(idx_df, aes(x = ind, y = reads)) +
      geom_col(fill = "lightblue", colour = "black") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8)) +
      xlab("Sample") +
      ylab("Total reads mapped") +
      ggtitle("Reads per individual")
    
    plot_reads_by_locus <- ggplot(idx_df, aes(x = loc, y = reads)) +
      geom_col(fill = "goldenrod", colour = "black") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
      xlab("Marker Name, Amplicon Name, or Chromosome") +
      ylab("Total reads mapped") +
      ggtitle("Reads per locus")
    
    # Save idxstats plots
    ggsave("reads_per_individual.pdf", plot_reads_by_ind, width = 12, height = 8)
    ggsave("reads_per_locus.pdf", plot_reads_by_locus, width = 12, height = 8)
    
    # Check for duplicate sample IDs before exporting
    duplicate_samples <- idx_df %>% 
      select(ind, loc) %>%
      group_by(ind, loc) %>%
      summarise(n = n(), .groups = "drop") %>%
      filter(n > 1) %>%
      arrange(ind, loc)
    
    # If duplicates found, write them to a file and throw an error
    if (nrow(duplicate_samples) > 0) {
      unique_duplicated_inds <- duplicate_samples %>%
        select(ind) %>%
        distinct() %>%
        pull(ind)
      
      # Write duplicates to a file for reference
      write.csv(duplicate_samples, "duplicate_samples.csv", row.names = FALSE)
      
      # Generate error message with duplicated sample IDs
      error_message <- paste0(
        "ERROR: Duplicate sample IDs detected for the following samples: ",
        paste(unique_duplicated_inds, collapse = ", "),
        ". See duplicate_samples.csv for details. This likely happened because multiple input files were reduced to the same sample ID during name extraction. Please use unique sample IDs or modify the sample name extraction function."
      )
      
      # Write error to file for Nextflow to capture
      sink("error_message.txt")
      cat(error_message)
      sink()
      
      # Stop execution with error code
      stop(error_message)
    }
    
    # Export matrix (only runs if no duplicates are found)
    matrix_df <- idx_df %>%
      select(-len, -unmapd) %>%
      spread(ind, reads)
    write.table(matrix_df,
      file = "reads_matrix.txt",
      quote = FALSE, sep = "\t", row.names = FALSE
    )
    
    # Read and process data for loci check
    data <- matrix_df
    
    # Convert to long format
    long_data <- data %>%
      pivot_longer(
        cols = -loc,
        names_to = "sample",
        values_to = "reads"
      ) %>%
      filter(loc != "*")
    
    # Calculate sample-level metrics
    sample_stats <- long_data %>%
      group_by(sample) %>%
      summarize(
        total_reads = sum(reads),
        mean_reads = mean(reads),
        median_reads = median(reads),
        sd_reads = sd(reads),
        cv = sd_reads / mean_reads * 100,
        non_zero_loci = sum(reads > MIN_READS),
        total_loci = TOTAL_PANEL_LOCI,
        pct_success = (non_zero_loci / TOTAL_PANEL_LOCI) * 100,
        shannon = diversity(reads, index = "shannon"),
        pielou = shannon / log(TOTAL_PANEL_LOCI),
        simpson = diversity(reads, index = "simpson")
      ) %>%
      arrange(desc(total_reads))
    
    # Calculate locus-level metrics
    locus_stats <- long_data %>%
      group_by(loc) %>%
      summarize(
        total_reads = sum(reads),
        mean_reads = mean(reads),
        median_reads = median(reads),
        sd_reads = sd(reads),
        cv = sd_reads / mean_reads * 100,
        non_zero_samples = sum(reads > MIN_READS),
        total_samples = n(),
        pct_success = (non_zero_samples / total_samples) * 100
      ) %>%
      arrange(desc(pct_success))
    
    # Identify problematic loci
    problem_loci <- locus_stats %>%
      filter(pct_success < (POOR_PERFORMANCE_THRESHOLD * 100)) %>%
      arrange(pct_success)
    
    # Create visualizations
    plot_total_reads_by_sample <- ggplot(sample_stats, aes(x = reorder(sample, -total_reads), y = total_reads)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        title = "Total Reads per Sample",
        x = "Sample",
        y = "Total Reads"
      )
    
    plot_success_vs_reads <- ggplot(sample_stats, aes(x = total_reads, y = pct_success)) +
      geom_point(aes(size = non_zero_loci), alpha = 0.6) +
      geom_smooth(method = "lm", alpha = 0.2) +
      theme_minimal() +
      labs(
        title = "Success Rate vs Total Reads",
        x = "Total Reads",
        y = "Percent Success",
        size = "Successful Loci"
      )
    
    # Create enhanced read depth heatmap (combining the best features of both prior heatmaps)
    heatmap_data <- long_data %>%
      mutate(read_category = cut(reads,
        breaks = c(-Inf, 0, 10, 100, 1000, Inf),
        labels = c("0", "1-10", "11-100", "101-1000", ">1000")
      ))
    
    plot_read_depth_heatmap <- ggplot(
      heatmap_data,
      aes(
        x = reorder(sample, -reads),
        y = reorder(loc, reads),
        fill = read_category
      )
    ) +
      geom_tile() +
      scale_fill_viridis_d(option = "viridis") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6)
      ) +
      labs(
        title = "Read Depth Distribution",
        x = "Sample",
        y = "Locus",
        fill = "Read Count"
      )
    
    plot_amplification_evenness <- ggplot(sample_stats, aes(x = reorder(sample, -pielou), y = pielou)) +
      geom_bar(stat = "identity", fill = "darkgreen") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        title = "Amplification Evenness by Sample",
        x = "Sample",
        y = "Pielou's Evenness"
      )
    
    # Save results
    write.csv(sample_stats,
      "sample_statistics.csv",
      row.names = FALSE
    )
    write.csv(locus_stats,
      "locus_statistics.csv",
      row.names = FALSE
    )
    write.csv(problem_loci,
      "problematic_loci.csv",
      row.names = FALSE
    )
    
    # Save plots
    combined_qc_plots <- plot_grid(
      plot_total_reads_by_sample, 
      plot_success_vs_reads, 
      plot_read_depth_heatmap, 
      plot_amplification_evenness, 
      ncol = 2
    )
    
    ggsave("loci_qc_plots.pdf",
      combined_qc_plots,
      width = 15,
      height = 12
    )
    
    # Generate summary report
    sink("loci_qc_summary.txt")
    cat("Amplicon Sequencing QC Analysis Summary\\n")
    cat("======================================\\n\\n")
    cat(sprintf("Total Samples: %d\\n", nrow(sample_stats)))
    cat(sprintf("Total Loci: %d\\n", TOTAL_PANEL_LOCI))
    cat(sprintf("Mean Reads per Sample: %.0f\\n", mean(sample_stats\$total_reads)))
    cat(sprintf("Median Reads per Sample: %.0f\\n", median(sample_stats\$total_reads)))
    cat(sprintf("Mean Success Rate: %.1f%%\\n", mean(sample_stats\$pct_success)))
    cat(sprintf("Number of Problematic Loci: %d\\n", nrow(problem_loci)))
    cat(sprintf("Mean Amplification Evenness: %.3f\\n", mean(sample_stats\$pielou)))
    sink()
    """
}
