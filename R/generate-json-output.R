#' Get gene-level changes in copy number from FACETS output.
#'
#' This maps protein-coding genes onto copy-number segmentation from FACETS output derives the copy-number status from that, based on the integer copy number.
#' Additionally, the copy-number log-ratio of the SNPs falling insided the gene boundaries are used to perform a Z-test against the dipLogR baseline.
#' 
#' @param hisens_output Full FACETS hisensitivity output from \code{run_facets}.
#' @param purity_output Full FACETS purity output from \code{run_facets}.
#' @param gene_level    Full gene level output from \code{gene-level-changes}.
#' @param arm_level     Full arm level output from \code{arm-level-changes}.
#' @return \code{json}
#' 
#' @import data.table
#' @import jsonlite
#' 
#' @examples
#' \dontrun{
#' generate-json-ouput(hisens_output, NULL, gene_level_ouput, arm_level_output, parameters)
#' }

#' @export generate_json
generate_json = function(hisens_output,
                            purity_output,
                            gene_level,
                            arm_level,
                            parameters) {
    # TODO: include arm_level output currently not included in json?

    # genomic cumsum coordinates for plotting
    #genome = get(parameters$genome)
    genome  = tibble::tribble(
    ~chrom,     ~size, ~centstart,  ~centend, ~centromere,
    1, 249250621,  121535434, 124535434,   121535434,
    2, 243199373,   92326171,  95326171,    92326171,
    3, 198022430,   90504854,  93504854,    90504854,
    4, 191154276,   49660117,  52660117,    49660117,
    5, 180915260,   46405641,  49405641,    46405641,
    6, 171115067,   58830166,  61830166,    58830166,
    7, 159138663,   58054331,  61054331,    58054331,
    8, 146364022,   43838887,  46838887,    43838887,
    9, 141213431,   47367679,  50367679,    47367679,
    10, 135534747,   39254935,  42254935,    39254935,
    11, 135006516,   51644205,  54644205,    51644205,
    12, 133851895,   34856694,  37856694,    34856694,
    13, 115169878,    1.6e+07,   1.9e+07,     1.6e+07,
    14, 107349540,    1.6e+07,   1.9e+07,     1.6e+07,
    15, 102531392,    1.7e+07,     2e+07,     1.7e+07,
    16,  90354753,   35335801,  38335801,    35335801,
    17,  79759049,   22263006,  25263006,    22263006,
    18,  78077248,   15460898,  18460898,    15460898,
    19,  59128983,   24681782,  27681782,    24681782,
    20,  63025520,   26369569,  29369569,    26369569,
    21,  48129895,   11288129,  14288129,    11288129,
    22,  51304566,    1.3e+07,   1.6e+07,     1.3e+07,
    23, 155270560,   58632012,  61632012,    58632012,
    24,  59373566,   10104553,  13104553,    10104553
)
    VERTICAL_LINE_X_COORD_LIST = c(0, cumsum(genome$size))

    # x-axis tick values for chromosomes
    Tick_Values = floor(head(VERTICAL_LINE_X_COORD_LIST, -1) +
                          (tail(VERTICAL_LINE_X_COORD_LIST, -1) - 
                            head(VERTICAL_LINE_X_COORD_LIST, -1))/2)


    # seg data for line plots
    segs = as.data.table(hisens_output$segs)

    # gene data for line plots
    genes = as.data.table(gene_level)

    # snp level data for scatter plots
    snps = as.data.table(hisens_output$snps)

    # assign dummy end location, which will be
    #  same as start location. Both start and end 
    #  are required for foverlaps function
    snps$maploc_end = snps$maploc

    # assing snp to a gene base on genomics coordinates
    setkey(genes, chrom, gene_start, gene_end)
    setkey(snps, chrom, maploc, maploc_end)
    snps_gene_mapping = foverlaps(snps, genes, type="within")

    # seg data for line plots mapped to all the genes
    #  within a given segment
    hisens = as.data.table(hisens_output$segs)
    hisens[aggregate(genes[,c(2)], 
                    by = genes[,c(3,10,11)], function(genes) paste(genes, collapse = ",")), 
          on=c(chrom = "chrom", 
                start = "seg_start", 
                end = "seg_end"), gene_list := i.gene]


    # log ratio plot line data
    Line_Data_CNLR = list()
    Segment_mean = list()
    for (i in seq(nrow(hisens))) {
      Line_Data_CNLR[[i]] = list(
        "id" = unbox(paste0("Chromosome ", hisens[i,]$chrom, ": Copy #", i)),
        "data" = list(
          list(
            "x" = unbox(VERTICAL_LINE_X_COORD_LIST[hisens[i,]$chrom] + hisens[i,]$start), 
            "y" = unbox(hisens[i,]$cnlr.median)),
          list(
            "x" = unbox(VERTICAL_LINE_X_COORD_LIST[hisens[i,]$chrom] + hisens[i,]$end), 
            "y" = unbox(hisens[i,]$cnlr.median))))
      
      #seg_key = paste0("Chromosome ", as.character(hisens[i,]$chrom), " Copy #", i)
      Segment_mean[[i]] = unbox(hisens[i,]$cnlr.median)
      names(Segment_mean)[i] = paste0("Chromosome ", as.character(hisens[i,]$chrom), ": Copy #", i)
    }

    # log ratio and log odds ration plot scatter data
    snps = as.data.table(hisens_output[[1]]$snps)
    CNLR_Scatter_Data = VALOR_Scatter_Data = list()
    for (chrom_var in unique(snps$chrom)) {
      cnlr_data_points = valor_data_points = list()
      snps_per_chrom = snps[snps$chrom == chrom_var,]
      # additional counters are being defined here instead of
      #  using the i variable in the for loop below as a counter
      #  because certain snps might have NA as value for valor
      #  and they will be skipped. But the indices of the 
      #  valor_data_points list should be continunous
      cnlr_data_point_counter = valor_data_point_counter = 1
      for (i in seq(1, nrow(snps_per_chrom))) {
        cnlr_data_points[[cnlr_data_point_counter]] = list(
          "x" = unbox(VERTICAL_LINE_X_COORD_LIST[chrom_var] + snps_per_chrom[i,]$maploc), 
          "y" = unbox(snps_per_chrom[i,]$cnlr))
        if (!is.na(snps_per_chrom[i,]$valor)) {
          valor_data_points[[valor_data_point_counter]] = list(
            "x" = unbox(VERTICAL_LINE_X_COORD_LIST[chrom_var] + snps_per_chrom[i,]$maploc), 
            "y" = unbox(snps_per_chrom[i,]$valor))
          # add +1 to valor_data_point_counter only when a 
          #  snp with a floating point value for valor (i.e., hetSNP) is
          #  encountered.
          valor_data_point_counter = valor_data_point_counter + 1
        }
        # add +1 to cnlr_data_point_counter for every SNP
        cnlr_data_point_counter = cnlr_data_point_counter + 1
      }
      CNLR_Scatter_Data[[chrom_var]] = list(
        "id" = unbox(paste0("Chromosome ", chrom_var)),
        "data" = cnlr_data_points
      )
      VALOR_Scatter_Data[[chrom_var]] = list(
        "id" = unbox(paste0("Chromosome ", chrom_var)),
        "data" = valor_data_points
      )
    }


    # log odds ratio plot line data
    Line_Data_VALOR = Log_DR = list()
    Log_DR_counter = 1
    for (i in seq(nrow(hisens))) {
      Line_Data_VALOR[[Log_DR_counter]] = list(
          "id" = unbox(paste0("Chromosome ", hisens[i,]$chrom, ": Copy #", i, "_Min")),
          "data" = list(
            list(
              "x" = unbox(VERTICAL_LINE_X_COORD_LIST[hisens[i,]$chrom] + hisens[i,]$start), 
              "y" = unbox(-sqrt(abs(hisens[i,]$mafR)))),
            list(
              "x" = unbox(VERTICAL_LINE_X_COORD_LIST[hisens[i,]$chrom] + hisens[i,]$end), 
              "y" = unbox(-sqrt(abs(hisens[i,]$mafR))))))
      Line_Data_VALOR[[Log_DR_counter+1]] = list(
          "id" = unbox(paste0("Chromosome ", hisens[i,]$chrom, ": Copy #", i, "_Max")),
          "data" = list(
            list(
              "x" = unbox(VERTICAL_LINE_X_COORD_LIST[hisens[i,]$chrom] + hisens[i,]$start), 
              "y" = unbox(sqrt(abs(hisens[i,]$mafR)))),
            list(
              "x" = unbox(VERTICAL_LINE_X_COORD_LIST[hisens[i,]$chrom] + hisens[i,]$end), 
              "y" = unbox(sqrt(abs(hisens[i,]$mafR))))))

      Log_DR[[Log_DR_counter]] = unbox(-sqrt(abs(hisens[i,]$mafR)))
      Log_DR[[Log_DR_counter+1]] = unbox(sqrt(abs(hisens[i,]$mafR)))
      names(Log_DR)[Log_DR_counter] = paste0("Chromosome ", as.character(hisens[i,]$chrom), ": Copy #", i, "_Min")
      names(Log_DR)[Log_DR_counter+1] = paste0("Chromosome ", as.character(hisens[i,]$chrom), ": Copy #", i, "_Max")
      Log_DR_counter = Log_DR_counter + 2
    }

    # integer copy number line data
    CNCF_ICNP = EM_ICNP = list()
    ICNP_counter = 1
    for (i in seq(nrow(segs))) {
      CNCF_ICNP[[ICNP_counter]] = list(
        "id" = unbox(paste0("Chromosome ", segs[i,]$chrom, ": Copy #", i, "_TCN")),
        "data" = list(
          list(
            "x" = unbox(VERTICAL_LINE_X_COORD_LIST[segs[i,]$chrom] + segs[i,]$start), 
            "y" = unbox(segs[i,]$tcn)),
          list(
            "x" = unbox(VERTICAL_LINE_X_COORD_LIST[segs[i,]$chrom] + segs[i,]$end), 
            "y" = unbox(segs[i,]$tcn))))
      CNCF_ICNP[[ICNP_counter+1]] = list(
        "id" = unbox(paste0("Chromosome ", segs[i,]$chrom, ": Copy #", i, "_LCN")),
        "data" = list(
          list(
            "x" = unbox(VERTICAL_LINE_X_COORD_LIST[segs[i,]$chrom] + segs[i,]$start), 
            "y" = unbox(segs[i,]$lcn)),
          list(
            "x" = unbox(VERTICAL_LINE_X_COORD_LIST[segs[i,]$chrom] + segs[i,]$end), 
            "y" = unbox(segs[i,]$lcn))))
      EM_ICNP[[ICNP_counter]] = list(
        "id" = unbox(paste0("Chromosome ", segs[i,]$chrom, ": Copy #", i, "_TCN_EM")),
        "data" = list(
          list(
            "x" = unbox(VERTICAL_LINE_X_COORD_LIST[segs[i,]$chrom] + segs[i,]$start), 
            "y" = unbox(segs[i,]$tcn.em)),
          list(
            "x" = unbox(VERTICAL_LINE_X_COORD_LIST[segs[i,]$chrom] + segs[i,]$end), 
            "y" = unbox(segs[i,]$tcn.em))))
      EM_ICNP[[ICNP_counter+1]] = list(
        "id" = unbox(paste0("Chromosome ", segs[i,]$chrom, ": Copy #", i, "_LCN_EM")),
        "data" = list(
          list(
            "x" = unbox(VERTICAL_LINE_X_COORD_LIST[segs[i,]$chrom] + segs[i,]$start), 
            "y" = unbox(segs[i,]$lcn.em)),
          list(
            "x" = unbox(VERTICAL_LINE_X_COORD_LIST[segs[i,]$chrom] + segs[i,]$end), 
            "y" = unbox(segs[i,]$lcn.em))))
      ICNP_counter = ICNP_counter + 2
    }
      
    # Gene level data segment of the json
    gene_dict = list()
    gene_dict_counter = 1
    for (i in seq(nrow(hisens))) {
      gene_list = strsplit(hisens[i,]$gene_list, ",")[[1]]
      # a segment has no gene, skip and go to the next
      if(all(is.na(gene_list))) {
        next
      }
      gene_data = list()
      gene_counter = 1
      for (g in gene_list) {
        scatter_data_for_gene = snps_gene_mapping[snps_gene_mapping$gene == g,]
        
        # snp level data organization for each gene
        # scatter_data_for_gene should never be empty for a given segment
        scatter_data_counter = 1
        gene_scatter_data_cnlr = list()
        gene_scatter_data_valor = list()
        for (s in seq(1, nrow(scatter_data_for_gene))) {
          gene_scatter_data_cnlr[[scatter_data_counter]] = list(
            "x" = unbox(VERTICAL_LINE_X_COORD_LIST[hisens[i,]$chrom] + scatter_data_for_gene[s,]$maploc), 
            "y" = unbox(scatter_data_for_gene[s,]$cnlr))
          if (!is.na(scatter_data_for_gene[s,]$valor)) {
            gene_scatter_data_valor[[scatter_data_counter]] = list(
              "x" = unbox(VERTICAL_LINE_X_COORD_LIST[hisens[i,]$chrom] + scatter_data_for_gene[s,]$maploc), 
              "y" = unbox(scatter_data_for_gene[s,]$valor))
          }
          scatter_data_counter = scatter_data_counter + 1
        }
        
        # gene level data organization
        # gene_cnlr_counter and gene_valor_counter should, in theory, will always be 1
        #  in the current implementation. If facets were to support intragenic gene level
        #  calls in the future with multiple rows for a given gene in the gene level file,
        #  this block needs to be updated to indicate multiple segments within a gene. These
        #  counters will be server to distinguish the segments within a gene.
        gene_cnlr_counter = 1
        gene_valor_counter = 1
        gene_start = unbox(VERTICAL_LINE_X_COORD_LIST[hisens[i,]$chrom] + genes[gene == g,]$gene_start)
        gene_end = unbox(VERTICAL_LINE_X_COORD_LIST[hisens[i,]$chrom] + genes[gene == g,]$gene_end)
        gene_tcn = unbox(hisens[i,]$tcn)
        gene_lcn = unbox(hisens[i,]$lcn)
        gene_tcn.em = unbox(hisens[i,]$tcn.em)
        gene_lcn.em = unbox(hisens[i,]$lcn.em)
        gene_cnlr_median = unbox(hisens[i,]$cnlr.median)
        gene_valor_positive_y = unbox(sqrt(abs(hisens[i,]$mafR)))
        gene_valor_negative_y = unbox(-sqrt(abs(hisens[i,]$mafR)))
        gene_line_data_cnlr = list(
          "id" = unbox(g),
          "data" = list(
            list(
              "x" = gene_start, 
              "y" = gene_cnlr_median),
            list(
              "x" = gene_end, 
              "y" = gene_cnlr_median)))
        gene_line_data_valor_min = list(
          "id" = unbox(paste0(g, " Copy #", as.character(gene_valor_counter), "Min")),
          "data" = list(
            list(
              "x" = gene_start, 
              "y" = gene_valor_negative_y),
            list(
              "x" = gene_end, 
              "y" = gene_valor_negative_y)))
        gene_line_data_valor_max = list(
          "id" = unbox(paste0(g, " Copy #", as.character(gene_valor_counter), "_Max")),
          "data" = list(
            list(
              "x" = gene_start, 
              "y" = gene_valor_positive_y),
            list(
              "x" = gene_end, 
              "y" = gene_valor_positive_y)))
        gene_line_data_CNCF_TCN = list(
          "id" = unbox(paste0(g, " Copy #", as.character(gene_cnlr_counter), "_TCN")),
          "data" = list(
            list(
              "x" = gene_start, 
              "y" = gene_tcn,
              list(
                "x" = gene_end, 
                "y" = gene_tcn))))
        gene_line_data_CNCF_LCN = list(
          "id" = unbox(paste0(g, " Copy #", as.character(gene_cnlr_counter), "_LCN")),
          "data" = list(
            list(
              "x" = gene_start, 
              "y" = gene_lcn,
              list(
                "x" = gene_end, 
                "y" = gene_lcn))))
        gene_line_data_EM_TCN = list(
          "id" = unbox(paste0(g, " Copy #", as.character(gene_cnlr_counter), "_TCN_EM")),
          "data" = list(
            list(
              "x" = gene_start, 
              "y" = gene_tcn.em,
              list(
                "x" = gene_end, 
                "y" = gene_tcn.em))))
        gene_line_data_EM_LCN = list(
          "id" = unbox(paste0(g, " Copy #", as.character(gene_cnlr_counter), "_LCN_EM")),
          "data" = list(
            list(
              "x" = gene_start, 
              "y" = gene_lcn.em,
              list(
                "x" = gene_end, 
                "y" = gene_lcn.em))))
        
        # organize line and scatter data for the gene, g
        gene_data_current_gene = list(
          g = list(
            "CNLRPlot" = list(
              "Line_Data" = gene_line_data_cnlr,
              "Scatter_Data" = list(
                "id" = unbox(g),
                "data" = gene_scatter_data_cnlr)
            ),
            "VALORPlot" = list(
              "Line_Data" = list(
                gene_line_data_valor_min, 
                gene_line_data_valor_max),
              "Scatter_Data" = list(
                "id" = unbox(g),
                "data" = gene_scatter_data_valor)
            ),
            "ICNPlot" = list(
              "Line_Data" = list(
                "CNCF" = list(
                  gene_line_data_CNCF_TCN,
                  gene_line_data_CNCF_LCN),
                "EM" = list(
                  gene_line_data_EM_TCN,
                  gene_line_data_EM_LCN)
              ),
              "Scatter_Data" = NULL
            ),
            "Gene_Info" = unbox(scatter_data_for_gene$cn_state[1])
          )
        )
        gene_data[[gene_counter]] = gene_data_current_gene
        names(gene_data)[gene_counter] = unbox(g)
        gene_counter = gene_counter + 1
      }
      # organize data for all genes in a given segment
      segment_id = paste0("Chromosome ", as.character(hisens[i,]$chrom), ": Copy #", gene_dict_counter)
      gene_dict[[gene_dict_counter]] = list(
        "Gene_List" = gene_list,
        "Gene_Data" = gene_data)
      names(gene_dict)[gene_dict_counter] = unbox(segment_id)
      gene_dict_counter = gene_dict_counter + 1
    }

    jsonlite::write_json(
      list(
        "FACETS_PROCESSED_DATA_FOR_PLOT" = list(
          "VERTICAL_LINE_X_COORD_LIST" = VERTICAL_LINE_X_COORD_LIST, 
          "Tick_Values" = Tick_Values, 
          "CNLRPlot" = list("Line_Data" =  Line_Data_CNLR,
                            "Segment_mean" = Segment_mean,
                            "Scatter_Data" = CNLR_Scatter_Data),
          "VALORPlot" = list("Line_Data" = Line_Data_VALOR,
          # TODO: is there a better way to label this data than using "Log_DR"?
                            "Log_DR" = Log_DR,
                            "Scatter_Data" = VALOR_Scatter_Data),
          "ICNPlot" = list("CNCF" = CNCF_ICNP,
                          "EM" = EM_ICNP),
          "Gene_Dictionary" = gene_dict,
          "All_Genes_List" = genes$gene),

        # TODO: add arm level data to the json?
        # "FACETS_ARM_LEVEL_DATA" = list(),

        # TODO: include other parameters printed to run details file?
        "PARAMETERS" = list(
          # input parameters to facets
          "sample_id" = unbox(parameters$sample_id),
          "purity_cval" = unbox(parameters$purity_cval),
          "cval" = unbox(parameters$cval),
          "dipLogR" = unbox(hisens_output$dipLogR),
          "purity_min_het" = unbox(parameters$purity_min_nhet),
          "min_het" = unbox(parameters$min_nhet),
          "normal_depth" = unbox(parameters$ndepth),
          "snp_window_size" = unbox(parameters$snp_nbhd),

          # variables from facets fit
          "loglik" = unbox(purity_output$loglik),
          "alBalLogR" = as.vector(purity_output$alBalLogR),
          # TODO: print purity and ploidy from purity fit?
          "purity" = unbox(hisens_output$purity),
          "ploidy" = unbox(hisens_output$ploidy),

          # facets2N parameters
          "MandUnormal" = unbox(parameters$MandUnormal),
          "useMatchedX" = unbox(parameters$useMatchedX),
          "cbs" = unbox(parameters$cbs),
          "het_thresh" = unbox(parameters$het_thresh),
          # boolean variable to indicate if donor snp count data was used
          "donor_counts" = unbox(parameters$donor_counts),
          
          # general parameters/variables
          "genome" = unbox(parameters$genome),
          "facets_version" = unbox(parameters$facets_version),
          "facetsSuite_version" = unbox(parameters$facetsSuite_version),
          "unmatched" = unbox(parameters$unmatched),
          "seed" = unbox(parameters$seed)
        )
      ), 
      pretty=TRUE, 
      path = paste0(parameters$outdir, "/", parameters$sample_id, "_facets_output.json")
    )
}



