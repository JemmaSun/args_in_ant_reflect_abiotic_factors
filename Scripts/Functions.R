#  _____                 _   _                 
# |  ___|   _ _ __   ___| |_(_) ___  _ __  ___ 
# | |_ | | | | '_ \ / __| __| |/ _ \| '_ \/ __|
# |  _|| |_| | | | | (__| |_| | (_) | | | \__ \
# |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
#
#  Jiarui Sun
colors37 = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36",
             "#552095","#507f2d","#db37aa","#84b67c","#a06fda","#df462a","#5b83db","#c76c2d",
             "#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3",
             "#8c3025","#417d61","#862977","#bba672","#403367","#da8a6d","#a79cd4","#71482c",
             "#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")


summarise_table <- function(df, group_for_row, group_for_col, summarise_by_column = NULL, 
                            method = c('n','sum', 'max', 'min', 'mean')[1],
                            env.df = NULL, changeSampleName = FALSE, refColumn = NULL) {
  if (method == 'n') {
    Func <- 'n()'
  } else { # Note: summarise_by_column CANNOT be NULL!!!
    Func <- paste(method,'(',summarise_by_column,')', sep='')
  }
  tmp.df <- df %>% group_by(eval(parse(text=group_for_row)),eval(parse(text=group_for_col))) %>% 
    summarise(count=eval(parse(text=Func)))
  tmp.df <- as.data.frame(tmp.df)
  colnames(tmp.df) <- c('Row_name', 'Col_name', 'count') 
  new.df <- tmp.df %>% ungroup() %>%
    mutate(Row_name=Row_name,Col_name=Col_name) %>% 
    spread(Col_name,count,fill=0)
  new.df <- as.data.frame(new.df)
  if (changeSampleName) {
    if (length(refColumn) == 1) { # If just column name is given, e.g. 'sample_id.GSVeasy' instead of env.df$sample_id.GSVeasy
      refColumn <- env.df[[refColumn]]
    }
    rownames(new.df) <- rownames(env.df)[match(new.df[,1], refColumn)]
    if (nrow(new.df) < nrow(env.df)) {
      missing_sample <- nrow(env.df) - nrow(new.df)
      message(paste("Warning: ", missing_sample, "samples have 0 hits. Consider removing them in downsteam analyses. "))
      
    }
    new.df <- new.df[rownames(env.df),]
    new.df[is.na(new.df)] <- 0
    rownames(new.df) <- rownames(env.df)
  } else {
    rownames(new.df) <- new.df[,1]
  }
  new.df <- new.df[,-1]
  return(new.df)
}

draw_cmm_heatmap <- function(otu, taxonomy_df, env_df = NULL, abundance_cutoff = NULL) {
  # First check OTU table is in correct format and contains small number of OTUs:
  if (!is.null(abundance_cutoff)) {
    otu <- otu[,which(apply(otu, 2, max) >= abundance_cutoff)]
  }
  if (all(dim(otu)>0)) {
    if (ncol(otu) > 100) {
      stop("You have more than 100 OTUs, please OTUs in a reasonable number (preferably ~ 50)")
    }
  } else {
    stop("Data frame dimensions too small!")
  }
  # Sort out taxonomy metadata:
  tax.abd <- taxonomy_df[colnames(otu),]
  tax.abd <- tax.abd[order(tax.abd$Taxonomy),]
  # Remove 'Bacteria;unclassified' or 'Archaea;unclassified' --
  #tax.abd <- tax.abd[-which(is.na(tax.abd$Phylum)),]
  
  tax.abd$p_lab <- tax.abd$Phylum
  tax.abd$c_lab <- paste("c__", tax.abd$Class, sep='')
  tax.abd$o_lab <- paste("o__", tax.abd$Order, sep='')
  for (x in c("p_lab", 'c_lab', 'o_lab')) {tax.abd[[x]][duplicated(tax.abd[[x]])] <- ""}
  otu <- otu[,tax.abd$OTU]
  
  # Any annotations to be added to the heatmap:
  # a) Categorical -- also show gaps in heatmap!
  # To-do!
  top_ann <- NULL
  bot_ann <- NULL
  col_split <- NULL
  
  # b) Numerical --
  # To-do!
  
  # Build rowAnnotation for taxonomy labels:
  l_ann <- rowAnnotation(
    phylum = anno_text(tax.abd$p_lab, gp = gpar(fontface = 'bold', fontsize = 8)),
    class = anno_text(tax.abd$c_lab, gp = gpar(fontsize = 8)),
    order = anno_text(tax.abd$o_lab, gp = gpar(fontsize = 8)),
    show_legend = F,
    show_annotation_name = F
  )
  r_ann <- rowAnnotation(
    fgs = anno_text(tax.abd$hm_style, gp = gpar(fontsize = 8))
  )
  
  # Draw heatmap:
  ht <- Heatmap(mat = as.matrix(t(otu)), col = colorRampPalette(c("white","grey15"))(100),
                border = "dark grey", rect_gp = gpar(col="grey",lwd=.5),
                cluster_rows = F, cluster_columns = F,
                column_split = col_split, row_split = tax.abd$Phylum,
                show_column_names = T, column_names_gp = gpar(fontsize = 8),
                show_row_names = FALSE,
                top_annotation = top_ann, bottom_annotation = bot_ann,
                left_annotation = l_ann, right_annotation = r_ann,
                heatmap_legend_param = list(title = "Relative abundance (%)"),
                show_heatmap_legend = T, row_title = NULL,
                width = nrow(otu)*unit(4, "mm"),
                height = ncol(otu)*unit(3.8, "mm"),
  )
  ann_legends <- list() # To be complete for metadata annotation
  draw(ht, annotation_legend_list = ann_legends)
  message("Done. A ComplexHeatmap object and Legend object have been added to the returned list.\n")
  cat("To modify or re-draw this heatmap, use `$hm` to access the heatmap object, and `$hm_legends` to access the annotation legend list.\n")
  cat("Example code to re-draw the heatmap: \n")
  cat("\n    `draw(res$hm, annotation_legend_list = res$hm_legends)`")
  cat("     # Remember to replace object name 'res'.")
  return(list(hm=ht,hm_legends=ann_legends))
}


auto_lm_forward_new <- function(data, env_params, env_df, test_method = "Chisq",
                                hm_title="Linear model - t values",
                                corr.test = FALSE) {
  otu_df <- data
  if (!test_method %in% c("Chisq", "F")) {
    stop("Error:::Invalid test_method! Please use 'Chisq' or 'F'.")
  }
  # Step 0 - check if samples in otu_df are included in env_df:
  check_sample_matches <- sapply(rownames(otu_df), function(x) x%in%rownames(env_df))
  if (!all(check_sample_matches)) {
    stop("Missing sample(s) in metadata: ", paste(names(which(!check_sample_matches)),  collapse=', '))
  }
  check_params_matches <- sapply(env_params, function(x) x%in%colnames(env_df))
  if (!all(check_params_matches)) {
    stop("Missing env parameter(s) in metadata: ", paste(names(which(!check_params_matches)),  collapse=', '))
  }
  env_df <- as.data.frame(env_df[rownames(otu_df), env_params])
  rownames(env_df) <- rownames(otu_df)
  colnames(env_df) <- env_params
  # Step 1 - Initiate output objects:
  lm_results <- data.frame(Data = NA, Formula = NA, Parameter = NA, Df = NA,
                           F_value = NA, p_value = NA, Signif = NA) # Initiate output dataframe
  lm_matrix <- data.frame(matrix(0, nrow = length(env_params),
                                 ncol = length(colnames(otu_df))))
  rownames(lm_matrix) <- c(env_params)
  colnames(lm_matrix) <- colnames(otu_df)
  lm_matrix.p <- lm_matrix # Added here to plot corrplot
  if (test_method == "Chisq") {
    intermediate_results <- data.frame(Data = NA, Round = 0, Predictors = NA,
                                       Type = NA, Parameter = NA, Df = NA,
                                       AIC = NA, p = NA, Winner = NA)
  } else {
    intermediate_results <- data.frame(Data = NA, Round = 0, Predictors = NA,
                                       Type = NA, Parameter = NA, Df = NA,
                                       F_value = NA, p = NA, Winner = NA)
  }
  # Step 2a - If uni-param: run correlation test, draw a heat map and quit --
  if (corr.test) {
    message('Calculating correlations between data variables and metadata parameters ...')
    res_list <- calculate_uni_param_lm(otu_df, env_params, env_df)
    message('Done. Draw heat map ...')
    corr.lm_matrix <- res_list$Corr.Summary
    if (!all(is.na(corr.lm_matrix))) {
      my_color <- three_scale_color(min(corr.lm_matrix), max(corr.lm_matrix),
                                    colorRampPalette(c('blue','white','red'))(101))
      ht <- Heatmap(as.matrix(corr.lm_matrix), col = my_color, name = hm_title,
                    column_title = paste("(Uni-parameter Correlations) ", hm_title, sep = ""),
                    border = "black", rect_gp = gpar(col="black",lwd=.5),
                    heatmap_legend_param = list(title = "t values"),
                    column_names_gp = gpar(fontsize = 8), row_names_gp = gpar(fontsize = 8),
                    cluster_rows = F, cluster_columns = F,
                    width = ncol(corr.lm_matrix)*unit(4, "mm"),
                    height = nrow(corr.lm_matrix)*unit(3.8, "mm"))
      draw(ht)
    }
    # Return results and quit:
    #return(res_list)
  } else {
    res_list <- list()
  }
  # Step 2b - Loop across OTUs or Genes:
  for (var_name in colnames(otu_df)) {
    message(paste("Working on variable", var_name, "..."))
    env_df.tmp <- env_df[,env_params]
    Data <- otu_df[[var_name]]
    if (all(is.na(Data))) {
      message("Skipped as only NAs are found.")
      next
    }
    if (any(is.na(Data))) {
      env_df.tmp <- env_df.tmp[-which(is.na(Data)),]
      Data <- Data[-which(is.na(Data))]
    }
    if (sum(Data) == 0) {
      message("Skipped as gene count is 0.")
      next
    }
    # Step 2.1 - Initiate objects for while-loop:
    predictors <- list()
    mbig <- lm(Data ~ ., data = env_df.tmp)
    m0 <- lm(Data ~ 1, data = env_df.tmp)
    # Step 2.2 - Start while-loop:
    curr_round = 1
    while (curr_round > 0) {
      # Step 2.2.1 Add1:
      add.tmp <- data.frame(add1(m0, scope = formula(mbig), test = test_method))
      ##### Sort out results for output:
      if (test_method == "Chisq") {
        add.tmp <- add.tmp[-(which(is.na(add.tmp$Pr..Chi.))),] # Remove '<none>'
        add.tmp <- add.tmp[order(add.tmp$AIC),][,c(1,4,5)] # Df, AIC, Pr
        colnames(add.tmp) <- c("Df", "AIC", "p")
      } else {
        add.tmp <- add.tmp[-(which(is.na(add.tmp$F.value))),] # Remove '<none>'
        add.tmp <- add.tmp[order(-abs(add.tmp$F.value)),][,c(1,5,6)] # Df, F, Pr
        colnames(add.tmp) <- c("Df", "F_value", "p")
      }
      add.tmp$Winner = F
      add.tmp <- cbind(data.frame(Data = var_name, Round = curr_round,
                                  Predictors = paste(". ~ ",
                                                     paste(predictors,
                                                           collapse = "+")),
                                  Type = "add1", Parameter = rownames(add.tmp)),
                       add.tmp)
      # Step 2.2.2 Find winner:
      if (add.tmp$p[1] < 0.05){
        add.tmp$Winner[1] = T
        added_param <- add.tmp$Parameter[1]
      } else {
        curr_round = 0
        message("Nothing can be added in this round.")
        next # Nothing can be added.
      }
      intermediate_results <- rbind(intermediate_results, add.tmp)
      # Step 2.2.3 Update model with Winner:
      formula.temp <- as.formula(paste(". ~ .", "+", added_param))
      m0 <- update(m0, formula.temp)
      # Step 2.2.4 Drop1 to test the Winner:
      drop.tmp <- data.frame(drop1(m0, test = test_method))
      ##### Store results for output:
      if (test_method == "Chisq") {
        drop.tmp <- drop.tmp[-(which(is.na(drop.tmp$Pr..Chi.))),] # Remove '<none>'
        drop.tmp <- drop.tmp[order(drop.tmp$AIC),][,c(1,4,5)] # Df, AIC, Pr
        colnames(drop.tmp) <- c("Df", "AIC", "p")
      } else {
        drop.tmp <- drop.tmp[-(which(is.na(drop.tmp$F.value))),] # Remove '<none>'
        drop.tmp <- drop.tmp[order(-abs(drop.tmp$F.value)),][,c(1,5,6)] # Df, F, Pr
        colnames(drop.tmp) <- c("Df", "F_value", "p")
      }
      drop.tmp$Winner = NA
      drop.tmp <- cbind(data.frame(Data = var_name, Round = curr_round,
                                   Predictors = paste(". ~ ",
                                                      paste(predictors, collapse = "+")),
                                   Type = "drop1",
                                   Parameter = rownames(drop.tmp)), drop.tmp)
      intermediate_results <- rbind(intermediate_results, drop.tmp)
      # Step 2.2.5 Test if all predictors are remain significant:
      temp.mod.pvals <- drop.tmp$p
      if ((length(predictors)+1) != length(temp.mod.pvals)) {
        curr_round = 0
        stop("Error! drop1() results missing for predictors.")
      } # There maybe sth. wrong with drop1 results.
      if (length(temp.mod.pvals[which(temp.mod.pvals >= 0.05)])>0){
        curr_round = 0
        message(paste(added_param, "is not significant in parsimonious model so not included."))
        formula.temp <- as.formula(paste(". ~ .", "-", added_param))
        m0 <- update(m0, formula.temp)
      } else {
        curr_round = curr_round + 1
        predictors <- append(predictors, added_param)
        message(paste(added_param, "is significant in parsimonious model and thereby included."))
      }
      if ((length(env_params) - length(predictors)) < 1){
        curr_round = 0
        message("No more params to add! Quit while loop!")
      } # No more params to add!
    }# End while loop
    final_model <- paste(". ~ ", paste(predictors, collapse = "+"))
    m.final <- m0
    message(paste("Final model for", var_name, "is: ", final_model))
    print(summary(m.final))
    aov.tmp <- anova(m.final)
    # Store ANOVA format results to data frame:
    aov.tmp.df <- as.data.frame(aov.tmp[,c(1,4,5)])
    colnames(aov.tmp.df) <- c("Df", "F_value",	"p_value")
    aov.tmp.df <- aov.tmp.df %>%
      mutate(Signif = case_when(p_value < 0.001 ~ "***",
                                p_value < 0.01 ~  "**",
                                p_value < 0.05 ~ "*",
                                TRUE ~ " "))
    aov.tmp.df <- cbind(data.frame(Data = var_name,
                                   Formula = paste(". ~ ", paste(predictors,
                                                                 collapse = "+")),
                                   Parameter = rownames(aov.tmp.df)),
                        aov.tmp.df)
    lm_results <- rbind(lm_results, aov.tmp.df)
    # Store Summary format results to data matrix:
    res.tmp <- summary(m.final)
    res.tmp.df <- as.data.frame(res.tmp$coefficients)
    colnames(res.tmp.df) <- c("Estimate", "SE", "t_value", "p_value")
    res.tmp.df <- res.tmp.df %>%
      mutate(Signif = case_when(p_value < 0.001 ~ "***",
                                p_value < 0.01 ~  "**",
                                p_value < 0.05 ~ "*",
                                TRUE ~ " "))
    for (env_param in predictors) {
      lm_matrix[env_param, var_name] <- res.tmp.df[env_param,'t_value']
      lm_matrix.p[env_param, var_name] <- res.tmp.df[env_param,'p_value']
    }
  } # End For-loop
  if (nrow(lm_results)>1) {
    lm_results <- lm_results[-1,]
    rownames(lm_results) <- 1:nrow(lm_results)
  }
  if (nrow(intermediate_results)>1) {
    intermediate_results <- intermediate_results[-1,]
    rownames(intermediate_results) <- 1:nrow(intermediate_results)
  }
  # Skip drawing heatmap of t values:
  message("Finished auto lm modelling via step-wise forward selection. The summary of lm() results for each otu/gene/variable were stored in $lm_summary, and intermediate forward selection records were saved in $intermediate_values.")
  res_list$Results = lm_results
  res_list$intermediate_values = intermediate_results
  res_list$t_values = lm_matrix
  res_list$p_values = lm_matrix.p
  return(res_list)
}


# This function draws ordination plot (tested on PCA and RDA)
# @param datapca The RDA object list by `rda()`
# @param env.df A data.frame stores all environmental parameters to fit on the plot. Please make sure this data frame contains only the parameters you want to add. Default is NULL and will not run envfit()
# @param envfit If true (then `env.df` must be true), add envfit() to plot and highlight those are significant. Default is False. Recommended for PCA (unconstraint) but not RDA (constraint) as RDA will show biplot arrows
# @param scaling The 'scaling' argument in plot() or rda() that determines how the scores should be scaled. Default is 3 = 'symmetric'
# @param sd.val This is a magic variable from Paul's func, which determine how diverse a species should be highlighted on plot.
# @param text.cex Font size for text on data points (not the axis or label size)
# @param sp.pch,sp.col,sp.cex Point shape, color and size for Speices points (e.g., OTUs or Genes)
# @param site.pch,site.col,site.cex Point shape, color and size for Sites points (e.g., Samples)
draw_rda_plot <- function(datapca, env.df = NULL, envfit = FALSE,
                          scaling = 3, sd.val = 3, text.cex = 1,
                          sp.pch = 3, sp.col = 'red', sp.cex = 0.6,
                          site.pch = 21, site.col = 'black', site.cex = 0.6){
  axis.perc_val <- round((100*eigenvals(datapca)[1:2]/datapca$tot.chi[[1]]), digits=2)
  plot(datapca, scaling = scaling, type = "n",
       xlab = paste(names(axis.perc_val)[1], " (", axis.perc_val[1],"%)", sep = ""),
       ylab = paste(names(axis.perc_val)[2], " (", axis.perc_val[2],"%)", sep = ""))
  # Add all points to Species but colored in grey (i.e., OTU/genes) --
  points(datapca, dis='sp', pch=sp.pch, col='grey', cex = sp.cex, scaling = 3) # Add points ('x' shape) for species (e.g., OTUs or genes)
  
  # Add points to Sites (i.e., samples) --
  points(datapca, scaling = scaling, display ="sites", pch = site.pch, bg=site.col, cex=site.cex, col="grey15")
  text(datapca, scaling = scaling, display ="sites", cex = text.cex)
  
  # Highlight species if they are important --
  species_val <- data.frame(scores(datapca,scaling=scaling)$sp)
  species_to_highlight <- c(which(species_val[,1] > sd.val * sd(scores(datapca)$sp[,1])),
                            which(species_val[,1] < 0 - (sd.val * sd(scores(datapca)$sp[,1]))),
                            which(species_val[,2] > sd.val * sd(scores(datapca)$sp[,2])),
                            which(species_val[,2] < 0 - (sd.val * sd(scores(datapca)$sp[,2])))) # Not know why the threshold is based on non-scaled values but this is what from Paul's function.
  # species_to_highlight <- c(which(species_val[,1] > sd.val * sd(species_val[,1])),
  #                           which(species_val[,1] < 0 - sd.val * sd(species_val[,1])),
  #                           which(species_val[,2] > sd.val * sd(species_val[,2])),
  #                           which(species_val[,2] < 0 - sd.val * sd(species_val[,2]))) # This is what I think that truely makes sense
  if (length(species_to_highlight) > 0) {
    species_to_highlight_val <- species_val[unique(sort(species_to_highlight)),]
    points(species_to_highlight_val, pch = sp.pch, col = sp.col, cex = sp.cex)
    text(species_to_highlight_val, labels = rownames(species_to_highlight_val), cex = text.cex, col = sp.col)
  }
  
  # Add arrow biplot for regression coefficients, either via envfit() or RDA constraintion
  if (envfit){
    if (!is.null(env.df)){
      ev <- envfit(datapca, env.df, choices = 1:2, scaling = scaling, na.rm = TRUE)
      plot(ev, add = TRUE, cex = .6, col = "grey40")
      ev_sig <- envfit(datapca, env.df[,names(which(ev$vectors$pvals < 0.05))], choices = 1:2, scaling = scaling, na.rm = TRUE)
      plot(ev_sig, add = TRUE, cex = .8, col = "blue")
    } else {message("Warning: You asked to fit env factors but metadata (`env.df`) are not provided.")}
  } else {
    text(datapca, dis="cn", scaling=scaling, col = "blue")
  }
  cat("Done. You may want to add legends yourself. Example legend command for samples (i.e., Sites): \n")
  print("legend('topleft', legend = levels(your_sample_group), pch = site.pch, pt.cex = site.cex, pt.bg = 'grey15')")
}

