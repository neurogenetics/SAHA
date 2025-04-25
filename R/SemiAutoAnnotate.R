#' Semi-automated annotation of query clusters
#'
#' This function provides a semi-automated approach for annotating query clusters.
#' Users are prompted to manually assign names to clusters based on visualizations
#' of marker-based or marker-free data.
#'
#' @param ann An object containing SAHA analysis results.
#' @param data_type The type of data to use for annotation ("Markers", "AvgExp", or "Both").
#' @param refine A data frame containing previously refined annotations (optional).
#' @param existing A data frame containing a saved, but incomplete semi_auto result.
#'
#' @return A data frame containing the original and new cluster names.
#'
#' @importFrom dplyr filter, arrange
#' @importFrom ggpubr ggarrange
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom grid grid.grabExpr
#'
#' @export
SemiAutoAnnotate = function(ann, data_type = NULL, refine = NULL, existing = NULL) {
  if (is.null(data_type)) {
    cat("Thank you for choosing the semi-automated approach. If you wish to proceed, please select which data_type to auto-annotate using: Markers, AvgExp, or Both.")
    return(NULL)
  }
  if (is.null(existing)) {
    cat("Thank you for choosing the semi-automated approach. To progress through each annotation, use the enter key. If you want to save your annotations, use the Enter key on a blank entry. If you want to quit the function, use your ESC key.")
  }
  
  if (data_type == "Markers") {
    temp = ann@results$marker_based$dotplot_all
    hand_names = data.frame(unique(temp$data$cluster))
    colnames(hand_names)[1] = "old_names"
    hand_names$new_names = ""
    
    if (!is.null(existing)) {
      hand_names$new_names = existing$new_names
      
      if (!is.null(refine)) {
        for (row in 1:nrow(hand_names)) {
          clust = hand_names$old_names[row]
          if (hand_names$new_names[row] == "") {
            match = refine$best_match[refine$cluster == clust]
            if (length(match) == 1 && match != "INCONCLUSIVE") {
              hand_names$new_names[row] = match
            }
          }
        }
        todo <- hand_names[hand_names$new_names == "", ]
        tokeep <- hand_names[hand_names$new_names != "", ]
      } else {
        todo = hand_names[hand_names$new_names == "", ]
      }
    } else if (!is.null(refine)) {
      todo <- hand_names %>% filter(old_names %in% unique(refine$cluster[refine$best_match == "INCONCLUSIVE"]))
      tokeep <- hand_names %>% filter(old_names %in% unique(refine$cluster[refine$best_match != "INCONCLUSIVE"]))
      tokeep$new_names = refine[refine$best_match != "INCONCLUSIVE", "best_match"]
    } else {
      todo = hand_names
    }
    
    for (i in todo$old_names) {
      temp2 <- ggplot(subset(temp$data, cluster == i), aes(x = celltype, y = cluster, alpha = as.numeric(-log(pvalue,10)))) +
        geom_point(aes(size = as.numeric(prop), color = sig)) +
        scale_size_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
        labs(title = "Every Cell Type by Cluster", x = " ", y = "Cluster") +
        scale_color_manual(values = c("black", "red")) +
        theme_bw() +
        theme(legend.position = "none", axis.text.x = element_text(angle=90, hjust=1))
      print(temp2)
      x = readline(paste0("What would you like to name cluster ", i, ": "))
      
      if (x == "") {
        message("Annotation paused. Progress saved.")
        hand_names$new_names[hand_names$old_names %in% todo$old_names] <- todo$new_names
        return(hand_names)
      }
      
      todo[todo$old_names == i, "new_names"] = x
    }
    
    if (!is.null(refine)) {
      hand_names = rbind(todo, tokeep)
      hand_names = hand_names %>% arrange(old_names)
    } else if (!is.null(existing)) {
      hand_names = existing
      hand_names$new_names[hand_names$old_names %in% todo$old_names] <- todo$new_names
    } else {
      hand_names = todo
    }
    return(hand_names)
  }
  
  if (data_type == "AvgExp") {
    ann3 = ann@results$marker_free$corr
    rownames(ann3) = gsub("^query\\.", "", rownames(ann3))
    colnames(ann3) = gsub("^db\\.", "", colnames(ann3))
    hand_names = data.frame(rownames(ann3))
    colnames(hand_names)[1] = "old_names"
    hand_names$new_names = ""
    
    if (!is.null(existing)) {
      hand_names$new_names = existing$new_names
      
      if (!is.null(refine)) {
        for (row in 1:nrow(hand_names)) {
          clust = hand_names$old_names[row]
          if (hand_names$new_names[row] == "") {
            match = refine$best_match[refine$cluster == clust]
            if (length(match) == 1 && match != "INCONCLUSIVE") {
              hand_names$new_names[row] = match
            }
          }
        }
        todo <- hand_names[hand_names$new_names == "", ]
        tokeep <- hand_names[hand_names$new_names != "", ]
      } else {
        todo = hand_names[hand_names$new_names == "", ]
      }
    } else if (!is.null(refine)) {
      todo <- hand_names %>% filter(old_names %in% unique(refine$cluster[refine$best_match == "INCONCLUSIVE"]))
      tokeep <- hand_names %>% filter(old_names %in% unique(refine$cluster[refine$best_match != "INCONCLUSIVE"]))
      tokeep$new_names = refine[refine$best_match != "INCONCLUSIVE", "best_match"]
    } else {
      todo = hand_names
    }
    
    for (i in todo$old_names) {
      temp = ann3[rownames(ann3) == i, ]
      mat = data.matrix(temp)
      print(Heatmap(mat,
                    name = "Pearson correlation",
                    cluster_columns = FALSE,
                    row_names_side = "left",
                    row_names_gp = gpar(fontsize = 11),
                    column_names_side = "bottom",
                    column_names_gp = gpar(fontsize = 11),
                    row_title = "query cluster",
                    row_title_rot = 90,
                    row_title_side = "left",
                    row_title_gp = gpar(fontface = "bold"),
                    column_title = "db",
                    column_title_rot = 0,
                    column_title_side = "bottom",
                    column_title_gp = gpar(fontface = "bold")
      ))
      x = readline(paste0("What would you like to name cluster ", i, ": "))
      
      if (x == "") {
        message("Annotation paused. Progress saved.")
        hand_names$new_names[hand_names$old_names %in% todo$old_names] <- todo$new_names
        return(hand_names)
      }
      
      todo[todo$old_names == i, "new_names"] = x
    }
    
    if (!is.null(refine)) {
      hand_names = rbind(todo, tokeep)
      hand_names = hand_names %>% arrange(old_names)
    } else if (!is.null(existing)) {
      hand_names = existing
      hand_names$new_names[hand_names$old_names %in% todo$old_names] <- todo$new_names
    } else {
      hand_names = todo
    }
    return(hand_names)
  }
  
  if (data_type == "Both") {
    ann2 = ann@results$marker_based$dotplot_all
    hand_names1 = data.frame(unique(ann2$data$cluster))
    colnames(hand_names1)[1] = "old_names"
    ann3 = ann@results$marker_free$corr
    rownames(ann3) = gsub("^query\\.", "", rownames(ann3))
    colnames(ann3) = gsub("^db\\.", "", colnames(ann3))
    hand_names2 = data.frame(rownames(ann3))
    colnames(hand_names2)[1] = "old_names"
    hand_names = merge(hand_names1, hand_names2, by = "old_names")
    
    if (nrow(hand_names) == 0) {
      warning("No matching cluster names found in marker-based and marker-free analysis. Consider running separately or renaming clusters.")
    }
    
    hand_names = hand_names %>% arrange(old_names)
    hand_names$new_names = ""
    
    if (!is.null(existing)) {
      hand_names$new_names = existing$new_names
      
      if (!is.null(refine)) {
        for (row in 1:nrow(hand_names)) {
          clust = hand_names$old_names[row]
          if (hand_names$new_names[row] == "") {
            match = refine$best_match[refine$cluster == clust]
            if (length(match) == 1 && match != "INCONCLUSIVE") {
              hand_names$new_names[row] = match
            }
          }
        }
        todo <- hand_names[hand_names$new_names == "", ]
        tokeep <- hand_names[hand_names$new_names != "", ]
      } else {
        todo = hand_names[hand_names$new_names == "", ]
      }
    } else if (!is.null(refine)) {
      todo <- hand_names %>% filter(old_names %in% unique(refine$cluster[refine$best_match == "INCONCLUSIVE"]))
      tokeep <- hand_names %>% filter(old_names %in% unique(refine$cluster[refine$best_match != "INCONCLUSIVE"]))
      tokeep$new_names = refine[refine$best_match != "INCONCLUSIVE", "best_match"]
    } else {
      todo = hand_names
    }
    
    for (i in todo$old_names) {
      temp2 <- ggplot(subset(ann2$data, cluster == i), aes(x = celltype, y = cluster, alpha = as.numeric(-log(pvalue,10)))) +
        geom_point(aes(size = as.numeric(prop), color = sig)) +
        scale_size_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
        labs(title = "Every Cell Type by Cluster", x = " ", y = "Cluster") +
        scale_color_manual(values = c("black", "red")) +
        theme_bw() +
        theme(legend.position = "none", axis.text.x = element_text(angle=90, hjust=1))
      
      temp = ann3[rownames(ann3) == i, ]
      mat = data.matrix(temp)
      p2 = Heatmap(mat,
                   name = "Pearson correlation",
                   cluster_columns = FALSE,
                   row_names_side = "left",
                   row_names_gp = gpar(fontsize = 11),
                   column_names_side = "bottom",
                   column_names_gp = gpar(fontsize = 11),
                   row_title = "query cluster",
                   row_title_rot = 90,
                   row_title_side = "left",
                   row_title_gp = gpar(fontface = "bold"),
                   column_title = "db",
                   column_title_rot = 0,
                   column_title_side = "bottom",
                   column_title_gp = gpar(fontface = "bold"))
      p2_grob <- grid.grabExpr(draw(p2))
      print(ggarrange(temp2, p2_grob, nrow = 2))
      
      x = readline(paste0("What would you like to name cluster ", i, ": "))
      
      if (x == "") {
        message("Annotation paused. Progress saved.")
        hand_names$new_names[hand_names$old_names %in% todo$old_names] <- todo$new_names
        return(hand_names)
      }
      
      todo[todo$old_names == i, "new_names"] = x
    }
    
    if (!is.null(refine)) {
      hand_names = rbind(todo, tokeep)
      hand_names = hand_names %>% arrange(old_names)
    } else if (!is.null(existing)) {
      hand_names = existing
      hand_names$new_names[hand_names$old_names %in% todo$old_names] <- todo$new_names
    } else {
      hand_names = todo
    }
    return(hand_names)
  }
}