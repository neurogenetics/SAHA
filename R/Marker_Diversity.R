#' Marker Diversity Function
#'
#' This function calculates the Shannon diversity index for markers within clusters and categorizes the diversity into "High", "Within 95% CI", and "Low" categories.
#'
#' @param ann An object containing annotation data.
#' 
#' @return A ggplot object showing the Shannon diversity by cluster.
#' @importFrom dplyr case_when
#' @importFrom tidyr spread
#' @importFrom vegan diversity
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual geom_hline theme_bw labs
#' @importFrom tibble column_to_rownames
#' @importFrom stats lm confint
#' @export

Marker_Diversity = function(ann){
   ann1=ann@ann1

   df = data.frame(ann1$query[,c("gene","cluster")])

   df$incidence = 1

   # Reshape using spread (assuming values for each id and time combination)
   wide_data <- df %>%
      spread(gene, incidence)  # Spreads "value" based on "time" for each "id"

   # Print the resulting wide dataframe
   wide_data <- wide_data %>%
      column_to_rownames("cluster")

   # Print the resulting dataframe
   #print(wide_data)



   # Calculate Shannon diversity for each cluster
   wide_data <- replace(wide_data, is.na(wide_data), 0)
   shannon_diversity <- diversity(wide_data, index = "shannon")

   shannon_diversity=data.frame(shannon_diversity)
   shannon_diversity$cluster=rownames(shannon_diversity)

   # # 99% sure you don't need this. When swithced to custom with pre-named clusters doesn't work... deal if we need to later.
   # shannon_diversity$cluster=as.numeric(shannon_diversity$cluster)


   model = lm(shannon_diversity ~ 1, data = shannon_diversity)
   CI=confint(model, level = 0.95)

   shannon_diversity$category = factor(
      case_when(
         shannon_diversity$shannon_diversity > CI[2] ~ "High Marker Diversity",
         shannon_diversity$shannon_diversity > CI[1] ~ "Within 95% CI",
         TRUE ~ "Low Marker Diversity"
      ))
   print(shannon_diversity)

   # Print or further analyze the shannon_diversity object (numeric vector)
   return(ggplot(shannon_diversity,aes(x=cluster,y=shannon_diversity))+
             geom_point(aes(color = category))+
             scale_color_manual(values = c("firebrick4","dodgerblue4","gray75"))+
             geom_hline(yintercept = mean(shannon_diversity$shannon_diversity))+
             geom_hline(yintercept = CI[1],linetype = 2)+
             geom_hline(yintercept = CI[2],linetype = 2)+
             theme_bw()+  # Threshold line 2
             labs(title = "Marker Diversity by Cluster", y = "Shannon Diveristy", x = "Cluster ID",
                  color = "Marker Diversity") + theme(axis.text.x = element_text(angle = 90)) + theme(axis.text = element_text(hjust = 0.95,
    vjust = 0.5)))  # Rename legend title

   }
