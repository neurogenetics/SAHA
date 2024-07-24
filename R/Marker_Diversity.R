# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

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
   print(wide_data)



   # Calculate Shannon diversity for each cluster
   wide_data <- replace(wide_data, is.na(wide_data), 0)
   shannon_diversity <- diversity(wide_data, index = "shannon")

   shannon_diversity=data.frame(shannon_diversity)
   shannon_diversity$cluster=rownames(shannon_diversity)

   shannon_diversity$cluster=as.numeric(shannon_diversity$cluster)

   model = lm(shannon_diversity ~ 1, data = shannon_diversity)
   CI=confint(model, level = 0.95)

   shannon_diversity$category = factor(
      case_when(
         shannon_diversity$shannon_diversity > CI[2] ~ "High Marker Diversity",
         shannon_diversity$shannon_diversity > CI[1] ~ "Within 95% CI",
         TRUE ~ "Low Marker Diversity"
      ))

   # Print or further analyze the shannon_diversity object (numeric vector)
   return(ggplot(shannon_diversity,aes(x=cluster,y=shannon_diversity))+
             geom_point(aes(color = category))+
             scale_color_manual(values = c("firebrick4","dodgerblue4","gray75"))+
             geom_hline(yintercept = mean(shannon_diversity$shannon_diversity))+
             geom_hline(yintercept = CI[1],linetype = 2)+
             geom_hline(yintercept = CI[2],linetype = 2)+
             theme_bw()+  # Threshold line 2
             labs(title = "Marker Diversity by Cluster", y = "Shannon Diveristy", x = "Cluster ID",
                  color = "Marker Diversity"))  # Rename legend title

}
