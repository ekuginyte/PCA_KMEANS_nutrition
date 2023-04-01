##### SET UP #####

# Libraries required
library_list <- c("knitr", "class", "cluster", "colorspace", "corrplot", "dendextend", "kableExtra", "dplyr", "factoextra", "GGally", "ggcorrplot", "mclust", "readr", "scatterplot3d", "tidyverse", "viridis", "wesanderson", "cowplot", "gtsummary", "cellranger")

# Load the libraries
for (i in library_list) {
  library(i, character.only = TRUE)
}

# Set working directory
setwd("~/Documents/MT5758 Multivariate Analysis/Assignments/Nutrition_data/Nutrition")

# Load data
data <- read_csv("data.csv")





##### DATA WRANGLING #####

# Check for missing values in each column
missing <- colSums(is.na(data))

# Not numeric types
non_numeric_cols <- sapply(data, function(x) !is.numeric(x))
non_numeric_col_names <- names(data)[non_numeric_cols]

# Find number of categories of food groups
n <- unique(data$FoodGroup) %>% length()

# Convert the food groups variable to factors and save it as a separate vector
labels <- data$FoodGroup

# Delete the unwanted columns with character variables, delete USRDA variables
#   as they are basically duplicates of the nutrition values
data <- data %>% dplyr::select(-c(1:7, 31:45)) 

# Convert all the other variables to numeric
data <- lapply(data, function(x) if(is.numeric(x)) x else as.numeric(as.character(x)))

# Convert to data frame
data <- as.data.frame(data)




#### FREQUENCY PLOT

# Save the data frame
data_full <- cbind(data, labels)

# Find number of entries in each category
group_freq <- table(data_full$labels)

# Sorted food labels
fg <- c("Dairy, Eggs", "Spices", "Baby Foods", "Fats, Oils", "Poultry", 
        "Soups, Sauces", "Sausages, Luncheon", "Br Cereals", 
        "Snacks", "Fruit, Fr Juice", "Pork", "Vegetables",
        "Nuts, Seeds", "Beef", "Beverages", "Finfish, Shellf.", 
        "Legumes", "Lamb, Veal, Game", "Baked Products", 
        "Sweets", "Cereal Gr, Pasta", "Fast Foods", "Meals, Entr, Side",
        "American Native", "Restaurant F") %>% 
  sort()

# Set colours to be different for each cluster
my_col1 <- wesanderson::wes_palette(n = 1, "Zissou1")

# Convert the table to a data frame and abbreviate food group names
group_freq <- data.frame(Food_Group = fg,
                         Frequency = as.numeric(group_freq))

# Plot the frequencies
p1 <- ggplot(group_freq, aes(x = Food_Group, y = Frequency)) +
  geom_bar(stat = "identity", color = "white", alpha = 0.9, 
           fill = my_col1, alpha = 0.8) +
  labs(x = "Food Group", y = "Frequency") +
  ylim(c = 0, 1000) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.background = element_rect(fill = "#F5F5F5", color = NA)) +
  guides(fill = FALSE) +
  geom_rect(aes(xmin = as.numeric(Food_Group) - 0.5, 
                xmax = as.numeric(Food_Group) + 0.5, 
                ymin = 0, ymax = Frequency), 
            fill = NA, color = "black")





##### CORRELATION #####

# Reset the grid
par(mfrow = c(1, 1)) 

# Get correlation matrix and plot correlogram
cor_mat <- cor(data)

# Plot the matrix
corrplot.mixed(cor_mat, upper = "circle", addgrid = c("upper"), 
               order = "hclust", tl.col = "black", tl.pos = "lt", diag = "l", 
               number.font = 0.2, tl.cex = 0.5, number.cex = 0.55)






##### DATA SCALING AND CENTRE-ING #####

# Standardize the data to prepare for Principal Component Analysis
scaled_data <- scale(data, center = TRUE, scale = TRUE)

# Add the IDs and the food groups back to the scaled data set
foodGroup <- data.frame(foodGroup = labels)
dataf <- cbind(foodGroup, scaled_data)






##### PCA SELECTION #####

# Perform PCA on the scaled data
pca <- prcomp(scaled_data)

# Extract the principal component scores
PCA_components <- pca$x %>% 
  data.frame()

### PCA variance
# Compute scores and get eigenvalues
eig <- eigen(cov(scaled_data))
Z <- scaled_data %*% eig$vectors

# Compute cumulative proportion of variance explained
pc_cumul_prop_var <- cumsum(eig$values/sum(eig$values))


### Plot both plots

# Set up the plot grid
par(mfrow = c(1, 2))


# Plot the PCA
plot(pca, type = "l", 
     xlab = "Principal Component Index", 
     ylab = "Proportion of Variance Explained", 
     col = "brown2", lwd = 2)


# Plot cumulative proportion of variance explained
plot(pc_cumul_prop_var, xlab = "Principal Component Index", 
     ylab = "Cum. Prop. of Var. Explained", 
     ylim = c(0, 1),
     col = "brown2",
     lwd = 2)
abline(h = 0.9, lwd = 1.5, col = "azure4", lty = 2)
abline(v = 14, lwd = 1.5, col = "azure4", lty = 2)
text(14, 0.75, "14 PCA C.", col = "gray60", adj = c(0, -.1))

# Extract the first 14 principal components
pca_components <- pca$x[, 1:14]

# Reset the plot grid to default
par(mfrow = c(1, 1))

y <- foodGroup %>% as.vector()

# Plot biplot of the first two PCA components
fviz_pca_biplot(pca,
                label = "var",
                habillage = y)





### K-MEANS CLUSTERING

# Perform K-means clustering on the scaled data
kmeans_fit <- kmeans(scaled_data, centers = 9, nstart = 50, iter.max = 50)

# Get the cluster assignments for each observation
clustersAssigned <- kmeans_fit$cluster

# Add the cluster assignments to the data frame WITH PCA COMPONENTS, overwrite 
# the scaled data
dataf <- cbind(foodGroup, scaled_data, clustersAssigned)

### EVALUATE CLUSTERS

# Calculate ARI and FMI
ari <- mclust::adjustedRandIndex(dataf$foodGroup, dataf$clustersAssigned)
fmi <- FM_index_R(dataf$foodGroup, dataf$clustersAssigned,  assume_sorted_vectors = FALSE)

# Print the results
cat("Adjusted Rand Index (ARI):", ari, "\n")
cat("Fowlkes-Mallows Index (FMI):", fmi, "\n")


### PLOT K MEANS

# Colour palette
pal <- c("red", "yellow")

# Scatterplot matrices coloured by clusters
pairs(pca_components[, 11:14], pch = 20, cex = 0.8, col = pal[clustersAssigned])


# Silhouette
silkmeans_p <- fviz_nbclust(scaled_data, kmeans, linecolor = "brown2", method = "s")

silkmeans_p2 <- silkmeans_p +
  ggtitle(" ") +
  theme(plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))

# Total within sum of squares 
wsskmeans <- fviz_nbclust(scaled_data, kmeans, linecolor = "brown2", method = "wss") 

# Visualise silhouhette information
sil_p <- fviz_silhouette(silhouette(kmeans_fit$cluster, dist(scaled_data)),
                         print.summary = FALSE)

sil_p2 <- sil_p +
  ggtitle("Average Silhouette Width = 0.33") + 
  theme(plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))

# Set up the plot grid
par(mfrow = c(1, 1))
# Suggests optimal k = 9
# Suggests k = +- 4

# Together plots
cowplot::plot_grid(silkmeans_p2, sil_p2, ncol = 2, height = 150)


# Visualise kmeans clustering
kmeans_p <- fviz_cluster(kmeans_fit, scaled_data, outlier.color = "black", geom = "point") +
  theme_minimal() +
  ggtitle(NULL)

# Count the number of entries in each cluster
cluster_counts <- table(kmeans_fit$cluster)

# Create a data frame with the cluster counts and percentages
pie_data <- data.frame(cluster = names(cluster_counts),
                       count = as.numeric(cluster_counts),
                       percent = round(as.numeric(cluster_counts) / 
                                         sum(cluster_counts) * 100, 1))

# Plot the pie chart
kmeans_pie_p <- ggplot(pie_data, aes(x = "", y = count, fill = cluster)) +
  geom_bar(stat = "identity", width = 1) +
  geom_bar(stat = "identity", width = 1, 
           color = "white", alpha = 0, lwd = 0.5) +
  coord_polar(theta = "y") +
  theme_void() +
  ggtitle("Cluster Counts") +
  ggtitle(NULL) +
  labs(fill = "Cluster") +
  geom_text(aes(label = paste0(round(count/sum(count)*100), "%")), 
            position = position_stack(vjust = 0.5), 
            color = "grey30", fontface = "bold")

# Plot together
cowplot::plot_grid(kmeans_p, kmeans_pie_p, ncol = 2, height = 150)


### PLOT THE COMPOSITION OF CLUSTERS

# Recode foodGroup factor levels with shorter labels
dataf$foodGroup_short <- dplyr::recode(dataf$foodGroup,
                                       "American Indian/Alaska Native Foods" = "Native Foods",
                                       "Baby Foods" = "Baby Food",
                                       "Baked Products" = "Baked",
                                       "Beef Products" = "Beef",
                                       "Beverages" = "Drinks",
                                       "Breakfast Cereals" = "Cereals",
                                       "Cereal Grains and Pasta" = "Grains & Pasta",
                                       "Dairy, Eggs" = "Dairy & Eggs",
                                       "Fast Foods" = "Fast Food",
                                       "Fats, Oils" = "Fats & Oils",
                                       "Finfish and Shellfish Products" = "Fish & Shellfish",
                                       "Fruits and Fruit Juices" = "Fruits & Juices",
                                       "Lamb, Veal, and Game Products" = "Lamb & Game",
                                       "Legumes and Legume Products" = "Legumes",
                                       "Meals, Entrees, and Side Dishes" = "Meals & Sides",
                                       "Nut and Seed Products" = "Nuts & Seeds",
                                       "Pork Products" = "Pork",
                                       "Poultry Products" = "Poultry",
                                       "Restaurant Foods" = "Restaurants",
                                       "Sausages and Luncheon Meats" = "Sausages",
                                       "Snacks" = "Snacks",
                                       "Soups, Sauces, and Gravies" = "Soups & Sauces",
                                       "Spices" = "Spices",
                                       "Sweets" = "Sweets",
                                       "Vegetables and Vegetable Products" = "Vegetables")

# Summarise the percentages of each food group in clusters
summaryTable <- dataf %>%
  group_by(clustersAssigned, foodGroup) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  group_by(clustersAssigned) %>%
  # Shorten the titles of food groups to they fit in the table
  mutate(foodGroup_short = dplyr::recode(foodGroup,
                                         "clustersAssigned" = "Cluster",
                                         "American Indian/Alaska Native Foods" = "Native Am",
                                         "Baby Foods" = "Baby",
                                         "Baked Products" = "Baked",
                                         "Beef Products" = "Beef",
                                         "Beverages" = "Drinks",
                                         "Breakfast Cereals" = "Cereals",
                                         "Cereal Grains and Pasta" = "Grains",
                                         "Dairy and Egg Products" = "Dairy, Eggs ",
                                         "Fast Foods" = "Fast Food",
                                         "Fats and Oils" = "Fats, Oils",
                                         "Finfish and Shellfish Products" = "Fish",
                                         "Fruits and Fruit Juices" = "Fruit, Juice ",
                                         "Lamb, Veal, and Game Products" = "Lamb, Veal. ",
                                         "Legumes and Legume Products" = "Legumes",
                                         "Meals, Entrees, and Side Dishes" = "Meals",
                                         "Nut and Seed Products" = "Nuts, Seeds",
                                         "Pork Products" = "Pork",
                                         "Poultry Products" = "Poultry",
                                         "Restaurant Foods" = "Restaurant",
                                         "Sausages and Luncheon Meats" = "Sausages",
                                         "Snacks" = "Snacks",
                                         "Soups, Sauces, and Gravies" = "Soups",
                                         "Spices and Herbs" = "Spices",
                                         "Sweets" = "Sweets",
                                         "Vegetables and Vegetable Products" = "Vegetables")) %>%
  mutate(prop = round(n * 100 / sum(n), 2)) %>%
  dplyr::select(-n, -foodGroup) %>%
  spread(foodGroup_short, prop) %>%
  rename(Clusters = clustersAssigned) %>% 
  # Replace "NA" with "-"X
  mutate_all(~ifelse(. == "NA", "-", as.numeric(.))) %>% 
  # Convert summary table from list to data frame
  as.data.frame() 

# Get rid of NA's and add "-"
summaryTable[is.na(summaryTable)] <- "-"

# Check if the proportions have been calculated well
sum(summaryTable[2, 2:26], na.rm = TRUE)

cluster_majority <- character(nrow(summaryTable))
for (i in 1:nrow(summaryTable)) {
  cluster_majority[i] <- colnames(summaryTable)[which.max(summaryTable[i, ])]
}

# Change colours of the values for the summary table
for (i in i:ncol(summaryTable)) {
  summaryTable[, i] <- cell_spec(summaryTable[, i], color = ifelse(summaryTable[, i] > 20, "tomato1", "black"))
}


# Cut the tables into two as they are too wide
summary_t1 <- summaryTable[, 1:13]
summary_t2 <- summaryTable[, c(1, 14:26)]

# Print the tables
kable(summary_t1, align = "c", caption = "Cluster Composition. Proportion of Each Food Group in a Cluster (Part 1)") %>% 
  row_spec(0, angle = 90) %>% 
  kableExtra::kable_styling(
    latex_options = c("striped", "repeat_header"),
    stripe_color = "gray!7",
    font_size = 8,
    full_width = TRUE) %>% 
  kable_styling(latex_options = "HOLD_position")

kable(summary_t2, align = "c", caption = "Cluster Composition. Proportion of Each Food Group in a Cluster (Part 2)") %>% 
  row_spec(0, angle = 90)  %>% 
  kableExtra::kable_styling(
    latex_options = c("striped", "repeat_header"),
    stripe_color = "gray!7",
    font_size = 8,
    full_width = TRUE)  %>% 
  kable_styling(latex_options = "HOLD_position")


kable(summary_t2, align = "c", caption = "Summary Table") %>% 
  row_spec(0, angle = 90)  %>% 
  kableExtra::kable_styling(
    latex_options = c("striped", "repeat_header"),
    stripe_color = "gray!7",
    font_size = 8,
    full_width = TRUE) 





##### RESULTS OF THE CLUSTERING #####

# Group the data by cluster and foodGroup, count the frequency of each group, and summarize by the total count within each cluster
df_summary <- dataf %>%
  group_by(clustersAssigned, foodGroup) %>%
  summarize(count = n()) %>%
  mutate(total_count = sum(count)) %>%
  mutate(percentage = count / total_count * 100) %>%
  dplyr::select(clustersAssigned, foodGroup, percentage) %>%
  pivot_wider(names_from = clustersAssigned, 
              values_from = percentage, 
              values_fill = 0) %>%
  arrange(foodGroup) %>%
  pivot_longer(cols = -foodGroup, names_to = "Cluster", values_to = "Percentage")


### PLOT RESULTS

# Set colours to be different for each cluster
my_col1 <- wesanderson::wes_palette(n = 3, "Zissou1")
my_col2 <- wesanderson::wes_palette(n = 3, "FantasticFox1")
my_col3 <- wesanderson::wes_palette(n = 3, "Darjeeling1")

# Recode the names so they fit in the table
df_summary$foodGroup <- dplyr::recode(df_summary$foodGroup,
                                      "American Indian/Alaska Native Foods" = "Native Foods",
                                      "Baby Foods" = "Baby Food",
                                      "Baked Products" = "Baked",
                                      "Beef Products" = "Beef",
                                      "Beverages" = "Drinks",
                                      "Breakfast Cereals" = "Cereals",
                                      "Cereal Grains and Pasta" = "Grains & Pasta",
                                      "Dairy, Eggs" = "Dairy & Eggs",
                                      "Fast Foods" = "Fast Food",
                                      "Fats, Oils" = "Fats & Oils",
                                      "Finfish and Shellfish Products" = "Fish & Shellfish",
                                      "Fruits and Fruit Juices" = "Fruits & Juices",
                                      "Lamb, Veal, and Game Products" = "Lamb & Game",
                                      "Legumes and Legume Products" = "Legumes",
                                      "Meals, Entrees, and Side Dishes" = "Meals & Sides",
                                      "Nut and Seed Products" = "Nuts & Seeds",
                                      "Pork Products" = "Pork",
                                      "Poultry Products" = "Poultry",
                                      "Restaurant Foods" = "Restaurants",
                                      "Sausages and Luncheon Meats" = "Sausages",
                                      "Snacks" = "Snacks",
                                      "Soups, Sauces, and Gravies" = "Soups & Sauces",
                                      "Spices" = "Spices",
                                      "Sweets" = "Sweets",
                                      "Vegetables and Vegetable Products" = "Vegetables")

# Filter the data frame to include only clusters 1, 2, and 3
df1 <- df_summary %>% 
  filter(Cluster %in% c("1", "2", "3"))

# Create the side-by-side bar chart # 1-3 clusters
ggplot(df1, aes(x = Percentage, y = foodGroup, fill = Cluster)) +
  geom_bar(position = position_dodge(width = 0.9), alpha = 0.95, stat = "identity") +
  # set the color palette
  scale_fill_manual(values = my_col1) +  
  labs(x = "Percentage of Food Group Elements in Respective Clusters", 
       y = "Food Group (25 Total)", fill = "Cluster Number") +
  theme_minimal() +
  #  move legend to bottom
  theme(legend.position = "bottom",
        # add border to legend box
        legend.box.background = element_rect(color = "black"),  
        # increase size of legend text
        legend.text = element_text(size = 8), 
        # increase size of axis titles
        axis.title = element_text(size = 8))

# Filter the data frame to include only clusters 1, 2, and 3
df2 <- df_summary %>% 
  filter(Cluster %in% c("4", "5", "8"))

# Create the side-by-side bar chart # 4-6 clusters
ggplot(df2, aes(x = Percentage, y = foodGroup, fill = Cluster)) +
  geom_bar(position = position_dodge(width = 0.9), alpha = 0.95, stat = "identity") +
  scale_fill_manual(values = my_col2) + 
  labs(x = "Percentage of Food Group Elements in Respective Clusters", 
       y = "Food Group (25 Total)", fill = "Cluster Number") +
  theme_minimal() +
  theme(legend.position = "bottom",   
        legend.box.background = element_rect(color = "black"),
        legend.text = element_text(size = 8),  
        axis.title = element_text(size = 8)) 



# Filter the data frame to include only clusters 1, 2, and 3
df3 <- df_summary %>% 
  filter(Cluster %in% c("6", "7", "9"))

# Create the side-by-side bar chart # 7-9 clusters
ggplot(df3, aes(x = Percentage, y = foodGroup, fill = Cluster)) +
  geom_bar(position = position_dodge(width = 0.9), alpha = 0.95, stat = "identity") +
  scale_fill_manual(values = my_col3) + 
  labs(x = "Percentage of Food Group Elements in Respective Clusters", 
       y = "Food Group (25 Total)", fill = "Cluster Number") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.box.background = element_rect(color = "black"),
        legend.text = element_text(size = 8),
        axis.title = element_text(size = 8)) 






##### CLUSTER EXPLORATION #####

# Create a data frame for the main nutrition components to be visualised
cluster_nutr <- dataf[, c(2:7, 25)]

# Create an empty list to store the cluster means
cluster_means <- list()
# Loop through the clusters and calculate the mean values of each variable
for (i in 1:9){
  cluster_data <- cluster_nutr[dataf$clustersAssigned == i, ]
  cluster_mean <- colMeans(cluster_data)
  cluster_means[[i]] <- cluster_mean
}

# Combine the cluster means into a data frame
cluster_means_df <- do.call(rbind, cluster_means) %>% data.frame()

# Cluster column must be a factor
cluster_means_df$clustersAssigned <- as.factor(cluster_means_df$clustersAssigned)

# Create a parallel coordinates plot for the selected variables,
# with lines colored by cluster

# Pick colours
my_col4 <- viridis_pal(option = "magma")(10)

# Plot the parallel coordinates
ggparcoord(cluster_means_df, columns = 1:6, groupColumn = 7, order = "anyClass", 
           showPoints = TRUE, scale = "globalminmax") +
  scale_color_manual(values = my_col4) +
  ylab("Value") +
  xlab("Nutrient") +
  scale_fill_discrete(name = "Clusters") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        fill = "Clusters") +
  geom_line(size = 0.7) +
  theme_minimal() +
  labs(color = "Clusters")