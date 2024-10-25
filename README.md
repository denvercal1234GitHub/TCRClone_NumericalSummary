
# AIMS

To be able to generate tables that summarize the numerical values of clonal expansion and sharing after performing TCR clonal analysis using scRepertoire or other packages on the Seurat object. 

# PREPARE INPUT FROM SEURAT OBJECT WITH TCR ASSIGNMENT 

For 2 aims (assess TCR clonal expansion and clonal sharing numerically), it starts with a Seurat object after `scRepertoire()`. 

```{r}
## Prepare input for expansion_summary() function and analyze_Ctaa_sharing() from a Seurat object that has CTaa info attached from scRepertoire package with cells (rows) without TCR info removed
## This input requires the presence of 2 columns in the metadata of Seurat: CTaa is the clone sequences; cluster_ID is cluster assignment from Seurat for each clone

## SeuratObj@meta.data -> SeuratObj_metaData


####EXAMPLE TO ILLUSTRATE THE FUNCTIONS

# Sample data frame as example for a metadata of Seurat Object 
SeuratObj_metaData <- data.frame(
  cluster_ID = c("T25", "T25", "T25","T25","T25", "T25","T25","T25","T25","T26","T26", "T26","T26", "T27", "T27", "T27", "T27"),
  CTaa = c("K","K","A", "A", "A","B", "B","B","B","Q", "D","B","B", "C", "E", "F", "A"),
  stringsAsFactors = FALSE
)

#### This example has cluster T25 having clone "A" and clone "B" that expanded at least once! 
#### T26 has clone "B" expanded twice (which is counted as 1 clone being expanded). 
#### T27 clones are all single clones, i.e., 0 clone expanded. 
SeuratObj_metaData
```


# AIM 1: SUMMARIZE CLONAL EXPANSION FOR EACH CLUSTER

`expansion_summary()` will summarise for each cluster, how many clones a cluster has that expanded at least X times. EX: If Expanded_Count = 1, it means this cluster has 1 clone that expanded at least once. If Expanded_Count = 2, it means this cluster has 2 clones that each of them had expanded at least once. And so on. Thus, the higher the Expanded_Count, the more expanded the cluster is as it means the cluster has more clones that at least expanded X times. 

`expansion_summary()` need an input as a dataframe that counts the occurrence of each unique clone. 

```{r}
## Step 1: Obtain frequency from SeuratObj_metaData
SeuratObj_metaData %>% group_by(CTaa, cluster_ID) %>% summarise(count = n()) %>% arrange(cluster_ID) -> SeuratObj_countOccurenceOfEachUniqueClone_df

## Step 2: Get unique cluster IDs
unique_clusters <- unique(SeuratObj_countOccurenceOfEachUniqueClone_df$cluster_ID)

## Step 3: Run expansion_summary() 
expansion_summary <- lapply(unique_clusters, function(cluster) {
  cluster_data <- SeuratObj_countOccurenceOfEachUniqueClone_df %>%
    filter(cluster_ID == cluster)
  
  # Count how many times a clone expanded at least once (count > 1). Change per your own threshold
  expanded_count <- sum(cluster_data$count > 1)
  
  # Return a summary for this cluster
  data.frame(
    Cluster = cluster,
    Expanded_Count = expanded_count
  )
})

## Step 4: Combine the summaries into a single data frame
expansion_summary_df <- do.call(rbind, expansion_summary)

## Step 5: Print the summary table. 
#### This example has cluster T25 having clone "A" and clone "B" and clone "K" that expanded at least once! 
#### T26 has clone "B" expanded twice (which is counted as 1 clone being expanded). 
#### T27 clones are all single clones, i.e., 0 clone expanded. 
expansion_summary_df
```


# AIM 2.  SUMMARIZE WHICH CLONES FROM WHICH CLUSTERS SHARED WITH WHICH OTHER CLUSTERS AND WHICH CLONES NOT SHARED WITH ANY CLUSTERS 

`analyze_Ctaa_sharing()` will export a list of 2 dataframes as results (`ctaa_analysis_result`)

***To report how many clones, i.e. the `nrow` that are unique (not shared with anyone) from a cluster, e.g., T1, and their clone sequences*** 

`ctaa_analysis_result$Unique %>% dplyr::filter(Cluster == "T1")`

***To report how many clone(s), i.e. the `nrow` that are shared from this cluster with any other clusters, and the sequences of the clone(s)***

`ctaa_analysis_result$Shared %>% dplyr::filter(Cluster == "T1")`


```{r}
# Load necessary library
library(dplyr)

# Function to analyze CTaa shared values and specify which clusters share them
analyze_Ctaa_sharing <- function(data) {
  # Step 1: Get unique clusters
  unique_clusters <- unique(data$cluster_ID)

  # Step 2: Initialize lists to store results
  shared_results <- list()
  unique_results <- list()

  # Step 3: Loop through each unique cluster
  for (cluster in unique_clusters) {
    # Get the CTaa values for the current cluster
    current_ctaa <- data %>%
      filter(cluster_ID == cluster) %>%
      pull(CTaa)

    # Initialize variables to track shared CTaa values
    shared_ctaa <- NULL

    # Step 4: Check against all other clusters
    for (other_cluster in unique_clusters) {
      if (cluster != other_cluster) {
        other_ctaa <- data %>%
          filter(cluster_ID == other_cluster) %>%
          pull(CTaa)

        # Find shared CTaa values
        shared_current <- intersect(current_ctaa, other_ctaa)
        if (length(shared_current) > 0) {
          shared_ctaa <- unique(c(shared_ctaa, shared_current))
        }
      }
    }

    # Find unique CTaa values for the current cluster
    unique_current <- setdiff(current_ctaa, shared_ctaa)

    # Store results for shared CTaa values
    if (length(shared_ctaa) > 0) {
      for (value in shared_ctaa) {
        # Create a list to store the names of clusters sharing the current CTaa value
        clusters_shared <- c(cluster)  # Start with the current cluster

        # Count occurrences across all clusters sharing the CTaa
        clone_count_in_receiving_end <- 0

        for (other_cluster in unique_clusters) {
          if (other_cluster != cluster) {
            count_in_other_cluster <- sum(data$CTaa == value & data$cluster_ID == other_cluster)
            if (count_in_other_cluster > 0) {
              clusters_shared <- c(clusters_shared, other_cluster)
              clone_count_in_receiving_end <- clone_count_in_receiving_end + count_in_other_cluster
            }
          }
        }

        # Create a data frame entry for the shared CTaa value
        shared_entry <- data.frame(
          Cluster = cluster,
          Ctaa = value,
          Clusters_Shared = paste(clusters_shared[clusters_shared != cluster], collapse = ", "),
          Clone_Count_in_receivingEnd = clone_count_in_receiving_end,  # Updated column name
          stringsAsFactors = FALSE
        )

        # Append to the list of shared results
        shared_results[[length(shared_results) + 1]] <- shared_entry
      }
    }

    # Store unique CTaa values for the current cluster
    if (length(unique_current) > 0) {
      for (value in unique_current) {
        # Count occurrences of the unique CTaa in the cluster
        count_unique <- sum(data$CTaa == value & data$cluster_ID == cluster)

        unique_entry <- data.frame(
          Cluster = cluster,
          Ctaa = value,
          Clone_Count_inClusterOfInterest = count_unique,
          stringsAsFactors = FALSE
        )

        unique_results[[length(unique_results) + 1]] <- unique_entry
      }
    }
  }

  # Step 5: Combine the shared results into a data frame
  shared_df <- do.call(rbind, shared_results)

  # Combine unique results into a data frame
  unique_df <- do.call(rbind, unique_results)

  return(list(Shared = shared_df, Unique = unique_df))
}

# Run the analysis
ctaa_analysis_result <- analyze_Ctaa_sharing(SeuratObj_metaData)

  
# View the shared results with cluster information (DOESN'T INCLUDE UNIQUE CLONES)
### Note:
### In forward direction: T25 share clone B with T26 and clone B appear in T26 twice. That is clone B expanded twice in T26. That is T25 only has 1 clone (clone B) shared and it share that with T26.  
### In reverse direction: T26 share clone B with T25 and clone B appear in T25 four times. That is clone B expanded 4 times in T25. 
ctaa_analysis_result$Shared 

### Can just do ctaa_analysis_result$Shared %>% dplyr::filter(Cluster == "T25") to show that how many clones T25 shared with other clusters. 
ctaa_analysis_result$Shared %>% dplyr::filter(Cluster == "T25")

```


```{r}
########### ONLY SHOW UNIQUE CLONES FOR EACH CLUSTER 
## Note: T25 has clone K that is not shared with anyone, and clone K appears twice in T25!
## Note: T26 has clone Q and D are not shared with any one, and clone Q appears once in T26, and clone D also appears once in T26. Hence, they are listed and counted as 1 each 
ctaa_analysis_result$Unique
```

