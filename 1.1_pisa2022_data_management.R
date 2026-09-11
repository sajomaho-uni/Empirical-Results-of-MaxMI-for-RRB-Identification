# START OF SCRIPT ----

# ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- #

# This script contains R code for the initial data management of the PISA 2022
# math item response data

# Author:
# Santeri Holopainen
# Turku Research Institute for Learning Analytics (TRILA)
# University of Turku

# Date: 11.9.2026

# R version: 4.5.0

# ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- #

# 1 Setup ----

# Clear the environment
rm(list = ls())

# Load the following packages
library(tidyverse)
library(tidylog)
library(crayon)
library(haven)

# Settings for tidylog
crayon = function(x) cat(green(x), sep = "\n")
options("tidylog.display" = list(crayon))
rm(crayon)

# 2 Initial data management ----

# The PISA 2022 cognitive item data is available at the OECD website by clicking the link
# "Download the datasets" under "PISA data > Public Use Files (PUFs)".
# The dataset's name in SPSS format is "CY08MSP_STU_COG.SAV".
# Link to the website: https://www.oecd.org/en/data/datasets/pisa-2022-database.html

# The meta data regarding the variables is stored in the Codebook excel file, 
# which can be accessed by clicking the link "Codebook" under  "Codebook and compendia". 
# The Codebook excel file is named "CY08MSP_CODEBOOK_27thJune24.xlsx".
# The excel file has many Sheets. The correct sheet is named "CY08MSP_STU_COG".

# The variable naming convention was retrieved from the technical report,
# which is available at https://doi.org/10.1787/01820d6d-en
# Especially Annex Table 12.A.2 (p. 266) was useful for this.

# Import the cognitive item data (this might take a while)
d <- read_sav("C:/Users/sjholo/Documents/Research/PhD dissertation/Article 3/Data analyses/Real-world example (PISA 2022)/CY08MSP_STU_COG.SAV")
d_backup <- d # backup

# Select only Finnish test-takers
d <- d %>% filter(CNT == "FIN")

# . 2.1 Construct the item map ----

# The item map will be used for extracting correct data for analysis from the main data file.

# The technical report is useful here:
# Annex Table 12.A.2 (p. 266) contains the naming convention of the variables.

# We are interested in the computer-based math items.
# These variables' names start with "CM", "CMA", "DM", or "DMA".
# The remaining characters are (in this order)
# - three numbers (unique item identifier), 
# - "Q",
# - a two-digit numeric item part code, and
# - additional information of the variable (e.g., "S", "SA", "SB", etc. for scored responses or "TT" for total timing)
# Example: "MA105Q01S" indicates the scored response for the first part of mathematics item 105.
item_map <- 
  data.frame(col_name = colnames(d)[str_detect(colnames(d), "^[DC]M.*Q0[1-5].*")]) %>% # get all variable names that match the description above
  mutate(item = str_remove_all(str_extract(col_name, "^[DC]M.*Q0[1-5]"), "^[DC]MA|^[DC]M|Q0[1-5]"), # extract the unique item identifier
         item_part = str_extract(col_name, "Q0[1-5]"), # extract the item part
         measure = str_remove(col_name, "^[DC]M.*Q0[1-5]"), # extract the type of measure (e.g., is it scored response or total time or something else?)
         full_name = str_extract(col_name, "M.*Q0[1-5]")) %>% # extract full name of the item (unique identifier + item part)
  mutate(id = 1) %>%
  pivot_wider(id_cols = c(item, item_part, full_name), names_from = measure, values_from = id) # finally, we pivot wide w.r.t. to the measure

# Set missing values to zeroes
item_map[is.na(item_map)] <- 0

# From the resulting data frame we know which items have which measures in the whole data file.
# The data set's dimensions should be 234 rows and 28 columns.

# We are interested only in items which have...
# - an automatically scored response,
# - the actual raw response, and
# - the total timing.
# => Leave only those items which have "S", "R", and "TT" in the item map
item_map <- 
  item_map %>% filter(S == 1 & R == 1 & TT == 1)

# The resulting data set should have 99 rows.

# Note!
# Item 112 part Q02 is the only item that has the following measures:
# HD
# N
# P
# Q
# WD
# This item is weird. According to the PISA 2022 code book, it does not have any data in the main file.
# => We remove it.
item_map <- item_map %>% filter(item != "112")

# The resulting data set should have 98 rows.

# At this point, let's do a check that all the selected items indeed are scored either correct or incorrect.
# => Find the min and max score for each item
item_map$min_score <- numeric(nrow(item_map))
item_map$max_score <- numeric(nrow(item_map))
for(i in 1:nrow(item_map)){
  vec_i <- na.omit(unname(unlist(d[,str_detect(colnames(d), paste0("^C", item_map$full_name[i], "S$"))])))
  if(length(vec_i) > 0){
    item_map$min_score[i] <- min(vec_i)
    item_map$max_score[i] <- max(vec_i)
  }else{
    item_map$min_score[i] <- 0
    item_map$max_score[i] <- 0
  }
}
sum(item_map$max_score != 1)
sum(item_map$min_score != 0)
# For all of the remaining items, min score is 0 and max score is 1

# Find the response options for each item
item_map$response_options <- numeric(nrow(item_map))
item_map$n_options <- numeric(nrow(item_map))
for(i in 1:nrow(item_map)){
  vec_i <- sort(unique(na.omit(unname(unlist(d[,str_detect(colnames(d), paste0("^D", item_map$full_name[i], "R$"))])))))
  if(length(vec_i) > 0){
    item_map$response_options[i] <- paste(sort(unique(vec_i)), collapse = "|")
    item_map$n_options[i] <- length(vec_i)
  }else{
    item_map$response_options[i] <- NA
    item_map$n_options[i] <- NA
  }
}

# Finally, remove irrelevant variables.
item_map <- item_map %>% select(item:full_name, response_options:n_options)

# The final item map should be a data frame of 98 rows and 5 columns

# Save the item map
save(item_map, file = "item_map_pisa2022.rda")

# . 2.2 Extract relevant data for analysis from the main file ----

# Extract data for the student IDs, scored responses, raw responses and 
# total timing from each of the 98 items.
temp_data <- list()
for(i in 1:nrow(item_map)){
  item_i <- item_map$item[i]
  name_i <- item_map$full_name[i]
  d_i <- d[,str_detect(colnames(d), paste0("CNTSTUID|", "^[DC]", name_i, "S$|", "^[DC]", name_i, "R$|", "^[DC]", name_i, "TT$"))]
  colnames(d_i)[str_detect(colnames(d_i), ".*S$")] <- "Correctness"
  colnames(d_i)[str_detect(colnames(d_i), ".*R$")] <- "Response"
  colnames(d_i)[str_detect(colnames(d_i), ".*TT$")] <- "Time"
  d_i <- na.omit(d_i)
  d_i$Correctness <- as.numeric(d_i$Correctness)
  d_i$Response <- as.numeric(d_i$Response)
  d_i$Time <- as.numeric(d_i$Time)
  temp_data[[name_i]] <- d_i
}
# The resulting "temp_data" object is a list of data frames each of which contains...
# - student IDs ("CNTSTUID"),
# - scored responses ("Correctness"),
# - raw responses ("Response"), and
# - total timing ("Time")
# for each item.

# Save the "temp_data" object on top of "d".
d <- temp_data

# Remove temporary objects from the environment
rm(temp_data, d_i, vec_i, item_i, name_i, i)

# If all went well, remove the backup data file:
rm(d_backup)

# Save the final data for analysis
save(d, file = "d_pisa2022.rda")

# END OF SCRIPT ----