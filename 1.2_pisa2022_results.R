# START OF SCRIPT ----

# ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- # ---- #

# This script contains R code for reproducing the results of employing the MaxMI 
# methods to the PISA 2022 math item response data

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
library(ggpubr)
library(mixtools)

# Settings for tidylog
crayon = function(x) cat(green(x), sep = "\n")
options("tidylog.display" = list(crayon))
rm(crayon)

# 2 Employ the RRB identification methods across all items ----

# . 2.1 Setup ----

# Clear the environment
rm(list = ls())

# Read the user-defined functions from "1.0_pisa2022_functions.R"
source("1.0_pisa2022_functions.R")

# Import the item map
load("item_map_pisa2022.rda")

# Import the cleaned PISA 2022 math item response data
load("d_pisa2022.rda")

# Bind the data
d <- bind_rows(d, .id = "item")
# The resulting data frame should have 110,521 rows and 5 columns

# Mutate response time from milliseconds to seconds
d <- d %>% mutate(Time = Time/1000)

# . 2.2 MaxMI-22 ----

# Loop through each item and employ the function "thr_MaxMI_x2y2_f"
thresholds_pisa2022_MaxMI22 <- list()
for(i in 1:length(item_map$item)){
  item_i <- item_map$full_name[i]
  d_i <- d %>% filter(item == item_i)
  thresholds_pisa2022_MaxMI22[[i]] <- thr_MaxMI_x2y2_f(d_k = d_i, n_cut = 10, limit = T)
}

# save
save(thresholds_pisa2022_MaxMI22, file = "thresholds_pisa2022_MaxMI22.rda")

# . 2.3 MaxMI-2Q ----

# Loop through each item and employ the function "thr_MaxMI_x2yq_f"
thresholds_pisa2022_MaxMI2Q <- list()
for(i in 1:length(item_map$item)){
  item_i <- item_map$full_name[i]
  d_i <- d %>% filter(item == item_i)
  thresholds_pisa2022_MaxMI2Q[[i]] <- thr_MaxMI_x2yq_f(d_k = d_i, n_cut = 10)
}

# save
save(thresholds_pisa2022_MaxMI2Q, file = "thresholds_pisa2022_MaxMI2Q.rda")

# . 2.4 MaxMI-32 ----

# Loop through each item and employ the function "thr_MaxMI_x3y2_f"
thresholds_pisa2022_MaxMI32 <- list()
for(i in 1:length(item_map$item)){
  item_i <- item_map$full_name[i]
  d_i <- d %>% filter(item == item_i)
  thresholds_pisa2022_MaxMI32[[i]] <- thr_MaxMI_x3y2_f(d_k = d_i)
}

# save
save(thresholds_pisa2022_MaxMI32, file = "thresholds_pisa2022_MaxMI32.rda")

# . 2.5 CUMP ----

# Loop through each item and employ the function "thr_CUMP_f"
thresholds_pisa2022_CUMP <- list()
for(i in 1:length(item_map$item)){
  item_i <- item_map$full_name[i]
  g_i <- 1/item_map$n_options[i]
  d_i <- d %>% filter(item == item_i)
  thresholds_pisa2022_CUMP[[i]] <- thr_CUMP_f(d_k = d_i, g_k = g_i, n_cut = 10)
}

# save
save(thresholds_pisa2022_CUMP, file = "thresholds_pisa2022_CUMP.rda")

# . 2.6 MLN ----

# Loop through each item and employ the function "thr_MLN_f"
thresholds_pisa2022_MLN <- list()
for(i in 1:length(item_map$item)){
  item_i <- item_map$full_name[i]
  d_i <- d %>% filter(item == item_i)
  thresholds_pisa2022_MLN[[i]] <- thr_MLN_f(d_k = d_i, seed_k = i)
}

# save
save(thresholds_pisa2022_MLN, file = "thresholds_pisa2022_MLN.rda")

rm(d_i, i, item_i, g_i)

# 3 Results presented in the main text ----

# . 3.1 Setup ----

# Clear the environment
rm(list = ls())

# Import the thresholds
load("thresholds_pisa2022_MaxMI22.rda")
load("thresholds_pisa2022_MaxMI32.rda")
load("thresholds_pisa2022_MaxMI2Q.rda")
load("thresholds_pisa2022_CUMP.rda")
load("thresholds_pisa2022_MLN.rda")

# Import the item map
load("item_map_pisa2022.rda")

# Import the cleaned PISA 2022 math item response data
load("d_pisa2022.rda")

# Bind the data
d <- bind_rows(d, .id = "item")
# The resulting data frame should have 110,521 rows and 5 columns

# Mutate response time from milliseconds to seconds
d <- d %>% mutate(Time = Time/1000)

# Extract the thresholds
temp1 <- 
  bind_rows(lapply(thresholds_pisa2022_MaxMI22, function(x) data.frame(thr = x$thr))) %>%
  mutate(item = item_map$full_name,
         method = "MaxMI-22")

temp2 <- 
  bind_rows(lapply(thresholds_pisa2022_MaxMI2Q, function(x) data.frame(thr = x))) %>%
  mutate(item = item_map$full_name,
         method = "MaxMI-2Q")

temp3 <- 
  bind_rows(lapply(thresholds_pisa2022_MaxMI32, function(x) data.frame(thr = x$thr1))) %>%
  mutate(item = item_map$full_name,
         method = "MaxMI-32 (1)")

temp4 <- 
  bind_rows(lapply(thresholds_pisa2022_MaxMI32, function(x) data.frame(thr = x$thr2))) %>%
  mutate(item = item_map$full_name,
         method = "MaxMI-32 (2)")

temp5 <- 
  bind_rows(lapply(thresholds_pisa2022_CUMP, function(x) data.frame(thr = x))) %>%
  mutate(item = item_map$full_name,
         method = "CUMP")

temp6 <- 
  bind_rows(lapply(thresholds_pisa2022_MLN, function(x) data.frame(thr = x$thr))) %>%
  mutate(item = item_map$full_name,
         method = "MLN")

# Join the thresholds to one data frame
thresholds <- bind_rows(temp1, temp2, temp3, temp4, temp5, temp6)
rm(temp1, temp2, temp3, temp4, temp5, temp6,
   thresholds_pisa2022_MaxMI22, 
   thresholds_pisa2022_MaxMI2Q,
   thresholds_pisa2022_MaxMI32,
   thresholds_pisa2022_CUMP, 
   thresholds_pisa2022_MLN)

# Determine the proportion of observations below the thresholds for each item
thresholds$p_x1 <- NA

# loop through each item
for(i in 1:length(item_map$item)){
  item_i <- item_map$full_name[i]
  d_i <- d %>% filter(item == item_i) %>% arrange(Time) %>% mutate(p_x1 = (1:n())/n())
  
  # loop through each method
  for(j in 1:length(unique(thresholds$method))){
    method_j <- unique(thresholds$method)[j]
    thr_ij <- thresholds$thr[thresholds$item == item_i & thresholds$method == method_j]
    p_x1_ij <- tail(d_i$p_x1[d_i$Time <= thr_ij], 1)
    thresholds$p_x1[thresholds$item == item_i & thresholds$method == method_j] <- p_x1_ij
  }
}
rm(i, d_i, item_i, j, method_j, thr_ij, p_x1_ij)

# . 3.2 Initial checks ----

# . . 3.2.1 Sample size range across the items ----

d %>%
  group_by(item) %>%
  summarize(n = n(), .groups = "drop") %>%
  summarize(min = min(n),
            max = max(n))

# . . 3.2.2 Proportion of computed thresholds ----

thresholds %>%
  group_by(method) %>%
  summarize(N = length(na.omit(thr)), .groups = "drop") %>%
  mutate(P = 100*N/98)

# . 3.3 Figure 1 in the main text ----

# . . 3.3.1 Setup ----

# Compute response accuracy, i.e., the proportion of correct responses, for each item
RAs <- 
  d %>%
  group_by(item) %>%
  summarize(accuracy = mean(Correctness),
            .groups = "drop")

# . . 3.3.2 The figure ----

# Form a data frame for drawing the figure
plot_data <- 
  thresholds %>% 
  # Filter out the second threshold of the MaxMI-32 method
  filter(method != "MaxMI-32 (2)") %>%
  # Compute and mutate necessary variables:
  mutate(# Compute a grouping variable which tells whether the threshold could be computed or not
         group = ifelse(!is.na(thr), 1, 2),
         group = factor(group, labels = c("Threshold was computed", "Threshold was not computed")),
         # Mutate p_x1 so that missing values are fixed to 0
         p_x1 = ifelse(is.na(p_x1), 0, p_x1),
         ) %>%
  # Join the thresholds with the response accuracies
  full_join(RAs, by = "item") %>%
  mutate(method = factor(method, levels = c("MaxMI-22", "MaxMI-2Q", "MaxMI-32 (1)", "CUMP", "MLN")))


# The figure
(p <- 
  ggplot(plot_data, mapping = aes(x = accuracy, y = p_x1, fill = group, shape = group))+
  geom_point(size = 2)+
  scale_shape_manual(values = c("Threshold was computed" = 21, "Threshold was not computed" = 22))+
  facet_wrap(~method)+
  labs(x = "Proportion of correct responses", 
       y = "Proportion of observations of response time\nthat are less than or equal to the threshold")+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 14)))

# save
ggsave(filename = "fig1.jpg",
       plot = p,
       width = 10000,
       height = 7000,
       units = "px",
       dpi = 1200)

rm(plot_data, RAs, p)


# . 3.4 Determine which three items could be used to demonstrate the behaviour of the MaxMI methods ----

# We want to select three items based on their ability to demonstrate the strengths and weaknesses of the different MaxMI methods.
# Ideally, in each item, one of the three methods successfully estimated a suitable threshold, while the others did not.
temp_table <- 
  thresholds %>%
  filter(str_detect(method, "MaxMI")) %>%
  pivot_wider(id_cols = item, names_from = method, values_from = p_x1) %>%
  # We do not need the second threshold of MaxMI-32
  select(-`MaxMI-32 (2)`) %>%
  # let's use 0.25 as a proportion threshold for determining whether the proportion 
  # of observations below the threshold is too large
  # => We create a grouping variable "cond":
  # cond = 0 if none of the methods produced a proportion less than 0.25
  # cond = 1 if MaxMI-22 succeeded
  # cond = 2 if MaxMI-2Q succeeded
  # cond = 3 if MaxMI-32 succeeded
  mutate(cond = ifelse(`MaxMI-22` < .25 & `MaxMI-2Q` > .25 & `MaxMI-32 (1)` > .25, 1, 0),
         cond = ifelse(`MaxMI-22` > .25 & `MaxMI-2Q` < .25 & `MaxMI-32 (1)` > .25, 2, cond),
         cond = ifelse(`MaxMI-22` > .25 & `MaxMI-2Q` > .25 & `MaxMI-32 (1)` < .25, 3, cond)) %>%
  filter(cond > 0)

# Visualize
plot_data <- 
  temp_table %>%
  pivot_longer(cols = `MaxMI-22`:`MaxMI-32 (1)`, values_to = "p_x1", names_to = "method") %>%
  mutate(cond = factor(cond,
                       labels = c("MaxMI-22 succeeded",
                                  "MaxMI-2Q succeeded",
                                  "MaxMI-32 (1) succeeded")))
ggplot(plot_data, mapping = aes(x = item, y = p_x1, fill = method))+
  geom_col(position = position_dodge())+
  facet_wrap(~cond, scales = "free_x")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "bottom")

# Let's choose MA131Q02, MA144Q01, and MA147Q04

rm(temp_table, plot_data)

# . 3.5 Table 1 in the main text ----

# Extract data for the three items
d_items <- d %>% filter(item %in% c("MA131Q02", "MA144Q01", "MA147Q04"))

# Extract thresholds for the three items
thresholds_items <-
  thresholds %>%
  filter(item %in% c("MA131Q02", "MA144Q01", "MA147Q04")) %>%
  mutate(thr = round(thr, 0)) %>%
  pivot_wider(id_cols = item, names_from = method, values_from = thr)

# The table:
d_items %>%
  group_by(item) %>%
  summarize(N = n(),
            Mean_Acc = round(mean(Correctness), 2),
            RT_Q25 = round(quantile(Time, .25), 0),
            RT_Md = round(median(Time), 0)) %>%
  inner_join(thresholds_items, by = "item")

rm(d_items, thresholds_items)

# . 3.6 Figure 2 in the main text ----

# Read the user-defined functions from "1.0_pisa2022_functions.R"
source("1.0_pisa2022_functions.R")

# . . 3.6.1 Plot A (item MA131Q02) ----

# Extract data for item MA131Q02
d_k <- 
  d %>% 
  filter(item == "MA131Q02") %>% 
  arrange(Time) %>% 
  mutate(id = 1:n()) %>%
  relocate(id, .before = item)

# Compute empirical MIRT and conditional probabilities, P(Y=1|X=x), x=0,1, 
# when X and Y are both binary:
d_k <-
  d_k %>%
  mutate(mirt_x2y2 = emp_mirt_x2y2_f(Correctness), # empirical MIRT
         p_x0 = 1 - (1:n())/n(), # P(X=0)=1-P(X=1), needed for q = P(Y = 1 | X = 0)
         p_x0y1 = mean(Correctness)-cumsum(Correctness)/n(), # P(X=0,Y=1)=P(Y=1)-P(Y=1,X=1), needed for q = P(Y = 1 | X = 0)
         p_y1_cond_x1 = emp_cump_f(Correctness), # r = P(Y=1|X=1), i.e., the CUMP curve
         p_y1_cond_x0 = p_x0y1 / pmax(p_x0, .Machine$double.xmin)) # q = P(Y=1|X=0)=P(X=0,Y=1)/P(X=0)
  
# Compute empirical MIRT and empirical conditional probabilities,
# P(Y=y|X=1), y = 1,2,3,4, when X is binary and Y is the raw response
d_k_x2yq <- 
  # Use the "emp_mirt_x2yq_f" function
  emp_mirt_x2yq_f(t = d_k$Time, y = d_k$Response) %>%
  # Compute P(Y=y|X=1) for each y = 1,2,3,4
  mutate(p_yq1_cond_x1 = p_x1y1 / p_x1,
         p_yq2_cond_x1 = p_x1y2 / p_x1,
         p_yq3_cond_x1 = p_x1y3 / p_x1,
         p_yq4_cond_x1 = p_x1y4 / p_x1) %>%
  rename(mirt_x2yq = mirt) %>%
  select(mirt_x2yq, p_yq1_cond_x1, p_yq2_cond_x1, p_yq3_cond_x1, p_yq4_cond_x1)

# Join
d_k <- 
  d_k %>%
  bind_cols(d_k_x2yq) %>%
  select(id, Time, Correctness, Response, 
         mirt_x2y2, p_y1_cond_x1, p_y1_cond_x0,
         mirt_x2yq, p_yq1_cond_x1, p_yq2_cond_x1, p_yq3_cond_x1, p_yq4_cond_x1)

# . . . 3.6.1.1 Plot A1 (results when both X and Y are binary) ----

# Form data frames for drawing the plot:

# 1) Data frame with empirical MIRT and conditional probabilities P(Y=1|X=x), x = 0,1
d_1 <- 
  d_k %>% 
  # Filter out first 10 observations to reduce spurious noise
  filter(id > 10) %>%
  select(Time, mirt_x2y2)

# 2) Data frame with frequencies in five-second time segments
d_2 <-
  d_k %>%
  # Compute the frequencies
  mutate(Time = cut(Time, seq(0, max(Time), 5),
                    labels = as.character(seq(5/2, max(Time)-5/2, 5))),
         Time = as.numeric(as.character(Time))) %>%
  group_by(Time) %>%
  summarize(n = n(), .groups = "drop") %>%
  # Scale the frequencies so that they can be plotted to same figure with MIRT
  mutate(ylim.prim.low = 0,
         ylim.prim.upp = max(d_1$mirt_x2y2),
         ylim.sec.low = 0,
         ylim.sec.upp = max(n),
         b = (ylim.prim.upp - ylim.prim.low)/(ylim.sec.upp - ylim.sec.low),
         a = ylim.prim.low - b*ylim.sec.low,
         n_s = a + b*n) %>%
  select(Time, n_s)

# 3) Data frame with empirical conditional probabilities P(Y=1|X=1) and P(Y=1|X=0)
d_3 <-
  d_k %>%
  # Filter out first 10 observations to reduce spurious noise
  filter(id > 10) %>%
  # Scale the probs so that they can be plotted to same figure with MIRT
  mutate(ylim.prim.low = 0,
         ylim.prim.upp = max(d_1$mirt_x2y2),
         ylim.sec.low = 0,
         ylim.sec.upp = 1,
         b = (ylim.prim.upp - ylim.prim.low)/(ylim.sec.upp - ylim.sec.low),
         a = ylim.prim.low - b*ylim.sec.low,
         p_y1_cond_x1_s = a + b*p_y1_cond_x1,
         p_y1_cond_x0_s = a + b*p_y1_cond_x0) %>%
  select(Time, p_y1_cond_x1_s, p_y1_cond_x0_s)

# 4) Data frame with the MaxMI-22 threshold
d_4 <- 
  thresholds %>%
  filter(item == "MA131Q02" & method == "MaxMI-22") %>%
  select(item, thr)

# Draw the plot
(p_A1 <- 
    ggplot()+
    geom_col(d_2, mapping = aes(x = Time, y = n_s), fill = "lightgrey", alpha = .5)+
    geom_line(d_1, mapping = aes(x = Time, y = mirt_x2y2), linewidth = .5, color = "orange")+
    geom_area(d_1, mapping = aes(x = Time, y = mirt_x2y2, fill = "mirt"), alpha = .2)+
    geom_line(d_3, mapping = aes(x = Time, y = p_y1_cond_x1_s), linewidth = 3, color = "white")+
    geom_line(d_3, mapping = aes(x = Time, y = p_y1_cond_x1_s, color = "F1 p_y1_cond_x1", linetype = "F1 p_y1_cond_x1"))+
    geom_line(d_3, mapping = aes(x = Time, y = p_y1_cond_x0_s, color = "F2 p_y1_cond_x0", linetype = "F2 p_y1_cond_x0"))+
    geom_point(d_4, mapping = aes(x = thr, y = 0),
               fill = "orange", color = "black", shape = 21, stroke = 1.5, size = 3)+
    scale_y_continuous(sec.axis = sec_axis(~ ./max(.), 
                                           breaks = c(0, 0.25, 0.5, 0.75, 1),
                                           labels = c("0.00", "0.25\n(=g)", "0.50", "0.75", "1.00"),
                                           name = expression(P[italic(N)](italic(U) == 1 ~ "|" ~ italic(X) == italic(x)))))+
    scale_fill_manual(values = c("orange"),
                      labels = c("Empirical MIRT"))+
    scale_color_manual(name = "curve",
                       values = c("blue", "darkblue"),
                       labels = c(expression(P[italic(N)](italic(U) == 1 ~ "|" ~ italic(X) == 1),
                                             P[italic(N)](italic(U) == 1 ~ "|" ~ italic(X) == 0))))+
    scale_linetype_manual(name = "curve",
                          values = c("solid", "longdash"),
                          labels = c(expression(P[italic(N)](italic(U) == 1 ~ "|" ~ italic(X) == 1),
                                                P[italic(N)](italic(U) == 1 ~ "|" ~ italic(X) == 0))))+
    labs(x = "RT or RT threshold (sec.)", y = "Empirical MIRT", 
         title = "A1")+
    theme_bw()+
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          legend.key.width = unit(3, "line"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          title = element_text(size = 16))+
    guides(color = guide_legend(nrow = 2))+
    coord_cartesian(xlim = c(0, 200)))

# . . . 3.6.1.2 Plot A2 (results when X is binary and Y is the raw response) ----

# Form data frames for drawing the plot:

# 1) Data frame with empirical MIRT
d_1 <-
  d_k %>%
  # Filter out first 10 observations to reduce spurious noise
  filter(id > 10) %>%
  select(Time, mirt_x2yq)

# 2) Data frame with frequencies in five-second time segments
d_2 <-
  d_k %>%
  # Compute the frequencies
  mutate(Time = cut(Time, seq(0, max(Time), 5),
                    labels = as.character(seq(5/2, max(Time)-5/2, 5))),
         Time = as.numeric(as.character(Time))) %>%
  group_by(Time) %>%
  summarize(n = n(), .groups = "drop") %>%
  # Scale the frequencies so that they can be plotted to same figure with MIRT
  mutate(ylim.prim.low = 0,
         ylim.prim.upp = max(d_1$mirt_x2yq),
         ylim.sec.low = 0,
         ylim.sec.upp = max(n),
         b = (ylim.prim.upp - ylim.prim.low)/(ylim.sec.upp - ylim.sec.low),
         a = ylim.prim.low - b*ylim.sec.low,
         n_s = a + b*n) %>%
  select(Time, n_s)

# 3) Data frame with empirical conditional probabilities P(Y=y|X=1)
d_3 <-
  d_k %>%
  # Filter out first 10 observations to reduce spurious noise
  filter(id > 10) %>%
  # Scale the probs so that they can be plotted to same figure with MIRT
  mutate(ylim.prim.low = 0,
         ylim.prim.upp = max(d_1$mirt_x2yq),
         ylim.sec.low = 0,
         ylim.sec.upp = 1,
         b = (ylim.prim.upp - ylim.prim.low)/(ylim.sec.upp - ylim.sec.low),
         a = ylim.prim.low - b*ylim.sec.low,
         p_y1_cond_x1_s = a+b*p_y1_cond_x1,
         p_yq1_cond_x1_s = a + b*p_yq1_cond_x1,
         p_yq2_cond_x1_s = a + b*p_yq2_cond_x1,
         p_yq3_cond_x1_s = a + b*p_yq3_cond_x1,
         p_yq4_cond_x1_s = a + b*p_yq4_cond_x1) %>%
  select(Time, p_y1_cond_x1_s, p_yq1_cond_x1_s, p_yq2_cond_x1_s, p_yq3_cond_x1_s, p_yq4_cond_x1_s)

# 4) Data frame with the MaxMI-2Q threshold
d_4 <- 
  thresholds %>%
  filter(item == "MA131Q02" & method == "MaxMI-2Q") %>%
  select(item, thr)

# Plot
(p_A2 <-
    ggplot()+
    geom_col(d_2, mapping = aes(x = Time, y = n_s), fill = "lightgrey", alpha = .5)+
    geom_line(d_1, mapping = aes(x = Time, y = mirt_x2yq), linewidth = .5, color = "orange")+
    geom_area(d_1, mapping = aes(x = Time, y = mirt_x2yq), fill = "orange", alpha = .2)+
    geom_line(d_3, mapping = aes(x = Time, y = p_y1_cond_x1_s), linewidth = 3, color = "white")+
    geom_line(d_3, mapping = aes(x = Time, y = p_yq1_cond_x1_s, color = "p_y1_cond_x1", linetype = "p_y1_cond_x1"))+
    geom_line(d_3, mapping = aes(x = Time, y = p_yq2_cond_x1_s, color = "p_y2_cond_x1", linetype = "p_y2_cond_x1"))+
    geom_line(d_3, mapping = aes(x = Time, y = p_yq3_cond_x1_s, color = "p_y3_cond_x1", linetype = "p_y3_cond_x1"))+
    geom_line(d_3, mapping = aes(x = Time, y = p_yq4_cond_x1_s, color = "p_y4_cond_x1", linetype = "p_y4_cond_x1"))+
    geom_point(d_4, mapping = aes(x = thr, y = 0), 
               fill = "orange", color = "black", shape = 21, stroke = 1.5, size = 3)+
    scale_y_continuous(sec.axis = sec_axis(~ ./max(.), 
                                           breaks = c(0, 0.25, 0.5, 0.75, 1),
                                           labels = c("0.00", "0.25\n(=g)", "0.50", "0.75", "1.00"),
                                           name = expression(P[italic(N)](italic(Y) == italic(y) ~ "|" ~ italic(X) == 1))))+
    scale_color_manual(name = "curve",
                       values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),
                       labels = c(expression(P[italic(N)](italic(Y) == 1 ~ "|" ~ italic(X) == 1), 
                                             P[italic(N)](italic(Y) == 2 ~ "|" ~ italic(X) == 1),
                                             P[italic(N)](italic(Y) == 3 ~ "|" ~ italic(X) == 1),
                                             P[italic(N)](italic(Y) == 4 ~ "|" ~ italic(X) == 1))))+
    scale_linetype_manual(name = "curve",
                          values = c("solid", "longdash", "dashed", "dotdash"),
                          labels = c(expression(P[italic(N)](italic(Y) == 1 ~ "|" ~ italic(X) == 1), 
                                                P[italic(N)](italic(Y) == 2 ~ "|" ~ italic(X) == 1),
                                                P[italic(N)](italic(Y) == 3 ~ "|" ~ italic(X) == 1),
                                                P[italic(N)](italic(Y) == 4 ~ "|" ~ italic(X) == 1))))+
    labs(x = "RT or RT threshold (sec.)", y = "Empirical MIRT", 
         shape = "Method", fill = "Method",
         title = "A2")+
    theme_bw()+
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          legend.key.width = unit(3, "line"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          title = element_text(size = 16))+
    guides(color = guide_legend(nrow = 2))+
    coord_cartesian(xlim = c(0, 200)))


# . . 3.6.2 Plot B (item MA144Q01) ----

# Extract data for item MA144Q01
d_k <- 
  d %>% 
  filter(item == "MA144Q01") %>% 
  arrange(Time) %>% 
  mutate(id = 1:n()) %>%
  relocate(id, .before = item)

# Compute empirical MIRT and conditional probabilities, P(Y=1|X=x), x=0,1, 
# when X and Y are both binary:
d_k <-
  d_k %>%
  mutate(mirt_x2y2 = emp_mirt_x2y2_f(Correctness), # empirical MIRT
         p_x0 = 1 - (1:n())/n(), # P(X=0)=1-P(X=1), needed for q = P(Y = 1 | X = 0)
         p_x0y1 = mean(Correctness)-cumsum(Correctness)/n(), # P(X=0,Y=1)=P(Y=1)-P(Y=1,X=1), needed for q = P(Y = 1 | X = 0)
         p_y1_cond_x1 = emp_cump_f(Correctness), # r = P(Y=1|X=1), i.e., the CUMP curve
         p_y1_cond_x0 = p_x0y1 / pmax(p_x0, .Machine$double.xmin)) # q = P(Y=1|X=0)=P(X=0,Y=1)/P(X=0)

# Compute empirical MIRT and empirical conditional probabilities,
# P(Y=y|X=1), y = 1,2,3,4, when X is binary and Y is the raw response
d_k_x2yq <- 
  # Use the "emp_mirt_x2yq_f" function
  emp_mirt_x2yq_f(t = d_k$Time, y = d_k$Response) %>%
  # Compute P(Y=y|X=1) for each y = 1,2,3,4
  mutate(p_yq1_cond_x1 = p_x1y1 / p_x1,
         p_yq2_cond_x1 = p_x1y2 / p_x1,
         p_yq3_cond_x1 = p_x1y3 / p_x1,
         p_yq4_cond_x1 = p_x1y4 / p_x1) %>%
  rename(mirt_x2yq = mirt) %>%
  select(mirt_x2yq, p_yq1_cond_x1, p_yq2_cond_x1, p_yq3_cond_x1, p_yq4_cond_x1)

# Join
d_k <- 
  d_k %>%
  bind_cols(d_k_x2yq) %>%
  select(id, Time, Correctness, Response, 
         mirt_x2y2, p_y1_cond_x1, p_y1_cond_x0,
         mirt_x2yq, p_yq1_cond_x1, p_yq2_cond_x1, p_yq3_cond_x1, p_yq4_cond_x1)

# . . . 3.6.2.1 Plot B1 (results when both X and Y are binary) ----

# Form data frames for drawing the plot:

# 1) Data frame with empirical MIRT and conditional probabilities P(Y=1|X=x), x = 0,1
d_1 <- 
  d_k %>% 
  # Filter out first 10 observations to reduce spurious noise
  filter(id > 10) %>%
  select(Time, mirt_x2y2)

# 2) Data frame with frequencies in five-second time segments
d_2 <-
  d_k %>%
  # Compute the frequencies
  mutate(Time = cut(Time, seq(0, max(Time), 5),
                    labels = as.character(seq(5/2, max(Time)-5/2, 5))),
         Time = as.numeric(as.character(Time))) %>%
  group_by(Time) %>%
  summarize(n = n(), .groups = "drop") %>%
  # Scale the frequencies so that they can be plotted to same figure with MIRT
  mutate(ylim.prim.low = 0,
         ylim.prim.upp = max(d_1$mirt_x2y2),
         ylim.sec.low = 0,
         ylim.sec.upp = max(n),
         b = (ylim.prim.upp - ylim.prim.low)/(ylim.sec.upp - ylim.sec.low),
         a = ylim.prim.low - b*ylim.sec.low,
         n_s = a + b*n) %>%
  select(Time, n_s)

# 3) Data frame with empirical conditional probabilities P(Y=1|X=1) and P(Y=1|X=0)
d_3 <-
  d_k %>%
  # Filter out first 10 observations to reduce spurious noise
  filter(id > 10) %>%
  # Scale the probs so that they can be plotted to same figure with MIRT
  mutate(ylim.prim.low = 0,
         ylim.prim.upp = max(d_1$mirt_x2y2),
         ylim.sec.low = 0,
         ylim.sec.upp = 1,
         b = (ylim.prim.upp - ylim.prim.low)/(ylim.sec.upp - ylim.sec.low),
         a = ylim.prim.low - b*ylim.sec.low,
         p_y1_cond_x1_s = a + b*p_y1_cond_x1,
         p_y1_cond_x0_s = a + b*p_y1_cond_x0) %>%
  select(Time, p_y1_cond_x1_s, p_y1_cond_x0_s)

# 4) Data frame with the MaxMI-22 threshold
d_4 <- 
  thresholds %>%
  filter(item == "MA144Q01" & method == "MaxMI-22") %>%
  select(item, thr)

# Draw the plot
(p_B1 <- 
    ggplot()+
    geom_col(d_2, mapping = aes(x = Time, y = n_s), fill = "lightgrey", alpha = .5)+
    geom_line(d_1, mapping = aes(x = Time, y = mirt_x2y2), linewidth = .5, color = "orange")+
    geom_area(d_1, mapping = aes(x = Time, y = mirt_x2y2, fill = "mirt"), alpha = .2)+
    geom_line(d_3, mapping = aes(x = Time, y = p_y1_cond_x1_s), linewidth = 3, color = "white")+
    geom_line(d_3, mapping = aes(x = Time, y = p_y1_cond_x1_s, color = "F1 p_y1_cond_x1", linetype = "F1 p_y1_cond_x1"))+
    geom_line(d_3, mapping = aes(x = Time, y = p_y1_cond_x0_s, color = "F2 p_y1_cond_x0", linetype = "F2 p_y1_cond_x0"))+
    geom_point(d_4, mapping = aes(x = thr, y = 0),
               fill = "orange", color = "black", shape = 21, stroke = 1.5, size = 3)+
    scale_y_continuous(sec.axis = sec_axis(~ ./max(.),  
                                           breaks = c(0, 0.25, 0.5, 0.75, 1),
                                           labels = c("0.00", "0.25\n(=g)", "0.50", "0.75", "1.00"),
                                           name = expression(P[italic(N)](italic(U) == 1 ~ "|" ~ italic(X) == italic(x)))))+
    scale_fill_manual(values = c("orange"),
                      labels = c("Empirical MIRT"))+
    scale_color_manual(name = "curve",
                       values = c("blue", "darkblue"),
                       labels = c(expression(P[italic(N)](italic(U) == 1 ~ "|" ~ italic(X) == 1),
                                             P[italic(N)](italic(U) == 1 ~ "|" ~ italic(X) == 0))))+
    scale_linetype_manual(name = "curve",
                          values = c("solid", "longdash"),
                          labels = c(expression(P[italic(N)](italic(U) == 1 ~ "|" ~ italic(X) == 1),
                                                P[italic(N)](italic(U) == 1 ~ "|" ~ italic(X) == 0))))+
    labs(x = "RT or RT threshold (sec.)", y = "Empirical MIRT", 
         title = "B1")+
    theme_bw()+
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          legend.key.width = unit(3, "line"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          title = element_text(size = 16))+
    guides(color = guide_legend(nrow = 2))+
    coord_cartesian(xlim = c(0, 200)))

# . . . 3.6.2.2 Plot B2 (results when X is binary and Y is the raw response) ----

# Form data frames for drawing the plot:

# 1) Data frame with empirical MIRT
d_1 <-
  d_k %>%
  # Filter out first 10 observations to reduce spurious noise
  filter(id > 10) %>%
  select(Time, mirt_x2yq)

# 2) Data frame with frequencies in five-second time segments
d_2 <-
  d_k %>%
  # Compute the frequencies
  mutate(Time = cut(Time, seq(0, max(Time), 5),
                    labels = as.character(seq(5/2, max(Time)-5/2, 5))),
         Time = as.numeric(as.character(Time))) %>%
  group_by(Time) %>%
  summarize(n = n(), .groups = "drop") %>%
  # Scale the frequencies so that they can be plotted to same figure with MIRT
  mutate(ylim.prim.low = 0,
         ylim.prim.upp = max(d_1$mirt_x2yq),
         ylim.sec.low = 0,
         ylim.sec.upp = max(n),
         b = (ylim.prim.upp - ylim.prim.low)/(ylim.sec.upp - ylim.sec.low),
         a = ylim.prim.low - b*ylim.sec.low,
         n_s = a + b*n) %>%
  select(Time, n_s)

# 3) Data frame with empirical conditional probabilities P(Y=y|X=1)
d_3 <-
  d_k %>%
  # Filter out first 10 observations to reduce spurious noise
  filter(id > 10) %>%
  # Scale the probs so that they can be plotted to same figure with MIRT
  mutate(ylim.prim.low = 0,
         ylim.prim.upp = max(d_1$mirt_x2yq),
         ylim.sec.low = 0,
         ylim.sec.upp = 1,
         b = (ylim.prim.upp - ylim.prim.low)/(ylim.sec.upp - ylim.sec.low),
         a = ylim.prim.low - b*ylim.sec.low,
         p_y1_cond_x1_s = a+b*p_y1_cond_x1,
         p_yq1_cond_x1_s = a + b*p_yq1_cond_x1,
         p_yq2_cond_x1_s = a + b*p_yq2_cond_x1,
         p_yq3_cond_x1_s = a + b*p_yq3_cond_x1,
         p_yq4_cond_x1_s = a + b*p_yq4_cond_x1) %>%
  select(Time, p_y1_cond_x1_s, p_yq1_cond_x1_s, p_yq2_cond_x1_s, p_yq3_cond_x1_s, p_yq4_cond_x1_s)

# 4) Data frame with the MaxMI-2Q threshold
d_4 <- 
  thresholds %>%
  filter(item == "MA144Q01" & method == "MaxMI-2Q") %>%
  select(item, thr)

# Plot
(p_B2 <-
    ggplot()+
    geom_col(d_2, mapping = aes(x = Time, y = n_s), fill = "lightgrey", alpha = .5)+
    geom_line(d_1, mapping = aes(x = Time, y = mirt_x2yq), linewidth = .5, color = "orange")+
    geom_area(d_1, mapping = aes(x = Time, y = mirt_x2yq), fill = "orange", alpha = .2)+
    geom_line(d_3, mapping = aes(x = Time, y = p_y1_cond_x1_s), linewidth = 3, color = "white")+
    geom_line(d_3, mapping = aes(x = Time, y = p_yq1_cond_x1_s, color = "p_y1_cond_x1", linetype = "p_y1_cond_x1"))+
    geom_line(d_3, mapping = aes(x = Time, y = p_yq2_cond_x1_s, color = "p_y2_cond_x1", linetype = "p_y2_cond_x1"))+
    geom_line(d_3, mapping = aes(x = Time, y = p_yq3_cond_x1_s, color = "p_y3_cond_x1", linetype = "p_y3_cond_x1"))+
    geom_line(d_3, mapping = aes(x = Time, y = p_yq4_cond_x1_s, color = "p_y4_cond_x1", linetype = "p_y4_cond_x1"))+
    geom_point(d_4, mapping = aes(x = thr, y = 0), 
               fill = "orange", color = "black", shape = 21, stroke = 1.5, size = 3)+
    scale_y_continuous(sec.axis = sec_axis(~ ./max(.), 
                                           breaks = c(0, 0.25, 0.5, 0.75, 1),
                                           labels = c("0.00", "0.25\n(=g)", "0.50", "0.75", "1.00"),
                                           name = expression(P[italic(N)](italic(Y) == italic(y) ~ "|" ~ italic(X) == 1))))+
    scale_color_manual(name = "curve",
                       values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),
                       labels = c(expression(P[italic(N)](italic(Y) == 1 ~ "|" ~ italic(X) == 1), 
                                             P[italic(N)](italic(Y) == 2 ~ "|" ~ italic(X) == 1),
                                             P[italic(N)](italic(Y) == 3 ~ "|" ~ italic(X) == 1),
                                             P[italic(N)](italic(Y) == 4 ~ "|" ~ italic(X) == 1))))+
    scale_linetype_manual(name = "curve",
                          values = c("solid", "longdash", "dashed", "dotdash"),
                          labels = c(expression(P[italic(N)](italic(Y) == 1 ~ "|" ~ italic(X) == 1), 
                                                P[italic(N)](italic(Y) == 2 ~ "|" ~ italic(X) == 1),
                                                P[italic(N)](italic(Y) == 3 ~ "|" ~ italic(X) == 1),
                                                P[italic(N)](italic(Y) == 4 ~ "|" ~ italic(X) == 1))))+
    labs(x = "RT or RT threshold (sec.)", y = "Empirical MIRT", 
         shape = "Method", fill = "Method",
         title = "B2")+
    theme_bw()+
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          legend.key.width = unit(3, "line"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          title = element_text(size = 16))+
    guides(color = guide_legend(nrow = 2))+
    coord_cartesian(xlim = c(0, 200)))

# . . 3.6.3 Plot C (item MA147Q04) ----

# Extract data for item MA147Q04
d_k <- 
  d %>% 
  filter(item == "MA147Q04") %>% 
  arrange(Time) %>% 
  mutate(id = 1:n()) %>%
  relocate(id, .before = item)

# Compute empirical MIRT and conditional probabilities, P(Y=1|X=x), x=0,1, 
# when X and Y are both binary:
d_k <-
  d_k %>%
  mutate(mirt_x2y2 = emp_mirt_x2y2_f(Correctness), # empirical MIRT
         p_x0 = 1 - (1:n())/n(), # P(X=0)=1-P(X=1), needed for q = P(Y = 1 | X = 0)
         p_x0y1 = mean(Correctness)-cumsum(Correctness)/n(), # P(X=0,Y=1)=P(Y=1)-P(Y=1,X=1), needed for q = P(Y = 1 | X = 0)
         p_y1_cond_x1 = emp_cump_f(Correctness), # r = P(Y=1|X=1), i.e., the CUMP curve
         p_y1_cond_x0 = p_x0y1 / pmax(p_x0, .Machine$double.xmin)) # q = P(Y=1|X=0)=P(X=0,Y=1)/P(X=0)

# Compute empirical MIRT and empirical conditional probabilities,
# P(Y=y|X=1), y = 1,2,3,4, when X is binary and Y is the raw response
d_k_x2yq <- 
  # Use the "emp_mirt_x2yq_f" function
  emp_mirt_x2yq_f(t = d_k$Time, y = d_k$Response) %>%
  # Compute P(Y=y|X=1) for each y = 1,2,3,4
  mutate(p_yq1_cond_x1 = p_x1y1 / p_x1,
         p_yq2_cond_x1 = p_x1y2 / p_x1,
         p_yq3_cond_x1 = p_x1y3 / p_x1,
         p_yq4_cond_x1 = p_x1y4 / p_x1) %>%
  rename(mirt_x2yq = mirt) %>%
  select(mirt_x2yq, p_yq1_cond_x1, p_yq2_cond_x1, p_yq3_cond_x1, p_yq4_cond_x1)

# Join
d_k <- 
  d_k %>%
  bind_cols(d_k_x2yq) %>%
  select(id, Time, Correctness, Response, 
         mirt_x2y2, p_y1_cond_x1, p_y1_cond_x0,
         mirt_x2yq, p_yq1_cond_x1, p_yq2_cond_x1, p_yq3_cond_x1, p_yq4_cond_x1)

# . . . 3.6.3.1 Plot C1 (results when both X and Y are binary) ----

# Form data frames for drawing the plot:

# 1) Data frame with empirical MIRT and conditional probabilities P(Y=1|X=x), x = 0,1
d_1 <- 
  d_k %>% 
  # Filter out first 10 observations to reduce spurious noise
  filter(id > 10) %>%
  select(Time, mirt_x2y2)

# 2) Data frame with frequencies in five-second time segments
d_2 <-
  d_k %>%
  # Compute the frequencies
  mutate(Time = cut(Time, seq(0, max(Time), 5),
                    labels = as.character(seq(5/2, max(Time)-5/2, 5))),
         Time = as.numeric(as.character(Time))) %>%
  group_by(Time) %>%
  summarize(n = n(), .groups = "drop") %>%
  # Scale the frequencies so that they can be plotted to same figure with MIRT
  mutate(ylim.prim.low = 0,
         ylim.prim.upp = max(d_1$mirt_x2y2),
         ylim.sec.low = 0,
         ylim.sec.upp = max(n),
         b = (ylim.prim.upp - ylim.prim.low)/(ylim.sec.upp - ylim.sec.low),
         a = ylim.prim.low - b*ylim.sec.low,
         n_s = a + b*n) %>%
  select(Time, n_s)

# 3) Data frame with empirical conditional probabilities P(Y=1|X=1) and P(Y=1|X=0)
d_3 <-
  d_k %>%
  # Filter out first 10 observations to reduce spurious noise
  filter(id > 10) %>%
  # Scale the probs so that they can be plotted to same figure with MIRT
  mutate(ylim.prim.low = 0,
         ylim.prim.upp = max(d_1$mirt_x2y2),
         ylim.sec.low = 0,
         ylim.sec.upp = 1,
         b = (ylim.prim.upp - ylim.prim.low)/(ylim.sec.upp - ylim.sec.low),
         a = ylim.prim.low - b*ylim.sec.low,
         p_y1_cond_x1_s = a + b*p_y1_cond_x1,
         p_y1_cond_x0_s = a + b*p_y1_cond_x0) %>%
  select(Time, p_y1_cond_x1_s, p_y1_cond_x0_s)

# 4) Data frame with the MaxMI-22 threshold
d_4 <- 
  thresholds %>%
  filter(item == "MA147Q04" & method == "MaxMI-22") %>%
  select(item, thr)

# Draw the plot
(p_C1 <- 
    ggplot()+
    geom_col(d_2, mapping = aes(x = Time, y = n_s), fill = "lightgrey", alpha = .5)+
    geom_line(d_1, mapping = aes(x = Time, y = mirt_x2y2), linewidth = .5, color = "orange")+
    geom_area(d_1, mapping = aes(x = Time, y = mirt_x2y2, fill = "mirt"), alpha = .2)+
    geom_line(d_3, mapping = aes(x = Time, y = p_y1_cond_x1_s), linewidth = 3, color = "white")+
    geom_line(d_3, mapping = aes(x = Time, y = p_y1_cond_x1_s, color = "F1 p_y1_cond_x1", linetype = "F1 p_y1_cond_x1"))+
    geom_line(d_3, mapping = aes(x = Time, y = p_y1_cond_x0_s, color = "F2 p_y1_cond_x0", linetype = "F2 p_y1_cond_x0"))+
    geom_point(d_4, mapping = aes(x = thr, y = 0),
               fill = "orange", color = "black", shape = 21, stroke = 1.5, size = 3)+
    scale_y_continuous(sec.axis = sec_axis(~ ./max(.), 
                                           breaks = c(0, 0.25, 0.5, 0.75, 1),
                                           labels = c("0.00", "0.25\n(=g)", "0.50", "0.75", "1.00"),
                                           name = expression(P[italic(N)](italic(U) == 1 ~ "|" ~ italic(X) == italic(x)))))+
    scale_fill_manual(values = c("orange"),
                      labels = c("Empirical MIRT"))+
    scale_color_manual(name = "curve",
                       values = c("blue", "darkblue"),
                       labels = c(expression(P[italic(N)](italic(U) == 1 ~ "|" ~ italic(X) == 1),
                                             P[italic(N)](italic(U) == 1 ~ "|" ~ italic(X) == 0))))+
    scale_linetype_manual(name = "curve",
                          values = c("solid", "longdash"),
                          labels = c(expression(P[italic(N)](italic(U) == 1 ~ "|" ~ italic(X) == 1),
                                                P[italic(N)](italic(U) == 1 ~ "|" ~ italic(X) == 0))))+
    labs(x = "RT or RT threshold (sec.)", y = "Empirical MIRT", 
         title = "C1")+
    theme_bw()+
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          legend.key.width = unit(3, "line"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          title = element_text(size = 16))+
    guides(color = guide_legend(nrow = 2))+
    coord_cartesian(xlim = c(0, 200)))

# . . . 3.6.3.2 Plot C2 (results when X is binary and Y is the raw response) ----

# Form data frames for drawing the plot:

# 1) Data frame with empirical MIRT
d_1 <-
  d_k %>%
  # Filter out first 10 observations to reduce spurious noise
  filter(id > 10) %>%
  select(Time, mirt_x2yq)

# 2) Data frame with frequencies in five-second time segments
d_2 <-
  d_k %>%
  # Compute the frequencies
  mutate(Time = cut(Time, seq(0, max(Time), 5),
                    labels = as.character(seq(5/2, max(Time)-5/2, 5))),
         Time = as.numeric(as.character(Time))) %>%
  group_by(Time) %>%
  summarize(n = n(), .groups = "drop") %>%
  # Scale the frequencies so that they can be plotted to same figure with MIRT
  mutate(ylim.prim.low = 0,
         ylim.prim.upp = max(d_1$mirt_x2yq),
         ylim.sec.low = 0,
         ylim.sec.upp = max(n),
         b = (ylim.prim.upp - ylim.prim.low)/(ylim.sec.upp - ylim.sec.low),
         a = ylim.prim.low - b*ylim.sec.low,
         n_s = a + b*n) %>%
  select(Time, n_s)

# 3) Data frame with empirical conditional probabilities P(Y=y|X=1)
d_3 <-
  d_k %>%
  # Filter out first 10 observations to reduce spurious noise
  filter(id > 10) %>%
  # Scale the probs so that they can be plotted to same figure with MIRT
  mutate(ylim.prim.low = 0,
         ylim.prim.upp = max(d_1$mirt_x2yq),
         ylim.sec.low = 0,
         ylim.sec.upp = 1,
         b = (ylim.prim.upp - ylim.prim.low)/(ylim.sec.upp - ylim.sec.low),
         a = ylim.prim.low - b*ylim.sec.low,
         p_y1_cond_x1_s = a+b*p_y1_cond_x1,
         p_yq1_cond_x1_s = a + b*p_yq1_cond_x1,
         p_yq2_cond_x1_s = a + b*p_yq2_cond_x1,
         p_yq3_cond_x1_s = a + b*p_yq3_cond_x1,
         p_yq4_cond_x1_s = a + b*p_yq4_cond_x1) %>%
  select(Time, p_y1_cond_x1_s, p_yq1_cond_x1_s, p_yq2_cond_x1_s, p_yq3_cond_x1_s, p_yq4_cond_x1_s)

# 4) Data frame with the MaxMI-22 threshold
d_4 <- 
  thresholds %>%
  filter(item == "MA147Q04" & method == "MaxMI-2Q") %>%
  select(item, thr)

# Plot
(p_C2 <-
    ggplot()+
    geom_col(d_2, mapping = aes(x = Time, y = n_s), fill = "lightgrey", alpha = .5)+
    geom_line(d_1, mapping = aes(x = Time, y = mirt_x2yq), linewidth = .5, color = "orange")+
    geom_area(d_1, mapping = aes(x = Time, y = mirt_x2yq), fill = "orange", alpha = .2)+
    geom_line(d_3, mapping = aes(x = Time, y = p_y1_cond_x1_s), linewidth = 3, color = "white")+
    geom_line(d_3, mapping = aes(x = Time, y = p_yq1_cond_x1_s, color = "p_y1_cond_x1", linetype = "p_y1_cond_x1"))+
    geom_line(d_3, mapping = aes(x = Time, y = p_yq2_cond_x1_s, color = "p_y2_cond_x1", linetype = "p_y2_cond_x1"))+
    geom_line(d_3, mapping = aes(x = Time, y = p_yq3_cond_x1_s, color = "p_y3_cond_x1", linetype = "p_y3_cond_x1"))+
    geom_line(d_3, mapping = aes(x = Time, y = p_yq4_cond_x1_s, color = "p_y4_cond_x1", linetype = "p_y4_cond_x1"))+
    geom_point(d_4, mapping = aes(x = thr, y = 0), 
               fill = "orange", color = "black", shape = 21, stroke = 1.5, size = 3)+
    scale_y_continuous(sec.axis = sec_axis(~ ./max(.), 
                                           breaks = c(0, 0.25, 0.5, 0.75, 1),
                                           labels = c("0.00", "0.25\n(=g)", "0.50", "0.75", "1.00"),
                                           name = expression(P[italic(N)](italic(Y) == italic(y) ~ "|" ~ italic(X) == 1))))+
    scale_color_manual(name = "curve",
                       values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),
                       labels = c(expression(P[italic(N)](italic(Y) == 1 ~ "|" ~ italic(X) == 1), 
                                             P[italic(N)](italic(Y) == 2 ~ "|" ~ italic(X) == 1),
                                             P[italic(N)](italic(Y) == 3 ~ "|" ~ italic(X) == 1),
                                             P[italic(N)](italic(Y) == 4 ~ "|" ~ italic(X) == 1))))+
    scale_linetype_manual(name = "curve",
                          values = c("solid", "longdash", "dashed", "dotdash"),
                          labels = c(expression(P[italic(N)](italic(Y) == 1 ~ "|" ~ italic(X) == 1), 
                                                P[italic(N)](italic(Y) == 2 ~ "|" ~ italic(X) == 1),
                                                P[italic(N)](italic(Y) == 3 ~ "|" ~ italic(X) == 1),
                                                P[italic(N)](italic(Y) == 4 ~ "|" ~ italic(X) == 1))))+
    labs(x = "RT or RT threshold (sec.)", y = "Empirical MIRT", 
         shape = "Method", fill = "Method",
         title = "C2")+
    theme_bw()+
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          legend.key.width = unit(3, "line"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          title = element_text(size = 16))+
    guides(color = guide_legend(nrow = 2))+
    coord_cartesian(xlim = c(0, 200)))

# . . 3.6.4 Combine the plots and draw the final figure ----

# Combine A1, B1, and C1
(p_1 <- ggarrange(p_A1, p_B1, p_C1, nrow = 3, legend = "bottom", common.legend = T))

# Combine A2, B2, and C2
(p_2 <- ggarrange(p_A2, p_B2, p_C2, nrow = 3, legend = "bottom", common.legend = T))

# The final figure
(p <- ggarrange(p_1, p_2, ncol = 2, legend = "bottom"))

# Save the figure
ggsave(filename = "fig2.jpg",
       plot = p,
       width = 10000,
       height = 11000,
       units = "px",
       dpi = 1200)

rm(p, p_1, p_2)


# . 3.7 Figure 3 in the main text ----

# Extract data for item MA147Q04
d_k <- 
  d %>% 
  filter(item == "MA147Q04") %>% 
  arrange(Time) %>% 
  mutate(id = 1:n()) %>%
  relocate(id, .before = item)

# Compute empirical MIRT when X is governed by two time thresholds and Y denotes correctness
d_k_x3y2 <- emp_mirt_x3y2_f(t = d_k$Time, y = d_k$Correctness)

# Form data frames for the plots

# 1) Data frame with time thresholds and empirical MIRT
d_1 <- d_k_x3y2 %>% select(t1,t2,mirt)

# 2) Data frame with the first MaxMI-32 threshold
d_2 <- 
  thresholds %>%
  filter(item == "MA147Q04" & str_detect(method, "MaxMI-32")) %>%
  pivot_wider(id_cols = item, names_from = method, values_from = thr) %>%
  rename(thr1 = `MaxMI-32 (1)`,
         thr2 = `MaxMI-32 (2)`) %>%
  select(item, thr1, thr2)

# Draw the figure
(p <- 
    ggplot()+
    geom_contour_filled(d_1, mapping = aes(x = t1, y = t2, z = mirt, fill = after_stat(level_mid)),
                        binwidth = .01)+
    scale_fill_viridis_c()+
    annotate(geom = "linerange", xmin = -100, xmax = d_2$thr1, y = d_2$thr2, linewidth = 1, color = "white")+
    annotate(geom = "linerange", ymin = -100, ymax = d_2$thr2, x = d_2$thr1, linewidth = 1, color = "white")+
    annotate(geom = "linerange", xmin = -100, xmax = d_2$thr1, y = d_2$thr2, linewidth = .5, linetype = "dashed")+
    annotate(geom = "linerange", ymin = -100, ymax = d_2$thr2, x = d_2$thr1, linewidth = .5, linetype = "dashed")+
    geom_point(d_2, mapping = aes(x = thr1, y = thr2), 
               shape = 21, color = "black", fill = "white", stroke = 1.2, size = 2.5)+
    labs(x = "First time threshold (sec.)", y = "Second time threshold (sec.)",
         fill = "Empirical\nMIRT")+
    theme_bw()+
    theme(legend.position = "inside",
          legend.position.inside = c(.85,.3),
          legend.background = element_rect(color = "black"),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12))+
    coord_cartesian(xlim = c(0, 200), ylim = c(0,200)))

# Save the figure
ggsave(filename = "fig3.jpg",
       plot = p,
       width = 6000,
       height = 5000,
       units = "px",
       dpi = 1200)

# END OF SCRIPT ----