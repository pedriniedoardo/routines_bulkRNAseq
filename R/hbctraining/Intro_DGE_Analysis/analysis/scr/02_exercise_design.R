# libraries ---------------------------------------------------------------
library(tidyverse)
library(rafalib)
library(caret)

# build the matrix --------------------------------------------------------
# Create the data frame
metadata <- data.frame(
  sample = paste0("sample", 1:12),
  treatment = rep(c("A", "B", "C"), each = 4),
  sex = rep(c("F", "F", "M", "M"), times = 3),
  replicate = rep(1:4, times = 3)
) %>%
  mutate(across(treatment:replicate,as.factor))

str(metadata)

# build the model matrix
X <- model.matrix(~ treatment + sex,data = metadata)
X
imagemat(X)

# test if the matrix is co-linear
dim(X)
qr(X)$rank

# notice that If we add also the replicate, the gender and the replicate are colinera
X2 <- model.matrix(~ treatment + sex + replicate,data = metadata)
X2
imagemat(X2)

# test if the matrix is co-linear
dim(X2)
qr(X2)$rank

# find where is colinear
# Find the linear combinations
combo_info <- findLinearCombos(X2)

# Print the results
print(combo_info)

X2[,combo_info$linearCombos[[1]]]

# Fill in the RNA isolation column of the metadata table. Since we can only prepare 2 samples at a time and we have 12 samples total, you will need to isolate RNA in 6 batches. In the RNA isolation column, enter one of the following values for each sample: group1, group2, group3, group4, group5, group6. Make sure to fill in the table so as to avoid confounding by batch of RNA isolation.
metadata_test01 <- metadata %>%
  mutate(batch = c(1,2,3,4,5,6,1,2,3,4,5,6)) %>%
  # mutate(batch = c(1,1,2,2,3,3,4,4,5,5,6,6)) %>%
  mutate(across(treatment:batch,as.factor))

X3 <- model.matrix(~ treatment + sex + batch,data = metadata_test01)
X3
imagemat(X3)

# test if the matrix is co-linear
dim(X3)
qr(X3)$rank

# check the distribution fo the treamtnets per batch and gender
metadata_test01 %>%
  group_by(treatment,batch) %>%
  summarise(n = n())

metadata_test01 %>%
  group_by(sex,batch) %>%
  summarise(n = n())

metadata_test01 %>%
  group_by(sex,treatment) %>%
  summarise(n = n())

# BONUS: To perform the RNA isolations more quickly, you devote two researchers to perform the RNA isolations. Create a researcher column and fill in the researchers’ initials for the samples they will prepare: use initials AB or CD.
metadata_test02 <- metadata_test01 %>%
  mutate(user = c(1,2,1,2,1,2,2,1,2,1,2,1)) %>%
  # mutate(user = rep(c(1,2),6)) %>%
  mutate(across(treatment:user,as.factor))

X4 <- model.matrix(~ treatment + sex + batch + user,data = metadata_test02)
X4
imagemat(X4)

# test if the matrix is co-linear
dim(X4)
qr(X4)$rank

# solution provided -------------------------------------------------------
# solution provided by the authors
metadata_solution <- data.frame(
  sample = paste0("sample", 1:12),
  treatment = rep(c("A", "B", "C"), each = 4),
  sex = rep(c("F", "M"), each = 2, times = 3),
  replicate = rep(1:4, times = 3),
  `batch` = c("group1", "group2", "group3", "group4", 
                      "group3", "group5", "group1", "group6", 
                      "group4", "group6", "group2", "group5"),
  user = c("AB", "CD", "AB", "CD", "AB", "CD", "AB", "CD", 
                 "CD", "AB", "CD", "AB"),
  stringsAsFactors = FALSE, # Good practice to avoid factors by default
  check.names = FALSE # Allows the space in "RNA isolation"
)

# Print the resulting table
metadata_solution

X5 <- model.matrix(~ treatment + sex + batch + user,data = metadata_solution)
X5
imagemat(X5)

# test if the matrix is co-linear
dim(X5)
qr(X5)$rank


