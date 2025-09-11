# AIM ---------------------------------------------------------------------
# sample pull genes belonging to GO terms


# libraries ---------------------------------------------------------------
library(tidyverse)


# read in the sample table with GO terms ----------------------------------

# read in the list of terms provided by Aletta
TOI <- read_csv("../data/GOterms_iron_myelin.csv")

# msigdbr approach --------------------------------------------------------
library(msigdbr)

# pull the collection
msigdbr_collections()

# pull the specific annotation
# the table contains the terms and the genes
gene_sets <- msigdbr(species = "Homo sapiens", collection = "C5",subcollection = "GO:BP")
head(gene_sets)

test01_01_summary <- gene_sets %>%
  group_by(gs_exact_source) %>%
  summarise(n = n())

test01_01_summary

# check that the terms provided are present in the current dataset
TOI %>%
  left_join(test01_01_summary,by = c("GO_id" = "gs_exact_source")) %>%
  mutate(test = !is.na(n)) %>%
  print(n = 30)

gene_sets %>%
  filter(gs_exact_source == "GO:0110076")

# org.Hs.eg.db approach ---------------------------------------------------
library(AnnotationDbi)
library(org.Hs.eg.db)

# check the package version
packageVersion("org.Hs.eg.db")

# Get all the available keytypes from the human database
keytypes(org.Hs.eg.db)

# Get the first 6 gene symbols from the database
head(keys(org.Hs.eg.db, keytype="SYMBOL"))

# Get the first 6 GO IDs from the database
head(keys(org.Hs.eg.db, keytype="GO"))

# what is the difference between GO and GOALL
# GO (Standard) Content	Only current, non-obsolete GO terms Size	Smaller, as it is a curated subset.
# GOALL (All) Includes all GO terms, including obsolete ones Larger, as it contains the full history of the GO.
test02_01 <- AnnotationDbi::select(org.Hs.eg.db,
                                   keys = unique(TOI$GO_id),
                                   columns = c('SYMBOL'),
                                   keytype = "GO")

head(test02_01)

test02_02 <- AnnotationDbi::select(org.Hs.eg.db,
                                   keys = unique(TOI$GO_id),
                                   columns = c('SYMBOL'),
                                   keytype = "GOALL")

head(test02_02)

test02_01_summary <- test02_01 %>%
  group_by(GO) %>%
  summarise(n = n())

test02_02_summary <- test02_02 %>%
  group_by(GOALL) %>%
  summarise(n = n())

# summarise the genes
TOI %>%
  left_join(test02_01_summary,by = c("GO_id" = "GO")) %>%
  mutate(test = !is.na(n)) %>%
  print(n = 30)

# sample pull the genes
test02_01 %>%
  filter(GO %in% "GO:0010446")

# summarise the genes
TOI %>%
  left_join(test02_02_summary,by = c("GO_id" = "GOALL")) %>%
  mutate(test = !is.na(n)) %>%
  print(n = 30)

# sample pull the genes
test02_02 %>%
  filter(GOALL %in% "GO:0010446")

# AnnotationHub approach --------------------------------------------------
library(AnnotationHub)
library(AnnotationDbi)

ah <- AnnotationHub()
human_orgdb_record <- query(ah, c("OrgDb", "Homo sapiens"))

# inspect the resutls for the query and pick the version
human_orgdb_record

# Retrieve the database object using its ID
orgdb <- human_orgdb_record[["AH117067"]]

# explore the keytypes for the database
keytypes(orgdb)

test03_01 <- AnnotationDbi::select(orgdb,
                                   keys = unique(TOI$GO_id),
                                   columns = c('SYMBOL'),
                                   keytype = "GO")

head(test03_01)

test03_02 <- AnnotationDbi::select(orgdb,
                                   keys = unique(TOI$GO_id),
                                   columns = c('SYMBOL'),
                                   keytype = "GOALL")

head(test03_02)

test03_01_summary <- test03_01 %>%
  group_by(GO) %>%
  summarise(n = n())

test03_02_summary <- test03_02 %>%
  group_by(GOALL) %>%
  summarise(n = n())

# summarise the genes
TOI %>%
  left_join(test03_01_summary,by = c("GO_id" = "GO")) %>%
  mutate(test = !is.na(n)) %>%
  print(n = 30)

# sample pull the genes
test03_01 %>%
  filter(GO %in% "GO:0010446")

# summarise the genes
TOI %>%
  left_join(test03_02_summary,by = c("GO_id" = "GOALL")) %>%
  mutate(test = !is.na(n)) %>%
  print(n = 30)

# sample pull the genes
test03_02 %>%
  filter(GOALL %in% "GO:0010446")

identical(test03_01,test02_01)
identical(test03_02,test02_02)
