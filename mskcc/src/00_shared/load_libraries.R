# ============================================================================
# Shared helper: sourced by the analysis scripts, never run on its own.
# ============================================================================
##################
# Author: Ethan Tse
# Date: Sept 30th, 2024
#
# Load recurring libraries, helper functions, and customizations
#
##################

library(dplyr)
library(tibble)
library(readr)
library(tidyr)
library(readxl)

# for plotting
library(cowplot)
library(ggplot2)
library(ggrepel)
library(ggpubr)

# Load helper functions
source("./src/00_shared/utils.R")
source("./src/00_shared/styles.R")

# short cancer type names
# Make cancer codes for each cancer type
cancer_codes <- c(
  "NSLSC",
  "LUAD", # Lung Adenocarcinoma
  "PRAD", # Prostate Adenocarcinoma
  "HGSOC", # High-Grade Serous Ovarian Cancer
  "BRCA", # Breast Invasive Ductal Carcinoma
  "UCEC", # Uterine Endometrioid Carcinoma
  "COAD", # Colon Adenocarcinoma
  "STAD", # Stomach Adenocarcinoma
  "PDAC", # Pancreatic Adenocarcinoma
  "KIRC", # Renal Clear Cell Carcinoma
  "READ", # Rectal Adenocarcinoma
  "USC", # Uterine Serous Carcinoma/Uterine Papillary Serous Carcinoma
  "BLCA", # Bladder Urothelial Carcinoma
  "COADREAD", # Colorectal Adenocarcinoma
  "MCC", # Merkel Cell Carcinoma
  "SKCM", # Cutaneous Melanoma
  "ESCA", # Esophageal Adenocarcinoma
  "NBL", # Neuroblastoma
  "UCS", # Uterine Carcinosarcoma/Uterine Malignant Mixed Mullerian Tumor
  "GEJ", # Adenocarcinoma of the Gastroesophageal Junction
  "CLL", # Chronic Lymphocytic Leukemia/Small Lymphocytic Lymphoma
  "DLBCL", # Diffuse Large B-Cell Lymphoma, NOS
  "FL", # Follicular Lymphoma
  "GIST", # Gastrointestinal Stromal Tumor
  "LMS", # Leiomyosarcoma
  "Other"
)
names(cancer_codes) <- c(
  "Non-Small Cell Lung Cancer",
  "Lung Adenocarcinoma",
  "Prostate Adenocarcinoma",
  "High-Grade Serous Ovarian Cancer",
  "Breast Invasive Ductal Carcinoma",
  "Uterine Endometrioid Carcinoma",
  "Colon Adenocarcinoma",
  "Stomach Adenocarcinoma",
  "Pancreatic Adenocarcinoma",
  "Renal Clear Cell Carcinoma",
  "Rectal Adenocarcinoma",
  "Uterine Serous Carcinoma/Uterine Papillary Serous Carcinoma",
  "Bladder Urothelial Carcinoma",
  "Colorectal Adenocarcinoma",
  "Merkel Cell Carcinoma",
  "Cutaneous Melanoma",
  "Esophageal Adenocarcinoma",
  "Neuroblastoma",
  "Uterine Carcinosarcoma/Uterine Malignant Mixed Mullerian Tumor",
  "Adenocarcinoma of the Gastroesophageal Junction",
  "Chronic Lymphocytic Leukemia/Small Lymphocytic Lymphoma",
  "Diffuse Large B-Cell Lymphoma, NOS",
  "Follicular Lymphoma",
  "Gastrointestinal Stromal Tumor",
  "Leiomyosarcoma",
  "Other"
)
