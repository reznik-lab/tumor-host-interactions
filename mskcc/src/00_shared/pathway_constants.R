# ============================================================================
# Shared helper: sourced by the analysis scripts, never run on its own.
# ============================================================================
# ============================================================================
# Pathway and Cancer Type Constants
# ============================================================================
#
# Purpose:
#   Define gene lists for canonical cancer signaling pathways (based on
#   Sanchez-Vega et al., Cell, 2018) and obesity-associated cancer types
#   for use in mutation matrix annotation and cohort selection.
#
# Contents:
#   1. CANCER TYPE CONSTANTS
#      - obesity_cancers: 9 cancer types with established obesity associations
#      - full_obesity_associated_cancers: expanded list of histological subtypes
#   2. SIGNALING PATHWAY GENE LISTS (11 pathways)
#      - Cell cycle, Hippo, Myc, Notch, NRF2, PI3K/AKT/mTOR, RTK/RAS/MAPK,
#        TGFβ, p53, Wnt
#
# Note: These constants are sourced by establish_msk_impact_gene_pathways_df_and_mtx.R
#       to annotate binary mutation matrices with pathway-level alterations.
#
# ============================================================================

# ---- CANCER TYPE CONSTANTS ----

# Primary cancer types with strong obesity associations
obesity_cancers <- c(
  "BRCA", # Breast invasive carcinoma
  "UCEC", # Uterine corpus endometrial carcinoma
  "USC", # Uterine serous carcinoma
  "HGSOC", # High-grade serous ovarian cancer
  "COADREAD", # Colorectal adenocarcinoma
  "STAD", # Stomach adenocarcinoma
  "PDAC", # Pancreatic ductal adenocarcinoma
  "KIRC", # Kidney renal clear cell carcinoma
  "ESCA" # Esophageal carcinoma
)

# Expanded list of obesity-associated cancer histological subtypes
# Used for sample classification in body composition analysis
full_obesity_associated_cancers <- c(
  # uterine
  "Endometrial Carcinoma",
  "Uterine Endometrioid Carcinoma",
  "Uterine Mixed Endometrial Carcinoma",
  "Uterine Clear Cell Carcinoma",
  "Uterine Serous Carcinoma/Uterine Papillary Serous Carcinoma",
  "Uterine Carcinosarcoma/Uterine Malignant Mixed Mullerian Tumor",
  # breast
  "Breast Invasive Ductal Carcinoma",
  "Invasive Breast Carcinoma",
  "Breast Invasive Cancer, NOS",
  "Breast Invasive Lobular Carcinoma",
  "Breast Invasive Carcinoma, NOS",
  "Breast Mixed Ductal and Lobular Carcinoma",
  # colorectal
  "Colorectal Adenocarcinoma",
  "Rectal Adenocarcinoma",
  "Colon Adenocarcinoma",
  "Mucinous Adenocarcinoma of the Colon and Rectum",
  "High-Grade Neuroendocrine Carcinoma of the Colon and Rectum",
  # esophageal
  "Esophageal Adenocarcinoma",
  "Adenocarcinoma of the Gastroesophageal Junction",
  "Esophagogastric Adenocarcinoma",
  # kidney
  "Renal Clear Cell Carcinoma",
  # "Papillary Renal Cell Carcinoma",
  # "Chromophobe Renal Cell Carcinoma",
  "Renal Cell Carcinoma",
  "Unclassified Renal Cell Carcinoma",
  # pancreatic
  "Pancreatic Adenocarcinoma",
  "Adenosquamous Carcinoma of the Pancreas",
  "Acinar Cell Carcinoma of the Pancreas",
  # liver
  "Hepatocellular Carcinoma",
  # gallbladder
  "Gallbladder Cancer",
  "Gallbladder Adenocarcinoma, NOS",
  # Stomach cancer (cardia predominance)
  "Stomach Adenocarcinoma",
  "Poorly Differentiated Carcinoma of the Stomach",
  "Signet Ring Cell Carcinoma of the Stomach",
  # Thyroid cancer
  "Papillary Thyroid Cancer",
  "Hurthle Cell Thyroid Cancer",
  "Poorly Differentiated / Anaplastic Thyroid Cancer",
  # Meningioma
  "Meningioma",
  # Multiple myeloma
  "Plasma Cell Myeloma",
  # Ovarian
  "High-Grade Serous Ovarian Cancer",
  "Low-Grade Serous Ovarian Cancer"
)

# ---- SIGNALING PATHWAY GENE LISTS ----
# Reference: Sanchez-Vega et al., Cell, 2018 (doi:10.1016/j.cell.2018.03.022)
# Note: Gene names are standardized to lowercase HUGO symbols
# These lists are used to annotate binary mutation matrices with pathway-level calls:
#   1 = at least one gene altered, 0 = genes covered but unaltered, NA = genes not covered

# Cell cycle pathway: RB1 inactivation and cyclin/CDK activation
cell_cycle_pathway <- c(
  "rb1", "cdkn2a", "cdkn2b", "cdkn2c",
  "cdkn1a", "cdkn1b",
  "ccnd1", "ccnd2", "ccnd3", "ccne1",
  "cdk4", "cdk6",
  "e2f3"
)

# Hippo pathway: YAP/WWTR1 activation pathway
hippo_pathway <- c(
  "nf2",
  "lats1",
  "lats2",
  "yap1",
  "wwtr1", # taz
  "fat1"
  # "stk11"
)

# Myc pathway: Myc oncogene and associated regulatory factors
myc_pathway <- c(
  "myc", # c-myc, mycc
  "mycn", # n-myc
  "mycl", # l-myc, mycl1
  "max", # myc-associated factor x
  "mga" # max gene associated
)

# Notch pathway: Notch receptor signaling and Myc-related transcriptional regulation
notch_pathway <- c(
  "notch1",
  "notch2",
  "notch3",
  "notch4",
  "fbxw7",
  "dtx1",
  "crebbp",
  "ep300",
  "kdm5a",
  "spen"
)

# NRF2 pathway: Oxidative stress response through Nrf2/KEAP1 axis
nrf2_pathway <- c(
  "nfe2l2", # nrf2
  "keap1", # kelch-like ech-associated protein 1
  "cul3" # cullin 3
)

# PI3K/AKT/mTOR pathway: Receptor tyrosine kinase signaling and metabolic control
pi3k_akt_mtor_pathway <- c(
  "pik3ca",
  "pik3cb",
  "pik3cd",
  "pik3cg",
  "pik3r1",
  "pik3r2",
  "pik3r3",
  "pten",
  "akt1",
  "akt2",
  "akt3",
  "mtor",
  "rptor",
  "rictor",
  "tsc1",
  "tsc2",
  "rheb",
  "pdpk1",
  "inpp4b",
  "inppl1", # ship2
  "ppp2r1a",
  "stk11"
)

# RTK/RAS/MAPK pathway: Receptor tyrosine kinase to mitogenic signaling cascade
rtk_ras_mapk_pathway <- c(
  # rtks
  "egfr", "erbb2", "erbb3", "erbb4",
  "fgfr1", "fgfr2", "fgfr3", "fgfr4",
  "pdgfra", "pdgfrb",
  "met", "alk", "ret", "ros1",
  "ntrk1", "ntrk2", "ntrk3",
  "flt1", "flt3", "flt4",
  "kit", "igf1r", "insr",
  "axl", "ddr2", "tek", "kdr",
  "csf1r", "csf3r", "mpl",
  "epha3", "epha5", "epha7", "ephb1", "jak2",
  # ras/mapk
  "kras", "nras", "hras", "rit1", "rras", "rras2",
  "braf", "raf1", "araf",
  "map2k1",
  "mapk1",
  "sos1",
  "cbl",
  "nf1",
  "rasa1",
  "abl1",
  "rac1", "ptpn11",
  # unsure
  "map2k2",
  "map2k4",
  "pak1",
  "mapk3",
  "map3k1", "map3k13", "map3k14",
  "shoc2", "spred1",
  "src",
  "yes1",
  "lyn",
  "crkl",
  "gab1", "gab2", "pak7", "cdc42"
)

# TGFβ pathway: Transforming growth factor-beta signaling (EMT and fibrosis)
tgfb_pathway <- c(
  "tgfbr1",
  "tgfbr2",
  "acvr1", # alk2
  "bmpr1a", # alk3
  "smad2",
  "smad3",
  "smad4"
)

# p53 pathway: Tumor suppressor and DNA damage response
p53_pathway <- c(
  "tp53", # p53, tumor protein p53
  "tp53bp1", # 53bp1
  "mdm2",
  "mdm4",
  "atm", # ataxia telangiectasia mutated
  "chek2", # chk2
  "cdkn2a"
)

# Wnt pathway: Wnt/β-catenin signaling and stemness
wnt_pathway <- c(
  "apc", # adenomatous polyposis coli
  "ctnnb1", # β-catenin
  "axin1", # axis inhibition protein 1
  "axin2", # axis inhibition protein 2, conductin
  "rnf43", # ring finger protein 43
  "znrf3", # zinc and ring finger 3
  "fam123b", # amer1
  "tcf7l2", # tcf4
  "gsk3b" # glycogen synthase kinase 3 beta
)

# ============================================================================
# End of pathway_constants.R
# ============================================================================
