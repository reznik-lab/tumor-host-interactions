# Shared definitions, sourced by every script.
#
# Produces no manuscript figure or table: the plotting theme, the BMI-usability
# rule and the minimum-cases-per-type floor that the reported analyses share.

# add conda environment to library path
.libPaths(c("/opt/conda/envs/packages/lib/R/library", .libPaths()))

## Shared library stack. Every script gets these by sourcing this file, so they
## should not be repeated in the scripts themselves. The tidyverse core
## packages are listed explicitly rather than left implicit in
## library(tidyverse), so scripts stay correct across tidyverse versions.
library(tidyverse)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(tibble)
library(stringr)
library(lubridate)
library(ggplot2)
library(cowplot)
library(tools)

theme_std <- function(base_size = 7) {
  theme_cowplot(font_size = base_size) %+replace%
    theme(
      line = element_line(colour = "black", linewidth = 0.1),
      axis.line.x = element_line(colour = "black", linewidth = 0.1),
      axis.line.y = element_line(colour = "black", linewidth = 0.1),
      axis.ticks = element_line(linewidth = 0.1),
      axis.ticks.length = unit(1, "pt")
    ) +
    theme(
      strip.background = element_blank()
      # strip.text = element_text(face = "bold")
    )
}


## ============================================================
## Which BMI measurements are usable -- ONE definition for the
## whole pipeline. Every script that reads measurementOccurrence
## must go through filter_usable_bmi() rather than writing its
## own filter, so the rows PSM matches on are the same rows the
## analyses run on.
##
## 1. NO CONCEPT FILTER. measurementOccurrence is already
##    BMI-only: t_disp_standard_concept_name has the single
##    value "Body mass index (BMI) [Ratio]". The
##    str_detect(..., "BMI") calls scattered downstream are
##    therefore no-ops, kept only as cheap assertions.
##
## 2. NO PROVENANCE FILTER. t_disp_measurement_type_concept_name
##    records where a measurement came from (EHR, From physical
##    examination, Patient self-report, No matching concept...).
##    Stage 01 used to drop "Patient self-report" (4,169),
##    "No matching concept" (273,651) and "No information" (3)
##    -- 277,823 of 3,257,958 rows, 8.5% -- while every later
##    stage kept them. So n_bmi_measurements, the PSM covariate,
##    counted a different set of measurements than the analyses
##    ever saw. All provenance types now count, everywhere.
##    ("Claim", also in that old exclusion list, never appears
##    in the data at all.)
##
## 3. VALUE RANGE IS THE ONLY GATE. [BMI_MIN, BMI_MAX] was
##    already applied in post_psm_filter_ctrls.R and
##    setup_me_modelling_df.R; it is now applied in Stage 01 and
##    Stage 05 too, which previously had no value gate. This
##    matters most for the weight-category cut points in
##    05_weight_class, whose case_when has an open-ended
##    `TRUE ~ "Obese"` tail: without the gate an implausibly
##    high BMI is silently classified rather than dropped.
## ============================================================
BMI_MIN <- 5
BMI_MAX <- 80

#' Keep only physiologically plausible BMI measurements.
#'
#' @param df  a measurement table
#' @param col name of the BMI value column. Defaults to the raw column name;
#'   setup_me_modelling_df.R renames it to "bmi" before filtering.
filter_usable_bmi <- function(df, col = "value_as_number") {
  df %>%
    filter(!is.na(.data[[col]]),
           .data[[col]] >= BMI_MIN,
           .data[[col]] <= BMI_MAX)
}


## ============================================================
## Minimum cancer-type size -- ONE definition for the whole
## pipeline. A cancer type is reported only if at least this
## many CASES of it are in the cohort.
##
## 1. COUNTED ON CASES ONLY, never on cases + controls pooled.
##    At 10:1 matching a pooled gate of 50 lets a type through
##    on ~5 cases (CODE_REVIEW_FINDINGS #17). Every caller
##    counts n_distinct(person_id[case == <case level>]); note
##    the case column is labelled "case"/"control" in Stage 04
##    and 0/1 in Stage 03, so the comparison differs by script
##    even though the threshold does not.
##
## 2. THIS IS A FLOOR, NOT THE ONLY GATE. Analyses that need
##    more data for their own reasons keep their extra, stricter
##    requirements on top of it -- the per-timebin minimum in
##    mixed_effects_modelling.R and the episode-count minimum
##    behind Figure 1B. Those gate different quantities (people
##    in a window, qualifying episodes) and are not case counts,
##    so they are NOT set from this constant.
## ============================================================
MIN_CASES_PER_TYPE <- 100


## ============================================================
## LANDMARK WINDOWS -- ONE definition for every script that
## fits a landmark model, for the same reason the BMI filter
## and the case-count floor live here. Two scripts now window
## the same measurements the same way; a second literal copy
## of this table is a drift waiting to happen.
##
## Windows are center +/- buffer, in DAYS from diagnosis
## (negative = before). In months, at 365.25/12 days/month:
##
##   5_years       [-1822.5, -1462.5]   48.05 - 59.88 before
##   3_years       [-1457.5, -1097.5]   36.06 - 47.89 before
##   1_year        [  -485,     -245]    8.05 - 15.93 before
##   6_months      [  -302.5,   -62.5]   2.05 -  9.94 before
##   at_diagnosis  [   -30,       30]    0.99 before to 0.99 after
##
## TWO THINGS THE KEYS DO NOT SAY:
##
## 1. "5_years" is the 48-60 month window and "3_years" is
##    36-48, not 60 and 36. That is deliberate -- the windows
##    are built to reproduce the trajectory panel's "48-60" and
##    "36-48" bins, and a true 5 years (-1826.25) would push
##    half the window past the days_from_dx >= -1825 cap in
##    setup_me_modelling_df.R. The manuscript facet already
##    reads "48-60 months"; the KEYS are what mislead, so do
##    not let "5 years" reach the Methods text.
##
## 2. "at_diagnosis" is SYMMETRIC, and the upstream filter
##    keeps measurements to +30 days, so this window really
##    does include up to a month AFTER diagnosis. It is
##    "within 30 days of diagnosis", not "at or before".
##
## The windows are unequal in width (7.9 months near diagnosis,
## 11.8 months out at 3-5 years, 2.0 at diagnosis) and they do
## not tile the interval -- nothing between 15.9 and 36.1
## months, nothing between 0.99 and 2.05. Each model sees only
## the measurements inside its own window.
## ============================================================
LANDMARK_TIMEPOINTS <- tibble::tibble(
  label  = c("6_months", "1_year", "3_years", "5_years", "at_diagnosis"),
  center = c(-182.5,     -365,     -1277.5,   -1642.5,   0),
  buffer = c(120,         120,      180,       180,      30)
) %>%
  dplyr::mutate(
    window_lo = center - buffer,
    window_hi = center + buffer
  )

## Display labels for the five, chronological order, furthest
## from diagnosis first -- the order a trajectory reads in.
LANDMARK_LABELS <- c(
  "5_years"      = "48-60 months",
  "3_years"      = "36-48 months",
  "1_year"       = "8-16 months",
  "6_months"     = "2-10 months",
  "at_diagnosis" = "At diagnosis"
)
LANDMARK_ORDER <- names(LANDMARK_LABELS)

stopifnot(setequal(LANDMARK_ORDER, LANDMARK_TIMEPOINTS$label))


## ============================================================
## Cancer type -> tissue-group mapping for pan-cancer analysis.
##
## assign_tissue_group(raw_label) takes a vector of raw
## cancer_type / condition labels and returns a vector of
## tissue-group labels, suitable for group_by() upstream of a
## per-group PSM matching run.
##
## Groups are split so that, within this ordering, subtype-
## specific hematologic branches and site-specific solid-tumor
## branches are checked before their broader/unspecified
## catch-alls -- see inline notes on ordering-sensitive branches.
## ============================================================

normalize_dx <- function(x) {
  x %>%
    str_to_lower() %>%
    str_replace_all("\\bprimary\\b", "") %>%
    str_replace_all("\\brecurrent\\b", "") %>%
    str_replace_all("\\b(left|right|bilateral)\\b", "") %>%
    str_squish()
}

assign_tissue_group <- function(raw_label) {
  dx <- normalize_dx(raw_label)
  ## The multi-primary guard has to see the label BEFORE normalize_dx() runs.
  ## normalize_dx() strips \bprimary\b -- correct for "Primary malignant
  ## neoplasm of stomach", where "primary" is noise -- but it also turns
  ## "multiple primary malignant neoplasms" into "multiple malignant
  ## neoplasms", so the "multiple primary" alternative below could never match.
  ##
  ## FIXED: that made this whole branch a no-op for multi-primary labels. An
  ## unqualified one still reached "Other" via the TRUE ~ fallback, by luck. A
  ## SITE-QUALIFIED one did not: "Multiple primary malignant neoplasms of
  ## colon" fell through to the colorectal branch and was filed as a genuine
  ## colorectal primary. Since cancer_type_grouped is what psm.R matches on and
  ## every per-type figure stratifies by, such a patient was counted as a
  ## single-site case carrying a diagnosis date and a "pre-diagnosis"
  ## trajectory that belong to a multi-primary history.
  raw <- str_to_lower(raw_label)
  case_when(
    str_detect(raw, "multiple primary") |
      str_detect(dx, "hereditary nonpolyposis|hnpcc|independent multiple sites")
      ~ "Other / multiple primary / unclear",
    ## --------------------------------------------------------
    ## Hematologic malignancies -- split into rough biological
    ## subgroups (Hodgkin lymphoma; non-Hodgkin lymphoma by
    ## B-cell/T-cell-NK/unspecified; leukemia by acute/chronic/
    ## other; MDS-MPN; multiple myeloma & plasma cell neoplasms)
    ## rather than one lumped "Hematologic" bucket. This whole
    ## section must run BEFORE the organ/site-specific branches
    ## below, since heme labels are frequently qualified by an
    ## anatomic site (e.g. "lymphoma of testis", "leukemic
    ## reticuloendotheliosis of lymph nodes of inguinal region
    ## and lower limb") that would otherwise be caught by those
    ## site-based regexes first. Within this section, subtype-
    ## specific branches (Hodgkin, T/NK-cell NHL, B-cell NHL,
    ## hairy cell leukemia via its "leukemic reticuloendotheliosis"
    ## synonym) run before the catch-all "unspecified" branches
    ## for lymphoma and leukemia, so subtype is preserved wherever
    ## the label states it.
    ## --------------------------------------------------------
    ## Hodgkin lymphoma (any site/subtype)
    str_detect(dx, "hodgkin")
      ~ "Lymphoma - Hodgkin",
    ## Non-Hodgkin lymphoma, T-cell / NK-cell subtypes.
    ## Must run before the B-cell and unspecified NHL branches.
    str_detect(dx, paste0(
      "t-cell lymphoma|t-lymphoblastic lymphoma|nk-cell lymphoma|",
      "mycosis fungoides|s(e|\u00e9)zary|",
      "angioimmunoblastic|hepatosplenic t-cell|enteropathy-associated|",
      "adult t-cell leukemia/lymphoma|peripheral t-cell lymphoma|",
      "cutaneous t-cell lymphoma"
    ))
      ~ "Lymphoma - non-Hodgkin (T-cell / NK-cell)",
    ## Non-Hodgkin lymphoma, B-cell subtypes.
    str_detect(dx, paste0(
      "b-cell lymphoma|b-lymphoblastic|diffuse large b-cell|",
      "follicular lymphoma|marginal zone lymphoma|",
      "mucosa-associated lymphoid tissue|malt lymphoma|",
      "burkitt|mantle cell lymphoma|small lymphocytic|",
      "lymphoplasmacytic|plasmablastic lymphoma|",
      "cutaneous b-cell lymphoma|precursor cell lymphoblastic lymphoma|",
      "diffuse follicle center lymphoma|high grade b-cell lymphoma|",
      "low grade b-cell lymphoma|nodular lymphoma"
    ))
      ~ "Lymphoma - non-Hodgkin (B-cell)",
    ## "Leukemic reticuloendotheliosis" is the historical name for
    ## hairy cell leukemia, not a lymphoma -- must be caught here,
    ## before the generic "reticuloendotheliosis" match below routes
    ## it to the NHL-unspecified bucket instead.
    str_detect(dx, "leukemic reticuloendotheliosis")
      ~ "Leukemia - chronic",
    ## Non-Hodgkin lymphoma / lymphoid neoplasm, no B/T/NK subtype
    ## specified in the label.
    str_detect(dx, paste0(
      "lymphoma|lymphosarcoma|reticulosarcoma|reticuloendotheliosis|",
      "non-hodgkin"
    ))
      ~ "Lymphoma - non-Hodgkin (unspecified/other)",
    ## Acute leukemias (AML and ALL, all FAB/molecular subtypes).
    ## Word-boundary anchors on "acute" are required so this does
    ## NOT false-match inside "subacute leukemia" (which has no
    ## space/boundary before "acute" and belongs in the
    ## other/unspecified leukemia bucket instead).
    str_detect(dx, paste0(
      "\\bacute myel|\\bacute lymphoid leukemia|\\bacute lymphoblastic leukemia|",
      "\\bacute promyelocytic|\\bacute leukemia|",
      "\\bacute myelomonocytic|b-cell acute lymphoblastic leukemia|",
      "t-lymphoblastic leukemia"
    ))
      ~ "Leukemia - acute",
    ## Chronic leukemias, including hairy cell leukemia.
    str_detect(dx, paste0(
      "chronic lymphoid leukemia|chronic lymphocytic leukemia|",
      "chronic myeloid leukemia|chronic leukemia|",
      "hairy cell leukemia|b-cell chronic lymphocytic leukemia"
    ))
      ~ "Leukemia - chronic",
    ## Leukemia, subtype/chronicity not otherwise specified above.
    str_detect(dx, paste0(
      "leukemi|leukaemi|myeloid leukemia|lymphoid leukemia|",
      "megakaryocytic leukemia|mast cell leukemia|",
      "large granular lymphocytic leukemia|subacute leukemia|",
      "subacute lymphoid leukemia"
    ))
      ~ "Leukemia - other/unspecified",
    ## Myelodysplastic / myeloproliferative neoplasms.
    str_detect(dx, "myelodysplastic|myeloproliferative|polycythemia vera")
      ~ "Myelodysplastic / myeloproliferative neoplasm (MDS/MPN)",
    ## Multiple myeloma and other plasma cell neoplasms.
    str_detect(dx, "myeloma|plasmacytoma|plasma cell leukemia")
      ~ "Multiple myeloma & plasma cell neoplasms",
    str_detect(dx, "melanoma")
      ~ "Skin - melanoma",
    str_detect(dx, "basal cell carcinoma|merkel cell|squamous cell carcinoma of skin|neoplasm of skin\\b|sebaceous")
      ~ "Skin - non-melanoma",
    str_detect(dx, "breast")
      ~ "Breast",
    str_detect(dx, "prostate")
      ~ "Prostate",
    ## --------------------------------------------------------
    ## Pancreatic neuroendocrine (incl. islets of Langerhans) --
    ## split out from exocrine pancreas. Must run BEFORE the
    ## generic pancreas branch below, since these labels also
    ## contain "pancrea" and would otherwise be caught there.
    ## --------------------------------------------------------
    str_detect(dx, "islets of langerhans|neuroendocrine (tumor|neoplasm) of pancreas|pancreatic neuroendocrine")
      ~ "Pancreatic neuroendocrine (incl. islets of Langerhans)",
    ## --------------------------------------------------------
    ## Ampulla of Vater -- split out from exocrine pancreas.
    ## --------------------------------------------------------
    str_detect(dx, "ampulla of vater")
      ~ "Ampulla of Vater",
    ## --------------------------------------------------------
    ## Anal canal -- split out from colorectal. Includes
    ## cloacogenic zone (anal transitional zone), embryologically
    ## part of the anal canal rather than the colon/rectum.
    ## Must run BEFORE the colorectal branch below.
    ## --------------------------------------------------------
    str_detect(dx, "anal canal|\\banus\\b|cloacogenic zone|perianal")
      ~ "Anal canal",
    str_detect(dx, "colon|rectosigmoid|rectum|cecum|appendix|hepatic flexure|splenic flexure|large intestine")
      ~ "Colorectal",
    str_detect(dx, "pancrea")
      ~ "Pancreas",
    str_detect(dx, "liver|gallbladder|bile duct|biliary tract|cholangiocarcinoma")
      ~ "Liver & biliary tract",
    str_detect(dx, "stomach|cardia|pyloric antrum|pylorus")
      ~ "Stomach",
    str_detect(dx, "esophag")
      ~ "Esophagus",
    str_detect(dx, "duodenum|jejunum|ileum|small intestine")
      ~ "Small intestine",
    str_detect(dx, "kidney|renal pelvis")
      ~ "Kidney",
    str_detect(dx, "urinary bladder|ureter|urethra|urachus")
      ~ "Bladder & lower urinary tract",
    str_detect(dx, "lung|bronchus|trachea|respiratory tract")
      ~ "Lung",
    ## --------------------------------------------------------
    ## Uterine / endometrial -- split out from the rest of the
    ## female genital tract (cervix, ovary, fallopian tube,
    ## vulva, vagina, etc.). Must run BEFORE the female genital
    ## tract branch below, since "uterine adnexa" refers to the
    ## ovaries/fallopian tubes (not the uterus itself) and is
    ## deliberately excluded here via negative lookahead so it
    ## falls through to the female genital tract group instead.
    ## --------------------------------------------------------
    str_detect(dx, "uterus|uterine(?! adnexa)|endometri|myometrium|corpus uteri|isthmus of uterus")
      ~ "Uterine & endometrial",
    str_detect(dx, "cervix|cervical(?! esophagus)|ovar|fallopian|vulva|vagina|broad ligament|clitoris|labia|parametrium|uterine adnexa|female genital")
      ~ "Female genital tract",
    str_detect(dx, "testis|testicular|penis|epididymis|spermatic cord|scrotum|male genital")
      ~ "Male genital tract",
    str_detect(dx, "thyroid")
      ~ "Thyroid",
    str_detect(dx, "adrenal|parathyroid|pituitary|pineal gland|paraganglion|carotid body|endocrine gland")
      ~ "Other endocrine",
    str_detect(dx, "tongue|tonsil|pharyn|laryng|glottis|epiglott|oral cavity|salivary|parotid|sublingual gland|submaxillary|nasal|sinus|lip\\b|palate|gum\\b|floor of mouth|buccal|mouth|uvula|vallecula|waldeyer|postcricoid|retromolar|mandible|\\bear\\b")
      ~ "Head & neck",
    str_detect(dx, "brain|cerebell|cerebr|meninges|spinal cord|cranial nerve|peripheral nerve|nervous system|frontal lobe|parietal lobe|occipital lobe|temporal lobe|craniopharyngeal|cauda equina")
      ~ "Central nervous system",
    str_detect(dx, "eye\\b|eyelid|conjunctiva|cornea|choroid|ciliary body|retina|orbit|lacrimal|optic nerve")
      ~ "Eye",
    str_detect(dx, "\\bbone\\b|skull|vertebral column|\\brib\\b|sternum|clavicle|long bone|short bone|articular cartilage|pelvic bones|scapula|limb bones")
      ~ "Bone",
    str_detect(dx, "soft tissue")
      ~ "Soft tissue",
    str_detect(dx, "thymus")
      ~ "Thymus",
    str_detect(dx, "heart")
      ~ "Heart",
    str_detect(dx, "peritoneum|retroperitoneum|abdomen|abdominal esophagus|intra-abdominal")
      ~ "Peritoneum / abdomen (unspecified)",
    str_detect(dx, "mediastinum|pleura")
      ~ "Mediastinum & pleura",
    str_detect(dx, "spleen")
      ~ "Spleen",
    str_detect(dx, "placenta")
      ~ "Placenta",
    str_detect(dx, "unknown site")
      ~ "Other / multiple primary / unclear",
    TRUE ~ "Other / multiple primary / unclear"
  )
}

## ============================================================
## P-VALUE TYPOGRAPHY -- ONE definition for every figure in the
## paper, for the same reason filter_usable_bmi() is shared: a
## convention that has to be identical across panels cannot be
## re-typed per script and stay identical.
##
## House style:
##   * capital P, ITALIC          -- italic(P)
##   * BH-adjusted values get a subscript  -- italic(P)[adj]
##   * the exact value, never a "< 0.001" threshold; below 0.001
##     it renders as m x 10^e rather than rounding to 0.000, so
##     3.10 x 10^-12 and 9.20 x 10^-4 stay distinguishable
##
## ggplot cannot italicise part of a string, so fmt_p() returns
## PLOTMATH source, not display text. EVERY layer that draws it
## must pass parse = TRUE:
##
##   geom_text(aes(label = label), parse = TRUE)
##
## and anywhere a legend title or axis label needs the symbol,
## use the expression directly:
##
##   scale_colour_manual(name = expression(italic(P)[adj] < 0.05))
##
## @param p         numeric vector of p-values. NA yields an empty
##                  label (or just the prefix), never "NA".
## @param adjusted  TRUE for BH-adjusted values -> italic(P)[adj].
## @param digits    decimals for values >= 0.001.
## @param prefix    optional plain text placed before the symbol,
##                  e.g. "LME estimate = 1.55,". Quoted for you;
##                  must not itself contain a single quote.
## ============================================================
fmt_p <- function(p, adjusted = FALSE, digits = 3, prefix = NULL) {
  sym <- if (adjusted) "italic(P)[adj]" else "italic(P)"
  pre <- if (is.null(prefix)) rep("", length(p)) else sprintf("'%s'", prefix)
  vapply(seq_along(p), function(i) {
    x <- p[i]
    if (is.na(x)) return(if (nzchar(pre[i])) pre[i] else "''")
    num <- if (x < 0.001) {
      e <- floor(log10(x))
      sprintf("%.2f %%*%% 10^%d", x / 10^e, e)
    } else {
      sprintf("'%s'", formatC(x, digits = digits, format = "f"))
    }
    body <- paste0(sym, "~'='~", num)
    if (nzchar(pre[i])) paste0(pre[i], "~", body) else body
  }, character(1))
}


## ============================================================
## CANCER-TYPE DISPLAY LABELS -- ONE mapping for every figure,
## for the same reason fmt_p() and filter_usable_bmi() live
## here: a convention that has to be identical across panels
## cannot be re-typed per script and stay identical. Three
## scripts previously carried their own short_type_label().
##
## PRESENTATION ONLY. The strings assign_tissue_group() returns
## are KEYS: psm.R matches on them, they are written into the
## stable_data TSVs, and every downstream join is on the raw
## value. Renaming them at the source would silently break
## joins against any file written before the rename. So nothing
## is ever STORED under a display label -- call this inside the
## ggplot() pipeline, never upstream of a join or a write_tsv().
##
## House style for the labels themselves:
##   * the distinguishing term first, not the genus
##     ("B-cell NHL", not "Lymphoma - non-Hodgkin (B-cell)"):
##     rows are ordered by prevalence, not alphabetically, so a
##     shared prefix buys nothing and costs axis width
##   * standard abbreviations where a clinician uses them
##     (NHL, MDS / MPN, NET, CNS)
##   * "(other)" as the one catch-all suffix
##   * "&" only for a genuine anatomic union
##
## The mapping MUST stay injective: the figures build factor
## levels from the mapped names (type_order_all in
## cachexia_frequency.R and its sensitivity twin), so two raw
## labels collapsing onto one display string would corrupt row
## order on every panel sharing that ordering. Asserted below at
## source() time rather than left to review.
##
## Types not listed here pass through unchanged -- Breast,
## Prostate, Lung, Colorectal, Stomach, Esophagus, Kidney,
## Thyroid, Head & neck, Liver & biliary tract, Small intestine,
## Anal canal, Ampulla of Vater, Eye, Bone, Soft tissue,
## Mediastinum & pleura, Thymus, Heart, Spleen, Placenta, and
## Pancreas. "Pancreas" is deliberately NOT relabelled
## "Pancreas (exocrine)" even though the neuroendocrine tumours
## are split out of it by assign_tissue_group().
## ============================================================
TYPE_DISPLAY_LABELS <- c(
  ## hematologic
  "Lymphoma - Hodgkin"                                     = "Hodgkin lymphoma",
  "Lymphoma - non-Hodgkin (B-cell)"                        = "B-cell NHL",
  "Lymphoma - non-Hodgkin (T-cell / NK-cell)"              = "T/NK-cell NHL",
  "Lymphoma - non-Hodgkin (unspecified/other)"             = "NHL (other)",
  "Leukemia - acute"                                       = "Acute leukemia",
  "Leukemia - chronic"                                     = "Chronic leukemia",
  "Leukemia - other/unspecified"                           = "Leukemia (other)",
  "Myelodysplastic / myeloproliferative neoplasm (MDS/MPN)" = "MDS / MPN",
  "Multiple myeloma & plasma cell neoplasms"               = "Multiple myeloma",
  ## solid tumours
  "Skin - melanoma"                                        = "Melanoma",
  "Skin - non-melanoma"                                    = "Skin (non-melanoma)",
  "Pancreatic neuroendocrine (incl. islets of Langerhans)" = "Pancreatic NET",
  "Bladder & lower urinary tract"                          = "Bladder & urinary tract",
  "Uterine & endometrial"                                  = "Uterus & endometrium",
  "Female genital tract"                                   = "Gynecologic (other)",
  "Male genital tract"                                     = "Testis & male genital",
  "Central nervous system"                                 = "Brain & CNS",
  "Other endocrine"                                        = "Endocrine (other)",
  "Peritoneum / abdomen (unspecified)"                     = "Peritoneum & abdomen",
  ## catch-all. NOT shortened to bare "Other": this bucket is
  ## where site-qualified multi-primary cases land, and after the
  ## assign_tissue_group() multi-primary fix it holds more of
  ## them than it used to. The axis should say so.
  "Other / multiple primary / unclear"                     = "Other / multiple primary"
)

stopifnot(
  ## no two raw types share a display label
  !anyDuplicated(unname(TYPE_DISPLAY_LABELS)),
  ## and no display label is itself a raw type, which would make
  ## the mapping read as a chain
  !any(unname(TYPE_DISPLAY_LABELS) %in% names(TYPE_DISPLAY_LABELS))
)

#' Raw cancer_type_grouped values -> the labels that go on a figure.
#'
#' @param x character/factor vector of raw group labels. Unmapped values are
#'   returned unchanged, so a new group added to assign_tissue_group() shows up
#'   under its raw name rather than as NA.
display_type_label <- function(x) {
  x   <- as.character(x)
  hit <- match(x, names(TYPE_DISPLAY_LABELS))
  unname(ifelse(is.na(hit), x, TYPE_DISPLAY_LABELS[hit]))
}


ggsave_wrapper <- function(filename, plot, extension = "pdf", width = 7, height = 7) {
  # Wrapper for ggsave that saves to the correct directory
  # can be used to force saved images to be of a specific extension
  # extracting extension from the file: https://stackoverflow.com/questions/7779037/extract-file-extension-from-file-path
  require(tools)

  if (!is.null(extension)) {
    # remove extension in filename
    filepath <- tools::file_path_sans_ext(filename)
    filepath <- paste0(filepath, ".", extension)
  }

  if (extension == "png") {
    # need to use bg = "white" for png because black and white
    # are reversed in the background
    ggsave(
      filename = filepath,
      plot = plot,
      device = extension,
      width = width,
      height = height,
      create.dir = TRUE,
      bg = "white",
      dpi = 600
    )
  } else if (extension == "svg") {
    ggsave(
      filename = filepath,
      plot = plot + theme(text = element_text(family = "sans")),
      device = svglite::svglite,
      width = width,
      height = height,
      create.dir = TRUE,
      bg = "white",
      dpi = 600
    )
  } else {
    ggsave(
      filename = filepath,
      plot = plot,
      device = extension,
      width = width,
      height = height,
      create.dir = TRUE
    )
  }

  # save rds of the ggplot so we can edit later if needed
  saveRDS(plot, file = paste0(tools::file_path_sans_ext(filepath), ".rds"))
}


## ============================================================
## Table output -- ONE writer, two formats.
##
## Every results table is written twice, as .tsv and as .xlsx,
## from a single call. The TSV stays the machine-readable copy
## the pipeline and git diffs read; the XLSX is the copy that
## goes into a supplement and gets opened by a human.
##
## Call it with the path you would have given write_tsv() --
## the extension is stripped and both are written beside each
## other, so a caller never names the same table twice and the
## two copies cannot drift apart or go out of sync.
##
## openxlsx rather than writexl: 02_matching/psm_quality_control.R
## already depends on it for the per-type balance workbooks, so
## this adds no new package to the Workbench image. It is loaded
## lazily -- a script that writes no tables does not need it.
##
## NOT for the intermediate files under stable_data/ and temp/.
## Those are pipeline plumbing, some of them one row per BMI
## measurement, and the XLSX format caps at 1,048,576 rows. The
## guard below refuses rather than truncating, because a silently
## short spreadsheet is worse than a missing one.
## ============================================================

## Excel's hard row cap, less one for the header row.
XLSX_MAX_ROWS <- 1048575L

#' Write a results table as both .tsv and .xlsx.
#'
#' @param x data frame to write.
#' @param path destination path. Any extension is stripped; ".tsv" and ".xlsx"
#'   are appended. Parent directories are created.
#' @param sheet worksheet name. Defaults to the file's base name. Excel caps
#'   sheet names at 31 characters and forbids []:*?/\, so the name is sanitised
#'   and truncated -- pass `sheet` explicitly if the truncation would be
#'   ambiguous between two tables in the same directory.
#' @return the base path, invisibly.
write_table <- function(x, path, sheet = NULL) {
  stopifnot(is.data.frame(x), length(path) == 1L, is.character(path))

  base <- tools::file_path_sans_ext(path)
  dir.create(dirname(base), recursive = TRUE, showWarnings = FALSE)

  ## TSV first and unconditionally: if the XLSX step fails -- a missing
  ## package, a list-column openxlsx cannot serialise -- the authoritative
  ## copy is already on disk rather than lost with the error.
  readr::write_tsv(x, paste0(base, ".tsv"))

  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("write_table() needs the openxlsx package to write ", basename(base),
         ".xlsx; the .tsv was written. install.packages(\"openxlsx\")")
  }

  if (nrow(x) > XLSX_MAX_ROWS) {
    warning(sprintf(
      "%s has %d rows, over Excel's %d-row cap; wrote the .tsv only. This table is too large to be a supplementary table -- check it is not an intermediate file.",
      basename(base), nrow(x), XLSX_MAX_ROWS), call. = FALSE)
    return(invisible(base))
  }

  if (is.null(sheet)) sheet <- basename(base)
  ## Excel forbids \ / ? * [ ] : in a sheet name and caps it at 31 characters.
  ## The ] must come FIRST inside the bracket expression -- anywhere else it
  ## closes the class early and the substitution silently does nothing.
  sheet <- substr(gsub("[]\\[*?:/]", "_", sheet), 1L, 31L)
  if (!nzchar(sheet)) sheet <- "Sheet1"

  ## The named-list form is what sets the sheet name in openxlsx; passing
  ## sheetName= to write.xlsx() is silently ignored for a bare data frame.
  openxlsx::write.xlsx(stats::setNames(list(x), sheet),
                       file = paste0(base, ".xlsx"),
                       overwrite = TRUE)

  invisible(base)
}
