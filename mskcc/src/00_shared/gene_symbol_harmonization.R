# ============================================================================
# Shared helper: sourced by the analysis scripts, never run on its own.
# ============================================================================
# Synonymous MSK-IMPACT gene symbols, collapsed to one name per gene before any
# mutation matrix is built, so a gene named differently across panel versions is
# not counted as two genes. Applied upstream of this bundle, at cohort-assembly
# time; the mapping is kept here because it is published as Supplementary
# Table 25 and is read out of this function rather than retyped.
#
# Direction is toward whichever symbol the pipeline collapses onto, which is not
# uniformly the newer one: MLL -> KMT2A moves to the IMPACT505 name, while
# AMER1 -> FAM123B moves to the older panels' name. This is not HUGO
# standardization.

rename_synonymous_genes <- function(df_with_hugo_symbol_col, column_name = "Hugo_Symbol") {
  impact_gene_mapping_vector <- c(
    # MLL in impact341, impact410, impact468 is KMT2A in impact505
    "MLL2" = "KMT2D",
    "MLL3" = "KMT2C",
    "MLL4" = "KMT2B",
    "MLL"  = "KMT2A"
    # unknown mapping for GTF2I
  )

  mskimpact_gene_mapping_vector <- c(
    # RAD51B in impact341, impact410, impact468 is RAD51L1 in impact505
    "RAD51B" = "RAD51L1",
    "RAD51D" = "RAD51L3",
    # AMER1 in impact505 is FAM123B in impact341, impact410, impact468
    "AMER1" = "FAM123B",
    "SETD8" = "KMT5A"
  )

  col <- column_name

  df_with_hugo_symbol_col[[col]] <- ifelse(
    df_with_hugo_symbol_col[[col]] %in% names(impact_gene_mapping_vector),
    impact_gene_mapping_vector[df_with_hugo_symbol_col[[col]]],
    df_with_hugo_symbol_col[[col]]
  )

  df_with_hugo_symbol_col[[col]] <- ifelse(
    df_with_hugo_symbol_col[[col]] %in% names(mskimpact_gene_mapping_vector),
    mskimpact_gene_mapping_vector[df_with_hugo_symbol_col[[col]]],
    df_with_hugo_symbol_col[[col]]
  )

  df_with_hugo_symbol_col
}
