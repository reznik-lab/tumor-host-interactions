# ============================================================================
# Shared helper: sourced by the analysis scripts, never run on its own.
# ============================================================================
# Gene set behind every PI3K/AKT/mTOR pathway-alteration call in this bundle.
# Symbols are lower-case: the mutation tables are matched on
# tolower(Hugo_Symbol). The same list, alongside the other oncogenic pathways,
# is published as Supplementary Table 21 and lives in pathway_constants.R --
# the two definitions are identical and were checked against each other.

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
