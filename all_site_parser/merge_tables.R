setwd("git_repos/nucleotide_to_AA_parser/all_site_parser")
require("data.table")
require("tidyverse")
df <- fread("SARS-CoV-2_hookup_table.csv")
sign_df <- fread("T_testMinReplicates_2020_3_MinOffspring_4_embedded_1_MinTipsOfEachAllele_2.csv")
df %>%
  filter(protein_name == "ORF3a" & codon_number == 251)

merged_df <- sign_df %>%
  filter(signif_t_test == "signif") %>%
  rename(nucleotide_pos = position)

merged_df <- df %>%
  right_join(merged_df, by = "nucleotide_pos")

fwrite(merged_df, "T_testMinReplicates_2020_3_MinOffspring_4_embedded_1_MinTipsOfEachAllele_2.hooked_up.csv")
