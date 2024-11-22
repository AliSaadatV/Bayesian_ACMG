library(tidyverse)

### Read data
clinvar <- read_delim("data/raw/ClinVar/clinvar_20220730_P_LP_B_LB_withchr_CLNSIG.tsv.gz", na=".",
                      col_names=c("CHROM", "POS", "REF", "ALT", "ESP", "EXAC", "TGP", "CLINSTAT", "CLNSIG"))
clinvar_annot <- read_delim("/mnt/data3/saadat/database/ClinVar/clinvar_20220730_P_LP_B_LB_withchr_vep_annot.tsv", 
                            na=".", col_names = c("Consequence",
                                                  "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature",
                                                  "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp",
                                                  "cDNA_position", "CDS_position", "Protein_position",
                                                  "Amino_acids", "Codons",
                                                  "Existing_variation", "DISTANCE", "STRAND", "FLAGS",
                                                  "SYMBOL_SOURCE", "HGNC_ID", "CANONICAL", "MANE_SELECT",
                                                  "MANE_PLUS_CLINICAL", "ENSP", "DOMAINS", "HGVS_OFFSET",
                                                  "LoF", "LoF_filter", "LoF_flags", "LoF_info"))
clinvar_annot_coord <- read_delim("data/raw/clinvar_20220730_P_LP_B_LB_withchr_vep_coord.tsv.gz", 
                                  na=".", col_names = c("CHROM", "POS", "REF", "ALT", "Allele"))
clinvar_annot <- cbind(clinvar_annot_coord, clinvar_annot)
clinvar_annot <- clinvar_annot %>% 
  left_join(clinvar) %>% 
  filter((is.na(ESP) | ESP<=0.02) & (is.na(EXAC) | EXAC<=0.02))
rm(clinvar, clinvar_annot_coord)

### count missense patho per gene
clinvar_patho_total <- clinvar_annot %>% 
  filter(str_detect(CLNSIG, "athogenic")) %>% 
  rowwise() %>%
  group_by(Gene) %>% 
  summarise(count_patho_total=n())

clinvar_missense_total <- clinvar_annot %>% 
  filter((str_detect(Consequence, "missense") | str_detect(Consequence, "start_lost") | str_detect(Consequence, "stop_lost") & (str_length(REF)==1 & str_length(ALT)==1))) %>%  
  rowwise() %>%
  group_by(Gene) %>% 
  summarise(count_missense_total=n())

clinvar_benign_missense <- clinvar_annot %>% 
  filter(str_detect(CLNSIG, "enign")) %>% 
  filter((str_detect(Consequence, "missense") | str_detect(Consequence, "start_lost") | str_detect(Consequence, "stop_lost") & (str_length(REF)==1 & str_length(ALT)==1))) %>% 
  rowwise() %>%
  group_by(Gene) %>% 
  summarise(count_benign_missense=n())

clinvar_patho_missense_noinframe <- clinvar_annot %>% 
  filter(str_detect(CLNSIG, "athogenic")) %>% 
  filter((str_detect(Consequence, "missense") | str_detect(Consequence, "start_lost") | str_detect(Consequence, "stop_lost") & (str_length(REF)==1 & str_length(ALT)==1))) %>% 
  rowwise() %>%
  group_by(Gene) %>% 
  summarise(count_patho_missense=n()) %>% 
  left_join(clinvar_patho_total) %>% 
  left_join(clinvar_missense_total) %>% 
  left_join(clinvar_benign_missense) %>% 
  replace_na(list(count_patho_total=0, count_missense_total=0, count_benign_missense=0)) %>% 
  mutate(ratio_missensepatho_to_totalpatho=count_patho_missense/count_patho_total,
         ratio_missensebenign_to_totalmissense=count_benign_missense/count_missense_total,
         ratio_missensepatho_to_totalmissense=count_patho_missense/count_missense_total,
         ratio_missensepatho_to_missensebenign=count_patho_missense/count_benign_missense) %>% 
  filter(ratio_missensepatho_to_totalpatho>=0.7 & ratio_missensebenign_to_totalmissense<=0.3 & count_patho_missense>=5)
write_delim(select(clinvar_patho_missense_noinframe, Gene), "data/processed/PP2_genes.txt.gz")

### count truncating variants
clinvar_patho_truncating <- clinvar_annot %>% 
  filter(str_detect(CLNSIG, "athogenic")) %>% 
  filter(str_detect(Consequence, "frameshift_variant") | str_detect(Consequence, "splice_acceptor_variant") | str_detect(Consequence, "splice_donor_variant") | str_detect(Consequence, "stop_gained")) %>% 
  rowwise() %>%
  group_by(Gene) %>% 
  summarise(count_patho_truncating=n()) %>% 
  left_join(clinvar_patho_total) %>% 
  mutate(ratio_truncatingpatho_to_totalpatho=count_patho_truncating/count_patho_total) %>% 
  filter(ratio_truncatingpatho_to_totalpatho>=0.7 & count_patho_truncating>=5)
write_delim(select(clinvar_patho_truncating, Gene), "data/processedBP1_genes.txt.gz")
