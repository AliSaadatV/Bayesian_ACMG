library(tidyverse)

### Read data
clinvar <- read_delim("data/raw/clinvar_20220730_P_LP_B_LB_withchr_CLNSIG.tsv.gz", na=".",
                      col_names=c("CHROM", "POS", "REF", "ALT", "ESP", "EXAC", "TGP", "CLINSTAT", "CLINSIG"))
clinvar_annot <- read_delim("data/raw/clinvar_20220730_P_LP_B_LB_withchr_vep_annot.tsv.gz", 
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

### clinvar
clinvar_annot <- clinvar_annot %>% 
  select(CHROM:ALT, Consequence, Gene, Feature, HGVSp, Protein_position, Amino_acids, CLINSTAT, CLINSIG)

write_delim(clinvar_annot, "data/processed/clinvar_PLP_BLB_CLNSIG_CLNSTAT.tsv.gz")

### uniprot 
unip <- read_delim("data/raw/humsvar_corrected.txt.gz", na = "-",
                   col_names = c("gene_name", "Swiss_Prot_AC", "FTID", "AAchange_uniprot", "Variant_category_uniprot", "dbSNP")) %>% 
  distinct(.keep_all = T)
gene_id <- gprofiler2::gconvert(unique(unip$gene_name)) %>% 
  select(input, target)
unip <- unip %>% 
  left_join(gene_id, by=c("gene_name"="input")) %>% 
  rename(Gene=target) %>% 
  distinct(.keep_all = T)

unip_na_genes <- unip %>% 
  filter(is.na(Gene)) #upload into https://www.uniprot.org/id-mapping

unip_na_genes_results <- read_delim("data/raw/uniprot_naIDs_convert.tsv.gz")
unip_na_genes_results <- unip_na_genes_results %>% 
  rename(Swiss_Prot_AC=From, Gene2=To) %>% 
  rowwise() %>% 
  mutate(Gene2=unlist(str_split(Gene2,"\\."))[1])

unip <- unip %>% 
  left_join(unip_na_genes_results) %>% 
  mutate(Gene=coalesce(Gene,Gene2)) %>% 
  select(-Gene2) %>% 
  filter(!is.na(Gene)) %>% 
  distinct()

write_delim(unip, "data/processed/humsvar_corrected_withID.txt.gz")
