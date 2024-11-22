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

### Read UCSC files (these files are 0-based. add 1 to Start coordinated to make them 1-based)
uniprot_zf <- read_delim("data/raw/UCSC_uniprotdomains.txt.gz", delim="\t",
                              skip=1, col_types = "cddcdcddcdccccccccccccccccccc") %>% 
  filter(unipDomain.annotationType=="zinc finger region") %>% 
  mutate(unipDomain.chromStart=unipDomain.chromStart+1) 

colnames(uniprot_zf) <-  unlist(str_split(colnames(uniprot_zf), "\\."))[seq(2, 2*ncol(uniprot_zf), by=2)]

uniprot_other <- read_delim("data/raw/UCSC_uniprotother.txt.gz", skip=1, delim="\t",
                            col_types = "cddcdcddcdccccccccccccccccccc") %>% 
  filter(unipOther.annotationType %in% c("active site", "binding site", "calcium-binding region", "DNA-binding region",
                                         "metal ion-binding site", "nucleotide phosphate-binding region", "site")) %>% 
  mutate(unipOther.chromStart=unipOther.chromStart+1)

colnames(uniprot_other) <-  unlist(str_split(colnames(uniprot_other), "\\."))[seq(2, 2*ncol(uniprot_other), by=2)]

### Break rows with multiple "blocks" into single blocks
uniprot <- rbind(uniprot_zf, uniprot_other) %>% 
  filter(chrom %in% c(str_c("chr", 1:22), "chrX", "chrY", "chrM"))

uniprot_corrected <- uniprot %>% 
  filter(blockCount==1) %>% 
  select(chrom, chromStart, chromEnd)

uniprot <- uniprot %>% 
  filter(blockCount>1)

for(i in 1:nrow(uniprot)){
  starts <- uniprot$chromStart[i] + as.numeric(unlist(str_split(uniprot$chromStarts[i], ",")))  
  ends <- starts + as.numeric(unlist(str_split(uniprot$blockSizes[i], ",")))
  temp <- data.frame(chrom=rep(uniprot$chrom[i], uniprot$blockCount[i]),
                     chromStart=starts,
                     chromEnd=ends)
  uniprot_corrected <- rbind(uniprot_corrected, temp)
}
uniprot_corrected <- uniprot_corrected %>% 
  rename(CHROM=chrom, Start=chromStart, End=chromEnd)

### find domain w/o benign variants (and those w/o benign and >=1 pathogenic)
clinvar_annot_benign <- clinvar_annot %>% 
  filter(str_detect(CLINSIG, "enign") & !str_detect(CLINSIG, "onflict") & !str_detect(CLINSIG, "ncertain") & !str_detect(CLINSIG, "risk_factor") & !str_detect(CLINSIG, "ffects") & !str_detect(CLINSIG, "drug_response")) %>% 
  filter((str_detect(Consequence, "missense") | str_detect(Consequence, "start_lost") | str_detect(Consequence, "stop_lost") & (str_length(REF)==1 & str_length(ALT)==1))) %>% 
  rowwise() %>%
  mutate(END=POS+str_length(ALT))
clinvar_annot_patho <- clinvar_annot %>% 
  filter(str_detect(CLINSIG, "athogenic")) %>% 
  filter((str_detect(Consequence, "missense") | str_detect(Consequence, "start_lost") | str_detect(Consequence, "stop_lost") & (str_length(REF)==1 & str_length(ALT)==1))) %>% 
  rowwise() %>%
  mutate(END=POS+str_length(ALT))
uniprot_domains_ranges <- GenomicRanges::GRanges(seqnames = uniprot_corrected$CHROM, ranges = IRanges::IRanges(start = uniprot_corrected$Start, end = uniprot_corrected$End))
clinvar_benign_ranges <- GenomicRanges::GRanges(seqnames = clinvar_annot_benign$CHROM, ranges = IRanges::IRanges(start = clinvar_annot_benign$POS, end = clinvar_annot_benign$END))
clinvar_patho_ranges <- GenomicRanges::GRanges(seqnames = clinvar_annot_patho$CHROM, ranges = IRanges::IRanges(start = clinvar_annot_patho$POS, end = clinvar_annot_patho$END))
uniprot_clinvar_overlap_noB <- IRanges::subsetByOverlaps(uniprot_domains_ranges, clinvar_benign_ranges, invert = T)
uniprot_clinvar_overlap_noB_1P <- IRanges::subsetByOverlaps(uniprot_clinvar_overlap_noB, clinvar_patho_ranges)
uniprot_clinvar_overlap_noB_1P <- GenomicRanges::reduce(uniprot_clinvar_overlap_noB_1P)
uniprot_clinvar_overlap_noB <- IRanges::subsetByOverlaps(uniprot_clinvar_overlap_noB, uniprot_clinvar_overlap_noB_1P, invert = T)
uniprot_clinvar_overlap_noB <- GenomicRanges::reduce(uniprot_clinvar_overlap_noB)
uniprot_clinvar_overlap_noB_1P <- data.frame(Chr=GenomicRanges::seqnames(uniprot_clinvar_overlap_noB_1P),
                                             Start=GenomicRanges::start(uniprot_clinvar_overlap_noB_1P),
                                             End=GenomicRanges::end(uniprot_clinvar_overlap_noB_1P))
uniprot_clinvar_overlap_noB <- data.frame(Chr=GenomicRanges::seqnames(uniprot_clinvar_overlap_noB),
                                          Start=GenomicRanges::start(uniprot_clinvar_overlap_noB),
                                          End=GenomicRanges::end(uniprot_clinvar_overlap_noB))

write_delim(uniprot_clinvar_overlap_noB, "data/processed/uniprot_functional_regions_no_benign.tsv.gz")
write_delim(uniprot_clinvar_overlap_noB_1P, "data/processed/uniprot_functional_regions_no_benign_1_patho.tsv.gz")


######## uniprot domains
uniprot_domain <- read_delim("data/raw/UCSC_uniprotdomains.txt.gz", delim="\t",
                         skip=1, col_types = "cddcdcddcdccccccccccccccccccc") %>% 
  filter(unipDomain.annotationType=="domain") %>% 
  mutate(unipDomain.chromStart=unipDomain.chromStart+1) 

colnames(uniprot_domain) <-  unlist(str_split(colnames(uniprot_domain), "\\."))[seq(2, 2*ncol(uniprot_domain), by=2)]

uniprot <- uniprot_domain %>% 
  filter(chrom %in% c(str_c("chr", 1:22), "chrX", "chrY", "chrM"))

uniprot_corrected <- uniprot %>% 
  filter(blockCount==1) %>% 
  select(chrom, chromStart, chromEnd)

uniprot <- uniprot %>% 
  filter(blockCount>1)

for(i in 1:nrow(uniprot)){
  starts <- uniprot$chromStart[i] + as.numeric(unlist(str_split(uniprot$chromStarts[i], ",")))  
  ends <- starts + as.numeric(unlist(str_split(uniprot$blockSizes[i], ",")))
  temp <- data.frame(chrom=rep(uniprot$chrom[i], uniprot$blockCount[i]),
                     chromStart=starts,
                     chromEnd=ends)
  uniprot_corrected <- rbind(uniprot_corrected, temp)
}
uniprot_corrected <- uniprot_corrected %>% 
  rename(CHROM=chrom, Start=chromStart, End=chromEnd)

### find domain w/o benign variants (and those w/o benign and >=1 pathogenic)
clinvar_annot_benign <- clinvar_annot %>% 
  filter(str_detect(CLINSIG, "enign") & !str_detect(CLINSIG, "onflict") & !str_detect(CLINSIG, "ncertain") & !str_detect(CLINSIG, "risk_factor") & !str_detect(CLINSIG, "ffects") & !str_detect(CLINSIG, "drug_response")) %>% 
  filter((str_detect(Consequence, "missense") | str_detect(Consequence, "start_lost") | str_detect(Consequence, "stop_lost") & (str_length(REF)==1 & str_length(ALT)==1))) %>% 
  rowwise() %>%
  mutate(END=POS+str_length(ALT))
clinvar_annot_patho <- clinvar_annot %>% 
  filter(str_detect(CLINSIG, "athogenic")) %>% 
  filter((str_detect(Consequence, "missense") | str_detect(Consequence, "start_lost") | str_detect(Consequence, "stop_lost") & (str_length(REF)==1 & str_length(ALT)==1))) %>% 
  rowwise() %>%
  mutate(END=POS+str_length(ALT))
uniprot_domains_ranges <- GenomicRanges::GRanges(seqnames = uniprot_corrected$CHROM, ranges = IRanges::IRanges(start = uniprot_corrected$Start, end = uniprot_corrected$End))
clinvar_benign_ranges <- GenomicRanges::GRanges(seqnames = clinvar_annot_benign$CHROM, ranges = IRanges::IRanges(start = clinvar_annot_benign$POS, end = clinvar_annot_benign$END))
clinvar_patho_ranges <- GenomicRanges::GRanges(seqnames = clinvar_annot_patho$CHROM, ranges = IRanges::IRanges(start = clinvar_annot_patho$POS, end = clinvar_annot_patho$END))
uniprot_clinvar_overlap_noB <- IRanges::subsetByOverlaps(uniprot_domains_ranges, clinvar_benign_ranges, invert = T)
uniprot_clinvar_overlap_noB_1P <- IRanges::subsetByOverlaps(uniprot_clinvar_overlap_noB, clinvar_patho_ranges)
uniprot_clinvar_overlap_noB_1P <- GenomicRanges::reduce(uniprot_clinvar_overlap_noB_1P)
uniprot_clinvar_overlap_noB <- IRanges::subsetByOverlaps(uniprot_clinvar_overlap_noB, uniprot_clinvar_overlap_noB_1P, invert = T)
uniprot_clinvar_overlap_noB <- GenomicRanges::reduce(uniprot_clinvar_overlap_noB)
uniprot_clinvar_overlap_noB_1P <- data.frame(Chr=GenomicRanges::seqnames(uniprot_clinvar_overlap_noB_1P),
                                             Start=GenomicRanges::start(uniprot_clinvar_overlap_noB_1P),
                                             End=GenomicRanges::end(uniprot_clinvar_overlap_noB_1P))
uniprot_clinvar_overlap_noB <- data.frame(Chr=GenomicRanges::seqnames(uniprot_clinvar_overlap_noB),
                                          Start=GenomicRanges::start(uniprot_clinvar_overlap_noB),
                                          End=GenomicRanges::end(uniprot_clinvar_overlap_noB))

write_delim(uniprot_clinvar_overlap_noB, "data/processed/uniprot_domains_no_benign.tsv.gz")
write_delim(uniprot_clinvar_overlap_noB_1P, "data/processed/uniprot_domains_no_benign_1_patho.tsv.gz")
