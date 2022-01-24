## get the files
files <- list.files(".","_single_clone.tsv", recursive=T, full.names=T)
names(files) <- basename(dirname(files))

## filter groups which are not relevant

# from P3, taking prior to vaccine
pat <- "P3_I[0-3]_S1"
files <- files[!grepl(pat, names(files))]

# from P9 taking naive and non naive, not combined
pat <- "P9_I[0-9]+_S3"
files <- files[!grepl(pat, names(files))]

# from P7 taking naive and non naive, not combined
pat <- "P7_I[0-9]+_S[2-3]"
files <- files[!grepl(pat, names(files))]

# not taking P10
pat <- "P10_"
files <- files[!grepl(pat, names(files))]

# read and combine files
dbs <- data.table::rbindlist(lapply(files,data.table::fread,stringsAsFactors=F,  select = c("sequence_id", "v_germline_end", "sequence_alignment", "germline_alignment_d_mask", "v_call", "j_call", "productive", "c_call") ),use.names=T, idcol="subject",fill=T)

# take only reads which start from the first nuc, filter V sequences wich has less than 312 nucs, has gaps or Ns

dbs_filter <- dbs %>% filter(!grepl("^[.]",sequence_alignment))%>% rowwise() %>% mutate(
  v_seq = substr(sequence_alignment, 1, v_germline_end),
  v_seq_germ = substr(germline_alignment_d_mask, 1, v_germline_end),
  project = strsplit(subject,"_")[[1]][1]) %>% filter(!grepl("N", v_seq), !grepl("-", v_seq), nchar(v_seq)>312)


dbs_filter <- shazam::observedMutations(dbs_filter, sequenceColumn="sequence_alignment",
                                 germlineColumn="germline_alignment_d_mask",
                                 regionDefinition=shazam::IMGT_V,
                                 frequency=TRUE,
                                 nproc=30)



allele_diff <- function(x,y) {
  germs <- c(x,y)
  germs <- lapply(germs, function(x)
    strsplit(x, '')[[1]])
  germs_m <- t(sapply(germs, `length<-`, max(lengths(germs))))
  setdiff_mat <- function(x) {
    sum(!unique(x) %in% c('.', NA, "N"))#
  }
  idx = which(apply(germs_m, 2, setdiff_mat) > 1)
  return(length(idx))
}

dbs_filter <- setDT(dbs_filter)
dbs_filter <- dbs_filter[,mut := allele_diff(v_seq, v_seq_germ), by = seq_len(nrow(dbs_filter))]
dbs_filtered_mut3 <- dbs_filter[mut<=3]

save(dbs_filtered_mut3, file = "data_p1_to_P9_filtered_mut3.rda")


data_ <- data.table(dbs_filtered_mut3)
data_[,v_call_new := stringi::stri_replace_all_regex(str = v_call, pattern = paste0(gsub("[*]","[*]",alleles_db$or_allele),"\\b"), replacement = alleles_db$new_allele, vectorize_all = F)]
data_[, v_gene := alakazam::getGene(v_call_new, first = F, collapse = T, strip_d = F)]
data_[, .(count = .N, multi = sum(grepl(",",v_gene))), by = v_gene]
data_ <- data_[!grepl(",",v_gene)]
data_[,v_allele := stringr::str_replace_all(gsub("(IG[HKL][VDJADEGMC]|TR[ABDG])[A-Z0-9\\(\\)]+[-/\\w]*[*]", "", v_call_new, perl = T), "([0-9][0-9]+),\\1", "\\1")]
data_[, .(count_allele = .N, multi_allele = sum(grepl(",",v_allele))), by = v_allele]
data_ <- data_[!(grepl(",",v_allele))]
data_ <- data_[!is.na(v_allele)]
data_[,v_allele2 := v_allele]

save(data_, file = "data_p1_to_P9_filtered_mut3_for_frac.rda")

cluster <-
  parallel::makeCluster(10, type = "FORK")
#Export cluster
parallel::clusterExport(cl = cluster,
                        c("allele_diff", "dbs_filter"),
                        envir =
                          environment())

# calculating mutation between IMGT sequence and the germline sequence, selecting a single sequence to each clone with the fewest mutations
dbs_filter$mut <-
  parSapply(cluster, 1:nrow(dbs_filter), function(i){
    length(allele_diff(c(dbs_filter$v_seq[i],dbs_filter$v_seq_germ[i])))
  })

dbs_filter <- setDT(dbs_filter)
dbs_filtered_mut3 <- dbs_filter[mut<=3]

save(dbs_filtered_mut3, file = "data_p1_to_P9_filtered_mut3.rda")

load(file = "alleles_dbs.rda")
load(file = "functional_groups.rda")
data_ <- data.table(dbs_filtered_mut3)
alleles_db <- alleles_dbs$IGH$functional$nonsingle$all$`318`$complete$`95`
data_[,v_call_new := stringi::stri_replace_all_regex(str = v_call, pattern = paste0(gsub("[*]","[*]",alleles_db$or_allele),"\\b"), replacement = alleles_db$new_allele, vectorize_all = F)]
data_[, v_gene := alakazam::getGene(v_call_new, first = F, collapse = T, strip_d = F)]

data_multi <- data_[,v_gene2 := v_gene]
data_multi <- data_multi[, .(v_gene2 = unlist(tstrsplit(v_gene2, ",", type.convert = FALSE))), by = setdiff(names(data_multi), "v_gene2")]
data_multi_sum <- data_multi[, .(count = .N, multi = sum(grepl(",",v_gene))), by = list(v_gene2,subject)]

data_filter <- data_[!grepl(",",v_gene)]
data_filter[,v_allele := stringr::str_replace_all(gsub("(IG[HKL][VDJADEGMC]|TR[ABDG])[A-Z0-9\\(\\)]+[-/\\w]*[*]", "", v_call_new, perl = T), "([0-9][0-9]+),\\1", "\\1")]

data_multi_a <-data_filter[,v_allele2 := v_allele]
data_multi_a <- data_multi_a[, .(v_allele2 = unlist(tstrsplit(v_allele2, ",", type.convert = FALSE))), by = setdiff(names(data_multi_a), "v_allele2")]
data_multi_a_sum <- data_multi_a[, .(count = .N, multi = sum(grepl(",",v_allele))), by = v_allele2]

data_filter_frac <- data_filter[!(grepl(",",v_allele))]
data_filter_frac <- data_filter_frac[!is.na(v_allele)]
data_filter_frac[,v_allele2 := v_allele]


save(data_multi_sum, data_multi_a_sum, file = "data_P1_to_P9_multi.rda")

## get the germlines
files_v <- list.files(".","V_extended_ref.fasta", recursive=T, full.names=T)
names(files_v) <- basename(dirname(files_v))

# from P3, taking prior to vaccine
pat <- "P3_I[0-3]_S1"
files_v <- files_v[!grepl(pat, names(files_v))]

# from P9 taking naive and non naive, not combined
pat <- "P9_I[0-9]+_S3"
files_v <- files_v[!grepl(pat, names(files_v))]

# from P7 taking naive and non naive, not combined
pat <- "P7_I[0-9]+_S[2-3]"
files_v <- files_v[!grepl(pat, names(files_v))]

# not taking P10
pat <- "P10_"
files_v <- files_v[!grepl(pat, names(files_v))]

vgerm <- unlist(unname(sapply(files_v,tigger::readIgFasta)))
vgerm <- vgerm[!duplicated(names(vgerm))]



## mat list


load("data.rda")
load("../vgerm_igh_with_novel.rda")
load("../functionality.rda")
library(DECIPHER)

vgerms$IGH <- vgerm

vgerm <- vgerm[!grepl(paste0(c("NL", "OR"), collapse = "|"), names(vgerm))]
vgerm <- vgerm[!names(vgerm) %in% functionality$IGH]
vgerm <- vgerm[!is.na(vgerm)]
vgerm <- gsub("[.]", "-", vgerm)
tmp_first <- stringi::stri_locate(vgerm, regex = "[A-Z|a-z]")     
for (i in 1:nrow(tmp_first)) {
  vgerm[[i]] <-
    paste0(gsub("[-]", "N", substr(vgerm[[i]], 1, tmp_first[i, 1] - 1)),
           substr(vgerm[[i]], tmp_first[i, 1], nchar(vgerm[[i]])))
}  

vgerm <-
  sapply(vgerm, function(x)
    ifelse(nchar(x) != max(nchar(vgerm)), paste0(x, paste0(
      rep("N", (max(nchar(
        vgerm
      )) - nchar(x))), collapse = ""
    )), x))

## cut to 3' end value
vgerm <-
  sapply(vgerm, function(x)
    substr(x, 1, 318)
  )

dna <- DNAStringSet(vgerm)

mat <-
  DistanceMatrix(
    dna,
    includeTerminalGaps = FALSE,
    penalizeGapGapMatches = FALSE,
    penalizeGapLetterMatches = T,
    verbose = F
  )

mat_list$IGH$functional$nonsingle$all$`318` <- mat


allele_groups <- function(tmp, genes, gr){
  new_alleles <- c()
  diag(tmp)=NA
  val <- tmp[rowSums(tmp==0,na.rm = T)!=0,colSums(tmp==0, na.rm = T)!=0]
  all_alleles <- rownames(tmp)
  gene <- getGene(all_alleles[1])
  similar <- which(val == 0, arr.ind = T)
  rm_a <- c()
  rm_s <- c()
  if(length(similar)!=0){
    similar[,1] <- rownames(val[similar[,1],])
    similar[,2] <- colnames(val[,as.numeric(similar[,2])])
    similar <- as.data.frame(similar, stringsAsFactors=F)
    similar$names <- sapply(1:nrow(similar), function(i) paste0(sort(similar[i,]), collapse = ","))
    idx_remove <- !duplicated(similar[,3])
    similar <- similar[idx_remove,1:2]
    for(ii in 1:nrow(similar)){
      a <- similar[ii,]
      a_n <- which(getGene(a,strip_d = F)!=gene)
      if(length(a_n)>1){
        rm_a <- c(rm_a,a[[2]])
        rm_s <- c(rm_s,a[[1]])
      }else{
        if(length(a_n)==1){
          rm_a <- c(rm_a,a[[a_n]])
          rm_s <- c(rm_s,a[[which(getGene(a,strip_d = F)==gene)]])
        }
      }
    }
    rm_a <- unique(rm_a)
    f_alleles <- all_alleles[which(!all_alleles %in% rm_a)]
  }else{
    f_alleles <- all_alleles
  }
  n <- length(f_alleles)
  for(x in 1:n){
    new <- paste0(gene,"G",gr,"*",ifelse(x<10,"0",""),x)
    new_alleles <- dplyr::bind_rows(data.frame(or_allele = f_alleles[x], new_allele = new, remove = F, stringsAsFactors = F), new_alleles)
  }
  if(length(rm_a)!=0) new_alleles <- dplyr::bind_rows(data.frame(or_allele = rm_a, new_allele = new_alleles$new_allele[new_alleles$or_allele %in% unique(rm_s)], remove = T, stringsAsFactors = F), new_alleles)
  
  return(new_alleles)
}

thresh = 95
method = "complete"
range_seq <- 1 - as.numeric(thresh) / 100
hc_sub <- hclust(as.dist(mat), method)
row_order <- labels(as.dendrogram(hc_sub, hang = -1))
cut <- data.frame(cluster = dendextend::cutree(as.dendrogram(hc_sub, hang = -1), h = range_seq))
hc_sub <- cut
colnames(hc_sub) <- paste0((1-range_seq)*100, "%")

library(alakazam)
imgt <- unique(getGene(rownames(hc_sub), strip_d = F, omit_nl = F))
imgt <- setNames(1:length(imgt), imgt)
hc_sub$IMGT <- imgt[getGene(rownames(hc_sub), strip_d = F, omit_nl = F)]

hc_sub$labels <- factor(rownames(hc_sub), levels = row_order)
hc_sub$labels_ind <- setNames(1:length(levels(hc_sub$labels)), levels(hc_sub$labels))[hc_sub$labels]
hc_plot <- reshape2::melt(hc_sub, id.vars = c("labels", "labels_ind"))

hc_plot$freq <- 1

hc_plot$variable <- factor(hc_plot$variable, levels = c("IMGT", paste0((1-range_seq)*100, "%")))

groups <-
  setNames(hc_plot$value[hc_plot$variable != "IMGT"], hc_plot$labels[hc_plot$variable !=
                                                                       "IMGT"])
library(dplyr)
alleles_db <- c()
for (gr in unique(groups)) {
  genes <- names(groups)[groups == gr]
  if (length(unique(getGene(genes, strip_d = F, omit_nl = F))) == 1){
    new <- sapply(genes, function(x) paste0(strsplit(x,"[*]")[[1]][1],"G",gr,"*",strsplit(x,"[*]")[[1]][2]))
    alleles_db <- bind_rows(alleles_db, data.frame(or_allele = genes, new_allele = new, remove = F))
    
  }else{
    prox_mat <- mat[genes, genes]
    
    alleles <- allele_groups(prox_mat, genes, gr = gr)
    alleles_db <- bind_rows(alleles_db, alleles)
  }
}

alleles_dbs$IGH$functional$nonsingle$all$`318`$complete$`95` <- alleles_db

alleles_db_l <-
  setNames(alleles_db$new_allele, alleles_db$or_allele)

groups_genes <- groups
names(groups_genes) <-
  sapply(names(groups_genes), function(x)
    if (x %in% names(alleles_db_l))
      alleles_db_l[x]
    else
      x)

tmp <-
  merge(
    data.frame(
      alleles = names(groups),
      cluster = groups,
      stringsAsFactors = F
    ),
    data.frame(
      alleles_new = names(groups_genes),
      cluster = groups,
      stringsAsFactors = F
    ),
    by = "cluster"
  )
tmp <-
  tmp %>% group_by(cluster) %>% dplyr::summarise(
    gene = paste0(unique(
      alakazam::getGene(alleles, strip_d = F)
    ), collapse = "/"),
    new_gene = paste0(unique(
      alakazam::getGene(alleles_new, strip_d = F)
    ), collapse = "/")
  )
genes_groups <- setNames(tmp$new_gene, tmp$gene)
functional_groups$IGH$functional$nonsingle$all$`318`$complete$`95` <- genes_groups



save(alleles_dbs, file = "../reference_book2/alleles_dbs.rda")
save(functional_groups, file = "../reference_book2/functional_groups.rda")

ids <- which(!names(vgerms$IGH) %in% alakazam::getAllele(names(vgerms_full$IGH), strip_d = F, omit_nl = F))
for(i in ids){
  tmp <- vgerms_full$IGH[420]
  names(tmp) <- names(vgerms$IGH[i])
  vgerms_full$IGH <- c(vgerms_full$IGH,tmp)
}

save(mat_list, vgerms, vgerms_msa, vgerms_full, functionality, color, col_vector, col_vec_dend, file = "../reference_book2/data.rda")

#### calc the fractions based on the new vgerm


data_frac <- sapply(c("IGH"), function(chain) {
  sapply(c("functional"), function(functional) {
    sapply(c("nonsingle"), function(single) {
      sapply(c("all"), function(rm_short) {
        sapply(as.character(318), function(seq_end) {
          sapply(c("complete"), function(method) {
            sapply(as.character(95), function(thresh) {
              
              tmp <- lapply(0:3, function(mut_val) {
                
                alleles_db <- alleles_dbs$IGH$functional$nonsingle$all$`318`$complete$`95`
                genotypes <- data_filter_frac[mut <= mut_val]
                
                genotypes[,n_row_sub := .N, by = subject]
                genotypes[,n_row := .N, by = list(subject,v_gene)]
                genotypes <- genotypes[, .(v_allele = unlist(tstrsplit(v_allele, ",", type.convert = FALSE))), by = setdiff(names(genotypes), "v_allele")]
                genotypes[,frac := 1/1]

                genotypes_fraction <- genotypes %>% ungroup() %>% dplyr::group_by(subject, v_gene, v_allele) %>%
                  dplyr::summarise(count = round(sum(frac),3), freq = round(sum(frac)/unique(n_row),8),
                                   freq2 = round(sum(frac)/unique(n_row_sub),8)) %>% ungroup() %>%  dplyr::arrange(subject, v_gene,desc(freq)) %>%
                  dplyr::group_by(subject, v_gene) %>% dplyr::mutate(loc = 1:dplyr::n(), project = strsplit(subject,"_")[[1]][1])
                
                hetro_samples <- genotypes %>% dplyr::filter(grepl("J6",j_call), !grepl(",",j_call)) %>% 
                  dplyr::select(subject, j_call) %>% dplyr::group_by(subject) %>% dplyr::mutate(nrow = dplyr::n()) %>% dplyr::ungroup() %>% 
                  dplyr::group_by(subject, j_call) %>% dplyr::summarise(frac = dplyr::n()/nrow) %>% dplyr::slice(1) %>%
                  dplyr::group_by(subject) %>% dplyr::summarise(frac_J6 = frac[which(frac>=0.2 & frac<=0.8)], 
                                                                A = sum(frac>=0.2 & frac<=0.8)>=2, 
                                                                J_A = paste0(j_call[which(frac>=0.2 & frac<=0.8)]))  %>% 
                  dplyr::filter(A)
                
                genotypes_fraction$J6 <- ifelse(genotypes_fraction$subject %in% hetro_samples$subject, 1, 2)
                
                genotypes_fraction$J6_TAG <- ""
                
                # style_str <- "white-space: nowrap; border: 1px solid white; background: steelblue; height: 15px;"
                
                for(sub in unique(hetro_samples$subject)){
                  no_check <- "❌ " 
                  check <- "✅"
                  j02 <- paste0("</b>J6*02: ",ifelse(hetro_samples %>% filter(subject==sub, J_A == "IGHJ6*02") %>% pull(frac_J6) %>% length(), check, no_check)) #hetro_samples %>% filter(subject==sub, J_A == "IGHJ6*02") %>% pull(frac_J6), 0)*1000,"px")
                  j03 <- paste0("</b>J6*03: ",ifelse(hetro_samples %>% filter(subject==sub, J_A == "IGHJ6*03") %>% pull(frac_J6) %>% length(), check, no_check)) #hetro_samples %>% filter(subject==sub, J_A == "IGHJ6*03") %>% pull(frac_J6), 0)*1000,"px")
                  j04 <- paste0("</b>J6*04: ",ifelse(hetro_samples %>% filter(subject==sub, J_A == "IGHJ6*04") %>% pull(frac_J6) %>% length(), check, no_check)) #hetro_samples %>% filter(subject==sub, J_A == "IGHJ6*04") %>% pull(frac_J6), 0)*1000,"px")
                  genotypes_fraction$J6_TAG[genotypes_fraction$subject == sub] <- paste0(j02,j03,j04)
                }

                genotypes_fraction_j6 <- genotypes %>% ungroup() %>% filter(subject %in% hetro_samples$subject) %>% 
                  group_by(subject) %>%
                  filter(j_call %in% hetro_samples$J_A[hetro_samples$subject %in% subject]) %>%
                  filter(grepl("J6",j_call), !grepl(",",j_call)) %>% ungroup() %>%
                  dplyr::group_by(subject, v_gene) %>%
                  dplyr::mutate(n_row = n()) %>% ungroup() %>%
                  dplyr::group_by(subject, v_gene, v_allele, j_call) %>%
                  dplyr::summarise(count = round(sum(frac),3), freq = round(sum(frac)/unique(n_row),8),
                                   freq2 = round(sum(frac)/unique(n_row_sub),8)) %>% ungroup() %>%  
                  dplyr::arrange(subject, v_gene,desc(freq)) %>%
                  dplyr::group_by(subject, v_gene) %>% 
                  dplyr::mutate(loc = 1:dplyr::n(), project = strsplit(subject,"_")[[1]][1])
                
                genotypes_fraction_comb <- rbind(genotypes_fraction_j6,genotypes_fraction)
                
                genotypes_fraction_comb$mut <- mut_val
                return(genotypes_fraction_comb)
              })
              genotypes_fraction_comb <- rbindlist(tmp)
              return(genotypes_fraction_comb)
            }, simplify = F)
          }, simplify = F)
        }, simplify = F)
      }, simplify = F)
    }, simplify = F)
  }, simplify = F)
}, simplify = F)


save(data_frac, file = "data_frac_new2.rda")


#### copy files from the reference book to the "new" site

## get the files

files <- list.files(".",".rmd")

## match between file and the new gene group

for(f in files){
  # f = files[1]
  
  lines <- readLines(f)
  
  ### check the group genes and the index
  
  id <- lines[1]
  
  g_id <- strsplit(id, " - ")[[1]][2]
  g_gene <- gsub("# ","",strsplit(id, " - ")[[1]][1])
  
  g_names <- names(genes_groups)
  g_names <- sapply(g_names, function(x) paste0(sort(strsplit(x,"/")[[1]]),collapse = "/"))
  if(grepl("/",g_gene)) g_gene <- paste0(sort(strsplit(g_gene,"/")[[1]]),collapse = "/")
  
  g_id_new <- genes_groups[g_names==g_gene]
  if(length(g_id_new)==0){
    print(f)
    next()
    }
  if(length(g_id_new)>1){
    print(f)
    next()
  }else{
    if(paste0(g_gene,g_id)==g_id_new){
      file.copy(f,"~/Dropbox (BIU)/genotype_sub_sampling/sister_genes/reference_book2/")
    }else{
      # change first line
      lines[1] <- paste0(strsplit(id, " - ")[[1]][1]," - ","G",strsplit(g_id_new,"[0-9]G")[[1]][2])
      # gerp line of g_group
      lines[grep("g_group =", lines)] <- gsub(paste0(g_gene,g_id),g_id_new,lines[grep("g_group =", lines)])
      # write file
      writeLines(lines, paste0("~/Dropbox (BIU)/genotype_sub_sampling/sister_genes/reference_book2/", ifelse(as.numeric(strsplit(g_id_new,"[0-9]G")[[1]][2])<10, paste0('0',as.numeric(strsplit(g_id_new,"[0-9]G")[[1]][2])), as.numeric(strsplit(g_id_new,"[0-9]G")[[1]][2])), "-G",as.numeric(strsplit(g_id_new,"[0-9]G")[[1]][2]),".rmd"))
    }
  }
  
  
}

files2 <- list.files("../reference_book2/",".rmd")
for(i in 1:51){
  
  if(any(grepl(paste0("G",i,"[.]"), files2))){
    next()
  }
  else{
    
    lines <- readLines(files[40])
    
    ### check the group genes and the index
    
    g_group = genes_groups[i]
    
    # change first line
    lines[1] <- paste0("# ", names(g_group), " - G",i)
    # gerp line of g_group
    lines[grep("g_group =", lines)] <- gsub("IGHV4-4G40",g_group,lines[grep("g_group =", lines)])
    # write file
    writeLines(lines, paste0("~/Dropbox (BIU)/genotype_sub_sampling/sister_genes/reference_book2/", ifelse(i<10, paste0('0',i), i), "-G",i,".rmd"))
    
  }
  
}



### sort absolute dict
alleles_db$gene <- getGene(alleles_db$or_allele, strip_d = F, omit_nl = F)
alleles_db$gene_group <- getGene(alleles_db$new_allele, strip_d = F, omit_nl = F)
rownames(alleles_db) <- NULL

for(i in 1:length(absolute_thresholds_dict)){
  print(i)
  g_old <- names(absolute_thresholds_dict)[i]
  alleles <- paste0("IGH",names(absolute_thresholds_dict[[i]]))
  
  g_new <- unique(alleles_db$gene_group[alleles_db$or_allele %in% alleles])
  print(g_new)
  print("")
  print("###############")
  if(length(g_new)>1){
    for(g in g_new){
      print("")
      print(g)
      missing_alleles <- alleles_db$or_allele[alleles_db$gene_group %in% g]
      #missing_alleles <- missing_alleles[!missing_alleles %in% alleles]
      cat(paste0("'",gsub("IGH","",missing_alleles),"'=",0.0001), sep = ",")
    }
  }else{
    missing_alleles <- alleles_db$or_allele[alleles_db$gene_group %in% g_new]
    #missing_alleles <- missing_alleles[!missing_alleles %in% alleles]
    cat(paste0("'",gsub("IGH","",missing_alleles),"'=",0.0001), sep = ",")
  }
  print("")
  print("###############")
}


#########################3


files2 <- list.files("../reference_book2/",".rmd")
files2 <- files2[-8]
for(f in files2){
  
  lines <- readLines(f)
  
  id <- grep("v_calls <- unique", lines)
  
  if(length(id)!=0){
    lines2 <- lines[1:(id-2)]
    lines3 <- c("```{r}", "v_calls <- unique(data[grepl(g_group,v_gene),v_call])", "seq_align(v_calls, allele_db, vgerms, chain, mat, g_group)","```")
    lines4 <- lines[(id-1):(id+2)]
    lines4[3] <- gsub("seq_align","seq_align2", lines4[3])
    lines5 <- lines[(id+3):length(lines)]
    
    lines <- c(lines2, lines3, lines4, lines5)
    
    writeLines(lines, f)
  }
  
}

########### thresholds
source("functions.R")
dfs <- data.table::rbindlist(lapply(absolute_thresholds_dict, function(l) data.frame(or_allele = paste0("IGH",names(l)), thresh = unlist(l), stringsAsFactors = F, row.names = NULL)), use.names = T, fill = T, idcol = "func_group")
dfs_f <- dfs[dfs$or_allele %in% alleles_db$or_allele,]
dfs_f <- dfs_f[!duplicated(dfs_f),]
ll <- setNames(alleles_db$new_allele, alleles_db$or_allele)
dfs_f$new_allele <- ll[dfs_f$or_allele]
alleles_db$group <- getGene(alleles_db$new_allele, first = F, strip_d = F, omit_nl = F)

ll <- setNames(alleles_db$group, alleles_db$or_allele)
dfs_f$group <- ll[dfs_f$or_allele]

remove
ll <- setNames(alleles_db$new_allele, alleles_db$or_allele)
dfs_f$rm <- ll[dfs_f$or_allele]
dfs_f <- dfs_f %>% dplyr::group_by(new_allele) %>% dplyr::summarise(func_group = func_group, or_allele = paste0(or_allele, collapse = "/"), thresh = unique(thresh))
