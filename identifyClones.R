# Kenneth B. Hoehn 
# 11/10/2022
# Combine BCR and seurat information
# Identify clonal clusters
# Calculate SHM
# Identify convergent sequence clusters
# Match sequences to Cov-AbDab

library(alakazam)
library(scoper)
library(dplyr)
library(Seurat)
library(shazam)
library(dowser)
library(ggpubr)
library(ggtree)

print(sessionInfo())

nproc=10

# timepoint names
ovn = c(
  "v1D14"="v1D14",
  "v2D6"="v2D6",
  "v2D9"="v2D9",
  "v2D14"="v2D14",
  "v2D28"="v2D28",
  "v2M6"="M6")

# TRUE if you want to perform data consolidation and clonal clustering
if(FALSE){

# processed Seurat object
s = readRDS("raw/Moderna.integrated.ken.rds")

# collapse plasmablast annotations
cluster_keys = as.character(unique(Idents(s)))
names(cluster_keys) = cluster_keys
cluster_keys[grepl("lasma",cluster_keys)] = "Plasmablast"
cluster_keys[grepl("PB",cluster_keys)] = "Plasmablast"
Idents(s) = cluster_keys[as.character(Idents(s))]

# housekeeping for sample ids
s$sample = gsub("sample","",unlist(lapply(strsplit(gsub(">","",rownames(s@meta.data)),split="<"),function(x)x[2])))
s$pool = gsub("pool","",unlist(lapply(strsplit(gsub(">","",rownames(s@meta.data)),split="<"),function(x)x[3])))
s[["celltype"]] = Idents(s)

# larger pool of BCRs, before GEX based filtering
barcodes = read.csv("raw/data_new.csv")
barcodes$sample = gsub("sample","",unlist(lapply(strsplit(gsub(">","",barcodes$id),split="<"),function(x)x[2])))
barcodes$pool = gsub("pool","",unlist(lapply(strsplit(gsub(">","",barcodes$id),split="<"),function(x)x[3])))

files = list.files("raw/vdj", recursive=TRUE, full.names=TRUE)
files = files[grepl("productive-T.tsv", files)]

# between 4% and 43% of BCR barcodes aren't found in combined seurat object
total = tibble()
for(file in files){
  print(file)
  
  fs = strsplit(file,split="/")[[1]]
  pool = gsub("Pool","",fs[3])
  sample = gsub("VDJ","",fs[4])

  if(grepl("3_D2M6", pool) || grepl("3_6mo", pool)){
    pool = "3"
    sample = gsub("VDJ","",fs[5])
  }

  t = readChangeoDb(file)
  t$new_barcode = paste0(gsub("-1","",t$cell_id),"-<sample",sample,"><pool",pool,">")

  m = match(t$new_barcode,rownames(s@meta.data))
  print(paste(pool, sample, mean(is.na(m))))

  mean(t$new_barcode == rownames(s@meta.data)[m],na.rm=TRUE)
  t = bind_cols(t, s@meta.data[m,])
  t$pool = pool
  t$sample = sample
  sort = sapply(strsplit(t$orig.ident, split="\\."),function(x)x[1])

  donor_m = match(t$new_barcode, barcodes$id)
  print(paste("unfiltered barcode:",mean(is.na(donor_m))))
  if(sum(t$Donor != barcodes$donor[donor_m],na.rm=TRUE) > 0){
      stop("Donor mismatch!")
  }
  if(sum(ovn[t$timepoint] != barcodes$timepoint[donor_m],na.rm=TRUE) > 0){
      stop("timepoint mismatch!")
  }
  if(sum(sort != barcodes$sort[donor_m],na.rm=TRUE) > 0){
      stop("sort mismatch!")
  }
  if(sum(!is.na(t$Donor) & is.na(barcodes$donor[donor_m]))){
    stop(paste("cells in seurat object but not barcode object (shouldn't happen)",
    sum(!is.na(t$Donor) & is.na(barcodes$donor[donor_m])))) 
  }
  print(paste("cells in barcode object but not in Seruat object",
    sum(is.na(donor_m))))
  t$Donor = barcodes$donor[donor_m]
  t$orig.ident = paste0(barcodes$sort[donor_m],".",barcodes$timepoint[donor_m])

  total = bind_rows(total, t)
}

# 23% of cells are not associated with a donor
# due to not matching
print(mean(is.na(total$Donor)))

ftotal = filter(total, !is.na(Donor) & productive)

# parse sort and day
identsplit = strsplit(ftotal$orig.ident, split="\\.")
ftotal$sort = unlist(lapply(identsplit, function(x)x[1]))
ftotal$timepoint = unlist(lapply(identsplit, function(x)x[2]))

# typo in metadata
ftotal$sort[ftotal$sort == "S2Ppos"] = "S2pos"
table(ftotal$Donor,ftotal$timepoint,ftotal$sort)

writeChangeoDb(ftotal, "processed/combined.tsv")
ftotal = readChangeoDb("processed/combined.tsv")

ftotal$group = ftotal$sort
ftotal$group[ftotal$group == "S2neg"] = "S-2P-"
ftotal$group[ftotal$group == "S2pos"] = "S-2P+"
ftotal$modality = "10X"

# Scoper dies if multiple heavy chains
multi_heavy = table(filter(ftotal, locus=="IGH" & !is.na(new_barcode))$new_barcode)
multi_heavy_cells = names(multi_heavy)[multi_heavy > 1]
ftotal = filter(ftotal, !new_barcode %in% multi_heavy_cells)

ftotal$sequence_id = paste0(ftotal$sequence_id,"-",ftotal$sample,"-",ftotal$pool)

print(table(ftotal$modality, ftotal$Donor))

# check clonal thresholds
dist_cross = distToNearest(filter(ftotal, locus=="IGH"),
        sequenceColumn="cdr3", 
        vCallColumn="v_call", jCallColumn="j_call",
        model="ham", normalize="len", nproc=nproc,
        cross="Donor")

pdf("results/crossDistance.pdf",height=20,width=8)
print(
  ggplot(subset(dist_cross, !is.na(cross_dist_nearest)), 
             aes(x=cross_dist_nearest)) + 
    theme_bw() + 
    xlab("Cross-sample_id Hamming distance") + 
    ylab("Count") +
    geom_histogram(color="white", binwidth=0.02) +
    geom_vline(xintercept=0.1, color="firebrick", linetype=2) +
    facet_grid(Donor ~ ., scales="free_y"))
dev.off()

print(table(dist_cross$modality, dist_cross$Donor))

plots = list()
thresholds = list()
pdf("results/dist_to_nearest.pdf",width=6,height=6)
for(Donor in unique(dist_cross$Donor)){
   print(Donor)
   temp = filter(dist_cross,!!Donor==Donor)
   print(table(temp$modality))
   dist_ham <- distToNearest(filter(temp,locus=="IGH"), sequenceColumn="cdr3", 
    vCallColumn="v_call", jCallColumn="j_call",
    model="ham", normalize="len", nproc=nproc)
   threshold_output <- findThreshold(dist_ham$dist_nearest,
    method = "gmm", model = "gamma-norm",
    cross=dist_cross$cross_dist_nearest,
    cutoff = "user", spc = 0.995)
    
    threshold <- threshold_output@threshold
    thresholds[[Donor]] = threshold
    print(threshold)
   g = ggplot(subset(dist_ham, !is.na(dist_nearest)),aes(x=dist_nearest,
    ,y = ..density..)) + 
     theme_bw() + 
     xlab("Hamming distance") + 
     ylab("Count") +
     scale_x_continuous(breaks=seq(0, 1, 0.1)) +
     geom_histogram(color="white", binwidth=0.02) +
     ggtitle(paste(Donor,threshold))+
     geom_histogram(
       aes(x=cross_dist_nearest,y = -..density..),
       color="white", binwidth=0.02,fill="black")+
     xlim(0,max(filter(dist_cross,
       !is.na(cross_dist_nearest))$cross_dist_nearest))+
     geom_vline(xintercept=0.1,color="grey")
   if(!is.na(threshold)){
       g = g + geom_vline(xintercept=threshold, color="firebrick", linetype=2)
   }
   plots[[Donor]] = g
   print(g)
}
dev.off()

clones = tibble()
for(donor in unique(ftotal$Donor)){
   temp = as.data.frame(hierarchicalClones(filter(ftotal, Donor == donor),
    threshold=thresholds[[donor]], cell_id = "new_barcode", locus = "locus",
    only_heavy = FALSE, split_light = TRUE, cdr3 = TRUE, nproc = nproc,
    verbose = FALSE, log = NULL, summarize_clones = TRUE))
  clones = bind_rows(clones, temp)
}
clones = filter(clones, !is.na(clone_id))
clones$clone_id = paste0(clones$Donor,"-",clones$clone_id)

print(table(clones$modality, clones$Donor))

saveRDS(clones,"processed/clones.rds")

# Create Germlines
references = readIMGT(dir = "~/share/germlines/imgt/human/vdj")

comb_germline = createGermlines(clones, references, nproc=nproc)

print(table(comb_germline$modality, comb_germline$Donor))

# calculate SHM frequency in the V gene
comb_germline <- observedMutations(comb_germline, 
      sequenceColumn="sequence_alignment",
      germlineColumn="germline_alignment_d_mask",
      regionDefinition=IMGT_V,
      frequency=TRUE,
      combine=TRUE, 
      nproc=nproc)

writeChangeoDb(comb_germline,"processed/all_cloned_data.tsv")

}# if(TRUE)

comb_germline = readChangeoDb("processed/all_cloned_data.tsv")

# Identify shared clones
pclones = hierarchicalClones(select(comb_germline, -clone_id), threshold=0.2, method="aa",
    cell_id = "new_barcode", locus = "locus",
    only_heavy = FALSE, split_light = TRUE, cdr3 = TRUE, nproc = nproc,
    verbose = FALSE, log = NULL, summarize_clones = TRUE)
pclones = ungroup(as.data.frame(pclones))

saveRDS(pclones, "processed/pclones.rds")

m = match(comb_germline$sequence_id, pclones$sequence_id)
comb_germline$convergent_cluster = pclones$clone_id[m]

number_donors = comb_germline %>%
  group_by(convergent_cluster) %>%
  summarize(convergent_cluster_donors = n_distinct(Donor))

m = match(comb_germline$convergent_cluster, number_donors$convergent_cluster)
comb_germline$convergent_cluster_donors = number_donors$convergent_cluster_donors[m]

writeChangeoDb(comb_germline,"processed/all_cloned_data.tsv")

# Return sequences that closely match to a particular sequence
# in the CoV-AbDab
match_aa = function(i, subject, query){
  if(i %% 10 == 0){print(paste(i, nrow(query)))}
  temp = query[i,]
  v = getGene(temp$Heavy.V.Gene)
  j = getGene(temp$Heavy.J.Gene)
  cdrh3 = gsub("\\s+","",temp$CDRH3)
  pmatch = filter(subject, locus == "IGH" &
    grepl(paste0(v,"\\*"), v_call) & 
    grepl(paste0(j,"\\*"), j_call) &
    cdr3_aal == nchar(cdrh3))
  if(nrow(pmatch) > 0){
    pmatch$aa_dist = unlist(lapply(pmatch$cdr3_aa, function(x)
      seqDist(cdrh3, x, getAAMatrix())))/nchar(cdrh3)
    pmatch$match = temp$Name
  }
  return(pmatch)
}

# Map to publicly characterized covid antibodies
comb_germline$cdr3_aa = translateDNA(comb_germline$cdr3)
comb_germline$cdr3_aal = nchar(comb_germline$cdr3_aa)
covid = read.csv("CoV-AbDab_031022.csv")
covid = filter(covid, grepl("uman",Heavy.V.Gene))
names(covid)[1] = "Name"

# find matches
matches = bind_rows(parallel::mclapply(1:nrow(covid),
  function(x)match_aa(x, comb_germline, covid),
  mc.cores=nproc))

saveRDS(matches, "processed/matches.rds")

# process match information
epitope = covid$Protein...Epitope
names(epitope) = covid$Name

neut = covid$Neutralising.Vs
names(neut) = covid$Name

clone_size = table(filter(comb_germline, locus=="IGH")$clone_id)

matches_long = matches %>%
  filter(aa_dist <= 0.2) %>%
  select(sequence_id, Donor, timepoint, group, c_call, celltype, clone_id, match, aa_dist) %>%
  mutate(clone_size=as.numeric(clone_size[clone_id]), 
    epitope = epitope[match], neutralizing=neut[match]) %>%
  arrange(desc(clone_size), clone_id)

matches_seq = matches_long %>%
  group_by(sequence_id, Donor, timepoint, group, c_call, celltype, clone_id, clone_size) %>%
  summarize(
    epitopes = paste(unique(epitope),collapse=","),
    matches = paste(unique(match),collapse=","),
    neutralizing = paste(unique(neutralizing),collapse=",")) %>%
  arrange(desc(clone_size), clone_id)

write.csv(matches_long, file="results/matches_long.csv",row.names=FALSE)
write.csv(matches_seq,  file="results/matches.csv",row.names=FALSE)
matches_long = read.csv("results/matches_long.csv")
matches_seq =  read.csv("results/matches.csv")

# add match information to BCR information
m = match(comb_germline$sequence_id, matches_seq$sequence_id)
comb_germline$public = !is.na(matches_seq$matches[m])
comb_germline$matches = matches_seq$matches[m]
comb_germline$epitopes = matches_seq$epitopes[m]
comb_germline$neutralizing = matches_seq$neutralizing[m]

# save annotated cloned data
saveRDS(comb_germline, "processed/all_public_cloned_data.rds")

