# Kenneth B. Hoehn 
# 7/Feb/2023, revised 10/May/2023 
# Perform clonal diversity and overlap analyses
# as well as SHM analyses


library(alakazam)
library(scoper)
library(dplyr)
library(tidyr)
library(Seurat)
library(shazam)
library(dowser)
library(ggpubr)
library(ggtree)
library(gridExtra)
library(RColorBrewer)
library(tidyr)
library(UpSetR)

print(sessionInfo())

nproc = 10

reannotate = FALSE

# new new cluster names
new_names = c(
"C1 - Naive B cell 1"="Naive-C1",
"C2 - Naive B cell 2"="Naive-C2",
"C3 - unswitched MBC"="Unswitched MBC-C3",
"C4 - Memory B cell 1"="Switched MBC-C4",
"C5.1 - Mixed MBC"="MBC-SC5.1",
"C5.2 - RM"="MBC-SC5.2",
"C5.3 - CD38+ AM"="MBC-SC5.3",
"C5.4 - TLM"="MBC-SC5.4",
"C5.5 - CD38- AM"="MBC-SC5.5",
"C6 - Plasmablast"="PB-C6",
"C7 - Naive B cell 3"="Naive-C7")

# new timepoint names
ovn = c(
  "D1D14"="v1D14",
  "D2D6"="v2D6",
  "D2D9"="v2D9",
  "D2D14"="v2D14",
  "D2D28"="v2D28",
  "6mo"="M6")

# timepoint colors
timepoint_palette = c(
  "v1D14"="#313695",
  "v2D6"="#4575b4",
  "v2D9"="#abd9e9",
  "v2D14"="#ffffbf",
  "v2D28"="#f46d43",
  "M6"="#a50026")

# set order of clusters
level.order = c("C1 - Naive B cell 1", "C2 - Naive B cell 2", "C3 - unswitched MBC",
 "C4 - Memory B cell 1", "C5.1 - Mixed MBC", "C5.2 - RM",
  "C5.3 - CD38+ AM", "C5.4 - TLM", "C5.5 - CD38- AM", "C6 - Plasmablast", "C7 - Naive B cell 3")
level.order = new_names[level.order]

# Set up palettes and plotting objects
subisotype_levels = c("IGHM", "IGHD", "IGHG3", "IGHG1", "IGHA1", "IGHG2", "IGHG4", "IGHE", "IGHA2")
subisotype_palette = brewer.pal(length(subisotype_levels), "Set1")
names(subisotype_palette) = subisotype_levels

# cluster info and colors
cinfo = read.csv("cluster_key.csv")

old2new = cinfo$New
names(old2new) = cinfo$Old

palette = cinfo$Color
names(palette) = new_names[cinfo$New]

if(reannotate){
  #comb_germline = readChangeoDb("processed/all_cloned_data.tsv")
  comb_germline = readRDS("processed/all_public_cloned_data.rds")
  
  # convert to new cell type names
  comb_germline$old_celltype = comb_germline$celltype
  comb_germline$celltype = old2new[comb_germline$old_celltype]
  
  print(table(filter(comb_germline, is.na(celltype) & !is.na(old_celltype))$old_celltype))
  
  comb_germline$new_celltype = new_names[comb_germline$celltype]
  
  print(unique(comb_germline$new_celltype)[!unique(comb_germline$new_celltype) %in% level.order])
  comb_germline$new_celltype = factor(comb_germline$new_celltype, levels=level.order)
  table(comb_germline$new_celltype)
  
  # assumes D2 28 days after D1
  comb_germline$time = NA
  comb_germline$time[comb_germline$timepoint == "v1D14"] = 14
  comb_germline$time[comb_germline$timepoint == "v2D6"]  = 28+6
  comb_germline$time[comb_germline$timepoint == "v2D9"]  = 28+9
  comb_germline$time[comb_germline$timepoint == "v2D14"] = 28+14
  comb_germline$time[comb_germline$timepoint == "v2D28"] = 28+28
  comb_germline$time[comb_germline$timepoint == "M6"]   = 6*30
  
  print(sum(is.na(comb_germline$time)))
  
  # set timepoint names and order
  comb_germline$timepoint = factor(comb_germline$timepoint, levels=ovn)
  
  # Compare SHM frequency of each cell type, per donor, per timepoint
  sclones = unique(filter(comb_germline, sort=="S2pos" & locus=="IGH")$clone_id)
  comb_germline$s2clone = comb_germline$clone_id %in% sclones
  
  comb_germline$naive = grepl("Naive", comb_germline$new_celltype)
  comb_germline$c_call = factor(comb_germline$c_call, 
    levels=subisotype_levels)
  
  saveRDS(comb_germline, "processed/all_public_analyzed_cloned_data.rds")
  comb_germline = readRDS("processed/all_public_analyzed_cloned_data.rds")
  
  # remove cells with only light chains
  # this can happen if heavy chain germline reconstruction fails for heavy chains
  hc = unique(filter(comb_germline,locus=="IGH")$new_barcode)
  lc = filter(comb_germline, !new_barcode %in% hc)
  comb_germline = filter(comb_germline, new_barcode %in% hc)
  
  table(lc$convergent_cluster)
  
  saveRDS(comb_germline, "processed/all_public_analyzed_filtered_cloned_data.rds")
}
comb_germline = readRDS("processed/all_public_analyzed_filtered_cloned_data.rds")

# S2+ clones have any S2+ cells
comb_germline$S2 = ""
comb_germline$S2[comb_germline$s2clone] = "S-2P+"
comb_germline$S2[!comb_germline$s2clone] = "S-2P-"

# for each clone, tally number of cells and cell types
clone_celltype_counts = comb_germline %>% 
  filter(locus=="IGH") %>% 
  group_by(Donor, clone_id) %>% 
  summarize(n=n_distinct(new_barcode), 
   ntypes=n_distinct(new_celltype, na.rm=TRUE),
   nsort=n_distinct(group, na.rm=TRUE))

# number of clones
print(n_distinct(clone_celltype_counts$clone_id))
print(n_distinct(comb_germline$clone_id))
# clones with > 1 cell
print(sum(clone_celltype_counts$n > 1))
# clones with >= 10 cells
print(sum(clone_celltype_counts$n >= 10))
# clones with > 1 cell type
print(mean(clone_celltype_counts$ntypes > 1))
# clones with > 1 cell type among clones with > 1 cell
print(mean(filter(clone_celltype_counts, n > 1)$ntypes > 1))
# clones with > 1 cell type among clones with > 1 sort
print(mean(filter(clone_celltype_counts, n > 1)$nsort > 1))

# for each sort and clone, tally number of cells
comb_germline %>% 
  filter(locus=="IGH") %>% 
  group_by(Donor, group, clone_id) %>% 
  summarize(n=n_distinct(new_barcode)) %>%
  group_by(group) %>%
  summarize(n_expanded=sum(n >= 10), n_clones=n(),
    percent_expanded=mean(n >= 10)*100)

# what % of PBs are in S-2P+ clones?
# only include cells from the final seurat object
#s2pclones = unique(filter(comb_germline, sort=="S2pos" & !is.na(celltype) &
# locus=="IGH")$clone_id)
comb_germline %>%
  filter(locus=="IGH" & celltype=="C6 - Plasmablast" & !is.na(celltype)) %>%
  group_by(s2clone) %>%
  summarize(ncells = n_distinct(new_barcode)) %>%
  ungroup() %>%
  summarize(number = ncells[s2clone], total=ncells[s2clone] + ncells[!s2clone],
    mean=number/(ncells[s2clone] + ncells[!s2clone]))

# Table S4 double check
triple_donors = comb_germline %>%
  group_by(convergent_cluster) %>%
  filter(convergent_cluster_donors == 3) %>%
  summarize(n_donor = n_distinct(Donor),
    n_cells = n_distinct(new_barcode),
    ighv=paste(unique(alakazam::getGene(v_call[locus == "IGH"])),collapse=","),
    ighkl=paste(unique(alakazam::getGene(v_call[locus != "IGH"])),collapse=",")) %>%
  arrange(desc(n_cells)) %>% data.frame()

# Table S4 double check
multi_donors = comb_germline %>%
  group_by(convergent_cluster) %>%
  summarize(n_donor = n_distinct(Donor),
    n_cells = n_distinct(new_barcode),
    ighv=paste(unique(alakazam::getGene(v_call[locus == "IGH"])),collapse=","),
    ighkl=paste(unique(alakazam::getGene(v_call[locus != "IGH"])),collapse=",")) %>%
  filter(n_donor == 2 & n_cells >=5) %>%
  arrange(desc(n_cells))

write.table(triple_donors,file="results/triple_donors.tsv",sep="\t",quote=FALSE,row.names=FALSE)
write.table(multi_donors, file="results/multi_donors.tsv", sep="\t",quote=FALSE,row.names=FALSE)

# calculate clonal diversity of each cell type, per donor, per timepoint
comb_germline$population = paste(comb_germline$Donor, 
  comb_germline$S2, comb_germline$timepoint, sep="|")

# Clonal diversity analysis of all cells vs covid
sample_curve = tibble()
for(d in unique(comb_germline$Donor)){
    sample_c <- alphaDiversity(
      filter(comb_germline, locus=="IGH" & Donor == d & !naive &
        sort != "Plasmablast"),
      group="population", clone="clone_id",
      min_q=2.0, max_q=2.0, step_q=0.0,min_n=50,
      ci=0.95, nboot=500, uniform=TRUE)
    aq = filter(sample_c@diversity,q==2)
    sample_curve = bind_rows(sample_curve, aq)
}

sample_curve$Donor = unlist(lapply(strsplit(sample_curve$population,
  split="\\|"), function(x)x[1]))
sample_curve$S2 = unlist(lapply(strsplit(sample_curve$population,
  split="\\|"), function(x)x[2]))
sample_curve$timepoint = unlist(lapply(strsplit(sample_curve$population,
  split="\\|"), function(x)x[3]))
sample_curve$timepoint = factor(sample_curve$timepoint, levels=ovn)

jsize = 1
pdf("results/diversity_non-naive.pdf",width=6,height=2,useDingbats=FALSE)
print(ggplot(sample_curve, aes(x=timepoint, y=d, group=S2, color=S2)) + 
    geom_line() + geom_point() + 
    geom_errorbar(aes(ymin=d_lower,ymax=d_upper),width=0.2) +
    facet_wrap(. ~ Donor, scales="free_y")+theme_bw() +
    ylab("Clonal diversity") + xlab("Timepoint")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=1, 
    color = "black"), axis.text.y = element_text(color = "black")) +
    theme(strip.background =element_rect(fill="white"))
    )
dev.off()

# calculate clonal overlap among each timepoint
clustercounts = comb_germline %>%
    filter(!naive & grepl("S2", sort)) %>%
    group_by(Donor, S2, clone_id, timepoint) %>%
    summarize(n = n_distinct(new_barcode))
timepoints = levels(comb_germline$timepoint)
results = tibble()
for(d in unique(comb_germline$Donor)){
  print(d)
  for(s in unique(comb_germline$S2)){
    for(a in timepoints){
        for(b in timepoints){
          clone_counts = clustercounts %>%
            filter(Donor == d & S2 == s) %>%
            filter(timepoint %in% c(a,b)) %>%
            group_by(clone_id) %>%
            summarize(distinct = n_distinct(timepoint)) 
    
          intersect_clones = clone_counts %>%
            filter(distinct == 2) %>%
            pull(clone_id)
    
          if(a == b){
            intersect_clones = clone_counts %>%
            pull(clone_id)
          }
    
          union_clones = clone_counts %>%
            pull(clone_id)
    
          results = bind_rows(results,
            tibble(a=a,b=b,Donor=d,S2=s,
              intersect=length(intersect_clones),
              union=length(union_clones),
              jaccard=intersect/union))
      }
    }
  }
}
results$intersect = as.character(results$intersect)
results$union = as.character(results$union)

# make combined matrix
combined_matrix = tibble()
for(d in unique(comb_germline$Donor)){
  print(d)
  temp = filter(results, Donor==d)
  for(a in timepoints){
    for(b in timepoints){
      if(which(timepoints==a) > which(timepoints==b)){
        values = filter(temp, !!a==a & !!b==b & S2=="S-2P+")
      }else if(which(timepoints == a) < which(timepoints == b)){
        values = filter(temp, !!a==a & !!b==b & S2=="S-2P-")
      }else{
        sp = filter(temp, !!a==a & !!b==b & S2=="S-2P+")
        sn = filter(temp, !!a==a & !!b==b & S2=="S-2P-")
        sp$intersect = paste0(sn$intersect, "\n", sp$intersect)
        sp$union = paste0(sn$union, "\n", sp$union)
        values = sp
      }
      combined_matrix = bind_rows(combined_matrix, values)
    }
  }
}
combined_matrix[combined_matrix$a == combined_matrix$b,]$jaccard = NA
combined_matrix$a = factor(combined_matrix$a, levels=ovn)
combined_matrix$b = factor(combined_matrix$b, levels=ovn)

# Make condensed clonal overlap figure
pdf("results/clonal_overlap_timepoint_non-naive.pdf",width=5.25,
  height=2.25,useDingbats=FALSE)
print(
ggplot(combined_matrix,aes(x=a,y=b,fill=jaccard))+geom_tile()+
  facet_grid(.~Donor)+scale_fill_distiller(palette="RdYlBu")+
  theme_bw()+geom_text(aes(label=(intersect)),size=2)+
  labs(fill="Jaccard\nindex") + 
  xlab("S-2P+ Clones")+ylab("S-2P- Clones")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(strip.background =element_rect(fill="white")))
dev.off()

# calculate clonal overlap among each new_celltype
clustercounts = comb_germline %>%
    filter(grepl("S2",sort)) %>%
    group_by(Donor, S2, clone_id, new_celltype) %>%
    summarize(n = n_distinct(new_barcode))
celltypes = levels(comb_germline$new_celltype)
results = tibble()
for(d in unique(comb_germline$Donor)){
  print(d)
  for(s in unique(comb_germline$S2)){
    for(a in celltypes){
        for(b in celltypes){
          clone_counts = clustercounts %>%
            filter(Donor == d & S2 == s) %>%
            filter(new_celltype %in% c(a,b)) %>%
            group_by(clone_id) %>%
            summarize(distinct = n_distinct(new_celltype)) 
    
          intersect_clones = clone_counts %>%
            filter(distinct == 2) %>%
            pull(clone_id)
    
          if(a == b){
            intersect_clones = clone_counts %>%
            pull(clone_id)
          }
    
          union_clones = clone_counts %>%
            pull(clone_id)
    
          results = bind_rows(results,
            tibble(a=a,b=b,Donor=d,S2=s,
              intersect=length(intersect_clones),
              union=length(union_clones),
              jaccard=intersect/union))
        }
    }
  }
}
results$intersect = as.character(results$intersect)
results$union = as.character(results$union)

d52 = filter(results, S2=="S-2P+" & grepl("5\\.2",b)) %>% data.frame()
print(d52 %>% group_by(Donor) %>% arrange(jaccard) %>% data.frame())

# make combined matrix
combined_matrix = tibble()
for(d in unique(comb_germline$Donor)){
  print(d)
  temp = filter(results, Donor==d)
  for(a in celltypes){
    for(b in celltypes){
      if(which(celltypes == a) > which(celltypes == b)){
        values = filter(temp, !!a == a & !!b==b & S2=="S-2P+")
      }else if(which(celltypes == a) < which(celltypes == b)){
        values = filter(temp, !!a == a & !!b==b & S2=="S-2P-")
      }else{
        sp = filter(temp, !!a == a & !!b==b & S2=="S-2P+")
        sn = filter(temp, !!a == a & !!b==b & S2=="S-2P-")
        sp$intersect = paste0(sn$intersect, "\n", sp$intersect)
        sp$union = paste0(sn$union, "\n", sp$union)
        values = sp
      }
      combined_matrix = bind_rows(combined_matrix, values)
    }
  }
}
combined_matrix[combined_matrix$a == combined_matrix$b,]$jaccard = NA
combined_matrix$a = factor(combined_matrix$a, levels=celltypes)
combined_matrix$b = factor(combined_matrix$b, levels=celltypes)

# Make condensed clonal overlap figure
pdf("results/clonal_overlap_celltype.pdf",width=9,height=4,useDingbats=FALSE)
print(
ggplot(combined_matrix,aes(x=a,y=b,fill=jaccard))+geom_tile()+
  facet_grid(.~Donor)+scale_fill_distiller(palette="RdYlBu")+
  theme_bw()+geom_text(aes(label=(intersect)),size=2)+
  labs(fill="Jaccard index") + 
  xlab("S-2P+ Clones")+ylab("S-2P- Clones")+
  theme(axis.text.x = element_text(angle = 45, vjust =1, hjust=1)) +
  theme(strip.background =element_rect(fill="white")))
dev.off()

# C5.2 - RM clones at 6 months
mem6_clones = comb_germline %>%
  filter(timepoint %in% c("M6") & celltype == "C5.2 - RM") %>%
  pull(clone_id)

# clones with resting memory at M6
mem6_clones = unique(mem6_clones)

# multi timepoint clones with defined celltypes
multi_time = comb_germline %>%
  filter(!is.na(new_celltype)) %>%
  group_by(clone_id) %>%
  summarize(ntime=n_distinct(timepoint)) %>%
  filter(ntime > 1) %>%
  pull(clone_id)

# multi celltype clones
multi_celltype = comb_germline %>%
  filter(!is.na(new_celltype)) %>%
  group_by(clone_id) %>%
  summarize(ntype=n_distinct(celltype)) %>%
  filter(ntype > 1) %>%
  pull(clone_id)

# multi timepoint clone data, labeled as C5.2, M6
mtime = filter(comb_germline, clone_id %in% multi_time)
mtime$overlap = "No C5.2, M6"
mtime$overlap[mtime$clone_id %in% mem6_clones] = "C5.2, M6"

# number of cells, clones, and multi-timepoint clones from
# C5.2 at M6
m6_clone_stats = comb_germline %>%
  filter(timepoint %in% c("M6") & celltype == "C5.2 - RM") %>%
  group_by(S2, Donor) %>%
  summarize(
    cells=n_distinct(new_barcode),
    clones=n_distinct(clone_id),
    multi_timepoint_clones=n_distinct(clone_id[clone_id %in% multi_time]),
    multi_celltype_clones=n_distinct(clone_id[clone_id %in% multi_celltype]))

write.csv(m6_clone_stats, file="results/m6_clone_stats.csv", 
  quote=FALSE, row.names=FALSE)

# calculate frequency of each celltype at each timepoint
cluster_time_frequency = mtime %>%
  filter(!is.na(new_celltype)) %>%
  group_by(Donor, S2, overlap, timepoint, new_celltype) %>%
  summarize(n=n_distinct(new_barcode)) %>%
  mutate(freq = n/sum(n))

cluster_time_frequency$timepoint = factor(cluster_time_frequency$timepoint,
 levels=ovn)
cluster_time_frequency$label = cluster_time_frequency$n
cluster_time_frequency$label[cluster_time_frequency$freq < 0.15] = ""
cluster_time_frequency$new_celltype = factor(cluster_time_frequency$new_celltype, 
  levels=level.order)

cluster_time_frequency$S2 = factor(cluster_time_frequency$S2, 
  levels=c("S-2P-", "S-2P+"))
cluster_time_frequency$overlap = factor(cluster_time_frequency$overlap, 
  levels=c("No C5.2, M6", "C5.2, M6"))
pdf("results/s2_6mo_memory_clones.pdf",width=8,height=4.5,useDingbats=FALSE)
g = ggplot(filter(cluster_time_frequency),
 aes(x=timepoint, y=freq, group=new_celltype, fill=new_celltype)) +
  geom_bar(stat='identity', color="black") +
 facet_grid(overlap+S2 ~ Donor) + theme_bw() +
 xlab("Timepoint") + ylab("Frequency") +
 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
 labs(fill="Annotation") +
 geom_text(position = position_stack(vjust = 0.5), size=3, aes(label=label)) +
 scale_fill_manual(values=palette[level.order]) +
 theme(strip.background =element_rect(fill="white"))
print(g)
dev.off()

# Check SHM of C5.2 - RM B cells at 6mo
mtime$period = "Early"
mtime$period[mtime$timepoint == "M6"] = "M6"
mufreq_time = mtime %>%
  group_by(Donor, S2, clone_id, period) %>%
  summarize(mu_freq = mean(mu_freq)) %>%
  spread(period, mu_freq) %>%
  filter(!is.na(`M6`)& !is.na(Early))

# Compare mutation frequency across donors, tissues, and cell types
boxsize=0.15
bracketsize=0.15
axissize=0.15
panelsize=0.25
sigsize=2
pointsize=0.1
fontsize=9
my_comparisons <- list(c("Early","M6"))
max_shm = max(mtime$mu_freq)
pdf("results/shm_time.pdf",width=3.5,height=3.5,useDingbats=FALSE)
print(ggpaired(filter(mufreq_time),
    cond1 = "Early", cond2 = "M6",
          outlier.shape=NA,line.size=0.05) +
  stat_compare_means(comparisons = my_comparisons,
   method="wilcox",method.args=list(correct=TRUE),
   size=3,paired=TRUE,linecolor="grey") +
  facet_grid(S2~Donor)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))+
  theme(legend.position="none") +
  ylab("Mean SHM/clone")+
  xlab("") +
  theme_bw() +
  theme(strip.background =element_rect(fill="white")))
dev.off()

# tally number of clones for each comparison
clone_count = mufreq_time %>% 
  group_by(S2, Donor) %>%
  summarize(nclones=n_distinct(clone_id))

write.csv(clone_count, file="results/shm_time_n.csv", 
  quote=FALSE, row.names=FALSE)

# Now check SHM of clones at early vs M6 by cluster
# Check SHM of C5.2 - RM B cells at 6mo
mtime$period_celltype = as.character(mtime$new_celltype)
mtime$period_celltype[mtime$period=="Early"] = "All"
mufreq_time_agg = mtime %>%
  filter(!is.na(period_celltype)) %>%
  group_by(S2, clone_id, period, period_celltype) %>%
  summarize(mu_freq = mean(mu_freq))

agg = tibble()
for(celltype in unique(mtime$new_celltype)){
  if(is.na(celltype)){next}
  m6clones = filter(mufreq_time_agg, period_celltype == celltype)$clone_id
  temp = filter(mufreq_time_agg, clone_id %in% m6clones)
  temp = filter(temp, period == "Early" | period_celltype == celltype)
  temp$celltype = celltype
  agg = bind_rows(agg, temp)
}

agg_spread = agg %>%
  select(-period_celltype) %>%
  group_by(S2, clone_id, celltype) %>%
  spread(period, mu_freq)

agg_spread$celltype = factor(agg_spread$celltype, levels=level.order)

# Compare mutation frequency across donors, tissues, and cell types
boxsize=0.15
bracketsize=0.15
axissize=0.15
panelsize=0.25
sigsize=1
pointsize=0.1
fontsize=7
my_comparisons <- list(c("Early","M6"))
max_shm = max(agg$mu_freq)
pdf("results/shm_time_celltype.pdf",width=3.5,height=7,useDingbats=FALSE)
print(ggpaired(filter(agg_spread),
    cond1 = "Early", cond2 = "M6",
          line.size=0.01) +
  stat_compare_means(comparisons = my_comparisons,
   method="wilcox",method.args=list(correct=TRUE),
   size=3,paired=TRUE,linecolor="grey") +
  facet_grid(celltype~S2)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)))+
  theme(legend.position="none") +
  ylab("Mean SHM/clone")+
  xlab("") +
  theme_bw() +
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text.y = element_text(size = 6)))
dev.off()

agg_spread$celltype = factor(agg_spread$celltype)
agg_spread_count = agg_spread %>%
  group_by(S2, celltype) %>% 
  summarize(nclones = n_distinct(clone_id))

write.csv(agg_spread_count, file="results/shm_time_celltype_n.csv", 
  quote=FALSE, row.names=FALSE)

# multi timepoint clones with defined celltypes
# tally frequency of mutli timepoint clones in S2+ vs S2-
multi_time_s2 = comb_germline %>%
  filter(!is.na(new_celltype)) %>%
  filter(!naive) %>%
  group_by(Donor, S2, clone_id) %>%
  summarize(ntime=n_distinct(timepoint),
    multi_time = ntime > 1, ncells = n_distinct(new_barcode)) %>%
  group_by(Donor, S2, multi_time) %>%
  summarize(clones = n_distinct(clone_id)) %>%
  mutate(freq = signif(clones/sum(clones), digits=2)) %>%
  mutate(total = sum(clones)) %>%
  filter(multi_time) %>%
  select(-multi_time) %>%
  mutate(multi_timepoint = clones) %>%
  select(Donor, S2, multi_timepoint, total, freq)

write.csv(multi_time_s2, file="results/multi_timepoint_clones.csv", 
  quote=FALSE, row.names=FALSE)

# Make upset plots of time-point clones
input2 = comb_germline %>%
  filter(!naive & locus=="IGH") %>%
  group_by(Donor, S2, clone_id, timepoint) %>%
  summarize(ncells = 1) %>%
  mutate(ntimepoints = n_distinct(timepoint)) %>%
  spread(timepoint, ncells, fill=0)

multi = filter(input2, ntimepoints > 1)

times = names(multi[,5:10])

plots = list()
for(donor in unique(input2$Donor)){
  for(s in c("S-2P-", "S-2P+")){
    temp = filter(multi, Donor == donor, S2 == s)
    pdf(paste0("results/",donor,"_",s,".pdf"), width=3, height=2.5)
    print(upset(data.frame(temp[,5:10]),  keep.order = T,
      sets = times, mb.ratio = c(0.6, 0.4),
      nsets=100, nintersects=10, text.scale = 1, 
       mainbar.y.label = paste(donor, s)))
    dev.off()
  }
}

# Compare convergent sequences celltypes
comb_germline$convergent = "Not convergent"
comb_germline$convergent[comb_germline$convergent_cluster_donors > 1] = "Convergent"

# set order of clusters so C5 subclusters are on the bottom
level.order2 = c("C1 - Naive B cell 1", "C2 - Naive B cell 2", "C3 - unswitched MBC",
 "C4 - Memory B cell 1", "C6 - Plasmablast", "C7 - Naive B cell 3", "C5.1 - Mixed MBC", 
 "C5.2 - RM", "C5.3 - CD38+ AM", "C5.4 - TLM", "C5.5 - CD38- AM")
level.order2 = new_names[level.order2]

convergent_counts = comb_germline %>%
  filter(!is.na(new_celltype)) %>%
  group_by(Donor, group, convergent, new_celltype) %>%
  summarize(cells=n_distinct(new_barcode)) %>%
  mutate(Freq = cells/sum(cells))
convergent_counts$label = convergent_counts$cells
convergent_counts$label[convergent_counts$Freq < 0.15] = ""
convergent_counts$new_celltype = factor(convergent_counts$new_celltype, 
  levels=level.order2)

g = ggplot(filter(convergent_counts,  group == "S-2P+"),
 aes(x=convergent, y=Freq, group=new_celltype, fill=new_celltype)) +
  geom_bar(stat='identity', color="black") +
 facet_grid(. ~ Donor) + theme_bw() +
 xlab("") + ylab("Frequency of convergent clusters") +
 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
 labs(fill="Annotation") +
 geom_text(position = position_stack(vjust = 0.5), size=3, aes(label=label)) +
 scale_fill_manual(values=palette[level.order2]) +
 theme(strip.background =element_rect(fill="white"))
pdf("results/convergent_celltypes.pdf", width=4.75,height=4)
print(g)
dev.off()

write.csv(convergent_counts, file="results/convergent_counts_celltype.csv",
  row.names=FALSE,quote=FALSE)

# are S2P+ convergent cells enriched for cluster 5?
contingency = convergent_counts %>%
  group_by(Donor) %>%
  filter(group == "S-2P+") %>%
  mutate(C5=grepl("SC5", new_celltype)) %>%
  select(Donor, convergent, cells, C5) %>%
  group_by(convergent, C5) %>%
  summarize(cells=sum(cells)) %>%
  tidyr::spread(convergent, cells)

chisq.test(contingency[,2:3])

# number of convergent clusters in each sort
print(comb_germline %>%
  group_by(group, convergent) %>%
  summarize(n=n_distinct(convergent_cluster)) %>%
  mutate(Freq = n/sum(n)*100))

print(comb_germline %>%
  group_by(convergent) %>%
  summarize(n_distinct(convergent_cluster)))

# Compare convergent sequences celltypes
convergent_counts = comb_germline %>%
  filter(!is.na(new_celltype)) %>%
  group_by(Donor, group, public, new_celltype) %>%
  summarize(cells=n_distinct(new_barcode)) %>%
  mutate(Freq = cells/sum(cells))
convergent_counts$label = convergent_counts$cells
convergent_counts$label[convergent_counts$Freq < 0.15] = ""

convergent_counts$public_antibody = "   Not public"
convergent_counts$public_antibody[convergent_counts$public] = "      Public"
convergent_counts$new_celltype = factor(convergent_counts$new_celltype, 
  levels=level.order2)
names(convergent_counts$new_celltype) = NULL

g = ggplot(filter(convergent_counts, group == "S-2P+"),
 aes(x=public_antibody, y=Freq, group=new_celltype, fill=new_celltype)) +
  geom_bar(stat='identity', color="black") +
 facet_grid(. ~ Donor) + theme_bw() +
 xlab("") + ylab("Frequency of public antibodies") +
 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
 labs(fill="Annotation") +
 geom_text(position = position_stack(vjust = 0.5), size=3, aes(label=label)) +
 scale_fill_manual(values=palette[level.order2]) +
 theme(strip.background =element_rect(fill="white"))
pdf("results/public_antibodies_celltype.pdf", width=4.75,height=4)
print(g)
dev.off()

write.csv(convergent_counts, file="results/convergent_counts_celltype.csv",
  row.names=FALSE,quote=FALSE)

# are S2P+ convergent cells enriched for cluster 5?
contingency = convergent_counts %>%
  group_by(Donor) %>%
  filter(group == "S-2P+") %>%
  mutate(C5=grepl("SC5", new_celltype)) %>%
  select(Donor, public_antibody, cells, C5) %>%
  group_by(public_antibody, C5) %>%
  summarize(cells=sum(cells)) %>%
  tidyr::spread(public_antibody, cells)

print(chisq.test(contingency[,2:3]))

covid = read.csv("CoV-AbDab_031022.csv")
covid = filter(covid, grepl("uman",Heavy.V.Gene))
names(covid)[1] = "Name"

# Gene usage plots
plots = list()
for(locus in unique(comb_germline$locus)){
  genes <- countGenes(filter(comb_germline, locus==!!locus), gene="v_call", 
    groups=c("locus","group","convergent"), mode="gene", fill=TRUE)

  ctemp = covid
  if(locus == "IGH"){
    clocus = substr(covid$Heavy.V.Gene,1,3)
    ctemp = ctemp[clocus == locus, ]
    ctemp$v_call = getGene(ctemp$Heavy.V.Gene)
  }else{
    clocus = substr(covid$Light.V.Gene,1,3)
    ctemp = ctemp[clocus == locus, ]
    ctemp$v_call = getGene(ctemp$Light.V.Gene)
  }
  cgenes = countGenes(ctemp, gene="v_call", fill=TRUE, mode="gene")
  cgenes$convergent = "CoV-AbDab"
  cgenes$locus = locus

  cgenes$group = "CoV-AbDab"
  genes = bind_rows(genes, cgenes)
  
  genes$group = factor(genes$group, levels=c("S-2P-","S-2P+","Plasmablast","CoV-AbDab"))
  
  label_names = as_labeller(c("S-2P-" = "S-2P-", "S-2P+" = "S-2P+", 
    "Plasmablast" = "PB","CoV-AbDab"="CoV-Ab"))

  discard = genes %>%
    group_by(gene) %>%
    summarize(max_freq = max(seq_freq)) %>%
    filter(max_freq < 0.05)

  colors = c("Not convergent"="#33a02c", "Convergent"="#ff7f00","CoV-AbDab"="#6a3d9a")
  genes$convergent = factor(genes$convergent, levels=names(colors))
  g = ggplot(filter(genes, !gene %in% discard$gene), aes(x=gene, y=seq_freq, 
    fill=convergent)) + 
  geom_bar(position="dodge", stat="identity") + 
  facet_grid(group ~ ., labeller=label_names) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Frequency") + xlab("") +
  labs(fill="") +
  scale_fill_manual(values=colors) +
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.key.size = unit(0.5, 'cm'),
    legend.title = element_text(size=8),
    legend.text = element_text(size=8))

  plots[[locus]] = g
}

pdf("results/gene_usage.pdf", width=6,height=10,useDingbats=FALSE)
gridExtra::grid.arrange(grobs=plots,nrow=3)
dev.off()

