# Kenneth B. Hoehn
# 11/16/22
# Tree analysis script for moderna vaccine project
# Build trees using Dowser and IgPhyML
# Perform measurable evolution test
# Perform SP test
# Make figures for tree-based analyses

library(alakazam)
library(scoper)
library(dplyr)
library(Seurat)
library(shazam)
library(dowser)
library(ggpubr)
library(ggtree)
library(tidyr)

# updated dowser::plotTrees function, common_scale bug fixed in development version
plotTrees <- function(trees, nodes=FALSE, tips=NULL, tipsize=NULL, 
    scale=0.01, node_palette="Dark2", tip_palette=node_palette, base=FALSE,
    layout="rectangular", node_nums=FALSE, tip_nums=FALSE, title=TRUE,
    labelsize=NULL, common_scale=FALSE, ambig="blend"){

    tiptype = "character"
    if(!base){
        cols <- c()
        # set up global tip and node palette
        if(!is.null(tips) && nodes && sum(tip_palette != node_palette) == 0){
            tipstates <- unique(c(unlist(lapply(trees$data,function(x)
                unique(x@data[[tips]])))))
            if(is.numeric(tipstates)){
                stop("Can't currently plot numeric tip values and node values")
            }
            tipstates = c(sort(tipstates),"Germline")
            nodestates <- sort(unique(unlist(lapply(trees$trees,function(x)
                    unique(unlist(strsplit(x$state,split=",")))
                    ))))
            combpalette <- getPalette(c(nodestates,tipstates),node_palette)
            trees$trees <- colorTrees(trees$trees,palette=combpalette,ambig=ambig)
            nodestates <- unlist(lapply(trees$trees,function(x){
                colors <- x$node.color
                names(colors) <- x$state
                colors
                }))
            nodepalette <- nodestates[unique(names(nodestates))]
            cols <- c(combpalette,nodepalette[!names(nodepalette) %in% names(combpalette)])
        }else{
            # set up global tip palette
            if(!is.null(tips)){
                tipstates <- unique(c(unlist(lapply(trees$data,function(x)
                    unique(x@data[[tips]])))))
                if(is.numeric(tipstates)){
                    tiptype <- "numeric"
                    cols <- range(tipstates)
                }else{
                    tipstates = c(sort(tipstates),"Germline")
                    if(is.null(names(tip_palette))){
                        tip_palette <- getPalette(tipstates,tip_palette)
                        tip_palette <- tip_palette[!is.na(names(tip_palette))]
                    }else{
                        nfound <- tipstates[!tipstates %in% names(tip_palette)]
                        if(length(nfound) > 0){
                            stop(paste(nfound,"not found in tip_palette"))
                        }
                    }
                    cols <- tip_palette
                }
            }
            # set up global node palette
            if(nodes){
                if(is.null(names(node_palette))){
                    nodestates <- unique(unlist(lapply(trees$trees,function(x)
                        unique(unlist(strsplit(x$state,split=",")))
                        )))
                    statepalette <- getPalette(sort(nodestates),node_palette)
                    statepalette <- statepalette[!is.na(names(statepalette))]
                }else{
                    statepalette <- node_palette
                }
                trees$trees <- colorTrees(trees$trees,palette=statepalette, ambig=ambig)
                
                nodestates <- unlist(lapply(trees$trees,function(x){
                    colors <- x$node.color
                    names(colors) <- x$state
                    colors
                    }))
                nodepalette <- nodestates[unique(names(nodestates))]
                cols <- c(tip_palette,nodestates)
            }
        }
        if(common_scale){
            # get maximum divergence value
            max_div <- max(unlist(lapply(trees$trees, function(x)max(getDivergence(x)))))
        }
        
        ps <- lapply(1:nrow(trees),function(x)plotTrees(trees[x,],
            nodes=nodes,tips=tips,tipsize=tipsize,scale=scale,node_palette=node_palette,
            tip_palette=tip_palette,base=TRUE,layout=layout,node_nums=node_nums,
            tip_nums=tip_nums,title=title,labelsize=labelsize, ambig=ambig))
        if(!is.null(tips) || nodes){
            ps  <- lapply(ps,function(x){
                    x <- x + theme(legend.position="right",
                    legend.box.margin=margin(0, -10, 0, 0))+
                    guides(color=guide_legend(title="State"))
                    if(tiptype == "character"){
                        x <- x + scale_color_manual(values=cols)
                    }else{
                        x <- x + scale_color_distiller(limits=cols,
                            palette=tip_palette)
                    }})
        }
        if(common_scale){
             ps  <- lapply(ps,function(x){
                x <- x + xlim(0, max_div*1.05)
            })
        }
        return(ps)
    }

    tree <- trees$trees[[1]]
    data <- trees$data[[1]]
    p <- ggtree::ggtree(tree,layout=layout)
    if(!is.null(data)){
        if(!is(data,"list")){
            data <- list(data)
        }
        index <- which(unlist(lapply(data,function(x)x@clone == tree$name)))
        if(length(index) == 0){
            stop("clone",tree$name," not found in list of clone objects")
        }
        if(length(index) > 1){
            stop("clone",tree$name," found more than once in list of clone objects")
        }
        data <- data[[index]]
        gl <- dplyr::tibble(sequence_id="Germline")
        for(n in names(data@data)){
            if(is(data@data[[n]],"numeric") || is(data@data[[n]], "integer")){
                gl[[n]] <- NA
            }else if(is(data@data[[n]], "character")){
                gl[[n]] <- "Germline"
            }else{
                gl[[n]] <- "Germline"
            }
        }
        data@data <- rbind(data@data,gl)
        p <- p %<+% data@data
    }
    if(!is.null(tree$pars_recon)){
        if(nodes){
            p <- p + aes(color=tree$state)
        }
    }
    if(!is.null(tips)){
        if(is.null(data)){
            stop("dataframe must be provided when tip trait specified")
        }
        if(!is.null(tipsize)){
            if(is(tipsize, "numeric")){
                p <- p + ggtree::geom_tippoint(aes(color=!!rlang::sym(tips)),size=tipsize)
            }else if(is(tipsize, "character")){
                p <- p + ggtree::geom_tippoint(aes(color=!!rlang::sym(tips),
                    size=!!rlang::sym(tipsize)))
            }
        }else{
            p <- p + ggtree::geom_tippoint(aes(color=!!rlang::sym(tips)))
        }
    }
    if(scale != FALSE){
        p <- p + ggtree::geom_treescale(width=scale)
    }
    if(title){
        p <- p + ggtitle(data@clone)
    }
    if(node_nums){
        if(is.null(labelsize)){
            p <- p + ggtree::geom_label(data=p$data[!p$data$isTip,],
                aes(label=!!rlang::sym("node")),label.padding = unit(0.1, "lines"),
                label.size=0.1)
        }else{
            p <- p + ggtree::geom_label(data=p$data[!p$data$isTip,],
                aes(label=!!rlang::sym("node")),label.padding = unit(0.1, "lines"),
                label.size=0.1,size=labelsize)
        }
    }
    if(tip_nums){
        if(is.null(labelsize)){
            p <- p + ggtree::geom_label(data=p$data[p$data$isTip,],
                aes(label=!!rlang::sym("node")),label.padding = unit(0.1, "lines"),
                label.size=0.1)
        }else{
            p <- p + ggtree::geom_label(data=p$data[p$data$isTip,],
                aes(label=!!rlang::sym("node")),label.padding = unit(0.1, "lines"),
                label.size=0.1,size=labelsize)
        }
    }
    p
}

# new cluster names
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

# set order of clusters
level.order = c("C1 - Naive B cell 1", "C2 - Naive B cell 2", "C3 - unswitched MBC",
 "C4 - Memory B cell 1", "C5.1 - Mixed MBC", "C5.2 - RM",
  "C5.3 - CD38+ AM", "C5.4 - TLM", "C5.5 - CD38- AM", "C6 - Plasmablast", "C7 - Naive B cell 3")
level.order = new_names[level.order]

# cluster info and colors
cinfo = read.csv("cluster_key.csv")
palette = cinfo$Color
names(palette) = new_names[cinfo$New]

# timepoint colors
timepoint_palette = c(
  "v1D14"="#313695",
  "v2D6"="#4575b4",
  "v2D9"="#abd9e9",
  "v2D14"="#ffffbf",
  "v2D28"="#f46d43",
  "M6"="#a50026")

# timepoint colors
time_palette = c(
  "14"="#313695",
  "34"="#4575b4",
  "37"="#abd9e9",
  "42"="#ffffbf",
  "56"="#f46d43",
  "180"="#a50026")

sessionInfo()

igphyml="/ysm-gpfs/home/kbh35/project/programs/igphyml/src/igphyml"
nproc=10

seq_cutoff = 3

comb_germline = readRDS("processed/all_public_analyzed_cloned_data.rds")

if(sum(is.na(comb_germline$time)) != 0){
  stop("NA time values")
}
print(mean(is.na(comb_germline$new_celltype)))
comb_germline$new_celltype = as.character(comb_germline$new_celltype)
comb_germline$new_celltype[is.na(comb_germline$new_celltype)] = "Missing"

# Remove sequences with PTCs
if(FALSE){
print("removing ptcs")
  comb_germline$aa = translateDNA(comb_germline$sequence_alignment)
  remove = grepl("\\*",comb_germline$aa)
  print(paste("removing",sum(remove),"sequences"))
  comb_germline = comb_germline[!remove,]
  saveRDS(comb_germline, "processed/all_public_analyzed_cloned_data_noPTC.rds")
}
comb_germline = readRDS("processed/all_public_analyzed_cloned_data_noPTC.rds")

hc = unique(filter(comb_germline,locus=="IGH")$new_barcode)
comb_germline = filter(comb_germline, new_barcode %in% hc)

# build phylogenetic trees, test for measurable evolution
if(FALSE){
  f = formatClones(filter(comb_germline,locus=="IGH"), 
      columns=c("Donor", "s2clone", "timepoint","modality"), 
      traits=c("time", "timepoint", "sort", "new_celltype", "modality", 
        "convergent_cluster_donors", "public"))

  trees = tibble()
  for(i in unique(f$Donor)){
      cf = filter(f, Donor == i & grepl(",",timepoint) & seqs >= seq_cutoff)
      tr = getTrees(cf, build="igphyml", partition="cf", nproc=nproc,
          exec=igphyml)
      trees = bind_rows(trees, tr)
  }   
  ot = trees

  trees = correlationTest(filter(trees,grepl(",",timepoint)), time="time", 
      permutations=1000000, nproc=nproc, perm_type="uniform", polyresolve=FALSE)
  trees = trees[order(trees$p),]

  saveRDS(trees, "results/trees_new2.rds")
}else{
  print("skipping tree building")
}

# SP test with all timepoints and cell types
if(FALSE){
  f = formatClones(filter(comb_germline, locus=="IGH" 
      & new_celltype != "Missing"), 
      columns=c("Donor", "s2clone", "timepoint", "new_celltype"), 
      traits=c("new_celltype"))

  results_list = list()
  for(i in unique(f$Donor)){
      print(i)
      cf = filter(f, Donor == i & grepl(",",new_celltype) & s2clone)

      results = findSwitches(cf,trait="new_celltype",
        id=paste0(i,"allmem"),igphyml=igphyml,dir="temp",
        permutations=1000,nproc=nproc,quiet=0,resolve=2,fixtrees=FALSE,
        downsample=TRUE,tip_switch=20,force_resolve=TRUE, rm_temp=TRUE)

      results_list[[i]] = results
  }
  saveRDS(results_list, "results/switches_all_bootstrap_new2.rds")
}

# SP test from earlier timepoints to M6 MBC-SC5.2
if(FALSE){
  comb_germline$new_celltype2 = comb_germline$new_celltype
  comb_germline$new_celltype2[comb_germline$new_celltype == "MBC-SC5.2" & 
    comb_germline$timepoint == "M6"] = "Rmem-6m"
  temp = filter(comb_germline, timepoint != "M6" | new_celltype2 == "Rmem-6m")

  f = formatClones(filter(temp, locus=="IGH" & new_celltype != "Missing" & modality != "bulk"), 
      columns=c("Donor", "s2clone", "timepoint", "new_celltype"), 
      traits=c("new_celltype2"))

  results_list = list()
  for(i in unique(f$Donor)){
      print(i)
      cf = filter(f, Donor == i & grepl(",",new_celltype) & s2clone)

      results = findSwitches(cf,trait="new_celltype2",
        id=paste0(i,"rmem6"),igphyml=igphyml,dir="temp2",
        permutations=1000,nproc=nproc,quiet=0,resolve=2,fixtrees=FALSE,
        downsample=TRUE,tip_switch=20,force_resolve=TRUE, rm_temp=TRUE)

      results_list[[i]] = results
  }   
  saveRDS(results_list, "results/switches_M6_bootstrap_new2.rds")
}

# SP test results considering M6 MBC-SC5.2 B cells
results_list = readRDS("results/switches_M6_bootstrap_new2.rds")
switching_6m = bind_rows(lapply(results_list, function(x)x$switches))

# SP test results considering all memory B cells
results_list = readRDS("results/switches_all_bootstrap_new2.rds")
switching_mem = bind_rows(lapply(results_list, function(x)x$switches))

res = tibble()
reps = tibble()
for(id in unique(switching_6m$ID)){
  temp = filter(switching_6m, ID==id)
  temp = tryCatch(testSP(temp, to="Rmem-6m",permuteAll=TRUE),error=function(e)e)
  if("error" %in% class(temp)){next}
  temp$means$ID = id
  temp$reps$ID = id
  temp$ID = id
  if(length(temp) == 0){next}
  res = bind_rows(res, temp$means)
  print(temp$means)
}
filter(res, PGT < 0.1)

res$text = ""
res$text[res$PGT < 0.1] = ""
res$text[res$PGT < 0.05] = "X"
res$ID = substr(res$ID, 1, 6)
res$FROM = gsub("-MBC"," MBC", res$FROM)
res$FROM = factor(res$FROM, levels=level.order)
pdf("results/SPtest_M6_RM.pdf",width=3.25,height=4,useDingbats=FALSE)
ggplot(res, aes(y=FROM, x=ID, fill=PGT, 
  label=text)) + geom_tile() +
  scale_fill_gradient(low="red",high="white",limits=c(0,1)) +
  geom_text(color="white", size=5) +
  theme_bw() + 
  ylab(paste0("Cell type to MBC-SC5.2, M6")) +
  xlab("Subject") +
  guides(fill = guide_legend(title = "SP test\nP value")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

res = tibble()
reps = tibble()
for(id in unique(switching_mem$ID)){
  temp = filter(switching_mem, ID==id)
  temp = testSP(temp, to="MBC-SC5.2",permuteAll=TRUE)
  temp$means$ID = id
  temp$reps$ID = id
  if(length(temp) == 0){next}
  res = bind_rows(res, temp$means)
  reps = bind_rows(reps, temp$reps)
  print(temp$means)
}
filter(res, PGT < 0.1)

res$text = ""
res$text[res$PGT < 0.1] = ""
res$text[res$PGT < 0.05] = "X"
res$ID = substr(res$ID, 1, 6)
res$FROM = gsub("-MBC"," MBC", res$FROM)
res$FROM = factor(res$FROM, levels=level.order)
pdf("results/SPtest_all_to_RM.pdf",width=3.5,height=4,useDingbats=FALSE)
ggplot(res, aes(y=FROM, x=ID, fill=PGT, 
  label=text)) + geom_tile() +
  scale_fill_gradient(low="red",high="white",limits=c(0,1)) +
  geom_text(color="white", size=5) +
  theme_bw() + 
  ylab(paste0("Cell type to MBC-SC5.2")) +
  xlab("Subject") +
  guides(fill = guide_legend(title = "SP test\nP value")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# Make tree plots
trees = readRDS("results/trees_new2.rds")
trees = trees[order(trees$seqs, decreasing=TRUE),]
trees$public = unlist(lapply(trees$data, function(x)sum(x@data$public) > 0))
trees$convergent = unlist(lapply(trees$data, function(x)max(x@data$convergent_cluster_donors)))

table(comb_germline$new_celltype)
table(comb_germline$timepoint)

# Plot all trees
shapes = c("S2pos"=21, "S2neg"=22, "Plasmablast"=23)
p = plotTrees(trees,scale=0.1)
for(i in 1:length(p)){
  p[[i]] = p[[i]] + 
    geom_tippoint(aes(fill=timepoint,shape=sort),stroke=0.4,size=2) + 
    scale_shape_manual(values=shapes, name="Annotation") +
    scale_fill_manual(values=timepoint_palette) + 
    guides(fill = guide_legend(override.aes = list(pch = 21) )) +
    geom_tiplab(aes(label=new_celltype,subset=new_celltype!="Missing")) +
    ggtitle(paste(trees$clone_id[[i]],"p =",signif(trees$p[[i]],digits=2),
                  "S2pos:",trees$s2clone[[i]], "public:",trees$public[[i]],
                  "donors:",trees$convergent[[i]]))+
    scale_x_continuous(expand=c(0.15,0))+
    scale_y_continuous(expand=c(0.05,0))+
    geom_tippoint(aes(subset=sort=="Germline"))
}
treesToPDF(p, file="results/all_trees.pdf",ncol=1,nrow=1)

# Plot all trees with MBC-SC5.2 at 6 months
m6memtrees = unlist(lapply(trees$data, function(x)
  sum(x@data$timepoint == "M6" & x@data$new_celltype=="MBC-SC5.2")))
m6trees = trees[m6memtrees > 0,]
p = plotTrees(m6trees,scale=0.1,common_scale=TRUE)
for(i in 1:length(p)){
  p[[i]] = p[[i]] + 
    geom_tippoint(aes(fill=timepoint,shape=sort),stroke=0.4,size=2) + 
    scale_shape_manual(values=shapes, name="Annotation") +
    scale_fill_manual(values=timepoint_palette) + 
    guides(fill = guide_legend(override.aes = list(pch = 21) )) +
    geom_tiplab(aes(label=new_celltype,subset=new_celltype!="Missing")) +
    ggtitle(paste(m6trees$clone_id[[i]],"p =",signif(m6trees$p[[i]],digits=2),
                  "S2pos:",m6trees$s2clone[[i]], "public:",m6trees$public[[i]],
                  "donors:",m6trees$convergent[[i]]))+
    scale_x_continuous(expand=c(0.15,0))+
    scale_y_continuous(expand=c(0.05,0))+
    geom_tippoint(aes(subset=sort=="Germline"))
}
treesToPDF(p, file="results/M6_restingmem_trees.pdf",ncol=1,nrow=1)

# make pretty trees figure
p = plotTrees(m6trees[1:3,],scale=0.1,common_scale=TRUE)
for(i in 1:length(p)){
  p[[i]] = p[[i]] + 
    geom_tippoint(aes(fill=timepoint,shape=sort),stroke=0.4,size=3) + 
    scale_shape_manual(values=shapes, name="Annotation") +
    scale_fill_manual(values=timepoint_palette) + 
    guides(fill = guide_legend(override.aes = list(pch = 21) )) +
    geom_tiplab(aes(label=new_celltype,subset=new_celltype!="Missing"),offset=0.005) +
    ggtitle(paste(m6trees$clone_id[[i]],"p =",signif(m6trees$p[[i]],digits=2),
                  "S2pos:",m6trees$s2clone[[i]]))+
    scale_x_continuous(expand=c(0.4,0))+
    scale_y_continuous(expand=c(0.05,0))+
    geom_tippoint(aes(subset=sort=="Germline"))
}
pdf("results/mem_trees.pdf",width=14,height=3,useDingbats=FALSE)
gridExtra::grid.arrange(grobs=p, ncol=3)
dev.off()

# Compare evolving vs non-evolving clones
comb_germline$ME = "Not evolving"
comb_germline$ME[comb_germline$clone_id %in% filter(trees, p < 0.05)$clone_id] = "Evolving"
comb_germline$S2 = "S-2P-"
comb_germline$S2[comb_germline$s2clone] = "S-2P+"

comb_germline$ME = factor(comb_germline$ME, levels=c("Not evolving","Evolving"))

cluster_time_frequency = comb_germline %>%
  filter(!is.na(seurat_clusters) & !naive) %>%
  group_by(Donor, S2, ME, timepoint, new_celltype) %>%
  summarize(n=n_distinct(new_barcode)) %>%
  mutate(freq = n/sum(n))

cluster_time_frequency$timepoint = factor(cluster_time_frequency$timepoint, levels=names(timepoint_palette))
cluster_time_frequency$label = cluster_time_frequency$n
cluster_time_frequency$label[cluster_time_frequency$freq < 0.15] = ""
cluster_time_frequency$new_celltype = factor(cluster_time_frequency$new_celltype, levels=level.order)

pdf("results/s2evolving_vs_nonevolving_clones.pdf",width=8.5,height=4.5,useDingbats=FALSE)
g = ggplot(filter(cluster_time_frequency),
 aes(x=timepoint, y=freq, group=new_celltype, fill=new_celltype)) +
  geom_bar(stat='identity', color="black") +
 facet_grid(S2+ME ~ Donor) + theme_bw() +xlab("Timepoint") + ylab("Frequency") +
 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
 labs(fill="Annotation") +
 geom_text(position = position_stack(vjust = 0.5), size=3, aes(label=label)) +
 scale_fill_manual(values=palette[level.order])+
 theme(strip.background =element_rect(fill="white"))
print(g)
dev.off()

cluster_time_frequency = comb_germline %>%
  filter(!is.na(seurat_clusters) & !naive) %>%
  group_by(Donor, S2, ME, timepoint, new_celltype, c_call) %>%
  summarize(n=n_distinct(new_barcode)) %>%
  mutate(freq = n/sum(n))

print(filter(cluster_time_frequency, ME=="Evolving" & grepl("Plasma", new_celltype)))

# Plot all ME trees with scatterplots
mtrees = filter(trees, p < 0.05)
p = plotTrees(mtrees,scale="none",common_scale=TRUE)
for(i in 1:length(p)){
  p[[i]] = p[[i]] + 
    geom_tippoint(aes(fill=timepoint,shape=sort),stroke=0.4,size=2,hjust=0.1) + 
    scale_shape_manual(values=shapes, name="Annotation") +
    scale_fill_manual(values=timepoint_palette) + 
    guides(shape="none",fill="none") +
    geom_tiplab(aes(label=new_celltype,subset=new_celltype!="Missing"),size=2) +
    ggtitle(paste(mtrees$clone_id[[i]],"p =",signif(mtrees$p[[i]],digits=2),
                  "S2pos:",mtrees$s2clone[[i]], "public:",mtrees$public[[i]],
                  "donors:",mtrees$convergent[[i]]))+
    theme(text = element_text(size = 3))+
    theme(plot.title = element_text(size = 8)) +
    geom_treescale(linesize=0.25,fontsize=3,width=0.1) +
    geom_tippoint(aes(subset=sort=="Germline"))
}

plots = list()
counts = 1
for(i in 1:length(mtrees$trees)){
    tree = mtrees$trees[[i]]
    data = mtrees$data[[i]]@data

    l = summary(lm(divergence ~ time,data=data))$coefficients[,1]

    # plot tree
    plots[[counts]] = p[[i]]

    # plot divergence scatterplot
    plots[[counts+1]] = ggplot(data,aes(x=time,y=divergence,
        shape=sort,fill=timepoint))+
        geom_point(size=2) + 
        scale_shape_manual(values=shapes) +
        scale_fill_manual(values=timepoint_palette)+
        labs(color="Time")+theme_bw()+xlab("Day")+
        ylab("Divergence")+geom_abline(intercept=l[1],slope=l[2])+
        theme(text = element_text(size = 8))+
        guides(fill=FALSE,shape=FALSE)

    counts = counts + 2
}
pdf("results/me_1.pdf",width=7.5,height=6)
gridExtra::grid.arrange(grobs=plots[1:12],ncol=4)
dev.off()

pdf("results/me_2.pdf",width=7.5,height=6)
gridExtra::grid.arrange(grobs=plots[13:24],ncol=4)
dev.off()

# compare to bar plot, just to double check
check = filter(trees, p < 0.05) %>%
group_by(Donor) %>%
do(bind_rows(lapply(.$data, function(x)x@data)) %>%
 group_by(timepoint, new_celltype) %>%
  summarize(n=sum(collapse_count)) %>%
  filter(new_celltype != "Missing")) 
print(data.frame(check))

# Get sequences of biggest clone in m6 trees
top_clone = m6trees$clone_id[1]
top_seqs = filter(comb_germline, clone_id == top_clone)

#Get corresponding convergent sequence cluster
top_cluster = unique(filter(comb_germline, clone_id == top_clone)$convergent_cluster)
cluster_seqs = filter(comb_germline, convergent_cluster == top_cluster)
matches = top_seqs$matches

table(cluster_seqs$timepoint, cluster_seqs$Donor)

# plot convergent sequence clones
# AA palette for plotting
aapalette = c(
"H" = "#a6cee3",
"K" = "#a6cee3",
"R" = "#a6cee3",
"D" = "#fb9a99",
"E" = "#fb9a99",
"S" = "#1f78b4",
"T" = "#1f78b4",
"N" = "#1f78b4",
"Q" = "#1f78b4",
"A" = "white",
"V" = "white",
"L" = "white",
"I" = "white",
"M" = "white",
"F" = "#cab2d6",
"Y" = "#cab2d6",
"W" = "#cab2d6",
"P" = "#fb9a99",
"G" = "#fb9a99",
"C" = "#fdbf6f",
"B" = "grey",
"Z" = "grey",
"X" = "grey")

#space the time periods out
cluster_seqs$new_celltype[cluster_seqs$new_celltype == "Missing"] = ""
cluster_seqs$ttime = as.character(cluster_seqs$timepoint)
cluster_seqs$ttime[nchar(as.character(cluster_seqs$timepoint)) == 4] = 
 paste0(cluster_seqs$timepoint[nchar(as.character(cluster_seqs$timepoint)) == 4],"  ")
cluster_seqs$ttime[nchar(as.character(cluster_seqs$timepoint)) == 2] = 
 paste0(cluster_seqs$timepoint[nchar(as.character(cluster_seqs$timepoint)) == 2],"     ")
cluster_seqs$id = paste(as.numeric(factor(cluster_seqs$new_barcode)),
  cluster_seqs$new_celltype, cluster_seqs$ttime, cluster_seqs$Donor, sep="|")
fd = filter(cluster_seqs, locus=="IGH")
fd = fd %>% arrange(Donor, timepoint, new_celltype)

# germline AA
glines = table(substr(fd$germline_alignment,1,312))
germline = names(glines[which.max(glines)])
AA = strsplit(translateDNA(gsub("\\.","",germline)),split="")
heavy_aa = bind_rows(lapply(1:length(AA),function(x){
  tibble(sequence_id="GermlineV",position=1:length(AA[[x]]),AA=AA[[x]])
}))
# sequence AA
AA = strsplit(translateDNA(gsub("\\.","",fd$sequence_alignment)),split="")
heavy_aa = bind_rows(heavy_aa,lapply(1:length(AA),function(x){
  tibble(sequence_id=fd$id[x],position=1:length(AA[[x]]),AA=AA[[x]])
}))

# get AA for AbDab
# matches from AbDab
ids = strsplit(fd$matches,split=",")[[1]]
covid = read.csv("CoV-AbDab_031022.csv")
covid = filter(covid, grepl("uman",Heavy.V.Gene))
names(covid)[1] = "Name"
fc = filter(covid, Name %in% ids & VHorVHH != "ND")
AA = strsplit(fc$VHorVHH,split="")
heavy_aa = bind_rows(heavy_aa,lapply(1:length(AA),function(x){
  tibble(sequence_id=fc$Name[x],position=1:length(AA[[x]]),AA=AA[[x]])
}))

heavy_aa$sequence_id = factor(heavy_aa$sequence_id,
  levels=unique(heavy_aa$sequence_id))

pdf("results/convergent_clone_aa.pdf",width=9.5,height=5,useDingbats=FALSE)
print(
  ggplot(heavy_aa, aes(x=position,y=sequence_id, fill=AA, label=AA)) + geom_tile() + 
  geom_text(size=2) + 
  theme_bw() + theme(legend.position = "none") + xlab("Position") + ylab("") +
  scale_x_continuous(expand=c(0,0.1)) +
  theme(text = element_text(size = 9), axis.text =element_text(size = 9)) +
  scale_fill_manual(values=aapalette) +
  scale_y_discrete(limits=rev)+
  ggtitle(paste0("Convergent cluster ",unique(cluster_seqs$convergent_cluster),
    " heavy chain matches"))
  )
dev.off()
