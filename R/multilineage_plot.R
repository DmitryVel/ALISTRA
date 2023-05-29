monocle_theme_opts <- function()
{
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(size=0.25, color="black")) +
    theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
}

theme_opts <- function()
{
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(plot.title = element_blank()) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_blank()) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.text.x = element_blank()) +
    theme(axis.title.x = element_blank()) +
    theme(axis.text.y = element_blank()) +
    theme(axis.title.y = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(axis.line.y = element_line(size=1, color="black")) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
}

#' @export
compress_UMAP <- function(cds, lineage, start, N = 2500){
  cds_name = deparse(substitute(cds))
  input = paste0("sel.cells = ",cds_name,"@lineages$", lineage)
  eval(parse(text=input))
  RNA_umap = reducedDims(cds)$"UMAP"[sel.cells,]
  window = nrow(RNA_umap)/N
  step = ((nrow(RNA_umap)-window)/N)
  #use sliding window to compress expression values and pseudotime
  print(paste0("Window: ", window))
  print(paste0("Step: ", step))
  umap_X = SlidingWindow("mean", RNA_umap[,1], window, step)
  umap_Y = SlidingWindow("mean", RNA_umap[,2], window, step)
  umap = as.data.frame(cbind(umap_X, umap_Y))
  umap
}

#' @export
common <- function(df){
  names(which.max(table(df)))
}

compress_meta <- function(meta, N){
  window = nrow(meta)/N
  step = ((nrow(meta)-window)/N)
  cluster = meta$cluster
  age = meta$age
  cluster.comp = SlidingWindow(common, cluster, window, step)
  age.comp = SlidingWindow(common, age, window, step)
  meta.comp = cbind(age.comp, cluster.comp)
  meta.comp = as.data.frame(meta.comp)
  colnames(meta.comp) <- c("age", "cluster")
  meta.comp
}

#' @export
isolate_lineage_ATAC <- function(cds, ATAC, lineage, sel_clusters = NULL, starting_clusters = NULL, subset = FALSE, N = 5, cl = 1){
  cds_name = deparse(substitute(cds))
  input = paste0("sub.graph = ",cds_name,"@graphs$", lineage)
  eval(parse(text=input))
  nodes_UMAP = cds@principal_graph_aux[["UMAP"]]$dp_mst
  if(subset == F){
    nodes_UMAP.sub = as.data.frame(t(nodes_UMAP[,names(V(sub.graph))]))
  }
  else{
    g = principal_graph(cds)[["UMAP"]]
    dd = degree(g)
    names1 = names(dd[dd > 2 | dd == 1])
    names2 = names(dd[dd == 2])
    names2 = sample(names2, length(names2)/subset, replace = F)
    names = c(names1, names2)
    names = intersect(names(V(sub.graph)), names)
    nodes_UMAP.sub = as.data.frame(t(nodes_UMAP[,names]))
  }
  #select cells along the graph
  mean.dist = path.distance(nodes_UMAP.sub)
  r = mean.dist*N
  cells_UMAP = as.data.frame(ATAC.integrated[['UMAP']][[]])
  colnames(cells_UMAP) = c("UMAP_1", "UMAP_2")
  sel.cells = cell.selector(nodes_UMAP.sub, cells_UMAP, r, cl = cl)
  #only keep cells in the progenitor and lineage-specific clusters
  sel.cells1 = c()
  sel.cells2 = sel.cells
  if(length(starting_clusters) > 0){
    sel.cells1 = names(Idents(ATAC.integrated)[Idents(ATAC.integrated) %in% starting_clusters])
  }
  if(length(sel_clusters) > 0){
    sel.cells2 = names(Idents(ATAC.integrated)[Idents(ATAC.integrated) %in% sel_clusters])
  }
  cells = unique(c(sel.cells1, sel.cells2))
  sel.cells = sel.cells[sel.cells %in% cells]
  return(sel.cells)
}

#' @export
import_ATAC <- function(cds_ref, ATAC_counts, UMAP, ATAC_meta, lineage){
  cds_name = deparse(substitute(cds_ref))
  input = paste0("sub.graph = ",cds_name,"@graphs$", lineage)
  eval(parse(text=input))
  nodes_UMAP = cds_ref@principal_graph_aux[["UMAP"]]$dp_mst
  #create monocle object with ATAC counts or gene activity counts
  gene_annotation = as.data.frame(rownames(ATAC_counts))
  rownames(gene_annotation) = rownames(ATAC_counts)
  colnames(gene_annotation) = 'gene_short_name'
  cds = new_cell_data_set(ATAC_counts, cell_metadata = ATAC_meta, gene_metadata = gene_annotation)
  reducedDims(cds)$"UMAP" <- UMAP
  s.clusters = as.character(ATAC_meta$cluster)
  names(s.clusters) <- rownames(ATAC_meta)
  cds@clusters$"UMAP"$"clusters" <- s.clusters
  cds@clusters$UMAP$partitions <- cds@clusters$UMAP$clusters
  cds@clusters$UMAP$partitions[cds@clusters$UMAP$partitions != "1"] <- "1"
  principal_graph(cds)[["UMAP"]] <- sub.graph
  cds@principal_graph_aux[["UMAP"]]$dp_mst <- nodes_UMAP[,names(V(sub.graph))]
  cds = import_monocle(cds)
  ref.umap = reducedDims(cds_ref)$"UMAP"
  atac.umap = reducedDims(cds)$"UMAP"
  n_neighbors = 10
  nn <- get.knnx(ref.umap, atac.umap, k=n_neighbors)
  pt <- cds_ref@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
  pt_atac = apply(nn$"nn.index", 1, function(x) mean(pt[x]))
  cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]] <- pt_atac
  cds@principal_graph_aux[["UMAP"]]$pseudotime <- pt_atac
  names(cds@principal_graph_aux[["UMAP"]]$pseudotime) <- row.names(colData(cds))
  sel.cells = colnames(ATAC_counts)
  cds_name = "cds"
  input = paste0(cds_name, "@graphs$", lineage, " <- sub.graph")
  eval(parse(text=input))
  input = paste0(cds_name, "@lineages$", lineage, " <- sel.cells")
  eval(parse(text=input))
  cds
}


#' @export
compress <- function(df, window = 3, step){
df.comp = SlidingWindow("mean", df, window, step)
}

#' @export
fit.m3 <- function(exp.sel, pt, max.pt, model = "expression ~ splines::ns(pseudotime, df=3)", N = 500){
require(speedglm)
family = stats::quasipoisson()
exp_data.sel = cbind(pt, exp.sel)
colnames(exp_data.sel) <- c("pseudotime","expression")
exp_data.sel = as.data.frame(exp_data.sel)
exp_data.sel$pseudotime <- as.numeric(as.character(exp_data.sel$pseudotime))
exp_data.sel$expression <- as.numeric(as.character(exp_data.sel$expression))
tryCatch({fit = speedglm(model, data = exp_data.sel, family = family, acc=1e-3, model=FALSE, y=FALSE)
d = as.data.frame(seq(from=0, to=max.pt, by = max.pt/(N-1)))
colnames(d) <- c("pseudotime")
fit = stats::predict(fit, newdata=d, type="response")
return(fit)
}, error=function(cond) {return(rep("NA", N))})
}

#' @export
as_matrix <- function(mat){

  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
    
  for (i in seq_along(val)){
      tmp[row_pos[i],col_pos[i]] <- val[i]
  }
    
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

#' @export
compress2 <- function(df, window, step){
  df.comp = SlidingWindow("mean", df, window, step)
}

#' @export
pull <- function(df, window, step){
  df.comp = SlidingWindow("sum", df, window, step)
}

#' @export
compress_lineages_v2 <- function(cds, start, window = F, N = 500, cores = F){
  lineages = names(cds@lineages)
  for(lineage in lineages){
    print(lineage)
    cds = compress_lineage_v2(cds, lineage = lineage, start = start, window = window, gene = FALSE, N = N, cores = cores)
    gc()
  }
  return(cds)
}

#' @export
compress_lineage_v2 <- function(cds, lineage, start, window = F, gene = FALSE, N = 500, cores = F, cells = FALSE){
  cds_name = deparse(substitute(cds))
  if(gene == FALSE){
    input = paste0("compress_expression_v2(",cds_name,", lineage = '", lineage, "', start = ", start, ", window = ", window, ", gene = ", gene, ", N = ", N, ", cores = ", cores, ")")
  }
  else{
    input = paste0("compress_expression_v2(",cds_name,", lineage = '", lineage, "', start = ", start, ", window = ", window, ", gene = '", gene, "', N = ", N, ", cores = ", cores, ")")
  }
  exp = eval(parse(text=input))
  input = paste0(cds_name, "@expression$", lineage, " <- exp$expression")
  eval(parse(text=input))
  input = paste0(cds_name, "@expectation$", lineage, " <- exp$expectation")
  eval(parse(text=input))
  input = paste0(cds_name, "@pseudotime$", lineage, " <- exp$pseudotime")
  eval(parse(text=input))
  eval(parse(text=paste0("return(",cds_name, ")")))
}

#' @export
compress_expression_v2 <- function(cds, lineage, start, window = F, gene = FALSE, N = 500, cores = F){
  cds_name = deparse(substitute(cds))
  if(cores != F){
    cl <- makeCluster(cores)
    clusterEvalQ(cl, c(library(evobiR)))
  }
  input = paste0("get_lineage_object(",cds_name,", lineage = '", lineage, "', start = ", start, ")")
  cds_subset = eval(parse(text=input))
  family = stats::quasipoisson()
  model = "expression ~ splines::ns(pseudotime, df=3)"
  names(cds_subset) <- rowData(cds_subset)$gene_short_name
  exp = as_matrix(exprs(cds_subset))
  exp = t(exp) /  pData(cds_subset)[, 'Size_Factor']
  pt <- cds_subset@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
  pt = pt[order(pt)]
  exp = exp[names(pt),]
  if(window == FALSE){
      window = nrow(exp)/N
  }
  step = ((nrow(exp)-window)/N)
  #use sliding window to compress expression values and pseudotime
  print(paste0("Window: ", window))
  print(paste0("Step: ", step))
  pt.comp = SlidingWindow("mean", pt, window, step)
  max.pt = max(pt.comp)
  if(gene != F){
    exp.comp = compress2(exp[,gene], window = window, step = step)
  }
  else{
    print(paste0("Compressing lineage ", lineage, " and fitting curves"))
    step = ((nrow(exp)-window)/N)
    exp.comp = pbapply(exp, 2, compress2, window = window, step = step)
  }
  if(gene != F){
    exp_data.sel = cbind(pt.comp, exp.comp)
    exp_data.sel = as.data.frame(exp_data.sel)
    colnames(exp_data.sel) <- c("pseudotime", "expression")
    exp_data.sel$pseudotime <- as.numeric(as.character(exp_data.sel$pseudotime))
    exp_data.sel$expression <- as.numeric(as.character(exp_data.sel$expression))
    d = as.data.frame(seq(from=0, to=max.pt, by = max.pt/(N-1)))
    colnames(d) <- c("pseudotime")
    tryCatch({fit = speedglm(model, data = exp_data.sel, family = family, acc=1e-3, model=FALSE, y=FALSE)
    fit = predict(fit, newdata = d, type='response')
    }, error=function(cond) {fit = as.data.frame(rep(0, N))})
    exp = as.data.frame(cbind(exp.comp, fit, d))
    colnames(exp) <- c("expression", "expectation", "pseudotime")
    exp$expression[exp$expression < 0] <- 0
    exp$expectation[exp$expectation < 0] <- 0
  }
  else{
    d = as.data.frame(seq(from=0, to=max.pt, by = max.pt/(N-1)))
    fit = pbapply(exp.comp, 2, fit.m3, pt = d, max.pt = max(d), N = N)
    fit = apply(fit, 2, as.numeric)
    return(list("expression" = exp.comp, "expectation" = fit, "pseudotime" = d))
  }
  exp$expression[exp$expression < 0] <- 0
  exp$expectation[exp$expectation < 0] <- 0
  if(cores != F){
    stopCluster(cl)
  }
  return(exp)
}

#' @export
compress_lineages <- function(cds, start, N = 500, cores = F){
lineages = names(cds@lineages)
for(lineage in lineages){
print(lineage)
cds = compress_lineage(cds, lineage, start, gene = FALSE, N = N, cores = cores)
}
return(cds)
}

#' @export
compress_ATAC_expression <- function(cds, lineage, start, N = 500, cores = F){
  cds_name = deparse(substitute(cds))
  if(cores != F){
    cl <- makeCluster(cores)
    clusterEvalQ(cl, c(library(evobiR)))
  }
  input = paste0("get_lineage_object(",cds_name,", lineage = '", lineage, "', start = ", start, ")")
  cds_subset = eval(parse(text=input))
  family = stats::quasipoisson()
  model = "expression ~ splines::ns(pseudotime, df=3)"
  names(cds_subset) <- rowData(cds_subset)$gene_short_name
  exp = as_matrix(exprs(cds_subset))
  pt <- cds_subset@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
  pt = pt[order(pt)]
  exp = exp[,names(pt)]
  window = ncol(exp)/N
  step = ((ncol(exp)-window)/N)
  #use sliding window to compress expression values and pseudotime
  print(paste0("Window: ", window))
  print(paste0("Step: ", step))
  pt.comp = SlidingWindow("mean", pt, window, step)
  max.pt = max(pt.comp)
  print(paste0("Compressing lineage ", lineage, " and fitting curves"))
  step = ((ncol(exp)-window)/N)
  print("pulling counts into metacells")
  if(cores > 1){
    print("multicore processing")
    exp.comp = pbapply(exp, 1, pull, window = window, step = step, cl = cl)
  }
  else{
    exp.comp = pbapply(exp, 1, pull, window = window, step = step)
  }
  exp.comp = round(t(exp.comp))
  return(list(exp.comp, pt.comp))
}

#' @export
compress_lineage <- function(cds, lineage, start, gene = FALSE, N = 500, cores = F, cells = FALSE){
cds_name = deparse(substitute(cds))
if(gene == FALSE){
input = paste0("compress_expression(",cds_name,", lineage = '", lineage, "', start = ", start, ", gene = ", gene, ", N = ", N, ", cores = ", cores, ")")
}
else{
input = paste0("compress_expression(",cds_name,", lineage = '", lineage, "', start = ", start, ", gene = '", gene, "', N = ", N, ", cores = ", cores, ")")
}
exp = eval(parse(text=input))
input = paste0(cds_name, "@expression$", lineage, " <- exp$expression")
eval(parse(text=input))
input = paste0(cds_name, "@expectation$", lineage, " <- exp$expectation")
eval(parse(text=input))
input = paste0(cds_name, "@pseudotime$", lineage, " <- exp$pseudotime")
eval(parse(text=input))
eval(parse(text=paste0("return(",cds_name, ")")))
}

#' @export
compress_expression <- function(cds, lineage, start, window = 3, gene = FALSE, N = 500, cores = F){
cds_name = deparse(substitute(cds))
if(cores != F){
cl <- makeCluster(cores)
clusterEvalQ(cl, c(library(evobiR)))
}
input = paste0("get_lineage_object(",cds_name,", lineage = '", lineage, "', start = ", start, ")")
cds_subset = eval(parse(text=input))
family = stats::quasipoisson()
model = "expression ~ splines::ns(pseudotime, df=3)"
names(cds_subset) <- rowData(cds_subset)$gene_short_name
exp = as.data.frame(as_matrix(exprs(cds_subset)))
exp = t(t(exp) /  pData(cds_subset)[, 'Size_Factor'])
pt <- cds_subset@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
pt <- as.data.frame(pt)
colnames(pt) <- c("pseudotime")
#exp = cbind(pt, round(t(exp)))
exp = cbind(pt, t(exp))
exp = exp[order(exp$pseudotime),]
step = ((nrow(exp)-window)/N)
pt = exp[,"pseudotime"]
#use sliding window to compress expression values and pseudotime
pt.comp = SlidingWindow("mean", pt, window, step)
max.pt = max(pt.comp)
if(gene != F){
exp = exp[,c("pseudotime", gene)]
exp.comp = compress(exp[,gene], window, step = step)
}
else{
print(paste0("Compressing lineage ", lineage, " and fitting curves"))
mat <- as.data.frame(exp[,2:ncol(exp)])
if(cores != F){
exp.comp = pbsapply(mat, compress, step = step, cl = cl)
}
else{
exp.comp = pbsapply(mat, compress, step = step)
}
}
if(gene != F){
exp_data.sel = cbind(pt.comp, exp.comp)
exp_data.sel = as.data.frame(exp_data.sel)
colnames(exp_data.sel) <- c("pseudotime", "expression")
exp_data.sel$pseudotime <- as.numeric(as.character(exp_data.sel$pseudotime))
exp_data.sel$expression <- as.numeric(as.character(exp_data.sel$expression))
d = as.data.frame(seq(from=0, to=max.pt, by = max.pt/(N-1)))
colnames(d) <- c("pseudotime")
tryCatch({fit = speedglm(model, data = exp_data.sel, family = family, acc=1e-3, model=FALSE, y=FALSE)
fit = predict(fit, newdata = d, type='response')
}, error=function(cond) {fit = as.data.frame(rep(0, N))})
exp = as.data.frame(cbind(exp.comp, fit, d))
colnames(exp) <- c("expression", "expectation", "pseudotime")
exp$expression[exp$expression < 0] <- 0
exp$expectation[exp$expectation < 0] <- 0
}
else{
mat <- as.data.frame(exp.comp)
d = as.data.frame(seq(from=0, to=max.pt, by = max.pt/(N-1)))
if(cores != F){
fit = pbsapply(mat, fit.m3, pt = d, max.pt = max(d), N = N, cl = cl)
}
else{
fit = pbsapply(mat, fit.m3, pt = d, max.pt = max(d), N = N)
}
fit = apply(fit, 2, as.numeric)
return(list("expression" = exp.comp, "expectation" = fit, "pseudotime" = d))
}
exp$expression[exp$expression < 0] <- 0
exp$expectation[exp$expectation < 0] <- 0
if(cores != F){
stopCluster(cl)
}
return(exp)
}

#' @export
plot_ridge <- function(cds, gene, lineages, scale_factor = FALSE, alpha = 0.6, text.size = 18, plot.title.size = 24, legend.key.size = 2, legend.text.size = 10, colors = c("red", "blue", "green", "cyan", "magenta", "purple", "orange", "black", "yellow", "tan"), N = 500, legend_position = "right"){
  cds_name = deparse(substitute(cds))
  input = paste0(cds_name,"@expression$", lineages[1])
  M = nrow(eval(parse(text = input)))
  pts = c()
  for(lineage in lineages){
    input = paste0(cds_name,"@pseudotime$", lineage)
    pt = eval(parse(text = input))[,1]
    pts = c(pts, pt)
  }
  max.pt = max(pts)
  dd = data.frame("1"=0,"2"=0,"3"=0)[FALSE,]
  pt = seq(from=0, to=max.pt, by = max.pt/(M-1))
  fits = c()
  ys = c()
  i = 0
  for(lineage in lineages){
    input = paste0("exp = ",cds_name,"@expression$", lineage)
    eval(parse(text=input))
    if(gene %in% colnames(exp)){
      input = paste0("fit = ",cds_name,"@expectation$", lineage,"[,'",gene,"']")
      eval(parse(text=input))
    }
    else{
      fit = rep(0, N)
    }
    fits = c(fits, fit)
  }
  for(lineage in lineages){
    input = paste0("exp = ",cds_name,"@expression$", lineage)
    eval(parse(text=input))
    if(gene %in% colnames(exp)){
      input = paste0("fit = ",cds_name,"@expectation$", lineage,"[,'",gene,"']")
      eval(parse(text=input))
    }
    else{
      fit = rep(0, N)
    }
    dd = rbind(dd, cbind(fit,pt,rep(lineage, M)))
    if(scale_factor == FALSE){
    	ys = c(ys, rep(0,M))
    }
    else{
    	ys = c(ys, rep((max(fits)/scale_factor)*i,M))
    }
    i = i+1
    }
  colnames(dd) <- c("expression", "pseudotime", "lineage")
  dd$pseudotime <- as.numeric(dd$pseudotime)
  dd$expression <- as.numeric(dd$expression)
  ymax = max(fits)
  colors.adj = c()
  for(color in colors){
    colors.adj = append(colors.adj, adjustcolor(color, alpha.f = alpha))
  }
  q <- ggplot(dd, aes(pseudotime, ys, height = expression, group = lineage, fill = lineage), fit = lineage) + geom_ridgeline() + scale_fill_manual(values = colors.adj)
  #q <- q + scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))
  q <- q + scale_y_log10()
  q <- q + ylim(y = c(0,ymax))
  q <- q + monocle_theme_opts() + ylab("Expression") + xlab("Pseudotime") + ggtitle(gene) + theme(legend.key.size = unit(legend.key.size, 'cm'), plot.title = element_text(size = plot.title.size, face="bold", hjust = 0.5), axis.text=element_text(size=text.size), axis.title=element_blank(), legend.text=element_text(size=legend.text.size), legend.title=element_text(size=text.size, face = "bold"), legend.position = legend_position)
  q
}

#' @export
get_pt <- function(lineage){
load(file = paste0(lineage, ".R"))
pt <- cds_subset@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
pt = pt[order(pt)]
return(pt)
}

#' @export
common <- function(df){
  names(which.max(table(df)))
}

#' @export
compress_meta <- function(meta, N){
  window = nrow(meta)/N
  step = ((nrow(meta)-window)/N)
  cluster = meta$cluster
  age = meta$age
  cluster.comp = SlidingWindow(common, cluster, window, step)
  age.comp = SlidingWindow(common, age, window, step)
  meta.comp = cbind(age.comp, cluster.comp)
  meta.comp = as.data.frame(meta.comp)
  colnames(meta.comp) <- c("age", "cluster")
  meta.comp
}

#' @export
plot_multiple <- function(cds, gene, lineages, meta = NULL, points = T, age.scale = T, scale.lineage = NULL, age.points = c("3rd trimester", "0-1 years", "2-4 years", "4-10 years"), breaks.labels = c("2nd", "3rd", "birth", "1y", "4y"), point_size = 0.1, line_size = 1, text.size = 14, plot.title.size = 36, legend.key.size = 0.5, legend.text.size = 10, colors = c("red", "blue", "green", "cyan", "magenta", "purple", "orange", "black", "yellow", "tan"), N = 500, legend_position = "none"){
  cds_name = deparse(substitute(cds))
  input = paste0(cds_name,"@expression$", lineages[1])
  N = nrow(eval(parse(text = input)))
  pts = c()
  if(length(scale.lineage) == 0){
  for(lineage in lineages){
    input = paste0(cds_name,"@pseudotime$", lineage)
    pt = eval(parse(text = input))[,1]
    pts = c(pts, pt)
  }
  max.pt = max(pts)
  }
  else{
  input = paste0(cds_name,"@pseudotime$", scale.lineage)
  pt = eval(parse(text = input))[,1]
  max.pt = max(pt)
  }
  print(max.pt)
  if(points == T){
    dd = as.data.frame(seq(from=0, to=max.pt, by = max.pt/(N-1)))
    cols = c("pseudotime")
    fits = c()
    exps = c()
    for(lineage in lineages){
      input = paste0("exp = ",cds_name,"@expression$", lineage)
      eval(parse(text=input))
      if(gene %in% colnames(exp)){
        input = paste0("exp = ",cds_name,"@expression$", lineage,"[,'",gene,"']")
        eval(parse(text=input))
        input = paste0("fit = ",cds_name,"@expectation$", lineage,"[,'",gene,"']")
        eval(parse(text=input))
      }
      else{
        exp = rep(0, N)
        fit = rep(0, N)
      }
      dd = cbind(dd, exp, fit)
      cols = append(cols, paste0("exp_", lineage))
      cols = append(cols, paste0("fit_", lineage))
      fits = c(fits, fit)
      exps = c(exps, exp)
    }
    colnames(dd) <- cols
    ymax = max(fits)
  }
  else{
    fits = c()
    dd = matrix(ncol = 3, nrow = 0,)
    for(lineage in lineages){
      input = paste0("exp = ",cds_name,"@expression$", lineage)
      eval(parse(text=input))
      if(gene %in% colnames(exp)){
        input = paste0("fit = ",cds_name,"@expectation$", lineage,"[,'",gene,"']")
        eval(parse(text=input))
      }
      else{
        fit = rep(0, N)
      }
      fits = c(fits, fit)
      dd = rbind(dd, cbind(seq(from=0, to=max.pt, by = max.pt/(N-1)), fit, rep(lineage, length(fit))))
    }
    ymax = max(fits)
    colnames(dd) <- c("pseudotime", "fit", "lineage")
    dd = as.data.frame(dd)
    dd$pseudotime <- as.numeric(dd$pseudotime)
    dd$fit <- as.numeric(dd$fit)
    dd$lineage <- factor(dd$lineage, levels = lineages)
  }
  q <- ggplot(data = dd)
  if(points == T){
    for(M in 1:length(lineages)){
      loop_input1 = paste0("geom_point(aes_string(x='pseudotime',y = '", paste0('exp_', lineages[M]), "',color='pseudotime'), size=I(", point_size, "))")
      loop_input2 = paste0("scale_color_gradient2(lineages[M],low='grey', ", "high='",colors[M],"')")
      loop_input3 = "new_scale_color()"
      loop_input4 = paste0("geom_line(aes_string(x='pseudotime', y = '", paste0('fit_', lineages[M]), "',size = I(", line_size, ")), color = '", colors[M],"')")
      q <- q + eval(parse(text=loop_input1)) + eval(parse(text=loop_input2)) + eval(parse(text=loop_input3)) + eval(parse(text=loop_input4))
    }
  }
  else{
    q <- q + geom_line(aes(x = pseudotime, y = fit, color = lineage), size = I(line_size)) + scale_color_manual(values = colors)
  }
  if(length(scale.lineage) == 1){
  input = paste0(cds_name,"@lineages$", scale.lineage)
  cells = eval(parse(text = input))
  age = meta[cells,c("age_num", "age_range")]
  }
  else{
  age = meta[,c("age_num", "age_range")]
  }
  age = age[order(age$age_num),]
  window = nrow(age)/N
  step = ((nrow(age)-window)/N)
  age.comp = SlidingWindow("mean", age$age_num, window, step)
  d = seq(from=0, to=max.pt, by = max.pt/(N-1))
  d = cbind(as.data.frame(d), age.comp)
  q <- q + scale_y_log10()
  if(age.scale == T){
    breaks.list = c(0)
    for(age.point in age.points){
      age.break = quantile(age[age$age_range == age.point,]$age_num, 0.95)
      age.break = d[which.min(abs(d[,2]-age.break)),1]
      breaks.list = append(breaks.list, age.break)
    }
    q <- q + scale_x_continuous(breaks = breaks.list, labels = breaks.labels)
  }
  q <- q + ylim(y = c(0,ymax))
  q <- q + monocle_theme_opts() + ylab("Expression") + xlab("Pseudotime") + ggtitle(gene) + theme(legend.key.size = unit(legend.key.size, 'cm'), plot.title = element_text(size = plot.title.size, face="bold", hjust = 0.5), axis.text=element_text(size=text.size), axis.text.x=element_text(angle = 60, hjust=1), axis.title=element_blank(), legend.text=element_text(size=legend.text.size), legend.title=element_text(size=text.size, face = "bold"), legend.position = legend_position)
  q
}



