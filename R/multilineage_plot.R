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
plot_multiple <- function(cds, gene, lineages, meta = NULL, points = T, age.scale = T, point_size = 0.1, line_size = 1, text.size = 14, plot.title.size = 36, legend.key.size = 0.5, legend.text.size = 10, colors = c("red", "blue", "green", "cyan", "magenta", "purple", "orange", "black", "yellow", "tan"), N = 500, legend_position = "none"){
  cds_name = deparse(substitute(cds))
  input = paste0(cds_name,"@expression$", lineages[1])
  N = nrow(eval(parse(text = input)))
  pts = c()
  for(lineage in lineages){
    input = paste0(cds_name,"@pseudotime$", lineage)
    pt = eval(parse(text = input))[,1]
    pts = c(pts, pt)
  }
  max.pt = max(pts)
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
    for(N in 1:length(lineages)){
      loop_input1 = paste0("geom_point(aes_string(x='pseudotime',y = '", paste0('exp_', lineages[N]), "',color='pseudotime'), size=I(", point_size, "))")
      loop_input2 = paste0("scale_color_gradient2(lineages[N],low='grey', ", "high='",colors[N],"')")
      loop_input3 = "new_scale_color()"
      loop_input4 = paste0("geom_line(aes_string(x='pseudotime', y = '", paste0('fit_', lineages[N]), "',size = I(", line_size, ")), color = '", colors[N],"')")
      q <- q + eval(parse(text=loop_input1)) + eval(parse(text=loop_input2)) + eval(parse(text=loop_input3)) + eval(parse(text=loop_input4))
    }
  }
  else{
    q <- q + geom_line(aes(x = pseudotime, y = fit, color = lineage), size = I(line_size)) + scale_color_manual(values = colors)
  }
  age = meta[,c("age_num", "age_range")]
  age = age[order(age$age_num),]
  window = nrow(meta)/N
  step = ((nrow(meta)-window)/N)
  age.comp = SlidingWindow("mean", age$age_num, window, step)
  d = cbind(as.data.frame(seq(from=0, to=max.pt, by = max.pt/(N-1))), age.comp)
  second_third = quantile(age[age$age_range == "2nd trimester",]$age_num, 0.95)
  second_third = d[which.min(abs(d[,2]-second_third)),1]
  birth = quantile(age[age$age_range == "3rd trimester",]$age_num, 0.95)
  birth = d[which.min(abs(d[,2]-birth)),1]
  infant = quantile(age[age$age_range == "0-1 years",]$age_num, 0.95)
  infant = d[which.min(abs(d[,2]-infant)),1]
  four = quantile(age[age$age_range == "2-4 years",]$age_num, 0.95)
  four = d[which.min(abs(d[,2]-four)),1]
  ten = quantile(age[age$age_range == "4-10 years",]$age_num, 0.95)
  ten = d[which.min(abs(d[,2]-ten)),1]
  q <- q + scale_y_log10()
    if(age.scale == T){
  q <- q + scale_x_continuous(breaks = c(0, second_third, birth, infant, four), labels = c("2nd", "3rd", "birth", "1y", "4y"))
  }
  q <- q + ylim(y = c(0,ymax))
  q <- q + monocle_theme_opts() + ylab("Expression") + xlab("Pseudotime") + ggtitle(gene) + theme(legend.key.size = unit(legend.key.size, 'cm'), plot.title = element_text(size = plot.title.size, face="bold", hjust = 0.5), axis.text=element_text(size=text.size), axis.text.x=element_text(angle = 60, hjust=1), axis.title=element_blank(), legend.text=element_text(size=legend.text.size), legend.title=element_text(size=text.size, face = "bold"), legend.position = legend_position)
  q
}


