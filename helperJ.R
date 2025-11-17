################################################################################
#source(file.path(dirname(wd), "helperJ.R"))
library(stringr)
library(tidyverse)
#cName
#df$cName <- apply(df[,c('Exp','PCR','Enz','Temp',)], 1, paste0, collapse = ".")
#sessionInfo()
#qOrder(index_table,'cName')
#mutate(FP = factor(FP, levels = FP_order))
################################################################################
#s=paste0("dollarDV=c('df$",paste0(cdv,collapse="','df$"),"')")
#eval(parse(t=s))
#s=paste0("df$cName[rowCount]=paste0(",paste0(dollarDV,collapse='[rowCount],'),"[rowCount])")
#eval(parse(t=s))
################################################################################
#cdf$Temp = factor(cdf$Temp, levels = c('55C','45C')
#<- tidyr::gather(mds_data, key = "Variable", value = "Value", -Temp, -PCR, -Exp)
#################################################################################
################################################################################
get_last_var <- function() {
  # Check the calling environment first
  env <- parent.frame()
  objects <- ls(envir = env)
  vars <- objects[sapply(objects, function(x) !is.function(get(x, envir = env)))]

  # Check if there are any variables
  if (length(vars) > 0) {
    last_var_name <- tail(vars, 1)
    return(last_var_name)
  } else {
    # If no variables found, check the global environment
    env <- .GlobalEnv
    objects <- ls(envir = env)
    vars <- objects[sapply(objects, function(x) !is.function(get(x, envir = env)))]

    if (length(vars) > 0) {
      last_var_name <- tail(vars, 1)
      return(last_var_name)
    } else {
      cat("No variables found in the global environment.\n")
      return(NULL)  # Return NULL if no variables are found
    }
  }
}
#get_last_var()
##################################################################################
# dataframe,vector of column names to be factorized levels will be in same order
#
qMatrix=function(df,x,y,d,f='mean'){
  df=df
  ux=unique(df[,x])
  uy=unique(df[,y])
  m=matrix(0,length(uy),length(ux))
  row.names(m)=uy
  colnames(m)=ux
  s=paste0("ag=aggregate(",d,"~",x,"+",y,",data=df,",f,")")
  print(s)
  eval(parse(t=s))
  for (c in ux){
    for(r in uy){
      m[r,c]=ag[ag[,y]==r & ag[,x]==c,d]

    }}
  m
}
######### Factorize a df ####################################
################################################################################
# dataframe,vector of column names to be factorized levels will be in same order
#
qFactordf=function(df,factors){
  for(f in factors){
    s=paste0("df$",f,"=factor(df$",f,",levels=unique(df$",f,"))")
    eval(parse(t=s))
  }
  df
}
################################################################################
######### simple aggregate will add more cName, factorize ####################################
################################################################################
# dataframe,named vector or list?
#qAg=function(df='df',iv='iv',dv='dv',fun='mean'){

sAg=function(df,values,categories,fn="mean"){
  #df=df;av=values;avs=categories;fn=fn;
  s=paste0("aggregate(cbind(",paste(values,collapse=","),")~",paste(categories,collapse="+"),",data=df,",fn,")")
  print(s)
  eval(parse(t=s))


}

################################################################################
######### quick aggregate will add more cName, factorize ####################################
################################################################################
# dataframe,named vector or list?
#qAg=function(df='df',iv='iv',dv='dv',fun='mean'){

qAg=function(df,av,avs,fn="mean"){
  df=df;av=av;avs=avs;fn=fn
  s=paste0("aGdf=aggregate(cbind(",paste(av,collapse=","),")~",paste(avs,collapse="+"),",data=df,",fn,")")
  print(s)
  eval(parse(t=s))
  for(d in dv){
    s=paste0("aGdf$",d,"=as.character(aGdf$",d,")")
    eval(parse(t=s))
  }

  for (r in 1:nrow(aGdf)){
    aGdf$cNames[r]=paste(aGdf[r,1:length(dv)],collapse=".")
    aGdf$Jitx[r]=1-(jitter(1,30))
    aGdf$Jity[r]=1-(jitter(1,30))
  }
  dv=c(dv,"cNames")
  aGdf=qFactordf(aGdf,dv)
  aGdf
}
################################################################################
#log 2
# natural_log((count+1)/(geometric_mean (count+1)))
################################################################################
qLog2=function(x,s='off'){
  if(s=='on'){print('qLog2 vector');print(x)}

  #Test for negatives
  if (any(x < 0)) {
    stop("Data frame contains negative values. CLR transformation requires non-negative values.")
  }
  if(sum(x)==0){return(NA)}
  log_base_2=log2(x/mean(x[x > 0]))

  # print(paste('mean(x)',mean(x)))
  print(paste('mean(x[x > 0]))',mean(x[x > 0])))
  print(paste('log_base_2=log2(x/mean(x[x > 0]))',log_base_2=log2(x/mean(x[x > 0]))))

  return(log_base_2)
  if(sum(x)==0){return(NAN)}
}
################################################################################
################################################################################
#centered log transformation
# natural_log((count+1)/(geometric_mean (count+1)))
################################################################################
# ignores zeros
qRCLR=function(x,s='off'){
  if(s=='on'){print('qCLR vector');print(x)}

  #Test for negatives
  if (any(x < 0)) {
    stop("Data frame contains negative values. CLR transformation requires non-negative values.")
  }
  centered_log=x
  centered_log[centered_log>0] <- log(x[x>0]) - mean(log(x[x>0]))
  return(centered_log)

}
################################################################################
# ignores zeros
qCLR=function(x,s='off'){
  if(s=='on'){print('qCLR vector');print(x)}

  #Test for negatives
  if (any(x < 0)) {
    stop("Data frame contains negative values. CLR transformation requires non-negative values.")
  }
  centered_log=x
  centered_log[centered_log>0] <- log(x[x>0]) - mean(log(x[x>0]))
  return(centered_log)

}

################################################################################
################################################################################
# PermANOVA Test
# perm_ANOVA_p=qPermANOVA(.,sample_name,FP,ASV,counts,s='on',r='p')
################################################################################
qPermANOVA <- function(df, row_names, grouping, features, values, r = 'p', s = 'off') {

  if (s == 'on') {
    print('PermANOVA row_names')
    print(row_names)
    print('PermANOVA grouping')
    print(grouping)
    print('PermANOVA features')
    print(features)
    print('PermANOVA values')
    print(values)
  }

  raw_df <- df

  df <- as.data.frame(cbind(row_names, grouping, features, values))
  colnames(df) <- c('row_name', 'grp', 'feature', 'value')

  pre_matrix <- df %>%
    pivot_wider(names_from = feature, values_from = value)

  distance_matrix <- pre_matrix %>%
    select(-grp) %>%
    mutate_at(-1, as.numeric) %>%
    column_to_rownames(var = 'row_name')

  tryCatch({
    result <- adonis2(distance_matrix ~ pre_matrix$grp, permutations = 999)
    if (r == 'p') {
      return(result$Pr[1])
    }
    if (r == 'R2') {
      return(result$R2[1])
    }
    if (r == 'Df') {
      return(result$Df[1])
    }
  }, error = function(e) {
    if(s=='on'){print(paste("Error:", e$message))}
    return(NA)
  })
}

################################################################################
# Kruskal-Wallace Test
# qKW(df, df$clr_counts, df$FP)
# df%>%
#   group_by(ASV)%>%
#   summarize(KW=qKW(.,clr_counts,FP))
################################################################################
qKW=function(df,x,grouping,s='off'){
  if(s=='on'){print('Kruskal-Wallace vector');print(x)}
  #Test for negatives
  KW=kruskal.test(x ~ grouping, data = df)
  return(KW$p.value)
}
################################################################################
#return index from population vector
#input a vector representing populations
################################################################################
# population vector
qShannon=function(x,s='off'){
  if(s=='on'){print('shannon vector');print(x)}
  x=x[x!=0];s=sum(x);-sum(x/s*log(x/s))}

qEvenness=function(x,s='off'){
  if(s=='on'){print('evenness vector');print(x)}
  x=x[x!=0];s=sum(x);-sum(x/s*log(x/s))/log(length(x))}

# y = reference population ratios. default is even distribution
qIdealScore=function(x,y=c(rep(1/length(x),length(x))),s='off'){
  if(s=='on'){print('Ideal Score vector');print(x)}
  sum(abs((x/sum(x))-(y/sum(y))))*100}

qBray <- function(x, s='off') {
  if (s == 'on') {
    print('Input vector:')
    print(x)
  }
  x = x[x != 0]
  s = sum(x)
  sum(abs(x/s - x/s))^2 / sum(x/s + x/s)
}
qEU <- function(x, y, s='off') {
  if (s == 'on') {
    print('Input vectors:')
    print(x)
    print(y)
  }
  sum((x - y)^2)^0.5
}


################################################################################
#y=y/sum(y)
#x=x/sum(x)
#v=x-y
#IdS=sum(abs(v))*100
#IdS
################################################################################
#########return factorized and leveled data frame of given vector column names ####################################
################################################################################
# dataframe,named vector or list?
qLeveldf=function(x,y){
  df=x
  dv=y
  for(d in dv){
    s=paste0("df$",d,"=factor(df$",d,",levels=unique(df$",d,"))")
    print(s)
    eval(parse(t=s))
  }
  df
}
################################################################################
################################################################################
#  for(f in 1:length(l)){assign(names(l[f]),unlist(l[f]))}
#qw[order(qw$Time,qw$Temp),]




dss2df <- function(dss) {data.frame(width=width(dss), seq=as.character(dss), names=names(dss))}
ev=function(x){eval(parse(text=x))}

pasteP=function(x){paste(x,collapse="+")}
pasteC=function(x){paste(x,collapse=",")}
pasteQ=function(x){noquote(paste0(x))}
nQ=function(x){noquote(x)}
pasteQC=function(x){noquote(paste(x,collapse=","))}
evP=function(x){eval(parse(text=paste0(x)))}


# Create a data frame of naming elements from a named list of vectors to make each row a unique name
#pasteP=function(x){paste(x,collapse="+")}
nameDF=function(x){
  df=expand.grid(x)
  # df=data.frame(df)
  dN=colnames(df)
  names(dN)=dN
  for (d in colnames(df)){
    # eval(parse(text=noquote(paste0("df$",d,"=as.factor(df$",d,")"))))
    #eval(parse(text=noquote(paste0("levels(df$",d,")=",dv))))
    dN[d]=paste0("df$",d)
    print(d)
  }
  s=paste0("df=df[order(",nQ(pasteC(dN)),"),]")
  eval(parse(text=pasteQ(s)))
  for (d in dN){
    eval(parse(text=noquote(paste0(d,"=as.character(",d,")"))))
  }
  df=as.data.frame(df)
  row.names(df)= 1:nrow(df)
  df
}
################################################################################
# make a factor and add levels
# input dataframe
# input vector (name= df$name) of ordered levels
################################################################################
# dataframe,'vector name',order vector
qFactor=function(x,y,z){
  df=x
  eval(parse(text=paste0("df$",y,"=factor(df$",y,",levels=z)")))
  df
}
################################################################################
################################################################################
# make a factor and add levels
# input dataframe
# input vector (name= df$name) of ordered levels
################################################################################
# dataframe,named vector or list?
qLevel=function(df,levelVector,vectorName){
  #for(n in 1:length(nL)){assign(names(nL[n]),unlist(nL[n]))}
  df=df
  lv=levelVector
  vn=vectorName
  eval(parse(text=paste0("df$",vn,"=as.factor(df$",vn,")")))
  eval(parse(text=paste0("df$",vn,"=factor(df$",vn,",levels=",lv,")")))
}
################################################################################
############################## reorder ####################################
# dataframe,ordering vector
qOrder=function(df,columns){

  #yv=y
  # names(yv)=yv
  # for (c in columns){
  #   eval(parse(text=paste0("df$",c,"=as.factor(df$",c,")")))
  #   eval(parse(text=paste0("df$",c,"=factor(df$",c,",levels=",c,")")))
  # }
  if(length(columns)>1){
    for (c in columns){
      s=paste0("df=df[order(",(pasteC(c)),"),]")
      eval(parse(t=pasteQ(s)))
    }}else{
      s=paste0("df=df[order(df$",columns,"),]")
      eval(parse(t=s))}
  df
}
################################################################################
#combineClm ####################################
# df, named List of vector columns. at least 2
# first is what will be combined
# the rest identify the original columns
# all vectors must be in alignment
#cNames="your vectors"
#yourVectors="column names to be combined under vector name"
#cNames=c("prReads","prNames")
#s=paste0("cList=list(",pasteC(cNames),")")
#(eval(parse(text=s)))
#names(cList)=cNames
#cdf=combineClm(df,cList)

#bRead=prPosBaseName
#fPRname=as.character(prPosBaseNamedf$Var1)
#position=as.character(prPosBaseNamedf$Var2)
#base=as.character(prPosBaseNamedf$Var3)
#cNames=c('bRead','prPosBaseName','fPRname','position','base')
#s=paste0("cList=list(",pasteC(cNames),")")
#(eval(parse(text=s)))
#names(cList)=cNames
#cdf=combineClm(df,cList)
############################## combineCol ####################################
combineClm=function(df,cList){




  #       "m13","n14","o15","p16","q17","r18","s19","t20","u21","v22","w23",
  #      "x24","y25","z26","aa27","bb28","cc29","dd30","ee31")
  cList
  unlist(cList[2])
  for(w in 1:length(cList)){assign(names(cList[w]),unlist(cList[w]))}
  increaseBy=length(unlist(cList[1]))

  for(z in names(cList)){
    s=paste0("df$",z,"=NA")
    print(s)
    eval(parse(text=s))
  }

  mLength=increaseBy+length(cList)
  abc=paste0("a",1:mLength)
  abc=paste0("a",1:increaseBy)
  s=paste0(paste0(abc,collapse="="),"=df")
  eval(parse(t=s))

  #abc=c("a1,a2,...

  #a1=a2=...=df



  #### first vector set to specific columns
  for (v in 1:increaseBy) {
    w=unlist(cList[1])
    s=paste0(abc[v],"[,names(cList[1])]=df[,w[v]]")
    print(s)
    print(names(cList[1]))
    print(w[v])
    eval(parse(t=s))
  }
  #### other vectors set to first vector
  wVector=vector(length=nrow(df))
  for(u in 2:length(cList)){
    w=unlist(cList[u])
    for (v in 1:increaseBy) {
      wVector[]=w[v]
      s= paste0(abc[v],"$",names(cList[u]),"=wVector")
      print(s)
      eval(parse(t=s))
    }}

  #combine rows

  nR=(nrow(df)-1)*increaseBy
  for(q in 1:increaseBy){
    iR=nR+q
    s=paste0(abc[q],"$rn=seq(q,iR,increaseBy)")
    print(s)
    eval(parse(t=s))
    s=paste0("row.names(",abc[q],")=seq(q,iR,increaseBy)")
    print(s)
    eval(parse(t=s))}

  s=paste0("df=rbind(",pasteC(abc[1:increaseBy]),")")
  print(s)
  eval(parse(t=s))
  df
}
###############################################################################
# qPrint <- function(var=get_last_var()) {
#   # Get the name of the variable
#   var_name <- deparse(substitute(var))
#
#   # Convert the variable to a string, handling vectors appropriately
#   if(is.vector(var)) {
#     value_str <- paste(var, collapse=", ")
#   } else {
#     value_str <- toString(var)
#   }
#
#   # Remove slashes from the string
#   value_str <- gsub("/", "", value_str)
#
#   # Create the output string
#
#   output_string <- sprintf('%s = %s', var_name, value_str)
#
#   # Print the output
#   print(output_string)
# }
qPrint <- function(var = NULL) {
  # Get the name of the variable from the environment if not provided
  if (is.null(var)) {
    env <- parent.frame()
    objects <- ls(envir = env)
    vars <- objects[sapply(objects, function(x) !is.function(get(x, envir = env)))]

    # Check if there are any variables
    if (length(vars) > 0) {
      last_var_name <- tail(vars, 1)
      var_name=last_var_name
      var=get(var_name, envir = env)
    } else {
      # If no variables found, check the global environment
      env <- .GlobalEnv
      objects <- ls(envir = env)
      vars <- objects[sapply(objects, function(x) !is.function(get(x, envir = env)))]

      if (length(vars) > 0) {
        last_var_name <- tail(vars, 1)
        var_name=last_var_name
        var=get(var_name, envir = env)
      } else {
        cat("No variables found in the global environment.\n")
        return(NULL)  # Return NULL if no variables are found
      }
    }
  }else{
    var_name <- deparse(substitute(var))
  }

  # Convert the variable to a string, handling vectors appropriately
  if (is.vector(var)) {
    value_str <- paste(var, collapse = ", ")
  } else {
    value_str <- toString(var)
  }

  # Remove slashes from the string
  value_str <- gsub("/", "", value_str)

  # Create the output string
  output_string <- sprintf('%s = %s', var_name, value_str)

  # Print the output
  print(output_string)
}

#qPrint()  # Prints the last variable in the calling environment
################################################################################
# Define the tPrint function
tPrint <- function(message) {
  # Get the current system time
  current_time <- Sys.time()

  # Print the message followed by the current time
  cat(paste0(message, " at ", format(current_time, "%Y-%m-%d %H:%M:%S")), "\n")
}

# Example usage
tPrint("start")

################################################################################
# functions useful for building graphs
#
one_sample_hist <- function(x, mu, xlab='', ylab='Frequency', main='',
                            col='skyblue', border='navy', ...) {
  hist(x, xlab=xlab, ylab=ylab, col=col,las=1, border=border, main=main)

  # add a line at the sample mean
  abline(v=mean(x), lwd=2, lty=2, col='navy')

  # add a line at the hypothetical mean
  if(!missing(mu)) {
    abline(v=mu, lwd=4, col='orangered')

    t_results = t.test(x, mu=mu)

    # round the p-value to 3 places
    pval = round(t_results$p.value, 3)

    # draw some text on the graph
    text(par('usr')[1]*1.1, par('usr')[4]*.85,
         paste0('p=', pval), cex=1.5, col='gray30')
  }

  return(t_results)
}


.jitr = function(x, len=length(x), r=0.35) {
  x + runif(len, -r, r)
}

.add_points = function(x, y, jitr_x=.2, pch=19, ...) {
  points(.jitr(x, length(y), r=.2), y, pch=pch, ...)
}


barplot_with_scatter = function(y, xlab='', ylab='',
                                col='gray30', border=col, border_only=TRUE, pt.col=adjustcolor(col,.75),
                                ebar.col = pt.col,
                                ebar.lwd=2, pch=16, jitr_x, ebar.code=0, draw_mean=FALSE, data, ylim, ...) {

  if(inherits(y, 'formula')) {
    by_cond <- aggregate(y, list, data=data)
    if(dim(by_cond)[2L] > 2) {
      stop('Function only supports one X variable at present')
    }

    m <- do.call(cbind, by_cond[,ncol(by_cond)])
    colnames(m) = by_cond[,1]
    if(missing(xlab)) xlab <- attr(terms(y), 'term.labels')
    if(missing(ylab)) ylab <- colnames(model.frame(y, data=data))[1]

    y <- m
  }

  y <- as.matrix(y)

  if(missing(ylim)) {
    lim <- range(c(0, pretty(c(y))))
  } else {
    lim = ylim
  }

  msem = t(apply(y, 2, mean_se))

  bar.col = col
  if(border_only) {
    bar.col = NA
  }

  xp=barplot(msem[,1], las=1,
             ylim=lim,
             col=bar.col, border=border,
             ylab=ylab, names.arg=colnames(y),
             xlab=xlab, ...)

  if(missing(jitr_x)) {
    jitr_x = diff(xp[1:2])/3
  }

  pt.col = rep_len(pt.col, nrow(msem))
  pch = rep_len(pch, nrow(msem))

  for(ii in 1:nrow(msem)) {
    .add_points(xp[ii], y[,ii],
                col=pt.col[ii], pch=pch[ii], jitr_x=jitr_x)

    if(draw_mean) lines(xp[ii] + c(-jitr_x, jitr_x),
                        rep(msem[ii,1], 2), col=col[ii], lwd=ebar.lwd)
  }

  .error_bars(xp, msem, lwd=ebar.lwd, col=ebar.col, code=ebar.code)
}

# plot.two_groups = function(data, line.col='gray30', mean.col='dodgerblue3') {
#   means = colMeans(data)
#   sems = apply(data, 2, standard_error)
#
#   x = barplot(means, ylim=c(0,1.1*max(data)), las=1, border=line.col, col=NA)
#   apply(data, 1, lines, col=line.col, x=x)
#   lines(x, means, col=mean.col, lwd=3)
#   error_bars(x, means, sems, col=mean.col, lwd=3)
# }

## Making error bars from x, y, and a measure of error
.error_bars = function (x, y = NULL, err = NULL, length = 0.05, type = "n",
                        col = "black", pt.col = col, code = 0, ...) {
  if (is.null(y)) {
    if (is.matrix(x)) {
      y <- x[, 1]
      err <- x[, 2]
    }
    else {
      y <- x
    }
    x <- seq_along(y)
  }

  if (is.matrix(y)) {
    err <- y[, 2]
    y <- y[, 1]
  }

  if (is.null(err)) {
    err <- y
    y <- x
    x <- seq_along(y)
  }
  .ebars.y(x, y, err, length, code = code, col = col, ...)
  points(x, y, type = type, col = pt.col, ...)
}

.ebars.x = function(x, y, sem, length = 0.05, code=3, ...) {
  arrows(x - sem, y, x + sem, y, angle = 90, code = code, length = length, ...)
}

.ebars.y = function(x, y, sem, length = 0.05, code = 3, ...) {
  # if (up) {
  arrows(x0 = x, y0 = as.numeric(y-sem), y1 = as.numeric(y + sem), angle = 90, code = code, length = length, ...)
  # }
  # if (down) {
  # arrows(x0 = x, y0 = as.numeric(y), y1 = as.numeric(y - sem), angle = 90, code = code, length = length, ...)
  # }
}

panel.fixed <- function(x,y, method='pearson', ...) {

  .cor = cor.test(x,y, method=method)
  .col = 'gray60'
  if(.cor$p.value < 0.05) {
    .col = ifelse(.cor$estimate > 0, 'palevioletred', 'dodgerblue3')
  }

  legend('center', legend=sprintf('%+.2f',.cor$estimate), adj = (.3),
         bty='n',  text.col=.col)
}
################################################################################
panel.cor <- function(x,y, method='pearson', ...) {

  .cor = cor.test(x,y, method=method)
  .col = 'gray60'
  if(.cor$p.value < 0.05) {
    .col = ifelse(.cor$estimate > 0, 'palevioletred', 'dodgerblue3')
  }

  legend('center', legend=sprintf('%+.2f',.cor$estimate), adj = (.3),
         bty='n', cex=max(.5, .cor$estimate*4), text.col=.col)
}
##############################################################################
.se = function(x, na.rm=TRUE, do_round=FALSE) {
  v = sqrt(var(x, na.rm = na.rm)/ sum(!is.na(x)))

  if(do_round) {
    if(v > 1) {
      v = round(v, 1)
    } else {
      v = round(v, abs(floor(log10(v))-1))
    }
  }

  return (v)
}

se = .se


points.jitter <- function(x, y, ..., shake=0.35) {
  x = rep(x, length(y))
  x = x + runif(length(x), -shake, shake)

  points(x,y, ...)
}



draw.error_bars <- function (x, y = NULL, sem = NULL, col='black', length = 0.05, type = "n",
                             pt.col = col, code = 3, direction='y', ...)
{

  if(missing(x) && is.null(y)) {
    stop('One of x or y must be provided')
  }

  if(is.data.frame(y)) {
    y = as.matrix(y)
  }

  if(missing(x)) {
    x <- if(is.matrix(y)) {
      1:nrow(y)} else {
        seq_along(y)
      }
  }

  if (is.null(y)) {
    if (is.matrix(x)) {
      y <- x[, 1]
      sem <- x[, 2]
    }
    else {
      y <- x
    }
    x <- seq_along(y)
  }
  if (is.matrix(y)) {
    sem <- y[, 2]
    y <- y[, 1]
  }
  if (is.null(sem)) {
    sem <- y
    y <- x
    x <- seq_along(y)
  }


  if(direction=='y') {
    .ebars.y(x, y, sem, length, code = code, col = col, ...)
  } else {
    .ebars.x(x,y,sem,length,code,col=col, ...)
  }
  points(x, y, type = type, col = pt.col, ...)
}


mean_se = function(x, na.rm=TRUE) {
  c('mean'=mean(x, na.rm=na.rm), 'SE'=.se(x, na.rm))
}


describe <- function(x, na.rm=TRUE, do_round=TRUE) {
  msd <- mean_sd(x, na.rm)

  c('N' = sum(!is.na(x)), msd, 'SE'=.se(x, na.rm=na.rm, do_round = do_round), 'isNA'=sum(is.na(x)),
    'min' = min(x, na.rm=na.rm),
    'mid' = median(x, na.rm=na.rm),
    'max' = max(x, na.rm=na.rm)
  )
}

mean_sd = function(x, na.rm=TRUE) {
  v = c('mean'=mean(x, na.rm=na.rm), 'SD'=sd(x, na.rm))

  if(!is.na(v[2]) && v[2] > 1) {
    v = round(v, 1)
  }

  return(v)
}


panel.viobox <- function(..., col='transparent', box.ratio) {
  panel.violin(..., col=col, varwidth = F, box.ratio = box.ratio)
  panel.bwplot(..., fill=NULL, box.ratio = .1)
}


panel.mse_plot <- function (x, y,
                            pch = if (is.null(groups)) dot.symbol$pch else sup.symbol$pch,
                            col = if (is.null(groups)) dot.symbol$col else sup.symbol$col,
                            lty = dot.line$lty, lwd = dot.line$lwd, col.line = dot.line$col,
                            groups = NULL, subscripts, pt.alpha = 1, ...)
{
  dot.line <- trellis.par.get("dot.line")
  dot.symbol <- trellis.par.get("dot.symbol")
  sup.symbol <- trellis.par.get("superpose.symbol")

  panel.xyplot(jitter(as.numeric(x)), y, subscripts=subscripts,
               groups=groups, col=col, pch=pch, lwd=lwd, ...)

  x <- as.numeric(x)
  y <- as.numeric(y)

  # remove missing values
  ind = ! (is.na(x) | is.na(y))
  x = x[ind]
  y = y[ind]


  #
  ux = as.numeric(sort(unique(x)))
  rx <- mean(diff(sort(ux)))
  #
  m = tapply(y[order(x)], sort(x), FUN = mean)
  se = tapply(y[order(x)], sort(x), FUN = .se)

  # assign('vals', list(
  #   'u'=ux, 'r'=rx, 'm'=m, 'x'=x, 'y'=y, 's'=subscripts
  # ), envir = globalenv())
  #
  #   print(ux)
  #   print(x)
  # print(subscripts)

  lsegments(
    x0 = ux-.2*rx,x1=ux+.2*rx,
    y0 = m, y1 = m,
    lwd=3, col=adjustcolor(col, alpha.f=100)
  )

  col = rep(col, length=length(m))

  mapply(
    function(.x, .y, .se, .col) {
      lpolygon(border = NA,
               x = c(.x-.1*rx, .x+.1*rx, .x+.1*rx, .x-.1*rx),
               y = c(.y-.se, .y-.se, .y+.se, .y+.se),
               col=adjustcolor(.col,alpha.f = 100), alpha=.2)
    }, ux, m, se, col)

}
panel.mse_dotplot <- panel.mse_plot

panel.se <- function(x, y, col=plot.line$col, line.col=col, pch=plot.symbol$pch, alpha.se=.25, ...) {
  v <- list(...)
  if('groups' %in% names(v)) {
    plot.line = trellis.par.get('superpose.line')
    plot.symbol = trellis.par.get('superpose.symbol')
    return(
      panel.superpose(x, y, v$subscripts, v$groups,
                      panel.groups=panel.se,
                      col=col, line.col=line.col, alpha.se=alpha.se, pch=pch, ...)
    )
  }

  plot.line <- trellis.par.get("plot.line")
  plot.symbol <- trellis.par.get("plot.symbol")
  xs <- if(is.factor(x)) {
    factor(c(levels(x) , rev(levels(x))), levels=levels(x))
  } else {
    xx <- sort(unique(x))
    c(xx, rev(xx))
  }
  means <- tapply(y, x, mean, na.rm = T)
  ses <- tapply(y, x, .se)

  panel.lines(head(xs, length(means)), means, col=col, type='o', pch=pch, cex=.8)
  panel.polygon(xs, c(means + ses, rev(means - ses)), col = col, alpha=alpha.se, border=NA)
}


# panel.mse_dotplot <- function(x, y, subscripts, horizontal=TRUE, col, ...) {
#
#   m = tapply(y, x, mean)
#   se = tapply(y, x, .se)
#
#   if(missing(col)) {
#     col = if('groups' %in% names(list(...))) {
#       trellis.par.get('superpose.symbol')$col
#     } else {
#       trellis.par.get('plot.symbol')$col
#     }
#   }
#
#   ux = as.numeric(unique(x))
#   rx <- mean(diff(sort(ux)))
#
#   lsegments(
#     x0 = ux-.2*rx,x1=ux+.2*rx,
#     y0 = m, y1 = m,
#     lwd=3, col=col
#   )
#
#   col = rep(col, length=length(m))
#
#   mapply(
#     function(.x, .y, .se, .col) {
#       # print(c(.x, .y, .se))
#       lpolygon(border = NA,
#         x = c(.x-.1*rx, .x+.1*rx, .x+.1*rx, .x-.1*rx),
#         y = c(.y-.se, .y-.se, .y+.se, .y+.se),
#         col=.col, alpha=.2)
#     }, ux, m, se, col)
#
#   panel.dotplot(jitter(as.numeric(x)),y,horizontal = horizontal, lwd = 0, col=col, subscripts=subscripts,...)
# }


table_as_counts_and_percentages <- function(tbl) {
  pct = sprintf('%2.1f', 100 * t(tbl) / colSums(tbl))

  t_tbl = t(tbl)

  m = matrix(mapply(function(p,c) {
    paste0(c, ' (', p, ')')
  }, c(pct), c(t_tbl)), nrow=nrow(t_tbl)
  )
  colnames(m) = colnames(t_tbl)
  rownames(m) = rownames(t_tbl)

  return(m)
}

logistic_summary <- function(glm_res, round=5) {
  stopifnot(
    family(glm_res)$family == 'binomial' &&
      family(glm_res)$link == 'logit'
  )

  summ = summary(glm_res)

  cf = summ$coefficients
  exp = matrix(round(digits = round, exp(cf[,1])))
  colnames(exp) = 'exp(Est)'

  cf = cbind(exp, cf)
  summ$coefficients = cf

  return (summ)
}


count_NA <- function(x) {
  sum(is.na(x))
}

mse_count <- function(x) {
  c(mean_se(x), 'N'=sum(!is.na(x)))
}


ebars.polygon = function (x, y, sem, alpha = 100/255, col = "black", fill = col,
                          stroke = col, border = NA, add_line = TRUE, lwd = 1, ...)
{

  if(missing(y)) {
    y = x
    x = seq_along(dim(y)[1])
  }

  if(missing(sem) && ncol(y) > 1) {
    sem = y[,2]
    y = y[,1]
  }

  is_finite = is.finite(y) & is.finite(sem)
  if (all(!is.finite(sem))) {
    is_finite <- is.finite(y)
    sem <- 0 * y
  }

  x = x[is_finite]
  y = y[is_finite]
  sem = sem[is_finite]
  sem = abs(sem)
  polygon(c(x, rev(x)), c(y + sem, rev(y - sem)), border = border,
          col = adjustcolor(fill, alpha))
  if (add_line)
    lines(x, y, col = stroke, lwd = lwd, ...)
}


#df <- df[order(levels(df$region)),]
############################################################
############################################################
#Alpha Diversity? Functions
# Tibble format?
#   columns, group,values
# values will be our x
############################################################
inverse_simpson <- function(x){
  n <- sum(x)
  1/(sum(x * (x-1) / (n * (n-1))))
}

qRichness<- function(x){
  sum(x>0)
} # compare to vegan specnumber(x) 1 - sum((x/n)^2)

shannon<-function(x){
  n <- sum(x)
  relative_abundance = ra=x[x>0]/n
  -sum(ra*log(ra))
} # compare to vegan diversity(x,index='shannon')


simpson<- function(x){
  n <- sum(x)
  sum(x * (x-1) / (n * (n-1)))
} # compare to vegan diversity(x,index='simpson') 1 - sum((x/n)^2
##################################################################
################################################################################
untried_perform_ttest <- function(df, group_var, count_var,s='off') {
  if(s=='on'){print('grouping variables');print(group_var);print('counting vector');print(count_var)}

  stats = c("sd", "mean")
  df %>%
    group_by({{group_var}},Type)  %>%
    summarize(
      ttest_p_value = tryCatch(t.test({{count_var}} ~ group_var)$p.value, error = function(e) NA),
      mean = mean(count_var, na.rm = TRUE),
      sd = sd(count_var, na.rm = TRUE),
      ttest_vector = toString(count_var),
      count_var = sum(count_var),
      .groups = 'drop'
    ) %>%
    mutate(across(all_of(stats), ~replace(., is.na(.), 0))) %>%
    mutate(adjusted_p_value = p.adjust(ttest_p_value, method = "fdr")) %>%
    mutate(significance = case_when(
      adjusted_p_value < 0.0001 ~ "****",
      adjusted_p_value < 0.001 ~ "***",
      adjusted_p_value < 0.01 ~ "**",
      adjusted_p_value < 0.05 ~ "*",
      adjusted_p_value >= 0.05 ~ "ns",
      TRUE ~ "unknown"
    )) %>%
    filter(!is.na(ttest_p_value))%>%
    select({{group_var}},ttest_p_value, adjusted_p_value, significance, mean, ttest_vector, everything()) %>%
    arrange(ttest_p_value) %>%
    mutate(
      across(all_of(stats), ~signif(., digits = 3)),
      adjusted_p_value = p.adjust(ttest_p_value, method = "fdr"),
      ttest_p_value = sprintf("%.3e", ttest_p_value),
      adjusted_p_value = sprintf("%.2e", adjusted_p_value)
    )
  return(df)
}
################################################################################

hash=paste(c(rep('#',80)),collapse="")

hash_line<-function(line_symbol='#'){
  hash_line=paste0('#',paste(c(rep(line_symbol,78)),collapse=""),'#');return(hash_line)}

hash_text <- function(text_string='hi', line_symbol="#", use_spaces = TRUE, width = 80) {
  # Example:
  # hash_text("Load Libraries", line_symbol = "=", use_spaces = TRUE)
  
  # Build padding text
  if (use_spaces) {
    text <- paste0(" ", text_string, " ")
  } else {
    text <- text_string
  }
  
  # Compute lengths
  text_len <- nchar(text)
  remaining <- width - text_len - 2  # subtract 2 for starting and ending #
  
  # Left/right padding
  left_pad  <- floor(remaining / 2)
  right_pad <- ceiling(remaining / 2)
  
  # Build line WITH ending #
  banner <- paste0(
    "#",
    paste0(rep(line_symbol, left_pad), collapse = ""),
    text,
    paste0(rep(line_symbol, right_pad), collapse = ""),
    "#"
  )
  
  return(banner)
}


################################################################################

################################################################################
log_runtime <- function(start_time, end_time, function_name, log_file = "runtime_log.txt") {
  runtime <- end_time - start_time
  log_entry <- paste(Sys.time(), "Function:", function_name, "Runtime:", runtime, "seconds", "\n")

  # Append the log entry to the specified log file
  write(log_entry, file = log_file, append = TRUE)
}
################################################################################
# write file to csv (names vector, filepath, prefix)
################################################################################
qWrite<- function(dataframe_names, output_data= NA, p_title=NA) {
  #qWrite('loadings_long_df',output_data,p_title)
  #qWrite('loadings_long_df')
  # Retrieve the data frames from the environment
  dataframes <- mget(dataframe_names, envir = .GlobalEnv)
  if(is.na(output_data)){output_data <- mget('output_data', envir = .GlobalEnv)}
  if(is.na(p_title)){p_title <- mget('p_title', envir = .GlobalEnv)}

  # Iterate through the list of dataframes and their names
  for (name in dataframe_names) {
    df <- dataframes[[name]]
    filename <- paste0(output_data, p_title,'_',name, '.csv')
    write.csv(df, filename, row.names = FALSE)
  }
}
################################################################################
# read file from csv (names vector, filepath, prefix)
#qRead(dataframe_names, output_data, p_title)
################################################################################
qRead <- function(dataframe_names, output_data= NA, p_title=NA) {

  #qRead(dataframe_names)
  #dataframes <- mget(dataframe_names, envir = .GlobalEnv)
  if(is.na(output_data)){output_data <- mget('output_data', envir = .GlobalEnv)}
  if(is.na(p_title)){p_title <- mget('p_title', envir = .GlobalEnv)}
  dataframes <- list()

  # Iterate over the vector of names
  for (name in dataframe_names) {
    # Construct the filename based on the output directory and title
    filename <- paste0(output_data, p_title,'_',name, '.csv')

    # Read the CSV file into a data frame
    df <- read.csv(filename, stringsAsFactors = FALSE)

    # Add the data frame to the list with the name as the key
    dataframes[[name]] <- df
  }

  # Assign the data frames to the global environment
  list2env(dataframes, envir = .GlobalEnv)
}
################################################################################
qHead <- function(data, n_cols = 6, n_rows = 3, max_colname_length = 20) {
  # Check if input is a data frame or matrix
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Input must be a data frame or matrix.")
  }

  # Truncate column names and ensure uniqueness
  colnames(data) <- str_trunc(colnames(data), width = max_colname_length, side = "right", ellipsis = "...")

  # Ensure unique column names by adding a prefix if needed
  colnames(data) <- make.unique(colnames(data), sep = "_")

  # Print a subset of the data (up to n_cols and n_rows)
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }

  # Print the data
  print(head(data %>% select(1:n_cols), n_rows))
}
################################################################################
# mean centering a data frame by category sub groups
################################################################################
mean_normalizing_by_category <- function(df, category_col) {
  # df with columns = features rows = observations row.names=identifier
  data_cols <- names(df)

  # Apply xPlode_sample_name to remove rownames and add meta data tied to them
  exploded_df <- df %>% xPlode_sample_name()

  # Identify the meta_data columns and data columns
  meta_data_cols <- setdiff(names(exploded_df), data_cols)

  mean_centered_df <- exploded_df %>%
    split(.[[category_col]]) %>%
    lapply(function(sub_df) {
      # Separate meta_data and data columns
      sub_meta_data <- sub_df[, meta_data_cols, drop = FALSE]
      sub_data <- sub_df[, data_cols, drop = FALSE]

      # Calculate the column-wise (feature) means
      sub_mean_standard <- colMeans(sub_data)

      sub_mean_centered_data <- sub_data %>%
        t() %>%
        `/`(sub_mean_standard) %>%
        t() %>%
        as.data.frame() %>%
        mutate(across(everything(), ~ replace_na(.x, 0))) # handle 0 count columns


      # Combine meta_data with normalized data columns
      sub_mean_centered_df <- cbind(sub_meta_data, sub_mean_centered_data)

      return(sub_mean_centered_df)
    }) %>%
    do.call(rbind, .) %>%
    {rownames(.) <- NULL; .} %>%  # Remove row names
    imPlode_sample_name() # Apply imPlode_sample_name to reverse xPlode

  return(mean_centered_df)
}
################################################################################
# Define a function to perform NMDS with increasing k
nmds_with_increasing_k <- function(data, max_tries = 100, max_iterations = 999) {

  #example: results <- nmds_with_increasing_k(data)

  # Compute distance matrix
  dist_matrix <- vegdist(data, method = "bray")

  # Initialize variables to store results
  stress_values <- c()
  k_values <- c()

  # Start with k = 1 and increase until convergence or max_tries reached
  k <- 1
  converged <- FALSE

  while (!converged && k <= max_tries) {
    # Perform PCA for initialization
    qPrint(k)
    set.seed(123)
    pca_result <- prcomp(dist_matrix, scale. = FALSE)
    init <- pca_result$x[, 1:k]

    # Perform NMDS with current k
    nmds_result <- metaMDS(dist_matrix, k = k,
                           maxit = 999,
                           trymax = 250,
                           autotransform = FALSE, init = init,wascores = TRUE)

    # Check convergence
    if (nmds_result$converged) {
      converged <- TRUE
      # Store stress value and k
      stress_values <- c(stress_values, nmds_result$stress)
      k_values <- c(k_values, k)
    } else {
      # Store non-converged stress (NA) and k
      stress_values <- c(stress_values, NA)
      k_values <- c(k_values, k)
    }

    # Increase k
    k <- k + 1
    qPrint(k)
  }

  # Return a data frame with k values and corresponding stress
  results_df <- data.frame(K = k_values, Stress = stress_values)
  qPrint(results_df)
  write.csv(results_df,'nmds_k_results.csv')
  return(nmds_result)
}
################################################################################
# xAgg aggregate an exploded matrix by subsampling
################################################################################

xAgg <- function(dfx, aggregate_columns) {

  # Example usage: result_df <- xAgg(dfx, c('Lat_Location','Type','Species'))

  # Check if all specified columns exist in the data frame
  missing_cols <- aggregate_columns[!aggregate_columns %in% colnames(dfx)]

  # If there are missing columns, return an error
  if (length(missing_cols) > 0) {
    stop(paste("The following columns are not found in the data frame:", paste(missing_cols, collapse = ", ")))
  }

  # Remove excess metadata columns that are not part of the aggregation
  dfx <- dfx %>% select(-any_of(setdiff(meta_data, aggregate_columns)))

  # Get unique combinations of the aggregation columns
  group_combinations <- dfx %>%
    select(any_of(aggregate_columns)) %>%
    distinct()

  # Initialize an empty data frame to hold the summarized results
  dfx_summarized <- data.frame(matrix(ncol = ncol(dfx) - length(aggregate_columns),
                                      nrow = nrow(group_combinations)))
  colnames(dfx_summarized) <- names(dfx)[!names(dfx) %in% aggregate_columns]  # Set column names to ASV names

  # Loop through each unique group combination to calculate sums
  for (combination_row in 1:nrow(group_combinations)) {
    group_name <- paste(group_combinations[combination_row, ], collapse = '_')
    tPrint(group_name)

    # Subset the data frame to include only the current group combination
    ag_group <- group_combinations[combination_row, ] %>%
      left_join(dfx) %>%
      select(-all_of(aggregate_columns))

    # Calculate column sums for the current group and store them in the summarized data frame
    dfx_summarized[combination_row, ] <- colSums(ag_group, na.rm = TRUE)
  }

  # Combine the group combinations with their corresponding summarized data
  dfx <- cbind(group_combinations, dfx_summarized)

  return(dfx)  # Return the resulting data frame
}
###############################################################################
# xAgg_mean aggregate an exploded matrix by subsampling by mean
################################################################################

xAgg_mean <- function(dfx, aggregate_columns) {

  # Example usage: result_df <- xAgg(dfx, c('Lat_Location','Type','Species'))

  # Check if all specified columns exist in the data frame
  missing_cols <- aggregate_columns[!aggregate_columns %in% colnames(dfx)]

  # If there are missing columns, return an error
  if (length(missing_cols) > 0) {
    stop(paste("The following columns are not found in the data frame:", paste(missing_cols, collapse = ", ")))
  }

  # Remove excess metadata columns that are not part of the aggregation
  dfx <- dfx %>% select(-any_of(setdiff(meta_data, aggregate_columns)))

  # Get unique combinations of the aggregation columns
  group_combinations <- dfx %>%
    select(aggregate_columns) %>%
    distinct()

  # Initialize an empty data frame to hold the summarized results
  dfx_summarized <- data.frame(matrix(ncol = ncol(dfx) - length(aggregate_columns),
                                      nrow = nrow(group_combinations)))
  colnames(dfx_summarized) <- names(dfx)[!names(dfx) %in% aggregate_columns]  # Set column names to ASV names

  # Loop through each unique group combination to calculate sums
  for (combination_row in 1:nrow(group_combinations)) {
    group_name <- paste(group_combinations[combination_row, ], collapse = '_')
    tPrint(group_name)

    # Subset the data frame to include only the current group combination
    ag_group <- group_combinations[combination_row, ] %>%
      left_join(dfx) %>%
      select(-all_of(aggregate_columns))

    # Calculate column sums for the current group and store them in the summarized data frame
    dfx_summarized[combination_row, ] <- colMeans(ag_group, na.rm = TRUE)
  }

  # Combine the group combinations with their corresponding summarized data
  dfx <- cbind(group_combinations, dfx_summarized)

  return(dfx)  # Return the resulting data frame
}
##############################################################################
# generate simplified taxon labels
##############################################################################
feature_labels <- function(group = 'feature', df = matrix_df) {
  # Ensure 'feature' is included in the group list
  feature_names <- names(df)
  if (!('feature' %in% group)) {
    group <- c(group, 'feature')
  }

  # Function to generate the new feature label based on taxonomy
  generate_new_feature_label <- function(feature) {
    # Split the feature taxonomic chain into segments
    feature_split <- strsplit(feature, ";")[[1]]

    # Find the first occurrence of "__" (double underscore) and the last occurrence of "uncultured"
    first_double_underscore <- match(TRUE, feature_split == "__", nomatch = NA)
    last_uncultured <- if (any(grepl("uncultured", feature_split))) {
      max(which(grepl("uncultured", feature_split)))
    } else {
      NA
    }

    # Determine the minimum occurrence (either the first "__" or last "uncultured")
    min_occurrence <- ifelse(!is.na(first_double_underscore) | !is.na(last_uncultured),
                             pmin(first_double_underscore, last_uncultured, na.rm = TRUE),
                             NA)

    # Get the number of segments
    num_segments <- length(feature_split)

    # Create the new feature_label based on the conditions
    if (is.na(min_occurrence)) {
      return(tail(feature_split, 1))  # If min_occurrence is NA, return the last segment
    } else if (min_occurrence == 1) {
      return(paste(feature_split, collapse = ";"))  # If min_occurrence == 1, return the whole feature_label
    } else {
      # Otherwise, return the segments from (min_occurrence - 1) through the end
      return(paste(feature_split[(min_occurrence - 1):num_segments], collapse = ";"))
    }
  }

  # Process the input data frame
  df_processed <- df %>%
    as.data.frame() %>%
    xPlode_sample_name() %>%
    pivot_longer(cols = all_of(feature_names), names_to = 'feature', values_to = 'counts')

  if (taxa_levs == 'ASV') {
    df_processed <- df_processed %>%
      mutate(ASV=feature)%>%
      left_join(ASV_Taxa %>% select(ASV, Taxa), by = "ASV")
  } else {
    df_processed <- df_processed %>%
      mutate(Taxa = feature)
  }

  df_processed<-df_processed%>%
    mutate(
      Domain = sapply(Taxa, extract_domain),  # Assuming extract_domain() is defined elsewhere
      feature_label = sapply(Taxa, generate_new_feature_label),
      feature_label = if_else(feature_label == '__', 'unclassified', feature_label),
      feature_label = sub("^p__", "", feature_label)  # Clean up 'p__' prefix if needed
    ) %>%
    group_by(sample_name) %>%
    mutate(sample_counts = sum(counts)) %>%
    ungroup() %>%
    mutate(rel_abun = counts / sample_counts) %>%
    group_by(feature, feature_label, Domain) %>%
    mutate(group_rel_abun = rel_abun) %>%
    group_by(across(all_of(group)), feature_label, Domain,Taxa) %>%
    summarize(group_rel_abun = mean(rel_abun), .groups = "drop") %>%
    ungroup() %>%
    mutate(
      rel_abun_label = paste0(signif(100 * group_rel_abun, 2), '%'),
      feature_label2 = paste0("<span style='color:", color_palette[Domain], "'><i>", feature_label, "</i> ", rel_abun_label, "</span>"),
      feature_label3 = paste0(feature_label, ' ', signif(100 * group_rel_abun, 2), '%')
    ) %>%
    arrange(desc(group_rel_abun)) %>%
    mutate(plot_order = factor(feature_label3, levels = rev(unique(feature_label3))))

  return(df_processed)
}
Kruskal_Mann_Whitney_Test <- function(df = matrix_df, testing_group = "primer_Label", category_group = "Type") {

  # example of use:  kw_group_results=Kruskal_Mann_Whitney_Test(matrix_df,plotting_set)

  pseudo_count <- min(df[df > 0], na.rm = TRUE) / 2
  pseudo_df = df+pseudo_count
  relative_abundance=pseudo_df/rowSums(pseudo_df)
  clr_df=clr(relative_abundance)

  kw_df <- clr_df %>%
    xPlode_sample_name()

  kw_group_results <- data.frame()

  # Separate data by the category group (e.g., "Type")
  unique_categories <- unique(kw_df[[category_group]])

  for (category in unique_categories) {

    category_data <- df %>%
      xPlode_sample_name() %>%
      filter(.data[[category_group]] == category) %>%
      imPlode_sample_name() %>%
      select(where(~ sum(.) != 0)) %>%
      t() %>% as.data.frame() %>%
      mutate(across(everything(), ~ . / sum(.))) %>%
      t() %>%
      clr() %>%
      as.data.frame() %>%
      xPlode_sample_name()

    # Check number of unique categories in the testing group
    check <- category_data %>%
      select(.data[[testing_group]]) %>%
      distinct() %>%
      pull()

    cat("Number of unique categories in", testing_group, ":", length(check), "\n")
    print(check)

    # Determine the test to use based on the number of categories in testing_group
    if (length(check) == 2) {
      test_type <- "Mann-Whitney"
      test_func <- wilcox.test
    } else if (length(check) > 2) {
      test_type <- "Kruskal-Wallis"
      test_func <- kruskal.test
    } else {
      warning(paste("Skipping category", category, "due to insufficient unique categories"))
      next
    }

    # Perform the chosen test for each group
    for (group in testing_group) {
      group_var <- category_data[[group]]  # Extract the grouping variable for the category

      # Skip if the grouping variable has less than two unique values
      if (length(unique(group_var)) > 1) {
        test_results <- category_data %>%
          select(-any_of(meta_data)) %>%  # Exclude the grouping variable(s)
          map_df(~ broom::tidy(test_func(.x ~ group_var)), .id = "feature") %>%
          mutate(testing = group,
                 category = category,
                 test_type = test_type) %>%
          ungroup()

        # Adjust p-values regardless of the test type
        test_results <- test_results %>%
          mutate(adjusted_pvalue = p.adjust(p.value, method = "BH"))  # Apply BH adjustment to p-value

        # Add significance for both tests
        test_results <- test_results %>%
          mutate(significance = case_when(
            p.value < 0.001 ~ "***",  # Highly significant
            p.value < 0.01  ~ "**",   # Significant
            p.value < 0.05  ~ "*",    # Marginally significant
            p.value < 0.1   ~ ".",    # Suggestive
            TRUE            ~ "ns"    # Not significant
          )) %>%
          mutate(adjusted_significance = case_when(
            adjusted_pvalue < 0.001 ~ "***",  # Highly significant
            adjusted_pvalue < 0.01  ~ "**",   # Significant
            adjusted_pvalue < 0.05  ~ "*",    # Marginally significant
            adjusted_pvalue < 0.1   ~ ".",    # Suggestive
            TRUE            ~ "ns"    # Not significant
          )) %>%
          mutate(adj_check = adjusted_pvalue / p.value)%>%
          rename(!!category_group := category)

        # Add the results to the final dataframe
        kw_group_results <- bind_rows(kw_group_results, test_results)
      } else {
        warning(paste("Skipping group", group, "in category", category, "because it has only one unique category."))
      }
    }
  }

  return(kw_group_results)
}
################################################################################
# end?
################################################################################

################################################################################
################################################################################
#slice(chull(NMDS1, NMDS2))
