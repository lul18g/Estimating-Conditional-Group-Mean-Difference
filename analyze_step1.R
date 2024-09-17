################################################################################ 
# step 1: obtain composite scores, factor scores, residual variance based on original data # 
################################################################################ 
install.packages("lavaan")
install.packages("psych")
library(lavaan)
library(psych)

#### unidimensional measurement model ####
id = "uni" 
if (id == "uni"){
  ## specify simulation conditions
  ld = c(85,63,455,657) #loading pattern: 85-(.8,.5);63-(.6,.3);455-(.455);657-(.657)
  d = c(0,2,5) #conditional group mean difference
  nn = c(50,100,150) #sample size
  
  cond = expand.grid(ld,d,nn)
  condidx = as.numeric(commandArgs(trailingOnly=TRUE))
  
  ld=cond[condidx,1]
  d=cond[condidx,2]
  nn=cond[condidx,3]
  
  ## create data folder for subsequent analysis for each condition
  dir = paste0('/gpfs/home/CGMD/step1analysis/','L',ld,'_d',d,'_nn',nn,'/')
  dir.create(dir)
  
  ## set working directory
  datafolder=paste0('/gpfs/home/CGMD/datagen/','L',ld,'_d',d,'_nn',nn,'/')
  setwd(datafolder)
  
  results=as.data.frame(matrix(ncol=49,nrow=1))
  colnames(results) <- c('rep', paste0('L',1:6,'_x'),paste0('se.L',1:6,'_x'),
                         paste0('e',1:6,'_x'),paste0('se.e',1:6,'_x'),
                         paste0('L',1:6,'_y'),paste0('se.L',1:6,'_y'),
                         paste0('e',1:6,'_y'),paste0('se.e',1:6,'_y'))
  
  for (r in 1:2000) {
    dataname=paste0('rep',r,'.dat')
    dat = read.table(dataname)
    colnames(dat)=c("x1","x2","x3","x4","x5","x6","y1","y2","y3","y4","y5","y6","g")
    
    ## add unweighted mean scores to dat
    dat$uwmeanx=rowSums(dat[,c("x1","x2","x3","x4","x5","x6")])/6
    dat$uwmeany=rowSums(dat[,c("y1","y2","y3","y4","y5","y6")])/6
    
    
    ## fit CFA models for Fx and Fy, respectively, to obtain parameter estimates and factor scores needed for subsequent analyses
    if (ld == 85){
      #CFA models
      cfa_x='fx =~ start(.8)*x1+start(.5)*x2+start(.8)*x3+start(.5)*x4+start(.8)*x5+start(.5)*x6
         fx~~1*fx'#UVI
      cfa_y='fy =~ start(.8)*y1+start(.5)*y2+start(.8)*y3+start(.5)*y4+start(.8)*y5+start(.5)*y6
         fy~~1*fy'
    }
    if (ld == 63){
      #CFA models
      cfa_x='fx =~ start(.6)*x1+start(.3)*x2+start(.6)*x3+start(.3)*x4+start(.6)*x5+start(.3)*x6
         fx~~1*fx'#UVI
      cfa_y='fy =~ start(.6)*y1+start(.3)*y2+start(.6)*y3+start(.3)*y4+start(.6)*y5+start(.3)*y6
         fy~~1*fy'
    }
    if (ld == 455){
      #CFA models
      cfa_x='fx =~ start(.455)*x1+start(.455)*x2+start(.455)*x3+start(.455)*x4+start(.455)*x5+start(.455)*x6
         fx~~1*fx'#UVI
      cfa_y='fy =~ start(.455)*y1+start(.455)*y2+start(.455)*y3+start(.455)*y4+start(.455)*y5+start(.455)*y6
         fy~~1*fy'
    }
    if (ld == 657){
      #CFA models
      cfa_x='fx =~ start(.657)*x1+start(.657)*x2+start(.657)*x3+start(.657)*x4+start(.657)*x5+start(.657)*x6
         fx~~1*fx'#UVI
      cfa_y='fy =~ start(.657)*y1+start(.657)*y2+start(.657)*y3+start(.657)*y4+start(.657)*y5+start(.657)*y6
         fy~~1*fy'
    }
    
    fit_cfa_x = sem(cfa_x,data=dat,estimator="ML",std.lv=T,mimic="Mplus") #fit CFA
    fit_cfa_y = sem(cfa_y,data=dat,estimator="ML",std.lv=T,mimic="Mplus")
    
    
    ## obtain unstandardized solutions
    uslt_cfa_x=parameterEstimates(fit_cfa_x)
    uslt_cfa_y=parameterEstimates(fit_cfa_y)
    
    res = as.data.frame(matrix(ncol=49,nrow=1))
    res[1,1]=r
    
    ## extract estimates of interest
    # cfa_x
    ld_x=data.frame(t(uslt_cfa_x[1:6,c("est")])) #unstandardized factor loadings
    colnames(ld_x)=c("L1_x","L2_x","L3_x","L4_x","L5_x","L6_x")
    ldx=rowSums(ld_x)
    res[1,2:7]=ld_x[1,1:6]
    
    res[1,8:13]=uslt_cfa_x[1:6,c("se")] #standard error of unstandardized factor loading estimates
    
    rv_x=data.frame(t(uslt_cfa_x[8:13,c("est")])) # unstandardized residual variance
    rvx=rowSums(rv_x) #sum of residual variance across items 
    res[1,14:19]=rv_x[1,1:6]
    
    res[1,20:25]=uslt_cfa_x[8:13,c("se")] # standard error of unstandardized residual variance estimates
    
    #cfa_y
    ld_y=data.frame(t(uslt_cfa_y[1:6,c("est")]))
    colnames(ld_y)=c("L1_y","L2_y","L3_y","L4_y","L5_y","L6_y")
    ldy=rowSums(ld_y)
    res[1,26:31]=ld_y[1,1:6]
    
    res[1,32:37]=uslt_cfa_y[1:6,c("se")]
    
    rv_y=data.frame(t(uslt_cfa_y[8:13,c("est")])) 
    rvy=rowSums(rv_y) 
    res[1,38:43]=rv_y[1,1:6]
    
    res[1,44:49]=uslt_cfa_y[8:13,c("se")]
    
    results[r,]=res
    
    ## compute weighted mean scores for Fx and Fy, respectively (using unstandarsized loadings as the weights)
    dat$wmeanx=(dat$x1*ld_x$L1_x+dat$x2*ld_x$L2_x+dat$x3*ld_x$L3_x+dat$x4*ld_x$L4_x+dat$x5*ld_x$L5_x+dat$x6*ld_x$L6_x)/6
    dat$wmeany=(dat$y1*ld_y$L1_y+dat$y2*ld_y$L2_y+dat$y3*ld_y$L3_y+dat$y4*ld_y$L4_y+dat$y5*ld_y$L5_y+dat$y6*ld_y$L6_y)/6
    
    ## obtain factor scores
    fx=data.frame(lavPredict(fit_cfa_x,newdata=dat,type = "lv", method = "regression",append.data=F))
    colnames(fx)=c("fsx")
    
    fy=data.frame(lavPredict(fit_cfa_y,newdata=dat,type = "lv", method = "regression",append.data=F))
    colnames(fy)=c("fsy")
    
    ## add factor scores to dat
    dat=data.frame(dat,fx,fy)
    
    ## compute reliability
    # alpha
    dat_x=dat[,c("x1","x2","x3","x4","x5","x6")]
    dat_y=dat[,c("y1","y2","y3","y4","y5","y6")]
    
    alpx=alpha(dat_x)
    alpha_x=alpx$total$raw_alpha
    dat$alp_x=alpha_x
    
    alpy=alpha(dat_y)
    alpha_y=alpy$total$raw_alpha
    dat$alp_y=alpha_y
    
    ## obtain measurement error variance that is used to adjust unreliability for single indicator (i.e., unweighted mean scores)
    dat$evx_alp=var(dat$uwmeanx)*(1-alpha_x) #error variance of mean scores computed via coefficient alpha
    dat$evy_alp=var(dat$uwmeany)*(1-alpha_y)
    
    ## code proper solutions
    results$proper_x = ifelse(results$e1_x< 0|results$e2_x< 0|results$e3_x< 0|results$e4_x< 0|
                                results$e5_x< 0|results$e6_x< 0,0,1)
    
    results$proper_y = ifelse(results$e1_y< 0|results$e2_y< 0|results$e3_y< 0|results$e4_y< 0|
                                results$e5_y< 0|results$e6_y< 0,0,1)
    
    if (results$proper_x[r] == 0){
      dat$wmeanx = 9999
      dat$fsx = 9999
    }
    
    if (results$proper_y[r] == 0){
      dat$wmeany = 9999
      dat$fsy = 9999
    }
    
    data = dat[,c("g","uwmeanx","uwmeany","wmeanx","wmeany","fsx","fsy","alp_x","alp_y",
                  "evx_alp","evy_alp")]
    ## export data for subsequent analysis
    filename=paste0(dir,'rep',r,'.dat')
    write.table(data,file=filename,sep="\t",col.names = FALSE,row.names = FALSE)
  } 
  
  results$loading=ld
  results$effect=d
  results$nn=nn
  resname=paste0('/gpfs/home/CGMD/step1analysis/results_','L',ld,'_d',d,'_nn',nn,'.csv')
  write.csv(results,file=resname)
}


#### multidimensional measurement model ####
id = "multi"
if (id == "multi") {
  ## specify simulation conditions
  gld = c(543,765) #group-specific factor loadings: 543-(.5,.4,.3);765-(.7,.6,.5)
  ld = c(3456,5658) #general factor loadings: 3456-(.3,.45,.6);5658-(.5,.65,.8)
  d = c(0,2,5) #conditional group mean difference
  nn = c(50,100,150) #sample size
  
  cond = expand.grid(gld,ld,d,nn)
  condidx = as.numeric(commandArgs(trailingOnly=TRUE))
  
  gld = cond[condidx,1]
  ld = cond[condidx,2]
  d = cond[condidx,3]
  nn = cond[condidx,4]
  
  ## create data folder for subsequent analysis for each condition
  dir = paste0('/gpfs/home/CGMD/step1analysis/','gL',gld,'_L',ld,'_d',d,'_nn',nn,'/')
  dir.create(dir)
  
  ## set working directory
  datafolder=paste0('/gpfs/home/CGMD/datagen/','gL',gld,'_L',ld,'_d',d,'_nn',nn,'/')
  setwd(datafolder)
  
  
  results1=as.data.frame(matrix(ncol=51,nrow=1)) #results from one-factor CFA
  colnames(results1) <- c('rep', paste0('L',1:6,'_x'),paste0('se.L',1:6,'_x'),
                          paste0('e',1:6,'_x'),paste0('se.e',1:6,'_x'),
                          paste0('L',1:6,'_y'),paste0('se.L',1:6,'_y'),
                          paste0('e',1:6,'_y'),paste0('se.e',1:6,'_y'),
                          'proper_x','proper_y')
  results2=as.data.frame(matrix(ncol=63,nrow=1)) #results from bi-factor model
  colnames(results2) <- c('rep', paste0('fL',1:6,'_x'),paste0('gL',1:3,'_x'),
                          paste0('se.fL',1:6,'_x'),paste0('se.gL',1:3,'_x'),
                          paste0('e',1:6,'_x'),paste0('se.e',1:6,'_x'),
                          paste0('fL',1:6,'_y'),paste0('gL',1:3,'_y'),
                          paste0('se.fL',1:6,'_y'),paste0('se.gL',1:3,'_y'),
                          paste0('e',1:6,'_y'),paste0('se.e',1:6,'_y'),
                          "proper_x","proper_y")
  
  for (r in 1:2000) {
    dataname=paste0('rep',r,'.dat')
    dat = read.table(dataname)
    colnames(dat)=c("x1","x2","x3","x4","x5","x6","y1","y2","y3","y4","y5","y6","g")
    
    ## add unweighted mean scores to dat
    dat$uwmeanx=rowSums(dat[,c("x1","x2","x3","x4","x5","x6")])/6
    dat$uwmeany=rowSums(dat[,c("y1","y2","y3","y4","y5","y6")])/6
    
    
    ## fit one-factor CFA models for Fx and Fy, respectively, to obtain parameter estimates and factor scores needed for subsequent analyses
    if (ld == 3456){
      #one-factor CFA models
      cfa1_x='fx =~ start(.3)*x1+start(.45)*x2+start(.6)*x3+start(.3)*x4+start(.45)*x5+start(.6)*x6
         fx~~1*fx'#UVI
      cfa1_y='fy =~ start(.3)*y1+start(.45)*y2+start(.6)*y3+start(.3)*y4+start(.45)*y5+start(.6)*y6
         fy~~1*fy'
    }
    if (ld == 5658){
      #one-factor CFA models
      cfa1_x='fx =~ start(.5)*x1+start(.65)*x2+start(.8)*x3+start(.5)*x4+start(.65)*x5+start(.8)*x6
         fx~~1*fx'#UVI
      cfa1_y='fy =~ start(.5)*y1+start(.65)*y2+start(.8)*y3+start(.5)*y4+start(.65)*y5+start(.8)*y6
         fy~~1*fy'
    }
    
    fit_cfa1_x = sem(cfa1_x,data=dat,estimator="ML",std.lv=T,mimic="Mplus") #fit CFA
    fit_cfa1_y = sem(cfa1_y,data=dat,estimator="ML",std.lv=T,mimic="Mplus")
    
    
    ## obtain unstandardized solutions
    uslt_cfa1_x=parameterEstimates(fit_cfa1_x)
    uslt_cfa1_y=parameterEstimates(fit_cfa1_y)
    
    res1 = as.data.frame(matrix(ncol=49,nrow=1))
    colnames(res1) <- c('rep', paste0('L',1:6,'_x'),paste0('se.L',1:6,'_x'),
                        paste0('e',1:6,'_x'),paste0('se.e',1:6,'_x'),
                        paste0('L',1:6,'_y'),paste0('se.L',1:6,'_y'),
                        paste0('e',1:6,'_y'),paste0('se.e',1:6,'_y'))
    res1[1,1]=r
    
    ## extract estimates of interest
    # cfa_x
    ld1_x=data.frame(t(uslt_cfa1_x[1:6,c("est")])) #unstandardized factor loadings
    colnames(ld1_x)=c("L1_x","L2_x","L3_x","L4_x","L5_x","L6_x")
    ld1x=rowSums(ld1_x)
    res1[1,2:7]=ld1_x[1,1:6]
    
    res1[1,8:13]=uslt_cfa1_x[1:6,c("se")] #standard error of unstandardized factor loading estimates
    
    rv1_x=data.frame(t(uslt_cfa1_x[8:13,c("est")])) #unstandardized residual variance
    rv1x=rowSums(rv1_x) #sum of residual variance across items 
    res1[1,14:19]=rv1_x[1,1:6]
    
    res1[1,20:25]=uslt_cfa1_x[8:13,c("se")] #standard error of unstandardized residual variance estimates
    
    
    #cfa_y
    ld1_y=data.frame(t(uslt_cfa1_y[1:6,c("est")]))
    colnames(ld1_y)=c("L1_y","L2_y","L3_y","L4_y","L5_y","L6_y")
    ld1y=rowSums(ld1_y)
    res1[1,26:31]=ld1_y[1,1:6]
    
    res1[1,32:37]=uslt_cfa1_y[1:6,c("se")]
    
    rv1_y=data.frame(t(uslt_cfa1_y[8:13,c("est")])) 
    rv1y=rowSums(rv1_y) 
    res1[1,38:43]=rv1_y[1,1:6]
    
    res1[1,44:49]=uslt_cfa1_y[8:13,c("se")]
    
    ## fit bi-factor models for covariate and outcome, respectively
    if (ld == 3456 & gld == 543){
      cfa2_x='fx =~ start(.3)*x1+start(.45)*x2+start(.6)*x3+start(.3)*x4+start(.45)*x5+start(.6)*x6
         gx =~ start(.5)*x1+start(.4)*x2+start(.3)*x3
         fx~~1*fx
         gx~~1*gx
         fx~~0*gx'
      
      cfa2_y='fy =~ start(.3)*y1+start(.45)*y2+start(.6)*y3+start(.3)*y4+start(.45)*y5+start(.6)*y6
         gy =~ start(.5)*y1+start(.4)*y2+start(.3)*y3
         fy~~1*fy
         gy~~1*gy
         fy~~0*gy'
    }
    if (ld == 3456 & gld == 765){
      cfa2_x='fx =~ start(.3)*x1+start(.45)*x2+start(.6)*x3+start(.3)*x4+start(.45)*x5+start(.6)*x6
         gx =~ start(.7)*x1+start(.6)*x2+start(.5)*x3
         fx~~1*fx
         gx~~1*gx
         fx~~0*gx'
      
      cfa2_y='fy =~ start(.3)*y1+start(.45)*y2+start(.6)*y3+start(.3)*y4+start(.45)*y5+start(.6)*y6
         gy =~ start(.7)*y1+start(.6)*y2+start(.5)*y3
         fy~~1*fy
         gy~~1*gy
         fy~~0*gy'
    }
    if (ld == 5658 & gld == 543){
      cfa2_x='fx =~ start(.5)*x1+start(.65)*x2+start(.8)*x3+start(.5)*x4+start(.65)*x5+start(.8)*x6
         gx =~ start(.5)*x1+start(.4)*x2+start(.3)*x3
         fx~~1*fx
         gx~~1*gx
         fx~~0*gx'
      
      cfa2_y='fy =~ start(.5)*y1+start(.65)*y2+start(.8)*y3+start(.5)*y4+start(.65)*y5+start(.8)*y6
         gy =~ start(.5)*y1+start(.4)*y2+start(.3)*y3
         fy~~1*fy
         gy~~1*gy
         fy~~0*gy'
    }
    if (ld == 5658 & gld == 765){
      cfa2_x='fx =~ start(.5)*x1+start(.65)*x2+start(.8)*x3+start(.5)*x4+start(.65)*x5+start(.8)*x6
         gx =~ start(.7)*x1+start(.6)*x2+start(.5)*x3
         fx~~1*fx
         gx~~1*gx
         fx~~0*gx'
      
      cfa2_y='fy =~ start(.5)*y1+start(.65)*y2+start(.8)*y3+start(.5)*y4+start(.65)*y5+start(.8)*y6
         gy =~ start(.7)*y1+start(.6)*y2+start(.5)*y3
         fy~~1*fy
         gy~~1*gy
         fy~~0*gy'
    }
    
    fit_cfa2_x = sem(cfa2_x,data=dat,estimator="ML",std.lv=T,mimic="Mplus")
    fit_cfa2_y = sem(cfa2_y,data=dat,estimator="ML",std.lv=T,mimic="Mplus")
    
    
    ## obtain unstandardized solutions
    uslt_cfa2_x=parameterEstimates(fit_cfa2_x)
    uslt_cfa2_y=parameterEstimates(fit_cfa2_y)
    
    res2 = as.data.frame(matrix(ncol=61,nrow=1))
    colnames(res2) <- c('rep', paste0('fL',1:6,'_x'),paste0('gL',1:3,'_x'),
                        paste0('se.fL',1:6,'_x'),paste0('se.gL',1:3,'_x'),
                        paste0('e',1:6,'_x'),paste0('se.e',1:6,'_x'),
                        paste0('fL',1:6,'_y'),paste0('gL',1:3,'_y'),
                        paste0('se.fL',1:6,'_y'),paste0('se.gL',1:3,'_y'),
                        paste0('e',1:6,'_y'),paste0('se.e',1:6,'_y'))
    res2[1,1] = r
    ## extract estimates of interest
    #cfa2_x
    ld_fx=data.frame(t(uslt_cfa2_x[1:6,c("est")])) #unstandardized factor loadings for the general factors
    colnames(ld_fx)=c("L1_fx","L2_fx","L3_fx","L4_fx","L5_fx","L6_fx")
    ldfx=rowSums(ld_fx)
    res2[1,2:7]=ld_fx[1,1:6]
    
    ld_gx=data.frame(t(uslt_cfa2_x[7:9,c("est")])) #unstandardized factor loadings for the group-specific factors
    ldgx=rowSums(ld_gx)
    res2[1,8:10]=ld_gx[1,1:3]
    
    res2[1,11:19]=uslt_cfa2_x[1:9,c("se")] #standard error of unstandardized factor loading estimates
    
    
    rv2_x=data.frame(t(uslt_cfa2_x[13:18,c("est")])) #unstandardized residual variance
    colnames(rv2_x)=c('e1','e2','e3','e4','e5','e6')
    rv2x=rowSums(rv2_x) #sum of residual variance across items 
    res2[1,20:25] = rv2_x[1,1:6]
    
    res2[1,26:31]=uslt_cfa2_x[13:18,c("se")] #standard error of unstandardized residual variance estimates
    
    
    #cfa2_y
    ld_fy=data.frame(t(uslt_cfa2_y[1:6,c("est")]))
    colnames(ld_fy)=c("L1_fy","L2_fy","L3_fy","L4_fy","L5_fy","L6_fy")
    ldfy=rowSums(ld_fy)
    res2[1,32:37]=ld_fy[1,1:6]
    
    ld_gy=data.frame(t(uslt_cfa2_y[7:9,c("est")]))
    ldgy=rowSums(ld_gy)
    res2[1,38:40]=ld_gy[1,1:3]
    
    res2[1,41:49]=uslt_cfa2_y[1:9,c("se")]
    
    rv2_y=data.frame(t(uslt_cfa2_y[13:18,c("est")])) 
    colnames(rv2_y)=c('e1','e2','e3','e4','e5','e6')
    rv2y=rowSums(rv2_y)  
    res2[1,50:55] = rv2_y[1,1:6]
    
    res2[1,56:61]=uslt_cfa2_y[13:18,c("se")]
    
    
    ## compute weighted mean scores for Fx and Fy, respectively (using unstandardized loadings based on one-factor cfa as the weights)
    dat$wmeanx=(dat$x1*ld1_x$L1_x+dat$x2*ld1_x$L2_x+dat$x3*ld1_x$L3_x+dat$x4*ld1_x$L4_x+dat$x5*ld1_x$L5_x+dat$x6*ld1_x$L6_x)/6
    dat$wmeany=(dat$y1*ld1_y$L1_y+dat$y2*ld1_y$L2_y+dat$y3*ld1_y$L3_y+dat$y4*ld1_y$L4_y+dat$y5*ld1_y$L5_y+dat$y6*ld1_y$L6_y)/6
    
    ## obtain factor scores (based on one-factor cfa)
    fx=data.frame(lavPredict(fit_cfa1_x,newdata=dat,type = "lv", method = "regression",append.data=F))
    colnames(fx)=c("fsx")
    
    fy=data.frame(lavPredict(fit_cfa1_y,newdata=dat,type = "lv", method = "regression",append.data=F))
    colnames(fy)=c("fsy")
    
    ## add factor scores to dat
    dat=data.frame(dat,fx,fy)
    
    ## compute reliability
    # alpha
    dat_x=dat[,c("x1","x2","x3","x4","x5","x6")]
    dat_y=dat[,c("y1","y2","y3","y4","y5","y6")]
    
    alpx=alpha(dat_x)
    alpha_x=alpx$total$raw_alpha
    dat$alp_x=alpha_x
    
    alpy=alpha(dat_y)
    alpha_y=alpy$total$raw_alpha
    dat$alp_y=alpha_y
    
    ## omega hierarchical
    true_fx=ldfx^2
    true_gx=ldgx^2
    error_x=rv2x
    omegaH_x=true_fx/(true_fx+true_gx+error_x)
    dat$omgH_x=omegaH_x
    
    true_fy=ldfy^2
    true_gy=ldgy^2
    error_y=rv2y
    omegaH_y=true_fy/(true_fy+true_gy+error_y)
    dat$omgH_y=omegaH_y
    
    ## obtain measurement error variance that is used to adjust unreliability for single indicator (i.e., unweighted mean scores)
    dat$evx_alp=var(dat$uwmeanx)*(1-alpha_x) #error variance of mean scores computed via coefficient alpha
    dat$evy_alp=var(dat$uwmeany)*(1-alpha_y)
    
    dat$evx_omgH=var(dat$uwmeanx)*(1-omegaH_x) #error variance of mean scores computed via coefficient omegaH
    dat$evy_omgH=var(dat$uwmeany)*(1-omegaH_y)
    
    ## code proper solutions
    res1$proper_x = ifelse(res1$e1_x< 0|res1$e2_x< 0|res1$e3_x< 0|res1$e4_x< 0|
                             res1$e5_x< 0|res1$e6_x< 0,0,1)
    
    res1$proper_y = ifelse(res1$e1_y< 0|res1$e2_y< 0|res1$e3_y< 0|res1$e4_y< 0|
                             res1$e5_y< 0|res1$e6_y< 0,0,1)
    
    res2$proper_x = ifelse(res2$e1_x< 0|res2$e2_x< 0|res2$e3_x< 0|res2$e4_x< 0|
                             res2$e5_x< 0|res2$e6_x< 0,0,1)
    
    res2$proper_y = ifelse(res2$e1_y< 0|res2$e2_y< 0|res2$e3_y< 0|res2$e4_y< 0|
                             res2$e5_y< 0|res2$e6_y< 0,0,1)
    
    if (res1$proper_x == 0){
      dat$wmeanx = 9999
      dat$fsx = 9999
    }
    
    if (res1$proper_y == 0){
      dat$wmeany = 9999
      dat$fsy = 9999
    }
    
    if (res2$proper_x == 0){
      dat$omgH_x = 9999
      dat$evx_omgH = 9999
    }
    
    if (res2$proper_y == 0){
      dat$omgH_y = 9999
      dat$evy_omgH = 9999
    }
    
    results1[r,] = res1
    results2[r,] = res2
    
    data = dat[,c("g","uwmeanx","uwmeany","wmeanx","wmeany","fsx","fsy","alp_x","alp_y",
                  "omgH_x","omgH_y","evx_alp","evy_alp","evx_omgH","evy_omgH")]
    ## export data for subsequent analysis
    filename=paste0(dir,'rep',r,'.dat')
    write.table(data,file=filename,sep="\t",col.names = FALSE,row.names = FALSE)
  } 
  
  results1$gloading=gld
  results1$loading=ld
  results1$effect=d
  results1$nn=nn
  resname1=paste0('/gpfs/home/CGMD/step1analysis/results_1factor_','gL',gld,'_L',ld,'_d',d,'_nn',nn,'.csv')
  write.csv(results1,file=resname1)
  
  results2$gloading=gld
  results2$loading=ld
  results2$effect=d
  results2$nn=nn
  resname2=paste0('/gpfs/home/CGMD/step1analysis/results_bifactor_','gL',gld,'_L',ld,'_d',d,'_nn',nn,'.csv')
  write.csv(results2,file=resname2)
}
