# 03_hypothesis_testing.R — RUN ONE HYPOTHESIS SECTION AT A TIME
# Run after 01_data_preparation_revised.R in the folder containing the .rds files.
# Requires: glmnet, ranger, xgboost, sandwich, lmtest (dplyr for R5). Optional for RF SHAP: fastshap.
# Run SETUP once, then H1, H2, H3, H4, H5a and H5b.
# AFTER the main run, optionally run H2 INCREMENTAL and robustness R1-R5.
# Run on the generated RDS files; verify the outputs and package versions.
# Main hypotheses use data_main.rds; robustness scenarios are at the end.


setwd("~/Desktop/Thesis_ESG_data/Data")
required <- c('glmnet','ranger','xgboost','sandwich','lmtest')
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing_pkgs)) stop('Install required packages: ', paste(missing_pkgs,collapse=', '))
SEED <- 20261027L
set.seed(SEED)
B <- 1000L
OUT <- 'results_hypotheses'
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
ind <- readRDS('ind_class.rds')
vars <- ind$var
miss_vars <- paste0('M_', ind$var[ind$make_M])
score <- 'esg_score'
pillars <- c('e_pillar','s_pillar','g_pillar')
ctl <- c('Size','Trend','PriorControversy')
folds <- list(c(2015,2017,2018),c(2015,2018,2019),c(2015,2019,2020))
# AUC: rank-based ROC-AUC; PR-AUC: non-interpolated average precision (AP).
roc_auc <- function(y,p) {
  a <- sum(y==1); b <- sum(y==0)
  if (a==0 || b==0) return(NA_real_)
  (sum(rank(p,ties.method='average')[y==1])-a*(a+1)/2)/(a*b)
}
pr_auc <- function(y,p) {
  n <- sum(y==1)
  if (!n) return(NA_real_)
  o <- order(p,decreasing=TRUE); yy <- y[o]
  # group tied scores to avoid arbitrary order within ties
  ends <- which(!duplicated(p[o],fromLast=TRUE))
  tp <- cumsum(yy)[ends]; k <- ends
  sum((tp-c(0,head(tp,-1)))/n * (tp/k))
}
f1 <- function(y,p,cut) {
  pos <- p >= cut; tp <- sum(pos & y==1); fp <- sum(pos & y==0); fn <- sum(!pos & y==1)
  if ((2*tp+fp+fn)==0) return(0)
  2*tp/(2*tp+fp+fn)
}
choose_cut <- function(y,p) {
  candidates <- sort(unique(c(0,1,p)))
  vals <- vapply(candidates,function(cut) f1(y,p,cut),numeric(1))
  candidates[which.max(vals)]
}
# Sector encoding fixed across folds; one level is reference for logistic regression.
make_x <- function(d, spec, sector_levels, controls=ctl) {
  d$sector25 <- factor(d$sector25, levels=sector_levels)
  sec <- model.matrix(~sector25, d)[,-1,drop=FALSE]
  colnames(sec) <- make.names(colnames(sec), unique = TRUE)
  esg <- switch(spec, M0=character(), M1a=score, M1b=pillars,
                M2=vars, M3=miss_vars, M4=c(vars,miss_vars),
                M2a=c(score,vars), M2b=c(pillars,vars),
                stop('Unknown specification ',spec))
  x <- cbind(as.matrix(d[,controls,drop=FALSE]), sec,
             as.matrix(d[,esg,drop=FALSE]))
  storage.mode(x) <- 'double'
  list(x=x, esg=esg, controls=setdiff(colnames(x),esg))
}
# The ESG predictors (INCLUDING disclosure dummies) are standardised within
# each training fold. Controls and sector reference coding remain unchanged.
scale_esg <- function(train, val, esg) {
  if (!length(esg)) return(list(tr=train,va=val))
  means <- colMeans(train[,esg,drop=FALSE]); sds <- apply(train[,esg,drop=FALSE],2,sd)
  sds[!is.finite(sds) | sds==0] <- 1
  train[,esg] <- sweep(sweep(train[,esg,drop=FALSE],2,means,'-'),2,sds,'/')
  val[,esg]   <- sweep(sweep(val[,esg,drop=FALSE],2,means,'-'),2,sds,'/')
  list(tr=train,va=val)
}
fit_predict <- function(train, val, spec, method, pars=list(), levels_sector, controls=ctl) {
  xx <- make_x(train,spec,levels_sector,controls)
  vv <- make_x(val,spec,levels_sector,controls)
  x <- xx$x; newx <- vv$x
  if (method %in% c('LogR','EN')) {
    if (anyNA(x) || anyNA(newx)) stop('Linear models require imputed indicators')
    scaled <- scale_esg(x,newx,xx$esg); x <- scaled$tr; newx <- scaled$va
  }
  y <- as.integer(train$Y)
  if (method=='LogR') {
    # Explicit intercept, no sector interaction terms; singular columns handled by glm.
    dat <- as.data.frame(x, check.names=FALSE); dat$Y <- y
    fit <- suppressWarnings(glm(Y ~ .,data=dat,family=binomial()))
    p <- suppressWarnings(as.numeric(predict(fit,newdata=as.data.frame(newx,check.names=FALSE),type='response')))
  } else if (method=='EN') {
    pf <- as.numeric(colnames(x) %in% xx$esg) # 0 controls, 1 ESG predictors
    if (!any(pf > 0)) stop('M0 uses plain LogR as EN reference (no penalised ESG predictors)')
    fit <- glmnet::glmnet(x,y,family='binomial',alpha=pars$alpha,
                         lambda=pars$lambda,penalty.factor=pf,
                         standardize=FALSE,intercept=TRUE)
    p <- as.numeric(predict(fit,newx=newx,type='response',s=pars$lambda))
  } else if (method=='RF') {
    df <- data.frame(Y=factor(y,levels=c(0,1)),as.data.frame(x,check.names=FALSE),check.names=FALSE)
    fit <- ranger::ranger(Y ~ .,data=df,probability=TRUE,
                          num.trees=pars$num.trees,mtry=min(pars$mtry,ncol(x)),
                          min.node.size=pars$min.node.size,
                          na.action='na.learn',seed=SEED,num.threads=2)
    # Confirm ranger did not silently omit training observations with missing values.
    stopifnot(fit$num.samples == nrow(df))
    p <- as.numeric(predict(fit,data=as.data.frame(newx,check.names=FALSE))$predictions[,'1'])
  } else if (method=='GB') {
    fit <- xgboost::xgb.train(params=list(objective='binary:logistic',eval_metric='aucpr',
                              max_depth=pars$max_depth,eta=pars$eta,
                              min_child_weight=pars$min_child_weight,
                              subsample=1,colsample_bytree=1,nthread=2,seed=SEED),
                              data=xgboost::xgb.DMatrix(x,label=y,missing=NA),
                              nrounds=pars$nrounds,verbose=0)
    p <- as.numeric(predict(fit,xgboost::xgb.DMatrix(newx,missing=NA)))
  } else stop('Unknown method ', method)
  if (any(!is.finite(p))) stop('Nonfinite predictions: ',spec,' ',method)
  p <- pmax(pmin(p,1),0)
  list(fit=fit,p=p,feature_names=colnames(x),train_x=x,test_x=newx,esg=xx$esg)
}
# Derive the glmnet regularisation path for each alpha and model specification.
# The same candidate lambdas are used for all expanding-window folds.
# Generate paths from FY2015--2017 only, to avoid using validation-year outcomes.
# glmnet internally rescales penalty.factor; do not use arbitrary fixed lambdas.
lambda_path <- function(train,spec,alpha,sector_levels,controls=ctl) {
  xx <- make_x(train,spec,sector_levels,controls)
  x <- scale_esg(xx$x,xx$x,xx$esg)$tr
  if (anyNA(x)) stop('Elastic Net lambda path needs imputed indicators')
  pf <- as.numeric(colnames(x) %in% xx$esg)
  if (!any(pf > 0)) stop('M0 uses plain LogR')
  path <- glmnet::glmnet(x,as.integer(train$Y),family='binomial',
                         alpha=alpha,penalty.factor=pf,
                         standardize=FALSE,nlambda=50)
  unique(as.numeric(path$lambda))
}
# Fixed grids for tree models; EN's grid is generated from the training data.
grid_for <- function(method,n_features) {
  if (method=='RF') return(expand.grid(num.trees=500L,
                                       mtry=unique(pmax(1,round(c(sqrt(n_features)/2, sqrt(n_features), n_features/3)))),
                                       min.node.size=c(15L,30L,60L)))
  if (method=='GB') return(expand.grid(nrounds=c(100L,300L,600L),
                                       max_depth=c(1L,2L,3L,4L,6L), eta=c(.01,.03,.1),
                                       min_child_weight=c(1,5)))
  data.frame(dummy=1)
}
# Fit once per tuning candidate in every expanding window.
tune <- function(train,spec,method,sector_levels,controls=ctl) {
  n_features <- ncol(make_x(train,spec,sector_levels,controls)$x)
  if (method=='EN') {
    alpha_grid <- c(.1,.25,.5,.75,1)
    initial_train <- train[train$year >= 2015 & train$year <= 2017,,drop=FALSE]
    grid <- do.call(rbind,lapply(alpha_grid,function(a)
      data.frame(alpha=a,lambda=lambda_path(initial_train,spec,a,sector_levels,controls))))
    rownames(grid) <- NULL
  } else grid <- grid_for(method,n_features)
  scores <- rep(NA_real_,nrow(grid)); oofs <- vector('list',nrow(grid))
  for (j in seq_len(nrow(grid))) {
    pset <- as.list(grid[j,,drop=FALSE]); fold_predictions <- list(); fold_scores <- numeric(length(folds))
    for (k in seq_along(folds)) {
      f <- folds[[k]]; tr <- train[train$year>=f[1] & train$year<=f[2],,drop=FALSE]
      va <- train[train$year==f[3],,drop=FALSE]
      res <- fit_predict(tr,va,spec,method,pset,sector_levels,controls)
      fold_scores[k] <- pr_auc(va$Y,res$p)
      fold_predictions[[k]] <- data.frame(RIC=va$RIC,year=va$year,Y=va$Y,p=res$p)
    }
    scores[j] <- mean(fold_scores)
    oofs[[j]] <- do.call(rbind,fold_predictions)
  }
  best <- which.max(scores)
  list(pars=as.list(grid[best,,drop=FALSE]),cv_pr_auc=scores[best],
       validation=oofs[[best]],grid=cbind(grid,mean_pr_auc=scores))
}
fit_one <- function(d,spec,method,tag='main',controls=ctl) {
  tr <- d[d$split=='train',,drop=FALSE]; te <- d[d$split=='test',,drop=FALSE]
  levels_sector <- sort(unique(as.character(d$sector25)))
  tuned <- tune(tr,spec,method,levels_sector,controls)
  final <- fit_predict(tr,te,spec,method,tuned$pars,levels_sector,controls)
  cut <- choose_cut(tuned$validation$Y,tuned$validation$p)
  pred <- data.frame(RIC=te$RIC,year=te$year,Y=te$Y,
                     p=final$p,model=spec,estimator=method,scenario=tag)
  metr <- data.frame(scenario=tag,model=spec,estimator=method,
                     AUC=roc_auc(pred$Y,pred$p),PR_AUC=pr_auc(pred$Y,pred$p),
                     F1=f1(pred$Y,pred$p,cut),F1_cutoff=cut,
                     CV_PR_AUC=tuned$cv_pr_auc)
  list(metrics=metr,pred=pred,final=final,cv=tuned)
}
# Matched firm-year paired cluster bootstrap for differences (B-A).
paired_boot <- function(a,b,B=1000,seed=SEED) {
  stopifnot(!anyDuplicated(a[c('RIC','year')]),!anyDuplicated(b[c('RIC','year')]))
  z <- merge(a[,c('RIC','year','Y','p')],b[,c('RIC','year','Y','p')],
             by=c('RIC','year','Y'),suffixes=c('_a','_b'))
  if (nrow(z)!=nrow(a) || nrow(z)!=nrow(b)) stop('Comparison does not have identical test rows')
  firms <- unique(z$RIC)
  ix <- split(seq_len(nrow(z)),factor(z$RIC,levels=firms))
  set.seed(seed)
  vals <- replicate(B, {
    draws <- sample(seq_along(firms),length(firms),replace=TRUE)
    zz <- z[unlist(ix[draws],use.names=FALSE),,drop=FALSE]
    c(AUC=roc_auc(zz$Y,zz$p_b)-roc_auc(zz$Y,zz$p_a),
      PR_AUC=pr_auc(zz$Y,zz$p_b)-pr_auc(zz$Y,zz$p_a))
  })
  data.frame(metric=c('AUC','PR_AUC'),
    delta=c(roc_auc(z$Y,z$p_b)-roc_auc(z$Y,z$p_a),
            pr_auc(z$Y,z$p_b)-pr_auc(z$Y,z$p_a)),
    lower=apply(vals,1,quantile,.025,na.rm=TRUE),
    upper=apply(vals,1,quantile,.975,na.rm=TRUE))
}
# Convert predictors to a clear indicator/sign/polarity table.
interpret_en <- function(res,spec) {
  fit <- res$final$fit; features <- res$final$feature_names
  val <- as.numeric(as.matrix(coef(fit,s=res$cv$pars$lambda))[-1,1])
  stopifnot(length(val)==length(features))
  z <- data.frame(variable=features,coefficient=val)
  z <- merge(z,ind[,c('var','type','direction')],by.x='variable',by.y='var',all.x=TRUE,sort=FALSE)
  z$non_disclosure <- grepl('^M_',z$variable)
  z$direction[z$non_disclosure] <- 'missingness'
  z$contrary_to_expected <- ifelse(z$direction=='better',z$coefficient>0,
                                  ifelse(z$direction=='worse',z$coefficient<0,NA))
  z$contrary_to_expected[is.na(z$coefficient) | z$coefficient == 0] <- NA
  z$spec <- spec;z
}
# Native XGBoost TreeSHAP on test observations; log-odds units.
interpret_gb <- function(res,te,spec,tag='main') {
  shap <- predict(res$final$fit,xgboost::xgb.DMatrix(res$final$test_x,missing=NA),
                  predcontrib=TRUE)
  shap <- as.matrix(shap)
  # XGBoost returns one extra intercept column, normally named BIAS.
  # Drop the final column positionally, then verify the feature count.
  stopifnot(ncol(shap) == ncol(res$final$test_x) + 1L)
  shap <- shap[,seq_len(ncol(shap)-1L),drop=FALSE]
  stopifnot(ncol(shap) == ncol(res$final$test_x))
  colnames(shap) <- colnames(res$final$test_x)
  row_shap <- data.frame(RIC=te$RIC,year=te$year,Y=te$Y,shap,check.names=FALSE)
  write.csv(row_shap,file.path(OUT,paste0(tag,'_',spec,'_GB_SHAP_by_row.csv')),row.names=FALSE)
  x_te <- res$final$test_x
  stopifnot(nrow(shap)==nrow(x_te),ncol(shap)==ncol(x_te))
  shap_direction <- vapply(seq_len(ncol(shap)),function(j) {
    if (sd(x_te[,j],na.rm=TRUE)==0 || sd(shap[,j],na.rm=TRUE)==0) return(NA_real_)
    suppressWarnings(cor(x_te[,j],shap[,j],use='complete.obs'))
  },numeric(1))
  ans <- data.frame(variable=colnames(shap),mean_abs_SHAP=colMeans(abs(shap)),
                    mean_SHAP=colMeans(shap),shap_value_correlation=shap_direction,
                    spec=spec,estimator='GB',scale='log_odds')
  ans$direction <- ind$direction[match(ans$variable,ind$var)]
  ans$contrary_to_expected <- ifelse(ans$direction=='better',ans$shap_value_correlation>0,
                                   ifelse(ans$direction=='worse',ans$shap_value_correlation<0,NA))
  ans$raises_predicted_risk <- ifelse(grepl('^M_',ans$variable),
                                     ans$shap_value_correlation>0,NA)
  ans
}
# Approximate RF SHAP: optional; VERY expensive on all test observations.
# RF SHAP uses probability units and is Monte Carlo approximate (nsim=30).
# Do NOT compare RF SHAP magnitudes against GB's log-odds-scale values.
# Uses a training background sample and test observations, predicting probabilities.
interpret_rf <- function(res,spec,background_n=200L,nsim=30L) {
  if (!requireNamespace('fastshap',quietly=TRUE)) stop('install.packages("fastshap")')
  set.seed(SEED)
  x <- as.data.frame(res$final$train_x,check.names=FALSE)
  nx <- as.data.frame(res$final$test_x,check.names=FALSE)
  bg <- x[sample.int(nrow(x),min(nrow(x),background_n)),,drop=FALSE]
  # Handle native NA predictors in the background and explained rows.
  pw <- function(object,newdata) as.numeric(predict(object,data=newdata)$predictions[,'1'])
  sh <- fastshap::explain(res$final$fit,X=bg,newdata=nx,
                          pred_wrapper=pw,nsim=nsim,adjust=FALSE)
  sh <- as.matrix(sh)
  data.frame(variable=colnames(sh),mean_abs_SHAP=colMeans(abs(sh)),
             mean_SHAP=colMeans(sh),spec=spec,estimator='RF',scale='probability')
}


# =============================================================================
# 00. LOAD DATA + REUSABLE MODEL
# =============================================================================
# Set your working directory to the folder containing data_main.rds and ind_class.rds.
d <- readRDS('data_main.rds')
stopifnot(all(c('RIC','year','split','Y','sector25',ctl,score,pillars,vars,miss_vars) %in% names(d)))
stopifnot(!anyNA(d[,c('RIC','year','Y','sector25',ctl)]))
stopifnot(all(d$year[d$split=='train'] <= 2020), all(d$year[d$split=='test'] >= 2021))
models <- new.env(parent=emptyenv())
sector_levels <- sort(unique(as.character(d$sector25)))
get_model <- function(spec, method) {
  id <- paste(spec,method,sep='_')
  if (!exists(id,envir=models,inherits=FALSE)) {
    message('Training ',id,' ...')
    ans <- fit_one(d,spec,method)
    assign(id,ans,envir=models)
    saveRDS(ans,file.path(OUT,paste0('main_',id,'.rds')))
    write.csv(ans$cv$grid,file.path(OUT,paste0('main_',id,'_tuning.csv')),row.names=FALSE)
    write.csv(ans$pred,file.path(OUT,paste0('main_',id,'_predictions.csv')),row.names=FALSE)
  }
  get(id,envir=models,inherits=FALSE)
}
comparison <- function(h,a_spec,a_method,b_spec,b_method) {
  aa <- get_model(a_spec,a_method)
  bb <- get_model(b_spec,b_method)
  out <- paired_boot(aa$pred,bb$pred,B=B)
  out$hypothesis <- h
  out$reference <- paste(a_spec,a_method)
  out$candidate <- paste(b_spec,b_method)
  out <- label_evidence(out)
  out
}
save_test <- function(x, name) {
  print(x)
  write.csv(x,file.path(OUT,paste0(name,'.csv')),row.names=FALSE)
}
# For each model comparison, define a joint hypothesis conclusion from BOTH
# primary metrics, before looking at test results. Each metric's individual
# CI remains visible as its own row.
label_evidence <- function(x) {
  stopifnot(all(c('hypothesis','reference','candidate','metric','lower','upper') %in% names(x)))
  key <- interaction(x$hypothesis,x$reference,x$candidate,drop=TRUE)
  x$metric_support <- ifelse(x$lower > 0,'positive',
                             ifelse(x$upper < 0,'negative','inconclusive'))
  x$comparison_conclusion <- NA_character_
  for (g in levels(key)) {
    ix <- which(key == g)
    if (!setequal(x$metric[ix],c('AUC','PR_AUC'))) next
    pos <- sum(x$lower[ix] > 0 & x$upper[ix] > 0)
    neg <- sum(x$upper[ix] < 0)
    x$comparison_conclusion[ix] <- if (pos == 2) 'supported' else
      if (pos == 1 && neg == 1) 'mixed' else
        if (pos == 1) 'partly_supported' else 'not_supported'
  }
  x
}
# These are rules for each MODEL COMPARISON, not an automatic overall
# decision for H1/H3/H5a (which contain several comparisons).

# =============================================================================
# H1. DO ESG COMPOSITE / PILLAR SCORES IMPROVE PREDICTION?
# M1a vs M0 and M1b vs M0, using plain LogR and Elastic Net.
# =============================================================================
h1 <- do.call(rbind,list(
  comparison('H1','M0','LogR','M1a','LogR'),
  comparison('H1','M0','LogR','M1b','LogR'),
  comparison('H1','M0','LogR','M1a','EN'),
  comparison('H1','M0','LogR','M1b','EN')
))
save_test(h1,'H1_composite_pillars_vs_controls')

# Export plain-LogR coefficient tables for M1a and M1b (Chapter 4.4).
# ESG variables have been standardised using full training-period statistics.
# Coefficients are in log-odds units per one training standard deviation.
# Controls retain their original units; sector terms are relative to the reference sector.
for (s in c('M1a','M1b')) {
  fit <- get_model(s,'LogR')$final$fit
  tab <- as.data.frame(summary(fit)$coefficients, check.names=FALSE)
  tab$variable <- rownames(tab)
  tab <- tab[,c('variable',setdiff(names(tab),'variable')),drop=FALSE]
  rownames(tab) <- NULL
  write.csv(tab,file.path(OUT,paste0('H1_',s,'_LogR_coefficients.csv')),row.names=FALSE)
}

tr <- d[d$split == 'train', ]
for (s in c('M1a','M1b')) {
  fit <- get_model(s,'LogR')$final$fit
  ct <- lmtest::coeftest(fit, vcov = sandwich::vcovCL(fit, cluster = tr$RIC))
  print(ct)
  m <- matrix(as.numeric(ct), nrow = nrow(ct), dimnames = dimnames(ct))
  write.csv(data.frame(variable = rownames(m), m, check.names = FALSE),
            file.path(OUT, paste0('H1_', s, '_LogR_clustered_SE.csv')), row.names = FALSE)
}
sd(tr$esg_score)   # training SD of the ESG score (coefficients above are per one SD)

# =============================================================================
# H2. DO INDIVIDUAL ESG INDICATORS BEAT PILLAR SCORES?
# Elastic Net only: M2 (indicators) versus M1b (pillars).
# =============================================================================
h2 <- comparison('H2','M1b','EN','M2','EN')
save_test(h2,'H2_indicators_vs_pillars')

# =============================================================================
# H2 INCREMENTAL. DO SCORES ADD INFORMATION BEYOND INDICATORS, AND VICE VERSA?
# Run AFTER H3-H5b have finished. Elastic Net only; new specifications M2a/M2b.
# Each comparison uses the same firm-years and the existing paired bootstrap.
# =============================================================================
h2_inc <- do.call(rbind,list(
  comparison('H2inc','M1a','EN','M2a','EN'), # Indicators beyond composite
  comparison('H2inc','M1b','EN','M2b','EN'), # Indicators beyond pillars
  comparison('H2inc','M2','EN','M2a','EN'),  # Composite beyond indicators
  comparison('H2inc','M2','EN','M2b','EN')   # Pillars beyond indicators
))
save_test(h2_inc,'H2_incremental')

# =============================================================================
# H3. DO TREE MODELS BEAT ELASTIC NET WITH THE SAME INDICATORS?
# RF and GB versus EN, separately for M2 and M4.
# =============================================================================
h3 <- do.call(rbind,list(
  comparison('H3','M2','EN','M2','RF'),
  comparison('H3','M2','EN','M2','GB'),
  comparison('H3','M4','EN','M4','RF'),
  comparison('H3','M4','EN','M4','GB')
))
save_test(h3,'H3_tree_models_vs_elastic_net')

# =============================================================================
# H4. WHICH REPORTED INDICATORS MATTER, AND IN WHICH DIRECTION?
# M2 EN coefficients + M2 GB SHAP; read signs against indicator polarity.
# =============================================================================
h4_en <- interpret_en(get_model('M2','EN'),'M2')
h4_en <- subset(h4_en, variable %in% vars)
h4_en$nonzero <- h4_en$coefficient != 0
save_test(h4_en,'H4_M2_elastic_net_directions')
# SHAP values on TEST data (XGBoost contribution scale: log-odds).
h4_gb <- interpret_gb(get_model('M2','GB'),d[d$split=='test',],'M2')
save_test(h4_gb,'H4_M2_gradient_boosting_SHAP')
# RF SHAP (Monte Carlo approximation on probability scale). Set TRUE to run.
# If FALSE, thesis methods should state SHAP is reported for GB only.
RUN_RF_SHAP <- FALSE
if (RUN_RF_SHAP) {
  h4_rf <- interpret_rf(get_model('M2','RF'),'M2')
  save_test(h4_rf,'H4_M2_random_forest_SHAP')
}

# =============================================================================
# H5a. DOES NON-DISCLOSURE ADD PREDICTIVE VALUE?
# M3 (missingness only) vs M0; M4 (indicators+missingness) vs M2.
# Compare like-for-like estimator where it exists (M0 plain LogR is
# unpenalised-control reference for the EN variant of M3).
# =============================================================================
h5a <- do.call(rbind,list(
  comparison('H5a','M0','LogR','M3','EN'),
  comparison('H5a','M2','EN','M4','EN'),
  comparison('H5a','M0','RF','M3','RF'),
  comparison('H5a','M2','RF','M4','RF'),
  comparison('H5a','M0','GB','M3','GB'),
  comparison('H5a','M2','GB','M4','GB')
))
save_test(h5a,'H5a_missingness_incremental_prediction')

# =============================================================================
# H5b. WHICH MISSING INDICATORS MATTER? OUTCOME VERSUS INPUT?
# M3 and M4 EN coefficients + SHAP. The outcome/input comparison is
# DESCRIPTIVE because only 3 input indicators have missingness flags.
# =============================================================================
h5b_en <- do.call(rbind,lapply(c('M3','M4'),function(s) {
  x <- interpret_en(get_model(s,'EN'),s)
  x <- subset(x,non_disclosure)
  x$type <- ind$type[match(sub('^M_','',x$variable),ind$var)]
  x$raises_predicted_risk <- x$coefficient > 0
  x$raises_predicted_risk[is.na(x$coefficient) | x$coefficient == 0] <- NA
  x
}))
save_test(h5b_en,'H5b_elastic_net_missingness_signs')
for (s in c('M3','M4')) {
  sh <- interpret_gb(get_model(s,'GB'),d[d$split=='test',],s)
  sh <- subset(sh,grepl('^M_',variable))
  sh$type <- ind$type[match(sub('^M_','',sh$variable),ind$var)]
  save_test(sh,paste0('H5b_',s,'_GB_missingness_SHAP'))
}
# RF SHAP is approximate (Monte Carlo) and uses probability units,
# while XGBoost TreeSHAP uses log-odds. Compare rankings, not magnitudes.
if (RUN_RF_SHAP) for (s in c('M3','M4')) {
  save_test(interpret_rf(get_model(s,'RF'),s),paste0('H5b_',s,'_RF_SHAP'))
}


# SAVE SUMMARY (run after whichever sections have been completed)
# =============================================================================
all_metrics <- lapply(ls(models),function(n) get(n,envir=models)$metrics)
save_test(do.call(rbind,all_metrics),'MAIN_all_model_metrics')
capture.output(sessionInfo(),file=file.path(OUT,'sessionInfo.txt'))

save_test(comparison('H3ctrl','M0','LogR','M0','GB'), 'H3_context_M0_GB_vs_LogR')
save_test(comparison('H3ctrl','M0','LogR','M0','RF'), 'H3_context_M0_RF_vs_LogR')

for (id in c('M0_RF','M2_RF','M3_RF','M4_RF','M0_GB','M2_GB','M3_GB','M4_GB')) {
  g <- read.csv(file.path(OUT, paste0('main_', id, '_tuning.csv')))
  best <- g[which.max(g$mean_pr_auc), ]
  pars <- setdiff(names(g), c('mean_pr_auc','num.trees'))
  edge <- sapply(pars, function(p) if (best[[p]] == min(g[[p]])) ' (min)' else if (best[[p]] == max(g[[p]])) ' (max)' else '')
  cat(id, ': ', paste0(pars, '=', unlist(best[pars]), edge, collapse = ', '), '\n', sep = '')
}
prof <- function(id, p) {
  g <- read.csv(file.path(OUT, paste0('main_', id, '_tuning.csv')))
  cat(id, '\n'); print(aggregate(as.formula(paste('mean_pr_auc ~', p)), g, max))
}
prof('M2_GB','max_depth'); prof('M4_GB','max_depth'); prof('M2_RF','mtry'); prof('M0_RF','min.node.size')

# =============================================================================
# ROBUSTNESS R1--R5 (RUN ONLY AFTER MAIN HYPOTHESES)
# =============================================================================
# R1: drop Trend; R2: drop Size; R3: native missing values for RF/GB ONLY;
# R4a: all indicators imputed at median; R4b: alternative outcome imputation.
# R5: continuous indicators as within-sector-by-year percentile ranks.
# These runs are kept in SEPARATE caches and folders to avoid any accidental
# reuse of main-sample fits. Each scenario uses its own dataset and controls.
# Refit the model comparisons for H1/H2/H3/H5a as applicable. H4/H5b
# interpretation can be produced by calling interpret_en/interpret_gb on
# corresponding scenario models if requested.
# Create R5 only when needed, after the main analysis has finished.
# Sector (not finer industry group) is used to avoid very small peer groups.
# Binary indicators and missingness flags are left unchanged.
create_r5 <- function(path='data_r5.rds') {
  if (file.exists(path)) return(invisible(path))
  if (!requireNamespace('dplyr',quietly=TRUE))
    stop('Install dplyr to generate R5: install.packages("dplyr")')
  stopifnot(all(c('var','unit') %in% names(ind)))
  if (anyNA(ind$unit)) stop('Missing indicator unit classifications in ind_class.rds')
  cont <- ind$var[ind$unit != 'binary']
  stopifnot(length(cont) == 13L, all(cont %in% names(d)))
  d5 <- readRDS('data_main.rds')
  d5 <- dplyr::ungroup(dplyr::mutate(
    dplyr::group_by(d5, sector25, year),
    dplyr::across(dplyr::all_of(cont), dplyr::percent_rank)
  ))
  d5 <- as.data.frame(d5)
  if (anyNA(d5[, cont])) stop('R5: NA ranks - check sector-years with only one firm')
  if (nrow(d5) != nrow(d)) stop('R5 changed the number of rows')
  saveRDS(d5,path)
  invisible(path)
}
robustness_scenarios <- list(
  R1 = list(path='data_main.rds', controls=setdiff(ctl,'Trend')),
  R2 = list(path='data_main.rds', controls=setdiff(ctl,'Size')),
  R3 = list(path='data_nan.rds', controls=ctl),
  R4a = list(path='data_r4a.rds', controls=ctl),
  R4b = list(path='data_r4b.rds', controls=ctl),
  R5 = list(path='data_r5.rds', controls=ctl)
)
run_robustness <- function(scenario) {
  sc <- robustness_scenarios[[scenario]]
  if (is.null(sc)) stop('Unknown scenario: ',scenario)
  if (scenario=='R5') create_r5(sc$path)
  dd <- readRDS(sc$path)
  needed <- unique(c('RIC','year','split','Y','sector25',sc$controls,
                     score,pillars,vars,miss_vars))
  if (!all(needed %in% names(dd))) stop('Missing scenario columns: ',
                                      paste(setdiff(needed,names(dd)),collapse=', '))
  if (anyNA(dd[,c('RIC','year','split','Y','sector25',sc$controls)]))
    stop('Missing identifiers, outcomes or controls in ', scenario)
  memo <- new.env(parent=emptyenv())
  get_scenario <- function(spec,method) {
    if (scenario=='R3' && method %in% c('EN','LogR'))
      stop('R3 cannot run EN/LogR with native NA values')
    key <- paste(spec,method,sep='_')
    if (!exists(key,envir=memo,inherits=FALSE)) {
      message(scenario, ': fitting ',key)
      ans <- fit_one(dd,spec,method,tag=scenario,controls=sc$controls)
      assign(key,ans,envir=memo)
      saveRDS(ans,file.path(OUT,paste0(scenario,'_',key,'.rds')))
      write.csv(ans$pred,file.path(OUT,paste0(scenario,'_',key,'_predictions.csv')),row.names=FALSE)
      write.csv(ans$cv$grid,file.path(OUT,paste0(scenario,'_',key,'_tuning.csv')),row.names=FALSE)
    }
    get(key,envir=memo,inherits=FALSE)
  }
  cmp <- function(h,ra,rb) {
    a <- get_scenario(ra[1],ra[2]); b <- get_scenario(rb[1],rb[2])
    z <- paired_boot(a$pred,b$pred,B=B)
    z$hypothesis <- h
    z$reference <- paste(ra,collapse=' ')
    z$candidate <- paste(rb,collapse=' ')
    z$scenario <- scenario
    label_evidence(z)
  }
  if (scenario=='R3') {
    # Probe ranger native missingness support with probability forests.
    if (!requireNamespace('ranger',quietly=TRUE)) stop('ranger required')
    probe <- data.frame(Y=factor(c(0,1,0,1,0,1),levels=c(0,1)),
                        x=c(1,NA,2,3,NA,4))
    probe_fit <- tryCatch(ranger::ranger(Y~.,data=probe,probability=TRUE,num.trees=10L,
                            na.action='na.learn',seed=SEED),
             error=function(e) stop('Installed ranger does not support na.learn with probability=TRUE: ',conditionMessage(e)))
    stopifnot(probe_fit$num.samples == nrow(probe))
    comparisons <- list(
      cmp('R3_RF_vs_GB_M2',c('M2','RF'),c('M2','GB')),
      cmp('R3_RF_vs_GB_M4',c('M4','RF'),c('M4','GB')),
      cmp('H5a',c('M0','RF'),c('M3','RF')),
      cmp('H5a',c('M2','RF'),c('M4','RF')),
      cmp('H5a',c('M0','GB'),c('M3','GB')),
      cmp('H5a',c('M2','GB'),c('M4','GB')))
    # R3 has NO comparison against Elastic Net because EN cannot use native NA.
  } else if (scenario %in% c('R4a','R4b','R5')) {
    # Only imputed-indicator models can differ from main. Skip M0/M1a/M1b/M3.
    comparisons <- list(
      cmp('H2',c('M1b','EN'),c('M2','EN')),
      cmp('H3',c('M2','EN'),c('M2','RF')),
      cmp('H3',c('M2','EN'),c('M2','GB')),
      cmp('H3',c('M4','EN'),c('M4','RF')),
      cmp('H3',c('M4','EN'),c('M4','GB')),
      cmp('H5a',c('M2','EN'),c('M4','EN')),
      cmp('H5a',c('M2','RF'),c('M4','RF')),
      cmp('H5a',c('M2','GB'),c('M4','GB')))
  } else {
    comparisons <- list(
      cmp('H1',c('M0','LogR'),c('M1a','LogR')),
      cmp('H1',c('M0','LogR'),c('M1b','LogR')),
      cmp('H1',c('M0','LogR'),c('M1a','EN')),
      cmp('H1',c('M0','LogR'),c('M1b','EN')),
      cmp('H2',c('M1b','EN'),c('M2','EN')),
      cmp('H3',c('M2','EN'),c('M2','RF')),
      cmp('H3',c('M2','EN'),c('M2','GB')),
      cmp('H3',c('M4','EN'),c('M4','RF')),
      cmp('H3',c('M4','EN'),c('M4','GB')),
      cmp('H5a',c('M0','LogR'),c('M3','EN')),
      cmp('H5a',c('M2','EN'),c('M4','EN')),
      cmp('H5a',c('M0','RF'),c('M3','RF')),
      cmp('H5a',c('M2','RF'),c('M4','RF')),
      cmp('H5a',c('M0','GB'),c('M3','GB')),
      cmp('H5a',c('M2','GB'),c('M4','GB')))
  }
  tab <- do.call(rbind,comparisons)
  write.csv(tab,file.path(OUT,paste0(scenario,'_comparisons.csv')),row.names=FALSE)
  print(tab)
  invisible(list(comparisons=tab,models=memo))
}
# Run individually because each scenario refits many tuned models:
# r1 <- run_robustness('R1')
# r2 <- run_robustness('R2')
# r3 <- run_robustness('R3')
# r4a <- run_robustness('R4a')
# r4b <- run_robustness('R4b')

# Delete any old R5 file so create_r5() rebuilds it from scratch
if (file.exists('data_r5.rds')) file.remove('data_r5.rds')
# r5 <- run_robustness('R5')

# =============================================================================
# ROBUSTNESS: SAVE MODEL METRICS PER SCENARIO
# Run each block after the matching run_robustness() call above.
# =============================================================================
## for r1
r1_metrics <- do.call(rbind, lapply(ls(r1$models), function(n) get(n, envir = r1$models)$metrics))
print(r1_metrics); write.csv(r1_metrics, file.path(OUT, 'R1_all_model_metrics.csv'), row.names = FALSE)
## for r2
r2_metrics <- do.call(rbind, lapply(ls(r2$models), function(n) get(n, envir = r2$models)$metrics))
print(r2_metrics); write.csv(r2_metrics, file.path(OUT, 'R2_all_model_metrics.csv'), row.names = FALSE)
x <- paired_boot(get_model('M0','LogR')$pred, r2$models$M2_GB$pred, B = B)
x$hypothesis <- 'R2ctx'; x$reference <- 'M0 LogR (main)'; x$candidate <- 'M2 GB (R2)'
save_test(x, 'R2_context_M0main_vs_M2GB_R2')
## for r3
r3_metrics <- do.call(rbind, lapply(ls(r3$models), function(n) get(n, envir = r3$models)$metrics))
print(r3_metrics); write.csv(r3_metrics, file.path(OUT, 'R3_all_model_metrics.csv'), row.names = FALSE)
identical(r3$models$M2_GB$pred$p, r3$models$M4_GB$pred$p)
c(ncol(r3$models$M2_GB$final$train_x), ncol(r3$models$M4_GB$final$train_x))
## for r4a
r4a_metrics <- do.call(rbind, lapply(ls(r4a$models), function(n) get(n, envir = r4a$models)$metrics))
print(r4a_metrics); write.csv(r4a_metrics, file.path(OUT, 'R4a_all_model_metrics.csv'), row.names = FALSE)

## for r4b
r4b_metrics <- do.call(rbind, lapply(ls(r4b$models), function(n) get(n, envir = r4b$models)$metrics))
print(r4b_metrics); write.csv(r4b_metrics, file.path(OUT, 'R4b_all_model_metrics.csv'), row.names = FALSE)

## for r5
r5_metrics <- do.call(rbind, lapply(ls(r5$models), function(n) get(n, envir = r5$models)$metrics))
print(r5_metrics); write.csv(r5_metrics, file.path(OUT, 'R5_all_model_metrics.csv'), row.names = FALSE)

# =============================================================================
# APPENDIX H: ROBUSTNESS TABLES -> LATEX ROWS
# =============================================================================
scen <- c('R1','R2','R3','R4a','R4b','R5')
titles <- c(R1='R1: without the time trend', R2='R2: without firm size',
            R3='R3: tree-based models with native missing values',
            R4a='R4a: median imputation for all indicators',
            R4b='R4b: 10th/90th percentile imputation for outcome metrics',
            R5='R5: continuous indicators as sector-year percentile ranks')
fmt <- function(v) ifelse(round(v,4) < 0, paste0('$-$', sprintf('%.4f', abs(v))), sprintf('%.4f', v))
ci  <- function(d,l,u) paste0(fmt(d), ' [', fmt(l), ', ', fmt(u), ']')
concl <- c(supported='supported', partly_supported='partly supported',
           mixed='mixed', not_supported='not supported')
ctx_file <- file.path(OUT, 'R2_context_M0main_vs_M2GB_R2.csv')

# Table H1: differences with CIs
rows <- character(0)
for (i in seq_along(scen)) {
  sc <- scen[i]
  x <- read.csv(file.path(OUT, paste0(sc, '_comparisons.csv')))
  x <- x[x$hypothesis %in% c('H1','H2','H3','H5a'), ]
  key <- unique(x[, c('hypothesis','reference','candidate')])
  if (i > 1) rows <- c(rows, '\\midrule')
  rows <- c(rows, sprintf('\\multicolumn{5}{l}{\\textit{Panel %s. %s}} \\\\', LETTERS[i], titles[sc]))
  for (k in seq_len(nrow(key))) {
    z <- x[x$hypothesis==key$hypothesis[k] & x$reference==key$reference[k] & x$candidate==key$candidate[k], ]
    a <- z[z$metric=='AUC', ]; p <- z[z$metric=='PR_AUC', ]
    rows <- c(rows, sprintf('%s & %s vs.\\ %s & %s & %s & %s \\\\',
                            key$hypothesis[k], key$candidate[k], key$reference[k],
                            ci(a$delta,a$lower,a$upper), ci(p$delta,p$lower,p$upper), concl[a$comparison_conclusion]))
  }
  if (sc == 'R2' && file.exists(ctx_file)) {
    cx <- read.csv(ctx_file)
    stopifnot(all(c('metric','delta','lower','upper') %in% names(cx)))
    a <- cx[cx$metric=='AUC', ]; p <- cx[cx$metric=='PR_AUC', ]
    rows <- c(rows, sprintf('Context & M2 GB (R2) vs.\\ M0 LogR (main) & %s & %s & -- \\\\',
                            ci(a$delta,a$lower,a$upper), ci(p$delta,p$lower,p$upper)))
  }
}
writeLines(rows, file.path(OUT, 'APPX_H_differences_rows.tex'))

# Table H2: AUC and PR-AUC per check
get_levels <- function(sc) {
  f <- file.path(OUT, paste0(sc, '_all_model_metrics.csv'))
  if (file.exists(f)) return(read.csv(f))
  rds <- list.files(OUT, pattern=paste0('^', sc, '_M[0-9]+[ab]?_(LogR|EN|RF|GB)\\.rds$'), full.names=TRUE)
  do.call(rbind, lapply(rds, function(r) readRDS(r)$metrics))
}
lev <- c(list(main=read.csv(file.path(OUT, 'MAIN_all_model_metrics.csv'))),
         setNames(lapply(scen, get_levels), scen))
ord <- unique(lev$main[, c('model','estimator')])
ord <- ord[order(ord$model, match(ord$estimator, c('LogR','EN','RF','GB'))), ]
cell <- function(x, m, e, col) {
  v <- x[x$model==m & x$estimator==e, col]
  if (length(v)==0) '--' else sprintf('%.3f', v[1])
}
lev_rows <- character(0)
for (j in 1:2) {
  col <- c('AUC','PR_AUC')[j]
  if (j == 2) lev_rows <- c(lev_rows, '\\midrule')
  lev_rows <- c(lev_rows, sprintf('\\multicolumn{8}{l}{\\textit{Panel %s. %s}} \\\\', LETTERS[j], c('AUC','PR-AUC')[j]))
  for (r in seq_len(nrow(ord))) {
    m <- ord$model[r]; e <- ord$estimator[r]
    cells <- vapply(names(lev), function(s) cell(lev[[s]], m, e, col), character(1))
    lev_rows <- c(lev_rows, paste0(m, ' ', e, ' & ', paste(cells, collapse=' & '), ' \\\\'))
  }
}
writeLines(lev_rows, file.path(OUT, 'APPX_H_levels_rows.tex'))
cat(rows, sep='\n'); cat(lev_rows, sep='\n')

# =============================================================================
# INDUSTRIES (EXPLANATORY): CO2 SHAP BY SECTOR, M2 GB
# =============================================================================

sh <- read.csv(file.path(OUT, 'main_M2_GB_SHAP_by_row.csv'), check.names = FALSE)
sh <- merge(sh, d[d$split == 'test', c('RIC', 'year', 'sector25')], by = c('RIC', 'year'))
co2 <- aggregate(co2_equivalent_emissions_total ~ sector25, data = sh,
                 FUN = function(x) c(mean_abs_SHAP = mean(abs(x)), n = length(x)))
co2 <- do.call(data.frame, co2)
print(co2[order(-co2[, 2]), ], digits = 3)
write.csv(co2, file.path(OUT, 'EXPL_CO2_SHAP_by_sector.csv'), row.names = FALSE)

# =============================================================================
# CHAPTER 5 MAIN TABLE: DESCRIPTIVE STATISTICS
# =============================================================================
dv <- c(Y = 'Controversy in $t+1$', PriorControversy = 'Prior controversy',
        Size = 'Size (log market cap)', esg_score = 'ESG score',
        e_pillar = 'Environmental pillar', s_pillar = 'Social pillar', g_pillar = 'Governance pillar')
stat <- function(x) c(N = length(x), Mean = mean(x), SD = sd(x), Min = min(x),
                      Median = median(x), Max = max(x))
desc <- t(sapply(names(dv), function(v) stat(d[[v]])))
desc <- cbind(desc,
              Mean_train = sapply(names(dv), function(v) mean(d[[v]][d$split == 'train'])),
              Mean_test  = sapply(names(dv), function(v) mean(d[[v]][d$split == 'test'])))
print(round(desc, 3))
write.csv(desc, file.path(OUT, 'DESC_main_variables.csv'))

fmt <- function(m, digits = 2) apply(m, 1, function(r) paste(formatC(r, format = 'f', digits = digits), collapse = ' & '))
cat(paste0(dv, ' & ', fmt(desc[, -1]), ' \\\\'), sep = '\n')

# =============================================================================
# APPENDIX D: LOGR COEFFICIENTS WITH CLUSTERED SEs
# =============================================================================

tr <- d[d$split == 'train', ]
get_ct <- function(s) {
  fit <- get_model(s, 'LogR')$final$fit
  ct  <- lmtest::coeftest(fit, vcov = sandwich::vcovCL(fit, cluster = tr$RIC))
  matrix(as.numeric(ct), nrow = nrow(ct), dimnames = dimnames(ct))
}
a <- get_ct('M1a'); b <- get_ct('M1b')

sec_raw  <- make.names(paste0('sector25', sector_levels[-1]), unique = TRUE)
sec_nice <- paste('Sector:', gsub('&', '\\\\&', sector_levels[-1]))
lab <- c(esg_score = 'ESG score', e_pillar = 'E pillar', s_pillar = 'S pillar', g_pillar = 'G pillar',
         Size = 'Size (log market cap)', PriorControversy = 'Prior controversy', Trend = 'Trend',
         '(Intercept)' = 'Intercept', setNames(sec_nice, sec_raw))
vars_order <- names(lab)
stopifnot(all(rownames(a) %in% vars_order), all(rownames(b) %in% vars_order))

cell <- function(m, v) if (v %in% rownames(m))
  sprintf('%.3f & %.3f & %s', m[v, 1], m[v, 2],
          ifelse(m[v, 4] < 0.001, '$<0.001$', sprintf('%.3f', m[v, 4]))) else ' & & '
rows_tex <- sprintf('%s & %s & %s \\\\', lab[vars_order],
                    sapply(vars_order, cell, m = a), sapply(vars_order, cell, m = b))
writeLines(rows_tex, file.path(OUT, 'APPX_LogR_coefficients_rows.tex'))
cat(rows_tex, sep = '\n')
cat('Reference sector:', sector_levels[1], '\n')
