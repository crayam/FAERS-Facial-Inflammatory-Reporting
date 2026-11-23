args <- commandArgs(trailingOnly=TRUE)
out <- args[1]
library(jsonlite)
a <- read.csv(file.path(out,'grouped_design_counts.csv'),check.names=FALSE)
X <- as.matrix(a[,1:9]); y <- a$positive; n <- a$N
colnames(X) <- c('Intercept','drug_tralokinumab','drug_upadacitinib','sex_Male','sex_Missing','HCP_NonHCP','HCP_Missing','geo_NonUS','quarter_index')
P <- ncol(X)
evaluate <- function(b) {
 eta <- drop(X %*% b); prob <- plogis(eta); w <- n*prob*(1-prob)
 info <- crossprod(X, X*w); ch <- try(chol(info),silent=TRUE)
 if(inherits(ch,'try-error'))return(list(ll=-Inf,score=rep(NA_real_,P)))
 logdet <- 2*sum(log(diag(ch))); inv <- chol2inv(ch)
 h <- w*rowSums((X %*% inv)*X)
 ll <- sum(y*eta - n*(pmax(eta,0)+log1p(exp(-abs(eta))))) + logdet/2
 score <- drop(crossprod(X,y-n*prob+h*(0.5-prob)))
 list(ll=ll,score=score,info=info,prob=prob,hat=h,logdet=logdet)
}
fit <- function(init,fixed=integer(),values=numeric()) {
 free <- setdiff(seq_len(P),fixed); base <- init; base[fixed] <- values
 expand <- function(par){b<-base;b[free]<-par;b}
 fn <- function(par)-evaluate(expand(par))$ll
 gr <- function(par)-evaluate(expand(par))$score[free]
 z <- optim(base[free],fn,gr,method='BFGS',control=list(maxit=5000,reltol=1e-14))
 par <- z$par; hist <- list(); last_change <- NA_real_; last_step <- NA_real_
 for(i in seq_len(40)) {
   ev <- evaluate(expand(par)); g <- ev$score[free]
   if(max(abs(g))<1e-8)break
   hess <- optimHess(par,fn,gr); step <- solve(hess,g)
   scale <- 1
   for(k in seq_len(50)){
     candidate<-par+scale*step; ce<-evaluate(expand(candidate))
     if(ce$ll>=ev$ll-1e-10 && (ce$ll>ev$ll || max(abs(ce$score[free]))<max(abs(g))))break
     scale<-scale/2
   }
   last_change<-abs(ce$ll-ev$ll);last_step<-max(abs(candidate-par));par<-candidate
   hist[[length(hist)+1]]<-c(ll=ce$ll,score=max(abs(ce$score[free])),step=last_step)
   if(max(abs(ce$score[free]))<1e-8)break
 }
 b<-expand(par); ev<-evaluate(b); score<-max(abs(ev$score[free]))
 if(!is.finite(ev$ll)||score>1e-6)stop(sprintf('Fit failed, score %.10g',score))
 list(beta=b,ll=ev$ll,score=score,optim_code=z$convergence,optim_counts=z$counts,refinement_steps=length(hist),last_ll_change=last_change,last_beta_step=last_step,converged=TRUE,free=free)
}
init<-rep(0,P);init[1]<-qlogis(sum(y)/sum(n))
full<-fit(init); ev<-evaluate(full$beta)
print(full)
se<-sqrt(diag(solve(ev$info)))
crit<-qchisq(.95,1); profiled<-list(); results<-list()
for(j in c(2,3)){
  cat('Profiling',colnames(X)[j],'\n')
  cache<-new.env(parent=emptyenv())
  profile<-function(theta){
    key<-sprintf('%.16g',theta)
    if(exists(key,cache,inherits=FALSE))return(get(key,cache))
    z<-fit(full$beta,j,theta);assign(key,z,cache);z
  }
  residual<-function(theta)2*(full$ll-profile(theta)$ll)-crit
  roots<-sapply(c(-1,1),function(direction){
    step<-2*se[j];edge<-full$beta[j]+direction*step
    for(k in seq_len(30)){if(residual(edge)>0)break;step<-step*1.6;edge<-full$beta[j]+direction*step}
    if(residual(edge)<=0)stop('Profile root could not be bracketed')
    uniroot(residual,interval=sort(c(full$beta[j],edge)),tol=1e-10)$root
  })
  null<-profile(0);lrt<-2*(full$ll-null$ll)
  ci_fits<-lapply(roots,profile)
  comparator<-if(j==2)'tralokinumab' else 'upadacitinib'
  results[[length(results)+1]]<-list(comparison=paste('dupilumab vs',comparator),comparator=comparator,log_OR=-full$beta[j],adjusted_OR=exp(-full$beta[j]),CI_lower=exp(-roots[2]),CI_upper=exp(-roots[1]),P_value=pchisq(lrt,1,lower.tail=FALSE),PLR_statistic=lrt,beta_index=j,beta_profile_lower=roots[1],beta_profile_upper=roots[2],null_fit=null,lower_profile_fit=ci_fits[[1]],upper_profile_fit=ci_fits[[2]],lower_LR_residual=residual(roots[1]),upper_LR_residual=residual(roots[2]),profile_evaluations=length(ls(cache)))
}
coef_table<-data.frame(term=colnames(X),beta=full$beta,information_SE=se,exp_beta=exp(full$beta),adjusted_score=ev$score)
write.csv(coef_table,file.path(out,'model_coefficients.csv'),row.names=FALSE)
write.csv(solve(ev$info),file.path(out,'information_covariance.csv'),row.names=TRUE)
H<-optimHess(full$beta,function(b)-evaluate(b)$ll,function(b)-evaluate(b)$score)
info<-list(software=R.version.string,implementation='R base optim BFGS + analytic Firth adjusted score + numerical-Newton refinement; independent NumPy scoring validation',package_logistf_used=FALSE,objective='sum(s_i*eta_i - n_i*log(1+exp(eta_i))) + 0.5*log(det(X^T W X)); no binomial combinatorial constant',grouping='Identical 9-column design patterns aggregated with exact positive/negative counts; equivalent to individual Bernoulli likelihood and information',CI_method='95% profile penalized likelihood; chi-square(1) cutoff',P_method='Two-sided penalized likelihood-ratio test of log OR = 0; chi-square(1)',penalty_under_constraint='Full 9-column information determinant retained for all constrained fits',N=sum(n),positive=sum(y),parameters=P,full=full,coefficients=coef_table,information_eigenvalues=eigen(ev$info,symmetric=TRUE,only.values=TRUE)$values,penalized_negative_Hessian_eigenvalues=eigen(H,symmetric=TRUE,only.values=TRUE)$values,probability_range=range(ev$prob),maximum_individual_hat=max(ev$hat/n),results=results)
write_json(info,file.path(out,'fit_results.json'),pretty=TRUE,auto_unbox=TRUE,digits=16,na='null')
saveRDS(info,file.path(out,'fit_results.rds'))
print(lapply(results,function(z)z[c('comparison','adjusted_OR','CI_lower','CI_upper','P_value','PLR_statistic','lower_LR_residual','upper_LR_residual')]))
