args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[1] else getwd(),winslash="/",mustWork=TRUE)
res <- file.path(root,"analysis/corrected_results")
validation_dir <- file.path(root,"validation")
dir.create(validation_dir,recursive=TRUE,showWarnings=FALSE)
checks <- list()
record <- function(label, pass, detail="") {checks[[length(checks)+1]] <<- data.frame(check=label, passed=pass, detail=detail)}
rd <- function(x) read.csv(file.path(res,x),check.names=FALSE)
near <- function(x,y) isTRUE(all.equal(as.numeric(x),as.numeric(y),tolerance=1e-7))
d <- rd("UC_reproducibility_gene_level.csv"); s <- rd("UC_within_center_reproducibility.csv")
record("UC rank correlation",near(cor(d$t_discovery,d$t_validation,method="spearman"),s$spearman_t))
same <- sign(d$t_discovery)==sign(d$t_validation); sel <- d$q_discovery<.05
record("UC directional counts",sum(sel)==s$discovery_significant_present && sum(sel & same)==s$direction_concordant_n && sum(sel & same & d$q_validation<.05)==s$both_significant_same_direction)
for (id in c("GSE16879","GSE36807")) {
 d <- rd(paste0("EoG_signature_",id,"_scores.csv")); x<-d$score[d$group=="UC"];y<-d$score[d$group=="CTRL"]
 n<-length(x);m<-length(y);g<-(mean(x)-mean(y))/sqrt(((n-1)*var(x)+(m-1)*var(y))/(n+m-2))*(1-3/(4*(n+m)-9))
 se<-sqrt((n+m)/(n*m)+g^2/(2*(n+m-2)))
 s<-rd(paste0("EoG_signature_",id,if(id=="GSE16879") "_sample_validation.csv" else "_conditional_sensitivity.csv"))
 record(paste(id,"score effect and approximate CI"),near(c(g,g-1.96*se,g+1.96*se,mean(x)-mean(y)),c(s$hedges_g,s$ci_low,s$ci_high,s$mean_difference)),paste("n",n,m,"g",g))
}
d<-rd("EoG_signature_GSE109142_exploratory_scores.csv");d<-d[complete.cases(d),]
fit<-lm(score~eos32+histology+age+sex,data=d);X<-model.matrix(fit);A<-solve(crossprod(X));h<-rowSums((X%*%A)*X);w<-(residuals(fit)/(1-h))^2
V<-A%*%crossprod(X,X*w)%*%A;se<-sqrt(diag(V));b<-coef(fit);p<-2*pt(-abs(b/se),df.residual(fit));s<-rd("EoG_signature_GSE109142_exploratory_model.csv")
record("HC3 model from 198 complete score records",nrow(d)==198 && near(b,s$estimate) && near(se,s$robust_se) && near(p,s$p_value) && near(b-1.96*se,s$ci_low) && near(b+1.96*se,s$ci_high),paste("df",df.residual(fit),"independent sandwich calculation"))
for (f in list.files(res,pattern="^ORA_.*csv$")) {
 d<-rd(f);p<-phyper(d$overlap-1,d$set_size,12236-d$set_size,209,lower.tail=FALSE)
 record(paste(f,"hypergeometric and complete-family BH"),near(p,d$p_value) && near(p.adjust(p,"BH"),d$q_value),paste(nrow(d),"tests;",sum(d$q_value<.05),"significant"))
}
for (f in c("Hallmark_primary_discovery_validation.csv","Hallmark_exploratory_EoG_UC.csv")) {
 d<-rd(f);a<-if(grepl("primary",f)) "EoC" else "EoG"
 ok<-TRUE
 for(k in c(a,"UC_discovery","UC_validation")) ok<-ok && near(p.adjust(d[[paste0("PValue_",k)]],"BH"),d[[paste0("FDR_",k)]])
 candidate<-d[[paste0("FDR_",a)]]<.05 & d$FDR_UC_discovery<.05 & d[[paste0("Direction_",a)]]==d$Direction_UC_discovery
 holm<-rep(NA_real_,nrow(d));holm[candidate]<-p.adjust(d$PValue_UC_validation[candidate],"holm")
 rep<-candidate & d[[paste0("Direction_",a)]]==d$Direction_UC_validation & holm<.05
 record(paste(f,"BH, candidate Holm, direction"),ok && near(holm,d$validation_holm) && identical(rep,d$replicated),paste(sum(rep),"replicated"))
}
mr<-file.path(root,"analysis/matched_null_2026-08-28/results")
d<-read.csv(file.path(mr,"null_endpoint_statistics.csv"));s<-read.csv(file.path(mr,"competitive_test_results.csv"))
for (i in seq_len(nrow(s))) {
 z<-d$statistic[d$null_model==s$null_model[i] & d$endpoint==s$endpoint[i]]
 p<-(1+sum(z>=s$observed[i]))/(length(z)+1)
 record(paste("Matched null",s$null_model[i],s$endpoint[i]),near(p,s$p_one_sided[i]) && length(z)==10000,paste("P",p))
}
out<-do.call(rbind,checks);write.csv(out,file.path(validation_dir,"derived_statistics_checks.csv"),row.names=FALSE)
print(out);if(any(!out$passed)) quit(status=2)
