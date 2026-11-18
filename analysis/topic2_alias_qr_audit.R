# Review-driven 26-record comparison of alias(), QR and expression-aware audit.
# Uses one common early-period NHANES sample and the current package sources.
root <- "E:/KimiData/kimi/Workspaces/论文写作"
pkg <- file.path(root, "课题2_v38_BMC_MIDM投稿包/04_代码与GitHub仓库/composite-index-audit")
Sys.setenv(CIA_CODE_DIR=file.path(pkg,"analysis"),
           CIA_DATA_ROOT=file.path(root,"work_v36/nhanes_data"),
           CIA_OUT=file.path(root,"work_revision_20261117_1831/topic2_results"),
           CIA_BOOTSTRAP="2")
dir.create(Sys.getenv("CIA_OUT"),recursive=TRUE,showWarnings=FALSE)
build <- readLines(file.path(pkg,"analysis/build_results.R"),warn=FALSE)
prefix <- paste(build[seq_len(which(grepl("^Ld <- build_period",build))[1])],collapse="\n")
eval(parse(text=prefix),envir=globalenv())
score_vector <- function(d,index,scale_override=NULL) {
  score <- v28_score_is_valid(d[[index]],index)
  scale <- if (is.null(scale_override)) registry$model_scale[match(index,registry$index)] else scale_override
  if (identical(scale,"log") && index!="HOMA_IR_log") score <- ifelse(is.finite(score)&score>0,log(score),NA_real_)
  score
}
ids <- c(names(reference_audit_classes()),"GLM7-b")
base_ids <- sub("-b$","",ids)
for (idx in unique(base_ids)) Pd[[paste0("C_",idx)]] <- score_vector(Pd,idx)
all_bases <- unique(unlist(lapply(unique(base_ids),function(idx) v28_component_model_terms(idx,registry))))
all_bases <- unique(c(all_bases,"metab"))
score_cols <- paste0("C_",unique(base_ids))
req <- unique(c("dx_diab",all_bases,score_cols))
num_req <- setdiff(req,"sex")
ok <- complete.cases(Pd[,req]) & apply(as.data.frame(lapply(Pd[,num_req,drop=FALSE],function(x) is.finite(x))),1,all)
d <- Pd[ok,,drop=FALSE]
stopifnot(nrow(d)>1000L,!anyDuplicated(d$SEQN))
rows <- list()
for (i in seq_along(ids)) {
  rec <- ids[i]; idx <- base_ids[i]
  mt <- v28_component_model_terms(idx,registry)
  if (rec=="GLM7-b") mt <- unique(c(mt,"metab"))
  work <- d
  work$C <- work[[paste0("C_",idx)]]
  f0 <- reformulate(mt,response="dx_diab")
  f1 <- update(f0,~ . + C)
  fit <- lm(f1,data=work)
  al <- alias(fit)
  dropped <- names(coef(fit))[is.na(coef(fit))]
  complete <- al$Complete
  alias_C <- if (!is.null(complete) && "C" %in% rownames(complete)) paste(sprintf("%.6g",complete["C",]),colnames(complete),sep="*",collapse=";") else ""
  xb <- model.matrix(f0,data=work)
  xa <- model.matrix(f1,data=work)
  rb <- qr(xb,tol=1e-10)$rank; ra <- qr(xa,tol=1e-10)$rank
  ranks <- vapply(c(1e-8,1e-10,1e-12),function(tol) qr(xa,tol=tol)$rank-qr(xb,tol=tol)$rank,integer(1))
  aa <- audit_one_record(rec,data=work,tol=1e-10)
  model_scale <- registry$model_scale[match(idx,registry$index)]
  audited <- if (model_scale=="log" && idx!="HOMA_IR_log") paste0("log(",idx,")") else idx
  z <- xa[,colnames(xa)!="(Intercept)",drop=FALSE]
  sds <- apply(z,2,sd); usable <- is.finite(sds)&sds>0
  kval <- tryCatch(kappa(scale(z[,usable,drop=FALSE]),exact=TRUE),error=function(e) NA_real_)
  rows[[rec]] <- data.frame(record=rec,index=idx,n=nrow(work),events=sum(work$dx_diab),
    registered_base=paste(mt,collapse=" + "),native_definition=registry$formula[match(idx,registry$index)],
    audited_composite_column=audited,alias_C=ifelse("C" %in% dropped,"aliased","estimable"),
    alias_dropped_columns=paste(dropped,collapse=";"),alias_C_dependence=alias_C,
    rank_base=rb,rank_augmented=ra,rank_gain=ra-rb,
    rank_gain_tol_1e8=ranks[1],rank_gain_tol_1e10=ranks[2],rank_gain_tol_1e12=ranks[3],
    standardized_condition_number=kval,expression_class=aa$actual_class,
    derivation_path=aa$path,blocking_reason=ifelse(is.na(aa$blocking_reason),"",aa$blocking_reason),
    stringsAsFactors=FALSE)
}
out <- do.call(rbind,rows); rownames(out)<-NULL
stopifnot(nrow(out)==26L)
write.csv(out,file.path(Sys.getenv("CIA_OUT"),"v38_alias_QR_expression_26_records.csv"),row.names=FALSE,na="")
summary <- as.data.frame.matrix(table(out$expression_class,out$alias_C))
summary$expression_class <- rownames(summary); rownames(summary)<-NULL
write.csv(summary,file.path(Sys.getenv("CIA_OUT"),"v38_alias_QR_expression_summary.csv"),row.names=FALSE)
meta <- data.frame(item=c("common_sample","period","QR tolerance","alias outcome","software source","scope"),
  value=c(nrow(d),"2017-March 2020","1e-10; sensitivity 1e-8 and 1e-12","self-reported diabetes; alias geometry does not depend on outcome",pkg,"Numerical comparison on realized columns; no independent source-formula certification"))
write.csv(meta,file.path(Sys.getenv("CIA_OUT"),"v38_alias_QR_metadata.csv"),row.names=FALSE)
print(out[,c("record","n","alias_C","alias_dropped_columns","rank_base","rank_augmented","rank_gain","expression_class")])
print(summary)
cat("TOPIC2_ALIAS_QR_26_COMPLETE\n")
