# Independent reference fit for the audit: factorial coding, uncentered M-values.
# This does not call the production four-cell R script or Python moderation code.
.libPaths(c(file.path(getwd(), ".tools/r-library"), .libPaths()))
suppressPackageStartupMessages(library(limma))
stopifnot(as.character(packageVersion("limma")) == "3.54.0")
args <- commandArgs(trailingOnly=TRUE)
assignments <- read.csv(args[1], check.names=FALSE, colClasses="character")
out <- args[2]
samples <- read.csv("data/processed/samples.csv", colClasses=c(participant="character", diet_code="character"))
cpg <- readLines("data/processed/cpg_ids.txt")
beta <- matrix(readBin("data/processed/beta.bin", "double", n=length(cpg)*nrow(samples),
                       size=8, endian="little"), nrow=length(cpg), byrow=TRUE)
M <- log2(beta/(1-beta))
rm(beta)
visit <- factor(samples$time, levels=c(0,1))
observed <- lapply(c("before", "after"), function(time)
    read.csv(gzfile(file.path("results/diet_code", paste0(time,"_limma.csv.gz")))))
ans <- list()
for (i in seq_len(nrow(assignments))) {
    diet <- factor(as.character(unlist(assignments[i, match(samples$participant, names(assignments))])),
                   levels=c("0","1"))
    stopifnot(!anyNA(diet))
    design <- model.matrix(~diet*visit)
    stopifnot(identical(colnames(design), c("(Intercept)","diet1","visit1","diet1:visit1")))
    contrasts <- cbind(before=c(0,-1,0,0), after=c(0,-1,0,-1))
    fit <- eBayes(contrasts.fit(lmFit(M,design),contrasts), trend=FALSE, robust=FALSE)
    row <- data.frame(assignment=assignments$assignment[i],
        before=sum(fit$p.value[,1]<0.05 & abs(fit$coefficients[,1])>0.5),
        after=sum(fit$p.value[,2]<0.05 & abs(fit$coefficients[,2])>0.5),
        df_prior=fit$df.prior, s2_prior=fit$s2.prior)
    if (assignments$assignment[i]=="observed") {
        for (j in seq_len(2)) {
            stopifnot(identical(cpg, observed[[j]]$cpg))
            row[[paste0("max_abs_effect_error_",j)]] <- max(abs(fit$coefficients[,j]-observed[[j]]$logFC))
            row[[paste0("max_abs_t_error_",j)]] <- max(abs(fit$t[,j]-observed[[j]]$t))
            row[[paste0("max_abs_p_error_",j)]] <- max(abs(fit$p.value[,j]-observed[[j]]$P.Value))
            row[[paste0("max_abs_bh_error_",j)]] <- max(abs(p.adjust(fit$p.value[,j],"BH")-observed[[j]]$adj.P.Val))
        }
    }
    ans[[i]] <- row
    # Preserve observed-only error columns when combining the rows.
    columns <- unique(unlist(lapply(ans,names)))
    combined <- do.call(rbind,lapply(ans,function(x) {
        for (column in setdiff(columns,names(x))) x[[column]] <- NA_real_
        x[columns]
    }))
    write.csv(combined,out,row.names=FALSE)
    cat(sprintf("Independent R fit %d/%d: %s, %d/%d red sites\n",
                i,nrow(assignments),assignments$assignment[i],row$before,row$after))
    flush.console()
}
