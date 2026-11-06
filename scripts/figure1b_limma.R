# Run from the project root. This is a candidate reconstruction, not author code.
.libPaths(c(file.path(getwd(), ".tools/r-library"), .libPaths()))
suppressPackageStartupMessages(library(limma))
stopifnot(as.character(packageVersion("limma")) == "3.54.0")
args <- commandArgs(trailingOnly=TRUE)
field <- if (length(args)) args[1] else "diet_code"
stopifnot(field %in% c("diet_code", "diet_named"))
scenario <- if (length(args)>=3) args[3] else "blocked_randomization"
stopifnot(scenario %in% c("blocked_randomization","corrected_visits","deposited"))
if (scenario=="blocked_randomization") stopifnot(field=="diet_code")
base <- if (scenario!="deposited") file.path("results",scenario) else "results"
sample_file <- if (scenario!="deposited") "samples_corrected_visits.csv" else "samples.csv"
out <- file.path(base, field)
shuffled_example <- length(args) >= 2 && nzchar(args[2])
if (shuffled_example) out <- file.path(out, "shuffled_example")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
samples <- read.csv(file.path("data/processed",sample_file), colClasses=c(participant="character", diet_code="character"))
if (shuffled_example) {
    assignment <- read.csv(args[2], colClasses="character")
    stopifnot(!anyDuplicated(assignment$participant),
              setequal(assignment$participant, unique(samples$participant)))
    # Match participants explicitly: both visits receive the same shuffled label.
    samples[[field]] <- assignment$shuffled_group[match(samples$participant, assignment$participant)]
}
cpg <- readLines("data/processed/cpg_ids.txt")
beta <- matrix(readBin("data/processed/beta.bin", "double", n=length(cpg)*nrow(samples), size=8, endian="little"), nrow=length(cpg), byrow=TRUE)
stopifnot(ncol(beta) == nrow(samples), all(is.finite(beta)), all(beta > 0 & beta < 1))
rownames(beta) <- cpg
# M = log2(beta/(1-beta)). The coefficient is a log2 methylation-odds ratio,
# not log2(mean(beta_group0)/mean(beta_group1)).
M <- log2(beta / (1-beta))
rm(beta)
groups <- if (field == "diet_code") c("0", "1") else c("Vegan Diet", "Meat-rich Diet")
group <- as.integer(samples[[field]] == groups[1])
stopifnot(all(samples[[field]] %in% groups))
cell <- factor(paste0("g", group, "t", samples$time), levels=c("g0t0", "g1t0", "g0t1", "g1t1"))
design <- model.matrix(~0+cell)
colnames(design) <- levels(cell)
# The five contrasts supplied by the author, retaining numeric group names.
contrasts <- makeContrasts(change_numerator=g1t1-g1t0, change_denominator=g0t1-g0t0,
    before=g1t0-g0t0, after=g1t1-g0t1,
    interaction=(g1t1-g1t0)-(g0t1-g0t0), levels=design)
fit <- eBayes(contrasts.fit(lmFit(M, design), contrasts), trend=FALSE, robust=FALSE)
decisions <- decideTests(fit, adjust.method="none", p.value=0.05, lfc=0.5)
summary <- list()
for (j in seq_len(2)) {
    time <- c("before", "after")[j]
    tab <- topTable(fit, coef=time, number=Inf, sort.by="none")
    tab$cpg <- rownames(tab)
    tab$highlighted <- tab$P.Value < 0.05 & abs(tab$logFC) > 0.5
    stopifnot(identical(as.logical(tab$highlighted),as.logical(decisions[,time]!=0)))
    con <- gzfile(file.path(out, paste0(time, "_limma.csv.gz")), "wt")
    write.csv(tab, con, row.names=FALSE)
    close(con)
    summary[[j]] <- data.frame(time=time, numerator=groups[1], denominator=groups[2],
        n_numerator=sum(group==1 & samples$time==j-1), n_denominator=sum(group==0 & samples$time==j-1),
        highlighted=sum(tab$highlighted), bh_below_005=sum(tab$adj.P.Val < 0.05),
        df_residual=fit$df.residual[1], df_prior=fit$df.prior,
        s2_prior=fit$s2.prior, df_total=fit$df.total[1])
}
write.csv(do.call(rbind, summary), file.path(out, "limma_summary.csv"), row.names=FALSE)
# A fit at each time point is a prespecified sensitivity to the unspecified design.
if (!shuffled_example) {
sensitivity <- list()
for (time in 0:1) {
    idx <- which(samples$time==time)
    d <- model.matrix(~group[idx])
    f <- eBayes(lmFit(M[,idx], d), trend=FALSE, robust=FALSE)
    sensitivity[[time+1]] <- data.frame(time=time, highlighted=sum(f$p.value[,2]<0.05 & abs(f$coefficients[,2])>0.5), bh_below_005=sum(p.adjust(f$p.value[,2], "BH")<0.05))
}
write.csv(do.call(rbind, sensitivity), file.path(out, "separate_time_sensitivity.csv"), row.names=FALSE)
}
writeLines(sub("[[:space:]]+$", "", capture.output(sessionInfo())), file.path(out, "R_session.txt"))
print(do.call(rbind, summary))
