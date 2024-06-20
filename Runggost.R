### Given input and background tables, run g:gost from g:Profiler


fore <- read.delim("~/Code/R/scripts/GOplotTools/input/Input_genesWithChangingExons.txt", header=F)[,1]
back <- read.delim("~/Code/R/scripts/GOplotTools/input/Input_background.txt", header=F)[,1]

rungGost <- function(
    fore, back, species, 
    sources = c("GO","KEGG","REAC","TF","CORUM","HPA"),
    minIntSize = 2,
    maxTermSize = 1000
    ) {

    libMissing <- !require("gprofiler2", quietly=T) && stop("Failed to load R package 'gprofiler2'")

    go <- gost(
        fore,
        custom_bg     = back,
        domain_scope  = ifelse(is.null(back), "annotated", "custom"),
        organism      = species,
        ordered_query = FALSE,
        multi_query   = FALSE,
        sources=c("GO", "REAC", "KEGG", "TF", "CORUM", "HPA")[1:2],
        highlight = T
    )

    gres <- go$result
    go.pvals    <- t(sapply(go$result$p_values, c))
    go.sig      <- t(sapply(go$result$significant, c))
    go.qsz      <- t(sapply(go$result$query_sizes, c))
    go.intersz  <- t(sapply(go$result$intersection_sizes, c))
    go.hilight  <- t(sapply(go$result$highlighted, c))

    go.enr     <- log2((go.intersz / go.qsz) / (go$result$term_size / go$result$effective_domain_size))
    go.enr[go.intersz < minIntSize] <- 0
    go.enr[go.enr < 0] <- 0


}