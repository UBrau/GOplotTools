### Tools to conduct GO analysis and/or parse and plot results from g:Profiler, FuncAssociate and DAVID
###
### User functions:
### - runGprofiler():          Run GO analysis locally with the g:GOSt function of g:Profiler
### - plotGprofilerDots():     Generate lollipop plot of results generated with runGprofiler()
###                            or from web interface. If the latter, add column 'log2Enr' with 
###                            log2 (enrichment).
### - plotFuncAssDots():       Generate lollipop plot of results obtained with FuncAssociate
### - plotFuncAssCategories(): Plot FuncAssociate enrichment scores as bar graph
### - plotDAVIDclusters():     Plot DAVID cluster scores as bar graph
### To use, source() this this file and invoke one of the functions
###
### U. Braunschweig, 2016-2024

runGprofiler <- function(
    fore, back = NULL, species, 
    sources = c("GO", "REAC", "KEGG", "TF", "CORUM", "HPA", "MIRNA", "HP", "HPA", "WP")[1:6],
    ordered = FALSE,
    exclude_iea = FALSE,
    measure_underrepresentation = FALSE,
    user_threshold = 0.05,
    outBase = "gProfilerOut"
    ) {
### Run g:GOSt from the g:Profiler suite locally, based on a list with ENSENBL gene IDs
### and optional custom background.
### 
### fore:           Vector of ENSEMBL gene IDs (in which to look for enriched categories)
###                 or list of such vectors, in case of multi-query
### back:           (Optional) custum background
### species:        Species identifier, e.g. 'mmusculus', 'hsapiens'
### sources:        Annotation sources [default: GO, KEGG]
### ordered:        Ordered query
### exclude_iea: Exclude electronically inferred annotations
### measure_underrepresentation: As the name suggests
### user_threshold: Threshold for p-value considered significant [default: 0.05]
### outBase:        Base name for saving results. Set to NA to skip saving. [default: gProfilerOut]
###
### Value: List of length two: 
###        $result: data.frame with categories, enriched with log2-fold enrichment
###        $meta: metadata
### Optionally, a CSV file with enriched categories and an R archive
###        with metadata are saved
print(sources)
    libMissing <- !require("gprofiler2", quietly=T) && stop("Failed to load R package 'gprofiler2'")

    go <- gost(
        fore,
        custom_bg      = back,
        domain_scope   = ifelse(is.null(back), "annotated", "custom"),
        organism       = species,
        exclude_iea    = exclude_iea,
        measure_underrepresentation = measure_underrepresentation,
        user_threshold = user_threshold,
        ordered_query  = ordered,
        multi_query    = is.list(fore),
        sources        = sources,
        highlight      = !ordered
    )

    if (is.null(go)) {
        write.csv(
            "No significant results", row.names = FALSE,
            file = paste0(outBase, "_NOresults.csv")
        )
        return(NULL)
    }
    go.log2Enr     <- log2(
        (go$result$intersection_size / go$result$query_size) / (go$result$term_size / go$result$effective_domain_size)
    )

    res <- data.frame(
        go$result,
        log2Enr = go.log2Enr
    )

    res <- res[order(res$significant, res$log2Enr, res$intersection_size, decreasing = TRUE),]
    res$parents <- sapply(res$parents, paste, collapse=", ")
    if (ordered) {
        res <- res[c(1, 15, 2:14)]
    } else {
        res <- res[,c(1, 16, 2:15)]
    }

    if (!is.na(outBase)) {
        write.csv(
            res, row.names = FALSE,
            file = paste0(outBase, "_results.csv")
        )
        tmp <- go$meta
        save(
            tmp,
            file = paste0(outBase, "_metadata.Rdata")
        )
    }
    ### This does not allow multi-queries, in which case the structure is different

    list(
        result = res,
        meta   = go$meta
    )
}



plotGprofilerDots <- function(over, under=NULL, outName=NA, main="",
                            minLog2Enr=log2(3), maxX=1000, maxCat=NA,
                            onlyHighlighted=TRUE, 
                            scaleNmax=NA, scalePmin=0.001, scalePmax=0.1,
                            simplePcol=TRUE, minPcol="brown1", maxPcol="indianred4",
                            sourceCol = c(
                                GO="black", KEGG="cadetblue4", REAC="bisque4", TF="coral4", CORUM="darkorchid4",
                                HPA="goldenrod4", MIRNA="darkolivegreen", HP="deepskyblue3", HPA="khaki3", WP="salmon3"
                            ),
                            circleScale=0.25, legend=TRUE, sourceLegendPos="bottomleft",
                            minXlim=c(-10,4), wid=7, hei=5) {
### Plot the log2-fold enrichment and category names from g:Profiler/g:GOSt analysis as dots,
### with location indicating log2-enrichment, size indicating number of genes, and shading adjusted p-value.
### Currently, only output with only overrepresented OR underrepresented terms is supported.
###
### over:             Over-represented terms from g:GOSt query (single query) with added 'log2Enr' column (data.frame)
### under:            (Optional) like over, for under-represented terms
### outName:          (Optional) Name of plot file (will save a pdf)
### main:             Figure main header
### minLog2Enr:       Minimum |log2(enrichment)| to report categories [default: log2(3)]
### maxX:             Maximum total number of genes associated with a term in the whole gene space for it to be reported
###                   (to remove too broad categories)
### maxCat:           Report up to this number of (over- or underrepresented) categories for each file
###                   [default:all categories]
### onlyHighlighted:  Plot only driver terms 'highlighted' by g:GOSt [default]
### scaleNmax:        Number of genes in group corresponding to largest possible circle
### scalePmin:        P value corresponding to most saturated color
### scalePmax:        P value corresponding to white
### simplePcol:       (Logical) Should a simplified P-value color scheme with two cutoffs be used rather than a scale
### minPcol:          Color to display p-values. If simplePcol=TRUE, the lower threshold.
### maxPcol:          Color to display v-values lower than the higher threshold; ingored if simplePcolor=FALSE.
### circleScale:      Scale all circles by this factor (1=100%). Useful because circle size depends on x axis.
### legend:           Plot a legend?
### sourceLegendPos:  Position of the sources legend [default: bottomleft]
### minXlim:          Minimum coordintes on x-axis. Change if text is cut off. If legend is requested, will me made at
###                   least 2 more than the largest LOD.
### wid, hei:         Width and height of the figure (in inches).

    libMissing <- !require("plotrix", quietly=T) && stop("Failed to load R package 'plotrix'")

    noCategs <- FALSE

    if (simplePcol & is.na(maxPcol)) {stop("maxPcol must be provided if simplePcol is TRUE")}
    if (!simplePcol & !is.na(maxPcol)) {warning("maxPcol is ignored if simplePcol is FALSE")}

    ## If no underrepresentation table, add dummy
    if (is.null(under)) {
        fa <- list(over, over[c(),])
    } else {
        fa <- list(over, under)
    }

    ## Parse input and restrict to plottable categories
    if (onlyHighlighted) {
        if (!all(sapply(fa, FUN=function(x) {"highlighted" %in% names(x)}))) {
            stop("No information about highlighted categories. Did you submit an ordered query?")
        }
        fa <- lapply(fa, FUN=function(x) {x[x$highlighted,]})
    }
    fa <- lapply(fa, FUN=function(x) {
                      if (nrow(x) == 0) {
                          return(x)
                      } else {
                          x <- x[abs(x$log2Enr) >= minLog2Enr,]
                          if (nrow(x) > 0) {
                              x$term_name <- .properlyUppercase(x$term_name)
                              return(x)
                          } else {
                              return(x)
                          }
                      }
    })
    names(fa) <- c("over", "under")

    exceedMaxX <- lapply(fa, FUN=function(x) {which(x$term_size > maxX)}) 
    if (sum(sapply(exceedMaxX, length)) > 0) {
        warning(paste(sum(sapply(exceedMaxX, length)), "categories with more than", maxX, "genes removed:\n"), 
            paste(fa$over$term_name[exceedMaxX[[1]]], fa$under$term_name[exceedMaxX[[2]]], collapse="\n"),
            sep="\n")
        fa <- lapply(1:length(fa), FUN=function(x) {fa[[x]][fa[[x]]$term_size <= maxX,]})
        names(fa) <- c("over", "under")
    }
    if (nrow(fa$over) > 0 & nrow(fa$under) > 0) {
        stop("Output with both over- and unerrepresented categories currently not supported")
    }

    if (!is.na(maxCat) && any(sapply(fa, nrow) > maxCat)) {
        warning("Truncating to top ", maxCat, " categories")
        fa <- lapply(fa, FUN=head, n=maxCat)
    }
    fa <- rbind(fa$over, fa$under)

    ## Calculate dot color and radius, category colors
    noCategs <- nrow(fa) < 1
    if (!noCategs) {
        fa$rad <- sqrt(fa$intersection_size)
        if (is.na(scaleNmax)) {scaleNmax <- max(fa$intersection_size)}
        fa$rad <- fa$rad  / sqrt(scaleNmax)
        
        if (simplePcol) {
            tmp <- fa$p_value
            fa$col <- "grey70"
            fa$col[tmp < 0.05] <- maxPcol
            fa$col[tmp < 0.001] <- minPcol
        } else {
            fa$col <- -log10(fa$p_value)
            if (is.na(scalePmin)) {scalePmin <- 10 ^ -max(fa$col)}
            if (is.na(scalePmax)) {scalePmax <- 10 ^ -min(fa$col)}
            palette <- colorRampPalette(c("white", minPcol))(10)
            fa$col <- (fa$col + log10(scalePmax)) / (-log10(scalePmin) + log10(scalePmax))
            fa$col <-  palette[1 + round(fa$col * (length(palette) - 1))]
        }

        fa$textCol <- "black"
        fa$textCol <- sapply(sub("(GO):.{2}", "\\1", fa$source), FUN=function(x) {sourceCol[names(sourceCol) == x]})
        sourceCol <- sourceCol[names(sourceCol) %in% sub("(GO):.{2}", "\\1", fa$source)]

        xlim.fig  <- c(min(minXlim[1], floor(min(fa$log2Enr, na.rm=T))),
                       max(minXlim[2], ceiling(max(fa$log2Enr, na.rm=T))))
    } else {
        xlim.fig <- c(-1,1)
    }

    ## Other plot preparations
    xlim.plot <- c(floor(min(c(0, fa$log2Enr), na.rm=T)),
                   ceiling(max(c(0, fa$log2Enr), na.rm=T)))
    if (legend & xlim.fig[2] < xlim.plot[2] + 2) {xlim.fig[2] <- xlim.plot[2] + 2}
    ylim <- c(0, nrow(fa) + 1)
    
    
    ## Plot
    if (!is.na(outName)) pdf(outName, wid=wid, hei=hei)
    oldpar <- par(no.readonly=T)
    
    par(mar=c(5,1,5,1))
    plot(1,1,type="n", bty="n", xlim=xlim.fig,
         ylim=ylim, yaxs="i",
         xaxt="n",yaxt="n",
         xlab="", ylab="", main = main)
    if (noCategs) {
        text(0, 0.5, "No enriched categories", col=fa$textCol)
    } else {
        for (i in 1:nrow(fa)) {
            draw.circle(x=fa$log2Enr[nrow(fa):1][i], y=i,
                        radius=circleScale * fa$rad[nrow(fa):1][i], col=fa$col[nrow(fa):1][i])
        }
        segments(x0=0, y0=(nrow(fa):1), x1=fa$log2Enr, lty=3)
        text(0, nrow(fa):1, fa$term_name, pos=ifelse(fa$log2Enr < 0, 4, 2), col = fa$textCol)
    }
    abline(v=0)
    axis(1, at=xlim.plot[1]:xlim.plot[2])
    axis(3, at=xlim.plot[1]:xlim.plot[2])
    title(xlab="log2 (enrichment)")
    
    ## Legend
    if (legend & !noCategs) {
        par(xpd=NA)
        if (nrow(fa) < 10) {
            yleg <- seq(ylim[2], ylim[1], length.out=10)
        } else {
            yleg <- ylim[2]:(ylim[2]-9)
        }

        xleg <- xlim.fig[2] - circleScale * c(2.1, 1)
        if (simplePcol) {
            legCol <- c(maxPcol, minPcol)
            legP <- c("<0.05","<0.001")
            for (i in c(1,2)) {
                draw.circle(x=xleg[2], y=yleg[2:3][i], radius=circleScale * 0.4, col=legCol[i])
            }
            text(xleg[1], yleg[2:3], adj=c(1, 0.5), labels=legP)
        } else {
            legCol <- colorRampPalette(c("white", minPcol))(4)
            legP   <- signif(10 ^ -seq(-log10(scalePmax), -log10(scalePmin), length.out=4), 1)
            if (lessFlag) legP[length(legP)] <- paste("<", legP[length(legP)], sep="")
            for (i in 1:4) {
                draw.circle(x=xleg[2], y=yleg[2:5][i], radius=circleScale * 0.4, col=legCol[i])
            }
            text(xleg[1], yleg[2:5], adj=c(1, 0.5), labels=legP)
        }

        nSizeCirc <- min(4, length(unique(fa$intersection_size)))
        legN <- round(seq(min(fa$intersection_size), scaleNmax, length.out=nSizeCirc))
        legRad <- sqrt(legN)  / sqrt(scaleNmax)
        for (i in 1:nSizeCirc) {
            draw.circle(x=xleg[2], y=yleg[7:(6 + nSizeCirc)][i], radius=circleScale * legRad[i], col=NA)
        }
        text(xleg[1], yleg[7:(6 + nSizeCirc)], adj=c(1, 0.5), labels=legN)
        text(mean(xleg),  0.95*(yleg[9] - yleg[10]) + yleg[c(2,7)], adj=c(0.5, 0.5), c("P (adj.)", "# genes"))

        if (!all(names(sourceCol) == "GO")) {
            legend(sourceLegendPos, text.col=sourceCol, names(sourceCol))
        }
        par(xpd=FALSE)        
    }
       
    par(oldpar)
    if (!is.na(outName)) dev.off()
}


plotFuncAssDots <- function(file, outName=NA, main="", minLOD=log10(3), maxX=1000, maxCat=NA,
                            mergeOverlapping=TRUE, inputGenes=NA, attrEntList=NA, mergeOl=0.7,
                            scaleNmax=NA, scalePmin=0.001, scalePmax=0.1,
                            simplePcol=TRUE, minPcol="brown1", maxPcol="indianred4",
                            circleScale=0.25, legend=T,
                            minXlim=c(-10,4), wid=7, hei=5, cores=1) {
### Plot the FuncAssociate log-odds and category names as dots,
### with location indicating LOD, size indicating number of genes, and shading adjusted p-value.
### Currently, only output with only overrepresented OR underrepresented terms is supported.
###
### file:             Path of file with enriched terms downloaded from FuncAssociate (saved 'Results')
### outName:          Optional name of plot file (will save a pdf)
### main:             Figure main header
### minLOD:           Minimum |LOD| to report categories [default: 5]
### maxX:             Maximum total number of genes associated with a term in the whole gene space for it to be reported
###                   (to remove too broad categories)
### maxCat:           Report up to this number of (over- or underrepresented) categories for each file
###                   [default:all categories]
### mergeOverlapping: Remove categories due to overlap. If there is mutual overlap of at least mergeOl, only the category
###                   with the highest log-odds will be kept. Requires inputGenes and attrEntList.
### inputGenes:       The original list of genes submitted to Funcassociate. Only required if mergeOverlapping=TRUE.       
### attrEntList:      An 'Attribute/Entity List' containing mappings of genes to GO terms, which can be downloaded from
###                   the FuncAssociate results page. Only required if mergeOverlapping=TRUE.
### mergeOl:          Minimum fraction of overlap for mergeOverlapping. If TRUE, if categories have >= minOl of the 
###		      enriched genes mutually in common, only the category with strongest enrichment is shown.
### scaleNmax:        Number of genes in group corresponding to largest possible circle
### scalePmin:        P value corresponding to most saturated color
### scalePmax:        P value corresponding to white
### simplePcol:       (Logical) Should a simplified P-value color scheme with two cutoffs be used rather than a scale
### minPcol:          Color to display p-values. If simplePcol=TRUE, the lower threshold.
### maxPcol:          Color to display v-values lower than the higher threshold; ingored if simplePcolor=FALSE.
### circleScale:      Scale all circles by this factor (1=100%). Useful because circle size depends on x axis.
### legend:           Plot a legend?
### minXlim:          Minimum coordintes on x-axis. Change if text is cut off. If legend is requested, will me made at
###                   least 2 more than the largest LOD.
### wid, hei:         Width and height of the figure (in inches).

    library(plotrix)
    library(parallel)

    noCategs <- FALSE

    if (simplePcol & is.na(maxPcol)) {stop("maxPcol must be provided if simplePcol is TRUE")}
    if (!simplePcol & !is.na(maxPcol)) {warning("maxPcol is ignored if simplePcol is FALSE")}

    if (mergeOverlapping & (is.na(inputGenes) | is.na(attrEntList))) {
        stop("If mergeOverlapping=TRUE, inputGenes and attrEntList must be provided")
    }

    ## Parse input and restrict to plottable categories
    fa <- .extractFuncAssCategories(file, min.LOD=minLOD)
    
    exceedMaxX <- lapply(fa, FUN=function(x) {which(x$X > maxX)})
    if (sum(sapply(exceedMaxX, length)) > 0) {
        warning(paste(sum(sapply(exceedMaxX, length)), "categories with more than", maxX, "genes removed:\n"), 
            paste(fa$over$attrib.name[exceedMaxX[[1]]], fa$under$attrib.name[exceedMaxX[[2]]], collapse="\n"),
            sep="\n")
        fa <- lapply(1:length(fa), FUN=function(x) {fa[[x]][fa[[x]]$X <= maxX,]})
        names(fa) <- c("over", "under")
    }
    if (nrow(fa$over) > 0 & nrow(fa$under) > 0) {
        stop("Output with both over- and unerrepresented categories currently not supported")
    }
    
    ## Get a list of genes in the input that were associated with each GO term
    inpGenes <- read.delim(inputGenes, header=T, as.is=T)[,1]
    catGenes <- read.csv(attrEntList, sep="\t", as.is=T)
    catGenes <- catGenes[match(c(as.character(fa$over$attrib.ID), as.character(fa$under$attrib.ID)),
                               catGenes$Significant.attribute),]
    if (class(catGenes) != "data.frame" ||
        ncol(catGenes)  != 3 ||
        !all(c("Significant.attribute", "Significant.attribute.description", "Entity.list") == names(catGenes)))
    {
        stop("Are you sure you are using the Attribute/Entity List (and not e.g. the Entity/Attribute List) from FuncAssociate?")
    }
    
    catGenes <- list(cat   = catGenes[,1],
                     name  = catGenes[,2],
                     genes = strsplit(as.character(catGenes[,3]), split=" ")
                     )            
    catGenes$genes <- lapply(catGenes$genes, FUN=function(x) {intersect(x, inpGenes)})

    
    ## If categories mutually overlap by more than set fraction, remove less significant one
    if (nrow(fa$over) + nrow(fa$under) > 0) {
        fa <- .reduceOlCategories(fa, inputGenes, catGenes, mergeOl, maxCat, mergeOverlapping, cores)
    }
    if (nrow(fa) < 1) noCategs <- TRUE

    ## Calculate dot color and radius
    if (!noCategs) {
        fa$rad <- sqrt(fa$N)
        if (is.na(scaleNmax)) {scaleNmax <- max(fa$N)}
        fa$rad <- fa$rad  / sqrt(scaleNmax)
        
        if (simplePcol) {
            tmp <- as.character(fa$P_adj)
            tmp[tmp == "<0.001"] <- "0"
            tmp <- as.numeric(as.character(tmp))
            fa$col <- "grey70"
            fa$col[tmp < 0.05] <- maxPcol
            fa$col[tmp < 0.001] <- minPcol
        } else {
            fa$col <- as.character(fa$P_adj)
            lessFlag <- ifelse(any(fa$col %in% c("0", "<0.001")), TRUE, FALSE)
            fa$col[fa$col == "<0.001"] <- "0.001"
            fa$col[fa$col == 0] <- "0.001"
            fa$col <- -log10(as.numeric(as.character(fa$col)))
            if (is.na(scalePmin)) {scalePmin <- 10 ^ -max(fa$col)}
            if (is.na(scalePmax)) {scalePmax <- 10 ^ -min(fa$col)}
            palette <- colorRampPalette(c("white", minPcol))(10)
            fa$col <- (fa$col + log10(scalePmax)) / (-log10(scalePmin) + log10(scalePmax))
            fa$col <-  palette[1 + round(fa$col * (length(palette) - 1))]
        }

        xlim.fig  <- c(min(minXlim[1], floor(min(fa$LOD, na.rm=T))),
                       max(minXlim[2], ceiling(max(fa$LOD, na.rm=T))))

    } else {
        xlim.fig=c(-1,1)
    }

    ## Other plot preparations
    xlim.plot <- c(floor(min(c(0, fa$LOD), na.rm=T)),
                   ceiling(max(c(0, fa$LOD), na.rm=T)))
    if (legend & xlim.fig[2] < xlim.plot[2] + 2) {xlim.fig[2] <- xlim.plot[2] + 2}
    ylim <- c(0, nrow(fa) + 1)
    
    
    ## Plot
    if (!is.na(outName)) pdf(outName, wid=wid, hei=hei)
    oldpar <- par(no.readonly=T)
    
    par(mar=c(5,1,5,1))
    plot(1,1,type="n", bty="n", xlim=xlim.fig,
         ylim=ylim, yaxs="i",
         xaxt="n",yaxt="n",
         xlab="", ylab="", main = main)
    if (noCategs) {
        text(0, 0.5, "No enriched categories")
    } else {
        for (i in 1:nrow(fa)) {
            draw.circle(x=fa$LOD[nrow(fa):1][i], y=i,
                        radius=circleScale * fa$rad[nrow(fa):1][i], col=fa$col[nrow(fa):1][i])
        }
        segments(x0=0, y0=(nrow(fa):1), x1=fa$LOD, lty=3)
        text(0, nrow(fa):1, fa$attrib.name, pos=ifelse(fa$LOD < 0, 4, 2))
    }
    abline(v=0)
    axis(1, at=xlim.plot[1]:xlim.plot[2])
    axis(3, at=xlim.plot[1]:xlim.plot[2])
    title(xlab="log10 (odds ratio)", adj=abs(xlim.fig[1]) / sum(abs(xlim.fig)))

    ## Legend
    if (legend & !noCategs) {
        par(xpd=NA)
        if (nrow(fa) < 10) {
            yleg <- seq(ylim[2], ylim[1], length.out=10)
        } else {
            yleg <- ylim[2]:(ylim[2]-9)
        }

        xleg <- xlim.fig[2] - circleScale * c(2.1, 1)
        if (simplePcol) {
            legCol <- c(maxPcol, minPcol)
            legP <- c("<0.05","<0.001")
            for (i in c(1,2)) {
                draw.circle(x=xleg[2], y=yleg[2:3][i], radius=circleScale * 0.4, col=legCol[i])
            }
            text(xleg[1], yleg[2:3], adj=c(1, 0.5), labels=legP)
        } else {
            legCol <- colorRampPalette(c("white", minPcol))(4)
            legP   <- signif(10 ^ -seq(-log10(scalePmax), -log10(scalePmin), length.out=4), 1)
            if (lessFlag) legP[length(legP)] <- paste("<", legP[length(legP)], sep="")
            for (i in 1:4) {
                draw.circle(x=xleg[2], y=yleg[2:5][i], radius=circleScale * 0.4, col=legCol[i])
            }
            text(xleg[1], yleg[2:5], adj=c(1, 0.5), labels=legP)
        }

        nSizeCirc <- min(4, length(unique(fa$N)))
        legN <- round(seq(min(fa$N), scaleNmax, length.out=nSizeCirc))
        legRad <- sqrt(legN)  / sqrt(scaleNmax)
        for (i in 1:nSizeCirc) {
            draw.circle(x=xleg[2], y=yleg[7:(6 + nSizeCirc)][i], radius=circleScale * legRad[i], col=NA)
        }
        text(xleg[1], yleg[7:(6 + nSizeCirc)], adj=c(1, 0.5), labels=legN)
        text(mean(xleg),  0.95*(yleg[9] - yleg[10]) + yleg[c(2,7)], adj=c(0.5, 0.5), c("P (adj.)", "# genes"))
        par(xpd=FALSE)        
    }
       
    par(oldpar)
    if (!is.na(outName)) dev.off()
}



plotFuncAssCategories <- function(fileUp=NA, fileDown=NA, outName=NA, main="", min.LOD=log10(5), maxCat=NA,
                                  invertNegOrder=TRUE,
                                  minXlim=c(-7,7), wid=7, hei=5,
                                  colUp="brown1", colDown="dodgerblue") {
### Plot the FuncAssociate log-odds and category names as bar graphs
###
### fileUp:         Path of file with "upregulated" terms
### fileDown:       Path of file with "downregulated" terms
### outName:        Optional name of plot file (will save a pdf)
### main:           Figure main header
### min.LOD:        Minimum |LOD| to report categories [default: 5]
### maxCat:         Report up to this number of (over- or underrepresented) categories for each file [default:all categories]
### invertNegOrder: Invert the order of categories in fileDown (only if fileUp is provided); logical
### minXlim:        Minimum coordintes on x-axis. Change if text is cut off.
### wid, hei:       Width and height of the figure (in inches).

    if (is.na(fileUp) & is.na(fileDown)) {stop("Need at least one input file!")}
    noCategs <- FALSE

    if (is.na(fileUp)) {
        n.fa.plus <- c(NA)
        fa.plus <- c()
    } else {
        fa.plus <- .extractFuncAssCategories(fileUp, min.LOD=min.LOD, maxCat=maxCat)
        fa.plus <- lapply(fa.plus, FUN=function(x) {matrix(x$LOD, nrow=1, dimnames=list(c(), x$attrib.name))})
        n.fa.plus <- sapply(fa.plus, ncol)
        if (n.fa.plus[2] == 0) {            
            fa.plus <- fa.plus[[1]]
        } else {
            fa.plus <- cbind(fa.plus[[1]], t(as.matrix(fa.plus[[2]][,ncol(fa.plus[[2]]):1])))
        }
    }

    if (is.na(fileDown)) {
        n.fa.minus <- c(NA)
        fa.minus <- c()
    } else {
        fa.minus <- .extractFuncAssCategories(fileDown, min.LOD=min.LOD, maxCat=maxCat)
        fa.minus <- lapply(fa.minus, FUN=function(x) {matrix(x$LOD, nrow=1, dimnames=list(c(), x$attrib.name))})
        n.fa.minus <- sapply(fa.minus, ncol)
        if (n.fa.minus[2] == 0) {            
            fa.minus <- fa.minus[[1]]
        } else {
            fa.minus <- cbind(fa.minus[[1]], t(as.matrix(fa.minus[[2]][,ncol(fa.minus[[2]]):1])))
        }
        if (!is.na(fileUp) & invertNegOrder & sum(n.fa.minus, na.rm=T) > 0) {fa.minus <- t(as.matrix(fa.minus[,ncol(fa.minus):1]))}
    }


    if (!is.na(outName)) pdf(outName, wid=wid, hei=hei)
    xlim.fig  <- c(min(minXlim[1], floor(min(c(fa.plus, fa.minus), na.rm=T))),
                   max(minXlim[2], ceiling(max(c(fa.plus, fa.minus), na.rm=T))))
    xlim.plot <- c(floor(min(c(0, fa.plus, fa.minus), na.rm=T)),
                   ceiling(max(c(0, fa.plus, fa.minus), na.rm=T)))
    oldpar <- par(no.readonly=T)
    upperOffset <- ifelse(is.na(fileUp) | is.na(fileDown), 0, sum(n.fa.plus, n.fa.minus, na.rm=T)/20)

    par(mar=c(5,1,4,1))
    plot(1,1,type="n", bty="n", xlim=xlim.fig,
         ylim=c(0.2, length(fa.plus) + length(fa.minus) + upperOffset + 0.6), yaxs="i",
         xaxt="n",yaxt="n",
         xlab="", ylab="", main = main)
    axis(1, at=xlim.plot[1]:xlim.plot[2])
    title(xlab="log10 (odds ratio)", adj=abs(xlim.fig[1]) / sum(abs(xlim.fig)))
    abline(v=0)
    
    if (!is.na(fileUp)) {
        if (sum(n.fa.plus, na.rm=T) > 0) {
            rect(0,       sum(n.fa.plus, n.fa.minus, na.rm=T):(1 + sum(n.fa.minus, na.rm=T)) - 0.4 + upperOffset,
                 fa.plus, sum(n.fa.plus, n.fa.minus, na.rm=T):(1 + sum(n.fa.minus, na.rm=T)) + 0.4 + upperOffset, col=colUp)
            if (n.fa.plus[1] > 0) {
                text(-0.1, (sum(n.fa.plus, n.fa.minus, na.rm=T):(1 + sum(n.fa.minus, na.rm=T)) + upperOffset)[fa.plus > 0],
                     dimnames(fa.plus)[[2]][fa.plus > 0], pos=2, cex=1)
            }
            if (n.fa.plus[2] > 0) {
                text(0.1, (sum(n.fa.plus, n.fa.minus, na.rm=T):(1 + sum(n.fa.minus, na.rm=T)) + upperOffset)[fa.plus < 0],
                     dimnames(fa.plus)[[2]][fa.plus < 0], pos=4, cex=1)
            }
            axis(3, at=xlim.plot[1]:xlim.plot[2], labels=NA)
        } else {
            noCategs <- TRUE
        }
    }
    
    if (!is.na(fileDown)) {
        if (sum(n.fa.minus, na.rm=T) > 0) {
            rect(0, sum(n.fa.minus):1 - 0.4, fa.minus, sum(n.fa.minus):1 + 0.4, col=colDown)
            if (n.fa.minus[1] > 0) {
                text(-0.1, (sum(n.fa.minus):1)[fa.minus > 0], dimnames(fa.minus)[[2]][fa.minus > 0], pos=2, cex=1)
            }
            if (n.fa.minus[2] > 0) {
                text(0.1, (sum(n.fa.minus):1)[fa.minus < 0], dimnames(fa.minus)[[2]][fa.minus < 0], pos=4, cex=1)
            }
            segments(x0=xlim.plot[1], y0=sum(n.fa.minus) + 0.5 + 0.5 * upperOffset, x1=xlim.plot[2])
            noCategs <- FALSE
        } else {
            noCategs <- noCategs & TRUE
        }
    }

    if (noCategs) {text(0.1, 0.5, "No siginficant categories", pos=4)}

    par(oldpar)
    if (!is.na(outName)) dev.off()
}


.extractFuncAssCategories <- function(file, min.LOD=log10(2), maxCat=NA) {
### Called by plotFuncAssCategories() and plotFuncAssDots()
### Given the name of a file that contains FuncAssociate output (looking for over- and underrepresentation),
### return categories with at least +/- min.LOD
###
### Parameters:
###   file:    Name of the file
###   min.LOD: Minimum |LOD| to report categories (+LOD means overrepresented, -LOD underrepresented) [default: 5]
###   maxCat:  Report up to this number of categories of both over- and underrepresented [default:all categories]
### Value:
###   list of two data.frame() with FuncAssociate tables
    
    dat    <- scan(file, what="character", sep="\n", blank.lines.skip=F)
    overHeadInd  <- grep("OVERREPRESENTED ATTRIBUTES", dat)
    underHeadInd <- grep("UNDERREPRESENTED ATTRIBUTES", dat)
    if (length(underHeadInd) == 0) underHeadInd <- NA
    overEnd  <- ifelse(!is.na(underHeadInd), underHeadInd - 2, length(dat))
    underEnd <- ifelse(!is.na(underHeadInd), length(dat), NA)

    if (overEnd - overHeadInd > 1) {
        over <- read.delim(file, skip=overHeadInd, nrows=overEnd - overHeadInd - 1)
    } else {
        over <- read.delim(file, skip=overHeadInd, nrows=1, blank.lines.skip=F)[c(),]
    }

    if (!is.na(underHeadInd) && underEnd - underHeadInd > 1) {
        under <- read.delim(file, skip=underHeadInd, nrows=underEnd - overHeadInd - 1)
    } else {
        under <- read.delim(file, skip=overHeadInd, nrows=1, blank.lines.skip=F)[c(),]
    }
    
    dat <- list(over=over , under=under)
    
    dat <- lapply(dat, FUN=function(x) {
                      if (length(x) == 1 && is.na(x)) {
                          return(x)
                      } else {
                          x <- x[abs(x$LOD) >= min.LOD,]
                          if (nrow(x) > 0) {
                              x$attrib.name <- .properlyUppercase(x$attrib.name)
                              return(head(x, ifelse(is.na(maxCat), nrow(x), maxCat)))
                          } else {
                              return(x)
                          }
                      }
    })
    if (all(sapply(dat, nrow) == 0)) {
        stop("No categories above minLOD")
    }
    dat
}



.reduceOlCategories <- function(fa, inpGenes, catGenes, mergeOl, maxCat, mergeOverlapping, cores=1) {
### Called by plotFuncAssCategories/plotFuncAssDots to merge and reduce the number of categories for plotting
    tmp <- rbind(fa$over, fa$under)

    if (mergeOverlapping && nrow(fa$over) + nrow(fa$under) > 1) {
        catOl <- mcmapply(1:length(catGenes$cat), FUN=function(x) {sapply(1:length(catGenes$cat), FUN=function(y) {
            length(which(catGenes$genes[[x]] %in% catGenes$genes[[y]]))
        })}, mc.cores=cores)
        catLink <- catOl / sapply(1:length(catGenes$genes), FUN=function(x) {sapply(catGenes$genes, length)}) >= mergeOl
        diag(catLink) <- FALSE
        nets <- lapply(1:nrow(catLink), FUN=function(x) {which((catLink & t(catLink))[x,])})  # mutual overlap
        nets <- lapply(1:length(nets), FUN=function(x) {
            if (length(nets[[x]] > 0)) {
                return(sort(c(x, nets[[x]])))
            } else {
                return(NA)
            }
        })
        nets <- nets[sapply(nets, FUN=function(x) {!all(is.na(x))})]
        nets <- nets[!duplicated(sapply(nets, FUN=function(x) {paste(x, collapse=" ")}))]
        remove <- data.frame(keep   = 1:nrow(catOl) %in% unlist(lapply(nets, FUN=function(x) {x[which.max(tmp$LOD[x])]})),
                             remove = 1:nrow(catOl) %in% unlist(lapply(nets, FUN=function(x) {setdiff(x, x[which.max(tmp$LOD[x])])}))
                             )
        remove <- remove$remove & !remove$keep  # remove categories only if they are not supposed to represent another cluster
        
        if (any(remove)) {
            warning(paste(length(which(remove)), "categories removed due to overlap:\n"), 
                    paste(tmp$attrib.name[remove], collapse="\n"),
                    sep="\n")
        }
        
        fa$over  <- fa$over[!(fa$over$attrib.ID %in% tmp$attrib.ID[remove]),]
        fa$under <- fa$under[!(fa$under$attrib.ID %in% tmp$attrib.ID[remove]),]
    }

    exceedMaxCat <- sapply(fa, nrow) - maxCat
    if (!is.na(exceedMaxCat[1]) && exceedMaxCat[1] > 0) {
        warning(exceedMaxCat[1], " overrepresented categories exceeded maxCat\n")
        fa$over <- fa$over[1:maxCat,]
    }
    if (!is.na(exceedMaxCat[2]) && exceedMaxCat[2] > 0) {
        warning(exceedMaxCat[2], " underrepresented categories exceeded maxCat\n")
        fa$under <- fa$under[1:maxCat,]
    }

    if (nrow(fa$under) > 1) {fa$under <- fa$under[nrow(fa$under):1,]}
    rbind(fa$over, fa$under)
}



plotDAVIDclusters <- function(fileUp, fileDown, outName=NA, main="", minXlim=c(-7,7), wid=7, hei=5, min.p=0.01,
                              colUp="brown1", colDown="dodgerblue") {
### Plot the DAVID cluster scores and category names extracted with extractCategories()
###
### fileUp:   Path of file with "upregulated" terms
### fileDown: Path of file with "downregulated" terms
### outName:  Optional name of plot file (will save a pdf)
### main:     Figure main header
### minXlim:  Minimum coordintes on x-axis. Change if text is cut off.
### wid, hei: width and height of the figure (in inches).
### min.p:    Select clusters in which at least one category has a min. Benjamini-p-value lower than this
    
    david.plus <- .extractDAVIDcategories(fileUp, min.p=min.p)
    david.plus <- matrix(david.plus$score, nrow=1, dimnames=list(c(), david.plus$category))
    if (ncol(david.plus) == 0) {david.plus <- matrix(0)}

    david.minus <- .extractDAVIDcategories(fileDown, min.p=min.p)
    david.minus <- matrix(-1 * david.minus$score, nrow=1, dimnames=list(c(), david.minus$category))
    if (ncol(david.minus) == 0) {david.minus <- matrix(0)}
    
    david <- cbind(david.plus,david.minus)
    
    xlim.fig  <- c(min(minXlim[1],floor(min(david, na.rm=T))), max(minXlim[2],ceiling(max(david, na.rm=T))))
    xlim.plot <- c(floor(min(0, david, na.rm=T)), ceiling(max(0, david, na.rm=T)))
    
    if (!is.na(outName)) pdf(outName, wid=wid, hei=hei)
    oldpar <- par(no.readonly=T)
    par(mar=c(5,1,4,1))
    plot(1,1,type="n", bty="n", xlim=xlim.fig, ylim=c(0.6, ncol(david) + 1), xaxt="n",yaxt="n",
         xlab="DAVID cluster score", ylab="", main = main)
    axis(1, at=xlim.plot[1]:xlim.plot[2], labels=abs(c(xlim.plot[1]:xlim.plot[2])))
    abline(v=0)
    rect(0, ncol(david):(ncol(david.minus) + 1) - 0.1, david.plus, ncol(david):(ncol(david.minus) + 1) + 0.7,
         col=rep(colUp, ncol(david.plus)))
    rect(0, 1:ncol(david.minus) - 0.4, david.minus, 1:ncol(david.minus) + 0.4,
         col=rep(colDown, ncol(david.minus)))
    text(-0.1, ncol(david):(ncol(david.minus) + 1) + 0.3, dimnames(david.plus)[[2]], pos=2, cex=1.1)
    text(0.1, 1:ncol(david.minus), dimnames(david.minus)[[2]], pos=4, cex=1.1)
    par(oldpar)
    if (!is.na(outName)) dev.off()
}


.extractDAVIDcategories <- function(file, min.p=0.01) {
### Called by plotDAVIDclusters()
### Given the name of a file that contains DAVID output (clustering tool),
### return only clusters with at least one category below a certain Benjamini p-value.
### The category within the cluster with the lowest (Bonferroni corrected) p-value
### is chosen as the label.
###
### Value: data.frame() with slots "category", "score" (DAVID cluster score)

    dat    <- scan(file, what="character", sep="\n")
    headl  <- grep("Annotation Cluster [0-9]+", dat)
    headl <- c(headl, length(dat)+1) # add one more to have an end for the last record
    score  <- as.numeric(sub(".*Score: ([0-9.E-]+)$", "\\1", dat[headl[-length(headl)]]))
    catlines <- lapply(1:(length(headl) - 1), FUN=function(x) {(headl[x]+2):(headl[x + 1] - 1)})
    topCat  <- sapply(catlines, FUN=function(x) {
        tmp <- strsplit(dat[x], split="\t")
        min.pBonf <- which.min(sapply(tmp, FUN=function(x) {x[11]}))
        label <- sub(".+[~:]","", sapply(tmp, FUN=function(x) {x[2]}))[min.pBonf]
        .properlyUppercase(label)
    })
    min.pBenj <- sapply(catlines, FUN=function(x) {
        tmp <- strsplit(dat[x], split="\t")
        min(sapply(tmp, FUN=function(x) {as.numeric(x[12])}))
    })

    data.frame(category = topCat,
               score    = score
               )[min.pBenj < min.p,]
}


.properlyUppercase <- function(x) {
### Called by other functions; creates plot-ready capitalization of GO terms
    x <- as.character(x)
    firstWord <- unlist(sapply(strsplit(x, split="[ ]+"), FUN="[[", 1))
    anyUpper <- firstWord != tolower(firstWord)
    ifelse(anyUpper, x, paste(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)), sep=""))
}

