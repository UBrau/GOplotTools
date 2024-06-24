# GOplotTools
## Tools to conduct functional enrichment analysis, parse and plot output from gene ontology (GO) analysis (web) tools

### Dependencies
- Any version of R
- CRAN packags _gprofiler2_ (for local g:Profiler analysis), _plotrix_ and optionally _parallel_

### Workflow to conduct and plot GO analysis using g:Profiler or FuncAssociate

#### 1. Generate a list of (foreground) genes and a suitable background

A suitable background often is not the whole genome but should represent the set of genes that *could* be found in the analysis, e.g. because they contain alternative splicing events or are expressed in the control sample. It is better to use a unique gene identifier such as the ENCODE gene ID rather than the gene name. Gene names change.

#### 2a. Upload to g:Profiler, FuncAssociate (or DAVID) and run appropriate analysis** 

Find [g:Profiler](https://biit.cs.ut.ee/gprofiler/gost), [FuncAssociate](http://llama.mshri.on.ca/funcassociate/) and [DAVID](https://david.ncifcrf.gov/).

If using the web interface of g:Profiler, add a column _log2Enr_ with the log2 enrichment.
If using FuncAssociate, save the results table as well as the 'Attribute/Entity List'.  
Examples can be found in the *input* folder. 

#### OR

#### 2b. Run g:Profiler analysis locally

using runGprofiler(), which uses the _gProfiler2_ package (Kolberg et al., _F1000Research_ 2020).

#### 3. Plot results using one of the plotting functions

For examples of 'lollipop' plots based on g:Profiler and FuncAssociate see _output_ folder.

**g:Profiler**
Huge categories will be removed and only 'highlighted' driver categories will be shown by default. See options. Remaining categories will be ploted such that log2-enrichment is on the x-axis, dot size represents the number of genes from the category that were in the foreground, and color reflects p-value. Sources are indicated by text color.

~~~~
source("GOplotTools.R")

foreground <- read.delim("input/Input_genesWithChangingExons.txt", header=FALSE)[,1]
background <- read.delim("input/Input_background.txt", header=FALSE)[,1]

enriched <- runGprofiler(fore = foreground, back = background, species = "mmusculus", outBase = NA)

plotGprofilerDots(over    = enriched$result,
                  main    = "Genes enriched in input over background",
                  outName = "output/gProfiler_plotGprofilerDots.pdf", 
                  wid = 7, hei = 3.5)
~~~~

**FuncAssociate:**
Huge categories will be removed by default. If categories overlap more than a threshold, only the more significant one will be kept. See options. Remaining categories will be plotted such that the LOD is on the x-axis, dot size represents the number of genes from the category that were in the foreground, and color reflects p-value:

~~~~
source("GOplotTools.R")

plotFuncAssDots(file        = "input/FuncAssociate_results.tsv",
                outName     = "output/FuncAssociate_plotFuncAssDots.pdf",
                inputGenes  = "input/Input_genesWithChangingExons.txt",
                attrEntList = "input/FuncAssociate_attrEntList.xls",
                main = "Genes enriched in input over background",
                wid = 9, hei = 4.5)
~~~~


#### _Disclaimer:_
Input checking is imperfect. If it doesn't work, check that you are submitting the correct files.
