# GOplotTools
### Tools to parse and plot output from gene ontology (GO) analysis web tools

#### Dependencies
- Any version of R
- CRAN package _plotrix_ and optionally _parallel_

#### Workflow to conduct and plot GO analysis using FuncAssociate

**1. Generate a list of (foreground) genes and a suitable background**

A suitable background often is not the whole genome but should represent the set of genes that *could* be found in the analysis, e.g. because they contain alternative splicing events or are expressed in the control sample. It is better to use a unique gene identifier such as the ENCODE gene ID rather than the gene name. Gene names change.

**2. Upload to FuncAssociate (or DAVID) and run appropriate analysis** 

Examples can be found in the *input* folder. Find [FuncAssociate](http://llama.mshri.on.ca/funcassociate/) and [DAVID](https://david.ncifcrf.gov/).

**3. Save the results table as well as the 'Attribute/Entity List'** 

if using FuncAssociate. Examples in the *input* folder.

**4. Plot results using *plotFuncAssDots()***

Huge categories will be removed. If categories overlap more than a threshold, only the more significant one will be kept. Remaining categories will be plotted such that the LOD is on the x-axis, dot size represents the number of genes from the category that were in the foreground, and color reflects p-value:

~~~~
source("GOplotTools.R")

plotFuncAssDots(file="input/FuncAssociate_results.tsv",
                outName="output/FuncAssociate_plotFuncAssDots.pdf",
                inputGenes="input/Input_genesWithChangingExons.txt",
                attrEntList="input/FuncAssociate_attrEntList.xls",
                main="Genes enriched in input over background",
                wid=9, hei=4.5)
~~~~

See plot in *output* folder for result.

#### _Disclaimer:_
Input checking is imperfect. If it doesn't work, check that you are submitting the correct files.
