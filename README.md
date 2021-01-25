# GOplotTools
### Tools to parse and plot output from gene ontology (GO) analysis web tools

#### Dependencies
- Any version of R
- CRAN package _plotrix_ and optionally _parallel_

#### Workflow to conduct and plot GO analysis using FuncAssociate

1. Generate a list of (foreground) genes and a suitable background
2. Upload to DAVID or FuncAssociate and run appropriate analysis
3. Save the result table as well as the 'entity attribute list'
4. Plot results using plotFuncAssDots()

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

See plot in _output_ folder for result.

#### _Disclaimer:_
Input checking is imperfect. If it doesn't work, check that you are submitting the correct files.
