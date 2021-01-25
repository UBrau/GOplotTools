source("GOplotTools.R")

plotFuncAssDots(file="input/FuncAssociate_results.tsv",
                outName="output/FuncAssociate_plotFuncAssDots.pdf",
                inputGenes="input/Input_genesWithChangingExons.txt",
                attrEntList="input/FuncAssociate_attrEntList.xls",
                main="Genes enriched in input over background",
                wid=9, hei=4.5)
