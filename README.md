## SwitchSeq
When working with RNA-seq data, several tools exist to quantify differences in splicing across conditions and to address the significance of those changes (e.g. [DEXSeq](http://www.bioconductor.org/packages/release/bioc/html/DEXSeq.html), [MMDIFF](http://www.ncbi.nlm.nih.gov/pubmed/24281695)). Quiet often though, these tools result in a long list of genes that is difficult to interpret. By relying on transcript level quantifications, SwitchSeq provides a simple (yet powerful) approach to identify, annotate and visualise the most extreme changes in splicing across two different conditions, namely **switch events**. In brief, switch events are defined as those cases where, for a given gene, the identity of the most abundant transcript changes across conditions:

![SwitchSeq overview](/doc/fig1.png)


## Quick navigation
* [Setup guide](https://github.com/mgonzalezporta/switchseq/wiki/Setup-guide)
* [Tutorial](https://github.com/mgonzalezporta/switchseq/wiki/Tutorial)
* [Software usage](https://github.com/mgonzalezporta/switchseq/wiki/Software-usage)

## License
SwitchSeq is licensed under GPLv2.

## Contact
Mar Gonzàlez-Porta:
<mar@ebi.ac.uk>
