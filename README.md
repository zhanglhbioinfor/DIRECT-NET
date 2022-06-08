# DIRECT-NET: Inferring CREs and constructing regulatory networks from single cell multi-omics data

## Capabilities
DIRECT-NET is an R toolkit for inferring CREs and constructing regulatory networks from parallel single cell RNA-seq and ATAC-seq data or only single cell ATAC-seq data. 

## Installation
To make it easy to run DIRECT-NET in most common scRNA-seq and scATAC-seq data analysis pipelines, DIRECT-NET is now implemented within Seurat V4/Signac workflow. Please first **[install Seurat R pacakge (>= 4.0)](https://satijalab.org/seurat/install.html)** via ```install.packages('Seurat')```. 

DIRECT-NET R package can then be easily installed from Github using devtools:  

```
devtools::install_github("zhanglhbioinfor/DIRECT-NET")
```
 
### Installation of other dependencies
- Install Signac pacakge : ```devtools::install_github("timoast/signac", ref = "develop")```. Please check [here](https://satijalab.org/signac/articles/install.html#development-version-1) if you encounter any issue.
- Install Cicero package: ```devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")```. Please check [here](https://cole-trapnell-lab.github.io/cicero-release/docs_m3/#installing-cicero) if you encounter any issue.

## Tutorials
The implementent of DIRECT-NET is now seamlessly compatible with the workflow of Seurat V4 package. 

## Examples

- [Demo of DIRECT-NET on parallel scRNA-seq and scATAC-seq data of PBMC](https://htmlpreview.github.io/?https://github.com/zhanglhbioinfor/DIRECT-NET/blob/main/tutorial/demo_DIRECTNET_PBMC.html)

- [Demo of DIRECT-NET on parallel scRNA-seq and scATAC-seq data of A549](https://htmlpreview.github.io/?https://github.com/zhanglhbioinfor/DIRECT-NET/blob/main/tutorial/demo_DIRECTNET_A549.html)

- [Demo of DIRECT-NET on scATAC-seq data of GM12878](https://htmlpreview.github.io/?https://github.com/zhanglhbioinfor/DIRECT-NET/blob/main/tutorial/demo_DIRECTNET_GM12878.html)

- [Demo of DIRECT-NET on scATAC-seq data of human Brain](https://htmlpreview.github.io/?https://github.com/zhanglhbioinfor/DIRECT-NET/blob/main/tutorial/demo_DIRECTNET_Brain.html)

## Help
If you have any problems, comments or suggestions, please contact us at Lihua Zhang (zhanglh@whu.edu.cn).


