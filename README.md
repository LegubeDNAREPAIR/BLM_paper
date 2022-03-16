# A POLD3/BLM dependent pathway handles DSBs in transcribed chromatin upon excessive RNA:DNA hybrids accumulation.

Cohen S, Guenolé A, Lazar I, Marnef A, Clouaire T, Vernekar DV, Puget N, Rocher V, Arnould C, Aguirrebengoa M, Genais M, Firmin N, Shamanna RA, Mourad R, Bohr VA, Borde V and Legube G.


## Overview

Scripts used to generate all figures from our article and some raw data.

## Data accession

Data will be available at [??](https://www.ebi.ac.uk/arrayexpress/experiments/??/)

## System requirements

### Alignment

see https://github.com/LegubeDNAREPAIR/HistoneMapping

### Figure generation

the scripts were written with `R`, and need some packages : 

* `library(rtracklayer)`.
* `library(BSgenome.Hsapiens.UCSC.hg19)`.
* `library(reshape2)`.
* `library(tidyverse)`.
* `library(Homo.sapiens)`.

Some data (profile) were generated using `deeptools`, as you can see [here](script/snakefile)

| Scripts           | Description                                                                                       | Figures               |
|-------------------|---------------------------------------------------------------------------------------------------|-----------------------|
| src/Functions.R   | Usefull functions for some R scripts                                                              |                       |
| BLM_BLESS_CAT.R   | Generate DSB Categories for BLM values                                                            |                       |
| profile_boxplot.R | Some code to generate average profile / boxplots with genomic position and coverage files(bigWig) | Figures 3D,3E,S4D     |
| Edu_plots.R       | Manipulate and plot Edu-seq data                                                                  | Figures 4D 4E S1G S6B |
| FigG4.R           | Manipulate published G4 ChIP-seq data (GSE99205)                                                  | Figure S1E            |
| Fig1B.R           |                                                                                                   | Figure 1B             |
| Fig1C_S1C.R       |                                                                                                   | Figures 1C, S1C       |
| FIG1D_S1F.R       |                                                                                                   | Figures 1D, S1F       |
| Fig2B_3B.R        |                                                                                                   | Figures 2B, 3B        |
| Fig2D.R           |                                                                                                   | Figure 2D             |
| FIGS1A.R          |                                                                                                   | Figure S1A            |
| FigS1D.R          |                                                                                                   | Figure S1D            |