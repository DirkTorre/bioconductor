# RNA-seq analysis recap

This document contains background information on the experiment and used methods.

VERDER MET: Data pre-processing

## links:

* [RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR](https://f1000research.com/articles/5-1408/v3)
* [RNAseq123](https://bioconductor.org/packages/release/workflows/html/RNAseq123.html)
* [A guide to creating design matrices for gene expression experiments](https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html)
* []()


## background

### experiment

The dataset comes from mouse mammary gland tissue.
The workflow uses this dataset to illustrate:
- differences between biological groups (e.g., different cell populations or conditions),
- how to detect differentially expressed genes,
- how to interpret results at the pathway level.
The biological details are not the focus—the dataset is simply a realistic example that allows the authors to demonstrate the analysis pipeline.

### code

#### experiment walkthrough

short summary:

1. acquire gene level counts.
2. preprocessing of the counts.
3. exploratory data analysis.
4. differential expression testing
5. pathway analysis

More detailed summary:

- Data packaging
  + Read count matrix
  + Organising sample metadata: Attach metadata to the samples. So we can make a design matrix for the experiment. The design matrix is needed to model batch effects and group differences.
  + Add gene annotations: Making gene id's human readable.
- Data pre-processing
  + Transform raw counts: samples have different sequencing depths, so we can't compare raw count's fairly. we use log-cpm.
  + Removing genes that are lowly expressed: because they are not statistically meaningful.
  + Normalising gene expression distributions: correct for over representation of sequences and global differences between samples. we use TMM.
  + Unsupervised clustering of samples: to explore sample relationships and detect batch effects etc. If samples do not cluster by biological group, differential expression will be weak. MDS reveals major sources of variation (e.g., cell type vs. batch). Helps decide which factors to include in the design matrix.
- Differential expression analysis
  + Creating a design matrix and contrasts: Specify which comparisons (contrasts) you want to test.
  + Removing heteroscedascity from count data: linear models assume constant variance. But the data has variance that depended on the mean. We use voom to correct this, which leads to more accurate p-values and fold change estimates.
  + Fitting linear models for comparisons of interest: estimate expression diffenreces between groups. linear modeling is used to stabilize variance estimates.
  + Examining the number of DE genes: Assess if the comparison a meaningfull result. Too few DE genes indicates poor clustering, low power or technical results. too many indicate batch effects or normaization problems.
  + Examining individual DE genes from top to bottom: validate the top hit's and ensure the results make sense.
  + Useful graphical representations of differential expression results: because plot's are easier to interpret.
- Gene set testing with camera: Test whether predefined gene sets (pathways, GO terms) show coordinated changes. camera accounts for inter‑gene correlation, improving statistical validity.


AI summary:

Below is a clear, accurate summary of the **research study linked to the workflow** and an explanation of **how the RNA‑seq analysis fits into that research**. All factual statements are grounded in the retrieved sources.   [F1000Research](https://f1000research.com/articles/5-1408)  [ResearchGate](https://www.researchgate.net/publication/304065881_RNA-seq_analysis_is_easy_as_1-2-3_with_limma_Glimma_and_edgeR/fulltext/5d2f71caa6fdcc2462e86ebb/RNA-seq-analysis-is-easy-as-1-2-3-with-limma-Glimma-and-edgeR.pdf)

---

#####  🧬 **Summary of the Research Behind the Workflow**

The article *“RNA‑seq analysis is easy as 1‑2‑3 with limma, Glimma and edgeR”* is **not a biological research study** in the traditional sense. Instead, it is a **workflow / software tool article** published in F1000Research. Its purpose is to **demonstrate a complete, practical RNA‑seq analysis pipeline** using Bioconductor packages—specifically **edgeR**, **limma (with voom)**, and **Glimma**.

To make the workflow concrete, the authors analyze a **real RNA‑seq dataset from the mouse mammary gland**. This dataset provides a biological context, but the *goal of the article is methodological*, not to report new biological discoveries.

######  What the article aims to do  
- Show how to import, organize, filter, normalize, and explore RNA‑seq count data.  
- Demonstrate differential expression analysis using **limma‑voom**.  
- Demonstrate interactive visualization using **Glimma**.  
- Show how to perform pathway‑level analysis using **camera**.  
- Provide a reproducible, end‑to‑end workflow for beginners and experienced users.

###### Why this matters  
The authors emphasize that Bioconductor makes RNA‑seq analysis **efficient, reproducible, and accessible**, and they illustrate this by walking through a complete analysis of a real dataset.   [F1000Research](https://f1000research.com/articles/5-1408)

---

###### 🐭 **What is the underlying biological dataset?**

The dataset comes from **mouse mammary gland tissue**.  
The workflow uses this dataset to illustrate:

- differences between biological groups (e.g., different cell populations or conditions),  
- how to detect differentially expressed genes,  
- how to interpret results at the pathway level.

The biological details are not the focus—the dataset is simply a **realistic example** that allows the authors to demonstrate the analysis pipeline.   [ResearchGate](https://www.researchgate.net/publication/304065881_RNA-seq_analysis_is_easy_as_1-2-3_with_limma_Glimma_and_edgeR/fulltext/5d2f71caa6fdcc2462e86ebb/RNA-seq-analysis-is-easy-as-1-2-3-with-limma-Glimma-and-edgeR.pdf)

---

###### 🔬 **How the RNA‑seq Analysis Fits the Research**

Although the article is methodological, the RNA‑seq analysis steps are chosen to reflect **best practices** for analyzing real biological data. Here’s how each major stage fits into the purpose of the article:

---

####### 1. **Data Packaging**  
The authors show how to structure the dataset into Bioconductor objects (e.g., DGEList) so that:

- counts, sample metadata, and gene annotations stay synchronized,  
- downstream functions can operate cleanly.

This mirrors what any real RNA‑seq project must do before analysis.

---

####### 2. **Data Pre‑processing**  
The workflow demonstrates essential steps required for valid statistical inference:

- **Transformations** (CPM, log‑CPM)  
- **Filtering lowly expressed genes**  
- **Normalization** (TMM)  
- **Unsupervised clustering** (MDS)

These steps ensure that the dataset is clean, comparable, and biologically interpretable before modeling.

---

####### 3. **Differential Expression Analysis**  
The authors use **limma‑voom**, which models the mean–variance relationship in RNA‑seq data and fits linear models.

This part of the workflow shows:

- how to build a design matrix,  
- how to specify contrasts,  
- how to identify differentially expressed genes,  
- how to visualize and interpret results.

This is the core of most RNA‑seq studies, and the article uses the mouse dataset to demonstrate it in a realistic setting.

---

####### 4. **Gene Set Testing (camera)**  
The workflow concludes with pathway‑level analysis using **camera**, which accounts for inter‑gene correlation.

This step shows how to move from individual genes to **biological processes**, which is essential for interpreting RNA‑seq results in a meaningful biological context.

---


### research data

data comes from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63310

Title: Transcriptome profiling of purified mouse mammary stem, progenitor and mature cell populations

Purpose: The aim of this study is to determine the absolute and relative expression levels of mRNA transcripts across multiple flow cytometrically sorted epithelial cell types including freshly isolated CD24+CD29hi mammary stem cell-enriched basal cells (MaSC/basal), CD24+CD29loCD61+ luminal progenitor-enriched (LP) and the CD24+CD29loCD61- mature luminal-enriched (ML) cell populations. Additionally, a comparison between these primary cell types and cultured MaSC/Basal-derived mammosphere cells (mammosphere) and the CommaD-βGeo (CommaDβ) cell line was performed.

Methods: Total RNA was extracted and purified from sorted luminal or basal populations from the mammary glands of female virgin 8- to 10-week-old FVB/N mice (3 independent samples per population), MaSC/Basal cells cultured for 1 week under mammosphere conditions and CommaDβ cells grown under maintenance conditions (Deugnier et al. 2006). Total RNA (100 ng) was used to generate sequencing libraries for whole transcriptome analysis following Illumina’s TruSeq RNA v2 sample preparation protocol. Completed libraries were sequenced on HiSeq 2000 with TruSeq SBS Kit v3- HS reagents (Illumina) as 100 bp single-end reads at the Australian Genome Research Facility (AGRF), Melbourne. Approximately 30 million 100 bp single-end reads were obtained for each sample. Reads were aligned to the mouse reference genome mm10 and mapped to known genomic features at the gene level using the Rsubread package (version 1.14.1) (Liao et al. 2013). Single reads were then summarized into gene-level counts using FeatureCounts (Liao et al. 2014).
  	
Overall design 	Total RNA was extracted and purified from sorted luminal or basal populations from the mammary glands of female virgin 8- to 10-week-old FVB/N mice (3 independent samples for population), MaSC/Basal cells cultured for 1 week under mammsophere conditions and CommaDβ cells grown under maintenance conditions (Deugnier et al. 2006) and their transcriptomes analysed by RNA-Seq.


Here is an AI summary in plain english, it's better than my initial summary.

Here’s a **plain‑English explanation** of the research behind **GSE63310**, based directly on the information retrieved from the web   [Github](https://github.com/UCSF-DSOS/bulk_RNA_seq/blob/master/notebooks/2_downloadRNAseq.Rmd)  [Enrichr](https://maayanlab.cloud/datasets2tools/landing/dataset/GSE63310).

---

#### 🧬 **GSE63310 — Explained in Plain English**

##### 🌟 What the researchers wanted to understand  
The study looked at **different types of cells in the mouse mammary gland** (the tissue that produces milk). The mammary gland contains several kinds of epithelial cells, each with a different role:

- **Mammary stem cells (MaSC/basal)** — the “mother” cells that can create all other mammary cell types.  
- **Luminal progenitor cells (LP)** — cells partway along the developmental path; they can still divide and specialize.  
- **Mature luminal cells (ML)** — fully specialized cells that line the ducts and help produce milk.

The researchers wanted to **compare the gene activity (which genes are turned on/off)** in these different cell types to understand:

- how they differ,  
- how they relate to each other,  
- and what makes stem cells unique compared to more mature cells.

They also compared these primary cells to:

- **mammosphere cultures** (stem‑cell–enriched cells grown in a dish), and  
- **Comma‑D cells**, a commonly used mammary epithelial cell line.

---

#### 🔬 **What they did**

###### 1. **Collected cells from mouse mammary glands**  
They used flow cytometry (cell sorting) to isolate very pure populations of:

- stem‑cell–enriched basal cells,  
- luminal progenitors,  
- mature luminal cells.

###### 2. **Extracted RNA and performed RNA‑seq**  
For each cell type, they:

- extracted RNA,  
- built sequencing libraries,  
- sequenced them on an Illumina HiSeq 2000 machine,  
- generated ~30 million reads per sample.

###### 3. **Mapped the reads to the mouse genome**  
They aligned the reads to the mouse genome (mm10) and counted how many reads mapped to each gene.

This produced the **gene‑level count data** used in the workflow you’re studying.

---

#### 🧠 **What the research is really about**

This dataset helps answer questions like:

- Which genes define mammary stem cells?  
- How do progenitor cells differ from mature cells?  
- Are cultured cells (mammospheres, Comma‑D) similar to real tissue cells?  
- What pathways are active in each cell type?

It’s a **cell‑type comparison study** aimed at understanding mammary gland biology and stem cell differentiation.

---

#### 🔗 **How the RNA‑seq analysis (from the F1000 workflow) fits into this research**

The workflow you’re studying (limma/voom + edgeR + Glimma) uses **this exact dataset** as its example. Here’s how each analysis step connects to the biological study:

##### ✔ Data packaging  
The dataset contains multiple cell types. Organizing sample metadata is essential to compare them correctly.

##### ✔ Pre‑processing  
Filtering, normalization, and clustering help ensure:

- low‑quality genes are removed,  
- samples behave as expected (e.g., stem cells cluster together),  
- differences reflect biology, not sequencing artifacts.

##### ✔ Differential expression  
This is the heart of the biological question:

- Which genes are more active in stem cells vs. progenitors?  
- Which genes define mature luminal cells?  
- How do cultured cells differ from real tissue?

The design matrix and contrasts directly encode these comparisons.

##### ✔ Gene set testing (camera)  
Instead of looking at single genes, this step asks:

- Which **biological pathways** differ between cell types?  
- Are stem cells enriched for self‑renewal pathways?  
- Are luminal cells enriched for milk‑production pathways?

This is crucial for interpreting the biology.

---

#### 🎯 **In short**

**GSE63310** is a study that compares gene expression across different mammary epithelial cell types in mice to understand how stem cells differ from more mature cells.

The RNA‑seq workflow you’re studying uses this dataset as a **real‑world example** to demonstrate:

- how to clean the data,  
- how to compare groups,  
- how to find differentially expressed genes,  
- and how to interpret biological pathways.

It’s a perfect teaching dataset because it has clear biological differences and multiple well‑defined cell types.

---







## definitions

### R packages

* **limma**: Data analysis, linear models and differential expression for omics data. 
* **Glimma**: This package produces interactive visualizations for RNA-seq data analysis, utilizing output from limma, edgeR, or DESeq2. It produces interactive htmlwidgets versions of popular RNA-seq analysis plots to enhance the exploration of analysis results by overlaying interactive features. The plots can be viewed in a web browser or embedded in notebook documents. 
* **edgeR**: Differential expression analysis of sequence count data. Implements a range of statistical methodology based on the negative binomial distributions, including empirical Bayes estimation, exact tests, generalized linear models, quasi-likelihood, and gene set enrichment. Can perform differential analyses of any type of omics data that produces read counts, including RNA-seq, ChIP-seq, ATAC-seq, Bisulfite-seq, SAGE, CAGE, metabolomics, or proteomics spectral counts. RNA-seq analyses can be conducted at the gene or isoform level, and tests can be conducted for differential exon or transcript usage. 

### life science glossery

mammosphere: sphere of cells that is formed from a single mammalian stem cell.

CD: antibodies.

luminal progenitor-enriched: a cell that can differentiate into different cell types cells from the inside of an organ that stores or transports stuff.

cytromery:
https://en.wikipedia.org/wiki/Cytometry
Cytometry is the measurement of number and characteristics of cells. Variables that can be measured by cytometric methods include cell size, cell count, cell morphology (shape and structure), cell cycle phase, DNA content, and the existence or absence of specific proteins on the cell surface or in the cytoplasm.[

CD24:
https://pubmed.ncbi.nlm.nih.gov/36013184/
Cluster of differentiation 24 (CD24) is a small, highly glycosylated cell adhesion protein that is normally expressed by immune as well as epithelial, neural, and muscle cells. Tumor CD24 expression has been linked with alterations in several oncogenic signaling pathways. In addition, the CD24/Siglec-10 interaction has been implicated in tumor immune evasion, inhibiting macrophage-mediated phagocytosis as well as natural killer (NK) cell cytotoxicity.


CD29: 
https://en.wikipedia.org/wiki/Integrin_beta_1
Integrin beta-1 (ITGB1), also known as CD29, is a cell surface receptor that in humans is encoded by the ITGB1 gene.[5] This integrin associates with integrin alpha 1 and integrin alpha 2 to form integrin complexes which function as collagen receptors. It also forms dimers with integrin alpha 3 to form integrin receptors for netrin 1 and reelin. These and other integrin beta 1 complexes have been historically known as very late activation (VLA) antigens.

Integrin beta 1 is expressed as at least four different isoforms. In cardiac muscle and skeletal muscle, the integrin beta-1D isoform is specifically expressed, and localizes to costameres, where it aids in the lateral force transmission from the Z-discs to the extracellular matrix. Abnormal levels of integrin beta-1D have been found in limb girdle muscular dystrophy and polyneuropathy.
  
CD61:
https://en.wikipedia.org/wiki/Integrin_beta_3
Integrin beta-3 (β3) or CD61 is a protein that in humans is encoded by the ITGB3 gene.[5] CD61 is a cluster of differentiation found on thrombocytes.[6]
The ITGB3 protein product is the integrin beta chain beta 3. Integrins are integral cell-surface proteins composed of an alpha chain and a beta chain. A given chain may combine with multiple partners resulting in different integrins. Integrin beta 3 is found along with the alpha IIb chain in platelets. Integrins are known to participate in cell adhesion as well as cell-surface-mediated signaling.[7]
Defectively expressed β3 integrin subunit has been correlated with presence of endometriosis, and has been suggested as a putative marker of this condition.[8]
CD61 has been shown to interact with PTK2,[9][10] ITGB3BP,[11][12] TLN1[13][14] and CIB1.[15]

### stats glossery

heteroscedascity: Heteroscedasticity refers to a condition in statistics where the variability of a dependent variable is not constant across all levels of an independent variable. So when there is difference in variance.
https://en.wikipedia.org/wiki/Homoscedasticity_and_heteroscedasticity
homoscedastic: relation between x and y is pretty much linear, hetrosecdastic: looks linear at the start, but it diverges.


Negative Binomial distribution: The negative binomial distribution is a probability distribution that models the number of failures in a series of independent trials before a specified number of successes occurs. It is often used in scenarios where the variance of the data exceeds the mean, making it a useful alternative to the Poisson distribution.

# other notes

[Camera: a competitive gene set test accounting for inter-gene correlation](https://pmc.ncbi.nlm.nih.gov/articles/PMC3458527/)
"Competitive gene set tests are commonly used in molecular pathway analysis to test for enrichment of a particular gene annotation category amongst the differential expression results from a microarray experiment."
I'm not fully sure if this is the same as [GSEA](https://en.wikipedia.org/wiki/Gene_set_enrichment_analysis).

