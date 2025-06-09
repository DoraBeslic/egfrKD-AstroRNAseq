# 🧬 EGFR-KD Astrocyte RNA-Seq Analysis (Replication of Fig. 1J–K)

This repository contains a reproducible pipeline and supporting files for replicating **Figure 1J and 1K** from the paper:

> *"Psychedelic control of neuroimmune interactions governing fear"*  
> *(Analysis of amygdala astrocyte responses during chronic stress by EGFR, [NCBI GEO DataSets: GSE262981])*

We focus on replicating the **differential gene expression** and **gene set enrichment** analyses of **astrocytes with EGFR knockdown (EGFR-KD)**, which demonstrated that **EGFR signaling in the amygdala limits stress-induced fear behavior**.

---

## 📄 Summary

- **Organism:** *Mus musculus* (GRCm39)
- **Tissue:** Amygdala
- **Cell type:** Astrocytes
- **Perturbation:** Transduced with sgEgfr (n = 3) or sgRosa26 (n = 4; control)
- **Analysis Goal:** Identify EGFR-regulated genes and pathways relevant to neuroimmune signaling during stress-induced fear behavior
- **Figures Replicated:**  
  - **Figure 1J:** Volcano plot of differentially expressed genes  
  - **Figure 1K:** Gene set enrichment plots

### Key Findings from Analysis

- Bulk RNA-seq analysis of EGFR-KD amygdala astrocytes revealed **increased activation of inflammatory pathways**, consistent with EGFR’s known anti-inflammatory role in astrocytes.  
- Expression of **Nptx1** and **Fos** was increased; these genes are linked to stabilization of fear-memory formation by amygdala astrocytes.  
- Loss of EGFR function also led to elevated transcription of genes related to **receptor protein tyrosine phosphatases**, which are downstream targets of **NF-κB signaling** and involved in heterotypic cell–cell interactions.  

These molecular changes underline the role of astrocyte EGFR in limiting pro-inflammatory signaling and fear behavior in the amygdala.

---

## 🧪 Biological Context

Chronic stress induces recruitment of meningeal monocytes and pro-inflammatory signaling in the brain. EGFR expression in amygdala astrocytes suppresses this response and limits fear behavior. In this analysis, we replicate findings that:
- EGFR-KD upregulates pro-inflammatory gene expression.
- These changes are associated with stress-induced behavior.
- Psychedelics may reverse these molecular and behavioral effects.

---

## ⚙️ Workflow Overview

This project was implemented using **Nextflow** with execution on **AWS Batch** for scalability and reproducibility.

| Step                     | Tool            | Version    |
|--------------------------|-----------------|------------|
| Data download            | `sra-tools`     | 3.2.1
| Quality control          | `fastp`, `FastQC` | 0.23.4, 0.12.1 |
| Multi-sample QC report   | `MultiQC`       | 1.14       |
| Alignment                | `STAR`          | 2.7.10b     |
| Gene quantification      | `featureCounts` | 2.0.1      |
| Gene filtering           | Custom (avg count > 0.5) | –  |
| Differential expression  | `edgeR`         | 4.6.2     |

---

## 🚀 How to Reproduce

1. **Configure AWS Batch** credentials and compute environment  
2. **Run pipeline**:  
   ```bash
   nextflow run main.nf -profile awsbatch
   
## Citations
Chung, E.N., Lee, J., Polonio, C.M. et al. Psychedelic control of neuroimmune interactions governing fear. Nature 641, 1276–1286 (2025). https://doi-org.ezproxy.neu.edu/10.1038/s41586-025-08880-9

Guangchuang Yu, Li-Gen Wang, Yanyan Han and Qing-Yu He. clusterProfiler: an R package for
comparing biological themes among gene clusters. OMICS: A Journal of Integrative Biology.
2012, 16(5):284-287

T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo,
and G Yu. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The
Innovation. 2021, 2(3):100141



