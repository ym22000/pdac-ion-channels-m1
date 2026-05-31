# Identification and bioinformatic analysis of ion channels that may modulate immunity and extracellular matrix in PDAC stroma

Document prepared for the internship project carried out by **MISSOUM Youssra**, Master 1 student in **Bioinformatics and Computational Biology**, Universite Cote d'Azur, Nice.
The work was performed at the **Institut de Biologie Valrose (IBV)** in the team **Regulation of ion channel in cancer** led by Olivier Soriani.
Team page: [http://ibv.unice.fr/research-team/soriani/](http://ibv.unice.fr/research-team/soriani/)

---

## Executive summary

The main objective of the project was:

> **to identify and characterize, by bioinformatics, ion channels and especially KCN genes that could modulate immunity and/or extracellular matrix in pancreatic ductal adenocarcinoma through stromal territories**

The project did not identify one universal "stromal immune channel". Instead, it defined a structured KCN landscape across several complementary layers:

1. **single-cell RNA-seq** to define a stromal KCN panel and place genes in fibroblast states
2. **multicohort bulk survival** to test whether the selected KCN genes also have clinical relevance
3. **local bulk perturbation datasets** to reconnect some candidates to experimental contexts from the laboratory
4. **spatial transcriptomics** to determine where these KCN genes actually fall in the tissue, with which programs, and in which niches

The strongest integrated message is that the KCN panel separates into at least three major biological worlds:

1. **tumor ductal / epithelial interface**
   - `KCNK1`, `KCNN4`, `KCNK6`
2. **stromal contractile / ECM / myofibroblastic world**
   - `KCNMA1`, `KCNMB1`, `KCNE4`, `KCND2`
3. **resident / perivascular / immune-border world**
   - `KCNJ8`, `KCNK3`, `KCND3`, `KCNAB2`, `KCNN3`

This matters for the internship question because stromal immunity in PDAC is often modulated indirectly through:

- fibroblast contractility
- ECM deposition and remodeling
- tumor-stroma interfaces
- vascular and perivascular territories
- inflammatory and complement-rich border zones

In other words, the most relevant ion channels in this project are not only "immune genes". They are also candidate regulators or markers of stromal niches that organize matrix, signaling, exclusion, and local tissue architecture.

---

## 1. Initial question and biological rationale

### 1.1 Core question

The project started from a broader question than a simple differential expression analysis:

- which ion channels, especially `KCN` genes, are associated with the stromal compartment of PDAC?
- are these channels linked to specific stromal states?
- are they connected to programs compatible with modulation of:
  - immunity
  - extracellular matrix
  - fibroblast activation
  - tumor-stroma interfaces?

### 1.2 Biological hypothesis

The project is based on a biologically plausible hypothesis:

- ion channels regulate membrane potential, cell volume, and ionic fluxes such as `K+`, `Ca2+`, `Na+`, and `Cl-`
- these properties can influence:
  - fibroblast contractility
  - migration and adhesion
  - cytokine secretion
  - ECM production
  - spatial organization of stromal territories
- in PDAC, such effects can indirectly modulate immune cell distribution and function

The goal was therefore to transform a plausible mechanistic idea into a structured bioinformatic prioritization.

---

## 2. Datasets used

## 2.1 Global overview

| Analytical block | Dataset | Source | Role in the project |
|---|---|---|---|
| Stromal scRNA-seq global atlas | `data/Stroma_Subset2021.rds` | Oh et al., *Nature Communications* 2023 | global stromal exploration |
| Fibroblast-focused atlas | `data/Stroma_Subset2021_fibroblast_focus.rds` | curated subset from the same atlas | KCN prioritization and subtype interpretation |
| External fibroblast validation | USER local inputs | USER / humanpdac | independent validation of stromal state framework |
| Spatial transcriptomics | `data/PDAC_Updated_ST.rds` | Khaliq et al., *Nature Genetics* 2024 | tissue localization and niche interpretation |
| Local bulk perturbation 1 | `code/bulk/BKCA_shBKCA` | IBV laboratory | `BKCa / KCNMA1` context |
| Local bulk perturbation 2 | `code/bulk/SIGMAR1_shSIGMAR1` | IBV laboratory | `SIGMAR1` context |
| Bulk survival cohorts | `code/bulk_survival_analysis/` | GSE183795, GSE62452, GSE28735, Bailey 2016 | clinical relevance |

## 2.2 Why this multimodal combination matters

Each layer answers a different question:

- **single-cell RNA-seq** tells us which fibroblast states are involved
- **bulk survival** tells us whether some genes have clinical signal at patient level
- **local bulk perturbations** reconnect the bioinformatic signals to experimental laboratory axes
- **spatial transcriptomics** tells us where these genes fall in tissue and with which programs they spatially co-organize

The strength of the project comes from the convergence of these layers rather than from one isolated test.

---

## 3. Global analytical chronology

The project was built in the following order:

1. stromal single-cell exploration
2. fibroblast-focused KCN prioritization
3. state and subtype interpretation in scRNA-seq
4. `myCAF` / `iCAF` program-oriented interpretation
5. external validation in USER
6. multicohort bulk survival on the final KCN panel
7. local bulk perturbation analyses from the laboratory
8. spatial transcriptomics:
   - maps
   - ecotypes
   - niche tests
   - Hotspot
   - SPARK-X
   - local correlations
   - territory GSEA
   - neighborhood analysis
   - pathway overlays
   - normal pancreas comparison
9. integrated synthesis

This document follows this chronology.

---

## 4. Single-cell RNA-seq: stromal and fibroblast KCN prioritization

## 4.1 Why start with single-cell RNA-seq?

This was the right entry point because the central biological question was stromal. A bulk starting point would have mixed:

- tumor epithelial signal
- stromal signal
- immune signal
- vascular signal

The single-cell atlas allows us to:

- isolate the stromal world
- focus on fibroblast states
- define a biologically coherent KCN panel before moving to spatial and clinical layers

## 4.2 Global stromal atlas used first

Source object:

- `data/Stroma_Subset2021.rds`

Documented global metrics:

- `21469` cells
- `16972` genes
- `40` samples
- `3` datasets
- `8` `CellType2` labels
- `77.32%` `Activated`
- `22.68%` `Normal`

Main datasets represented in the stromal object:

- `Powers`
- `Peng`
- `Qadir`

This first object was used to understand the overall stromal landscape before narrowing the analysis to fibroblast-focused questions.

## 4.3 Fibroblast-focused object used for KCN prioritization

Source object:

- `data/Stroma_Subset2021_fibroblast_focus.rds`

The main fibroblast states considered in the project were:

- `qPSC`
- `smPSC`
- `myCAF`
- `csCAF`
- `iCAF`

Documented subtype counts in the AUCell / fibroblast-focused analyses:

- `qPSC = 5743`
- `smPSC = 5521`
- `myCAF = 4809`
- `csCAF = 4242`
- `iCAF = 338`

This distribution already contains an important caveat:

- the `iCAF` compartment is much smaller than the others
- therefore all `iCAF` conclusions must be read with more caution

## 4.4 Activated vs Normal stromal states

The project used the Moffitt stromal framework to define a simple global polarity:

`StromaScore = Moffitt.F5_ActivatedStroma.top25 - Moffitt.F13_NormalStroma.top25`

Rule:

- `Activated` if `F5 > F13`
- `Normal` otherwise

This step did not yet select KCN genes. It established a biologically interpretable backbone:

- a global activated stromal world
- a global normal-like stromal world

This was useful because later KCN genes could be interpreted not only as subtype markers but also as genes linked to a broad stromal polarity.

## 4.5 KCN search space and prefiltering

The initial KCN search space was defined from gene symbols matching the `^KCN` pattern.

Then a permissive but biologically motivated prefilter was applied:

- detected in at least `1%` of cells in the fibroblast object
- detected at at least `1%` in at least `2` stromal subtypes

Why use a low threshold?

- many ion channels are biologically meaningful despite sparse detection
- too strong a threshold would have removed niche-specific candidates

Why not keep every KCN without filtering?

- to reduce instability from genes supported by only a handful of cells
- to stabilize downstream mixed models and pseudobulk tests

## 4.6 Why four DEG strategies were used

The central methodological choice of the single-cell block was:

> **do not trust a single DEG framework for sparse KCN genes**

This is particularly important for ion channels because they are often:

- weakly detected
- zero-inflated
- heterogeneous across patients
- informative either through detection rate or through intensity among positives

Therefore, the final KCN panel was built by convergence across four different analytical strategies.

## 4.7 Test 1: pseudobulk DESeq2 Activated vs Normal

### Goal

Identify KCN genes distinguishing the global `Activated` stromal state from the `Normal` stromal state.

### Unit of analysis

Pseudobulk aggregated by:

- `Patient`
- `StromaState`

Retention rules:

- `n_cells >= 20` per pseudo-sample
- keep only patients represented in both states

### DESeq2 design

`~ Patient + n_cells + condition`

where:

- `condition = Activated / Normal`

### DESeq2 settings

- genes kept if `rowSums(counts) >= 10`
- `estimateSizeFactors(type = "poscounts")`
- `DESeq(sfType = "poscounts", fitType = "local")`

Tested coefficient:

- `condition_Activated_vs_Normal`

Decision threshold:

- `BH FDR < 0.05`

### Why this model is useful

This pseudobulk framework:

- reduces cell-level pseudoreplication
- treats the patient as the biological unit
- gives a robust global stromal contrast

### Main result pattern

Representative genes enriched in `Activated`:

- `KCNK6`
- `KCNT2`
- `KCND3`
- `KCNN3`
- `KCNJ8`

Representative genes enriched in `Normal`:

- `KCNAB1`
- `KCNA5`
- `KCNMA1`
- `KCNK17`
- `KCNC4`
- `KCNQ4`
- `KCNS3`
- `KCNE4`
- `KCNB1`

Documented number of significant KCN genes in this test:

- `14`

### Limits

- loss of cell-level resolution
- a rare but real subcluster signal can be diluted
- this is a state-level rather than niche-level analysis

## 4.8 Test 2: pseudobulk DESeq2 subtype vs rest

### Goal

Identify KCN genes enriched in a specific fibroblast subtype.

### Logic

For each subtype, build a contrast:

- `target` vs `rest`

### Design

`~ Patient + n_cells + condition`

where:

- `condition = target / rest`

### Decision thresholds

- `BH FDR < 0.05`
- `|log2FC| >= 0.25`

### Why this test matters

A gene can be very informative for one fibroblast state without being a strong global `Activated vs Normal` marker.

This test was therefore necessary to capture subtype-specific KCN patterns.

### Global result pattern

Approximate documented counts of significant KCN genes by subtype:

- `qPSC`: `8`
- `smPSC`: `9`
- `myCAF`: `7`
- `csCAF`: `1`
- `iCAF`: almost no robust signal

Examples:

- `qPSC`: `KCNJ8`, `KCNK17`, `KCNE4`, `KCNK6`
- `smPSC`: `KCNA5`, `KCNAB1`, `KCNMB1`, `KCNMA1`
- `myCAF`: `KCNMA1`, `KCNMB1`, `KCND3`, `KCNK3`
- `csCAF`: mainly `KCNK17`

### Limits

- the `rest` group is biologically heterogeneous
- close subtypes can partially mask each other

## 4.9 Test 3: hurdle mixed model

### Goal

Model sparse KCN genes more appropriately by separating:

- **detection**
- **positive-cell expression intensity**

### Detection model

`glmer(detected ~ is_target + log10_nCount_RNA + Dataset + (1 | Patient), family = binomial)`

### Positive expression model

`lmer(expr_pos ~ is_target + log10_nCount_RNA + Dataset + (1 | Patient))`

### Eligibility rules

Detection part:

- at least `20` positive cells overall
- at least `3` positive patients in target and rest
- both detection classes present

Positive-expression part:

- at least `20` positive cells overall
- at least `10` positive cells in target
- at least `10` positive cells in rest
- at least `3` positive patients in target and rest
- more than one distinct positive value

### Combination

- combine p-values with Fisher's method
- support classified as:
  - `both_parts`
  - `detection_only`
  - `positive_only`
  - `none`

### Why this model is especially relevant for KCN genes

Some channels are informative because:

- they are detected in more cells
- or because among positive cells they are stronger

A hurdle framework captures both types of biology.

### Documented results

Number of KCN genes with support by subtype:

- `qPSC`: `18`
- `smPSC`: `17`
- `myCAF`: `16`
- `csCAF`: `16`
- `iCAF`: `10`

This was one of the most helpful tests for sparse and niche-like KCN genes.

### Limits

- more complex model
- less intuitive than a plain log2FC
- still relies on good subgroup support

## 4.10 Test 4: MAST

### Goal

Add a standard single-cell DEG method adapted to zero-inflated expression.

### Main settings

Using `FindMarkers(..., test.use = "MAST", latent.vars = "nCount_RNA")`

Decision threshold:

- `padj_BH < 0.05`

Direction interpreted with:

- `avg_log2FC >= 0.25`

### Why MAST was useful

MAST provides a well-established single-cell DEG perspective complementary to pseudobulk and mixed models.

### Documented results

Approximate counts of significant KCN genes by subtype:

- `qPSC`: `19`
- `smPSC`: `18`
- `myCAF`: `18`
- `csCAF`: `21`
- `iCAF`: `10`

Examples:

- `qPSC`: `KCNJ8`, `KCNK17`, `KCNE4`
- `smPSC`: `KCNA5`, `KCNAB1`, `KCNMB1`, `KCNMA1`
- `myCAF`: `KCNK17`, `KCNJ8`, `KCNA5`, `KCNE4`, `KCNMB1`, `KCND2`
- `iCAF`: weaker but `KCNJ8`, `KCNS3`, `KCNK1`, `KCNK6` were among the most visible

### Limits

- still sensitive to single-cell structure
- patient is less explicitly handled than in pseudobulk

## 4.11 Final union panel: 25 KCN genes

### Inclusion rule

A gene entered the final panel if it was retained by at least one of the four DEG strategies:

- Test 1: significant by BH
- Test 2: significant by BH and `|log2FC| >= 0.25`
- Test 3: `support_pattern != none`
- Test 4: significant by BH

### Final panel

The final union contains **25 KCN genes**:

`KCNA5`, `KCNA7`, `KCNAB1`, `KCNAB2`, `KCNB1`, `KCNC4`, `KCND2`, `KCND3`, `KCNE4`, `KCNJ15`, `KCNJ16`, `KCNJ8`, `KCNK1`, `KCNK15`, `KCNK17`, `KCNK3`, `KCNK6`, `KCNMA1`, `KCNMB1`, `KCNMB4`, `KCNN3`, `KCNN4`, `KCNQ4`, `KCNS3`, `KCNT2`

### Inter-method robustness table

| Number of supporting methods | Genes |
|---|---|
| `4/4` | `KCNA5`, `KCNAB1`, `KCNC4`, `KCND3`, `KCNE4`, `KCNJ8`, `KCNK17`, `KCNK6`, `KCNMA1` |
| `3/4` | `KCNK3`, `KCNMB1`, `KCNN3`, `KCNS3`, `KCNT2` |
| `2/4` | `KCNA7`, `KCND2`, `KCNJ15`, `KCNK1`, `KCNK15`, `KCNMB4`, `KCNN4`, `KCNQ4` |
| `1/4` | `KCNAB2`, `KCNB1`, `KCNJ16` |

### How to read this panel

- `4/4`: highly robust candidates
- `3/4`: strong but more method-dependent
- `2/4`: biologically interesting, moderate statistical robustness
- `1/4`: kept by biological openness, should be interpreted more carefully

This final union was crucial because it defined the shared candidate space for all downstream analyses.

## 4.12 Polarity of the 25 KCN genes in the Normal vs Activated axis

Main table:

- `code/moffitt_kcn_deg_4_methods/tables/KCN_normal_vs_activated_poles_membership.tsv`

Documented categories:

| Category | Genes |
|---|---|
| `Activated_only` | `KCNAB2`, `KCND2`, `KCND3`, `KCNJ15`, `KCNK1`, `KCNK15`, `KCNMB4`, `KCNN3`, `KCNN4` |
| `Normal_only` | `KCNA5`, `KCNAB1`, `KCNB1`, `KCNE4`, `KCNJ16`, `KCNK17`, `KCNK3`, `KCNQ4`, `KCNS3` |
| `Both` | `KCNA7`, `KCNC4`, `KCNJ8`, `KCNK6`, `KCNMA1`, `KCNMB1`, `KCNT2` |

This result is conceptually important because it shows early that:

- the KCN panel is not a single activated-stroma signature
- some genes look more resident or normal-like
- some genes lie in transitional or context-dependent axes

## 4.13 GSVA Activated vs Normal: why it was included

### Goal

GSVA was not used to select KCN genes. It was used to describe the biological landscape of the stromal states in which KCN genes are embedded.

### Main choices

- pseudobulk by `Patient x StromaState`
- transformation:
  - `log_cpm = log1p((counts / library_size) * 1e6)`
- Hallmark collection from `msigdbr`
- keep pathways with `10-500` overlapping genes
- `kcdf = "Gaussian"`
- paired Wilcoxon by patient
- BH correction

### Why this was useful

This step answers:

> what kind of stromal biological world is represented by Activated stroma?

That is directly relevant for the internship because the topic is immunity plus ECM, not only ion channel expression.

### Main result pattern

Pathways higher in `Activated` stroma:

- `EPITHELIAL_MESENCHYMAL_TRANSITION`
- `INTERFERON_ALPHA_RESPONSE`
- `INTERFERON_GAMMA_RESPONSE`
- `INFLAMMATORY_RESPONSE`
- `COMPLEMENT`
- `ANGIOGENESIS`
- `COAGULATION`
- `KRAS_SIGNALING_UP`
- `HEDGEHOG_SIGNALING`
- `IL2_STAT5_SIGNALING`
- `IL6_JAK_STAT3_SIGNALING`
- `GLYCOLYSIS`
- `APICAL_JUNCTION`

Pathways higher in `Normal` stroma:

- `OXIDATIVE_PHOSPHORYLATION`
- `ADIPOGENESIS`
- `MYC_TARGETS_V1`
- `DNA_REPAIR`
- `MYOGENESIS`
- `FATTY_ACID_METABOLISM`

### Interpretation

Even before assigning each KCN gene, the activated stroma is already strongly compatible with:

- inflammation
- interferon signaling
- complement
- angiogenesis
- coagulation
- EMT-like remodeling

This provides the biological terrain into which KCN candidates are later interpreted.

### Why GSEA was not the main selection tool

Because the aim of the single-cell block was first to identify sparse KCN genes reliably. For that, DEG frameworks and cell-level scoring were more appropriate.

GSEA-like methods were more useful later for:

- pathway-level interpretation
- territory-level reading
- contextualization of already selected KCN genes

## 4.14 AUCell: why AUCell rather than classical GSEA at this stage

### Core logic

AUCell gives a **per-cell score** for a signature. This is much better suited than a classical bulk-style GSEA when asking:

> within a fibroblast subtype, are `KCN-high` cells functionally different from `KCN-low` cells?

### High vs Low definition

- if median expression = `0`:
  - `High = expression > 0`
  - `Low = expression == 0`
- otherwise:
  - `High = expression > median`
  - `Low = expression <= median`

Minimum requirements:

- at least `5` cells in each group per patient
- at least `20` cells total in each group
- at least `3` paired patients

### AUCell settings

- `AUCell_buildRankings(..., plotStats = FALSE)`
- `AUCell_calcAUC(..., aucMaxRank = ceiling(n_genes * 0.05))`

### What it brought

AUCell asked a finer question:

> inside a given stromal subtype, what biological program characterizes the `KCN-high` cells?

### Examples from the documented results

`KCNE4` in `csCAF`:

- higher `EMT`
- higher `TGF_BETA_SIGNALING`
- higher `ANGIOGENESIS`
- higher `COMPLEMENT`
- higher `HYPOXIA`
- higher `MYOGENESIS`

`KCNJ8` in `csCAF`:

- higher `MTORC1_SIGNALING`
- higher `ANGIOGENESIS`
- higher `OXIDATIVE_PHOSPHORYLATION`
- higher `MYC_TARGETS_V1`
- higher `COMPLEMENT`

`KCNK6` in `csCAF`:

- higher `MYOGENESIS`
- higher `CHOLESTEROL_HOMEOSTASIS`
- higher `OXIDATIVE_PHOSPHORYLATION`
- higher `EMT`

### Interpretation

This suggested that some KCN genes are not only markers of a subtype but markers of **functional sub-states within a subtype**.

## 4.15 `myCAF` / `iCAF` correlation analysis

### Goal

Place the 25 KCN genes relative to two major stromal axes that matter for the project:

- `myCAF`: contractile, ECM-rich, exclusion-oriented
- `iCAF`: inflammatory and secretory

### Inputs and method

Main script:

- `code/mycaf_icaf_correlations/run_kcn_mycaf_icaf_correlations.R`

Inputs:

- `kcn_union_25.tsv`
- `genes_top30_myCAF_human.txt`
- `genes_top30_iCAF_human.txt`

Main metric:

- `Spearman rho` between each KCN and each signature gene

Summary metrics:

- `mean_rho`
- `median_rho`
- `mean_abs_rho`
- `min_p_adj_bh`
- `n_fdr_lt_0_05`

### Main significant signals

Representative examples:

| Gene | Program | FDR |
|---|---|---:|
| `KCND2` | `myCAF` | `7.72e-11` |
| `KCNK17` | `myCAF` | `2.63e-09` |
| `KCNAB1` | `myCAF` | `9.22e-08` |
| `KCNC4` | `myCAF` | `2.02e-05` |
| `KCNN4` | `myCAF` | `4.39e-05` |
| `KCNS3` | `myCAF` | `4.99e-05` |
| `KCNE4` | `myCAF` | `5.58e-05` |
| `KCNK3` | `iCAF` | `6.96e-05` |
| `KCND3` | `iCAF` | `1.97e-04` |
| `KCNN3` | `iCAF` | `2.46e-04` |
| `KCNAB2` | `iCAF` | `8.86e-03` |

### Very important nuance on `KCNMA1`

`KCNMA1` later becomes one of the most interesting spatial genes, but in this pseudobulk single-cell correlation block it is **not** a very strong global `myCAF` signal:

- `median rho ~ 0.021`
- only `3/30` significant correlations with the `myCAF` top30 markers

This is a key lesson of the whole project:

> a gene may be biologically very meaningful as a spatial territory marker while remaining modest in a pseudobulk signature correlation framework

### Interpretation

This block was useful because it already separated:

- a more contractile / ECM-like group
- a more inflammatory / resident / border-like group

But it did not fully resolve tissue geography. That is why the project later needed spatial transcriptomics.

## 4.16 External validation in USER

### Goal

Avoid relying on a single stromal atlas.

### Main logic

Instead of blindly transferring labels, the USER fibroblasts were:

1. reconstructed from the local inputs
2. filtered for fibroblasts
3. reclustered
4. rescored using Moffitt-like stromal signatures
5. annotated using:
   - average signature score
   - overlap between cluster DEGs and the signatures

Final retained states:

- `myCAF`
- `iCAF`
- `csCAF`
- `PSC`

### Final USER fibroblast counts

- `myCAF = 11308` (`79.26%`)
- `csCAF = 1326` (`9.29%`)
- `PSC = 1281` (`8.98%`)
- `iCAF = 352` (`2.47%`)

### Marker sanity checks

- `myCAF`: `ACTA2` detected in `44.31%`, `FAP` in `58.57%`
- `csCAF`: `C7` in `68.55%`
- `iCAF`: `IL11` in `5.68%`, `FAP` in `72.44%`

### Why this matters

This validation does not prove each KCN mechanistically, but it strengthens the state framework underlying the rest of the project.

## 4.17 Main conclusions from the single-cell block

The single-cell part does not identify one homogeneous stromal channel program. Instead, it already suggests several biological worlds:

- a ductal / tumor-like world
- a contractile / ECM / myofibroblastic world
- a resident / border / perivascular / inflammatory world

Most importantly, it produced the **25-gene KCN panel** that became the backbone of the rest of the project.

## 4.18 Main limits of the single-cell block

- transcriptomic only, not functional
- `Activated vs Normal` is a useful backbone, not an absolute truth
- subtype-vs-rest contrasts are sensitive to rest-group heterogeneity
- `iCAF` is underpowered
- sparse KCN genes remain model-sensitive
- some genes, especially `KCNMA1`, are not fully captured by subtype-level pseudobulk frameworks

---

## 5. Multicohort bulk survival: clinical relevance of the 25 KCN genes

## 5.1 Why place survival immediately after scRNA-seq?

At this stage of the project, the 25-gene KCN panel has been biologically prioritized in fibroblasts. The next logical question is:

> do these candidates also have clinical significance at patient level?

This is why the multicohort survival block naturally comes **after** the single-cell selection and **before** the local bulk perturbation and spatial niche analyses.

## 5.2 Why survival is useful but not sufficient

Bulk survival does **not** tell us:

- which cell type expresses the gene
- which niche carries the signal
- whether the effect is stromal, tumoral, vascular, or immune

What it does provide is:

- clinical prioritization
- a patient-level relevance layer
- a way to distinguish candidates with possible translational value

## 5.3 Cohorts used

The integrated survival bundle combines four public PDAC cohorts:

- `PH_GSE183795`
- `Hussain_GSE62452`
- `Zhang_GSE28735`
- `Bailey_QCMG_UQ_2016`

Documented survival-eligible sample counts:

- `Bailey`: `96` eligible, `64` events, median OS about `14.98 months`
- `Hussain`: `65` eligible, `49` events, median OS about `14.9 months`
- `PH`: `134` eligible, `86` events, median OS about `11.31 months`
- `Zhang`: `42` eligible, `29` events, median OS about `14.5 months`

Total in the pooled analysis:

- `337` patients
- `228` events

## 5.4 Why a fused survival bundle was built

The cohorts come from different platforms and preprocessing histories:

- microarray
- RNA-seq

Therefore, the analysis did **not** concatenate raw intensities directly.

Instead, the pipeline:

1. loaded all four cohort bundles
2. retained tumor samples with valid survival data
3. defined the common KCN gene space across cohorts
4. standardized each gene **within each cohort**
5. pooled the samples afterward

This preserves within-cohort ranking while reducing platform effects.

## 5.5 Gene space analyzed

The fused bundle is built on a common space of `93` KCN genes.

However, the actual survival analyses focus on:

- the final **25-gene KCN panel** from the single-cell block

This is important conceptually:

- the survival block does not rediscover new KCN genes from scratch
- it tests the clinical weight of biologically curated candidates

## 5.6 Clinical harmonization

Clinical variables were harmonized into simplified factors:

- `stage_group = I / II / III / IV / Unknown`
- `grade_group = Low_1_2 / High_3_4 / Unknown`
- `margin_group = R0 / R1 / Unknown`

The main adjusted Cox model retained:

- `stage_group`
- `grade_group`
- `strata(cohort)`

Margin information was harmonized in the bundle but not kept in the final main model.

## 5.7 Expression standardization

Within each cohort, gene expression was standardized gene by gene:

`z = (x - mean_cohort) / sd_cohort`

Therefore, each survival coefficient is interpreted as:

> effect on risk for `+1 SD` expression increase within the gene's cohort of origin

This is a key methodological point because it avoids naive cross-platform comparison of absolute intensities.

## 5.8 Main adjusted Cox model

Main script:

- `code/bulk_survival_analysis/run_kcn25_multicohort_adjusted_cox_survival.R`

Model:

`Surv(os_months, event) ~ gene_value + stage_group + grade_group + strata(cohort)`

Interpretation:

- `gene_value` = within-cohort z-score expression
- `strata(cohort)` = cohort-specific baseline hazard

This is not a naive pooled Cox model. It allows each cohort to keep its own baseline risk structure.

## 5.9 Main metrics in the adjusted Cox results

The main output table reports:

- `beta`
- `HR = exp(beta)`
- `hr_low`, `hr_high`
- `SE`
- `z`
- `p_value`
- `p_adj_bh`
- `concordance`

Interpretation:

- `HR > 1`: high expression associated with poorer prognosis
- `HR < 1`: high expression associated with better prognosis
- `p_adj_bh`: multiple-testing corrected significance across the 25 genes

## 5.10 Significant genes in the adjusted Cox model

Significant after BH correction:

| Gene | HR | 95% CI | p-value | FDR | Interpretation |
|---|---:|---|---:|---:|---|
| `KCNMB4` | `1.283` | `[1.102 ; 1.494]` | `0.00130` | `0.0177` | unfavorable if high |
| `KCNC4` | `0.802` | `[0.699 ; 0.920]` | `0.00159` | `0.0177` | protective if high |
| `KCNAB1` | `0.807` | `[0.703 ; 0.925]` | `0.00212` | `0.0177` | protective |
| `KCNE4` | `0.808` | `[0.699 ; 0.933]` | `0.00379` | `0.0237` | protective |
| `KCNK6` | `1.202` | `[1.057 ; 1.367]` | `0.00509` | `0.0254` | unfavorable |
| `KCNT2` | `0.848` | `[0.744 ; 0.965]` | `0.0126` | `0.0478` | protective |
| `KCNN3` | `0.834` | `[0.723 ; 0.963]` | `0.0134` | `0.0478` | protective |

Examples of non-significant but still interesting trends:

- `KCNMA1`: `HR = 0.892`, `p = 0.077`, `FDR = 0.241`
- `KCNK17`: `HR = 0.890`, `p = 0.106`
- `KCNK1`: `HR = 0.899`, `p = 0.107`
- `KCNK3`: `HR = 0.902`, `p = 0.158`

## 5.11 Biological interpretation of the main survival results

The bulk survival layer does not mirror the spatial layer perfectly.

For example:

- `KCNMA1` becomes one of the most compelling stromal spatial candidates later
- but it is only a trend in the adjusted Cox model

Conversely:

- `KCNMB4` is the strongest bulk survival hit
- while spatially it is much less dominant than `KCNMA1`

This is a very important conceptual result:

> clinical bulk relevance and tissue niche identity are not the same thing

## 5.12 Proportional hazards checks

Additional check table:

- `code/bulk_survival_analysis/survival_kcn25_venn_adjusted_cox/tables/kcn25_multicohort_cox_zph_checks.tsv`

Main outputs:

- `gene_ph_p`
- `global_ph_p`
- `gene_ph_fdr`
- `global_ph_fdr`

Several genes violate proportional hazards assumptions, including:

- `KCNS3`
- `KCNJ8`
- `KCNB1`
- `KCNK1`
- `KCNK15`
- `KCNA7`
- `KCNQ4`

This is precisely why a time-varying Cox extension was added.

## 5.13 Exploratory KM cutoff scan

Main script:

- `code/bulk_survival_analysis/run_kcn25_multicohort_km_scan_min15_survival.R`

Goal:

- test whether some genes show a clearer `high vs low` survival separation at a specific cutoff

Constraints:

- `min_group_n = 15`
- `min_events_per_group = 5`

Corrections:

- BH within gene across cutoffs
- BH across genes on the best cutoff per gene

### Strongest KM scan results

- `KCNMB4`: cutoff `z = 1.013`, `49 high vs 288 low`, `raw p = 5.56e-11`, `BH across genes = 1.39e-09`, median OS `8.88 vs 20.14 months`
- `KCNAB1`: `raw p = 3.26e-06`, `BH across genes = 4.08e-05`
- `KCNA5`: `raw p = 9.55e-06`, `BH across genes = 7.96e-05`
- `KCNC4`: `raw p = 4.62e-05`, `BH across genes = 2.89e-04`
- `KCNE4`: `raw p = 1.14e-04`, `BH across genes = 5.69e-04`
- `KCNN3`: `raw p = 3.17e-04`, `BH across genes = 0.00132`
- `KCNK17`: `raw p = 4.29e-04`, `BH across genes = 0.00153`

Interpretation:

- this is a useful visual and exploratory layer
- but it remains more optimistic because it searches for a favorable cutoff

## 5.14 Time-varying Cox model

Main script:

- `code/bulk_survival_analysis/run_kcn25_multicohort_time_varying_cox_survival.R`

Model:

`Surv(os_months, event) ~ gene_value + tt(gene_value) + stage_group + grade_group + strata(cohort)`

with:

`tt(x, t) = x * (log(t) - mean(log(t)))`

Goal:

- test whether gene effect changes over follow-up time

Main outputs:

- `beta_main`
- `beta_time_interaction`
- `lrt_p_value`
- `lrt_p_adj_bh`
- time-specific HR trajectories at:
  - `6`
  - `12`
  - `24`
  - `36`
  - `48`
  - `60` months

### Significant time-varying effects after BH

- `KCNB1`: `LRT FDR = 0.00245`
- `KCNS3`: `LRT FDR = 0.00245`
- `KCNJ8`: `LRT FDR = 0.00271`
- `KCNK1`: `LRT FDR = 0.0148`
- `KCNK15`: `LRT FDR = 0.0148`
- `KCNA7`: `LRT FDR = 0.0440`

Examples:

- `KCNJ8`: HR about `1.56` at `6 months`, then about `0.82` at `60 months`
- `KCNB1`: HR about `1.45` at `6 months`, then about `0.83` at `60 months`
- `KCNK1`: HR about `0.60` at `6 months`, then about `0.99` at `60 months`

Interpretation:

- some KCN genes are not well summarized by a single constant hazard ratio
- their association with risk may weaken, reverse, or flatten over time

## 5.15 What the survival block adds to the project

This block does not tell us where a KCN gene lives in the tissue. It tells us:

- whether it has patient-level relevance
- whether the effect is continuous
- whether a visual `high vs low` split exists
- whether the effect is stable over time

This is why it belongs after the single-cell selection but before the spatial block.

## 5.16 Main limits of the survival block

- bulk mixes tumor, stroma, immune, and vascular signals
- within-cohort z-scoring reduces platform effects but does not remove every bias
- the main model adjusts only for `stage` and `grade`
- several genes violate proportional hazards assumptions
- KM scans are exploratory and can inflate apparent separation
- survival is a prioritization layer, not a mechanistic proof

---

## 6. Local bulk perturbation datasets from the laboratory

## 6.1 Why keep these local bulk analyses after the multicohort survival block?

The multicohort survival layer gives patient-level relevance. The local bulk perturbation datasets answer a different question:

> when a biologically meaningful laboratory axis is perturbed, do coherent transcriptomic programs emerge around the channel of interest?

These bulk datasets are not central for the 25-KCN selection itself, but they help connect the project to the biological context of the laboratory.

## 6.2 BKCa / shBKCa

### Dataset

- `8` samples
- `4 CTRL`
- `4 shBKCa`
- `17889` collapsed genes

### Main results

- `6734` DEGs at `padj < 0.05`
- top CAF signature by Wilcoxon: `ifCAF_manual`
- strongest positive correlation with `KCNMA1`:
  - `FN1`, `rho = 0.952`, `FDR = 0.0348`
- strongest negative correlation:
  - `MGST1`, `rho = -1`, `FDR = 0.0030`
- strongest Hallmark:
  - `INTERFERON_ALPHA_RESPONSE`, `NES = -2.539`, `FDR = 3.9e-11`

### Interpretation

This block is important because it gives an external coherence to the `KCNMA1 / BKCa` axis:

- strong link to `FN1`
- matrix-oriented signal
- immune-inflammatory pathway involvement

Even if the data are bulk and not stromally pure, they support the broader channel-ECM-immunity framework.

## 6.3 SIGMAR1 / shSIGMAR1

### Dataset

- `6` samples
- `3 CTRL`
- `3 shSIGMAR1`
- `14729` collapsed genes

### Main results

- `1546` DEGs at `padj < 0.05`
- top CAF signature: `myCAF_top30`, but not significant after FDR
- strongest positive correlation with `SIGMAR1`:
  - `EZR`, `rho = 0.943`
- strongest negative correlation:
  - `CXADR`, `rho = -1`
- strongest Hallmark:
  - `G2M_CHECKPOINT`, `NES = -2.047`, `FDR = 1.8e-10`

### Interpretation

This block is less directly tied to the 25 KCN panel, but it remains useful because it connects the project to a broader ion-channel signaling context relevant to the host laboratory.

## 6.4 Main limits of the local bulk perturbation block

- very small sample sizes
- one experimental context per axis
- bulk mixtures without pure stromal resolution
- good for coherence, weak for definitive prioritization alone

---

## 7. Spatial transcriptomics: tissue localization, niches, and programs

## 7.1 Why the spatial layer was necessary

After scRNA-seq and bulk, the missing question was:

> where do these KCN genes actually fall in PDAC tissue?

This is essential because in PDAC:

- immunity depends on tissue architecture
- ECM is spatially organized
- interfaces matter
- perivascular and border territories matter

The spatial layer is therefore where the project moves from **candidate genes** to **candidate niches**.

## 7.2 Spatial object and global metrics

Main object:

- `data/PDAC_Updated_ST.rds`

Documented global metrics:

- `91,496` spots
- `30` sections
- `8` assays
- default assay: `rctd_fullfinal`

Origin distribution:

- `Pancreas`: `35,458` spots (`38.75%`)
- `Liver`: `28,520` (`31.17%`)
- `Lymph node`: `17,698` (`19.34%`)
- `Normal Pancreas`: `9,820` (`10.73%`)

Sections used most heavily in this project:

- `10` primary PDAC sections:
  - `IU_PDA_T1`
  - `IU_PDA_T2`
  - `IU_PDA_T3`
  - `IU_PDA_T4`
  - `IU_PDA_T6`
  - `IU_PDA_T8`
  - `IU_PDA_T9`
  - `IU_PDA_T10`
  - `IU_PDA_T11`
  - `IU_PDA_T12`
- `3` normal pancreas sections:
  - `IU_PDA_NP2`
  - `IU_PDA_NP10`
  - `IU_PDA_NP11`

## 7.3 What kind of labels already exist in the spatial object

Main label layers used in the project:

- `cc_ischia_10` / `CompositionCluster_CC` for ecotypes
- `first_type` for dominant cell type per spot
- RCTD-related assays for composition information

The object is therefore already rich enough to support:

- ecotype-oriented analyses
- dominant-cell-type neighborhood analyses
- pathway scoring overlays

## 7.4 Global tissue composition and why it matters

Dominant `first_type` counts include:

- `myCAF`: `26,726` (`29.21%`)
- `Tumor Epithelial cells`: `24,627` (`26.92%`)
- `Normal Epithelial cells`: `17,997` (`19.67%`)
- `Hepatocytes`: `9,668` (`10.57%`)
- `B cells`: `4,407` (`4.82%`)
- `iCAF`: `2,246` (`2.45%`)

Main ecotypes include:

- `CC12`: `11.22%`
- `CC15`: `10.61%`
- `CC2`: `10.60%`
- `CC14`: `7.85%`
- `CC1`: `7.76%`
- `CC8`: `7.69%`
- `CC5`: `5.99%`
- `CC3`: `3.50%`

Spot class confidence:

- `doublet_certain`: `69.79%`
- `singlet`: `16.95%`
- `doublet_uncertain`: `12.12%`
- `reject`: `1.14%`

This is a crucial limitation and interpretation rule:

- Visium spots are mixed
- spatial conclusions should be read as **territory-level** conclusions, not strict one-cell-one-gene claims

## 7.5 First spatial level: maps and marker projections

The first step was descriptive but essential:

- ecotype maps
- `first_type` maps
- KCN projections
- marker projections

Reference markers used to interpret tissue structure:

- stromal contractile / ECM:
  - `ACTA2`
  - `CCN2`
  - `POSTN`
  - `TAGLN`
  - `COL12A1`
  - `THBS1`
- epithelial / ductal:
  - `EPCAM`
  - `KRT20`
- resident / immune-border:
  - `C7`
  - `C1R`
  - `CCL19`
  - `CCL21`
  - `MFAP4`
  - `CLU`
  - `DCN`

These maps do not prove a niche by themselves, but they are necessary to learn how the tissue is organized and to read later quantitative results correctly.

## 7.6 Ecotype detection by positivity

An early spatial analysis asked:

> in which ecotypes are KCN genes more often detected?

This analysis used:

- `Origin %in% c("Normal Pancreas", "Pancreas")`
- union of the 25 KCN genes
- expression in `Spatial/data`
- positivity defined as `expression > 0`

This first level measured:

- **percentage of positive spots**

This is informative for sparse KCN genes, but it does not yet capture intensity or section heterogeneity.

## 7.7 Fisher exact niche enrichment: useful first step, but later insufficient

### Initial niche definitions

- `CC1+CC5`
- `CC2+CC3`

### Initial question

For each `gene x niche` in primary pancreas:

> are KCN-positive spots overrepresented inside the niche compared with outside?

### Test

For each gene and niche:

- build a `2x2` table:
  - inside / outside niche
  - KCN+ / KCN-
- run `fisher.test()`

Main outputs:

- `odds ratio`
- `p-value`
- `BH FDR`
- `% positive inside`
- `% positive outside`

### Main early results

In `CC2+CC3`:

- `KCNK1`: `OR = 4.29`, `FDR = 0`
- `KCNN4`: `OR = 2.88`, `FDR = 1.09e-211`
- `KCNJ16`: `OR = 6.02`, `FDR = 9.26e-73`
- `KCNC4`: `OR = 2.06`, `FDR = 3.67e-58`
- `KCNQ4`: `OR = 2.24`, `FDR = 8.21e-29`
- `KCNK6`: `OR = 1.52`, `FDR = 3.73e-29`

In `CC1+CC5`:

- `KCNMA1`: `OR = 2.39`, `FDR = 8.20e-206`
- `KCNE4`: `OR = 2.32`, `FDR = 6.46e-113`
- `KCNMB1`: `OR = 2.00`, `FDR = 3.60e-63`
- `KCND2`: `OR = 1.83`, `FDR = 2.01e-27`

### Why this was a good first step

It gave a very simple and defensible first answer:

> positive spots for this KCN gene are more frequent in this niche than expected by chance

### Why this later became insufficient

The Fisher framework remains limited because:

- it binarizes expression
- it pools information too coarsely
- it ignores expression intensity
- it ignores section-to-section heterogeneity
- it ignores spatial dependence explicitly

This is why the project later replaced it with a more informative **per-section Wilcoxon / log2FC / AUC** framework.

## 7.8 Per-section ecotype vs rest analysis: Wilcoxon, AUC, log2FC, FDR

This is one of the major updates added after the first synthesis.

### Goal

For each primary section and each KCN gene:

- compare `CC1+CC5` vs the rest of the section
- compare `CC2+CC3` vs the rest of the section

### Main metrics

- `log2FC`
- `AUC`
- `%Pos in`
- `%Pos out`
- `dPct = %Pos in - %Pos out`
- Wilcoxon `p-value`
- `FDR`

### How to interpret these metrics

- `log2FC > 0`: higher expression in the tested niche
- `AUC > 0.5`: spots from the niche tend to have higher expression than the rest
- `FDR < 0.05`: statistically supported within the section
- number of significant sections: robustness across tissue slices

### Why this approach is better than Fisher

It retains:

- expression intensity
- section-level heterogeneity
- a more informative effect size
- a repeated structure across the 10 primary sections

### Main KCN results

Strong `CC2+CC3` block:

- `KCNK1`: `10/10` positive sections, `10/10` `FDR < 0.05`, `median log2FC = 0.414`, `median AUC = 0.610`
- `KCNN4`: `10/10` positive, `10/10` significant, `median log2FC = 0.257`, `median AUC = 0.581`
- `KCNK6`: `10/10` positive, `8/10` significant, `median AUC = 0.537`

More heterogeneous `CC1+CC5` signals:

- `KCNMA1`: `5/10` sections with positive `log2FC`, `3/10` with `FDR < 0.05`, `median log2FC ~ 0`, `median AUC = 0.503`
- `KCNMB1`: `4/10` positive, `3/10` significant, `median AUC = 0.499`
- `KCNE4`: `9/10` positive, `5/10` significant, `median log2FC = 0.033`, `median AUC = 0.512`

### Important interpretation

This analysis gave a decisive conceptual result:

- `KCNN4` and `KCNK1` behave like strong ecotype-associated genes
- `KCNMA1` does **not** behave like a clean ecotype marker

This is not a failure. It suggests that `KCNMA1` is better captured by:

- spatial patches
- local programs
- continuous stromal territories

rather than by a discrete ecotype label alone.

## 7.9 Marker controls in the same ecotype-vs-rest framework

A control analysis was run with canonical stromal markers:

- `ACTA2`
- `CCN2`
- `POSTN`
- `TAGLN`
- `COL12A1`

Main `CC1+CC5` summary results:

- `CCN2`: `median log2FC = 0.552`, `median AUC = 0.658`, `7/10` sections significant
- `ACTA2`: `0.434`, `0.630`, `10/10`
- `POSTN`: `0.429`, `0.606`, `10/10`
- `TAGLN`: `0.333`, `0.599`, `10/10`
- `COL12A1`: `0.300`, `0.574`, `8/10`

### Why this was very important

It showed two things at once:

1. the ecotype-vs-rest approach does capture a real stromal signal
2. even classical myCAF / ECM markers are not perfectly pure or binary in these sections

This directly supports the interpretation that the problem is **not** that `KCNMA1` is irrelevant. The problem is that a binary ecotype-vs-rest contrast is structurally too coarse for a continuous stromal territorial signal.

## 7.10 Hotspot: which KCN genes form true spatial patches?

### Main script

- `code/spatial_pdac_analysis/communication/hotspot/all_kcns_primary/run_hotspot_all_kcns_primary.py`

### Why Hotspot was needed

Wilcoxon ecotype contrasts ask:

- is the gene higher in one niche than in the rest?

Hotspot asks:

- does the gene form a true local spatial pattern in the tissue?

### Main parameters

- `10` primary sections
- model: `bernoulli`
- `n_neighbors = 30`
- `MIN_GENE_SPOTS = 10`
- color scale upper bound at quantile `0.95` for display

### Main metrics

- `C`: local spatial autocorrelation
- `Z`: standardized spatial signal
- `FDR`

Important:

- the map color represents raw counts
- not the `Z` score itself

### Main Hotspot summary

| Gene | Tested sections | Significant sections | median Z | Interpretation |
|---|---:|---:|---:|---|
| `KCNMA1` | 10 | 7 | `7.01` | very strong niche gene |
| `KCNK1` | 10 | 7 | `3.52` | strong spatial niche gene |
| `KCNMB1` | 10 | 7 | `3.46` | strong spatial niche gene |
| `KCNK3` | 10 | 6 | `2.63` | robust but more resident |
| `KCNN4` | 10 | 5 | `2.46` | robust, more tumor/ductal |
| `KCND3` | 10 | 5 | `2.00` | recurrent, weaker niche |

### Why this result was decisive

It clearly showed that:

- `KCNMA1` is one of the strongest spatial genes in primary PDAC
- even if it is not one of the best binary ecotype markers

This is one of the most important results of the whole project.

## 7.11 SPARK-X: independent spatial variability validation

### Main script

- `code/spatial_pdac_analysis/sparkx/run_sparkx_all_kcns_primary.R`

### Why add SPARK-X?

Hotspot is excellent for compact local patches. SPARK-X gives a complementary test of spatial dependence using mixed kernels and can capture smoother or more diffuse patterns.

### Main settings

- section-wise analysis
- raw counts plus coordinates
- `sparkx(..., option = "mixture")`
- multiple kernels
- multiple-testing correction by `Benjamini-Yekutieli`

### Main SPARK-X summary

- `KCNN4`: `10/10` sections significant, `median adjusted p = 2.97e-18`, `median best kernel stat = 87.7`
- `KCNMA1`: `10/10`, `5.21e-10`, `49.9`
- `KCNK6`: `10/10`
- `KCNJ8`: `10/10`
- `KCNK1`: `9/10`
- `KCNS3`: `9/10`
- `KCNAB2`: `9/10`

### Interpretation

SPARK-X strongly reinforces:

- the ductal/tumor block
- the contractile/ECM block
- part of the resident / border block

Convergence of Hotspot and SPARK-X is especially convincing.

## 7.12 Local correlations: with which genes does each KCN gene spatially co-vary?

### Main script

- `code/spatial_pdac_analysis/communication/hotspot/all_kcns_primary/local_correlations/run_kcn_local_correlations_primary.py`

### Logic

For each KCN gene:

- keep spatially informative genes (`FDR < 0.05`, `C > 0`)
- compute local correlations in the same Hotspot framework
- summarize across sections by:
  - `n_sections_positive_z`
  - `median_local_corr_z`

### Why this analysis is powerful

It answers:

> with which local molecular program does this KCN gene live in tissue?

This is not the same as a global correlation in the entire dataset. Two genes can share a local spatial program even if their overall global correlation is only moderate.

### Main examples

#### `KCNN4`

Top recurrent local partners include:

- `C19orf33`
- `KCNK1`
- `ELF3`
- `LCN2`
- `SLPI`
- `S100P`
- `KRT7`
- `MUC1`
- `ITGB4`
- `TSPAN8`
- `KRT18`

Interpretation:

- clear **tumor epithelial / ductal territory**

#### `KCNK1`

Strong recurrent pairs:

- `KCNK1-SLPI`: `10/10`, `median z = 4.62`
- `KCNK1-ELF3`: `10/10`, `4.49`
- `KCNK1-LCN2`: `10/10`, `3.80`

Interpretation:

- same ductal / epithelial world as `KCNN4`

#### `KCNMA1`

In the detailed per-section moffitt_stromal_exploration_outputs, especially in `T1`, top partners include:

- `GREM1`
- `SLC2A3`
- `SGK1`
- `CNN1`
- `FOSB`
- `KLF6`
- `SERPINE1`
- `THBS1`
- `LOXL1`
- `COL12A1`
- `COL11A1`
- `CTHRC1`
- `MMP2`
- `CCN2`

Interpretation:

- strongly compatible with a **contractile / myofibroblastic / ECM remodeling territory**

#### Resident / immune-border genes

Examples:

- `KCNK3-C7`: `10/10`, `median z = 2.24`
- `KCNJ8-C7`: `10/10`, `1.35`
- `KCND3-C7`: `9/10`, `3.10`
- `KCNK3-CLU`: `9/10`, `2.91`

Interpretation:

- resident / perivascular / immune-border-like programs

## 7.13 Territory GSEA on KCN-high spots

### Biological question

After identifying spatial patches, the next question was:

> what biological program characterizes the territory where the KCN signal is strongest?

### Main script

- `code/spatial_pdac_analysis/communication/hotspot/all_kcns_primary/territory_gsea/run_all_kcn_territory_gsea.R`

### Territory definition

For each primary section:

- positive spots: `KCN > 0` in `SCT/data`
- at least `30` positive spots required
- `KCN-high` = top `10%` of positive spots
- at least `20` high spots required

### Why this definition is useful

It avoids:

- collapsing the whole tissue into one contrast
- defining the territory only by presence/absence
- overinterpreting sections with very few positive spots

### GSEA framework

- assay used for ranking: `SCT`
- gene score = `mean(high) - mean(rest)`
- Hallmarks via `msigdbr`
- overlap minimum `10`
- `fgseaMultilevel(eps = 0)`

### Main results by KCN world

#### Stromal contractile / ECM world

`KCNMA1`:

- `10` usable sections
- `21` significant pathways total
- `3` positive and `18` negative
- top positive:
  - `HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION`
  - `NES = 1.63`
  - `FDR = 3.55e-08`
- top negative:
  - `HALLMARK_MYC_TARGETS_V1`
  - `NES = -1.70`
  - `FDR = 8.42e-15`

`KCNMB1`:

- top positive:
  - `HALLMARK_MYOGENESIS`
  - `NES = 2.13`

`KCNE4`:

- top positive:
  - `HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION`
  - `NES = 2.28`

Interpretation:

- strong compatibility with contractile / ECM / myofibroblastic territory biology

#### Tumor ductal / epithelial world

`KCNN4`:

- `16` significant pathways
- top negative:
  - `HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION`
  - `NES = -2.47`
  - `FDR = 1.35e-22`

`KCNK1`:

- top negative:
  - `EMT`
  - `NES = -2.46`
  - `FDR = 1.30e-24`

Interpretation:

- these territories oppose stromal contractile biology and fit much better with ductal tumor programs

#### Resident / border world

`KCND3`, `KCNJ8`, `KCNAB2`:

- more mixed and less massive signals
- consistent with resident / perivascular / border-like contexts

## 7.14 COMMOT

COMMOT was tested on selected candidates to ask whether KCN-high territories sat within particularly strong secreted signaling contexts.

Main result:

- no dominant secreted-signaling hyper-communicative niche emerged as a major explanation

Why this negative result still matters:

- it suggests the most informative biology of these KCN genes is more about:
  - spatial architecture
  - ECM
  - contractility
  - tissue interface organization

than about a broad increase in secreted-signaling communication per se.

## 7.15 Neighborhood analysis using `first_type`

### Main script

- `code/spatial_pdac_analysis/cell_type_neighborhood/run_kcn_high_first_type_knn_neighborhood.py`

### Question

Around `KCN-high` spots, which dominant spot types are more frequently found?

### Main settings

- high spots = top `10%` of positive spots
- at least `20` high spots
- at least `30` positive spots
- `K = 15` nearest neighbors
- compare neighborhood around high spots vs other spots

Main summary metrics:

- `median_high_neighbor_fraction`
- `median_other_neighbor_fraction`
- `median_delta_neighbor_fraction`
- `median_enrichment_ratio`
- Wilcoxon p-values and FDRs

### Main results

`KCNN4-high`:

- enriched neighborhood toward `Tumor Epithelial cells`
- `9/10` sections with positive delta
- `median delta neighbor fraction = 0.0323`

`KCNK1-high`:

- tumor epithelial enrichment as well
- `7/10` positive

`KCNMA1-high`:

- depleted in epithelial neighbors:
  - `Normal Epithelial cells`: `median delta = -0.058`
  - `Tumor Epithelial cells`: `-0.039`
- positive tendency toward `myCAF`:
  - `median delta = 0.0558`

`KCND3-high`:

- weaker and more mixed tendencies
- compatible with a resident / immune-border interpretation

### Interpretation

This does not prove that the neighboring `first_type` expresses the KCN gene. It shows the **local tissue context** of the KCN-high territory.

## 7.16 Continuous pathway overlays with `escape`: why they were added

The ecotype-vs-rest logic remained unsatisfying for `KCNMA1`. This motivated a new layer:

- left panel: the KCN Hotspot-style spatial projection
- right panel: a pathway score per spot calculated with `escape`

Then, for each section:

- pathway scores were standardized
- a `Spearman` correlation was computed between:
  - the raw counts of the KCN gene
  - the standardized pathway score

This is not a formal local co-localization metric, but it is a strong **spot-level continuous concordance** measure.

## 7.17 `KCNMA1` pathway overlays

Main pathways tested:

- `Cancer_associated_fibroblasts`
- `Matrix`
- `Matrix_remodeling`
- `myCAF_contractile`
- control: `NK_cells`

Main summary results:

- `Cancer_associated_fibroblasts`: `median rho = 0.154`, `10/10` sections `FDR < 0.05`
- `Matrix`: `0.112`, `10/10`
- `Matrix_remodeling`: `0.116`, `10/10`
- `myCAF_contractile`: `0.147`, `10/10`, `max rho = 0.402`
- `NK_cells`: `median rho = -0.002`, `5/10`

### Interpretation

This was a major conceptual advance:

- `KCNMA1` does not simply follow "any" pathway
- it aligns much better with a **continuous stromal contractile / CAF / ECM program**
- the best interpretation is not "KCNMA1 is a perfect ecotype marker"
- the best interpretation is:
  - `KCNMA1` marks a **continuous contractile stromal territory**

## 7.18 `KCNN4` pathway overlays

Main tumor-oriented modules tested:

- `ductal_tumoral_epithelial`
- `tumor_stress_secretory`

Main summary results:

- `ductal_tumoral_epithelial`: `median rho = 0.235`, `10/10` sections significant, `max rho = 0.423`
- `tumor_stress_secretory`: `median rho = 0.215`, `10/10`, `max rho = 0.375`

### Interpretation

This is one of the cleanest validations in the project:

- `KCNN4` strongly follows a continuous tumor ductal / epithelial program
- more than it follows a stromal one

## 7.19 `KCND3` pathway overlay

Main module tested:

- `iCAF_immune_border`

Main summary result:

- `median rho = 0.064`
- `8/10` sections with `FDR < 0.05`
- `max rho = 0.120`

### Interpretation

This is weaker than `KCNMA1` or `KCNN4`, but still informative:

- `KCND3` does not strongly overlap one sharply defined module
- instead it shows a weak but recurrent positive association
- this supports a looser resident / immune-border interpretation

## 7.20 Spatial marker panel rerun in Hotspot

A dedicated Hotspot run was also performed for reference markers such as:

- `ACTA2`
- `CCN2`
- `POSTN`
- `TAGLN`
- `COL12A1`

This was important because it showed that large canonical stromal markers can form smoother and visually more graded patterns than some KCN genes.

It helped explain why:

- `KCNMA1` may look visually more binary in some sections
- even when it is statistically a very real spatial patch

## 7.21 Normal pancreas spatial extension

Additional normal-pancreas analyses were added:

- Hotspot
- local correlations
- territory GSEA

Main normal sections:

- `IU_PDA_NP10`
- `IU_PDA_NP11`
- `IU_PDA_NP2`

### Main lessons from the normal-pancreas extension

- several KCN genes already have structured spatial patterns in the normal pancreas
- this prevents overinterpreting every signal as purely tumor-specific

Examples:

- `KCNJ8` and `KCNK3` support a resident / perivascular axis that exists outside tumor
- `KCNK1` keeps an epithelial anchoring in normal tissue, but not with the same program as in PDAC
- `KCNMA1` and `KCNE4` show that part of the stromal/ECM logic already exists in normal pancreas but is amplified or reorganized in PDAC

## 7.22 Targeted spatial extension: `SIGMAR1`

`SIGMAR1` was analyzed in the same spatial framework without artificially forcing it into the KCN panel.

Main branches:

- Hotspot
- local correlations
- neighborhood
- territory GSEA

Main cautious interpretation:

- `SIGMAR1` is weaker than the strongest KCN spatial genes
- it does not define a major robust niche comparable to `KCNMA1` or `KCNN4`
- its signal appears more diffuse and heterogeneous

## 7.23 Targeted spatial extension: `KCNMA1 + KCNMB1` and `KCNMA1 + KCNMB4`

Two pseudo-modules were built at spot level:

- `KCNMA1_KCNMB1 = KCNMA1 + KCNMB1`
- `KCNMA1_KCNMB4 = KCNMA1 + KCNMB4`

Goal:

- test whether the BK axis is more interpretable when `KCNMA1` is read together with a beta subunit

Main interpretation:

- both pseudo-modules clearly fall in a `myCAF` / contractile / ECM world
- `KCNMA1 + KCNMB1` is more spatially convincing than `KCNMA1 + KCNMB4`

This supports the idea that:

- the most coherent spatial BK partner for `KCNMA1` in this project is `KCNMB1`

## 7.24 Main conclusions from the spatial block

The spatial layer reorganized the whole biological interpretation.

It showed that:

- `KCNMA1` is best read as a contractile stromal territory gene
- `KCNN4` and `KCNK1` are best read as tumor ductal territory genes
- `KCNJ8`, `KCND3`, `KCNK3`, and `KCNAB2` support a resident / border / perivascular world

Most importantly, the spatial block explains why some genes are not well summarized by single-cell subtype labels or by bulk survival alone.

## 7.25 Main limits of the spatial block

- Visium remains spot-level, not pure single-cell
- doublet and mixed spots are common
- spatial co-organization does not prove cell-intrinsic co-expression
- Fisher, Hotspot, SPARK-X, local correlations, and overlays do not imply causality
- pathway overlays rely on curated modules and spot-level correlations, not on direct mechanistic validation

---

## 8. Integrated biological synthesis across all layers

## 8.1 The 25 KCN panel is not one single program

This is one of the strongest conclusions of the project.

The KCN panel separates into at least three biologically meaningful blocks.

## 8.2 Tumor ductal / epithelial block

Main genes:

- `KCNK1`
- `KCNN4`
- `KCNK6`
- partially `KCNC4`, `KCNQ4`, `KCNJ15`

Converging evidence:

- enrichment in `CC2+CC3`
- strong spatial variability
- local correlations with:
  - `ELF3`
  - `LCN2`
  - `SLPI`
  - `KRT7`
  - `KRT18`
  - `MUC1`
- territory GSEA opposing EMT / ECM-like programs
- tumor-oriented pathway overlays

Biological reading:

- these genes are not primarily fibroblast-only markers
- they are linked to tumor epithelial or ductal interface territories
- they remain highly relevant because tumor-stroma interfaces are crucial for local immune and ECM organization

## 8.3 Stromal contractile / ECM / myofibroblastic block

Main genes:

- `KCNMA1`
- `KCNMB1`
- `KCNE4`
- `KCND2`
- partially `KCNS3`

Converging evidence:

- `myCAF`-oriented interpretation in single-cell
- enrichment tendency in `CC1+CC5`
- strong Hotspot for `KCNMA1` and `KCNMB1`
- local correlations with:
  - `THBS1`
  - `CCN2`
  - `COL12A1`
  - `CTHRC1`
  - `ACTG2`
  - `CNN1`
- territory GSEA toward:
  - `EMT`
  - `MYOGENESIS`
  - ECM-rich biology
- pathway overlays toward:
  - CAF
  - Matrix
  - Matrix remodeling
  - `myCAF_contractile`

Biological reading:

- this is the block most directly relevant to **ECM plus indirect immune modulation**
- these genes are the strongest candidates for stromal mechanistic follow-up

## 8.4 Resident / perivascular / immune-border block

Main genes:

- `KCNJ8`
- `KCNK3`
- `KCND3`
- `KCNAB2`
- `KCNN3`

Converging evidence:

- more resident or mixed single-cell behavior
- local correlations with:
  - `C7`
  - `CCL19`
  - `CCL21`
  - `MFAP4`
  - `CLU`
- strong or moderate spatial variability
- weaker but recurrent `iCAF_immune_border`-like patterns for `KCND3`
- neighborhood patterns compatible with border, immune, or perivascular contexts

Biological reading:

- these genes are less purely contractile than the `KCNMA1` block
- they are particularly interesting for the question of stromal immune geography and border territories

---

## 9. Most convincing candidates in the project

If one combines:

- robustness across methods
- biological coherence
- spatial evidence
- relevance to the internship question

the strongest candidates are:

1. `KCNMA1`
2. `KCNMB1`
3. `KCNK1`
4. `KCNN4`
5. `KCNE4`
6. `KCND2`
7. `KCNJ8`
8. `KCNK3`

### Practical reading table

| Gene | Main reading | Link to immunity / ECM question |
|---|---|---|
| `KCNMA1` | contractile stromal territory | very strong |
| `KCNMB1` | contractile / perivascular stromal territory | very strong |
| `KCNE4` | activated fibro / ECM | strong |
| `KCND2` | myCAF / ECM / EMT-like | strong |
| `KCNK1` | tumor ductal / epithelial interface | strong but more indirect |
| `KCNN4` | tumor ductal / epithelial interface | strong but more indirect |
| `KCNJ8` | resident / perivascular | strong for border / territory logic |
| `KCNK3` | resident / iCAF-like / perivascular | strong for inflammatory-border logic |
| `KCND3` | weaker immune-border / resident context | moderate but coherent |
| `KCNAB2` | immune-border / resident context | promising but less robust |

---

## 10. Global robustness of the project

The project is robust as a multi-layer exploratory bioinformatic study because it combines:

- multiple data modalities
- multiple DEG frameworks
- external validation
- several independent spatial methods
- both biological and clinical interpretation layers

The most robust conclusions are:

- the existence of a final 25-gene KCN panel
- the separation into at least three major biological blocks
- the strong spatial role of:
  - `KCNMA1`
  - `KCNMB1`
  - `KCNK1`
  - `KCNN4`
- the strong contractile / ECM reading of:
  - `KCNMA1`
  - `KCNMB1`
  - `KCNE4`
  - `KCND2`

More exploratory elements remain:

- weaker 1/4 support genes
- some local bulk perturbation interpretations
- fine-grained mechanistic interpretation of immune-border genes

---

## 11. Main limitations

## 11.1 Conceptual limitations

- the project is transcriptomic, not functional
- it identifies associations and priorities, not causal mechanisms

## 11.2 Single-cell limitations

- subtype labels remain abstractions of biological continua
- `iCAF` is underrepresented
- sparse KCN genes remain model-sensitive

## 11.3 Bulk survival limitations

- signal is mixed across tumor, stroma, and immunity
- platform harmonization reduces but does not remove all biases
- some genes violate proportional hazards assumptions

## 11.4 Local bulk perturbation limitations

- small sample sizes
- no stromal purification
- laboratory-context dependent

## 11.5 Spatial limitations

- Visium spots are mixed
- spatial co-organization is not equivalent to same-cell co-expression
- pathway overlays are informative but not direct causal readouts

## 11.6 Biological limitation

The project does **not** directly test:

- immune-cell recruitment after channel perturbation
- ECM remodeling after channel perturbation
- electrophysiological channel activity
- direct co-culture or functional interaction assays

Therefore, the project prioritizes candidates and niches rather than closing the mechanistic question.

---

## 12. General conclusion

The project answers the internship question in a nuanced but biologically strong way.

It shows that:

- the stromal KCN landscape in PDAC is structured rather than uniform
- the most relevant candidates are distributed across distinct tissue territories
- the strongest axes for the internship topic are:
  - **contractile stromal / ECM / myofibroblastic**
  - **tumor ductal / epithelial interface**
  - **resident / perivascular / immune-border**

The most accurate final formulation is:

> The internship did not identify one single ion channel responsible for stromal immunomodulation in PDAC. Instead, it defined a structured landscape of KCN genes associated with distinct stromal states, tissue niches, and ECM / immunity-related programs, providing a strong basis for future functional validation.

---

## 13. Most logical experimental follow-ups

1. **stromal functional validation**
   - `KCNMA1`, `KCNMB1`, `KCNE4`, `KCND2`
   - contractility, ECM, migration, secretion

2. **tumor-stroma interface validation**
   - `KCNK1`, `KCNN4`, `KCNK6`
   - ductal plasticity, secretory programs, interface behavior

3. **resident / immune-border validation**
   - `KCNJ8`, `KCNK3`, `KCND3`, `KCNAB2`
   - inflammatory gradients, complement, vascular-border context

4. **immune-oriented assays**
   - CAF + immune co-culture
   - migration / infiltration assays after channel perturbation
   - cytokine and ECM readouts after channel modulation

---

## Main files supporting this synthesis

- `code/moffitt_stroma_exploration_pipeline.R`
- `code/moffitt_stroma_kcn_deg3_pipeline.R`
- `code/moffitt_stroma_kcn_add_test4_mast.R`
- `code/moffitt_stroma_gsva_activated_vs_normal.R`
- `code/mycaf_icaf_correlations/run_kcn_mycaf_icaf_correlations.R`
- `code/user_moffitt_validation/run_user_moffitt_validation.R`
- `code/bulk_survival_analysis/build_fused_survival_bundle.R`
- `code/bulk_survival_analysis/run_kcn25_multicohort_adjusted_cox_survival.R`
- `code/bulk_survival_analysis/run_kcn25_multicohort_km_scan_min15_survival.R`
- `code/bulk_survival_analysis/run_kcn25_multicohort_time_varying_cox_survival.R`
- `code/bulk/BKCA_shBKCA/run_bkca_shbkca_bulk_rnaseq_analysis.R`
- `code/bulk/SIGMAR1_shSIGMAR1/run_sigmar1_shsigmar1_bulk_rnaseq_analysis.R`
- `code/spatial_pdac_analysis/kcn/fisher_niche_enrichment_primary/run_kcn_niche_fisher_tests.R`
- `code/spatial_pdac_analysis/kcn/wilcoxon_ecotype_expression_contrasts_primary/run_kcn_ecotype_expression_contrasts_primary.R`
- `code/spatial_pdac_analysis/marker_ecotype_wilcoxon_primary/run_marker_ecotype_expression_contrasts_primary.R`
- `code/spatial_pdac_analysis/communication/hotspot/all_kcns_primary/run_hotspot_all_kcns_primary.py`
- `code/spatial_pdac_analysis/sparkx/run_sparkx_all_kcns_primary.R`
- `code/spatial_pdac_analysis/communication/hotspot/all_kcns_primary/local_correlations/run_kcn_local_correlations_primary.py`
- `code/spatial_pdac_analysis/communication/hotspot/all_kcns_primary/territory_gsea/run_all_kcn_territory_gsea.R`
- `code/spatial_pdac_analysis/cell_type_neighborhood/run_kcn_high_first_type_knn_neighborhood.py`
- `code/spatial_pdac_analysis/communication/hotspot/kcnma1_pathway_overlay_primary/run_kcnma1_pathway_overlay_primary.R`
- `code/spatial_pdac_analysis/communication/hotspot/kcnn4_pathway_overlay_primary/run_kcnn4_pathway_overlay_primary.R`
- `code/spatial_pdac_analysis/communication/hotspot/kcnd3_pathway_overlay_primary/run_kcnd3_pathway_overlay_primary.R`
- `code/spatial_pdac_analysis/communication/hotspot/marker_panel_primary/run_hotspot_marker_panel_primary.py`
- `code/spatial_pdac_analysis/communication/hotspot/sigmar1_primary/run_hotspot_sigmar1_primary.py`
- `code/spatial_pdac_analysis/communication/hotspot/kcnma1_subunits/run_hotspot_kcnma1_subunits_primary.py`
- `code/spatial_pdac_analysis/communication/hotspot/all_kcns_normal_pancreas/run_hotspot_all_kcns_normal_pancreas.py`
- `code/spatial_pdac_analysis/communication/hotspot/all_kcns_normal_pancreas/local_correlations/run_kcn_local_correlations_normal_pancreas.py`
- `code/spatial_pdac_analysis/communication/hotspot/all_kcns_normal_pancreas/kcn_territory_gsea/run_all_kcn_territory_gsea_normal_pancreas.R`

Last updated: 2026-05-24
