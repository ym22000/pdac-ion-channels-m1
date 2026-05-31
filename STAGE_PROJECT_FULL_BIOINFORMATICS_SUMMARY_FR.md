# Identification et analyse bio-informatique des canaux ioniques susceptibles de moduler l'immunité et la matrice extracellulaire dans le stroma du PDAC

Document préparé dans le cadre du projet de stage réalisé par **MISSOUM Youssra**, étudiante en **Master 1 Bio-informatique et Biologie Computationnelle** à l'Université Côte d'Azur, Nice.
Le travail a été mené à l'**Institut de Biologie Valrose (IBV)** dans l'équipe **Regulation of ion channel in cancer** dirigée par Olivier Soriani.
Page de l'équipe : [http://ibv.unice.fr/research-team/soriani/](http://ibv.unice.fr/research-team/soriani/)

---

## Résumé exécutif

L'objectif principal du projet était :

> **d'identifier puis de caractériser, par bio-informatique, des canaux ioniques et en particulier des gènes KCN susceptibles de moduler l'immunité et/ou la matrice extracellulaire dans l'adénocarcinome canalaire pancréatique à travers des territoires stromaux**

Le projet n'a pas identifié un unique "canal immunitaire stromal" universel. Il a plutôt défini un paysage KCN structuré à travers plusieurs couches complémentaires :

1. **single-cell RNA-seq** pour définir un panel stromal de gènes KCN et replacer ces gènes dans des états fibroblastiques
2. **survie bulk multicohorte** pour tester si les gènes KCN sélectionnés ont aussi une pertinence clinique
3. **jeux bulk locaux de perturbation** pour reconnecter certains candidats à des contextes expérimentaux du laboratoire
4. **spatial transcriptomics** pour déterminer où tombent réellement ces gènes KCN dans le tissu, avec quels programmes, et dans quelles niches

Le message intégré le plus fort est que le panel KCN se répartit en au moins trois grands mondes biologiques :

1. **interface tumorale ductale / épithéliale**
   - `KCNK1`, `KCNN4`, `KCNK6`
2. **monde stromal contractile / ECM / myofibroblastique**
   - `KCNMA1`, `KCNMB1`, `KCNE4`, `KCND2`
3. **monde résident / périvasculaire / immune-border**
   - `KCNJ8`, `KCNK3`, `KCND3`, `KCNAB2`, `KCNN3`

Ce point est important pour la question du stage car l'immunité stromale dans le PDAC est souvent modulée indirectement par :

- la contractilité fibroblastique
- le dépôt et le remodelage de la matrice extracellulaire
- les interfaces tumeur-stroma
- les territoires vasculaires et périvasculaires
- les zones de bordure inflammatoires et riches en complément

Autrement dit, les canaux ioniques les plus pertinents dans ce projet ne sont pas seulement des "gènes immunitaires". Ce sont aussi des régulateurs ou marqueurs potentiels de niches stromales qui organisent la matrice, la signalisation, l'exclusion et l'architecture locale du tissu.

---

## 1. Question initiale et justification biologique

## 1.1 Question centrale

Le projet partait d'une question plus large qu'une simple analyse d'expression différentielle :

- quels canaux ioniques, en particulier les gènes `KCN`, sont associés au compartiment stromal du PDAC ?
- ces canaux sont-ils liés à des états stromaux spécifiques ?
- sont-ils connectés à des programmes compatibles avec une modulation :
  - de l'immunité
  - de la matrice extracellulaire
  - de l'activation fibroblastique
  - des interfaces tumeur-stroma ?

## 1.2 Hypothèse biologique

Le projet repose sur une hypothèse biologiquement plausible :

- les canaux ioniques régulent le potentiel membranaire, le volume cellulaire, et les flux ioniques tels que `K+`, `Ca2+`, `Na+` et `Cl-`
- ces propriétés peuvent influencer :
  - la contractilité fibroblastique
  - la migration et l'adhérence
  - la sécrétion de cytokines
  - la production de matrice extracellulaire
  - l'organisation spatiale des territoires stromaux
- dans le PDAC, de tels effets peuvent indirectement moduler la distribution et la fonction des cellules immunes

L'objectif était donc de transformer une intuition mécanistique plausible en une priorisation bio-informatique structurée.

---

## 2. Jeux de données utilisés

## 2.1 Vue d'ensemble

| Bloc analytique | Jeu de données | Source | Rôle dans le projet |
|---|---|---|---|
| Atlas stromal scRNA-seq global | `data/Stroma_Subset2021.rds` | Oh et al., *Nature Communications* 2023 | exploration stromale globale |
| Atlas focalisé fibroblastes | `data/Stroma_Subset2021_fibroblast_focus.rds` | sous-ensemble curé du même atlas | priorisation des KCN et interprétation par sous-types |
| Validation externe fibroblastique | USER local inputs | USER / humanpdac | validation indépendante du cadre des états stromaux |
| Spatial transcriptomics | `data/PDAC_Updated_ST.rds` | Khaliq et al., *Nature Genetics* 2024 | localisation tissulaire et interprétation des niches |
| Bulk local de perturbation 1 | `code/bulk/BKCA_shBKCA` | laboratoire IBV | contexte `BKCa / KCNMA1` |
| Bulk local de perturbation 2 | `code/bulk/SIGMAR1_shSIGMAR1` | laboratoire IBV | contexte `SIGMAR1` |
| Cohortes bulk de survie | `code/bulk_survival_analysis/` | GSE183795, GSE62452, GSE28735, Bailey 2016 | pertinence clinique |

## 2.2 Pourquoi cette combinaison multimodale est importante

Chaque couche répond à une question différente :

- le **single-cell RNA-seq** dit quels états fibroblastiques sont impliqués
- la **survie bulk** dit si certains gènes ont un signal clinique à l'échelle patient
- les **bulk de perturbation** reconnectent les signaux bio-informatiques à des axes expérimentaux du laboratoire
- le **spatial transcriptomics** dit où tombent ces gènes dans le tissu et avec quels programmes ils se co-organisent spatialement

La force du projet vient de la convergence de ces couches plutôt que d'un test isolé.

---

## 3. Chronologie analytique globale

Le projet s'est construit dans l'ordre suivant :

1. exploration stromale single-cell
2. priorisation KCN centrée sur les fibroblastes
3. interprétation par états et sous-types en scRNA-seq
4. interprétation orientée `myCAF` / `iCAF`
5. validation externe dans USER
6. survie bulk multicohorte sur le panel final KCN
7. analyses bulk locales de perturbation issues du laboratoire
8. spatial transcriptomics :
   - cartes
   - écotypes
   - tests de niche
   - Hotspot
   - SPARK-X
   - local correlations
   - territory GSEA
   - analyses de voisinage
   - overlays de pathways
   - comparaison avec le pancréas normal
9. synthèse intégrée

Ce document suit cette chronologie.

---

## 4. Single-cell RNA-seq : priorisation stromale et fibroblastique des KCN

## 4.1 Pourquoi commencer par le single-cell RNA-seq ?

C'était le bon point d'entrée car la question biologique centrale était stromale. Un point de départ bulk aurait mélangé :

- signal tumoral épithélial
- signal stromal
- signal immunitaire
- signal vasculaire

L'atlas single-cell permet de :

- isoler le monde stromal
- se focaliser sur les états fibroblastiques
- définir un panel KCN biologiquement cohérent avant d'aller vers les couches spatiales et cliniques

## 4.2 Atlas stromal global utilisé d'abord

Objet source :

- `data/Stroma_Subset2021.rds`

Métriques globales documentées :

- `21469` cellules
- `16972` gènes
- `40` échantillons
- `3` datasets
- `8` labels `CellType2`
- `77.32%` `Activated`
- `22.68%` `Normal`

Principaux datasets représentés dans l'objet stromal :

- `Powers`
- `Peng`
- `Qadir`

Ce premier objet a servi à comprendre le paysage stromal global avant de restreindre l'analyse à des questions fibroblastiques.

## 4.3 Objet focalisé fibroblastes utilisé pour la priorisation KCN

Objet source :

- `data/Stroma_Subset2021_fibroblast_focus.rds`

Les principaux états fibroblastiques considérés dans le projet étaient :

- `qPSC`
- `smPSC`
- `myCAF`
- `csCAF`
- `iCAF`

Effectifs documentés des sous-types dans les analyses AUCell / fibroblastiques :

- `qPSC = 5743`
- `smPSC = 5521`
- `myCAF = 4809`
- `csCAF = 4242`
- `iCAF = 338`

Cette distribution contient déjà une limite importante :

- le compartiment `iCAF` est beaucoup plus petit que les autres
- donc toutes les conclusions sur `iCAF` doivent être lues avec davantage de prudence

## 4.4 États stromaux Activated vs Normal

Le projet a utilisé le cadre stromal de Moffitt pour définir une polarité globale simple :

`StromaScore = Moffitt.F5_ActivatedStroma.top25 - Moffitt.F13_NormalStroma.top25`

Règle :

- `Activated` si `F5 > F13`
- `Normal` sinon

Cette étape n'a pas encore sélectionné les gènes KCN. Elle a établi une colonne vertébrale biologiquement interprétable :

- un monde stromal global activé
- un monde stromal global normal-like

Ce point était utile car, ensuite, les gènes KCN pouvaient être interprétés non seulement comme des marqueurs de sous-type mais aussi comme des gènes liés à une polarité stromale globale.

## 4.5 Espace de recherche KCN et préfiltrage

L'espace initial des KCN a été défini à partir des symboles géniques correspondant au motif `^KCN`.

Ensuite, un préfiltrage permissif mais biologiquement motivé a été appliqué :

- détecté dans au moins `1%` des cellules de l'objet fibroblastique
- détecté à au moins `1%` dans au moins `2` sous-types stromaux

Pourquoi utiliser un seuil bas ?

- beaucoup de canaux ioniques restent biologiquement intéressants malgré une détection clairsemée
- un seuil trop strict aurait éliminé des candidats de niche

Pourquoi ne pas garder tous les KCN sans filtrage ?

- pour réduire l'instabilité provenant de gènes supportés par très peu de cellules
- pour stabiliser les modèles mixtes et les tests pseudobulk en aval

## 4.6 Pourquoi quatre stratégies DEG ont été utilisées

Le choix méthodologique central du bloc single-cell a été :

> **ne pas faire confiance à une seule méthode DEG pour des gènes KCN clairsemés**

Ce point est particulièrement important pour les canaux ioniques car ils sont souvent :

- faiblement détectés
- riches en zéros
- hétérogènes selon les patients
- informatifs soit par leur fréquence de détection, soit par leur intensité chez les positifs

Le panel final KCN a donc été construit par convergence entre quatre stratégies analytiques différentes.

## 4.7 Test 1 : pseudobulk DESeq2 Activated vs Normal

### Objectif

Identifier les gènes KCN qui distinguent l'état stromal global `Activated` de l'état stromal `Normal`.

### Unité d'analyse

Pseudobulk agrégé par :

- `Patient`
- `StromaState`

Règles de rétention :

- `n_cells >= 20` par pseudo-échantillon
- conserver uniquement les patients représentés dans les deux états

### Design DESeq2

`~ Patient + n_cells + condition`

où :

- `condition = Activated / Normal`

### Paramètres DESeq2

- gènes gardés si `rowSums(counts) >= 10`
- `estimateSizeFactors(type = "poscounts")`
- `DESeq(sfType = "poscounts", fitType = "local")`

Coefficient testé :

- `condition_Activated_vs_Normal`

Seuil de décision :

- `BH FDR < 0.05`

### Pourquoi ce modèle est utile

Ce cadre pseudobulk :

- réduit la pseudoréplication au niveau cellule
- traite le patient comme unité biologique
- fournit un contraste stromal global robuste

### Profil principal de résultats

Exemples de gènes enrichis en `Activated` :

- `KCNK6`
- `KCNT2`
- `KCND3`
- `KCNN3`
- `KCNJ8`

Exemples de gènes enrichis en `Normal` :

- `KCNAB1`
- `KCNA5`
- `KCNMA1`
- `KCNK17`
- `KCNC4`
- `KCNQ4`
- `KCNS3`
- `KCNE4`
- `KCNB1`

Nombre documenté de KCN significatifs dans ce test :

- `14`

### Limites

- perte de résolution cellulaire
- un signal réel mais porté par un sous-cluster rare peut être dilué
- c'est une analyse d'état global plutôt qu'une analyse de niche

## 4.8 Test 2 : pseudobulk DESeq2 sous-type vs reste

### Objectif

Identifier les gènes KCN enrichis dans un sous-type fibroblastique spécifique.

### Logique

Pour chaque sous-type, construire un contraste :

- `target` vs `rest`

### Design

`~ Patient + n_cells + condition`

où :

- `condition = target / rest`

### Seuils de décision

- `BH FDR < 0.05`
- `|log2FC| >= 0.25`

### Pourquoi ce test est important

Un gène peut être très informatif pour un état fibroblastique précis sans être un grand marqueur global `Activated vs Normal`.

Ce test était donc nécessaire pour capturer des motifs KCN spécifiques d'un sous-type.

### Profil global des résultats

Nombres approximatifs documentés de KCN significatifs par sous-type :

- `qPSC`: `8`
- `smPSC`: `9`
- `myCAF`: `7`
- `csCAF`: `1`
- `iCAF`: presque aucun signal robuste

Exemples :

- `qPSC`: `KCNJ8`, `KCNK17`, `KCNE4`, `KCNK6`
- `smPSC`: `KCNA5`, `KCNAB1`, `KCNMB1`, `KCNMA1`
- `myCAF`: `KCNMA1`, `KCNMB1`, `KCND3`, `KCNK3`
- `csCAF`: surtout `KCNK17`

### Limites

- le groupe `rest` est biologiquement hétérogène
- des sous-types proches peuvent partiellement se masquer entre eux

## 4.9 Test 3 : modèle hurdle mixte

### Objectif

Modéliser plus correctement les gènes KCN clairsemés en séparant :

- **la détection**
- **l'intensité d'expression chez les cellules positives**

### Modèle de détection

`glmer(detected ~ is_target + log10_nCount_RNA + Dataset + (1 | Patient), family = binomial)`

### Modèle d'expression positive

`lmer(expr_pos ~ is_target + log10_nCount_RNA + Dataset + (1 | Patient))`

### Règles d'éligibilité

Partie détection :

- au moins `20` cellules positives au total
- au moins `3` patients positifs dans target et rest
- présence des deux classes de détection

Partie expression positive :

- au moins `20` cellules positives au total
- au moins `10` cellules positives dans target
- au moins `10` cellules positives dans rest
- au moins `3` patients positifs dans target et rest
- plus d'une valeur positive distincte

### Combinaison

- combinaison des p-values par la méthode de Fisher
- support classé en :
  - `both_parts`
  - `detection_only`
  - `positive_only`
  - `none`

### Pourquoi ce modèle est particulièrement pertinent pour les KCN

Certains canaux sont informatifs parce que :

- ils sont détectés dans davantage de cellules
- ou parce que, parmi les cellules positives, ils sont plus élevés

Le cadre hurdle capture ces deux aspects de la biologie.

### Résultats documentés

Nombre de gènes KCN avec support par sous-type :

- `qPSC`: `18`
- `smPSC`: `17`
- `myCAF`: `16`
- `csCAF`: `16`
- `iCAF`: `10`

Ce test a été l'un des plus utiles pour les gènes KCN rares et de type niche.

### Limites

- modèle plus complexe
- moins intuitif qu'un simple log2FC
- dépend encore d'un bon support de sous-groupe

## 4.10 Test 4 : MAST

### Objectif

Ajouter une méthode DEG single-cell standard adaptée aux expressions riches en zéros.

### Paramètres principaux

Utilisation de `FindMarkers(..., test.use = "MAST", latent.vars = "nCount_RNA")`

Seuil de décision :

- `padj_BH < 0.05`

Direction interprétée avec :

- `avg_log2FC >= 0.25`

### Pourquoi MAST était utile

MAST fournit une perspective DEG single-cell bien établie, complémentaire du pseudobulk et des modèles mixtes.

### Résultats documentés

Nombres approximatifs de KCN significatifs par sous-type :

- `qPSC`: `19`
- `smPSC`: `18`
- `myCAF`: `18`
- `csCAF`: `21`
- `iCAF`: `10`

Exemples :

- `qPSC`: `KCNJ8`, `KCNK17`, `KCNE4`
- `smPSC`: `KCNA5`, `KCNAB1`, `KCNMB1`, `KCNMA1`
- `myCAF`: `KCNK17`, `KCNJ8`, `KCNA5`, `KCNE4`, `KCNMB1`, `KCND2`
- `iCAF`: plus faible, mais `KCNJ8`, `KCNS3`, `KCNK1`, `KCNK6` figuraient parmi les plus visibles

### Limites

- reste sensible à la structure single-cell
- le patient est moins explicitement modélisé que dans les pseudobulks

## 4.11 Panel final : 25 gènes KCN

### Règle d'inclusion

Un gène entre dans le panel final s'il est retenu par au moins une des quatre stratégies DEG :

- Test 1 : significatif après BH
- Test 2 : significatif après BH et `|log2FC| >= 0.25`
- Test 3 : `support_pattern != none`
- Test 4 : significatif après BH

### Panel final

L'union finale contient **25 gènes KCN** :

`KCNA5`, `KCNA7`, `KCNAB1`, `KCNAB2`, `KCNB1`, `KCNC4`, `KCND2`, `KCND3`, `KCNE4`, `KCNJ15`, `KCNJ16`, `KCNJ8`, `KCNK1`, `KCNK15`, `KCNK17`, `KCNK3`, `KCNK6`, `KCNMA1`, `KCNMB1`, `KCNMB4`, `KCNN3`, `KCNN4`, `KCNQ4`, `KCNS3`, `KCNT2`

### Tableau de robustesse inter-méthodes

| Nombre de méthodes supportant le gène | Gènes |
|---|---|
| `4/4` | `KCNA5`, `KCNAB1`, `KCNC4`, `KCND3`, `KCNE4`, `KCNJ8`, `KCNK17`, `KCNK6`, `KCNMA1` |
| `3/4` | `KCNK3`, `KCNMB1`, `KCNN3`, `KCNS3`, `KCNT2` |
| `2/4` | `KCNA7`, `KCND2`, `KCNJ15`, `KCNK1`, `KCNK15`, `KCNMB4`, `KCNN4`, `KCNQ4` |
| `1/4` | `KCNAB2`, `KCNB1`, `KCNJ16` |

### Comment lire ce panel

- `4/4`: candidats très robustes
- `3/4`: solides mais plus dépendants du cadre méthodologique
- `2/4`: biologiquement intéressants, robustesse statistique modérée
- `1/4`: conservés par ouverture biologique, à interpréter avec plus de prudence

Cette union finale a été cruciale car elle a défini l'espace de candidats partagé pour toutes les analyses suivantes.

## 4.12 Polarité des 25 KCN dans l'axe Normal vs Activated

Table principale :

- `code/moffitt_kcn_deg_4_methods/tables/KCN_normal_vs_activated_poles_membership.tsv`

Catégories documentées :

| Catégorie | Gènes |
|---|---|
| `Activated_only` | `KCNAB2`, `KCND2`, `KCND3`, `KCNJ15`, `KCNK1`, `KCNK15`, `KCNMB4`, `KCNN3`, `KCNN4` |
| `Normal_only` | `KCNA5`, `KCNAB1`, `KCNB1`, `KCNE4`, `KCNJ16`, `KCNK17`, `KCNK3`, `KCNQ4`, `KCNS3` |
| `Both` | `KCNA7`, `KCNC4`, `KCNJ8`, `KCNK6`, `KCNMA1`, `KCNMB1`, `KCNT2` |

Ce résultat est conceptuellement important car il montre très tôt que :

- le panel KCN n'est pas une signature unique de stroma activé
- certains gènes ressemblent plutôt à un monde résident ou normal-like
- certains gènes se situent dans des axes transitionnels ou dépendants du contexte

## 4.13 GSVA Activated vs Normal : pourquoi cette étape a été intégrée

### Objectif

GSVA n'a pas été utilisé pour sélectionner les gènes KCN. Il a été utilisé pour décrire le paysage biologique des états stromaux dans lesquels les gènes KCN s'inscrivent.

### Choix principaux

- pseudobulk par `Patient x StromaState`
- transformation :
  - `log_cpm = log1p((counts / library_size) * 1e6)`
- collection Hallmark depuis `msigdbr`
- garder les pathways avec `10-500` gènes en recouvrement
- `kcdf = "Gaussian"`
- Wilcoxon apparié par patient
- correction BH

### Pourquoi c'était utile

Cette étape répond à :

> quel type de monde biologique représente le stroma activé ?

C'est directement pertinent pour le stage car le sujet est immunité plus ECM, pas seulement l'expression de canaux ioniques.

### Profil principal des résultats

Pathways plus hauts dans le stroma `Activated` :

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

Pathways plus hauts dans le stroma `Normal` :

- `OXIDATIVE_PHOSPHORYLATION`
- `ADIPOGENESIS`
- `MYC_TARGETS_V1`
- `DNA_REPAIR`
- `MYOGENESIS`
- `FATTY_ACID_METABOLISM`

### Interprétation

Avant même de discuter chaque gène KCN, le stroma activé est déjà fortement compatible avec :

- inflammation
- signalisation interféron
- complément
- angiogenèse
- coagulation
- remodelage de type EMT

Cela définit le terrain biologique dans lequel les candidats KCN seront interprétés ensuite.

### Pourquoi la GSEA n'a pas été l'outil principal de sélection

Parce que l'objectif du bloc single-cell était d'abord d'identifier de manière fiable des gènes KCN clairsemés. Pour cela, des cadres DEG et des scores cellule par cellule étaient plus adaptés.

Les méthodes de type GSEA ont été plus utiles plus tard pour :

- l'interprétation orientée pathways
- la lecture de territoires
- la contextualisation de gènes KCN déjà sélectionnés

## 4.14 AUCell : pourquoi AUCell plutôt qu'une GSEA classique à ce stade

### Logique centrale

AUCell donne un **score par cellule** pour une signature. C'est beaucoup mieux adapté qu'une GSEA bulk classique lorsqu'on demande :

> à l'intérieur d'un sous-type fibroblastique, les cellules `KCN-high` sont-elles fonctionnellement différentes des cellules `KCN-low` ?

### Définition de High vs Low

- si la médiane d'expression = `0` :
  - `High = expression > 0`
  - `Low = expression == 0`
- sinon :
  - `High = expression > median`
  - `Low = expression <= median`

Conditions minimales :

- au moins `5` cellules dans chaque groupe par patient
- au moins `20` cellules totales dans chaque groupe
- au moins `3` patients appariés

### Paramètres AUCell

- `AUCell_buildRankings(..., plotStats = FALSE)`
- `AUCell_calcAUC(..., aucMaxRank = ceiling(n_genes * 0.05))`

### Ce que cela a apporté

AUCell posait une question plus fine :

> dans un sous-type stromal donné, quel programme biologique caractérise les cellules `KCN-high` ?

### Exemples issus des résultats documentés

`KCNE4` dans `csCAF` :

- augmentation de `EMT`
- augmentation de `TGF_BETA_SIGNALING`
- augmentation de `ANGIOGENESIS`
- augmentation de `COMPLEMENT`
- augmentation de `HYPOXIA`
- augmentation de `MYOGENESIS`

`KCNJ8` dans `csCAF` :

- augmentation de `MTORC1_SIGNALING`
- augmentation de `ANGIOGENESIS`
- augmentation de `OXIDATIVE_PHOSPHORYLATION`
- augmentation de `MYC_TARGETS_V1`
- augmentation de `COMPLEMENT`

`KCNK6` dans `csCAF` :

- augmentation de `MYOGENESIS`
- augmentation de `CHOLESTEROL_HOMEOSTASIS`
- augmentation de `OXIDATIVE_PHOSPHORYLATION`
- augmentation de `EMT`

### Interprétation

Cela suggère que certains gènes KCN ne sont pas seulement des marqueurs de sous-types, mais des marqueurs de **sous-états fonctionnels à l'intérieur d'un sous-type**.

## 4.15 Analyse de corrélation `myCAF` / `iCAF`

### Objectif

Replacer les 25 gènes KCN par rapport à deux axes stromaux majeurs pour le projet :

- `myCAF` : contractile, riche en ECM, orienté exclusion
- `iCAF` : inflammatoire et sécrétoire

### Inputs et méthode

Script principal :

- `code/mycaf_icaf_correlations/run_kcn_mycaf_icaf_correlations.R`

Inputs :

- `kcn_union_25.tsv`
- `genes_top30_myCAF_human.txt`
- `genes_top30_iCAF_human.txt`

Métrique principale :

- `Spearman rho` entre chaque gène KCN et chaque gène de signature

Métriques résumées :

- `mean_rho`
- `median_rho`
- `mean_abs_rho`
- `min_p_adj_bh`
- `n_fdr_lt_0_05`

### Principaux signaux significatifs

Exemples représentatifs :

| Gène | Programme | FDR |
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

### Nuance très importante sur `KCNMA1`

`KCNMA1` devient plus tard l'un des gènes spatiaux les plus intéressants, mais dans ce bloc de corrélation pseudobulk single-cell, ce n'est **pas** un grand signal global `myCAF` :

- `median rho ~ 0.021`
- seulement `3/30` corrélations significatives avec les marqueurs top30 `myCAF`

C'est une leçon clé de l'ensemble du projet :

> un gène peut être biologiquement très pertinent comme marqueur de territoire spatial tout en restant modeste dans un cadre de corrélation pseudobulk de signature

### Interprétation

Ce bloc a été utile car il séparait déjà :

- un groupe plus contractile / ECM-like
- un groupe plus inflammatoire / résident / border-like

Mais il ne résolvait pas complètement la géographie tissulaire. C'est pourquoi le projet a ensuite eu besoin du spatial transcriptomics.

## 4.16 Validation externe dans USER

### Objectif

Éviter de dépendre d'un seul atlas stromal.

### Logique principale

Au lieu de transférer aveuglément les labels, les fibroblastes USER ont été :

1. reconstruits à partir des inputs locaux
2. filtrés pour ne garder que les fibroblastes
3. re-clusterisés
4. rescorrés avec des signatures stromales de type Moffitt
5. annotés à partir :
   - du score moyen de signature
   - du recouvrement entre les DEGs du cluster et les signatures

États finaux retenus :

- `myCAF`
- `iCAF`
- `csCAF`
- `PSC`

### Effectifs finaux des fibroblastes USER

- `myCAF = 11308` (`79.26%`)
- `csCAF = 1326` (`9.29%`)
- `PSC = 1281` (`8.98%`)
- `iCAF = 352` (`2.47%`)

### Contrôles de cohérence marqueurs

- `myCAF`: `ACTA2` détecté dans `44.31%`, `FAP` dans `58.57%`
- `csCAF`: `C7` dans `68.55%`
- `iCAF`: `IL11` dans `5.68%`, `FAP` dans `72.44%`

### Pourquoi c'est important

Cette validation ne prouve pas chaque gène KCN mécaniquement, mais elle renforce le cadre des états stromaux qui sous-tend le reste du projet.

## 4.17 Principales conclusions du bloc single-cell

Le bloc single-cell n'identifie pas un programme stromal KCN homogène. Il suggère déjà plusieurs mondes biologiques :

- un monde ductal / tumor-like
- un monde contractile / ECM / myofibroblastique
- un monde résident / border / périvasculaire / inflammatoire

Surtout, il a produit le **panel de 25 gènes KCN** qui a servi de colonne vertébrale au reste du projet.

## 4.18 Principales limites du bloc single-cell

- transcriptomique uniquement, pas fonctionnel
- `Activated vs Normal` est une bonne colonne vertébrale, pas une vérité absolue
- les contrastes sous-type vs reste sont sensibles à l'hétérogénéité du groupe reste
- `iCAF` est sous-puissant
- les gènes KCN clairsemés restent sensibles au modèle
- certains gènes, en particulier `KCNMA1`, sont imparfaitement capturés par les cadres pseudobulk de sous-types

---

## 5. Survie bulk multicohorte : pertinence clinique des 25 gènes KCN

## 5.1 Pourquoi placer la survie juste après le scRNA-seq ?

À ce stade du projet, le panel KCN de 25 gènes a été priorisé biologiquement dans les fibroblastes. La question logique suivante est :

> ces candidats ont-ils aussi une signification clinique à l'échelle patient ?

C'est pourquoi le bloc de survie multicohorte vient naturellement **après** la sélection single-cell et **avant** les bulk de perturbation locaux et les analyses spatiales de niche.

## 5.2 Pourquoi la survie est utile mais insuffisante

La survie bulk ne dit **pas** :

- quel type cellulaire exprime le gène
- quelle niche porte le signal
- si l'effet est stromal, tumoral, vasculaire ou immunitaire

Ce qu'elle apporte en revanche :

- une priorisation clinique
- une couche de pertinence à l'échelle patient
- une façon de distinguer des candidats avec potentiel translationnel

## 5.3 Cohortes utilisées

Le bundle intégré de survie combine quatre cohortes publiques PDAC :

- `PH_GSE183795`
- `Hussain_GSE62452`
- `Zhang_GSE28735`
- `Bailey_QCMG_UQ_2016`

Effectifs documentés de patients éligibles à la survie :

- `Bailey`: `96` éligibles, `64` événements, médiane d'OS environ `14.98 mois`
- `Hussain`: `65` éligibles, `49` événements, médiane d'OS environ `14.9 mois`
- `PH`: `134` éligibles, `86` événements, médiane d'OS environ `11.31 mois`
- `Zhang`: `42` éligibles, `29` événements, médiane d'OS environ `14.5 mois`

Total dans l'analyse poolée :

- `337` patients
- `228` événements

## 5.4 Pourquoi un bundle de survie fusionné a été construit

Les cohortes proviennent de plateformes et historiques de preprocessing différents :

- microarray
- RNA-seq

L'analyse n'a donc **pas** concaténé directement les intensités brutes.

À la place, le pipeline :

1. a chargé les quatre bundles cohortes
2. a retenu les échantillons tumoraux avec survie valide
3. a défini l'espace commun des gènes KCN entre cohortes
4. a standardisé chaque gène **à l'intérieur de chaque cohorte**
5. a poolé les échantillons ensuite

Cela préserve le classement intra-cohorte tout en réduisant les effets de plateforme.

## 5.5 Espace génique analysé

Le bundle fusionné est construit sur un espace commun de `93` gènes KCN.

Cependant, les analyses de survie elles-mêmes portent sur :

- le **panel final de 25 gènes KCN** issu du bloc single-cell

Ce point est conceptuellement important :

- le bloc survie ne redécouvre pas de nouveaux gènes KCN à partir de zéro
- il teste le poids clinique de candidats déjà curés biologiquement

## 5.6 Harmonisation clinique

Les variables cliniques ont été harmonisées en facteurs simplifiés :

- `stage_group = I / II / III / IV / Unknown`
- `grade_group = Low_1_2 / High_3_4 / Unknown`
- `margin_group = R0 / R1 / Unknown`

Le modèle de Cox principal ajusté a retenu :

- `stage_group`
- `grade_group`
- `strata(cohort)`

L'information de marge a été harmonisée dans le bundle mais non retenue dans le modèle principal final.

## 5.7 Standardisation de l'expression

À l'intérieur de chaque cohorte, l'expression génique a été standardisée gène par gène :

`z = (x - mean_cohort) / sd_cohort`

Ainsi, chaque coefficient de survie s'interprète comme :

> effet sur le risque pour une augmentation de `+1 écart-type` d'expression du gène dans sa cohorte d'origine

C'est un point méthodologique clé car cela évite une comparaison naïve d'intensités absolues entre plateformes.

## 5.8 Modèle de Cox ajusté principal

Script principal :

- `code/bulk_survival_analysis/run_kcn25_multicohort_adjusted_cox_survival.R`

Modèle :

`Surv(os_months, event) ~ gene_value + stage_group + grade_group + strata(cohort)`

Interprétation :

- `gene_value` = expression z-score intra-cohorte
- `strata(cohort)` = risque de base spécifique à chaque cohorte

Ce n'est pas un modèle de Cox poolé naïf. Chaque cohorte conserve sa propre structure de risque de base.

## 5.9 Principales métriques dans les résultats de Cox ajusté

La table principale rapporte :

- `beta`
- `HR = exp(beta)`
- `hr_low`, `hr_high`
- `SE`
- `z`
- `p_value`
- `p_adj_bh`
- `concordance`

Interprétation :

- `HR > 1`: forte expression associée à un pronostic plus défavorable
- `HR < 1`: forte expression associée à un pronostic plus favorable
- `p_adj_bh`: significativité après correction multiple sur les 25 gènes

## 5.10 Gènes significatifs dans le modèle de Cox ajusté

Significatifs après correction BH :

| Gène | HR | IC95% | p-value | FDR | Interprétation |
|---|---:|---|---:|---:|---|
| `KCNMB4` | `1.283` | `[1.102 ; 1.494]` | `0.00130` | `0.0177` | défavorable si haut |
| `KCNC4` | `0.802` | `[0.699 ; 0.920]` | `0.00159` | `0.0177` | protecteur si haut |
| `KCNAB1` | `0.807` | `[0.703 ; 0.925]` | `0.00212` | `0.0177` | protecteur |
| `KCNE4` | `0.808` | `[0.699 ; 0.933]` | `0.00379` | `0.0237` | protecteur |
| `KCNK6` | `1.202` | `[1.057 ; 1.367]` | `0.00509` | `0.0254` | défavorable |
| `KCNT2` | `0.848` | `[0.744 ; 0.965]` | `0.0126` | `0.0478` | protecteur |
| `KCNN3` | `0.834` | `[0.723 ; 0.963]` | `0.0134` | `0.0478` | protecteur |

Exemples de tendances intéressantes non significatives :

- `KCNMA1`: `HR = 0.892`, `p = 0.077`, `FDR = 0.241`
- `KCNK17`: `HR = 0.890`, `p = 0.106`
- `KCNK1`: `HR = 0.899`, `p = 0.107`
- `KCNK3`: `HR = 0.902`, `p = 0.158`

## 5.11 Interprétation biologique principale des résultats de survie

La couche bulk survie ne reflète pas parfaitement la couche spatiale.

Par exemple :

- `KCNMA1` devient ensuite l'un des candidats stromaux spatiaux les plus convaincants
- mais n'est qu'une tendance dans le modèle de Cox ajusté

À l'inverse :

- `KCNMB4` est le signal le plus fort en survie bulk
- alors que spatialement il est beaucoup moins dominant que `KCNMA1`

Il s'agit d'un résultat conceptuellement très important :

> la pertinence clinique bulk et l'identité de niche tissulaire ne sont pas la même chose

## 5.12 Vérification de l'hypothèse de risques proportionnels

Table complémentaire :

- `code/bulk_survival_analysis/survival_kcn25_venn_adjusted_cox/tables/kcn25_multicohort_cox_zph_checks.tsv`

Principaux moffitt_stromal_exploration_outputs :

- `gene_ph_p`
- `global_ph_p`
- `gene_ph_fdr`
- `global_ph_fdr`

Plusieurs gènes violent l'hypothèse de risques proportionnels, notamment :

- `KCNS3`
- `KCNJ8`
- `KCNB1`
- `KCNK1`
- `KCNK15`
- `KCNA7`
- `KCNQ4`

C'est précisément pour cela qu'une extension en Cox time-varying a été ajoutée.

## 5.13 Scan exploratoire de cutoffs en Kaplan-Meier

Script principal :

- `code/bulk_survival_analysis/run_kcn25_multicohort_km_scan_min15_survival.R`

Objectif :

- tester si certains gènes montrent une séparation de survie plus nette pour un certain cutoff `high vs low`

Contraintes :

- `min_group_n = 15`
- `min_events_per_group = 5`

Corrections :

- BH within gene sur les cutoffs
- BH across genes sur le meilleur cutoff par gène

### Résultats KM les plus forts

- `KCNMB4`: cutoff `z = 1.013`, `49 high vs 288 low`, `raw p = 5.56e-11`, `BH across genes = 1.39e-09`, médiane d'OS `8.88 vs 20.14 mois`
- `KCNAB1`: `raw p = 3.26e-06`, `BH across genes = 4.08e-05`
- `KCNA5`: `raw p = 9.55e-06`, `BH across genes = 7.96e-05`
- `KCNC4`: `raw p = 4.62e-05`, `BH across genes = 2.89e-04`
- `KCNE4`: `raw p = 1.14e-04`, `BH across genes = 5.69e-04`
- `KCNN3`: `raw p = 3.17e-04`, `BH across genes = 0.00132`
- `KCNK17`: `raw p = 4.29e-04`, `BH across genes = 0.00153`

Interprétation :

- c'est une couche utile pour l'exploration et la visualisation
- mais elle reste plus optimiste car elle recherche un cutoff favorable

## 5.14 Modèle de Cox time-varying

Script principal :

- `code/bulk_survival_analysis/run_kcn25_multicohort_time_varying_cox_survival.R`

Modèle :

`Surv(os_months, event) ~ gene_value + tt(gene_value) + stage_group + grade_group + strata(cohort)`

avec :

`tt(x, t) = x * (log(t) - mean(log(t)))`

Objectif :

- tester si l'effet du gène change au cours du suivi

Principaux moffitt_stromal_exploration_outputs :

- `beta_main`
- `beta_time_interaction`
- `lrt_p_value`
- `lrt_p_adj_bh`
- trajectoires de HR aux temps :
  - `6`
  - `12`
  - `24`
  - `36`
  - `48`
  - `60` mois

### Effets temporels significatifs après BH

- `KCNB1`: `LRT FDR = 0.00245`
- `KCNS3`: `LRT FDR = 0.00245`
- `KCNJ8`: `LRT FDR = 0.00271`
- `KCNK1`: `LRT FDR = 0.0148`
- `KCNK15`: `LRT FDR = 0.0148`
- `KCNA7`: `LRT FDR = 0.0440`

Exemples :

- `KCNJ8`: HR environ `1.56` à `6 mois`, puis environ `0.82` à `60 mois`
- `KCNB1`: HR environ `1.45` à `6 mois`, puis environ `0.83` à `60 mois`
- `KCNK1`: HR environ `0.60` à `6 mois`, puis environ `0.99` à `60 mois`

Interprétation :

- certains gènes KCN ne sont pas bien résumés par un HR constant unique
- leur association au risque peut s'atténuer, s'inverser ou s'aplatir au cours du temps

## 5.15 Ce que le bloc survie apporte au projet

Ce bloc ne dit pas où vit un gène KCN dans le tissu. Il dit :

- s'il a une pertinence à l'échelle patient
- si son effet est continu
- si une séparation visuelle `high vs low` existe
- si cet effet est stable dans le temps

C'est pour cela qu'il doit venir après la sélection single-cell mais avant le bloc spatial.

## 5.16 Principales limites du bloc survie

- le bulk mélange signaux tumoraux, stromaux, immunitaires et vasculaires
- le z-scoring intra-cohorte réduit les effets de plateforme sans les éliminer complètement
- le modèle principal n'ajuste que sur `stage` et `grade`
- plusieurs gènes violent l'hypothèse de risques proportionnels
- les scans KM sont exploratoires et peuvent amplifier artificiellement la séparation
- la survie est une couche de priorisation, pas une preuve mécanistique

---

## 6. Jeux bulk locaux de perturbation issus du laboratoire

## 6.1 Pourquoi garder ces analyses bulk locales après la survie multicohorte ?

La survie multicohorte apporte une pertinence clinique à l'échelle patient. Les jeux bulk locaux de perturbation répondent à une autre question :

> lorsqu'un axe biologique pertinent du laboratoire est perturbé, des programmes transcriptomiques cohérents émergent-ils autour du canal d'intérêt ?

Ces jeux bulk ne sont pas centraux pour la sélection initiale des 25 KCN, mais ils aident à reconnecter le projet au contexte biologique expérimental du laboratoire.

## 6.2 BKCa / shBKCa

### Jeu de données

- `8` échantillons
- `4 CTRL`
- `4 shBKCa`
- `17889` gènes collapseés

### Principaux résultats

- `6734` DEGs à `padj < 0.05`
- top signature CAF par Wilcoxon : `ifCAF_manual`
- corrélation positive la plus forte avec `KCNMA1` :
  - `FN1`, `rho = 0.952`, `FDR = 0.0348`
- corrélation négative la plus forte :
  - `MGST1`, `rho = -1`, `FDR = 0.0030`
- Hallmark principal :
  - `INTERFERON_ALPHA_RESPONSE`, `NES = -2.539`, `FDR = 3.9e-11`

### Interprétation

Ce bloc est important car il donne une cohérence externe à l'axe `KCNMA1 / BKCa` :

- lien fort à `FN1`
- signal orienté matrice
- implication de voies immuno-inflammatoires

Même si les données sont bulk et non purement stromales, elles soutiennent le cadre plus général canal-ECM-immunité.

## 6.3 SIGMAR1 / shSIGMAR1

### Jeu de données

- `6` échantillons
- `3 CTRL`
- `3 shSIGMAR1`
- `14729` gènes collapseés

### Principaux résultats

- `1546` DEGs à `padj < 0.05`
- top signature CAF : `myCAF_top30`, mais non significative après FDR
- corrélation positive la plus forte avec `SIGMAR1` :
  - `EZR`, `rho = 0.943`
- corrélation négative la plus forte :
  - `CXADR`, `rho = -1`
- Hallmark principal :
  - `G2M_CHECKPOINT`, `NES = -2.047`, `FDR = 1.8e-10`

### Interprétation

Ce bloc est moins directement lié au panel des 25 KCN, mais il reste utile car il connecte le projet à un contexte plus large de signalisation liée aux canaux ioniques, pertinent pour le laboratoire d'accueil.

## 6.4 Principales limites du bloc bulk local de perturbation

- effectifs très faibles
- un seul contexte expérimental par axe
- bulk sans purification stromale
- utile pour la cohérence, faible comme preuve définitive de priorisation

---

## 7. Spatial transcriptomics : localisation tissulaire, niches et programmes

## 7.1 Pourquoi la couche spatiale était nécessaire

Après le scRNA-seq et le bulk, la question manquante était :

> où tombent réellement ces gènes KCN dans le tissu PDAC ?

Cette question est essentielle car, dans le PDAC :

- l'immunité dépend de l'architecture tissulaire
- l'ECM est spatialement organisée
- les interfaces comptent
- les territoires périvasculaires et de bordure comptent

La couche spatiale est donc celle où le projet passe des **gènes candidats** aux **niches candidates**.

## 7.2 Objet spatial et métriques globales

Objet principal :

- `data/PDAC_Updated_ST.rds`

Métriques globales documentées :

- `91,496` spots
- `30` sections
- `8` assays
- assay par défaut : `rctd_fullfinal`

Distribution par origine :

- `Pancreas`: `35,458` spots (`38.75%`)
- `Liver`: `28,520` (`31.17%`)
- `Lymph node`: `17,698` (`19.34%`)
- `Normal Pancreas`: `9,820` (`10.73%`)

Sections les plus utilisées dans ce projet :

- `10` sections primaires PDAC :
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
- `3` sections pancréas normal :
  - `IU_PDA_NP2`
  - `IU_PDA_NP10`
  - `IU_PDA_NP11`

## 7.3 Quels types de labels existent déjà dans l'objet spatial

Principales couches d'annotation utilisées dans le projet :

- `cc_ischia_10` / `CompositionCluster_CC` pour les écotypes
- `first_type` pour le type cellulaire dominant par spot
- couches liées à RCTD pour les informations de composition

L'objet est donc déjà suffisamment riche pour supporter :

- des analyses orientées écotypes
- des analyses de voisinage fondées sur le type dominant
- des overlays de scores de pathways

## 7.4 Composition tissulaire globale et pourquoi cela compte

Comptes principaux de `first_type` :

- `myCAF`: `26,726` (`29.21%`)
- `Tumor Epithelial cells`: `24,627` (`26.92%`)
- `Normal Epithelial cells`: `17,997` (`19.67%`)
- `Hepatocytes`: `9,668` (`10.57%`)
- `B cells`: `4,407` (`4.82%`)
- `iCAF`: `2,246` (`2.45%`)

Principaux écotypes :

- `CC12`: `11.22%`
- `CC15`: `10.61%`
- `CC2`: `10.60%`
- `CC14`: `7.85%`
- `CC1`: `7.76%`
- `CC8`: `7.69%`
- `CC5`: `5.99%`
- `CC3`: `3.50%`

Classes de confiance des spots :

- `doublet_certain`: `69.79%`
- `singlet`: `16.95%`
- `doublet_uncertain`: `12.12%`
- `reject`: `1.14%`

Il s'agit d'une limite cruciale et d'une règle d'interprétation :

- les spots Visium sont mixtes
- les conclusions spatiales doivent être lues comme des conclusions **au niveau territoire**, et non comme des affirmations strictes une cellule - un gène

## 7.5 Premier niveau spatial : cartes et projections de marqueurs

La première étape a été descriptive mais essentielle :

- cartes d'écotypes
- cartes `first_type`
- projections des KCN
- projections de marqueurs

Marqueurs de référence utilisés pour interpréter la structure du tissu :

- stroma contractile / ECM :
  - `ACTA2`
  - `CCN2`
  - `POSTN`
  - `TAGLN`
  - `COL12A1`
  - `THBS1`
- épithélial / ductal :
  - `EPCAM`
  - `KRT20`
- résident / immune-border :
  - `C7`
  - `C1R`
  - `CCL19`
  - `CCL21`
  - `MFAP4`
  - `CLU`
  - `DCN`

Ces cartes ne prouvent pas une niche à elles seules, mais elles sont nécessaires pour apprendre à lire le tissu et interpréter correctement les résultats quantitatifs suivants.

## 7.6 Détection par écotype fondée sur la positivité

Une première analyse spatiale a demandé :

> dans quels écotypes les gènes KCN sont-ils plus souvent détectés ?

Cette analyse utilisait :

- `Origin %in% c("Normal Pancreas", "Pancreas")`
- l'union des 25 gènes KCN
- l'expression dans `Spatial/data`
- une positivité définie comme `expression > 0`

Ce premier niveau mesurait :

- le **pourcentage de spots positifs**

Cela est informatif pour des gènes KCN clairsemés, mais ne capture pas encore l'intensité ni l'hétérogénéité inter-sections.

## 7.7 Enrichissement de niche par Fisher exact : bon premier pas, mais ensuite insuffisant

### Définitions initiales des niches

- `CC1+CC5`
- `CC2+CC3`

### Question initiale

Pour chaque `gene x niche` dans le pancréas primaire :

> les spots KCN-positifs sont-ils surreprésentés à l'intérieur de la niche par rapport à l'extérieur ?

### Test

Pour chaque gène et chaque niche :

- construction d'une table `2x2` :
  - inside / outside niche
  - KCN+ / KCN-
- application de `fisher.test()`

Principaux moffitt_stromal_exploration_outputs :

- `odds ratio`
- `p-value`
- `BH FDR`
- `% positive inside`
- `% positive outside`

### Principaux premiers résultats

Dans `CC2+CC3` :

- `KCNK1`: `OR = 4.29`, `FDR = 0`
- `KCNN4`: `OR = 2.88`, `FDR = 1.09e-211`
- `KCNJ16`: `OR = 6.02`, `FDR = 9.26e-73`
- `KCNC4`: `OR = 2.06`, `FDR = 3.67e-58`
- `KCNQ4`: `OR = 2.24`, `FDR = 8.21e-29`
- `KCNK6`: `OR = 1.52`, `FDR = 3.73e-29`

Dans `CC1+CC5` :

- `KCNMA1`: `OR = 2.39`, `FDR = 8.20e-206`
- `KCNE4`: `OR = 2.32`, `FDR = 6.46e-113`
- `KCNMB1`: `OR = 2.00`, `FDR = 3.60e-63`
- `KCND2`: `OR = 1.83`, `FDR = 2.01e-27`

### Pourquoi c'était un bon premier pas

Cela donnait une première réponse très simple et défendable :

> les spots positifs pour ce gène KCN sont plus fréquents dans cette niche que prévu par hasard

### Pourquoi cela est ensuite devenu insuffisant

Le cadre Fisher reste limité parce qu'il :

- binarise l'expression
- pool trop grossièrement l'information
- ignore l'intensité d'expression
- ignore l'hétérogénéité entre sections
- n'intègre pas explicitement la dépendance spatiale

C'est pour cela que le projet l'a ensuite remplacé par un cadre **Wilcoxon / log2FC / AUC par coupe** plus informatif.

## 7.8 Analyse par coupe écotype vs reste : Wilcoxon, AUC, log2FC, FDR

Il s'agit de l'une des mises à jour majeures ajoutées après la première synthèse.

### Objectif

Pour chaque section primaire et chaque gène KCN :

- comparer `CC1+CC5` au reste de la section
- comparer `CC2+CC3` au reste de la section

### Principales métriques

- `log2FC`
- `AUC`
- `%Pos in`
- `%Pos out`
- `dPct = %Pos in - %Pos out`
- `p-value` de Wilcoxon
- `FDR`

### Comment interpréter ces métriques

- `log2FC > 0`: expression plus élevée dans la niche testée
- `AUC > 0.5`: les spots de la niche tendent à exprimer davantage le gène que le reste
- `FDR < 0.05`: effet soutenu statistiquement dans la section
- nombre de sections significatives : robustesse à travers les lames

### Pourquoi cette approche est meilleure que Fisher

Elle conserve :

- l'intensité d'expression
- l'hétérogénéité inter-sections
- une taille d'effet plus informative
- une structure répétée sur les 10 sections primaires

### Principaux résultats KCN

Bloc `CC2+CC3` fort :

- `KCNK1`: `10/10` sections positives, `10/10` avec `FDR < 0.05`, `median log2FC = 0.414`, `median AUC = 0.610`
- `KCNN4`: `10/10` positives, `10/10` significatives, `median log2FC = 0.257`, `median AUC = 0.581`
- `KCNK6`: `10/10` positives, `8/10` significatives, `median AUC = 0.537`

Signaux `CC1+CC5` plus hétérogènes :

- `KCNMA1`: `5/10` sections avec `log2FC` positif, `3/10` avec `FDR < 0.05`, `median log2FC ~ 0`, `median AUC = 0.503`
- `KCNMB1`: `4/10` positives, `3/10` significatives, `median AUC = 0.499`
- `KCNE4`: `9/10` positives, `5/10` significatives, `median log2FC = 0.033`, `median AUC = 0.512`

### Interprétation importante

Cette analyse a donné un résultat conceptuel décisif :

- `KCNN4` et `KCNK1` se comportent comme de vrais gènes associés à un écotype
- `KCNMA1` **ne** se comporte pas comme un marqueur d'écotype propre

Ce n'est pas un échec. Cela suggère que `KCNMA1` est mieux capturé par :

- des patches spatiaux
- des programmes locaux
- des territoires stromaux continus

plutôt que par un simple label d'écotype discret.

## 7.9 Contrôles marqueurs dans le même cadre écotype vs reste

Une analyse contrôle a été relancée avec des marqueurs stromaux canoniques :

- `ACTA2`
- `CCN2`
- `POSTN`
- `TAGLN`
- `COL12A1`

Principaux résultats résumés dans `CC1+CC5` :

- `CCN2`: `median log2FC = 0.552`, `median AUC = 0.658`, `7/10` sections significatives
- `ACTA2`: `0.434`, `0.630`, `10/10`
- `POSTN`: `0.429`, `0.606`, `10/10`
- `TAGLN`: `0.333`, `0.599`, `10/10`
- `COL12A1`: `0.300`, `0.574`, `8/10`

### Pourquoi c'était très important

Cela a montré deux choses en même temps :

1. l'approche écotype vs reste capte bien un vrai signal stromal
2. même les marqueurs myCAF / ECM classiques ne sont pas parfaitement purs ni binaires dans ces sections

Cela soutient directement l'idée que le problème n'est **pas** que `KCNMA1` soit non pertinent. Le problème est qu'un contraste binaire écotype vs reste est structurellement trop grossier pour un signal territorial stromal continu.

## 7.10 Hotspot : quels gènes KCN forment de vrais patches spatiaux ?

### Script principal

- `code/spatial_pdac_analysis/communication/hotspot/all_kcns_primary/run_hotspot_all_kcns_primary.py`

### Pourquoi Hotspot était nécessaire

Les contrastes Wilcoxon par écotype demandent :

- le gène est-il plus haut dans une niche que dans le reste ?

Hotspot demande :

- le gène forme-t-il un vrai motif spatial local dans le tissu ?

### Paramètres principaux

- `10` sections primaires
- modèle : `bernoulli`
- `n_neighbors = 30`
- `MIN_GENE_SPOTS = 10`
- borne supérieure de la palette de couleur au quantile `0.95`

### Principales métriques

- `C`: autocorrélation spatiale locale
- `Z`: force spatiale standardisée
- `FDR`

Important :

- la couleur de la carte représente les counts bruts
- pas directement le score `Z`

### Principale synthèse Hotspot

| Gène | Sections testées | Sections significatives | median Z | Interprétation |
|---|---:|---:|---:|---|
| `KCNMA1` | 10 | 7 | `7.01` | très fort gène de niche |
| `KCNK1` | 10 | 7 | `3.52` | fort gène de niche spatiale |
| `KCNMB1` | 10 | 7 | `3.46` | fort gène de niche spatiale |
| `KCNK3` | 10 | 6 | `2.63` | robuste mais plus résident |
| `KCNN4` | 10 | 5 | `2.46` | robuste, plus tumoral/ductal |
| `KCND3` | 10 | 5 | `2.00` | niche récurrente, plus faible |

### Pourquoi ce résultat a été décisif

Il a clairement montré que :

- `KCNMA1` est l'un des gènes spatiaux les plus forts du PDAC primaire
- même s'il n'est pas l'un des meilleurs marqueurs binaires d'écotype

Il s'agit de l'un des résultats les plus importants de tout le projet.

## 7.11 SPARK-X : validation indépendante de la variabilité spatiale

### Script principal

- `code/spatial_pdac_analysis/sparkx/run_sparkx_all_kcns_primary.R`

### Pourquoi ajouter SPARK-X ?

Hotspot est excellent pour les patches locaux compacts. SPARK-X fournit un test complémentaire de dépendance à l'espace à l'aide de kernels multiples et peut capturer des motifs plus lisses ou plus diffus.

### Paramètres principaux

- analyse section par section
- counts bruts plus coordonnées
- `sparkx(..., option = "mixture")`
- multiples kernels
- correction multiple par `Benjamini-Yekutieli`

### Principale synthèse SPARK-X

- `KCNN4`: `10/10` sections significatives, `median adjusted p = 2.97e-18`, `median best kernel stat = 87.7`
- `KCNMA1`: `10/10`, `5.21e-10`, `49.9`
- `KCNK6`: `10/10`
- `KCNJ8`: `10/10`
- `KCNK1`: `9/10`
- `KCNS3`: `9/10`
- `KCNAB2`: `9/10`

### Interprétation

SPARK-X renforce fortement :

- le bloc ductal / tumoral
- le bloc contractile / ECM
- une partie du bloc résident / border

La convergence entre Hotspot et SPARK-X est particulièrement convaincante.

## 7.12 Local correlations : avec quels gènes chaque KCN co-varie-t-il spatialement ?

### Script principal

- `code/spatial_pdac_analysis/communication/hotspot/all_kcns_primary/local_correlations/run_kcn_local_correlations_primary.py`

### Logique

Pour chaque gène KCN :

- garder les gènes spatialement informatifs (`FDR < 0.05`, `C > 0`)
- calculer les corrélations locales dans le même cadre Hotspot
- résumer entre sections avec :
  - `n_sections_positive_z`
  - `median_local_corr_z`

### Pourquoi cette analyse est puissante

Elle répond à :

> avec quel programme moléculaire local ce gène KCN vit-il dans le tissu ?

Ce n'est pas la même chose qu'une corrélation globale sur l'ensemble de l'objet. Deux gènes peuvent partager un programme spatial local même si leur corrélation moyenne globale reste modérée.

### Principaux exemples

#### `KCNN4`

Partenaires locaux récurrents :

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

Interprétation :

- territoire **tumoral épithélial / ductal** net

#### `KCNK1`

Paires récurrentes fortes :

- `KCNK1-SLPI`: `10/10`, `median z = 4.62`
- `KCNK1-ELF3`: `10/10`, `4.49`
- `KCNK1-LCN2`: `10/10`, `3.80`

Interprétation :

- même monde ductal / épithélial que `KCNN4`

#### `KCNMA1`

Dans les moffitt_stromal_exploration_outputs détaillés par section, notamment en `T1`, les top partenaires incluent :

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

Interprétation :

- très compatible avec un **territoire contractile / myofibroblastique / de remodelage ECM**

#### Gènes résidents / immune-border

Exemples :

- `KCNK3-C7`: `10/10`, `median z = 2.24`
- `KCNJ8-C7`: `10/10`, `1.35`
- `KCND3-C7`: `9/10`, `3.10`
- `KCNK3-CLU`: `9/10`, `2.91`

Interprétation :

- programmes de type résident / périvasculaire / immune-border

## 7.13 Territory GSEA sur les spots KCN-high

### Question biologique

Après avoir identifié des patches spatiaux, la question suivante était :

> quel programme biologique caractérise le territoire où le signal KCN est le plus fort ?

### Script principal

- `code/spatial_pdac_analysis/communication/hotspot/all_kcns_primary/territory_gsea/run_all_kcn_territory_gsea.R`

### Définition du territoire

Pour chaque section primaire :

- spots positifs : `KCN > 0` dans `SCT/data`
- au moins `30` spots positifs requis
- `KCN-high` = top `10%` des spots positifs
- au moins `20` spots high requis

### Pourquoi cette définition est utile

Elle évite :

- de réduire tout le tissu à un seul contraste
- de définir le territoire uniquement sur présence / absence
- de surinterpréter des sections avec très peu de spots positifs

### Cadre GSEA

- assay utilisé pour le ranking : `SCT`
- score par gène = `mean(high) - mean(rest)`
- Hallmarks via `msigdbr`
- overlap minimum `10`
- `fgseaMultilevel(eps = 0)`

### Principaux résultats par monde KCN

#### Monde stromal contractile / ECM

`KCNMA1` :

- `10` sections utilisables
- `21` pathways significatifs au total
- `3` positifs et `18` négatifs
- top positif :
  - `HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION`
  - `NES = 1.63`
  - `FDR = 3.55e-08`
- top négatif :
  - `HALLMARK_MYC_TARGETS_V1`
  - `NES = -1.70`
  - `FDR = 8.42e-15`

`KCNMB1` :

- top positif :
  - `HALLMARK_MYOGENESIS`
  - `NES = 2.13`

`KCNE4` :

- top positif :
  - `HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION`
  - `NES = 2.28`

Interprétation :

- forte compatibilité avec une biologie de territoire contractile / ECM / myofibroblastique

#### Monde tumoral ductal / épithélial

`KCNN4` :

- `16` pathways significatifs
- top négatif :
  - `HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION`
  - `NES = -2.47`
  - `FDR = 1.35e-22`

`KCNK1` :

- top négatif :
  - `EMT`
  - `NES = -2.46`
  - `FDR = 1.30e-24`

Interprétation :

- ces territoires s'opposent à la biologie stromale contractile et correspondent beaucoup mieux à des programmes tumoraux ductaux

#### Monde résident / border

`KCND3`, `KCNJ8`, `KCNAB2` :

- signaux plus mixtes et moins massifs
- cohérents avec des contextes résidents / périvasculaires / de bordure

## 7.14 COMMOT

COMMOT a été testé sur des candidats sélectionnés pour demander si les territoires `KCN-high` se situaient dans des contextes particulièrement forts de communication sécrétée.

Résultat principal :

- aucune grande niche hyper-communicante en signalisation sécrétée n'a émergé comme explication dominante

Pourquoi ce résultat négatif reste utile :

- il suggère que la biologie la plus informative de ces gènes KCN est davantage liée :
  - à l'architecture spatiale
  - à l'ECM
  - à la contractilité
  - à l'organisation des interfaces tissulaires

qu'à une augmentation globale d'une communication sécrétée.

## 7.15 Analyse de voisinage avec `first_type`

### Script principal

- `code/spatial_pdac_analysis/cell_type_neighborhood/run_kcn_high_first_type_knn_neighborhood.py`

### Question

Autour des spots `KCN-high`, quels types dominants de spots sont retrouvés plus fréquemment ?

### Paramètres principaux

- high spots = top `10%` des spots positifs
- au moins `20` spots high
- au moins `30` spots positifs
- `K = 15` plus proches voisins
- comparaison du voisinage autour des spots high vs autres spots

Principales métriques résumées :

- `median_high_neighbor_fraction`
- `median_other_neighbor_fraction`
- `median_delta_neighbor_fraction`
- `median_enrichment_ratio`
- p-values Wilcoxon et FDRs

### Principaux résultats

`KCNN4-high` :

- voisinage enrichi en `Tumor Epithelial cells`
- `9/10` sections avec delta positif
- `median delta neighbor fraction = 0.0323`

`KCNK1-high` :

- enrichissement tumoral épithélial également
- `7/10` positifs

`KCNMA1-high` :

- déplété en voisins épithéliaux :
  - `Normal Epithelial cells`: `median delta = -0.058`
  - `Tumor Epithelial cells`: `-0.039`
- tendance positive vers `myCAF` :
  - `median delta = 0.0558`

`KCND3-high` :

- tendances plus faibles et plus mixtes
- compatibles avec une lecture résident / immune-border

### Interprétation

Cela ne prouve pas que le `first_type` voisin exprime le gène KCN. Cela décrit le **contexte tissulaire local** du territoire KCN-high.

## 7.16 Overlays continus de pathways avec `escape` : pourquoi ils ont été ajoutés

La logique écotype vs reste restait peu satisfaisante pour `KCNMA1`. Cela a motivé l'ajout d'une nouvelle couche :

- panneau gauche : projection spatiale du gène façon Hotspot
- panneau droit : score de pathway par spot calculé avec `escape`

Ensuite, pour chaque section :

- les scores de pathway ont été standardisés
- une corrélation de `Spearman` a été calculée entre :
  - les counts bruts du gène KCN
  - le score standardisé du pathway

Ce n'est pas une métrique formelle de co-localisation locale, mais c'est une très bonne mesure de **concordance spot-level continue**.

## 7.17 Overlays de pathways pour `KCNMA1`

Principaux pathways testés :

- `Cancer_associated_fibroblasts`
- `Matrix`
- `Matrix_remodeling`
- `myCAF_contractile`
- contrôle : `NK_cells`

Principaux résultats résumés :

- `Cancer_associated_fibroblasts`: `median rho = 0.154`, `10/10` sections `FDR < 0.05`
- `Matrix`: `0.112`, `10/10`
- `Matrix_remodeling`: `0.116`, `10/10`
- `myCAF_contractile`: `0.147`, `10/10`, `max rho = 0.402`
- `NK_cells`: `median rho = -0.002`, `5/10`

### Interprétation

Il s'agit d'une avancée conceptuelle majeure :

- `KCNMA1` ne suit pas n'importe quel pathway
- il s'aligne bien mieux avec un **programme stromal contractile / CAF / ECM continu**
- la meilleure interprétation n'est pas "`KCNMA1` est un parfait marqueur d'écotype"
- la meilleure interprétation est :
  - `KCNMA1` marque un **territoire stromal contractile continu**

## 7.18 Overlays de pathways pour `KCNN4`

Principaux modules tumoraux testés :

- `ductal_tumoral_epithelial`
- `tumor_stress_secretory`

Principaux résultats résumés :

- `ductal_tumoral_epithelial`: `median rho = 0.235`, `10/10` sections significatives, `max rho = 0.423`
- `tumor_stress_secretory`: `median rho = 0.215`, `10/10`, `max rho = 0.375`

### Interprétation

Il s'agit de l'une des validations les plus nettes du projet :

- `KCNN4` suit fortement un programme tumoral ductal / épithélial continu
- davantage qu'un programme stromal

## 7.19 Overlay de pathway pour `KCND3`

Module principal testé :

- `iCAF_immune_border`

Résultat résumé principal :

- `median rho = 0.064`
- `8/10` sections avec `FDR < 0.05`
- `max rho = 0.120`

### Interprétation

Ce résultat est plus faible que pour `KCNMA1` ou `KCNN4`, mais reste informatif :

- `KCND3` ne recouvre pas fortement un module unique très défini
- il montre plutôt une association positive faible mais répétée
- cela soutient une interprétation plus lâche de type résident / immune-border

## 7.20 Panel spatial de marqueurs relancé avec Hotspot

Un run Hotspot dédié a aussi été réalisé pour des marqueurs de référence comme :

- `ACTA2`
- `CCN2`
- `POSTN`
- `TAGLN`
- `COL12A1`

Cela était important car cela a montré que de grands marqueurs stromaux canoniques peuvent former des motifs plus lisses et visuellement plus graduels que certains gènes KCN.

Cela a aidé à expliquer pourquoi :

- `KCNMA1` peut avoir un rendu visuel plus binaire dans certaines sections
- tout en restant statistiquement un vrai patch spatial

## 7.21 Extension spatiale au pancréas normal

Des analyses supplémentaires ont été ajoutées sur le pancréas normal :

- Hotspot
- local correlations
- territory GSEA

Sections normales principales :

- `IU_PDA_NP10`
- `IU_PDA_NP11`
- `IU_PDA_NP2`

### Principales leçons de l'extension pancréas normal

- plusieurs gènes KCN ont déjà des motifs spatiaux structurés dans le pancréas normal
- cela évite de surinterpréter tous les signaux comme purement spécifiques de la tumeur

Exemples :

- `KCNJ8` et `KCNK3` soutiennent un axe résident / périvasculaire déjà visible hors tumeur
- `KCNK1` garde un ancrage épithélial en tissu normal, mais pas avec le même programme qu'en PDAC
- `KCNMA1` et `KCNE4` montrent qu'une partie de la logique stromale / ECM existe déjà dans le pancréas normal mais est amplifiée ou réorganisée dans le PDAC

## 7.22 Extension spatiale ciblée : `SIGMAR1`

`SIGMAR1` a été analysé dans le même cadre spatial sans l'inclure artificiellement dans le panel KCN.

Principales branches :

- Hotspot
- local correlations
- neighborhood
- territory GSEA

Interprétation principale prudente :

- `SIGMAR1` est plus faible que les meilleurs gènes KCN spatiaux
- il ne définit pas une grande niche robuste comparable à `KCNMA1` ou `KCNN4`
- son signal semble plus diffus et plus hétérogène

## 7.23 Extension spatiale ciblée : `KCNMA1 + KCNMB1` et `KCNMA1 + KCNMB4`

Deux pseudo-modules ont été construits au niveau spot :

- `KCNMA1_KCNMB1 = KCNMA1 + KCNMB1`
- `KCNMA1_KCNMB4 = KCNMA1 + KCNMB4`

Objectif :

- tester si l'axe BK devient plus lisible lorsque `KCNMA1` est lu avec une sous-unité bêta

Interprétation principale :

- les deux pseudo-modules tombent clairement dans un monde `myCAF` / contractile / ECM
- `KCNMA1 + KCNMB1` est plus convaincant spatialement que `KCNMA1 + KCNMB4`

Cela soutient l'idée que :

- le partenaire BK spatialement le plus cohérent pour `KCNMA1` dans ce projet est `KCNMB1`

## 7.24 Principales conclusions du bloc spatial

La couche spatiale a complètement réorganisé l'interprétation biologique.

Elle a montré que :

- `KCNMA1` est mieux lu comme un gène de territoire stromal contractile
- `KCNN4` et `KCNK1` sont mieux lus comme des gènes de territoire tumoral ductal
- `KCNJ8`, `KCND3`, `KCNK3` et `KCNAB2` soutiennent un monde résident / border / périvasculaire

Surtout, le bloc spatial explique pourquoi certains gènes sont mal résumés par les labels de sous-types single-cell ou par la survie bulk seule.

## 7.25 Principales limites du bloc spatial

- Visium reste du spot-level, pas du single-cell pur
- les doublets et spots mixtes sont fréquents
- la co-organisation spatiale n'est pas équivalente à une co-expression dans la même cellule
- Fisher, Hotspot, SPARK-X, local correlations et overlays n'impliquent pas une causalité
- les overlays de pathways reposent sur des modules curés et des corrélations spot-level, pas sur une validation mécanistique directe

---

## 8. Synthèse biologique intégrée à travers toutes les couches

## 8.1 Le panel des 25 KCN n'est pas un programme unique

C'est l'une des conclusions les plus fortes du projet.

Le panel KCN se répartit en au moins trois blocs biologiquement cohérents.

## 8.2 Bloc tumoral ductal / épithélial

Gènes principaux :

- `KCNK1`
- `KCNN4`
- `KCNK6`
- partiellement `KCNC4`, `KCNQ4`, `KCNJ15`

Indices convergents :

- enrichissement dans `CC2+CC3`
- forte variabilité spatiale
- local correlations avec :
  - `ELF3`
  - `LCN2`
  - `SLPI`
  - `KRT7`
  - `KRT18`
  - `MUC1`
- territory GSEA opposée à `EMT` / ECM
- overlays de pathways tumoraux

Lecture biologique :

- ces gènes ne sont pas principalement des marqueurs purement fibroblastiques
- ils sont liés à des territoires tumoraux épithéliaux ou ductaux
- ils restent hautement pertinents car les interfaces tumeur-stroma sont cruciales pour l'organisation locale de l'immunité et de l'ECM

## 8.3 Bloc stromal contractile / ECM / myofibroblastique

Gènes principaux :

- `KCNMA1`
- `KCNMB1`
- `KCNE4`
- `KCND2`
- partiellement `KCNS3`

Indices convergents :

- interprétation orientée `myCAF` en single-cell
- tendance à l'enrichissement dans `CC1+CC5`
- Hotspot fort pour `KCNMA1` et `KCNMB1`
- local correlations avec :
  - `THBS1`
  - `CCN2`
  - `COL12A1`
  - `CTHRC1`
  - `ACTG2`
  - `CNN1`
- territory GSEA vers :
  - `EMT`
  - `MYOGENESIS`
  - une biologie riche en ECM
- overlays de pathways vers :
  - CAF
  - Matrix
  - Matrix remodeling
  - `myCAF_contractile`

Lecture biologique :

- c'est le bloc le plus directement pertinent pour **ECM plus modulation immunitaire indirecte**
- ces gènes sont les meilleurs candidats pour des validations mécanistiques stromales

## 8.4 Bloc résident / périvasculaire / immune-border

Gènes principaux :

- `KCNJ8`
- `KCNK3`
- `KCND3`
- `KCNAB2`
- `KCNN3`

Indices convergents :

- comportement plus résident ou mixte en single-cell
- local correlations avec :
  - `C7`
  - `CCL19`
  - `CCL21`
  - `MFAP4`
  - `CLU`
- variabilité spatiale forte ou modérée
- patterns `iCAF_immune_border` plus faibles mais récurrents pour `KCND3`
- voisinages compatibles avec des contextes de bordure, immunitaires ou périvasculaires

Lecture biologique :

- ces gènes sont moins purement contractiles que le bloc `KCNMA1`
- ils sont particulièrement intéressants pour la question de la géographie immunitaire stromale et des territoires de bordure

---

## 9. Candidats les plus convaincants du projet

Si l'on combine :

- robustesse inter-méthodes
- cohérence biologique
- signal spatial
- pertinence pour la question du stage

les candidats les plus forts sont :

1. `KCNMA1`
2. `KCNMB1`
3. `KCNK1`
4. `KCNN4`
5. `KCNE4`
6. `KCND2`
7. `KCNJ8`
8. `KCNK3`

### Tableau de lecture pratique

| Gène | Lecture principale | Lien avec la question immunité / ECM |
|---|---|---|
| `KCNMA1` | territoire stromal contractile | très fort |
| `KCNMB1` | territoire stromal contractile / périvasculaire | très fort |
| `KCNE4` | fibro activé / ECM | fort |
| `KCND2` | myCAF / ECM / EMT-like | fort |
| `KCNK1` | interface tumorale ductale / épithéliale | fort mais plus indirect |
| `KCNN4` | interface tumorale ductale / épithéliale | fort mais plus indirect |
| `KCNJ8` | résident / périvasculaire | fort pour la logique bordure / territoire |
| `KCNK3` | résident / iCAF-like / périvasculaire | fort pour la logique inflammatory-border |
| `KCND3` | contexte immune-border / résident plus faible | modéré mais cohérent |
| `KCNAB2` | contexte immune-border / résident | prometteur mais moins robuste |

---

## 10. Robustesse globale du projet

Le projet est robuste comme étude bio-informatique exploratoire multicouche car il combine :

- plusieurs modalités de données
- plusieurs cadres DEG
- une validation externe
- plusieurs méthodes spatiales indépendantes
- des couches d'interprétation à la fois biologiques et cliniques

Les conclusions les plus robustes sont :

- l'existence du panel final de 25 gènes KCN
- la séparation en au moins trois grands blocs biologiques
- le rôle spatial fort de :
  - `KCNMA1`
  - `KCNMB1`
  - `KCNK1`
  - `KCNN4`
- la lecture fortement contractile / ECM de :
  - `KCNMA1`
  - `KCNMB1`
  - `KCNE4`
  - `KCND2`

Des éléments plus exploratoires restent :

- les gènes à support `1/4`
- certaines interprétations issues des bulk de perturbation locaux
- l'interprétation mécanistique fine des gènes immune-border

---

## 11. Principales limites

## 11.1 Limites conceptuelles

- le projet est transcriptomique, pas fonctionnel
- il identifie des associations et des priorités, pas des mécanismes causaux

## 11.2 Limites single-cell

- les labels de sous-types restent des abstractions de continus biologiques
- `iCAF` est sous-représenté
- les gènes KCN clairsemés restent sensibles au modèle

## 11.3 Limites de survie bulk

- le signal mélange tumeur, stroma et immunité
- l'harmonisation intra-cohorte réduit les effets de plateforme sans tout supprimer
- certains gènes violent l'hypothèse de risques proportionnels

## 11.4 Limites des bulk locaux de perturbation

- faibles effectifs
- absence de purification stromale
- dépendance au contexte expérimental du laboratoire

## 11.5 Limites spatiales

- Visium est spot-level et non single-cell pur
- la co-organisation spatiale n'est pas équivalente à une co-expression dans la même cellule
- les overlays de pathways sont informatifs mais pas des readouts causaux directs

## 11.6 Limite biologique

Le projet ne teste **pas** directement :

- le recrutement de cellules immunes après perturbation du canal
- le remodelage ECM après perturbation du canal
- l'activité électrophysiologique réelle des canaux
- des co-cultures ou interactions fonctionnelles directes

Le projet priorise donc des candidats et des niches plutôt qu'il ne clôt la question mécanistique.

---

## 12. Conclusion générale

Le projet répond à la question du stage de manière nuancée mais biologiquement forte.

Il montre que :

- le paysage KCN stromal dans le PDAC est structuré et non uniforme
- les candidats les plus pertinents se distribuent dans des territoires tissulaires distincts
- les axes les plus forts pour la question du stage sont :
  - **stromal contractile / ECM / myofibroblastique**
  - **interface tumorale ductale / épithéliale**
  - **résident / périvasculaire / immune-border**

La formulation finale la plus juste est :

> Le stage n'a pas identifié un unique canal ionique responsable de l'immunomodulation stromale dans le PDAC. Il a plutôt défini un paysage structuré de gènes KCN associés à des états stromaux distincts, des niches tissulaires distinctes et des programmes liés à l'ECM et à l'immunité, fournissant une base solide pour les validations fonctionnelles futures.

---

## 13. Suites expérimentales les plus logiques

1. **validation fonctionnelle stromale**
   - `KCNMA1`, `KCNMB1`, `KCNE4`, `KCND2`
   - contractilité, ECM, migration, sécrétion

2. **validation des interfaces tumeur-stroma**
   - `KCNK1`, `KCNN4`, `KCNK6`
   - plasticité ductale, programmes sécrétoires, comportement d'interface

3. **validation des mondes résident / immune-border**
   - `KCNJ8`, `KCNK3`, `KCND3`, `KCNAB2`
   - gradients inflammatoires, complément, contexte de bordure vasculaire

4. **essais orientés immunité**
   - co-cultures CAF + cellules immunes
   - tests de migration / infiltration après perturbation des canaux
   - readouts cytokines et ECM après modulation des canaux

---

## Principaux fichiers soutenant cette synthèse

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
