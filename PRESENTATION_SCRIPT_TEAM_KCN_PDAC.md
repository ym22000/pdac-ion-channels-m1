# Script detaille de presentation orale

## Titre
**Identification et analyse par bio-informatique des canaux ioniques qui modulent l'immunite et la matrice extracellulaire dans le stroma du PDAC**

---

## Objectif de ce document

Ce document sert de base pour une presentation orale devant une equipe de biologistes.

L'idee est de :
- garder un fil conducteur clair ;
- expliquer les methodes bioinformatiques de facon accessible ;
- montrer pourquoi chaque analyse a ete faite ;
- insister sur la partie spatiale, qui a ete la plus informative ;
- replacer chaque resultat dans l'histoire biologique des `KCN` dans le `PDAC`.

Le niveau de detail est volontairement plus riche qu'un simple plan de slides, pour que ce document puisse aussi servir de support de revision avant l'oral.

---

## Duree conseillee

- version detaillee : `15-20 min`
- version plus courte : `10-12 min`

---

## Message central a faire passer

Le message principal du projet est le suivant :

**on n'a pas identifie un seul canal ionique stromal du PDAC, mais plusieurs familles de `KCN`, chacune associee a une niche biologique differente, avec des implications distinctes pour l'immunite, la matrice extracellulaire et l'interface tumeur-stroma.**

Les trois grands mondes qui ressortent sont :
- un bloc **tumoral epithelial/ductal**
- un bloc **stromal contractile / ECM / myCAF-like**
- un bloc **resident / perivasculaire / immune-border**

---

# Slide 1 - Introduction biologique

## Titre suggere
**Pourquoi s'interesser aux canaux ioniques dans le stroma du PDAC ?**

## Ce qu'on met sur la slide
- PDAC = tumeur tres riche en stroma
- le stroma controle :
  - la matrice extracellulaire
  - la rigidite tissulaire
  - l'inflammation
  - l'exclusion immune
- les canaux ioniques controlent :
  - le potentiel membranaire
  - les flux de `Ca2+`, `K+`, `Na+`, `Cl-`
  - le volume cellulaire
  - la migration
  - la secretion
- question : quels `KCN` sont associes a l'immunite et a l'ECM dans le stroma du PDAC ?

## Script oral
Le point de depart du stage, c'est que le PDAC est une tumeur tres stromale.  
Et ce stroma n'est pas juste un decor : il participe activement a la progression tumorale.

Il controle :
- la production et le remodelage de la matrice extracellulaire,
- la rigidite du tissu,
- les interactions avec les fibroblastes,
- et aussi la facon dont les cellules immunes accedent, ou n'accedent pas, a la tumeur.

En parallele, les canaux ioniques jouent un role tres large dans la biologie cellulaire :
- regulation du calcium,
- homeostasie du potassium,
- volume cellulaire,
- migration,
- secretion,
- et signalisation.

Donc la question generale du projet a ete :

**peut-on identifier des canaux ioniques, surtout de type `KCN`, associes a des etats stromaux capables de moduler l'immunite et la matrice extracellulaire dans le PDAC ?**

## Message a marteler
Le but n'etait pas juste de faire une liste de genes, mais de comprendre **dans quels mondes biologiques** vivent ces canaux.

---

# Slide 2 - Strategie generale

## Titre suggere
**Une approche integree : single-cell, bulk, survie et spatial**

## Ce qu'on met sur la slide
1. `single-cell RNA-seq stromal`
2. signatures `myCAF/iCAF`
3. `bulk RNA-seq`
4. survie multicohorte
5. `spatial transcriptomics`

## Script oral
J'ai construit le projet par couches successives.

D'abord, il fallait identifier un panel credible de `KCN` dans le stroma.

Ensuite, il fallait leur donner un contexte transcriptionnel :
- sont-ils plus proches de `myCAF` ?
- de `iCAF` ?
- d'etats plus residents ?

Puis j'ai ajoute :
- du `bulk RNA-seq`, pour relier certains KCN a des perturbations ou a des programmes plus larges ;
- une couche rapide de survie, pour voir s'il y avait un signal clinique ;
- et surtout du `spatial transcriptomics`, qui a ete la partie la plus informative biologiquement.

Pourquoi ?
Parce que tant qu'on reste en single-cell ou en bulk, on sait qu'un gene est present ou associe a un etat, mais on ne sait pas **ou il se place dans le tissu**.

Or, pour un projet sur l'immunite et l'ECM dans le PDAC, la localisation spatiale est centrale.

---

# Slide 3 - Jeux de donnees utilises

## Titre suggere
**Les donnees mobilisees**

## Ce qu'on met sur la slide
- atlas stromal `scRNA-seq`
- objet spatial `PDAC_Updated_ST.rds`
- bulk du labo :
  - `BKCa/shBKCa`
  - `SIGMAR1/shSIGMAR1`
- cohortes bulk publiques pour la survie

## Script oral
Pour la partie single-cell, j'ai utilise un atlas stromal PDAC deja structure autour du compartiment fibroblastique/stromal.

Pour la partie spatiale, j'ai utilise un objet `Visium` PDAC comprenant :
- des tumeurs primaires,
- des metastases,
- et aussi `3` coupes de pancreas sain.

J'ai egalement exploite deux jeux `bulk` du laboratoire :
- `BKCa/shBKCa`
- `SIGMAR1/shSIGMAR1`

Enfin, j'ai utilise plusieurs cohortes bulk publiques pour une analyse de survie multicohorte.

## Comment le dire simplement
Le single-cell a servi a **selectionner** les KCN.
Le spatial a servi a **comprendre leurs niches**.
Le bulk et la survie ont servi de **couches complementaires**.

---

# Slide 4 - Pourquoi ne pas selectionner les KCN "a l'oeil" ?

## Titre suggere
**Construire un panel KCN robuste**

## Ce qu'on met sur la slide
- eviter une liste arbitraire
- croiser plusieurs methodes
- panel final = `25 KCN`

## Script oral
Une difficulte importante, c'etait d'eviter une approche arbitraire du type :
"on prend les canaux les plus connus", ou "on prend ceux qui ont l'air de sortir".

J'ai donc choisi une logique simple :

**si un KCN revient avec plusieurs methodes statistiques independantes, il a plus de chances d'etre un vrai signal biologique.**

Le but etait donc de construire un panel KCN qui soit :
- defendable,
- reproductible,
- et biologiquement interpretable.

Au final, cette strategie a conduit a un **panel de `25 KCN`**.

---

# Slide 5 - Single-cell : methodes utilisees

## Titre suggere
**Comment les KCN ont ete identifies en single-cell**

## Ce qu'on met sur la slide
- `pseudobulk DESeq2`
- `DESeq2 subtype vs rest`
- `MAST`
- modeles de type `hurdle`

## Script oral
La premiere grande etape a ete l'identification des KCN dans le compartiment stromal.

Pour cela, j'ai utilise plusieurs approches.

### 1. Pseudobulk
Le principe du pseudobulk est le suivant :
- au lieu d'analyser chaque cellule individuellement,
- on regroupe des cellules d'un meme groupe,
- on somme les counts,
- puis on applique une methode de `bulk RNA-seq`.

Pourquoi faire ca ?
Parce que le single-cell cellule par cellule est tres bruité, avec beaucoup de zéros et une forte variabilite technique.

Le pseudobulk permet de recuperer un signal plus stable au niveau d'un groupe biologique.

### 2. DESeq2
`DESeq2` modelise les counts par une **loi binomiale negative**.

Pourquoi cette loi ?
Parce que les donnees de comptage ne suivent pas une simple loi de Poisson : leur variance est souvent plus grande que leur moyenne, ce qu'on appelle la **surdispersion**.

En pratique, `DESeq2` estime :
- un `log2 fold change`
- une p-value
- puis une p-value corrigee pour tests multiples

### 3. MAST et modeles hurdle
Ces methodes sont utiles en single-cell parce qu'elles gerent mieux la structure "beaucoup de zeros".

Intuitivement, elles se posent deux questions :
1. le gene est-il detecte ou non ?
2. quand il est detecte, est-il plus exprime ?

Pourquoi c'etait utile ici ?
Parce que plusieurs `KCN` sont peu ou moderement exprimes, avec une detection rare mais biologiquement informative.

## Pourquoi plusieurs methodes ?
Parce qu'aucune methode seule n'est parfaite.

Le croisement des approches permet de limiter :
- les faux positifs,
- les effets de methode,
- les genes qui sortent "par hasard" dans une seule analyse.

---

# Slide 6 - Resultats single-cell

## Titre suggere
**Premier resultat : les KCN ne forment pas un seul programme**

## Ce qu'on met sur la slide
- au moins deux axes :
  - `myCAF / ECM / contractile`
  - `resident / inflammatory / iCAF-like`
- panel final `25 KCN`

## Script oral
Le premier resultat important est qu'il n'existe pas un seul programme KCN uniforme dans le stroma.

On voit deja, des le single-cell, plusieurs tendances.

Un premier bloc de KCN est plus compatible avec :
- le stroma active,
- la contractilite,
- l'ECM,
- et des etats `myCAF-like`

Par exemple :
- `KCNMA1`
- `KCNMB1`
- `KCNE4`
- `KCND2`

Un autre bloc est plus compatible avec :
- des etats residents,
- des bordures inflammatoires,
- des territoires perivasculaires,
- ou des lectures `iCAF-like`

Par exemple :
- `KCNJ8`
- `KCNK3`
- `KCND3`
- `KCNAB2`

Donc, des cette etape, le message est deja :

**les KCN ne racontent pas une seule biologie stromale.**

---

# Slide 7 - Pourquoi relier les KCN aux signatures myCAF/iCAF ?

## Titre suggere
**Positionner les KCN sur les grands etats stromaux**

## Ce qu'on met sur la slide
- scores de signatures
- `GSVA`
- `AUCell`
- correlations avec `myCAF` et `iCAF`

## Script oral
Une fois le panel KCN etabli, il fallait savoir a quels grands programmes stromaux ces genes ressemblent.

Pour ca, j'ai utilise des approches de type scores de signatures :
- `GSVA`
- `AUCell`

### Explication simple
Une signature, c'est une liste de genes representant un etat biologique.

Exemple :
- une signature `myCAF` contient des genes contractiles et matriciels
- une signature `iCAF` contient des genes plus inflammatoires ou de dialogue avec l'immunite

L'idee est donc :
- si un `KCN` suit plutot les genes `ACTA2`, `TAGLN`, `MYL9`, `COL11A1`, il est plus proche de `myCAF`
- s'il suit plutot des genes de complement, de cytokines ou de bordure inflammatoire, il est plus compatible avec un axe `iCAF-like`

## Ce que cette etape a apporte
Elle a servi a situer les KCN sur un axe biologique plus interpretable :
- contractile / ECM
vs
- inflammatoire / resident / border-stroma

Mais a ce stade, il manquait encore une dimension cruciale :
**l'espace.**

---

# Slide 8 - Bulk et survie : pourquoi les avoir faits ?

## Titre suggere
**Des couches complementaires, mais secondaires**

## Ce qu'on met sur la slide
- bulk labo = contexte fonctionnel
- survie = contexte clinique
- mais le coeur du projet reste spatial

## Script oral
J'ai ensuite ajoute deux couches complementaires.

### Bulk du laboratoire
Le `bulk RNA-seq` a servi a relier certains KCN a des programmes plus larges dans des perturbations biologiques.

Par exemple, le bloc `BKCa/shBKCa` a ete utile pour replacer `KCNMA1` dans une logique :
- stromale,
- matricielle,
- et inflammatoire.

### Survie
La partie survie a ete faite sur plusieurs cohortes bulk.

La methode principale a ete le **modele de Cox**.

### Explication simple
Le modele de Cox evalue si l'expression d'un gene est associee au risque d'evenement dans le temps.

On obtient un **hazard ratio** :
- `> 1` : risque plus eleve
- `< 1` : risque plus faible

### Pourquoi je dis que c'est secondaire
Parce que le bulk melange :
- tumeur,
- stroma,
- immunite,
- et parfois differences techniques.

Donc c'est utile comme contexte, mais ce n'est pas la meilleure couche pour decrire finement des niches stromales.

---

# Slide 9 - Pourquoi le spatial a ete la partie centrale

## Titre suggere
**Du gene au territoire**

## Ce qu'on met sur la slide
- `Visium` = resolution spot, pas cellule unique pure
- mais acces au contexte tissulaire
- vrai enjeu : localiser les KCN dans le tissu

## Script oral
Le spatial a ete la partie la plus informative du projet.

Pourquoi ?
Parce qu'un meme `KCN` peut etre :
- stromal dans un contexte,
- ductal dans un autre,
- perivasculaire dans un autre encore.

Et ces differences sont invisibles si on ne regarde pas le tissu.

Le `Visium` n'est pas une technologie cellule unique pure : un spot peut contenir plusieurs cellules.
Mais il permet de voir des **territoires**, des **interfaces**, des **bordures**, et des **niches dominantes**.

Et pour une question autour de l'ECM et de l'immunite, c'est exactement ce qu'il fallait.

---

# Slide 10 - Premier niveau spatial : cartes et niches

## Titre suggere
**Ou tombent les KCN dans le tissu ?**

## Ce qu'on met sur la slide
- cartes d'expression
- niches `CC1+CC5`, `CC2+CC3`
- tests d'enrichissement de Fisher

## Script oral
Le premier niveau d'analyse spatiale a ete descriptif :
- cartes d'expression des KCN
- projection sur les sections
- comparaison avec les niches ecotypiques

J'ai ensuite utilise des tests de Fisher pour demander :

**est-ce qu'un KCN est enrichi dans telle niche plus que prevu ?**

### Explication simple
On construit un tableau `2 x 2` :
- spot dans la niche / hors niche
- spot positif pour le KCN / negatif

Puis on teste si la distribution observee differe de ce qu'on attendrait si tout etait aleatoire.

### Message principal
Ca a deja fait emerger trois mondes :
- `KCNK1`, `KCNN4`, `KCNK6` vers les territoires ductaux/tumoraux
- `KCNMA1`, `KCNMB1`, `KCNE4`, `KCND2` vers les territoires stromaux/ECM
- `KCNJ8`, `KCNK3`, `KCND3`, `KCNAB2` vers des territoires plus residents ou perivasculaires

---

# Slide 11 - Hotspot : comment prouver qu'un patch spatial est reel

## Titre suggere
**Hotspot : un gene forme-t-il une vraie structure spatiale ?**

## Ce qu'on met sur la slide
- expression visuelle != spatialite statistique
- `Hotspot`
- `SPARK-X`

## Script oral
Une carte d'expression seule ne suffit pas.

On peut avoir :
- quelques spots forts mais disperses,
- ou au contraire un motif spatial reel, mais pas tres impressionnant visuellement.

J'ai donc utilise `Hotspot`.

### Idee intuitive
`Hotspot` demande :

**si un spot exprime un gene, ses voisins proches ont-ils tendance a l'exprimer aussi ?**

Si oui, cela indique une **autocorrelation spatiale locale**.

### Lecture mathematique simple
On construit un graphe de voisins dans l'espace.
Puis on regarde si les spots voisins se ressemblent plus que prevu par hasard.

Les sorties principales sont :
- `C`
- `Z`
- `FDR`

En pratique :
- `Z` eleve + `FDR < 0.05` = vraie organisation spatiale

### Pourquoi j'ai ajoute `SPARK-X`
`Hotspot` est tres bon pour des motifs locaux et des patches.
`SPARK-X` apporte une validation complementaire en testant plusieurs formes de variabilite spatiale.

Donc :
- `Hotspot` = niche locale
- `SPARK-X` = validation de spatialite globale

### Resultats forts
Les plus robustes spatialement en PDAC ont ete :
- `KCNMA1`
- `KCNK1`
- `KCNMB1`
- `KCNN4`

Et `SPARK-X` a aussi consolide :
- `KCNJ8`
- `KCNK6`
- `KCNAB2`
- `KCNE4`

---

# Slide 12 - Local correlations : de quoi sont faits les patches ?

## Titre suggere
**Hotspot local correlations : qui vit avec chaque KCN ?**

## Ce qu'on met sur la slide
- pas seulement "ou"
- mais "avec quels genes"

## Script oral
Une fois qu'on sait qu'un KCN forme une niche spatiale, la question suivante est :

**quels genes suivent localement le meme motif spatial ?**

Autrement dit :
- si `KCNMA1` est fort dans un patch,
- quels genes montent dans ce meme patch ?

Cette analyse est tres utile parce qu'elle revele la **composition moleculaire locale** de la niche.

### Resultats tres parlants
#### Bloc ductal tumoral
- `KCNK1`, `KCNN4`, `KCNK6`
- partenaires :
  - `ELF3`
  - `LCN2`
  - `SLPI`
  - `MUC1`
  - `KRT7`
  - `KRT18`
  - `S100P`

#### Bloc stromal contractile / ECM
- `KCNMA1`, `KCNMB1`, `KCNE4`, `KCND2`
- partenaires :
  - `ACTG2`
  - `PDLIM3`
  - `CTHRC1`
  - `THBS1`
  - `ITGA11`
  - `COL11A1`
  - `CCN1`
  - `CCN2`

#### Bloc resident / perivasculaire / immune-border
- `KCNJ8`, `KCND3`, `KCNAB2`, `KCNK3`
- partenaires :
  - `C7`
  - `CCL19`
  - `CCL21`
  - `MFAP4`
  - `CLU`
  - parfois `HLA-DRA`, `APOC1`, `LYZ`

### Pourquoi c'est important
Cette etape a ete une des plus fortes biologiquement, parce qu'elle a donne une **identite** a chaque niche KCN.

---

# Slide 13 - Territory GSEA : quel programme biologique caracterise le territoire KCN-high ?

## Titre suggere
**GSEA de territoire : lire les programmes biologiques des spots KCN-high**

## Ce qu'on met sur la slide
- spots `KCN-high`
- comparaison `high vs rest`
- Hallmarks

## Script oral
Apres avoir regarde les genes partenaires, j'ai voulu une lecture plus globale au niveau de programmes biologiques.

J'ai donc fait une **GSEA de territoire**.

### Principe
Pour chaque coupe, je definis les spots `KCN-high` :
- le top `10%` des spots positifs du KCN,
- avec des seuils minimaux pour garder un nombre de spots suffisant.

Ensuite, pour chaque gene, je calcule :

`moyenne dans high - moyenne dans rest`

Donc :
- valeur positive = gene enrichi dans le territoire KCN-high
- valeur negative = gene appauvri

Puis je classe les genes selon ce score, et j'applique `fgsea` sur les Hallmarks.

### Ce que ca repond
Ca ne dit pas :
- quelle cellule exprime exactement le KCN

Ca dit :
- **quel programme biologique caracterise le territoire ou ce KCN est le plus fort**

### Resultats principaux
#### Ductal tumoral
- `KCNK1` et `KCNN4`
- programmes epitheliaux/secretory/ductal

#### Stromal ECM / contractile
- `KCNMA1`, `KCNMB1`
- programmes `EMT`, `myogenesis`, ECM, stroma contractile

#### Resident / perivasculaire / immune-border
- `KCNJ8`, `KCND3`, `KCNAB2`
- programmes plus residents, de bordure, ou perivasculaires

### Point important
J'ai aussi fait ce type de GSEA sur les `3` coupes de pancreas sain, pour distinguer ce qui existe deja hors tumeur de ce qui devient reellement PDAC-like.

---

# Slide 14 - Analyse de voisinage cellulaire

## Titre suggere
**Autour des spots KCN-high, qui retrouve-t-on ?**

## Ce qu'on met sur la slide
- voisinage `KNN`
- `K = 15`
- comparaison `KCN-high` vs reste

## Script oral
La derniere grande couche spatiale a ete une analyse de voisinage.

La question etait :

**autour des spots ou un KCN est fort, quels types de spots retrouve-t-on davantage ?**

J'ai pris les spots `KCN-high`, puis j'ai regarde les `15` voisins les plus proches.

Pour chaque type cellulaire `first_type`, j'ai compare :
- sa frequence autour des spots `KCN-high`
- versus sa frequence autour du reste des spots

### Ce que ca a confirme
#### Proches des `Tumor Epithelial cells`
- `KCNK1`
- `KCNN4`

#### Proches des `myCAF`
- `KCNMA1`
- `KCNMB1`
- `KCNE4`
- `KCND2`

#### Plus mixtes, perivasculaires ou immune-border
- `KCNJ8`
- `KCND3`
- `KCNAB2`

### Comment presenter cette partie
Je la presente comme une couche **descriptive mais utile** :
- moins forte statistiquement que Hotspot ou GSEA,
- mais tres parlante pour le contexte local.

---

# Slide 15 - Trois grands mondes KCN dans le PDAC

## Titre suggere
**Le resultat central du projet**

## Ce qu'on met sur la slide
### Bloc 1
`KCNK1`, `KCNN4`, `KCNK6`

### Bloc 2
`KCNMA1`, `KCNMB1`, `KCNE4`, `KCND2`, `KCNS3`

### Bloc 3
`KCNJ8`, `KCND3`, `KCNAB2`, `KCNK3`, `KCNN3`

## Script oral
Le grand resultat du projet, c'est qu'il n'y a pas "un" programme KCN du PDAC.

Il y a trois grands mondes.

### 1. Monde ductal tumoral
- `KCNK1`
- `KCNN4`
- `KCNK6`

Ils definissent un territoire :
- epithelial
- ductal
- secretory
- tumor-associated

### 2. Monde stromal contractile / ECM / myCAF-like
- `KCNMA1`
- `KCNMB1`
- `KCNE4`
- `KCND2`
- `KCNS3`

Ce bloc est le plus directement lie au theme :
- matrice extracellulaire
- rigidite
- contractilite
- exclusion immune indirecte

### 3. Monde resident / perivasculaire / immune-border
- `KCNJ8`
- `KCND3`
- `KCNAB2`
- `KCNK3`
- `KCNN3`

Ce bloc semble plus lie :
- aux bordures,
- aux zones de transition,
- a l'interface stroma-immunite,
- et aux territoires perivasculaires.

---

# Slide 16 - Focus sur quelques KCN majeurs

## Titre suggere
**Les meilleurs candidats du projet**

## Ce qu'on met sur la slide
- `KCNMA1`
- `KCNK1`
- `KCNN4`
- `KCNJ8`

## Script oral
Si je devais retenir quelques `KCN` principaux :

### `KCNMA1`
C'est probablement le meilleur marqueur du **stroma contractile / myofibroblastique / ECM**.

Il vit avec :
- `ACTG2`
- `PDLIM3`
- `CTHRC1`
- `THBS1`

Il est :
- fort en Hotspot,
- fort en GSEA,
- fort dans le voisinage `myCAF`.

Donc c'est un excellent candidat pour l'axe **matrice extracellulaire**.

### `KCNK1`
Grand marqueur du **territoire tumoral ductal**.

Il est relie a :
- `ELF3`
- `LCN2`
- `SLPI`
- `MUC1`

Et son territoire est :
- ductal,
- secretory,
- epithelial.

### `KCNN4`
Tres important, parce qu'il pouvait sembler ambigu selon certaines couches non spatiales.

Mais spatialement, il raconte surtout une niche :
- tumorale,
- ductale,
- epitheliale,
et pas un fibro `myCAF` pur.

### `KCNJ8`
Tres beau marqueur :
- resident
- perivasculaire
- immune-border

Il est deja tres structure en pancreas sain, puis recontextualise dans le PDAC.

---

# Slide 17 - Comparaison pancreas sain vs PDAC

## Titre suggere
**Le PDAC cree-t-il ces territoires, ou les reprogramme-t-il ?**

## Ce qu'on met sur la slide
- certains signaux existent deja en normal
- le PDAC amplifie, durcit, ou redirige
- exemples :
  - `KCNJ8`
  - `KCNJ15`
  - `KCNK17`
  - `KCNMA1`

## Script oral
Le fait d'avoir ajoute les `3` coupes de pancreas sain a ete tres utile.

Ca a montre qu'il existe plusieurs cas.

### Cas 1 : territoire deja present puis reprogramme
#### `KCNJ8`
En normal, il marque deja un monde resident/perivasculaire.
En PDAC, il est plus interprete comme un territoire immune-border/perivasculaire stromal.

#### `KCNJ15`
En normal, il est deja tres structure.
En PDAC, il bascule vers une interface tumeur-remodeling.

### Cas 2 : territoire resident conserve
#### `KCNK17`
Il garde une lecture endocrine/residente, peu franchement tumorale.

### Cas 3 : territoire deja organise mais fortement amplifie en PDAC
#### `KCNMA1`
Il existe deja spatialement en normal,
mais dans le PDAC il devient un grand marqueur de stroma contractile/ECM.

### Message global
Le PDAC ne cree pas tout de novo.
Il **amplifie**, **reoriente** ou **durcit** des architectures biologiques deja presentes.

---

# Slide 18 - Limites

## Titre suggere
**Ce que le projet montre, et ce qu'il ne montre pas encore**

## Ce qu'on met sur la slide
- transcriptomique
- `Visium` = spot-level
- associations spatiales != causalite
- survie bulk peu specifique

## Script oral
Il y a plusieurs limites importantes.

### 1. Projet majoritairement transcriptomique
On parle ici de :
- genes,
- signatures,
- programmes,
- voisinages,
mais pas encore de demonstration fonctionnelle directe.

### 2. Visium = resolution spot
Un spot peut contenir plusieurs cellules.

Donc on parle plutot de :
- territoires
- niches dominantes
- interfaces

et pas d'attribution cellule unique pure.

### 3. Spatialite != causalite
Voir deux signaux co-localiser ne prouve pas qu'ils se regulent directement.

Ca fournit :
- des hypotheses robustes,
- des priorites,
- mais pas encore le mecanisme final.

### 4. Survie bulk
La survie multicohorte est utile,
mais elle melange les signaux :
- tumoraux,
- stromaux,
- immunitaires.

Donc elle reste une couche secondaire par rapport au spatial.

---

# Slide 19 - Conclusion generale

## Titre suggere
**Conclusion**

## Ce qu'on met sur la slide
- pas un seul canal, mais plusieurs familles de `KCN`
- trois grands mondes biologiques
- meilleurs candidats :
  - `KCNMA1`, `KCNMB1`, `KCNE4`, `KCND2`
  - `KCNK1`, `KCNN4`, `KCNK6`
  - `KCNJ8`, `KCND3`, `KCNAB2`

## Script oral
Pour conclure, le projet montre qu'il n'existe pas un seul canal stromal du PDAC.

Il existe plusieurs familles de `KCN`, associees a des mondes biologiques differents.

### Monde 1 : tumoral ductal
- `KCNK1`
- `KCNN4`
- `KCNK6`

### Monde 2 : stromal contractile / ECM / myCAF-like
- `KCNMA1`
- `KCNMB1`
- `KCNE4`
- `KCND2`

### Monde 3 : resident / perivasculaire / immune-border
- `KCNJ8`
- `KCND3`
- `KCNAB2`

Et le message final est :

**les canaux ioniques ne semblent pas moduler l'immunite et la matrice par un seul mecanisme, mais via plusieurs niches du microenvironnement tumoral.**

---

# Slide 20 - Perspectives

## Titre suggere
**Perspectives**

## Ce qu'on met sur la slide
- validation spatiale
- validation fonctionnelle
- co-localisation
- perturbations ciblees

## Script oral
Les suites les plus logiques seraient :

### Validation spatiale
- `RNAscope`
- immunofluorescence
- co-localisation avec :
  - `myCAF`
  - cellules tumorales ductales
  - endothelium / perivasculaire
  - cellules immunes

### Validation fonctionnelle
Tester quelques candidats forts :
- `KCNMA1`
- `KCNK1`
- `KCNN4`
- `KCNJ8`

### Question biologique de fond
Est-ce que ces canaux modulent :
- la matrice,
- la contractilite,
- la secretion,
- ou l'accessibilite immune
de facon differente selon la niche ou ils vivent ?

---

## Resume ultra court a dire si on te coupe

Si je dois resumer en une phrase :

**j'ai identifie un panel de `KCN` associes au stroma du PDAC, puis j'ai montre, surtout grace au spatial transcriptomics, que ces KCN se distribuent entre une niche ductale tumorale, une niche stromale contractile/ECM, et une niche resident/perivasculaire/immune-border.**

---

## Resume final a connaitre presque par coeur

Le projet n'a pas mis en evidence un seul canal ionique central du stroma PDAC, mais plusieurs familles de `KCN` qui occupent des territoires differents du microenvironnement.

Le bloc le plus lie a la matrice extracellulaire et au stroma active comprend surtout :
- `KCNMA1`
- `KCNMB1`
- `KCNE4`
- `KCND2`

Le bloc le plus lie au territoire tumoral ductal comprend surtout :
- `KCNK1`
- `KCNN4`
- `KCNK6`

Et un troisieme bloc, plus resident, perivasculaire ou immune-border, comprend :
- `KCNJ8`
- `KCND3`
- `KCNAB2`
- `KCNK3`

L'apport principal du travail est donc d'avoir transforme une liste de canaux en une carte de niches biologiques.

---

Last updated: 2026-05-16
