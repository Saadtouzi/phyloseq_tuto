Dada2 tutorial
================

  - [I -Importation des données](#i--importation-des-données)
  - [II- Filtration des données](#ii--filtration-des-données)
      - [1- Elaboration des listes des fichiers fastq des reads (Forward
        ans
        reverse)](#elaboration-des-listes-des-fichiers-fastq-des-reads-forward-ans-reverse)
      - [2- Les profils de qualités des
        reads](#les-profils-de-qualités-des-reads)
      - [3- Filtration des données](#filtration-des-données)
  - [III- Modèle d’erreur](#iii--modèle-derreur)
  - [IV- Inférence d’échantillon](#iv--inférence-déchantillon)
  - [V- Construction de la table
    d’observation](#v--construction-de-la-table-dobservation)
      - [1- Construction des contigues](#construction-des-contigues)
      - [2- Construction de la table
        d’observations](#construction-de-la-table-dobservations)
      - [3- Elimination des chimères](#elimination-des-chimères)
      - [4- Suivi des reads dans le
        pipeline](#suivi-des-reads-dans-le-pipeline)
  - [V- Assignation taxonomique](#v--assignation-taxonomique)
  - [VI- Evaluation de l’efficacité](#vi--evaluation-de-lefficacité)

``` r
library(rmarkdown)
library(knitr)
```

``` r
library("dada2")
```

    ## Loading required package: Rcpp

``` r
library("Rcpp")
```

\#Dans ce travail les auteurs ont travaillé avec des métadonnées
générées à partir du séquençage des communautés microbiennes
présentes dans 360 échantillons fécaux provenant de 12 souris.Le
séquençage a été réalisé par Illumina Miseq 2x250 en ciblant la région
V4 du gène de l’ARNr 16S. Les données fastq générées sont téléchargées
(01\_data-import), et importées pour les analyses biostatiques et
bio-informatiques.

# I -Importation des données

\#Les données sont importées sous le nom MiSeq\_SOP en utilisant le code
suivant.

``` r
path <- "~/phyloseq_tuto/MiSeq_SOP" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)
```

    ##  [1] "F3D0_S188_L001_R1_001.fastq"   "F3D0_S188_L001_R2_001.fastq"  
    ##  [3] "F3D1_S189_L001_R1_001.fastq"   "F3D1_S189_L001_R2_001.fastq"  
    ##  [5] "F3D141_S207_L001_R1_001.fastq" "F3D141_S207_L001_R2_001.fastq"
    ##  [7] "F3D142_S208_L001_R1_001.fastq" "F3D142_S208_L001_R2_001.fastq"
    ##  [9] "F3D143_S209_L001_R1_001.fastq" "F3D143_S209_L001_R2_001.fastq"
    ## [11] "F3D144_S210_L001_R1_001.fastq" "F3D144_S210_L001_R2_001.fastq"
    ## [13] "F3D145_S211_L001_R1_001.fastq" "F3D145_S211_L001_R2_001.fastq"
    ## [15] "F3D146_S212_L001_R1_001.fastq" "F3D146_S212_L001_R2_001.fastq"
    ## [17] "F3D147_S213_L001_R1_001.fastq" "F3D147_S213_L001_R2_001.fastq"
    ## [19] "F3D148_S214_L001_R1_001.fastq" "F3D148_S214_L001_R2_001.fastq"
    ## [21] "F3D149_S215_L001_R1_001.fastq" "F3D149_S215_L001_R2_001.fastq"
    ## [23] "F3D150_S216_L001_R1_001.fastq" "F3D150_S216_L001_R2_001.fastq"
    ## [25] "F3D2_S190_L001_R1_001.fastq"   "F3D2_S190_L001_R2_001.fastq"  
    ## [27] "F3D3_S191_L001_R1_001.fastq"   "F3D3_S191_L001_R2_001.fastq"  
    ## [29] "F3D5_S193_L001_R1_001.fastq"   "F3D5_S193_L001_R2_001.fastq"  
    ## [31] "F3D6_S194_L001_R1_001.fastq"   "F3D6_S194_L001_R2_001.fastq"  
    ## [33] "F3D7_S195_L001_R1_001.fastq"   "F3D7_S195_L001_R2_001.fastq"  
    ## [35] "F3D8_S196_L001_R1_001.fastq"   "F3D8_S196_L001_R2_001.fastq"  
    ## [37] "F3D9_S197_L001_R1_001.fastq"   "F3D9_S197_L001_R2_001.fastq"  
    ## [39] "filtered"                      "HMP_MOCK.v35.fasta"           
    ## [41] "Mock_S280_L001_R1_001.fastq"   "Mock_S280_L001_R2_001.fastq"  
    ## [43] "mouse.dpw.metadata"            "mouse.time.design"            
    ## [45] "stability.batch"               "stability.files"

# II- Filtration des données

## 1- Elaboration des listes des fichiers fastq des reads (Forward ans reverse)

\#On cherche à ordonner dans la première partie du code les fichier
fastq en les mettant dans le même ordre. Dans la deuxième partie du code
on défini un modèle pour nommer les échantillons fastq.

``` r
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```

## 2- Les profils de qualités des reads

``` r
plotQualityProfile(fnFs[1:2])
```

![](02_dada-analysist_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->
\#Le plot montre les score de qualités des reads forward à chaque
position (nucléotide). La courbe en vert représente le score qualité
moyen à chaque point, en orange les quartiles de la distribution de ces
points et en gris une carte thermique de la fréquence de chaque score.
On remarque que les scores qualités diminuent considérablement dans les
dernières positions à cause des erreurs de séquençage. Donc ce qu’il
faut faire, c’est tronquer la première partie des reads en éliminant les
10 derniers nucléotides. Les reads feront une taille d’à peu près 240pb.

``` r
plotQualityProfile(fnRs[1:2])
```

![](02_dada-analysist_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->
\#Le deuxième plot montre le score qualité des reads reverse. On
remarque qu’ils sont moins bon, et donc la partie tronquée sera moins
longue. La partie prise en compte alors sera que les 160 premières pb.

## 3- Filtration des données

\#On construit deux nouveaux objets filtFs et filtRs qui sont les noms
des fichiers contenant les données filtrées

``` r
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

\#La filtration des données est basée sur les paramètres générées à
partir des scores de qualité en plus de paramétres de filtrage standards
de DADA2.

``` r
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
```

    ##                               reads.in reads.out
    ## F3D0_S188_L001_R1_001.fastq       7793      7113
    ## F3D1_S189_L001_R1_001.fastq       5869      5299
    ## F3D141_S207_L001_R1_001.fastq     5958      5463
    ## F3D142_S208_L001_R1_001.fastq     3183      2914
    ## F3D143_S209_L001_R1_001.fastq     3178      2941
    ## F3D144_S210_L001_R1_001.fastq     4827      4312

# III- Modèle d’erreur

\#DADA2 utilise un modèle d’erreur paramétrique pour différencier les
erreurs de séquençage des variations réelles de nucléotides dans les
séquences.Les paramètres de ce modèle peuvent être générés à partir des
données elles-mêmes en utilisant une forme d’appréhendement non
supervisée dans laquelle l’inférence d’échantillon est alternée avec
l’estimation des paramètres jusqu’à ce que les deux soient cohérents.

``` r
errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## 33514080 total bases in 139642 reads from 20 samples will be used for learning the error rates.

``` r
errR <- learnErrors(filtRs, multithread=TRUE)
```

    ## 22342720 total bases in 139642 reads from 20 samples will be used for learning the error rates.

``` r
plotErrors(errF, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis
    
    ## Warning: Transformation introduced infinite values in continuous y-axis

![](02_dada-analysist_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->
\#Les taux d’erreur pour chaque transition possible (A → C, A → G,…)
sont indiqués. Les points sont les taux d’erreur observés pour chaque
score de qualité consensuel. La ligne noire montre les taux d’erreur
estimés après convergence de l’algorithme. La ligne rouge montre les
taux d’erreur attendus selon la définition nominale du Q-score. Ici, les
taux d’erreur estimés (ligne noire) correspondent bien aux taux observés
(points), et les taux d’erreur diminuent avec une qualité accrue comme
prévu. Tout semble être raisonnable.

# IV- Inférence d’échantillon

\#Après avoir filtré les données, et construit le modèle d’erreurs les
codes qui suivent visent à appliquer le modèle d’erreur obtenu
précedemment afin de corriger les erreurs sur les séquences (Sens et
anti-sens). L’algorithme est issu de la fonction dada (d’où dadaFs et
dadaRs pour forward et reverse).

``` r
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
```

    ## Sample 1 - 7113 reads in 1979 unique sequences.
    ## Sample 2 - 5299 reads in 1639 unique sequences.
    ## Sample 3 - 5463 reads in 1477 unique sequences.
    ## Sample 4 - 2914 reads in 904 unique sequences.
    ## Sample 5 - 2941 reads in 939 unique sequences.
    ## Sample 6 - 4312 reads in 1267 unique sequences.
    ## Sample 7 - 6741 reads in 1756 unique sequences.
    ## Sample 8 - 4560 reads in 1438 unique sequences.
    ## Sample 9 - 15637 reads in 3590 unique sequences.
    ## Sample 10 - 11413 reads in 2762 unique sequences.
    ## Sample 11 - 12017 reads in 3021 unique sequences.
    ## Sample 12 - 5032 reads in 1566 unique sequences.
    ## Sample 13 - 18075 reads in 3707 unique sequences.
    ## Sample 14 - 6250 reads in 1479 unique sequences.
    ## Sample 15 - 4052 reads in 1195 unique sequences.
    ## Sample 16 - 7369 reads in 1832 unique sequences.
    ## Sample 17 - 4765 reads in 1183 unique sequences.
    ## Sample 18 - 4871 reads in 1382 unique sequences.
    ## Sample 19 - 6504 reads in 1709 unique sequences.
    ## Sample 20 - 4314 reads in 897 unique sequences.

``` r
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
```

    ## Sample 1 - 7113 reads in 1660 unique sequences.
    ## Sample 2 - 5299 reads in 1349 unique sequences.
    ## Sample 3 - 5463 reads in 1335 unique sequences.
    ## Sample 4 - 2914 reads in 853 unique sequences.
    ## Sample 5 - 2941 reads in 880 unique sequences.
    ## Sample 6 - 4312 reads in 1286 unique sequences.
    ## Sample 7 - 6741 reads in 1803 unique sequences.
    ## Sample 8 - 4560 reads in 1265 unique sequences.
    ## Sample 9 - 15637 reads in 3414 unique sequences.
    ## Sample 10 - 11413 reads in 2522 unique sequences.
    ## Sample 11 - 12017 reads in 2771 unique sequences.
    ## Sample 12 - 5032 reads in 1415 unique sequences.
    ## Sample 13 - 18075 reads in 3290 unique sequences.
    ## Sample 14 - 6250 reads in 1390 unique sequences.
    ## Sample 15 - 4052 reads in 1134 unique sequences.
    ## Sample 16 - 7369 reads in 1635 unique sequences.
    ## Sample 17 - 4765 reads in 1084 unique sequences.
    ## Sample 18 - 4871 reads in 1161 unique sequences.
    ## Sample 19 - 6504 reads in 1502 unique sequences.
    ## Sample 20 - 4314 reads in 732 unique sequences.

``` r
dadaFs[[1]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 128 sequence variants were inferred from 1979 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

``` r
dadaFs[[2]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 113 sequence variants were inferred from 1639 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

``` r
dadaFs[[3]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 97 sequence variants were inferred from 1477 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

\#Les exemples donnés ci-dessus sont pour les échantillons 1, 2 et 3
respectivement. Pour l’echantillon 3 par exemple, DADA2 a trouvé que
qu’il y’avait 97 variants dans 1477 séquences uniques. \#NB: 1477
séquences uniques ne veut pas forcément dire 1477 bactéries, car
certaines bactéries peuvent avoir plusieurs opéron du gène de l’ARN16s.

# V- Construction de la table d’observation

## 1- Construction des contigues

\#Il s’agit ici de construire des contigues à partir des reads R1 er R2
en faisant un over lap. Ceci est possible dans ce cas car la taille des
reads permet leur chevauchement (R1: 240pb et R2: 160pb). Cette
opération aurait été impossible si le séquençage pourtait sur les
régions V3 ou V5 qui font 500 à 600pb.

``` r
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```

    ## 6551 paired-reads (in 106 unique pairings) successfully merged out of 6907 (in 199 pairings) input.

    ## 5025 paired-reads (in 100 unique pairings) successfully merged out of 5188 (in 156 pairings) input.

    ## 4973 paired-reads (in 80 unique pairings) successfully merged out of 5268 (in 166 pairings) input.

    ## 2595 paired-reads (in 52 unique pairings) successfully merged out of 2756 (in 109 pairings) input.

    ## 2553 paired-reads (in 60 unique pairings) successfully merged out of 2785 (in 119 pairings) input.

    ## 3622 paired-reads (in 53 unique pairings) successfully merged out of 4103 (in 157 pairings) input.

    ## 6079 paired-reads (in 81 unique pairings) successfully merged out of 6515 (in 198 pairings) input.

    ## 3961 paired-reads (in 90 unique pairings) successfully merged out of 4384 (in 188 pairings) input.

    ## 14231 paired-reads (in 143 unique pairings) successfully merged out of 15358 (in 351 pairings) input.

    ## 10526 paired-reads (in 120 unique pairings) successfully merged out of 11166 (in 279 pairings) input.

    ## 11156 paired-reads (in 137 unique pairings) successfully merged out of 11799 (in 298 pairings) input.

    ## 4329 paired-reads (in 84 unique pairings) successfully merged out of 4788 (in 180 pairings) input.

    ## 17431 paired-reads (in 153 unique pairings) successfully merged out of 17812 (in 272 pairings) input.

    ## 5850 paired-reads (in 81 unique pairings) successfully merged out of 6095 (in 159 pairings) input.

    ## 3716 paired-reads (in 86 unique pairings) successfully merged out of 3894 (in 147 pairings) input.

    ## 6865 paired-reads (in 99 unique pairings) successfully merged out of 7193 (in 187 pairings) input.

    ## 4430 paired-reads (in 67 unique pairings) successfully merged out of 4605 (in 127 pairings) input.

    ## 4574 paired-reads (in 100 unique pairings) successfully merged out of 4736 (in 172 pairings) input.

    ## 6094 paired-reads (in 109 unique pairings) successfully merged out of 6314 (in 172 pairings) input.

    ## 4269 paired-reads (in 20 unique pairings) successfully merged out of 4281 (in 28 pairings) input.

\#Le tableau ci-dessus montre les séquences chevauchées (Sens et
anti-sens).

``` r
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
```

    ##                                                                                                                                                                                                                                                       sequence
    ## 1 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCAAGTCAGCGGTAAAATTGAGAGGCTCAACCTCTTCGAGCCGTTGAAACTGGTTTTCTTGAGTGAGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCTCAACTGACGCTCATGCACGAAAGTGTGGGTATCGAACAGG
    ## 2 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCCAAGTCAGCGGTAAAATTGCGGGGCTCAACCCCGTACAGCCGTTGAAACTGCCGGGCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGAAAGTGCGGGGATCAAACAGG
    ## 3 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGCTGTTAAGTCAGCGGTCAAATGTCGGGGCTCAACCCCGGCCTGCCGTTGAAACTGGCGGCCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCCCGACTGACGCTGAGGCACGAAAGCGTGGGTATCGAACAGG
    ## 4 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGCTTTTAAGTCAGCGGTAAAAATTCGGGGCTCAACCCCGTCCGGCCGTTGAAACTGGGGGCCTTGAGTGGGCGAGAAGAAGGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCCTTCCGGCGCCCTACTGACGCTGAGGCACGAAAGTGCGGGGATCGAACAGG
    ## 5 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGACTCTCAAGTCAGCGGTCAAATCGCGGGGCTCAACCCCGTTCCGCCGTTGAAACTGGGAGCCTTGAGTGCGCGAGAAGTAGGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCCTACCGGCGCGCAACTGACGCTCATGCACGAAAGCGTGGGTATCGAACAGG
    ## 6 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGATGCCAAGTCAGCGGTAAAAAAGCGGTGCTCAACGCCGTCGAGCCGTTGAAACTGGCGTTCTTGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGAAAGCGTGGGTATCGAACAGG
    ##   abundance forward reverse nmatch nmismatch nindel prefer accept
    ## 1       579       1       1    148         0      0      1   TRUE
    ## 2       470       2       2    148         0      0      2   TRUE
    ## 3       449       3       4    148         0      0      1   TRUE
    ## 4       430       4       3    148         0      0      2   TRUE
    ## 5       345       5       6    148         0      0      1   TRUE
    ## 6       282       6       5    148         0      0      2   TRUE

\#Dans le second tableau, la colonne nmatch montre la longueure de la
taille de l’overlaps. la colonne Prefer montre quelle séquence est prise
entre R1 et R2 s’il y a des mismatshing au niveau des chevauchements.
Pour choisir, l’ordinateur va regarder le meilleur Q score entre R1 et
R2. Et puis toutes les séquences incapables de former les contigues on
été supprimées.

## 2- Construction de la table d’observations

\#Il s’agit de construire la table d’observations des ASV par
échantillons à partir des contigues emergés.

``` r
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```

    ## [1]  20 293

``` r
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
```

    ## 
    ## 251 252 253 254 255 
    ##   1  88 196   6   2

\#La table contient 293 ASV, et les longueurs des contigues, se situent
toutes autour de 253pb., ça veut dire dans la plage attendue pour la
région V4.

## 3- Elimination des chimères

\#Les chimères sont produits pendant l’étape de la PCR. Mais
généralement, ce sont des éléments rares. DADA2 va aussi se baser sur
cette caractéristique pour les éliminer. Mais ce qui va lui permettre de
les détecter, c’est le fait qu’ils vont avoir un début de séquence qui
peut être apparié avec une séquence abondante et une fin de séquence qui
peut s’apparier avec une autre séquence abondante.

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
```

    ## Identified 61 bimeras out of 293 input sequences.

``` r
dim(seqtab.nochim)
```

    ## [1]  20 232

``` r
sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 0.964263

\#Il y a donc 3.5% de séquences chimériques dans le jeu de données.

## 4- Suivi des reads dans le pipeline

\#Il s’agit de faire une sorte de recap des filtrations appliquées en
créant une nouvelle fonction getN, afin de vérifier l’avancement du
travail avant procéder à l’assignation taxonomique.

``` r
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```

    ##        input filtered denoisedF denoisedR merged nonchim
    ## F3D0    7793     7113      6996      6978   6551    6539
    ## F3D1    5869     5299      5227      5239   5025    5014
    ## F3D141  5958     5463      5339      5351   4973    4850
    ## F3D142  3183     2914      2799      2833   2595    2521
    ## F3D143  3178     2941      2822      2868   2553    2519
    ## F3D144  4827     4312      4146      4224   3622    3483

# V- Assignation taxonomique

\#Pour fair l’assignation taxonomique, il va nous falloir une base de
données qui est dans notre SILVA, puis un algorithme d’assignation qui
va être NBG classifier, et qui va attribuer à chaque sequence le rang
taxonomique le plus loin possible.

``` r
taxa <- assignTaxonomy(seqtab.nochim,"~/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
```

``` r
taxa<- addSpecies(taxa, "~/silva_species_assignment_v138.fa.gz")
```

``` r
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
```

    ##      Kingdom    Phylum         Class         Order           Family          
    ## [1,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [2,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [3,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [4,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [5,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Bacteroidaceae"
    ## [6,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ##      Genus         Species
    ## [1,] NA            NA     
    ## [2,] NA            NA     
    ## [3,] NA            NA     
    ## [4,] NA            NA     
    ## [5,] "Bacteroides" NA     
    ## [6,] NA            NA

\#Comme prévu, les Bacteroidetes sont bien représentés parmi les taxons
les plus abondants dans ces échantillons fécaux. Peu d’attributions
d’espèces ont été faites. Cela est expliqué par d’une part la
difficulté de faire une assignation d’espèces de façon assez sure à
partir du gène de l’ARN16S, et d’autre part du fait que les données du
microbiote intestinal de souris indigène dans les bases de données de
référence sont peu présentes.

# VI- Evaluation de l’efficacité

\#Afin de contrôler la qualité du séquençage et des analyses effectuées,
une communauté fictive de 20 échantillons, dont les séquences sont
connues à priori, a été ajoutée au jeu de données. Nous revenons donc à
ces échantillons et comparons les variantes de séquence inférées par
DADA2 à la composition attendue de la communauté.

``` r
unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
```

    ## DADA2 inferred 20 sample sequences present in the Mock community.

``` r
mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
```

    ## Of those, 20 were exact matches to the expected reference sequences.

\#DADA2 a identifié 20 ASV qui correspondent exactement aux génomes de
référence des membres attendus de la communauté. Le taux d’erreur
résiduel après le pipeline DADA2 pour cet échantillon est de 0%.

\#NB: Lecode ci-dessous ne fait pas partie de l’analyse DADA2, mais il
sert pour enregistrer le travail qui a été fait pour DADA2 afin de le
retrouver et l’utiliser directement pour la partie phyloseq.

``` r
save.image("~/phyloseq_tuto/02_dada-analysist")
```
