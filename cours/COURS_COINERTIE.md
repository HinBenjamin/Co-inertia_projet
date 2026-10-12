# Cours : L'Analyse de Co-inertie (CoIA)

> **Auteurs de référence :** S. Dolédec & D. Chessel (1994), S. Dray, D. Chessel & J. Thioulouse (2003).
> **Cadre méthodologique :** Diagrammes de dualité, package `ade4` (R).

---

## 1. Contexte et Problématique

En écologie et en biologie (ex. diversité microbienne, écologie des communautés, métagénomique), on dispose très souvent de **deux tableaux de données mesurés sur les mêmes $n$ unités d'échantillonnage** (sites, individus, prélèvements) :

* Un tableau $X$ de dimensions $(n \times p)$ : ex. variables environnementales, physico-chimiques, physico-chimie du sol.
* Un tableau $Y$ de dimensions $(n \times q)$ : ex. abondances taxonomiques, OTU/ASV, fonctions écologiques.

### L'objectif du couplage :

Analyser simultanément les deux tableaux pour comprendre comment la structure de variabilité de $X$ est liée à celle de $Y$.

---

## 2. Fondements Théoriques : Les Triplets Statistiques et l'Inertie

Dans le formalisme des **diagrammes de dualité** (Escoufier 1987), chaque tableau est associé à un **triplet statistique** :

* $(X, Q_X, D_n)$ pour le premier tableau.
* $(Y, Q_Y, D_n)$ pour le second tableau.

Où :

* $D_n$ ($n \times n$) est la matrice diagonale des **poids des individus** (généralement $D_n = \frac{1}{n} I_n$, ou les poids marginaux d'une AFC). **Important :** les deux analyses doivent partager la même métrique des individus $D_n$.
* $Q_X$ ($p \times p$) et $Q_Y$ ($q \times q$) sont les **métriques associées aux variables** (ex. identité pour une ACP sur covariances, inverse des variances pour une ACP normée).

### 2.1. Inertie d'un nuage de points

L'inertie mesure la dispersion globale du nuage par rapport à son centre de gravité :

$$
I(X) = \operatorname{Tr}(X Q_X X^T D_n) = \sum_{k=1}^r \lambda_k
$$

> **Remarque clé :** L'inertie mesure la *taille globale* (la quantité d'information) du nuage, mais **pas sa forme** ni sa direction principale d'élongation.

### 2.2. Variance et Covariance Vectorielles

* **Variance vectorielle (Escoufier) :**

  $$
  \operatorname{Vv}(X) = \operatorname{Tr}\left((X Q_X X^T D_n)^2\right)
  $$

  Elle mesure la structure interne et la forme du nuage.
* **Covariance vectorielle :**

  $$
  \operatorname{Covv}(X, Y) = \operatorname{Tr}\left(X Q_X X^T D_n Y Q_Y Y^T D_n\right)
  $$

  Elle mesure la liaison globale entre les deux structures de données.

### 2.3. Le Coefficient RV (Escoufier 1973)

Le coefficient $RV$ est l'analogue multivarié d'un coefficient de corrélation au carré ($r^2$) :

$$
RV(X, Y) = \frac{\operatorname{Covv}(X, Y)}{\sqrt{\operatorname{Vv}(X) \cdot \operatorname{Vv}(Y)}} \in [0, 1]
$$

* $RV = 0$ : aucune co-structure entre les deux tableaux.
* $RV = 1$ : concordance géométrique parfaite entre les deux nuages.
* **Test de permutation (Monte-Carlo) :** Dans `ade4`, on permute aléatoirement les lignes de l'un des tableaux pour tester si la liaison observée est significativement supérieure au hasard (`randtest`).

---

## 3. Principe de l'Analyse de Co-inertie (CoIA)

### 3.1. Formulation du critère d'optimisation

L'analyse de co-inertie cherche simultanément deux axes normés $u \in \mathbb{R}^p$ ($u^T Q_X u = 1$) et $v \in \mathbb{R}^q$ ($v^T Q_Y v = 1$) tels que la **covariance au carré** entre les projections $f = X Q_X u$ et $g = Y Q_Y v$ soit **maximale** :

$$
\max_{u, v} \operatorname{Cov}_{D_n}^2(f, g)
$$

Or, la covariance se décompose en :

$$
\operatorname{Cov}^2(f, g) = \operatorname{Var}(f) \times \operatorname{Var}(g) \times \operatorname{Cor}^2(f, g)
$$

### 3.2. Le compromis fondamental de la co-inertie

La co-inertie réalise un compromis optimal à trois composantes :

1. **$\operatorname{Cor}^2(f, g)$ :** Les scores sur les deux tableaux doivent être fortement corrélés.
2. **$\operatorname{Var}(f)$ :** L'axe doit expliquer une part importante de l'inertie du tableau $X$.
3. **$\operatorname{Var}(g)$ :** L'axe doit expliquer une part importante de l'inertie du tableau $Y$.

---

## 4. Comparaison avec les autres méthodes de couplage

| Critère                                     | Analyse de Co-inertie (CoIA)                                       | Analyse Canonique des Corrélations (ACC) | Analyse Redondante (RDA / ACPVI)              |
| :------------------------------------------- | :----------------------------------------------------------------- | :---------------------------------------- | :-------------------------------------------- |
| **Objectif**                           | Maximiser la**covariance**                                   | Maximiser la**corrélation**        | Expliquer$Y$ par une régression sur $X$  |
| **Rôle des tableaux**                 | **Symétrique** (co-structure)                               | **Symétrique**                     | **Asymétrique** ($X$ explique $Y$) |
| **Contrainte de dimension**            | **Aucune** ($p+q \gg n$ autorisé) | Forte ($n > p + q$) | $n > p$ obligatoire                     |                                               |
| **Inversion de matrice**               | **Non** (décomposition spectrale directe)                   | **Oui** (instable si colinéarité) | **Oui** (inversion de $X^T D_n X$)    |
| **Flexibilité des types de données** | **Maximale** (mixe ACP, AFC, ACM, Hill-Smith...)             | Variables quantitatives continues         | Principalement quantitatif (ou contraint)     |

> **Force majeure :** La co-inertie est idéale lorsque $p \gg n$ ou $q \gg n$ (ex. centaines d'OTUs ou de gènes sur quelques dizaines de stations/échantillons).

---

## 5. Mise en œuvre pratique dans R (`ade4`)

### Exemple type : Couplage Environnement (ACP) et Espèces (AFC)

```r
library(ade4)

# 1. Analyse séparée du tableau X (ex: physico-chimie du sol)
# ACP centrée et réduite
pca_env <- dudi.pca(df_env, row.w = rep(1/nrow(df_env), nrow(df_env)), 
                    scannf = FALSE, nf = 2)

# 2. Analyse séparée du tableau Y (ex: abondances d'espèces/OTUs)
# AFC (Analyse Factorielle des Correspondances)
coa_spe <- dudi.coa(df_spe, scannf = FALSE, nf = 2)

# Important : pour coupler ACP et AFC, on utilise souvent le poids des stations issu de l'AFC :
pca_env <- dudi.pca(df_env, row.w = coa_spe$lw, scannf = FALSE, nf = 2)

# 3. Réalisation de la Co-inertie
coin <- coinertia(pca_env, coa_spe, scannf = FALSE, nf = 2)

# 4. Résumé et test de significativité (Monte-Carlo)
summary(coin)
test_rv <- randtest(coin, nrepet = 999)
plot(test_rv)

# 5. Visualisations
plot(coin)
```

### Éléments clés de l'objet `coinertia` :

* `$eig` : valeurs propres (co-inertie associée à chaque dimension).
* `$RV` : valeur observée du coefficient RV.
* `$l1`, `$lX` : coordonnées des individus selon le tableau $X$.
* `$m1`, `$mY` : coordonnées des individus selon le tableau $Y$.
* `$c1` : projections des variables du premier tableau.
* `$li1` : projections des variables du second tableau.

### Représentation graphique :

* **Graphe des flèches (`plot(coin)`) :** Chaque individu est représenté par deux points reliés par une flèche. La tête de la flèche correspond à la position selon un tableau, l'origine selon l'autre tableau.
  * *Flèches courtes* : forte concordance des deux tableaux pour ces sites.
  * *Flèches longues* : divergence / désaccord entre la composition biologique et les conditions environnementales.

---

## 6. Références Académiques Principales

1. **Dolédec, S., & Chessel, D. (1994).** *Co-inertia analysis: an alternative method for studying species-environment relationships.* Freshwater Biology, 31(3), 277-294.
2. **Dray, S., Chessel, D., & Thioulouse, J. (2003).** *Co-inertia analysis and the linking of ecological data tables.* Ecology, 84(11), 3078-3089.
3. **Chessel, D., & Dufour, A.-B.** *Notes de cours cssb8 : Co-inertie, co-structures et compromis.* Master Écologie, Évolution, Biométrie, Université Lyon 1.
4. **Culhane, A. C., Perrière, G., & Higgins, D. G. (2003).** *Cross-platform comparison and visualisation of gene expression data using co-inertia analysis.* BMC Bioinformatics, 4:59.
