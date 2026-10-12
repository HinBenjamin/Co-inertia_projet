# Manuel Complet et Exhaustif : L'Analyse de Co-inertie (CoIA)

> **Destiné à l'étude approfondie, la préparation d'examens et la recherche (M2 / Doctorat).**
> **Auteurs fondateurs :** Yves Escoufier (Diagrammes de dualité, 1973, 1987), Daniel Chessel & Sylvain Dolédec (1994), Stéphane Dray, Daniel Chessel & Jean Thioulouse (2003).

---

## Sommaire

1. [Introduction Intuitive &amp; Problématique](#1-introduction-intuitive--problématique)
2. [Cadre Mathématique : Les Diagrammes de Dualité d&#39;Escoufier](#2-cadre-mathématique--les-diagrammes-de-dualité-descoufier)
3. [Inertie, Variance Vectorielle et Covariance Vectorielle](#3-inertie-variance-vectorielle-et-covariance-vectorielle)
4. [Le Coefficient RV et le Test de Permutation de Monte-Carlo](#4-le-coefficient-rv-et-le-test-de-permutation-de-monte-carlo)
5. [La Théorie Algébrique de l&#39;Analyse de Co-inertie (CoIA)](#5-la-théorie-algébrique-de-lanalyse-de-co-inertie-coia)
6. [Le Compromis Fondamental de la Co-inertie](#6-le-compromis-fondamental-de-la-co-inertie)
7. [Comparatif Exhaustif : CoIA vs ACC vs RDA vs Procuste vs PLS](#7-comparatif-exhaustif--coia-vs-acc-vs-rda-vs-procuste-vs-pls)
8. [Typologie des Couplages et Gestion des Métriques](#8-typologie-des-couplages-et-gestion-des-métriques)
9. [Anatomie Complète d&#39;un Objet `coinertia` sous R (`ade4`)](#9-anatomie-complète-dun-objet-coinertia-sous-r-ade4)
10. [Interprétation des Sorties Graphiques](#10-interprétation-des-sorties-graphiques)
11. [Généralisation à $K$ Tableaux : Co-inertie Multiple &amp; STATIS](#11-généralisation-à-k-tableaux--co-inertie-multiple--statis)
12. [Script R Complet et Cas Pratique Écologique](#12-script-r-complet-et-cas-pratique-écologique)
13. [Foire Aux Questions, Pièges et Questions d&#39;Examen](#13-foire-aux-questions-pièges-et-questions-dexamen)

---

## 1. Introduction Intuitive & Problématique

### 1.1. Le problème des deux tableaux

En sciences expérimentales, et en particulier en écologie, en microbiologie et en génomique, nous sommes constamment confrontés au couplage de deux points de vue différents sur les mêmes unités d'échantillonnage :

* **Tableau $X$ ($n$ lignes $\times p$ colonnes) :** Mesures physico-chimiques, variables environnementales, paramètres de sol, pH, humidité, nutriments.
* **Tableau $Y$ ($n$ lignes $\times q$ colonnes) :** Abondances d'espèces, comptages d'OTUs/ASVs bactériens, fonctions métaboliques, expressions de gènes.

Les deux tableaux partagent obligatoirement **les mêmes $n$ individus ou stations d'échantillonnage**.

```
       Tableau X (Milieu)               Tableau Y (Communauté)
            p variables                       q espèces
        ┌─────────────────┐               ┌─────────────────┐
      1 │                 │             1 │                 │
      2 │                 │             2 │                 │
 n  ... │      X          │        n  ... │      Y          │
sites   │                 │      sites    │                 │
      n │                 │             n │                 │
        └─────────────────┘               └─────────────────┘
                 ▲                                 ▲
                 └─────────── Mêmes sites ─────────┘
```

### 1.2. Pourquoi les approches naïves échouent-elles ?

1. **Faire une corrélation variable par variable (Pearson/Spearman) :** Si $p = 20$ variables physico-chimiques et $q = 500$ espèces, on obtient $20 \times 500 = 10\,000$ corrélations ! C'est le piège des tests multiples, sans vision synthétique de la structure globale.
2. **Faire une régression multiple colonne par colonne :** En présence de colinéarité entre variables ou dès que $p > n$, l'inversion de matrice devient impossible (problème d'inversibilité de $X^T X$).
3. **Faire deux ACP indépendantes et comparer les plans à l'œil :** Rien ne garantit que le premier axe de l'ACP de $X$ (qui maximise la variance de $X$) soit lié au premier axe de l'ACP de $Y$. Ils peuvent capturer des phénomènes totalement indépendants.

**L'Analyse de Co-inertie (CoIA)** a été créée par Chessel, Dolédec, Dray et Thioulouse précisément pour résoudre cette question : **trouver un système d'axes coordonnés qui extrait simultanément la variabilité maximale de chaque tableau tout en maximisant leur corrélation.**

---

## 2. Cadre Mathématique : Les Diagrammes de Dualité d'Escoufier

Pour comprendre la co-inertie telle qu'elle est formalisée à l'Université de Lyon et dans le package `ade4`, il est indispensable de maîtriser le formalisme unifié des **triplets statistiques** développé par Yves Escoufier (1987).

Toute analyse factorielle linéaire (ACP, AFC, ACM...) est entièrement définie par un **triplet statistique** :

$$
(X, Q_X, D_n)
$$

1. **$X$ ($n \times p$) :** Le tableau des données (centré, réduit, transformé ou issu de profils).
2. **$Q_X$ ($p \times p$) :** La métrique des variables (définissant la distance dans $\mathbb{R}^n$).
   * Dans une ACP non normée : $Q_X = I_p$ (matrice identité).
   * Dans une ACP normée (sur corrélations) : $Q_X = \operatorname{diag}(1 / s_j^2)$.
3. **$D_n$ ($n \times n$) :** La matrice diagonale des **poids des individus** (définissant la distance dans $\mathbb{R}^p$).
   * Pondération uniforme standard : $D_n = \frac{1}{n} I_n$.
   * En AFC (Analyse Factorielle des Correspondances) : $D_n = \operatorname{diag}(p_{i\cdot})$ (poids marginaux des lignes).

### L'opérateur d'inertie

À tout triplet $(X, Q_X, D_n)$, on associe l'opérateur linéaire d'individus dans $\mathbb{R}^n$ :

$$
W_X = X Q_X X^T D_n
$$

Cet opérateur $W_X$ est symétrique par rapport à la métrique $D_n$. Il contient **toute l'information géométrique** du nuage des $n$ points dans l'espace des variables : distances inter-individus et orientations.

---

## 3. Inertie, Variance Vectorielle et Covariance Vectorielle

L'espace vectoriel des opérateurs d'individus $\mathcal{M}_n$ peut être muni du produit scalaire de Hilbert-Schmidt :

$$
\langle A, B \rangle = \operatorname{Tr}(A D_n B D_n)
$$

Ce produit scalaire permet de mesurer la taille, la structure et la ressemblance entre tableaux de manière rigoureuse :

### 3.1. Inertie totale d'un nuage

L'inertie globale est la trace de l'opérateur :

$$
\operatorname{Iner}(X) = \operatorname{Tr}(W_X) = \operatorname{Tr}(X Q_X X^T D_n) = \sum_{k=1}^{\operatorname{rang}(X)} \lambda_k
$$

* **Propriété géométrique :** L'inertie mesure la **dispersion totale** du nuage de points autour de son centre de gravité.
* **Limite fondamentale :** L'inertie mesure la *taille* du nuage, mais **absolument pas sa forme** (un nuage sphérique isotrope et un nuage très allongé sur une ligne peuvent avoir exactement la même inertie globale).

### 3.2. Variance vectorielle (forme du nuage)

La variance vectorielle est la norme au carré de l'opérateur $W_X$ dans l'espace des opérateurs :

$$
\operatorname{Vv}(X) = \|W_X\|^2 = \operatorname{Tr}(W_X^2) = \operatorname{Tr}\left((X Q_X X^T D_n)^2\right) = \sum_{k=1}^r \lambda_k^2
$$

* **Interprétation :** La variance vectorielle mesure l'anisotropie et l'intensité de la structure.
  * Si $\operatorname{Vv}(X) \approx (\operatorname{Iner}(X))^2$, un seul axe concentre toute l'information (nuage filiforme unidimensionnel).
  * Si $\operatorname{Vv}(X)$ est faible par rapport à l'inertie au carré, la dispersion est éparpillée sur de nombreuses dimensions équivalentes (bruit isotrope).

### 3.3. Covariance vectorielle (liaison entre deux tableaux)

Soient deux tableaux mesurés sur les mêmes $n$ individus, représentés par $(X, Q_X, D_n)$ et $(Y, Q_Y, D_n)$.
**Condition absolue :** Les deux analyses séparées doivent impérativement avoir la **même métrique des individus $D_n$**.

La covariance vectorielle est le produit scalaire entre les deux opérateurs d'inertie :

$$
\operatorname{Covv}(X, Y) = \langle W_X, W_Y \rangle = \operatorname{Tr}\left(W_X D_n W_Y D_n\right) = \operatorname{Tr}\left(X Q_X X^T D_n Y Q_Y Y^T D_n\right)
$$

* **Interprétation :** $\operatorname{Covv}(X, Y)$ quantifie le recouvrement spatial et la concordance globale des distances entre les individus mesurées dans $X$ et dans $Y$.

---

## 4. Le Coefficient RV et le Test de Permutation de Monte-Carlo

### 4.1. Définition du coefficient RV d'Escoufier (1973)

Par analogie avec le coefficient de corrélation au carré classique entre deux variables scalaires :

$$
\operatorname{Cor}^2(x, y) = \frac{\operatorname{Cov}^2(x, y)}{\operatorname{Var}(x) \operatorname{Var}(y)}
$$

Escoufier a défini le **coefficient RV** entre deux tableaux comme le cosinus au carré entre les opérateurs $W_X$ et $W_Y$ :

$$
RV(X, Y) = \frac{\operatorname{Covv}(X, Y)}{\sqrt{\operatorname{Vv}(X) \cdot \operatorname{Vv}(Y)}} = \frac{\operatorname{Tr}\left(X Q_X X^T D_n Y Q_Y Y^T D_n\right)}{\sqrt{\operatorname{Tr}\left((X Q_X X^T D_n)^2\right) \cdot \operatorname{Tr}\left((Y Q_Y Y^T D_n)^2\right)}}
$$

### 4.2. Propriétés mathématiques du coefficient RV

1. **Bornes :** $0 \le RV(X, Y) \le 1$.
2. **Invariance :** $RV$ est invariant par homothétie et rotation orthogonale des nuages.
3. **Cas extrêmes :**
   * $RV = 0$ : Aucune liaison linéaire entre les structures géométriques de $X$ et $Y$.
   * $RV = 1$ : Les deux nuages de points sont homotopes (à une rotation et un facteur d'échelle près, la configuration spatiale des $n$ individus est identique).

### 4.3. Test de permutation de Monte-Carlo (`randtest`)

Comme les variables sont multivariées et corrélées, on ne peut pas appliquer un test $t$ de Student ou de Fisher conventionnel. On utilise donc un **test de permutation non-paramétrique** :

* **Hypothèse nulle ($H_0$) :** Le tableau $X$ et le tableau $Y$ sont indépendants (les lignes de $X$ n'ont aucun lien spécifique avec les lignes de $Y$).
* **Procédure :**
  1. Calculer le $RV_{\text{obs}}$ réel.
  2. Permuter aléatoirement l'ordre des lignes de $Y$ (ou de $X$) tout en conservant $D_n$.
  3. Recalculer le coefficient $RV_{\text{sim}}$ sur les données permutées.
  4. Répéter l'opération $K$ fois (ex. $K = 999$ permutations).
  5. Calculer la p-value empirique :
     $$
     p = \frac{\#\{RV_{\text{sim}} \ge RV_{\text{obs}}\} + 1}{K + 1}
     $$
* Si $p < 0.05$, on rejette $H_0$ : la co-structure entre les deux tableaux est statistiquement significative.

---

## 5. La Théorie Algébrique de l'Analyse de Co-inertie (CoIA)

Alors que le coefficient $RV$ fournit un résumé scalaire global de la liaison, l'analyse de co-inertie calcule un **système complet d'axes de projection orthogonaux** pour visualiser et décomposer cette co-structure.

### 5.1. Formulation de l'optimisation

On cherche un vecteur unitaire $u_1 \in \mathbb{R}^p$ (axe pour les variables de $X$) et un vecteur unitaire $v_1 \in \mathbb{R}^q$ (axe pour les variables de $Y$) tels que :

$$
\max_{u, v} \operatorname{Cov}_{D_n}^2\left(X Q_X u, Y Q_Y v\right)
$$

Sous les contraintes de normalisation métrique :

$$
u^T Q_X u = 1 \quad \text{et} \quad v^T Q_Y v = 1
$$

### 5.2. Résolution algébrique (Lagrangien et SVD)

Posons la fonction de Lagrange :

$$
\mathcal{L}(u, v, \alpha, \beta) = \left(u^T Q_X X^T D_n Y Q_Y v\right)^2 - \alpha (u^T Q_X u - 1) - \beta (v^T Q_Y v - 1)
$$

En dérivant par rapport à $u$ et $v$ et en annulant les dérivées partielles, on montre que ce problème équivaut à la **Décomposition en Valeurs Singulières (SVD)** ou à la diagonalisation de l'opérateur de co-inertie :

$$
W_{XY} = X^T D_n Y Q_Y Y^T D_n X Q_X
$$

Pour le $k$-ième axe factoriel :

* Le vecteur propre $u_k$ vérifie :
  $$
  (X^T D_n Y Q_Y Y^T D_n X) Q_X u_k = \lambda_k u_k
  $$
* Le vecteur propre $v_k$ est obtenu par dualité :
  $$
  v_k = \frac{1}{\sqrt{\lambda_k}} Y^T D_n X Q_X u_k
  $$
* La valeur propre $\lambda_k$ est la **co-inertie** associée au $k$-ième axe :
  $$
  \lambda_k = \operatorname{Cov}_{D_n}^2(f_k, g_k)
  $$

  Où $f_k = X Q_X u_k$ est le score des individus selon $X$, et $g_k = Y Q_Y v_k$ est le score des individus selon $Y$.

### 5.3. Décomposition de la co-inertie totale

La somme de toutes les valeurs propres de co-inertie est égale à la covariance vectorielle totale :

$$
\sum_{k=1}^{\min(p, q)} \lambda_k = \operatorname{Covv}(X, Y) = \operatorname{Tr}\left(X Q_X X^T D_n Y Q_Y Y^T D_n\right)
$$

Le pourcentage d'inertie expliqué par le premier axe s'écrit donc :

$$
\% \operatorname{CoIA}_1 = \frac{\lambda_1}{\sum_k \lambda_k} \times 100
$$

---

## 6. Le Compromis Fondamental de la Co-inertie

C'est le point central de la méthode, souvent interrogé en examen.

Par définition même de la covariance statistique :

$$
\operatorname{Cov}(f_k, g_k) = \operatorname{Var}^{1/2}(f_k) \times \operatorname{Var}^{1/2}(g_k) \times \operatorname{Cor}(f_k, g_k)
$$

En élevant au carré, la valeur propre de co-inertie $\lambda_k$ se factorise rigoureusement en :

$$
\lambda_k = \underbrace{\operatorname{Var}_{D_n}(f_k)}_{\text{Inertie projetée de } X} \times \underbrace{\operatorname{Var}_{D_n}(g_k)}_{\text{Inertie projetée de } Y} \times \underbrace{\operatorname{Cor}_{D_n}^2(f_k, g_k)}_{\text{Corrélation entre les deux scores}}
$$

### Pourquoi cette formule est-elle remarquable ?

Pour qu'un axe de co-inertie obtienne une forte valeur propre $\lambda_k$, il faut impérativement que les **trois termes soient élevés simultanément** :

1. **$\operatorname{Cor}^2(f_k, g_k)$ doit être forte :** Les deux tableaux doivent raconter la même histoire sur cet axe.
2. **$\operatorname{Var}(f_k)$ doit être forte :** L'axe doit capturer une direction de forte variabilité dans le tableau $X$.
3. **$\operatorname{Var}(g_k)$ doit être forte :** L'axe doit capturer une direction de forte variabilité dans le tableau $Y$.

> ⚠️ **Le piège évité :** Si deux variables (une dans $X$ et une dans $Y$) sont très corrélées ($r = 0.99$) mais ne représentent que du bruit résiduel (variance de $0.001$), l'Analyse Canonique des Corrélations (ACC) va les placer en premier axe, alors que la Co-inertie va les rejeter vers les derniers axes car leur variance est négligeable !

---

## 7. Comparatif Exhaustif : CoIA vs ACC vs RDA vs Procuste vs PLS

Voici un tableau de synthèse comparatif complet :

| Caractéristique                           | Analyse de Co-inertie (CoIA)                | Analyse Canonique des Corrélations (ACC)         | Analyse Redondante (RDA / ACPVI)                                                                              | Analyse Procuste (GPA)                | Régression PLS (2B-PLS)          |
| :----------------------------------------- | :------------------------------------------ | :------------------------------------------------ | :------------------------------------------------------------------------------------------------------------ | :------------------------------------ | :-------------------------------- |
| **Auteurs / Origine**                | Chessel & Dolédec (1994)                   | Hotelling (1936)                                  | Rao (1964) / van den Wollenberg (1977)                                                                        | Gower (1975)                          | Wold (1982), Rohlf & Corti (2000) |
| **Objectif maximisé**               | $\operatorname{Cov}^2(X u, Y v)$          | $\operatorname{Cor}^2(X u, Y v)$                | $\operatorname{Var}(\widehat{Y}_{|X})$ (variance expliquée)                                                | Minimise$\|X_{\text{rot}} - Y\|^2$  | $\operatorname{Cov}(X u, Y v)$  |
| **Statut des tableaux**              | **Symétrique** (co-structure)        | **Symétrique**                             | **Asymétrique** ($X \rightarrow Y$)                                                                  | **Symétrique** (superposition) | Asymétrique ou symétrique       |
| **Contrainte de dimension**          | **Aucune** ($p, q \gg n$ autorisé) | **Très stricte :** $n > p + q$           | **Stricte :** $n > p$ | Aucune (après réduction dimensionnelle) | Aucune ($p, q \gg n$ autorisé) |                                       |                                   |
| **Inversion de matrice**             | **Aucune** (SVD spectrale pure)       | $(X^T X)^{-1}$ et $(Y^T Y)^{-1}$ nécessaires | $(X^T X)^{-1}$ nécessaire                                                                                  | Aucune                                | Aucune                            |
| **Sensibilité à la colinéarité** | **Insensible**                        | Extrêmement sensible                             | Sensible (inflation des variances)                                                                            | Faible                                | Insensible                        |
| **Types d'analyses couplées**       | ACP, AFC, ACM, Hill-Smith...                | Variables continues standard                      | Variables continues standard                                                                                  | Nuages de coordonnées quelconques    | Variables continues standard      |

---

## 8. Typologie des Couplages et Gestion des Métriques

Dans le package `ade4`, la fonction `coinertia()` peut coupler n'importe quelle paire d'analyses d'inertie élémentaires, à condition que **les poids des individus soient rigoureusement identiques**.

### 8.1. Couplage ACP - ACP (deux tableaux quantitatifs)

* **Contexte :** Deux tableaux de mesures quantitatives continues (ex: paramètres morphologiques vs paramètres physiologiques).
* **Métrique des lignes :** $D_n = \frac{1}{n} I_n$.
* **Métrique des colonnes :** $Q_X = \operatorname{diag}(1/s_X^2)$ et $Q_Y = \operatorname{diag}(1/s_Y^2)$ (ACP normées).

### 8.2. Couplage ACP - AFC (Environnement vs Abondances d'espèces) — *Le cas classique en écologie*

* **Contexte :** Tableau $X$ de physico-chimie continue ($n \times p$) et tableau $Y$ de comptages d'espèces ($n \times q$).
* **L'analyse de $Y$ (espèces) :** On effectue une **AFC (Analyse Factorielle des Correspondances)** car le tableau contient des zéros et des abondances. L'AFC pondère chaque station par sa richesse totale relative :
  $$
  p_{i\cdot} = \frac{\sum_j Y_{ij}}{\sum_{i,j} Y_{ij}}
  $$

  La métrique des lignes de l'AFC est donc $D_n = \operatorname{diag}(p_{i\cdot})$.
* **L'analyse de $X$ (milieu) :** Pour que le couplage soit mathématiquement valide, l'ACP du milieu doit **obligatoirement adopter la même métrique des lignes** :
  ```r
  # Dans R avec ade4 :
  coa_spe <- dudi.coa(df_spe, scannf = FALSE, nf = 2)
  pca_env <- dudi.pca(df_env, row.w = coa_spe$lw, scannf = FALSE, nf = 2) # << Utilisation de coa_spe$lw
  coi <- coinertia(pca_env, coa_spe, scannf = FALSE, nf = 2)
  ```

### 8.3. Couplage AFC - AFC (Comptages vs Comptages)

* Deux tableaux d'abondances ou de présences/absences sur les mêmes stations (ex: Macroinvertébrés benthiques vs Poissons).

### 8.4. Couplage ACM - ACP / Hill-Smith (Données qualitatives ou mixtes)

* Si le tableau $X$ contient des facteurs qualitatifs (substrat, type de sol, zone géographique), on utilise `dudi.acm()` ou `dudi.hillsmith()`.

---

## 9. Anatomie Complète d'un Objet `coinertia` sous R (`ade4`)

Lorsqu'on exécute `coi <- coinertia(dudiX, dudiY, scannf = FALSE, nf = 2)`, l'objet retourné contient les champs suivants :

| Champ                                                                    | Type / Dimension                                                                       | Description mathématique et signification                                     |
| :----------------------------------------------------------------------- | :------------------------------------------------------------------------------------- | :----------------------------------------------------------------------------- |
| **`$eig`** | Vecteur ($\min(p, q)$)                            | Les valeurs propres de co-inertie$\lambda_k = \operatorname{Cov}_{D}^2(f_k, g_k)$.   |                                                                                |
| **`$RV`** | Scalaire ($[0, 1]$)                                | Le coefficient RV d'Escoufier calculé entre les deux tableaux.                        |                                                                                |
| **`$rank`**                                                      | Entier                                                                                 | Le rang mathématique de l'analyse conjointe.                                  |
| **`$nf`**                                                        | Entier                                                                                 | Le nombre d'axes factoriels retenus par l'utilisateur (généralement 2 ou 3). |
| **`$tabX`** | Matrice ($n \times p$)                           | Le tableau$X$ transformé et pondéré tel qu'utilisé par l'analyse.                |                                                                                |
| **`$tabY`** | Matrice ($n \times q$)                           | Le tableau$Y$ transformé et pondéré tel qu'utilisé par l'analyse.                |                                                                                |
| **`$c1`** | Matrice ($p \times \text{nf}$)                     | Coordonnées canoniques des variables de$X$ (vecteurs propres $Q_X u_k$).          |                                                                                |
| **`$li1`** | Matrice ($q \times \text{nf}$)                    | Coordonnées canoniques des variables de$Y$ (vecteurs propres $Q_Y v_k$).          |                                                                                |
| **`$l1`** / **`$lX`** | Matrice ($n \times \text{nf}$) | **Scores des individus projetés selon le tableau $X$ :** $f_k = X Q_X u_k$. |                                                                                |
| **`$m1`** / **`$mY`** | Matrice ($n \times \text{nf}$) | **Scores des individus projetés selon le tableau $Y$ :** $g_k = Y Q_Y v_k$. |                                                                                |
| **`$aX`** | Matrice ($p \times \text{nf}$)                     | Poids canoniques des variables de$X$.                                                |                                                                                |
| **`$aY`** | Matrice ($q \times \text{nf}$)                     | Poids canoniques des variables de$Y$.                                                |                                                                                |

---

## 10. Interprétation des Sorties Graphiques

L'exécution de la commande standard `plot(coi)` génère une planche composite à 4 cadrans :

```
┌─────────────────────────────────┬─────────────────────────────────┐
│  1. Valeurs propres ($eig)      │  2. Projections variables X     │
│  Barplot décroissant des inerties│  Cercle des corrélations ($c1)  │
├─────────────────────────────────┼─────────────────────────────────┤
│  3. Graphe des flèches          │  4. Projections variables Y     │
│  Liaison stations X -> Y ($l1/$m1)│  Dispersion espèces ($li1)      │
└─────────────────────────────────┴─────────────────────────────────┘
```

### 10.1. Le graphe des flèches (Individus / Stations)

C'est le graphique diagnostique le plus important :

* Chaque station $i$ est représentée par **deux points reliés par une flèche** :
  * L'origine de la flèche correspond à la projection de la station selon le milieu ($X$).
  * La pointe de la flèche correspond à la projection de la station selon la biologie ($Y$).
* **Interprétation de la longueur de la flèche :**
  * **Flèche très courte :** Très forte adéquation entre l'environnement et la communauté biologique présente sur ce site.
  * **Flèche très longue :** Forte discordance / résidu important. Le peuplement biologique ne correspond pas à ce que l'environnement physico-chimique prédisait (ex: pollution récente non détectée par la chimie ponctuelle, effet d'historique, barrière de dispersion).

### 10.2. Le graphe des variables

* La superposition mentale des graphiques des variables de $X$ et de $Y$ permet de lire directement les **optima écologiques** :
  * Les espèces projetées dans la même direction qu'une variable environnementale (ex: nitrates, température) sont positivement associées à ce facteur.
  * Deux variables de tableaux différents pointant dans la même direction présentent une forte co-variation structurelle.

---

## 11. Généralisation à $K$ Tableaux : Co-inertie Multiple & STATIS

Lorsque l'on dispose de plus de deux tableaux mesurés sur les mêmes individus (ex: $K$ prélèvements temporels, ou Physico-chimie + Bactéries + Champignons + Faune du sol), la co-inertie se généralise sous deux formes majeures :

### 11.1. L'Inter-structure (Matrice des RV)

Soient $K$ tableaux $X_1, X_2, \dots, X_K$.
On calcule la matrice carrée $(K \times K)$ des coefficients de corrélation vectorielle :

$$
\mathbf{RV} = \begin{pmatrix} 1 & RV(X_1, X_2) & \dots & RV(X_1, X_K) \\ RV(X_2, X_1) & 1 & \dots & RV(X_2, X_K) \\ \vdots & \vdots & \ddots & \vdots \\ RV(X_K, X_1) & RV(X_K, X_2) & \dots & 1 \end{pmatrix}
$$

La première valeur propre et le premier vecteur propre de cette matrice $\mathbf{RV}$ (méthode **STATIS** de L'Hermier des Plantes 1976) définissent des poids optimaux $\alpha_k$ attribués à chaque tableau.

### 11.2. Le Compromis

On construit un opérateur de synthèse pondéré, appelé **compromis** :

$$
W_{\text{compromis}} = \sum_{k=1}^K \alpha_k W_{X_k}
$$

L'analyse spectrale de cet opérateur donne une vue consensuelle de l'ensemble des tableaux.

### 11.3. L'Analyse de Co-inertie Multiple (ACOM / MCOA)

Dans `ade4::mcoa()`, l'algorithme recherche un espace de référence synthétique commun où chaque tableau vient se projeter, permettant de comparer simultanément $K$ tableaux de nature hétérogène.

---

## 12. Script R Complet et Cas Pratique Écologique

Voici un script reproductible complet utilisant le jeu de données classique de la rivière Doubs (`data(doubs)`) disponible dans `ade4` :

```r
# ==============================================================================
# SCRIPT DE RÉFÉRENCE : ANALYSE DE CO-INERTIE SOUS R AVEC ADE4
# ==============================================================================

# 1. Chargement des bibliothèques
if (!requireNamespace("ade4", quietly = TRUE)) install.packages("ade4")
library(ade4)

# 2. Chargement des données écologiques de la rivière Doubs
# env : 30 sites x 11 variables environnementales (pente, débit, pH, oxygène...)
# fish : 30 sites x 27 espèces de poissons
data(doubs)
env  <- doubs$env
fish <- doubs$fish

cat("Dimension Milieu :", dim(env), "\n")
cat("Dimension Poissons :", dim(fish), "\n")

# 3. Étape 1 : Analyses préliminaires séparées
# AFC sur la communauté de poissons (gestion des zéros et abondances)
coa_fish <- dudi.coa(fish, scannf = FALSE, nf = 2)

# ACP sur l'environnement en utilisant IMPÉRATIVEMENT les poids de ligne de l'AFC !
pca_env <- dudi.pca(env, row.w = coa_fish$lw, scannf = FALSE, nf = 2)

# 4. Étape 2 : Analyse de Co-inertie
coi <- coinertia(pca_env, coa_fish, scannf = FALSE, nf = 2)

# 5. Résultats chiffrés fondamentaux
cat("\n--- VALEURS PROPRES DE CO-INERTIE ---\n")
print(coi$eig)
cat("Pourcentage axe 1 :", round(coi$eig[1] / sum(coi$eig) * 100, 2), "%\n")
cat("Coefficient RV observé :", round(coi$RV, 4), "\n")

# 6. Test de significativité par Monte-Carlo (randtest)
set.seed(42)
test_mc <- randtest(coi, nrepet = 999)
cat("\n--- TEST DE MONTE-CARLO SUR LE RV ---\n")
print(test_mc)

# 7. Sorties graphiques
# Planche complète
plot(coi)

# Zoom sur le graphe des flèches (concordance des sites)
s.match(coi$lX, coi$mY, clabel = 0.7, 
        sub = "Graphe des flèches : Milieu (rond) -> Poissons (flèche)")
```

---

## 13. Foire Aux Questions, Pièges et Questions d'Examen

### Q1 : Que signifie concrètement une longue flèche sur le graphe des individus ?

> **Réponse :** Une flèche longue indique une divergence (un résidu fort) entre les deux tableaux pour cette station. Sa composition biologique s'écarte significativement de ce que laissait prévoir son profil physico-chimique.

### Q2 : Peut-on réaliser une co-inertie quand $p > n$ (plus de variables que d'échantillons) ?

> **Réponse :** **Oui, absolument.** C'est l'un des avantages majeurs de la CoIA par rapport à la RDA ou l'ACC. La CoIA ne nécessite aucune inversion de matrice matricielle ($X^T X$ n'a pas besoin d'être inversé), elle repose sur une décomposition spectrale qui converge toujours, même en métagénomique avec $p = 10\,000$ gènes et $n = 15$ échantillons.

### Q3 : Pourquoi l'ACP de l'environnement doit-elle être pondérée par `coa$lw` ?

> **Réponse :** Pour que la covariance vectorielle $\operatorname{Tr}(W_X D_n W_Y D_n)$ et le produit scalaire entre individus aient un sens géométrique, l'espace des individus doit être doté d'une **unique et même métrique $D_n$**. Si l'AFC utilise des poids marginaux $p_{i\cdot}$, l'ACP doit obligatoirement utiliser ces mêmes poids de ligne.

### Q4 : Quelle est la différence entre le coefficient RV et la p-value du `randtest` ?

> **Réponse :**
>
> * Le **$RV$** est une mesure de la **taille d'effet** (l'intensité intrinsèque du lien, entre 0 et 1).
> * La **p-value** mesure la **significativité statistique** contre l'absence de lien ($H_0$).
> * *Exemple :* Avec un très petit échantillon ($n = 5$), on peut observer un $RV = 0.60$ non significatif ($p = 0.15$). À l'inverse, avec $n = 500$, un $RV = 0.12$ peut être hautement significatif ($p = 0.001$).

### Q5 : Pourquoi préférer la CoIA à l'Analyse Canonique des Corrélations (ACC) ?

> **Réponse :** L'ACC maximise la corrélation sans tenir compte de la variance expliquée, ce qui conduit à sur-interpréter des axes mathématiques artificiels sans réalité biologique. De plus, l'ACC est inapplicable dès que $p + q \ge n$ et s'effondre en cas de colinéarité. La CoIA maximise la covariance, garantissant un compromis optimal entre corrélation et inertie réelle.
