# Analyse de Co-inertie (CoIA) — Analyse Multitableaux

Ce dépôt rassemble les supports théoriques, méthodologiques et pratiques sur l'**Analyse de Co-inertie** (CoIA, package R `ade4`), ainsi qu'une présentation interactive au format Quarto.

---

## 📁 Structure du projet

```text
.
├── cours/                      # Cours et manuels théoriques complets
│   ├── COURS_COINERTIE.md      # Synthèse de cours (triplets statistiques, inertie, RV, dualité)
│   └── MANUEL_COMPLET_COINERTIE.md # Manuel détaillé de référence (mathématiques, code R, interprétations)
├── presentation/               # Diaporama Quarto (Reveal.js)
│   ├── Presentation.qmd        # Source de la présentation
│   └── Presentation.html       # Présentation compilée en HTML autonome (embed-resources)
├── docs/                       # Articles scientifiques et supports de référence (PDF)
│   ├── Chessel_Dufour_cssb8_Coinertie.pdf
│   ├── Course6_Coinertia_Analysis.pdf
│   ├── Culhane_etal_2003_Coinertia_Microarrays.pdf
│   ├── Dray_etal_2003_Coinertia.pdf
│   └── PBIL_tdr62_TD_Coinertie_R.pdf
├── scripts/                    # Scripts R pour les analyses et exemples
├── .gitignore                  # Fichiers et dossiers ignorés par Git
└── README.md                   # Présentation générale du projet
```

---

## 🎯 Objectifs & Principe de la méthode

L'analyse de co-inertie (Dolédec & Chessel 1994, Dray et al. 2003) est une méthode d'ordination symétrique permettant de coupler **deux tableaux de données mesurés sur les mêmes individus** :
- **Tableau $X$ ($n \times p$)** : Ex. variables environnementales / physico-chimiques
- **Tableau $Y$ ($n \times q$)** : Ex. abondances d'espèces ou diversité taxonomique / fonctionnelle

### Pourquoi utiliser la co-inertie ?
Contrairement à la régression ou aux analyses canoniques classiques (ACC, RDA), la co-inertie :
1. **Tolère un grand nombre de variables** ($p, q \gg n$) sans problème de multicolinéarité ou de sur-ajustement.
2. **Offre une grande flexibilité** : chaque tableau peut être analysé préalablement par la méthode la plus appropriée (ACP normée ou non, AFC, PCoA, ACM) via le formalisme des triplets statistiques $(X, Q_X, D_n)$ et $(Y, Q_Y, D_n)$.
3. **Maximise la covariance** : elle identifie des axes projetés qui maximisent simultanément la variance de chaque nuage et la corrélation entre les projections.

---

## 🚀 Utilisation & Rendu

### Prérequis R
Pour exécuter les analyses ou recompiler la présentation :
```r
install.packages(c("ade4", "adegraphics", "knitr"))
```

### Visualiser la présentation
- Ouvrez directement [presentation/Presentation.html](file:///Users/hin.benjamin/AMD_projet/presentation/Presentation.html) dans n'importe quel navigateur web.
- Pour recompiler le diaporama avec [Quarto](https://quarto.org/) :
  ```bash
  quarto render presentation/Presentation.qmd
  ```

---

## 👥 Auteurs
- Quentin
- Inès
- Benjamin
