# Analyse de co-inertie



## 1. Introduction & Problématique


### 1.1. 2 tableaux 

Il s'agit d'étude la structure de variabilité dans deux ou plusieurs tableaux, on veut chercher à discuter de la relations entres les tableaux au niveau des relations entre les structures. En faisant le parallèle avec qu'en analyse multi-tableaux on regarde avec importance comment se comporte les structures le comportement de chacune des analyses séparées. Par exemple, en science voire en écologie, microbiologie et en gémonique, il y'a des cas où l'on est confronté en couplage de deux points de vue différents sur les mêmes individus : 

- Tableau $X (lignes * p colonnes)$  : Méta-données
- Tableau $Y (lignes * p colonnes)$  : Table d'abaondance 

Les deux tableaux partagent obligatoirement les mêmes individus. 

### 1.2. Limites des approches hors co-inertie
Nous utilisons cette méthode car elles possèdent des aspects qui ne sont pas bien aborder dans les autres méthodes. Par exemple dans le cas de faire deux ACP sur les deux tablaux ne marche pas car rien les phénomènes capturés chez l'un et les autres axes de sont pas les mêmes et peuvent capturé des phénomènes indépendant. 

Ainsi, l'analyse de la co-inertie nous permet de trouver un ensemble d'axe qui extrait simultanément la variabilité maximale de chaque tableau tout en maximisant leur corrélation. 

## 2. Formalisation de la co-inertie

On va cherhcer à associé chaque tablaeu à un triplet statistique qui est chacun décomposée de la manière suivante avec : une matrice diagonale, et une matice Q associé à un ensemble de métrique associé au variable. 

### 2.1 Inertie et nuage de points 

L'inertie mesure la dispersion du nuage par rapport à son centre de gravité, elle mesure la taille globale du nuage mais pas sa forme ni sa direction principale d'élongation. Elle est de la forme : 

I(X) = \operatorname{Tr}(X Q_X X^T D_n) = \sum_{k=1}^r \lambda_k

### 2.2 Variance et Covariance vectorielle 

Ces métriques nous permetent de mesurer la structure interne et la forme du nuage ainsi que la liaison globale entre les deux structures de données. 


# Résumé de la méthode 

On va chercher à mettre en relation deuc tableaux de données distincts qui ont été mesuré sur les mêmes individus. On veut expliquer un tableau par un autre, on veut chercher à trouver un strcuture commune au deux tableaux. Autrement dit voir si les individus s'organisent de la même manière d'un point du tableau 1 et du tableau. 

On commence par analyser chaque tableau séparement avce la méthode la plus approprié. On va ensuite chercher l'axe d'ordinatioin pour lequel on a la le max de covariance. La co inertie peut être vue comme la corrélation de la variance du tableau 1 et celui du tableau 2. Cette méthode nous permet de d'anayser de jeu de données avec beaucoup de variable mais aussi elle est plus felxible vis-à-vis du type de tableau mais il faut faire la bonne analyse au début. 

Afin d'interpréter au mieux cette méthode il faut voir regarder le Rv qui mesure la corrélation entre les deux tableaux. On peut tester cette corrélation via le test de Monte-Carlo. Dans le graphique en flèche il faut regarder si les flèches sont courtes ou pas cela nous donne des informations sur la concordance entre les deux tableaux est fortes et si elles le sont toutes cela nous permet de dire que la co-struture du tableau est forte.  

