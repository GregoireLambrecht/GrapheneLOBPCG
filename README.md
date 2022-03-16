# GrapheneLOBPCG

Ces codes sont en liens avec le projet de département de première année de l'Ecole des Ponts et Chaussées.

PuitsPotentiel est relatif à la résolution de l'équation de Schrodinger indépendante du temps associée à un puits de potentiel. Le fichier LOBPCG est
fournit par nos encadrants et codés par des chercheurs du CERMICS (laboratoire de mathématiques appliquées de l'Ecole). Ce fichier comporte un algorithme de recherche de
valeurs propres, et de vecteurs propres. Le fichier PuitsPotentiel est la résolution de l'équation en 1D, 2D et 3D, pour un potentiel V donnée par l'utilisateur. Il trace la densité de Probabilité pour le mode demeandé. Il présente aussi un algorithme de comparaison des temps d'éxecution entre les fonctions de LinearAlgebra et de LOBPCG

Les fichiers DiagBand et DiagBandCalc sont nos traveaux sur la résolution de l'équation pour un potentiel en nid d'abeille en 2D. Le fichier DiagBand consiste dans le calcul des valeurs propres de Hk pour un k donné et un poteniel V donné, en 2D.
Le fichier DiagBandCalc trace un diagramme de bande pour un poteniel donné. Dans l'exemple, le poteniel n'est pas un potentiel en nid d'abeille, mais un potentiel choisit aléatoirement. Nous n'avons pas eu le temps de tester notre programme sur des situations physiques. Donc nous n'avons aucun moyen de savoir si notre code est juste.
