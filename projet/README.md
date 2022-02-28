Le code C++ a été développé sous macOS version 10.15.2 avec les outils
de développement de Xcode version 11.3.1, et l'environnement de
développement CLion version 2019.3.

La plupart du code C++ a été fournit par Professeur Hecht et fait partie
du cours 5MM30. Les améliorations principales sont liées à la résolution
des éléments finis en 3D avec des coefficients complexes (non-constants).

=========================================================================
Prérequis
=========================================================================

Pour compiler et exécuter le code il faut les outils informatiques
suivants :

    Compilateur (g++)
    cmake
    UMFPACK (via SuiteSparse)
    FreeFEM (pour visualiser les résultats)

Pour visualiser les résultats des calculs, le logiciel utilise FreeFEM.
Si FreeFEM n'est pas disponible, le logiciel échouera à la fin de son
exécution. Mais avant de se terminer, le logiciel écrit plusieurs fichiers
contenant des données qui peuvent être utilisés pour générer des
visualisations.

=========================================================================
Compilation
=========================================================================

Pour configurer et compiler le code, saisissez les commandes suivantes
depuis le répertoire racine du projet :

    cmake .
    cmake --build .

Si tout se passe bien le logiciel SolveOven devrait être présent
dans le répertoire racine du projet.

=========================================================================
Dépannage
=========================================================================

Si la compilation échoue, il est probable que le compilateur n'arrive
pas à trouver des fichiers en-tête liés á UMFPACK (umfpack.h).

Sous macOS avec SuiteSparse fournit par Homebrew, ces fichiers se trouvent
normalement sous /usr/local/include. Si ces fichiers se trouvent
ailleurs sur votre système il faut mettez à jour la ligne 8 du fichier
CMakeLists.txt :

    ...
    include_directories(/usr/local/include)
    ...

Après avoir changé CMakeLists.txt il faut lancer de nouveau les
2 commandes :

    cmake .
    cmake --build .

=========================================================================
Exécution
=========================================================================

Afin de lancer le logiciel, saisissez la commande suivante depuis le
répertoire racine du projet :

    ./SolveOven


SolveOven (src/apps/SolveOven.cpp) effectue les étapes suivantes :

    Résolution du problème de Helmoltz
    Résolution du problème de Poisson
    Génération des fichiers contenant les données des résultats
    Lancement de FreeFEM (le script freefem/view_oven_solution.edp)

Le script FreeFEM (freefem/view_oven_solution.edp) affiche 3 plots au total.
Pour procéder au prochain plot appuyez sur la touche entrée.

=========================================================================
Scripts FreeFEM
=========================================================================

Dans le répertoire freefem se trouve 3 scripts FreeFEM :

    view_oven_solution.edp : afficher des plots à partir
    des données générées par SolveOven
    generate_oven_mesh_3d.edp : générer un maillage 3D qui
    consiste en le four (région rectangulaire) et, facultativement, l'aliment
    à cuire (un sphère)
    solve_oven_freefem.edp : résoudre le problème du four avec
    FreeFEM (pour vérifier la solution C++)

De plus, le script view_oven_solution.edp génère des fichiers de données
du format VTK pour l'utilisation avec ParaView. La plupart des figures dans ce
rapport ont été créées avec ParaView à partir de ces fichiers VTK.

=========================================================================
Fichiers de données
=========================================================================

Durant l'exécution, le logiciel SolveOven et le script
view_oven_solution.edp génèrent des fichiers de données qui sont
sauvegardés dans le répertoire visualisation.

    freefem_helmholtz.txt contient les valeurs (complexes)
    de la solution de l'équation de Helmholtz au format convenable pour
    FreeFEM

    freefem_heat.txt contient les valeurs de la solution
    de l'équation de Poisson au format convenable pour FreeFEM

    helmholtz_real.vtu contient les valeurs réelles de la solution
    de l'équation de Helmholtz au format VTK (pour ParaView)

    helmholtz_imag.vtu contient les valeurs imaginaires de la
    solution de l'équation de Helmholtz au format VTK (pour ParaView)

    heat.vtu contient les valeurs de la solution de
    l'équation de Poisson au format VTK (pour ParaView)

    helmholtz.csv contient les valeurs de la solution de
    l'équation de Helmholtz au format CSV

    heat.csv contient les valeurs de la solution de
    l'équation de Poisson au format CSV

=========================================================================
Paramètres variables et pointeurs de fonction
=========================================================================

Afin de faciliter le calcul des paramètres variables (par exemple,
epsilon(x) et K(x)) le code C++ emploie des pointeurs
de fonction.

Pendant l'assemblage des matrices on calcule la quadrature des
équations sur chaque simplexe (élément). Cette quadrature requiert
le calcul des opérations différentielles et aussi le calcul des paramètres
variables qui apparaissent dans les équations. Chaque opération est
contenue dans la structure des données Operation. Cette
structure se trouve dans le fichier include/solver/Solver3D.hpp :

    template<class TypeScalar>
    struct Operation {
        function<TypeScalar(Simplex3&)> computeWeight;
        int operatorU, operatorV;
    };

La fonction computeWeight est utilisée pendant le processus
de l'assemblage pour calculer les paramètres variables. Dans le fichier
src/apps/SolveOven.cpp, par exemple, on attribue à
computeWeight la fonction suivante pour calculer la permittivité
relative :

    complex<double> computePermittivity(Simplex3& simplex) {
        if (inObject(simplex)) {
            return inversePermittivityObject;
        }
        return inversePermittivityAir;
    }

En redéfinissant la fonction inObject() on peut alors utiliser
cette méthode pour mettre toute sorte
d'objet dans le four sans besoin de créer des maillages supplémentaires.
On peut se servir uniquement du maillage du four vide pour effectuer
toutes les expériences numériques. Ces objets auront bien sur des
frontières plutôt irréguliers, mais pour vite tester la cuisson
de divers objets, c'est une méthode super efficace.

En outre, on peut faire varier partout dans l'enceinte (et l'aliment)
les valeurs de la permittivité et de la conductivité thermique.
