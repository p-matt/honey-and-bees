# Algo génétique:
> 100 abeilles doivent parcourir 50 positions en un minimum de temps sur un plan donné  
> Chaque abeille a sa propre course, une abeille correspond à une course  
> 3 type d'abeilles: parents, enfants, mutants  
> Chaque génération élimine les abeilles les moins efficaces dans le but d'améliorer la population  
> Point de départ: (500, 500)

Le programme fonctionne avec un fichier d'entrée qui inclus toutes les positions.  
En premier lieu l'utilisateur doit choisir de créer un nouveau essain (random) ou bien de repartir depuis une population déjà existante.  
Certains des paramètres sont ajustés au fil du programme (mutationEffect, mutationRate) et d'autres sont fixés avant le lancement (birthRate, nbParent ...)

Concept de l'algorithme génétique:

    La sélection : Choix des individus les mieux adaptés.
    Le croisement : Mélange par la reproduction des particularités des individus choisis.
    La mutation : Altération aléatoire des particularités d'un individu.

<img src="https://user-images.githubusercontent.com/73102263/100369230-f0082480-3004-11eb-9624-fd3eb0685e30.gif" align="center" width="500" height="500" />
