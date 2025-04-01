Scripts 

*Prep_forcage_met_SiSPAT_Brune
à partir d'un tableau de toutes les obs préparation du fichier de forçage climatique pour les dates souhaitées 

*Prep_CI_SiSPAT_v2_Brune
à partir de mesures des sondes d'humidité et température dans le sol prépération des fichiers de profil deconditions initales (plusieurs options de conditions à la limite) 
ATTENTION 
celui créé pas parfait -> il faut changer la première ligne et rajouter une ligne (mêmes valeurs) pour les temp et dupliquer la dernière humidité 

*Compar_obs_Brune_p1
ouverture des sorties de modélisation au bon format (3 ouverts en simulatnés) 
+tableaux des mesures faites à la station
calcul de scores, graphiques de visualisation, comparaison entre les 3 simus et les obs 
bilan radiatif, bilan d'énergie, evpotraspiration , bilan d'eau, températures et humidités dans le sol 

*Graphes_rapport
la même chose que celui d'avant mais graphes spécifiques créés pour le rapport de stage 

*retention 
tracer les courbes théoriques de Van genuchten et Brooks and corey  h(tetha) et K(tetha) 

*sols simu: visualisation des fichiers sols (temp, humidité, a tous les nœuds et tous les 30 min) 
*Gener_param_files: fichier de génération des fichiers de paramètres en série

* scores_ens     
calcul des scores de l'ensemble généré rangé sous forme de matrice 

* analyse scores ensembles
à partir du tableau des scores, analyse des simulations, test des critères pour discriminer un optimum de meilleures simulations 

