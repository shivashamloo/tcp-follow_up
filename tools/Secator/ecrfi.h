/**********************************************************************************/
/*                                                                                */
/*Procedure d'ecriture du fichier.dst des distances en vue du calcul de phylogenie*/
/*                                                                                */
/**********************************************************************************/

void ecritureFichierDistancesSequences(char *nomFichier,int nbIndividus,
				       individu_t **individus);

/******************************************************************/
/*                                                                */
/*Procedure d'ecriture du fichier resultat presentant les clusters*/
/*                                                                */
/******************************************************************/

void ecritureFichierClusters(int nbIndividus,int nbDimensions,individu_t **individus,
			     int nbClusters,char *nomFichier);

/**********************************************************/
/*                                                        */
/*Procedure d'ecriture du fichier de l'alignement multiple*/
/*avec les sequences reordonnees suivant leurs groupes    */
/*                                                        */
/**********************************************************/

void ecritureFichierClustersAlignement(char *nomFichier,int longueurAlignement,
				       int nbIndividus,individu_t **individus,
				       char **sequences,int nbClusters);

/***************************************************/
/*                                                 */
/*Procedure d'ecriture des etiquettes d'un ensemble*/
/*                                                 */
/***************************************************/

void ecritureFichierEtiquettes(int tailleTestSet,individu_t **testSet,char *nomFichier);
