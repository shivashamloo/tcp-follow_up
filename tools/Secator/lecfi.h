/***********************************************************/
/*                                                         */
/*Procedure de lecture du fichier d'un arbre phylogenetique*/
/*(au format bionj)                                        */
/*                                                         */
/***********************************************************/

void lectureArbrePhylogenetique(char *nomFichier,char **arbre);

/**********************************************************/
/*                                                        */
/*Procedure de lecture du fichier d'un alignement multiple*/
/*                                                        */
/**********************************************************/

void lectureFichierAlignement(char *nomAlignement,int *longueurAlignement,int *nbIndividus,
				individu_t ***individus,char ***sequences);

/***************************************/
/*                                     */
/*Procedure de lecture d'un fichier.tfa*/
/*                                     */
/***************************************/

void lectureFichierTfa(char *nomAlignement,int *longueurAlignement,
			 int *nbIndividus,char ***sequences,individu_t ***individus);

/***************************************/
/*                                     */
/*Procedure de lecture d'un fichier.msf*/
/*                                     */
/***************************************/

void lectureFichierMsf(char *nomAlignement,int *longueurAlignement,
			 int *nbIndividus,char ***sequences,individu_t ***individus);

/********************************************/
/*                                          */
/*Procedure de lecture du fichier contenant */
/*les coordonnees des individus a classifier*/
/*                                          */
/********************************************/

void lectureFichierCoordonnees(int *nbIndividus,int *nbDimensions,individu_t ***individus,
				 char *nomFichierEntree,int typeDonnees);

/************************************************/
/*                                              */
/*Procedure de lecture du fichier contenant     */
/*les distances entre les individus a classifier*/
/*Les individus sont stockes dans un tableau    */
/*de pointeurs sur structure                    */
/*                                              */
/************************************************/

void lectureFichierDistances(int *nbIndividus,int nbDimensions,individu_t ***individus,
			       char *nomFichierEntree,double ***distances);


/*lecture du fichier des etiquettes designant les clusters*/
void lectureEtiquettes(int nbIndividus,individu_t **individus,char *nomFichier);
