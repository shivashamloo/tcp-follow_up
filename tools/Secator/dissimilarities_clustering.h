#define EXPOSANT_HOLDER_SEUIL 0.5

/*************************************************************/
/*                                                           */
/*classification des dissimilarites pour trouver le nombre de*/
/*dissimilarites elevees et donc le nombre de groupes        */
/*(avec l'exposant de Holder)                                */
/*                                                           */
/*************************************************************/

void clusteringDissimilaritesHolder(int nbDissimilarites,double *historiqueDissimilarites,
				      int *nbDissimilaritesElevees,int resolution);

/******************************************************************************/
/*                                                                            */
/*Procedure de calcul des exposants d'Holder sur les valeurs de dissimilarites*/
/*                                                                            */
/******************************************************************************/

void calculExposantsHolder(int nbDissimilarites,double *historiqueDissimilarites,
			     int *nbExposantsHolder,double *exposantsHolder);

/**********************************************************************/
/*                                                                    */
/*Procedure de separation en deux groupes des valeurs de dissimilarite*/
/*par calcul du minimum d'inertie intraclasse                         */
/*                                                                    */
/**********************************************************************/

void separationDissimilaritesBete(int nbDissimilarites,double *historiqueDissimilarites,
				    int *nbDissimilaritesPremierGroupe);

/************************************************************************/
/*                                                                      */
/*Procedure de separation en deux groupes des valeurs de dissimilarite  */
/*par calcul du minimum d'inertie intraclasse mais uniquement aux points*/
/*d'exposants d'Holder faibles                                          */
/*                                                                      */
/************************************************************************/

void separationDissimilaritesEvolue(int nbExposantsHolder,int nbDissimilarites,
				      double *historiqueDissimilarites,
				      int *nbDissimilaritesPremierGroupe,
				      double *exposantsHolder);
