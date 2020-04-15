/*******************************************************************************/
/*                                                                             */
/*Procedure de decouverte du nombre de clusters par Secator a partir d'un arbre*/
/*                                                                             */
/*******************************************************************************/

void secatorTree(int nbIndividus,noeud_t *racine,double *seuilDissimilarite,int *nbClusters);

/********************************************************/
/*                                                      */
/*Procedure de decoupage de l'arbre en coupant quand la */
/*dissimilarite est superieure au seuil de dissimilarite*/
/*                                                      */
/********************************************************/

void decoupageArbreSeuilDissimilarite(int nbMotifs,noeud_t *noeud,double seuilDissimilarite,
				      int *nbClusters);

/******************************************************************/
/*                                                                */
/*Procedure de recherche des premiers noeuds valides d'une branche*/
/*i.e les premiers a avoir une perte d'inertie superieur au seuil */
/*                                                                */
/******************************************************************/

void cherchePremiersNoeudsValides(noeud_t *noeud,double seuilDissimilarite,
				  int *nbNoeudsTrouves,noeud_t **noeudsTrouves);

/*************************************************************************/
/*                                                                       */
/*Procedure de recherche de toutes les feuilles a partir d'un noeud donne*/
/*                                                                       */
/*************************************************************************/

void chercheFeuilles(noeud_t *noeud,int *nbFeuilles,noeud_t **feuilles);

/**********************************************************************************/
/*                                                                                */
/*Procedure de decouverte du seuil de dissimilarite a partir du nombre de clusters*/
/*                                                                                */
/**********************************************************************************/

void nbClustersToDissimilarityThreshold(int nbIndividus,int nbClusters,noeud_t *racine,
					double *seuilDissimilarite);

/*******************************************************************/
/*                                                                 */
/*Procedure recursive de recuperation des dissimilarites d'un arbre*/
/*                                                                 */
/*******************************************************************/

void recupereDissimilarites(noeud_t *noeud,int *compteur,double *dissimilarites);
