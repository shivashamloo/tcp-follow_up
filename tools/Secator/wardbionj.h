/******************************************/
/*                                        */
/*Procedure de classification hierarchique*/
/*                                        */
/******************************************/

void wardSurUnArbreBionj(int nbSequences,noeud_t **adressesNoeuds,
			 noeud_t *noeudFictif,int *compteurIndicesDissimilarite,
			 double *historiqueIndicesDissimilarite,
			 noeud_t **historiqueNoeuds);

/*********************************************************/
/*                                                       */
/*Procedure d'initialisation des indices de dissimilarite*/
/*                                                       */
/*********************************************************/

void initialisationIndicesDissimilarite(int nbIndividus,double **indicesDissimilariteAnciens,
					noeud_t **adressesNoeuds);

/***************************************************/
/*                                                 */
/*Procedure de mise a jour des indices de Ward     */
/*lors de la classification hierarchique ascendante*/
/*                                                 */
/***************************************************/

void miseAJourIndicesWard(int groupe1,int groupe2,noeud_t **adressesNoeuds,
			  int nbGroupesCourants,int *numerosGroupesCourants,
			  double **indicesWardAnciens,double **indicesWardNouveaux);
