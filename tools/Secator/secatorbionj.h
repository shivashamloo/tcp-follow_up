/****************************************/
/*                                      */
/*Procedure de tri par ordre decroissant*/
/*                                      */
/****************************************/

void triRapideLocal(double *valeurs,noeud_t **noeuds,int gauche,int droite);

/*******************************************/
/*                                         */
/*sous-procedure de triRapide             */ 
/*                                         */
/*******************************************/

void echangerLocal(double *valeurs,noeud_t **noeuds,int element1,int element2);

/*****************************************************************/
/*                                                               */
/*Procedure de classification hierarchique base sur l'arbre bionj*/
/*avec decouverte du nombre de groupes par Secator               */
/*                                                               */
/*****************************************************************/

void secatorBionj(int nbIndividus,noeud_t *racine,int *nbClusters,int weighting);

/**************************************************/
/*                                                */
/*Procedure de recherche des adresses des feuilles*/
/*                                                */
/**************************************************/

void chercheAdressesNoeuds(noeud_t *noeud,noeud_t **adressesNoeuds);

/**************************************************************/
/*                                                            */
/*Procedure de construction de la solution dans le cas general*/
/*                                                            */
/**************************************************************/

void clustersBuilding(int nbIndividus,double seuilPertesInertie,
		      noeud_t *racine,noeud_t **adressesNoeuds,int *nbClusters);
