#define TROUVE 1
#define PAS_TROUVE 0

/******************************************************************************************/
/*                                                                                        */
/*Procedure construisant l'arbre bionj a partir d'un alignement de sequence dans sequences*/
/*                                                                                        */
/******************************************************************************************/

void creationArbreBionj(int nbIndividus,individu_t **individus,
			  char *fichierEntree,noeud_t **racine);


/*******************************************************/
/*                                                     */
/*Procedure de calcul des distances a partir de l'arbre*/
/*                                                     */
/*******************************************************/

void calculDistances(int nbIndividus,individu_t **individus,noeud_t *racine);

/**************************************/
/*                                    */
/*Procedure de marquage a 0 des noeuds*/
/*                                    */
/**************************************/

void marquageNoeuds(noeud_t *noeud);

/**************************************************/
/*                                                */
/*Procedure d'association des individus aux noeuds*/
/*                                                */
/**************************************************/

void associationIndividusNoeuds(int nbIndividus,individu_t **individus,noeud_t *noeud);

/**************************************************/
/*                                                */
/*Procedure d'association des individus aux noeuds*/
/*                                                */
/**************************************************/

void associationIndividusNoeuds2(int *indiceIndividu,individu_t **individus,noeud_t *noeud);

/*************************************************************/
/*                                                           */
/*Procedure d'enracinement de l'arbre sur la sequence fictive*/
/*                                                           */
/*************************************************************/

void enracinementArbre(int nbIndividus,noeud_t **racine);

/**********************************************************************/
/*                                                                    */
/*Procedure de renommage des noeuds internes car apres restructuration*/
/*le nombre de noeuds n'est plus le meme                              */
/*                                                                    */
/**********************************************************************/

void renommageNoeudsInternes(int nbIndividus,noeud_t *noeud,int *nbNoeudsInternes,
			       int *nbFeuilles);

/****************************************************************************/
/*                                                                          */
/*Procedure de restructuration de l'arbre et de calcul de la nouvelle racine*/
/*                                                                          */
/****************************************************************************/

void restructurationArbre(noeud_t *noeud,noeud_t *nouveauPere,double nouvelleDistance);

/************************************************************/
/*                                                          */
/*Procedure de recherche de la sequence fictive dans l'arbre*/
/*                                                          */
/************************************************************/

int chercheNoeudSequenceFictive(int nbIndividus,noeud_t *noeud,
				   noeud_t **noeudSequenceFictive);

/*****************************************************/
/*                                                   */
/*Procedure de construction de l'arbre phylogenetique*/
/*                                                   */
/*****************************************************/

void constructionArbre(int *nbIndividus,char *arbre,noeud_t **racine);

/*******************************************************************************/
/*                                                                             */
/*Procedure de copie d'une portion de chaine de caracteres delimitee de maniere*/
/*stricte par deux pointeurs sur char dans une autre chaine de caracteres      */
/*                                                                             */
/*******************************************************************************/

void copiePartieString(char *debut,char *fin,char *destination);

/******************************************************************************/
/*                                                                            */
/*Procedure de remplacement d'une partie d'une chaine de caracteres, delimitee*/
/*de maniere large par debut et fin, par une autre chaine de caracteres       */
/*                                                                            */
/******************************************************************************/

void remplacePartieString(char **chaine,char *debut,char *fin,char *chaineDeRemplacement);

/*******************************************************************/
/*                                                                 */
/*Procedure d'affichage de l'arbre (essentiellement pour debuggage)*/
/*                                                                 */
/*******************************************************************/

void affichageArbre(noeud_t *noeud);


