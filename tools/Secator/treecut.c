#include "main.h"
#include "treecut.h"
#include "dissimilarities_clustering.h"
#include "tools.h"

/*******************************************************************************/
/*                                                                             */
/*Procedure de decouverte du nombre de clusters par Secator a partir d'un arbre*/
/*                                                                             */
/*******************************************************************************/

void secatorTree(int nbIndividus,noeud_t *racine,double *seuilDissimilarite,int *nbClusters)
{
  /*declaration de variables*/
  int compteur;  
  double *dissimilarites;
  /*fin declaration de variables*/

  /*allocation memoire*/
  dissimilarites=(double *)malloc(sizeof(double)*(nbIndividus-1));
  /*fin allocation memoire*/

  compteur=0;
  recupereDissimilarites(racine,&compteur,dissimilarites);

  /*tri des dissimilarites par ordre decroissant*/
  triRapide(dissimilarites,0,compteur-1); 

  /*decouverte du nombre de clusters*/
  clusteringDissimilaritesHolder(nbIndividus-1,dissimilarites,nbClusters,0);

  *seuilDissimilarite=dissimilarites[*nbClusters-2];


  /*desallocation memoire*/
  free(dissimilarites);
  /*fin desallocation memoire*/
}

/**********************************************************************************/
/*                                                                                */
/*Procedure de decouverte du seuil de dissimilarite a partir du nombre de clusters*/
/*                                                                                */
/**********************************************************************************/

void nbClustersToDissimilarityThreshold(int nbIndividus,int nbClusters,noeud_t *racine,
					double *seuilDissimilarite)
{
  /*declaration de variables*/
  int i,compteur;
  double *dissimilarites;
  /*fin declaration de variables*/

  /*allocation memoire*/
  dissimilarites=(double *)malloc(sizeof(double)*(nbIndividus-1));
  /*fin allocation memoire*/

  compteur=0;
  recupereDissimilarites(racine,&compteur,dissimilarites);

  /*tri des dissimilarites par ordre decroissant*/
  triRapide(dissimilarites,0,compteur-1);

  *seuilDissimilarite=dissimilarites[nbClusters-2];

  /*desallocation memoire*/
  free(dissimilarites);
  /*fin desallocation memoire*/
}

/*******************************************************************/
/*                                                                 */
/*Procedure recursive de recuperation des dissimilarites d'un arbre*/
/*                                                                 */
/*******************************************************************/

void recupereDissimilarites(noeud_t *noeud,int *compteur,double *dissimilarites)
{
  if(noeud->copain1!=NULL)
    {
      dissimilarites[*compteur]=noeud->dissimilarite;
      (*compteur)++;
      recupereDissimilarites(noeud->copain1,compteur,dissimilarites);
      recupereDissimilarites(noeud->copain2,compteur,dissimilarites);
    }
}

/********************************************************/
/*                                                      */
/*Procedure de decoupage de l'arbre en coupant quand la */
/*dissimilarite est superieure au seuil de dissimilarite*/
/*                                                      */
/********************************************************/

void decoupageArbreSeuilDissimilarite(int nbMotifs,noeud_t *noeud,
				      double seuilDissimilarite,int *nbClusters) 
{
  /*declaration des variables*/
  int i,j,nbAnciennesBranches,nbNouvellesBranches,nbNoeudsTrouves,changement;
  int nbFeuilles;
  noeud_t **anciennesBranches,**nouvellesBranches,**noeudsTrouves,**feuilles;
  /*fin declaration des variables*/
  
  /*allocation memoire*/
  anciennesBranches=(noeud_t **)malloc(sizeof(noeud_t *)*nbMotifs);
  nouvellesBranches=(noeud_t **)malloc(sizeof(noeud_t *)*nbMotifs);
  noeudsTrouves=(noeud_t **)malloc(sizeof(noeud_t *)*nbMotifs);
  feuilles=(noeud_t **)malloc(sizeof(noeud_t *)*nbMotifs);
  /*fin allocation memoire*/

  nbAnciennesBranches=1;
  anciennesBranches[0]=noeud;
  
  while(1)
    {
      nbNouvellesBranches=0;
      changement=NON;
      
      /*on cherche les premiers noeuds valides a partir de chaque branche*/
      for(i=0;i<nbAnciennesBranches;i++)
	{ 
	  nbNoeudsTrouves=0;

	  cherchePremiersNoeudsValides(anciennesBranches[i],seuilDissimilarite,
					  &nbNoeudsTrouves,noeudsTrouves);

	  if(nbNoeudsTrouves==0)
	    {
	      nouvellesBranches[nbNouvellesBranches]=anciennesBranches[i];
	      nbNouvellesBranches++;
	    }
	  else
	    {  
	      changement=OUI;
	      for(j=0;j<nbNoeudsTrouves;j++)
		{
		  nouvellesBranches[nbNouvellesBranches]=noeudsTrouves[j]->copain1;
		  nbNouvellesBranches++;
		  
		  nouvellesBranches[nbNouvellesBranches]=noeudsTrouves[j]->copain2;
		  nbNouvellesBranches++;
		}
	    }
	}
      
      if(changement==NON)
	{
	  break;
	}
      else
	{
	  nbAnciennesBranches=nbNouvellesBranches;
	  for(i=0;i<nbAnciennesBranches;i++)
	    {
	      anciennesBranches[i]=nouvellesBranches[i];
	    }
	}      
    }

  *nbClusters=0;
  /*on constitue les groupes a partir des feuilles de chaque branche*/
  for(i=0;i<nbNouvellesBranches;i++)
    {
      nbFeuilles=0;
      chercheFeuilles(nouvellesBranches[i],&nbFeuilles,feuilles);
      
      for(j=0;j<nbFeuilles;j++)
	{
	  feuilles[j]->individu->cluster=*nbClusters;
	}
      (*nbClusters)++;
    }
  
  /*desallocation memoire*/
  free(anciennesBranches);
  free(nouvellesBranches);
  free(noeudsTrouves);
  free(feuilles);
  /*fin desallocation memoire*/

}

/******************************************************************/
/*                                                                */
/*Procedure de recherche des premiers noeuds valides d'une branche*/
/*i.e les premiers a avoir une perte d'inertie superieur au seuil */
/*                                                                */
/******************************************************************/

void cherchePremiersNoeudsValides(noeud_t *noeud,double seuilDissimilarite,
				  int *nbNoeudsTrouves,noeud_t **noeudsTrouves)
{
  if(noeud->dissimilarite>seuilDissimilarite-DBL_MIN)
    {
      noeudsTrouves[*nbNoeudsTrouves]=noeud;
      (*nbNoeudsTrouves)++;
    }
  else
    {
      if(noeud->copain1!=NULL)
	{
	  cherchePremiersNoeudsValides(noeud->copain1,seuilDissimilarite,
				       nbNoeudsTrouves,noeudsTrouves);
	  cherchePremiersNoeudsValides(noeud->copain2,seuilDissimilarite,
				       nbNoeudsTrouves,noeudsTrouves);
	}
    }
}

/*************************************************************************/
/*                                                                       */
/*Procedure de recherche de toutes les feuilles a partir d'un noeud donne*/
/*                                                                       */
/*************************************************************************/

void chercheFeuilles(noeud_t *noeud,int *nbFeuilles,noeud_t **feuilles)
{
  if(noeud->copain1==NULL)
    {
      feuilles[*nbFeuilles]=noeud;
      (*nbFeuilles)++;
    }
  else
    {
      chercheFeuilles(noeud->copain1,nbFeuilles,feuilles);
      chercheFeuilles(noeud->copain2,nbFeuilles,feuilles);
    }
}


