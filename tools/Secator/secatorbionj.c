#include "main.h"
#include "secatorbionj.h"
#include "dissimilarities_clustering.h"
#include "weighting.h"
#include "treecut.h"

/*****************************************************************/
/*                                                               */
/*Procedure de classification hierarchique base sur l'arbre bionj*/
/*avec decouverte du nombre de groupes par Secator               */
/*                                                               */
/*****************************************************************/

void secatorBionj(int nbIndividus,noeud_t *racine,int *nbClusters,int weighting) {
  /*declaration des variables*/
  int i,j,k,l,compteurIndicesDissimilarite,nbExposantsHolder,resolution,nbFeuilles;
  double *historiqueIndicesDissimilarite,seuilPertesInertie,*RMSD_groupes;
  double *exposantsHolder,**pourcentagesIdentite;
  noeud_t **adressesNoeuds,**historiqueNoeuds,**feuilles;
  /*fin declaration des variables*/
  
  resolution=0;

  /*allocation memoire*/
  adressesNoeuds=(noeud_t **)malloc(sizeof(noeud_t)*(2*nbIndividus-1));
  historiqueIndicesDissimilarite=(double *)malloc(sizeof(double)*nbIndividus);
  historiqueNoeuds=(noeud_t **)malloc(sizeof(noeud_t *)*nbIndividus);
  feuilles=(noeud_t **)malloc(sizeof(noeud_t *)*nbIndividus);
  /*fin allocation memoire*/

  /*on recupere les adresses de tous les noeuds*/
  chercheAdressesNoeuds(racine,adressesNoeuds);

  if(weighting==OUI)
    {
      /*weighting des sequences s'il y a lieu*/
      calculPoidsSequences(racine,0);
    }
  else
    {
      for(i=0;i<nbIndividus;i++)
	{
	  adressesNoeuds[i]->poids=1.0;
	}
    }

  /*maintenant on peut deroote l'arbre comme on a les poids*/
  racine->copains[1]->copains[0]=racine->copains[2];
  racine->copains[2]->copains[0]=racine->copains[1];

  /*classification hierachique des sequences en respectant la topologie de l'arbre*/
  wardSurUnArbreBionj(nbIndividus,adressesNoeuds,racine,&compteurIndicesDissimilarite,
		      historiqueIndicesDissimilarite,historiqueNoeuds);


  /*affichage de l'arbre*/
  /*affichage_arbre(&noeud_fictif,OUI);*/

  /*tri par ordre decroissant des pertes d'inertie inter-classe*/
  triRapideLocal(historiqueIndicesDissimilarite,historiqueNoeuds,0,
		 compteurIndicesDissimilarite-1);


  if(fabs(historiqueIndicesDissimilarite[0]-
	  historiqueIndicesDissimilarite[nbIndividus-2])<DBL_MIN)
    {
      *nbClusters=1;
      nbFeuilles=0;
      chercheFeuilles(racine,&nbFeuilles,feuilles);
      for(i=0;i<nbIndividus;i++)
	{
	  feuilles[i]->individu->cluster=0;
	}
    }
  else
    {
	  /*dissimilarities clustering using Holder exponent*/
	  clusteringDissimilaritesHolder(nbIndividus-1,historiqueIndicesDissimilarite,
					 nbClusters,resolution);

      /*clusters's building*/
      /*here nbClusters can change !*/
      clustersBuilding(nbIndividus,historiqueIndicesDissimilarite[*nbClusters-2],
			racine,adressesNoeuds,nbClusters);
    }  

  /*desallocation memoire*/
  for(i=0;i<2*nbIndividus-1;i++)
    {
      free(adressesNoeuds[i]);
    }
  free(adressesNoeuds);
  free(historiqueIndicesDissimilarite);
  free(historiqueNoeuds);
  free(feuilles);
  /*fin desallocation memoire*/
  
}

/********************************************************/
/*                                                      */
/*Procedure de recherche des adresses de tous les noeuds*/
/*                                                      */
/********************************************************/

void chercheAdressesNoeuds(noeud_t *noeud,noeud_t **adressesNoeuds)
{
  adressesNoeuds[noeud->numero]=noeud;

  if(noeud->copains[1]!=NULL)
    {
      chercheAdressesNoeuds(noeud->copains[1],adressesNoeuds);
      chercheAdressesNoeuds(noeud->copains[2],adressesNoeuds);
    }
}

/**************************************************************/
/*                                                            */
/*Procedure de construction de la solution dans le cas general*/
/*                                                            */
/**************************************************************/

void clustersBuilding(int nbIndividus,double seuilPertesInertie,
		      noeud_t *racine,noeud_t **adressesNoeuds,int *nbClusters)
{
  /*declaration des variables*/
  int i,j,nbAnciennesBranches,nbNouvellesBranches,nbNoeudsTrouves,changement,nbFeuilles;
  noeud_t **anciennesBranches,**nouvellesBranches,**noeudsTrouves,**feuilles;
  /*fin declaration des variables*/
  
  /*allocation memoire*/
  anciennesBranches=(noeud_t **)malloc(sizeof(noeud_t *)*nbIndividus);
  nouvellesBranches=(noeud_t **)malloc(sizeof(noeud_t *)*nbIndividus);
  noeudsTrouves=(noeud_t **)malloc(sizeof(noeud_t *)*nbIndividus);
  feuilles=(noeud_t **)malloc(sizeof(noeud_t *)*nbIndividus);
  /*fin allocation memoire*/

  nbAnciennesBranches=1;
  anciennesBranches[0]=racine;
  
  while(1)
    {
      nbNouvellesBranches=0;
      changement=NON;
      
      /*on cherche les premiers noeuds valides a partir de chaque branche*/
      for(i=0;i<nbAnciennesBranches;i++)
	{ 
	  nbNoeudsTrouves=0;

	  cherchePremiersNoeudsValides(anciennesBranches[i],seuilPertesInertie,
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

  /*on each bough leaves receive their clusters numbers*/
  for(i=0;i<nbNouvellesBranches;i++)
    {
      nbFeuilles=0;
      chercheFeuilles(nouvellesBranches[i],&nbFeuilles,feuilles);
      
      for(j=0;j<nbFeuilles;j++)
	{
	  feuilles[j]->individu->cluster=i;
	}
    }
  *nbClusters=nbNouvellesBranches;
  
  /*desallocation memoire*/
  free(anciennesBranches);
  free(nouvellesBranches);
  free(noeudsTrouves);
  free(feuilles);
  /*fin desallocation memoire*/
}

/****************************************/
/*                                      */
/*Procedure de tri par ordre decroissant*/
/*                                      */
/****************************************/

void triRapideLocal(double *valeurs,noeud_t **noeuds,int gauche,int droite)
{ 
  /*declaration des variables*/
  int elementSuivant;
  int indiceSeparateur;
  /*fin de declaration des variables*/
  
  if(gauche>=droite)
    {
      return;
    }
  
  indiceSeparateur=gauche;
  elementSuivant=gauche+1;
  
  while(elementSuivant<=droite)
    {
      if(valeurs[elementSuivant]>valeurs[indiceSeparateur])
	{
	  echangerLocal(valeurs,noeuds,indiceSeparateur+1,elementSuivant);
	  echangerLocal(valeurs,noeuds,indiceSeparateur,indiceSeparateur+1);
	  indiceSeparateur++;
	}
      elementSuivant++;
    }
  
  triRapideLocal(valeurs,noeuds,gauche,indiceSeparateur-1);
  triRapideLocal(valeurs,noeuds,indiceSeparateur+1,droite);
}

/*******************************************/
/*                                         */
/*sous-procedure de triRapide             */ 
/*                                         */
/*******************************************/

void echangerLocal(double *valeurs,noeud_t **noeuds,int element1,int element2)
{
  /*declaration des variables*/
  double tempVal;
  noeud_t *tempNoeud;
  /*fin de declaration des variables*/

  tempVal=valeurs[element1];
  valeurs[element1]=valeurs[element2];
  valeurs[element2]=tempVal;

  tempNoeud=noeuds[element1];
  noeuds[element1]=noeuds[element2];
  noeuds[element2]=tempNoeud;
}



