#include "mainClusPack.h"
#include "divide.h"
#include "2means.h"
#include "dpc.h"
#include "Maths/space.h"
#include "Maths/tools.h"
#include "normalized_2cut.h"

/**************************************/
/*                                    */
/*Clustering par divisions successives*/
/*                                    */
/**************************************/

void divide(int nbIndividus,int nbDimensions,int *nbClusters,individu_t **individus,
	    int typeDonnees,int clusteringMethod,int nbClustersSelection,
	    int donneesNormalisees,int typeDensite,int nbVoisins)
{
  /*declaration des variables*/
  int i,j,k,resultatTestGroupe1,resultatTestGroupe2,compteur,nbFoisQuonSepare;
  int nbSimulations,succesClassification;
  cluster_t *premierGroupe,*groupeCourant,*groupe1,*groupe2,*nouveauGroupe;
  double **coordonneesIndividus,*cotesHyperpave,*valeursMinimum,*valeursMaximum;
  double qualite,risque;
  /*fin declaration des variables*/

  /********************/
  /*                  */
  /*Corps du programme*/
  /*                  */
  /********************/

  /*allocation memoire*/
  premierGroupe=(cluster_t *)malloc(sizeof(cluster_t));
  premierGroupe->individus=(individu_t **)malloc(sizeof(individu_t *)*nbIndividus);
  premierGroupe->precedent=NULL;
  premierGroupe->suivant=NULL;

  groupe1=(cluster_t *)malloc(sizeof(cluster_t));
  groupe1->individus=(individu_t **)malloc(sizeof(individu_t *)*nbIndividus);

  groupe2=(cluster_t *)malloc(sizeof(cluster_t));
  groupe2->individus=(individu_t **)malloc(sizeof(individu_t *)*nbIndividus);
  coordonneesIndividus=(double **)malloc(sizeof(double *)*nbIndividus);
  for(i=0;i<nbIndividus;i++)
    {
      coordonneesIndividus[i]=(double *)malloc(sizeof(double)*nbDimensions);
    }
  cotesHyperpave=(double *)malloc(sizeof(double)*nbDimensions);
  valeursMinimum=(double *)malloc(sizeof(double)*nbDimensions);
  valeursMaximum=(double *)malloc(sizeof(double)*nbDimensions);
  /*fin allocation memoire*/

  for(i=0;i<nbIndividus;i++)
    {
      for(j=0;j<nbDimensions;j++)
	{
	  coordonneesIndividus[i][j]=individus[i]->valeursTraitees[j];
	}
    }

  if(nbClustersSelection==DPC)
    {
      qualite=0;
      risque=0.01;
    
      if(typeDonnees==ALIGNEMENT)
	{
	  risque=0.01;
	}
      else if(typeDensite==DENSITE2)
	{
	  risque=0.01;
	}  
      
      nbSimulations=20;
      /*      if(typeDonnees==ALIGNEMENT)
	{ 
	  for(i=0;i<nbDimensions;i++)
	    {
	      valeursMinimum[i]=0;
	      cotesHyperpave[i]=1.0;
	    }
	}
      else
      {*/
	  /*on cherche les valeurs minimum et maximum pour chaque axe*/
	  for(i=0;i<nbDimensions;i++)
	    {
	      valeursMinimum[i]=coordonneesIndividus[0][i];
	      valeursMaximum[i]=coordonneesIndividus[0][i];
	    }
	  
	  for(i=0;i<nbIndividus;i++)
	    {
	      for(j=0;j<nbDimensions;j++)
		{
		  if(valeursMinimum[j]>coordonneesIndividus[i][j])
		    {
		      valeursMinimum[j]=coordonneesIndividus[i][j];
		    }
		  if(valeursMaximum[j]<coordonneesIndividus[i][j])
		    {
		      valeursMaximum[j]=coordonneesIndividus[i][j];
		    }
		}
	    }
     
	  /*on calcule la longueur des cotes de l'hyperpave*/
	  for(i=0;i<nbDimensions;i++)
	    {
	      cotesHyperpave[i]=valeursMaximum[i]-valeursMinimum[i];
	    }
	  /*}*/
    }

  /*on initialise le premier groupe*/
  premierGroupe->nbIndividus=nbIndividus;
  for(i=0;i<nbIndividus;i++)
    {
      premierGroupe->individus[i]=individus[i];
    }
  premierGroupe->insecable=NON;

  *nbClusters=1;
  while(1)
    {
      groupeCourant=premierGroupe;
      while(groupeCourant!=NULL)
	{
	  if(groupeCourant->insecable==NON)
	    {
	      break;
	    }
	  else
	    {
	      groupeCourant=groupeCourant->suivant;
	    }
	}
      if(groupeCourant==NULL)
	{
	  break;
	}

      /*on regarde si dans le groupe courant il y a au moins deux individus differents*/
      /*auquel cas on peut appliquer l'algo 2-means*/
      if(clusteringMethod==KMEANS)
	{
	  /*2-means*/
	  succesClassification=
	    twoMeans(nbDimensions,groupeCourant,groupe1,groupe2,typeDonnees);	  
	}
      
      if((groupe1->nbIndividus<0)||(groupe2->nbIndividus<0))
	{
	  groupeCourant->insecable=OUI;
	}
      else
	{
	  if(succesClassification==NON)
	    {
	      resultatTestGroupe1=NON;
	      resultatTestGroupe2=NON;
	    }
	  else if(nbClustersSelection==DPC)
	    {
	      nbFoisQuonSepare=0;
	      for(i=0;i<nbSimulations;i++)
		{
		  resultatTestGroupe1=
		    testDPC(nbIndividus,nbDimensions,groupeCourant,groupe1,
			     groupe2,coordonneesIndividus,cotesHyperpave,
			     valeursMinimum,risque,typeDonnees,typeDensite,
			     donneesNormalisees,1);
		  if(resultatTestGroupe1==OUI)
		    {
		      nbFoisQuonSepare++;
		    }
		}
	      
	      if(nbFoisQuonSepare>=(int)floor(0.6*(double)nbSimulations))
		{
		  resultatTestGroupe1=OUI;
		  qualite+=(double)nbFoisQuonSepare/(double)nbSimulations;
		}
	      else
		{
		  resultatTestGroupe1=NON;
		  qualite+=(double)(nbSimulations-nbFoisQuonSepare)/
		    (double)nbSimulations;
		}
	    }
	  
	  if(resultatTestGroupe1==OUI)
	    {
	      (*nbClusters)++;
	    }
	  printf("current number of classes : %d\n",*nbClusters);
	  
	  resultatTestGroupe2=resultatTestGroupe1;
	  
	  if((resultatTestGroupe1==OUI)||(resultatTestGroupe2==OUI))
	    { 
	      /*allocation memoire*/
	      nouveauGroupe=(cluster_t *)malloc(sizeof(cluster_t));
	      nouveauGroupe->individus=(individu_t **)malloc(sizeof(individu_t *)*
							      groupe1->nbIndividus);
	      /*fin allocation memoire*/
	      
	      nouveauGroupe->precedent=premierGroupe;
	      nouveauGroupe->suivant=premierGroupe->suivant;
	      if(premierGroupe->suivant!=NULL)
		{
		  premierGroupe->suivant->precedent=nouveauGroupe;
		}
	      premierGroupe->suivant=nouveauGroupe;
	      
	      nouveauGroupe->nbIndividus=groupe1->nbIndividus;
	      for(i=0;i<nouveauGroupe->nbIndividus;i++)
		{
		  nouveauGroupe->individus[i]=groupe1->individus[i];
		}
	      
	      groupeCourant->nbIndividus=groupe2->nbIndividus;
	      for(i=0;i<groupeCourant->nbIndividus;i++)
		{
		  groupeCourant->individus[i]=groupe2->individus[i];
		}
	      
	      if(nouveauGroupe->nbIndividus<3)
		{
		  nouveauGroupe->insecable=OUI;
		}
	      else
		{
		  nouveauGroupe->insecable=NON;
		}
	      
	      if(groupeCourant->nbIndividus<3)
		{
		  groupeCourant->insecable=OUI;
		}
	      else
		{
		  groupeCourant->insecable=NON;
		}
	    }
	  else
	    {
	      groupeCourant->insecable=OUI;
	    }
	}
    }

  /*om compte le nombre de groupes*/
  groupeCourant=premierGroupe;
  *nbClusters=0;
  
  while(groupeCourant!=NULL)
    {
      for(i=0;i<groupeCourant->nbIndividus;i++)
	{
	  groupeCourant->individus[i]->cluster=*nbClusters;
	}
      groupeCourant=groupeCourant->suivant;
      (*nbClusters)++;
    }

  /*desallocation memoire*/
  free(groupe1->individus);
  free(groupe1);
  free(groupe2->individus);
  free(groupe2);

  groupeCourant=premierGroupe->suivant;
  if(groupeCourant==NULL)
    {
      free(premierGroupe->individus);
      free(premierGroupe);
    }
  else
    {
      while(1)
	{
	  free(groupeCourant->precedent->individus);
	  free(groupeCourant->precedent);

	  if(groupeCourant->suivant!=NULL)
	    {
	      groupeCourant=groupeCourant->suivant;
	    }
	  else
	    { 
	      free(groupeCourant->individus);
	      free(groupeCourant);
	      break;
	    }
	}
    }
  for(i=0;i<nbIndividus;i++)
    {
      free(coordonneesIndividus[i]);
    }
  free(coordonneesIndividus);
  free(valeursMinimum); 
  free(valeursMaximum);
  free(cotesHyperpave);
  /*fin desallocation memoire*/
}

