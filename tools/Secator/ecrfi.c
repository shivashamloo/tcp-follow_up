#include "main.h"

/**********************************************************************************/
/*                                                                                */
/*Procedure d'ecriture du fichier.dst des distances en vue du calcul de phylogenie*/
/*                                                                                */
/**********************************************************************************/

void ecritureFichierDistancesSequences(char *nomFichier,int nbIndividus,
				       individu_t **individus)
{
  /*declaration des variables*/
  FILE *fichier;
  int i,j;
  /*fin declaration des variables*/

  fichier=fopen(nomFichier,"w+");

  fprintf(fichier,"%d\n",nbIndividus+1);
  for(i=0;i<nbIndividus;i++)
    {
      fprintf(fichier,"%d ",individus[i]->id);
      for(j=0;j<nbIndividus;j++)
	{
	  fprintf(fichier,"%f ",(float)individus[i]->valeurs[individus[j]->id]);
	  if((j+1)%8==0)
	    {
	      fprintf(fichier,"\n");
	    }
	}
      fprintf(fichier,"1.0 \n");
    }
  fprintf(fichier,"fictive_sequence ");
  for(j=0;j<nbIndividus;j++)
    {
      fprintf(fichier,"1.0 ");
      if((j+1)%8==0)
	{
	  fprintf(fichier,"\n");
	}
    }
  fprintf(fichier,"0.0 \n");
  fclose(fichier);
}

/******************************************************************/
/*                                                                */
/*Procedure d'ecriture du fichier resultat presentant les clusters*/
/*                                                                */
/******************************************************************/

void ecritureFichierClusters(int nbIndividus,int nbDimensions,individu_t **individus,
			       int nbClusters,char *nomFichier)
{
  /*declaration des variables*/
  int i,j,k,compteur,*taillesClusters,nbIndividusTotal,nbClustersRectifie;
  individu_t ***clusters;
  FILE *file;
  /*fin declaration des variables*/

  /*allocation memoire*/
  taillesClusters=(int *)malloc(sizeof(int)*(nbClusters+1));
  clusters=(individu_t ***)malloc(sizeof(individu_t *)*(nbClusters+1));
  /*fin allocation memoire*/

  nbIndividusTotal=nbIndividus; 
  
  for(i=0;i<nbIndividus;i++) {
      nbIndividusTotal+=individus[i]->nbIndividusSimilaires;
    }

 // printf("nbClusters dans ecrfi : %d\n",nbClusters);
  for(i=0;i<nbClusters+1;i++)
    {
      clusters[i]=(individu_t **)malloc(sizeof(individu_t *)*nbIndividusTotal);
    }
  for(i=0;i<nbClusters+1;i++)
    {
      taillesClusters[i]=0;
    }

  for(i=0;i<nbIndividus;i++)
    {
      if(individus[i]->cluster>-1)
	{
	  clusters[individus[i]->cluster][taillesClusters[individus[i]->cluster]]=individus[i];
	  (taillesClusters[individus[i]->cluster])++;
	
	}
      else
	{
	  clusters[nbClusters][taillesClusters[nbClusters]]=individus[i];
	  (taillesClusters[nbClusters])++;
	
	  
	}
    }

  /*calcul du vrai nombre de clusters, en supprimant les clusters vides*/
  nbClustersRectifie=0;
  for(i=0;i<nbClusters+1;i++) {
    if(taillesClusters[i]>0) {
      nbClustersRectifie++;
    }
  }

  file=fopen(nomFichier,"w");
  fprintf(file,"Number of clusters : %d\n",nbClustersRectifie);
  
  compteur=0;
  for(i=0;i<nbClusters;i++)
    {
      if(taillesClusters[i]>0) {
	fprintf(file,"\nCluster %d ; size=%d\n",compteur,taillesClusters[i]);
	for(j=0;j<taillesClusters[i];j++)
	  {
	    fprintf(file,"%s",clusters[i][j]->nom);
	    if(clusters[i][j]->description!=NULL)
	      {
		fprintf(file,"%s",clusters[i][j]->description);
	      }
	    else
	      {
		fprintf(file,"\n");
	      }
	  }
	compteur++;
      }
    }
  if(taillesClusters[nbClusters]>0)
    {
      fprintf(file,"\nunclustered points %d ; size=%d\n",i,taillesClusters[nbClusters]);
      for(j=0;j<taillesClusters[nbClusters];j++)
	{
	  fprintf(file,"%s",clusters[nbClusters][j]->nom);
	  if(clusters[i][j]->description!=NULL)
	    {
	      fprintf(file,"%s",clusters[i][j]->description);
	    }
	  else
	    {
	      fprintf(file,"\n");
	    }
	}
    }
  fclose(file);
 
  /*desallocation memoire*/
  free(taillesClusters);
  for(i=0;i<nbClusters;i++)
    {
      free(clusters[i]);
    }
  free(clusters);
  /*fin desallocation memoire*/
}

/**********************************************************/
/*                                                        */
/*Procedure d'ecriture du fichier de l'alignement multiple*/
/*avec les sequences reordonnees suivant leurs groupes    */
/*                                                        */
/**********************************************************/

void ecritureFichierClustersAlignement(char *nomFichier,int longueurAlignement,
				       int nbIndividus,individu_t **individus,
				       char **sequences,int nbClusters)
{
  /*declaration des variables*/
  int i,j,k,*taillesClusters,nbIndividusTotal;
  individu_t ***clusters;
  FILE *file;
  /*fin declaration des variables*/

  /*allocation memoire*/
  taillesClusters=(int *)malloc(sizeof(int)*(nbClusters+1));
  clusters=(individu_t ***)malloc(sizeof(individu_t *)*(nbClusters+1));
  /*fin allocation memoire*/

  nbIndividusTotal=nbIndividus;
  for(i=0;i<nbClusters+1;i++)
    {
      clusters[i]=(individu_t **)malloc(sizeof(individu_t *)*nbIndividusTotal);
    }
  for(i=0;i<nbClusters+1;i++)
    {
      taillesClusters[i]=0;
    }

  for(i=0;i<nbIndividus;i++)
    {
      if(individus[i]->cluster>-1)
	{
	  clusters[individus[i]->cluster][taillesClusters[individus[i]->cluster]]=individus[i];
	  (taillesClusters[individus[i]->cluster])++;
	}
      else
	{
	  clusters[nbClusters][taillesClusters[nbClusters]]=individus[i];
	  (taillesClusters[nbClusters])++;
	}
    }

  file=fopen(nomFichier,"w");
  for(i=0;i<nbClusters;i++)
    {
      fprintf(file,">GROUP_%d\n",i+1);
      for(j=0;j<longueurAlignement;j++)
	{
	  fprintf(file,"-");
	}
      fprintf(file,"\n");
      for(j=0;j<taillesClusters[i];j++)
	{
	  fprintf(file,">%s\n",clusters[i][j]->nom);
	  for(k=0;k<longueurAlignement;k++)
	    {
	      if(sequences[clusters[i][j]->id][k]=='O')
		{
		  fprintf(file,"-");
		}
	      else
		{
		  fprintf(file,"%c",sequences[clusters[i][j]->id][k]);
		}
	    }
	  fprintf(file,"\n");
	}
    }
  if(taillesClusters[nbClusters]>0)
    {
      fprintf(file,">UNCLUSTERED\n");
      for(j=0;j<taillesClusters[nbClusters];j++)
	{
	  fprintf(file,">%s\n",clusters[nbClusters][j]->nom);
	  for(k=0;k<longueurAlignement;k++)
	    {
	      if(sequences[clusters[nbClusters][j]->id][k]=='O')
		{
		  fprintf(file,"-");
		}
	      else
		{
		  fprintf(file,"%c",sequences[clusters[nbClusters][j]->id][k]);
		}
	    }
	  fprintf(file,"\n");
	}
    }
  fclose(file);
 
  /*desallocation memoire*/
  free(taillesClusters);
  for(i=0;i<nbClusters;i++)
    {
      free(clusters[i]);
    }
  free(clusters);
  /*fin desallocation memoire*/
}


/***************************************************/
/*                                                 */
/*Procedure d'ecriture des etiquettes d'un ensemble*/
/*                                                 */
/***************************************************/

void ecritureFichierEtiquettes(int tailleTestSet,individu_t **testSet,char *nomFichier) {
	/*declaration de variables*/
	int i;
	FILE *fichier;
	/*fin declaration de variables*/
	printf("fichier : %s\n",nomFichier);
	fichier=fopen(nomFichier,"w");
	for(i=0;i<tailleTestSet;i++) {
	  fprintf(fichier,"%s %d\n",testSet[i]->nom,testSet[i]->cluster);	
	}
	fclose(fichier);
}
