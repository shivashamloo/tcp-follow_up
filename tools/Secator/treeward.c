#include "main.h"
#include "treeward.h"

/********************************************************************/
/*                                                                  */
/*Procedure de classification par la methode des voisins reciproques*/ 
/*                                                                  */
/********************************************************************/

void classificationWard(int nbIndividus,int nbDimensions,individu_t **individus,
			 int typeDonnees,noeud_t **racine)
{
  /*declaration des variables*/
  int i,j,k,l,nbRacines,graine,plusProcheVoisin,ancienneGraine,compteurNumeros;
  double **dissimilaritesCourantes,dissimilariteMinimum;
  double *poids,*historiqueDissimilarites;
  double dissimilariteMinimumParia,seuilDissimilarite;
  couple_t meilleurRegroupement;
  noeud_t *feuilles,**noeuds,*nouveauNoeud;
  /*fin de declaration des variables*/

  /*allocation de la memoire*/
  dissimilaritesCourantes=(double **)malloc(sizeof(double *)*nbIndividus);
  for(i=0;i<nbIndividus;i++)
    {
      dissimilaritesCourantes[i]=(double *)malloc(sizeof(double)*nbIndividus);
    }
  poids=(double *)malloc(sizeof(double)*nbIndividus);
  feuilles=(noeud_t *)malloc(sizeof(noeud_t)*nbIndividus);
  noeuds=(noeud_t **)malloc(sizeof(noeud_t *)*nbIndividus);

  historiqueDissimilarites=(double *)malloc(sizeof(double)*nbIndividus);
  /*fin allocation de la memoire*/

  /*initialisation des dissimilarites*/
  for(i=0;i<nbIndividus;i++)
    {
      poids[i]=1.0;
      dissimilaritesCourantes[i][i]=0;
    }

  compteurNumeros=0;
  /*initialisation des dissimilarites*/
  for(i=0;i<nbIndividus-1;i++)
    {
      for(j=i+1;j<nbIndividus;j++)
	{       
	  if((nbDimensions<=0)||(typeDonnees==DISTANCES))
	    {
	      dissimilaritesCourantes[i][j]=pow((double)(individus[i]->valeurs[j]),
						 2.0)*poids[i]*poids[j]/(poids[i]+poids[j]);
	      dissimilaritesCourantes[j][i]=dissimilaritesCourantes[i][j];
	    }
	  else
	    {
	      dissimilaritesCourantes[i][j]=0.0;
	      for(k=0;k<nbDimensions;k++)
		{
		  dissimilaritesCourantes[i][j]+=
		    pow((double)(individus[i]->valeurs[k]-
				 individus[j]->valeurs[k]),2.0);
		}
	      
	      dissimilaritesCourantes[i][j]*=poids[i]*poids[j]/(poids[i]+poids[j]);
	      dissimilaritesCourantes[j][i]=dissimilaritesCourantes[i][j];
	    }
	}
    }
    
   /*initialisation de la construction de l'arbre*/
   for(i=0;i<nbIndividus;i++)
     {
       feuilles[i].copain1=NULL;
       feuilles[i].copain2=NULL;
       feuilles[i].numero=i;
       feuilles[i].qualite=FEUILLE;
       feuilles[i].distances[1]=0;
       feuilles[i].distances[2]=0;
       feuilles[i].dissimilarite=0;
       feuilles[i].individu=individus[i];
       noeuds[i]=&(feuilles[i]);
     }
  
   /*debut de la classification hierarchique proprement dite*/
   graine=0;
   for(nbRacines=nbIndividus;nbRacines>1;nbRacines--)
    {
      dissimilariteMinimum=DBL_MAX;
      dissimilariteMinimumParia=DBL_MAX;
      if(graine==0)
	{
	  plusProcheVoisin=1;
	}
      else
	{
	  plusProcheVoisin=0;
	}
      dissimilariteMinimum=dissimilaritesCourantes[graine][plusProcheVoisin];
      if(nbRacines==2)
	{
	  meilleurRegroupement.i=0;
	  meilleurRegroupement.j=1;
	}
      else
	{
	  for(i=0;i<nbRacines;i++)
	    {
	      if((i!=graine)&&(dissimilaritesCourantes[graine][i]<dissimilariteMinimum))
		{
		  plusProcheVoisin=i;
		  dissimilariteMinimum=dissimilaritesCourantes[graine][i];
		}
	    }
	  ancienneGraine=graine;
	  graine=plusProcheVoisin;
	  if(graine==0)
	    {
	      plusProcheVoisin=1;
	    }
	  else
	    {
	      plusProcheVoisin=0;
	    }
	  dissimilariteMinimum=dissimilaritesCourantes[plusProcheVoisin][graine];

	  while(1)
	    {
	      for(i=0;i<nbRacines;i++)
		{
		  if((i!=graine)&&(dissimilaritesCourantes[graine][i]<dissimilariteMinimum))
		    {
		      plusProcheVoisin=i;
		      dissimilariteMinimum=dissimilaritesCourantes[graine][i];
		    }
		}
	      if(plusProcheVoisin==ancienneGraine)
		{
		  break;
		}
	      else
		{
		  ancienneGraine=graine;
		  graine=plusProcheVoisin;
		  if(graine==0)
		    {
		      plusProcheVoisin=1;
		    }
		  else
		    {
		      plusProcheVoisin=0;
		    }
		  dissimilariteMinimum=dissimilaritesCourantes[plusProcheVoisin][graine];
		}
	    }
	
	  if(graine<plusProcheVoisin)
	    {
	      meilleurRegroupement.i=graine;
	      meilleurRegroupement.j=plusProcheVoisin;
	    }
	  else
	    {
	      meilleurRegroupement.i=plusProcheVoisin;
	      meilleurRegroupement.j=graine;
	      graine=plusProcheVoisin;
	    }
	}

      historiqueDissimilarites[nbRacines-2]=dissimilariteMinimum;
	
      /*construction de l'arbre*/
      nouveauNoeud=(noeud_t *)malloc(sizeof(noeud_t));
     
      nouveauNoeud->copain1=noeuds[meilleurRegroupement.i];
      nouveauNoeud->copain2=noeuds[meilleurRegroupement.j];
      nouveauNoeud->numero=compteurNumeros;
      compteurNumeros++;

      nouveauNoeud->qualite=NOEUD_ORDINAIRE;
      nouveauNoeud->dissimilarite=dissimilariteMinimum;
      nouveauNoeud->distances[1]=
	dissimilariteMinimum-noeuds[meilleurRegroupement.i]->dissimilarite;
      nouveauNoeud->distances[2]=
	dissimilariteMinimum-noeuds[meilleurRegroupement.j]->dissimilarite;
      noeuds[meilleurRegroupement.i]->copains[0]=nouveauNoeud;
      noeuds[meilleurRegroupement.j]->copains[0]=nouveauNoeud;
      noeuds[meilleurRegroupement.i]=nouveauNoeud;
      for(i=meilleurRegroupement.j;i<nbRacines-1;i++)
	{
	  noeuds[i]=noeuds[i+1];
	}
      
      /*mise a jour des dissimilarites*/
      for(i=0;i<meilleurRegroupement.i;i++)
	{
	  dissimilaritesCourantes[i][meilleurRegroupement.i]=
	    ((poids[meilleurRegroupement.i]+poids[i])*
	     dissimilaritesCourantes[i][meilleurRegroupement.i]+
	     (poids[meilleurRegroupement.j]+poids[i])*
	     dissimilaritesCourantes[i][meilleurRegroupement.j]-
	     (poids[i])*
	     dissimilaritesCourantes[meilleurRegroupement.i][meilleurRegroupement.j])/
	    (poids[i]+poids[meilleurRegroupement.i]+poids[meilleurRegroupement.j]);
	  dissimilaritesCourantes[meilleurRegroupement.i][i]=
	    dissimilaritesCourantes[i][meilleurRegroupement.i];
	}
      for(i=meilleurRegroupement.i+1;i<meilleurRegroupement.j;i++)
	{
	  dissimilaritesCourantes[i][meilleurRegroupement.i]=
	    ((poids[meilleurRegroupement.i]+poids[i])*
	     dissimilaritesCourantes[i][meilleurRegroupement.i]+
	     (poids[meilleurRegroupement.j]+poids[i])*
	     dissimilaritesCourantes[i][meilleurRegroupement.j]-
	     (poids[i])*
	     dissimilaritesCourantes[meilleurRegroupement.i][meilleurRegroupement.j])/
	    (poids[i]+poids[meilleurRegroupement.i]+
	     poids[meilleurRegroupement.j]);
	  dissimilaritesCourantes[meilleurRegroupement.i][i]=
	    dissimilaritesCourantes[i][meilleurRegroupement.i];
	}
      for(i=meilleurRegroupement.j+1;i<nbRacines;i++)
	{
	  dissimilaritesCourantes[i][meilleurRegroupement.i]=
	    ((poids[meilleurRegroupement.i]+poids[i])*
	     dissimilaritesCourantes[i][meilleurRegroupement.i]+
	     (poids[meilleurRegroupement.j]+poids[i])*
	     dissimilaritesCourantes[i][meilleurRegroupement.j]-
	     (poids[i])*
	     dissimilaritesCourantes[meilleurRegroupement.i][meilleurRegroupement.j])/
	    (poids[i]+poids[meilleurRegroupement.i]+
	     poids[meilleurRegroupement.j]);
	  dissimilaritesCourantes[meilleurRegroupement.i][i]=
	    dissimilaritesCourantes[i][meilleurRegroupement.i];
	}
      dissimilaritesCourantes[meilleurRegroupement.i][meilleurRegroupement.i]=0;
      for(i=meilleurRegroupement.j;i<nbRacines-1;i++)
	{
	  for(j=0;j<nbRacines;j++)
	    {
	      dissimilaritesCourantes[i][j]=dissimilaritesCourantes[i+1][j];
	    }
	}
      for(i=meilleurRegroupement.j;i<nbRacines-1;i++)
	{
	  for(j=0;j<nbRacines-1;j++)
	    {
	      dissimilaritesCourantes[j][i]=dissimilaritesCourantes[j][i+1];
	    }
	}
      
      /*fin de la mise a jour des dissimilarites*/
      poids[meilleurRegroupement.i]+=poids[meilleurRegroupement.j];
      for(i=meilleurRegroupement.j;i<nbRacines-1;i++)
	{
	  poids[i]=poids[i+1];
	}
    }
   /*fin de la classification hierarchique*/

   *racine=noeuds[0];

   /*desallocation de la memoire*/
   for(i=0;i<nbIndividus;i++)
     {
       free(dissimilaritesCourantes[i]);
     }
   free(dissimilaritesCourantes);
  
   free(poids);
   free(noeuds);
   free(historiqueDissimilarites);
   /*fin desallocation de la memoire*/
}
