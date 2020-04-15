#include "main.h"

/*********************************************************************/
/*                                                                   */
/*Procedure de calcul des pourcentages d'identite entre les sequences*/
/*                                                                   */
/*********************************************************************/

void calculIdentitesSequences(int nbIndividus,individu_t **individus,
				int longueurAlignement,char **sequences)
{
  /*declaration des variables*/
  int i,j,k,nbIdentites,nbPositionsCommunes;
  /*fin declaration des variables*/

  for(i=0;i<nbIndividus;i++)
    {
      individus[i]->valeurs[i]=1.0;
      for(j=i+1;j<nbIndividus;j++)
	{
	  nbPositionsCommunes=0;
	  nbIdentites=0;
	  for(k=0;k<longueurAlignement;k++)
	    {
	       if((sequences[i][k]!='O')&&(sequences[j][k]!='O')&&(sequences[i][k]!='X')&&(sequences[j][k]!='X'))
		{
		  nbPositionsCommunes++;
		  if(sequences[i][k]==sequences[j][k])
		    {
		      nbIdentites++;
		    }
		}
	    }
	  if(nbPositionsCommunes==0)
	    {
	      individus[i]->valeurs[j]=0.0;
	    }
	  else
	    {
	      individus[i]->valeurs[j]=
		(double)nbIdentites/(double)nbPositionsCommunes;
	    }
	  individus[j]->valeurs[i]=individus[i]->valeurs[j];
	}
    }

  /*printf("%d\n",nbIndividus);
    for(i=0;i<nbIndividus;i++)
    {   
    printf("%s ",individus[i]->nom);
    for(j=0;j<nbIndividus;j++)
    {
    printf("%.2f ",individus[i]->valeursTraitees[j]);
    }
    printf("\n");
    }  
    exit(0);*/
}

/*********************************************************************/
/*                                                                   */
/*Procedure de calcul des pourcentages d'identite entre les sequences*/
/*                                                                   */
/*********************************************************************/

void calculDistancesSequences(int nbIndividus,individu_t **individus,
			      int longueurAlignement,char **sequences)
{
  /*declaration des variables*/
  int i,j,k,nbIdentites,nbPositionsCommunes;
  /*fin declaration des variables*/

  for(i=0;i<nbIndividus;i++)
    {
      individus[i]->valeurs[i]=0.0;
      for(j=i+1;j<nbIndividus;j++)
	{
	  nbPositionsCommunes=0;
	  nbIdentites=0;
	  for(k=0;k<longueurAlignement;k++)
	    {
	       if((sequences[i][k]!='O')&&(sequences[j][k]!='O')&&(sequences[i][k]!='X')&&(sequences[j][k]!='X'))
		{
		  nbPositionsCommunes++;
		  if(sequences[i][k]==sequences[j][k])
		    {
		      nbIdentites++;
		    }
		}
	    }
	  if(nbPositionsCommunes==0)
	    {
	      individus[i]->valeurs[j]=1.0;
	    }
	  else
	    {
	      individus[i]->valeurs[j]=
		1.0-(double)nbIdentites/(double)nbPositionsCommunes;
	    }
	  individus[j]->valeurs[i]=individus[i]->valeurs[j];
	}
    }

  /*for(i=0;i<nbIndividus;i++) {
    for(j=i+1;j<nbIndividus;j++) {
      printf("%.4f\n",individus[j]->valeursTraitees[i]);
    }
    }*/
}

/******************************************************************************/
/*                                                                            */
/*Procedure de calcul des distances entre les sequences en les stockant a part*/
/*                                                                            */
/******************************************************************************/

void calculDistancesApartSequences(int nbIndividus,int longueurAlignement,
				     char **sequences,double **distances)
{
  /*declaration des variables*/
  int i,j,k,nbIdentites,nbPositionsCommunes;
  /*fin declaration des variables*/

  for(i=0;i<nbIndividus;i++)
    {
      distances[i][i]=0.0;
      for(j=i+1;j<nbIndividus;j++)
	{
	  nbPositionsCommunes=0;
	  nbIdentites=0;
	  for(k=0;k<longueurAlignement;k++)
	    {
	       if((sequences[i][k]!='O')&&(sequences[j][k]!='O')&&(sequences[i][k]!='X')&&(sequences[j][k]!='X'))
		{
		  nbPositionsCommunes++;
		  if(sequences[i][k]==sequences[j][k])
		    {
		      nbIdentites++;
		    }
		}
	    }
	  if(nbPositionsCommunes==0)
	    {
	      distances[i][j]=1.0;
	    }
	  else
	    {
	      distances[i][j]=1.0-(double)nbIdentites/(double)nbPositionsCommunes;
	    }
	  distances[j][i]=distances[i][j];
	}
    }
}
