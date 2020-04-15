#include "main.h"
#include "weighting.h"

/*********************************************/
/*                                           */
/*Procedure de calcul des poids des sequences*/
/*                                           */
/*********************************************/

void calculPoidsSequences(noeud_t *noeud,double poids)
{
  /*declaration des variables*/
  int nbFeuillesGauche,nbFeuillesDroite;
  /*fin declaration des variables*/

  if(noeud->copains[1]!=NULL)
    {
      nbFeuillesGauche=0;
      nbFeuillesDroite=0;
      calculNombreFeuilles(noeud->copains[1],&nbFeuillesGauche);
      calculNombreFeuilles(noeud->copains[2],&nbFeuillesDroite);
      if(noeud->copains[1]==NULL)
	{
	  /*on est a la racine*/
	  calculPoidsSequences(noeud->copains[1],0);
	  calculPoidsSequences(noeud->copains[2],0);
	}
      else
	{
	  calculPoidsSequences(noeud->copains[1],poids+noeud->distances[0]/
				 (double)nbFeuillesGauche);
	  calculPoidsSequences(noeud->copains[2],poids+noeud->distances[0]/
				 (double)nbFeuillesDroite);
	}
    }
  else
    {
      noeud->poids=poids+noeud->distances[0];
    }
}

/***************************************************************/
/*                                                             */
/*Procedure de calcul du nombre de feuilles a partir d'un noeud*/
/*                                                             */
/***************************************************************/

void calculNombreFeuilles(noeud_t *noeud,int *nbFeuilles)
{
  if(noeud->copains[1]==NULL)
    {
      (*nbFeuilles)++;
    }
  else
    {
      calculNombreFeuilles(noeud->copains[1],nbFeuilles);
      calculNombreFeuilles(noeud->copains[2],nbFeuilles);
    }
  
}



