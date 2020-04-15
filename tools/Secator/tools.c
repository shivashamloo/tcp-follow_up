#include "main.h"
#include "tools.h"


/******************************************************/
/*                                                    */
/*Procedure de calcul de la distance entre deux points*/
/*                                                    */
/******************************************************/

double calculDistance(int nbDimensions,double *point1,double *point2)
{
  /*declaration des variables*/
  int i;
  double distance,temp;
  /*fin declaration des variables*/

  distance=0;
  for(i=0;i<nbDimensions;i++)
    {
      temp=point1[i]-point2[i];
      distance+=temp*temp;
    }
  distance=sqrt((double)distance);

  return distance;
}


/***************************************************************************/
/*                                                                         */
/*Procedure de regression lineaire sur les mesures (renvoie le coefficient)*/
/*                                                                         */
/***************************************************************************/

double regressionLineaire(int n,double **mesures)
{
  /*declaration des variables*/
  int i;
  double moyennex,moyenney,covariance,variancex;
  /*fin declaration des variables*/

  /*calcul des moyennes*/
  moyennex=0;
  moyenney=0;

  for(i=0;i<n;i++)
    {
      moyennex+=mesures[i][0];
      moyenney+=mesures[i][1];
    }
  
  moyennex/=(double)n;
  moyenney/=(double)n;

  /*calcul de la variance de x*/
  variancex=0;

  for(i=0;i<n;i++)
    {
      variancex+=pow((double)(moyennex-mesures[i][0]),2.0);
    }
  variancex/=(double)n;

  /*calcul de la covariance entre x et y*/
  covariance=0;
  
  for(i=0;i<n;i++)
    {
      covariance+=(mesures[i][0]-moyennex)*(mesures[i][1]-moyenney);
    }
  covariance/=(double)n;

  return covariance/variancex;
}

/*****************************************************/
/*                                                   */
/*renvoie la valeur absolue de a                     */
/*                                                   */
/*****************************************************/
double valeurAbsolue(double a)
{
  if(a>0)
    {
      return a;
    }
  else
    {
      return -a;
    }
}

 
/*****************************************************/
/*                                                   */
/*renvoie le minimum de a et de b                    */
/*                                                   */
/*****************************************************/
int min(int a,int b)
{
  if(a<b)
    {
      return a;
    }
  else
    {
      return b;
    }
}

/****************************************/
/*                                      */
/*Procedure de tri par ordre decroissant*/
/*                                      */
/****************************************/

void triRapide(double *valeurs,int gauche,int droite)
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
      if(valeurs[elementSuivant]>=valeurs[indiceSeparateur])
	{
	  echangerValeursPourTri(valeurs,indiceSeparateur+1,elementSuivant);
	  echangerValeursPourTri(valeurs,indiceSeparateur,indiceSeparateur+1);
	  indiceSeparateur++;
	}
      elementSuivant++;
    }
  
  triRapide(valeurs,gauche,indiceSeparateur-1);
  triRapide(valeurs,indiceSeparateur+1,droite);
}

/******************************/
/*                            */
/*sous-procedure de tri_rapide*/ 
/*                            */
/******************************/

void echangerValeursPourTri(double *valeurs,int element1,int element2)
{
  /*declaration des variables*/
  double tempVal;
  /*fin de declaration des variables*/

  tempVal=valeurs[element1];
  valeurs[element1]=valeurs[element2];
  valeurs[element2]=tempVal;
}

/****************************************/
/*                                      */
/*Procedure de tri par ordre decroissant*/
/*                                      */
/****************************************/

void triRapideValeursPositions(double *valeurs,int *positions,int gauche,int droite)
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
      if(valeurs[elementSuivant]>=valeurs[indiceSeparateur])
	{
	  echangerValeursPositionsPourTri(valeurs,positions,indiceSeparateur+1,
					  elementSuivant);
	  echangerValeursPositionsPourTri(valeurs,positions,indiceSeparateur,
					  indiceSeparateur+1);
	  indiceSeparateur++;
	}
      elementSuivant++;
    }
  
  triRapideValeursPositions(valeurs,positions,gauche,indiceSeparateur-1);
  triRapideValeursPositions(valeurs,positions,indiceSeparateur+1,droite);
}

/******************************/
/*                            */
/*sous-procedure de tri_rapide*/ 
/*                            */
/******************************/

void echangerValeursPositionsPourTri(double *valeurs,int *positions,int element1,int element2)
{
  /*declaration des variables*/
  int tempPosition;
  double tempVal;
  /*fin de declaration des variables*/

  tempVal=valeurs[element1];
  valeurs[element1]=valeurs[element2];
  valeurs[element2]=tempVal;

  tempPosition=positions[element1];
  positions[element1]=positions[element2];
  positions[element2]=tempPosition;
}

/*****************************************************/
/*                                                   */
/*Procedure de calcul de la factorielle d'un entier n*/
/*                                                   */
/*****************************************************/

int calculFactorielle(int n)
{
  /*declaration de variables*/
  int i,resultat;
  /*fin declaration de variables*/

  resultat=1;
  for(i=2;i<=n;i++)
    {
      resultat*=i;
    }
  return resultat;
}

/***************************************************************/
/*                                                             */
/*Procedure de calcul de la factorielle "impaire" d'un entier n*/
/*                                                             */
/***************************************************************/

int calculFactorielleImpaire(int n)
{
  /*declaration de variables*/
  int i,resultat;
  /*fin declaration de variables*/

  resultat=1;
  for(i=2;i<=n;i++)
    {
      resultat*=2*i-1;
    }
  return resultat;
}
