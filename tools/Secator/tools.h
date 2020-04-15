/******************************************************/
/*                                                    */
/*Procedure de calcul de la distance entre deux points*/
/*                                                    */
/******************************************************/

static inline double calculDistanceCarre(int nbDimensions,double *point1,double *point2)
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

  return distance;
}

/******************************************************/
/*                                                    */
/*Procedure de calcul de la distance entre deux points*/
/*                                                    */
/******************************************************/

double calculDistance(int nbDimensions,double *point1,double *point2);


/***************************************************************************/
/*                                                                         */
/*Procedure de regression lineaire sur les mesures (renvoie le coefficient)*/
/*                                                                         */
/***************************************************************************/

double regressionLineaire(int n,double **mesures);

/*****************************************************/
/*                                                   */
/*renvoie la valeur absolue de a                     */
/*                                                   */
/*****************************************************/
double valeurAbsolue(double a);
 
/*****************************************************/
/*                                                   */
/*renvoie le minimum de a et de b                    */
/*                                                   */
/*****************************************************/
int min(int a,int b);

/****************************************/
/*                                      */
/*Procedure de tri par ordre decroissant*/
/*                                      */
/****************************************/

void triRapide(double *valeurs,int gauche,int droite);

/******************************/
/*                            */
/*sous-procedure de tri_rapide*/ 
/*                            */
/******************************/

void echangerValeursPourTri(double *valeurs,int element1,int element2);

/****************************************/
/*                                      */
/*Procedure de tri par ordre decroissant*/
/*                                      */
/****************************************/

void triRapideValeursPositions(double *valeurs,int *positions,int gauche,int droite);

/******************************/
/*                            */
/*sous-procedure de tri_rapide*/ 
/*                            */
/******************************/

void echangerValeursPositionsPourTri(double *valeurs,int *positions,int element1,int element2);


/*****************************************************/
/*                                                   */
/*Procedure de calcul de la factorielle d'un entier n*/
/*                                                   */
/*****************************************************/

int calculFactorielle(int n);

/***************************************************************/
/*                                                             */
/*Procedure de calcul de la factorielle "impaire" d'un entier n*/
/*                                                             */
/***************************************************************/

int calculFactorielleImpaire(int n);
