/*********************************************************************/
/*                                                                   */
/*Procedure de calcul des pourcentages d'identite entre les sequences*/
/*                                                                   */
/*********************************************************************/

void calculIdentitesSequences(int nbIndividus,individu_t **individus,
			      int longueurAlignement,char **sequences);

/*********************************************************************/
/*                                                                   */
/*Procedure de calcul des pourcentages d'identite entre les sequences*/
/*                                                                   */
/*********************************************************************/

void calculDistancesSequences(int nbIndividus,individu_t **individus,
				int longueurAlignement,char **sequences);

/******************************************************************************/
/*                                                                            */
/*Procedure de calcul des distances entre les sequences en les stockant a part*/
/*                                                                            */
/******************************************************************************/

void calculDistancesApartSequences(int nbIndividus,int longueurAlignement,
				   char **sequences,double **distances);
