#include "main.h"
#include "treebionj.h"
#include "bionj.h"

/******************************************************************************************/
/*                                                                                        */
/*Procedure construisant l'arbre bionj a partir d'un alignement de sequence dans sequences*/
/*                                                                                        */
/******************************************************************************************/

void creationArbreBionj(int nbIndividus,individu_t **individus,
			char *fichierEntree,noeud_t **racine)
{
  /*declaration de variables*/
  char nomFichierEntreeBionj[TAILLE_NOM],nomFichierSortieBionj[TAILLE_NOM],*arbre;
  int entierBidon;
  /*fin declaration de variables*/

  strcpy(nomFichierEntreeBionj,fichierEntree);
  if(strstr(nomFichierEntreeBionj,".")!=NULL)
    {
      sprintf(strstr(nomFichierEntreeBionj,"."),".dst");
    }
  else
    {
      strcat(nomFichierEntreeBionj,".dst");
    }
  
  strcpy(nomFichierSortieBionj,fichierEntree);
  if(strstr(nomFichierSortieBionj,".")!=NULL)
    {
      sprintf(strstr(nomFichierSortieBionj,"."),".nj");
    }
  else
    {
      strcat(nomFichierSortieBionj,".nj");
    }

  /*ecriture du fichier.dst des distances en vue du calcul de phylogenie*/
  ecritureFichierDistancesSequences(nomFichierEntreeBionj,nbIndividus,individus);

  /*calcul de l'arbre phylogenetique par bionj*/
  bionj(nomFichierEntreeBionj,nomFichierSortieBionj);

  /*lecture de l'arbre phylogenetique*/
  lectureArbrePhylogenetique(nomFichierSortieBionj,&arbre);

  /*(re)construction de l'arbre phylogenetique*/
  constructionArbre(&entierBidon,arbre,racine);
  
  /*enracinement de l'arbre*/ 
  enracinementArbre(nbIndividus,racine); 

  /*association des individus aux noeuds*/
  associationIndividusNoeuds(nbIndividus,individus,*racine); 
}

/*******************************************************/
/*                                                     */
/*Procedure de calcul des distances a partir de l'arbre*/
/*                                                     */
/*******************************************************/

void calculDistances(int nbIndividus,individu_t **individus,noeud_t *racine) {
  /*declaration de variables*/
  int i,j;
  noeud_t *noeud1,*noeud2,*noeudCourant,*premierAncetreCommun;
  double distance=0;
  /*fin declaration de variables*/

  for(i=0;i<nbIndividus;i++) {
    for(j=i+1;j<nbIndividus;j++) {
      /*on demarque tous les noeuds*/
      marquageNoeuds(racine);

      /*on cherche le premier ancetre commun*/
      noeud1=individus[i]->noeud;
      noeud2=individus[j]->noeud;

      noeudCourant=noeud1;
      noeudCourant->marquage=1;
      while(noeudCourant!=racine) {
	noeudCourant=noeudCourant->copains[0];
	noeudCourant->marquage=1;
      }

      noeudCourant=noeud2;
      while(noeudCourant->marquage==0) {
	noeudCourant=noeudCourant->copains[0];

      }
      premierAncetreCommun=noeudCourant;

      /*on recalcule la distance entre chaque paire d'individus*/
      distance=0;
      noeudCourant=noeud1;
      while(noeudCourant!=premierAncetreCommun) {
	distance+=noeudCourant->distances[0];
	noeudCourant=noeudCourant->copains[0];
      }

      noeudCourant=noeud2;
      while(noeudCourant!=premierAncetreCommun) {
	distance+=noeudCourant->distances[0];
	noeudCourant=noeudCourant->copains[0];
      }

      individus[i]->valeurs[individus[j]->id]=distance;
      individus[j]->valeurs[individus[i]->id]=distance;
    }
  }
}

/**************************************/
/*                                    */
/*Procedure de marquage a 0 des noeuds*/
/*                                    */
/**************************************/

void marquageNoeuds(noeud_t *noeud)
{
  /*declaration de variables*/
  int i;
  char indice[100];
  /*fin declaration de variables*/

  noeud->marquage=0;
	
  if(noeud->copains[1]!=NULL)
    {
      marquageNoeuds(noeud->copains[1]);
      marquageNoeuds(noeud->copains[2]);
    }
}


/**************************************************/
/*                                                */
/*Procedure d'association des individus aux noeuds*/
/*                                                */
/**************************************************/

void associationIndividusNoeuds(int nbIndividus,individu_t **individus,noeud_t *noeud)
{
  /*declaration de variables*/
  int i;
  char indice[100];
  /*fin declaration de variables*/

  for(i=0;i<nbIndividus;i++)
    {
      sprintf(indice,"%d",individus[i]->id);
      if(strcmp(noeud->etiquette,indice)==0)
	{
	  noeud->individu=individus[i];
	  individus[i]->noeud=noeud;
	  break;
	}
    }

  if(noeud->copains[1]!=NULL)
    {
      associationIndividusNoeuds(nbIndividus,individus,noeud->copains[1]);
      associationIndividusNoeuds(nbIndividus,individus,noeud->copains[2]);
    }
}

/**************************************************/
/*                                                */
/*Procedure d'association des individus aux noeuds*/
/*                                                */
/**************************************************/

void associationIndividusNoeuds2(int *indiceIndividu,individu_t **individus,noeud_t *noeud)
{
  /*declaration de variables*/
  int i;
  char indice[100];
  /*fin declaration de variables*/

  if(noeud->copains[1]!=NULL) {
    associationIndividusNoeuds2(indiceIndividu,individus,noeud->copains[1]);
    associationIndividusNoeuds2(indiceIndividu,individus,noeud->copains[2]);
  } else {
    individus[*indiceIndividu]->id=atoi(noeud->etiquette);
    individus[*indiceIndividu]->nom=strdup(noeud->etiquette);
    noeud->individu=individus[*indiceIndividu];
    individus[*indiceIndividu]->noeud=noeud;
    (*indiceIndividu)++;
  }
}


/*************************************************************/
/*                                                           */
/*Procedure d'enracinement de l'arbre sur la sequence fictive*/
/*                                                           */
/*************************************************************/

void enracinementArbre(int nbIndividus,noeud_t **racine)
{
  /*declaration des variables*/
  noeud_t *noeudSequenceFictive;
  int nbNoeudsInternes,nbFeuilles;
  /*fin declaration des variables*/

  /*procedure de recherche de la sequence fictive dans l'arbre*/
  if(chercheNoeudSequenceFictive(nbIndividus,(*racine)->copains[1],
				    &noeudSequenceFictive)==PAS_TROUVE)
    {
      chercheNoeudSequenceFictive(nbIndividus,(*racine)->copains[2],
				     &noeudSequenceFictive);
    }

  *racine=noeudSequenceFictive->copains[0];

  if((*racine)->copains[0]==NULL)
    {
      if((*racine)->copains[1]==NULL)
	{
	  *racine=(*racine)->copains[2];
	}
      else
	{
	  *racine=(*racine)->copains[1];
	}
    }
  else
    {
      if((*racine)->copains[1]==noeudSequenceFictive)
	{
	  (*racine)->copains[1]=NULL;
	}
      else
	{
	  (*racine)->copains[2]=NULL;
	}
     
      /*on restructure l'arbre en partant de la sequence fictive*/
      restructurationArbre(*racine,NULL,0);
    }

  nbNoeudsInternes=0;
  nbFeuilles=0;

  /*renommage des noeuds internes comme le nombre de noeuds a change*/
  renommageNoeudsInternes(nbIndividus,*racine,&nbNoeudsInternes,&nbFeuilles); 
}


/**********************************************************************/
/*                                                                    */
/*Procedure de renommage des noeuds internes car apres restructuration*/
/*le nombre de noeuds n'est plus le meme                              */
/*                                                                    */
/**********************************************************************/

void renommageNoeudsInternes(int nbIndividus,noeud_t *noeud,int *nbNoeudsInternes,
			     int *nbFeuilles)
{
  if(noeud->copains[1]!=NULL)
    {
      noeud->numero=*nbNoeudsInternes+nbIndividus;
      sprintf(noeud->etiquette,"etiquette%d",*nbNoeudsInternes+nbIndividus);
      (*nbNoeudsInternes)++;
      renommageNoeudsInternes(nbIndividus,noeud->copains[1],nbNoeudsInternes,nbFeuilles);
      renommageNoeudsInternes(nbIndividus,noeud->copains[2],nbNoeudsInternes,nbFeuilles);
    }
  else
    {
      noeud->numero=*nbFeuilles;
      (*nbFeuilles)++;
    }
}

/****************************************************************************/
/*                                                                          */
/*Procedure de restructuration de l'arbre et de calcul de la nouvelle racine*/
/*                                                                          */
/****************************************************************************/

void restructurationArbre(noeud_t *noeud,noeud_t *nouveauPere,double nouvelleDistance)
{
  /*declaration des variables*/
  noeud_t *ancienPere,*ancienFilsGauche,*ancienFilsDroite;
  double ancienneDistance;
  /*fin declaration des variables*/

  ancienPere=noeud->copains[0];
  ancienFilsGauche=noeud->copains[1];
  ancienFilsDroite=noeud->copains[2];
  ancienneDistance=noeud->distances[0];
  
  if(noeud->copains[0]->copains[0]!=NULL)
    {
      if(noeud->copains[1]==nouveauPere)
	{
	  noeud->copains[1]=ancienPere;
	}
      else
	{
	  noeud->copains[2]=ancienPere;
	}
      noeud->copains[0]=nouveauPere;
      noeud->distances[0]=nouvelleDistance;
      restructurationArbre(ancienPere,noeud,ancienneDistance);
    }
  /*on se trouve au niveau d'un fils de l'ancienne racine*/
  /*ancienPere pointant sur l'ancienne racine*/
  else
    { 
      if(noeud->copains[1]==nouveauPere)
	{
	  if(ancienPere->copains[1]==noeud)
	    {
	      noeud->copains[1]=ancienPere->copains[2];
	      ancienPere->copains[2]->copains[0]=noeud;
	      ancienPere->copains[2]->distances[0]+=ancienneDistance;
	    }
	  else
	    {
	      noeud->copains[1]=ancienPere->copains[1];
	      ancienPere->copains[1]->copains[0]=noeud;
		  ancienPere->copains[1]->distances[0]+=ancienneDistance;
	    }
	}
      else
	{
	  if(ancienPere->copains[1]==noeud)
	    {
	      noeud->copains[2]=ancienPere->copains[2];
	      ancienPere->copains[2]->copains[0]=noeud;
	      ancienPere->copains[2]->distances[0]+=ancienneDistance;
	      ancienPere->copains[1]->distances[0]+=ancienneDistance;
	    }
	  else
	    {
	      noeud->copains[2]=ancienPere->copains[1];
	      ancienPere->copains[1]->copains[0]=noeud;
	      
	    }
	}
      
      noeud->copains[0]=nouveauPere;
      noeud->distances[0]=nouvelleDistance;
    }
}

/*******************************************************************/
/*                                                                 */
/*Procedure d'affichage de l'arbre (essentiellement pour debuggage)*/
/*                                                                 */
/*******************************************************************/

void affichageArbre(noeud_t *noeud)
{
  printf("%s(%d) : \n",noeud->etiquette,noeud->numero);

  if(noeud->copains[1]!=NULL)
    {
      printf("fils1 : %s(%d) \n",noeud->copains[1]->etiquette,noeud->copains[1]->numero);
      printf("fils2 : %s(%d) \n",noeud->copains[2]->etiquette,noeud->copains[2]->numero);

      affichageArbre(noeud->copains[1]);
      affichageArbre(noeud->copains[2]);
    }
}

/************************************************************/
/*                                                          */
/*Procedure de recherche de la sequence fictive dans l'arbre*/
/*                                                          */
/************************************************************/

int chercheNoeudSequenceFictive(int nbIndividus,noeud_t *noeud,
				   noeud_t **noeudSequenceFictive)
{
  if(noeud->copains[1]==NULL)
    {
      if(strcmp(noeud->etiquette,"fictive_sequence")==0)
	{
	  *noeudSequenceFictive=noeud;
	  return TROUVE;
	}
      else
	{
	  return PAS_TROUVE;
	}
    }
  else
    {
      if(chercheNoeudSequenceFictive(nbIndividus,noeud->copains[1],
					noeudSequenceFictive)==TROUVE)
	{
	  return TROUVE;
	}
      else if(chercheNoeudSequenceFictive(nbIndividus,noeud->copains[2],
					     noeudSequenceFictive)==TROUVE)
	{
	  return TROUVE;
	}
      else
	{
	  return PAS_TROUVE;
	}
    }
}


/*****************************************************/
/*                                                   */
/*Procedure de construction de l'arbre phylogenetique*/
/*                                                   */
/*****************************************************/

void constructionArbre(int *nbIndividus,char *arbre,noeud_t **racine)
{
  /*declaration des variables*/
  int i,j,nombreNoeuds,compteurIndividus;
  char *positionParentheseOuvrante,*positionParentheseFermante;
  char *positionPremiersDoublePoints,*positionDeuxiemesDoublePoints;
  char *positionTroisiemesDoublePoints,*positionPremiereVirgule,*positionSecondeVirgule;
  char etiquetteGauche[200],etiquetteDroite[200],etiquetteMilieu[200];
  char nouvelleEtiquette[200],distanceGauche[200],distanceDroite[200],distanceMilieu[200];
  noeud_t *noeudGauche,*noeudMilieu,*noeudDroite,*noeudPere,**noeudsCourants;
  /*fin declaration des variables*/

  /*on cherche le nombre d'individus en cherchant le nombre de double-points 
    non precedes par une parenthese fermante*/
  *nbIndividus=0;
  
  /*on cherche les premiers double_points*/
  positionPremiersDoublePoints=arbre;
  while((positionPremiersDoublePoints=strstr(&(positionPremiersDoublePoints[1]),":"))!=NULL) {
    if(positionPremiersDoublePoints[-1]!=')') {
      (*nbIndividus)++;
    }
  }

  /*allocation de memoire*/
  noeudsCourants=(noeud_t **)malloc(sizeof(noeud_t *)*(*nbIndividus));
  /*fin allocation de memoire*/

  /*construction de l'arbre*/
  nombreNoeuds=0;
  compteurIndividus=0;
  for(i=0;i<*nbIndividus-3;i++)
    {
      /*on cherche la premiere parenthese fermante*/
      positionParentheseFermante=strstr(arbre,")");

      j=-1;
      while(positionParentheseFermante[j]!='(')
	{
	  j--;
	}
      positionParentheseOuvrante=&(positionParentheseFermante[j]);
      positionPremiersDoublePoints=strstr(positionParentheseOuvrante,":");
      positionPremiereVirgule=strstr(positionPremiersDoublePoints,",");
      positionDeuxiemesDoublePoints=strstr(positionPremiereVirgule,":");

      copiePartieString(positionParentheseOuvrante,positionPremiersDoublePoints,
			  etiquetteGauche);

      copiePartieString(positionPremiereVirgule,positionDeuxiemesDoublePoints,
			  etiquetteDroite);

      copiePartieString(positionPremiersDoublePoints,positionPremiereVirgule,
			  distanceGauche);
      
      copiePartieString(positionDeuxiemesDoublePoints,positionParentheseFermante,
			  distanceDroite);

      if(strstr(etiquetteGauche,"etiquette")==NULL)
	{
	  /*allocation memoire*/
	  noeudGauche=(noeud_t *)malloc(sizeof(noeud_t));
	  /*fin allocation memoire*/

	  noeudGauche->numero=compteurIndividus;compteurIndividus++;
	  strcpy(noeudGauche->etiquette,etiquetteGauche);
	  
	  noeudGauche->copains[1]=NULL;
	  noeudGauche->copains[2]=NULL;
	  noeudGauche->copain1=NULL;
	  noeudGauche->copain2=NULL;
	}
      else
	{
	  for(j=0;j<*nbIndividus;j++)
	    {
	      if(strcmp((noeudsCourants[j])->etiquette,etiquetteGauche)==0)
		{
		  noeudGauche=noeudsCourants[j];
		  break;
		}
	    }
	}

      if(strstr(etiquetteDroite,"etiquette")==NULL)
	{
	  /*allocation memoire*/
	  noeudDroite=(noeud_t *)malloc(sizeof(noeud_t));
	  /*fin allocation memoire*/

	  noeudDroite->numero=compteurIndividus;compteurIndividus++;
	  strcpy(noeudDroite->etiquette,etiquetteDroite);

	  noeudDroite->copains[1]=NULL;
	  noeudDroite->copains[2]=NULL;
	  noeudDroite->copain1=NULL;
	  noeudDroite->copain2=NULL;
	}
      else
	{
	  for(j=0;j<*nbIndividus;j++)
	    {
	      if(strcmp((noeudsCourants[j])->etiquette,etiquetteDroite)==0)
		{
		  noeudDroite=noeudsCourants[j];
		  break;
		}
	    }
	}
      
      noeudGauche->distances[0]=atof(distanceGauche);
      noeudDroite->distances[0]=atof(distanceDroite);
      
      sprintf(nouvelleEtiquette,"etiquette%d",*nbIndividus+nombreNoeuds);

      /*allocation memoire*/
      noeudPere=(noeud_t *)malloc(sizeof(noeud_t));
      /*fin allocation memoire*/

      strcpy(noeudPere->etiquette,nouvelleEtiquette);
      noeudPere->numero=*nbIndividus+nombreNoeuds;
      noeudPere->copains[1]=noeudGauche;
      noeudPere->copains[2]=noeudDroite;
      noeudPere->distances[1]=atof(distanceGauche);
      noeudPere->distances[2]=atof(distanceDroite);
      noeudGauche->copains[0]=noeudPere;
      noeudDroite->copains[0]=noeudPere;
      
      noeudsCourants[nombreNoeuds]=noeudPere;

      /*on remplace les deux noeuds par un nouveau noeud dans la chaine de caracteres arbre*/
      remplacePartieString(&arbre,positionParentheseOuvrante,positionParentheseFermante,
			     nouvelleEtiquette);

      nombreNoeuds++;
    }

  positionPremiersDoublePoints=strstr(arbre,":");
  positionPremiereVirgule=strstr(positionPremiersDoublePoints,",");
  positionDeuxiemesDoublePoints=strstr(positionPremiereVirgule,":");
  positionSecondeVirgule=strstr(positionDeuxiemesDoublePoints,",");
  positionTroisiemesDoublePoints=strstr(positionSecondeVirgule,":");
  positionParentheseFermante=strstr(arbre,")");

  copiePartieString(arbre,positionPremiersDoublePoints,etiquetteGauche);
  copiePartieString(positionPremiereVirgule,positionDeuxiemesDoublePoints,etiquetteMilieu);
  copiePartieString(positionSecondeVirgule,positionTroisiemesDoublePoints,etiquetteDroite);
  copiePartieString(positionPremiersDoublePoints,positionPremiereVirgule,
		      distanceGauche);
  copiePartieString(positionDeuxiemesDoublePoints,positionSecondeVirgule,
		      distanceMilieu);
  copiePartieString(positionTroisiemesDoublePoints,positionParentheseFermante,
		      distanceDroite);

  if(strstr(etiquetteGauche,"etiquette")==NULL)
    {
      /*allocation memoire*/
      noeudGauche=(noeud_t *)malloc(sizeof(noeud_t));
      /*fin allocation memoire*/
      
      noeudGauche->numero=compteurIndividus;compteurIndividus++;
      strcpy(noeudGauche->etiquette,etiquetteGauche);

      noeudGauche->copains[1]=NULL;
      noeudGauche->copains[2]=NULL;
      noeudGauche->copain1=NULL;
      noeudGauche->copain2=NULL;
    }
  else
    {
      for(j=0;j<*nbIndividus;j++)
	{
	  if(strcmp((noeudsCourants[j])->etiquette,etiquetteGauche)==0)
	    {
	      noeudGauche=noeudsCourants[j];
	      break;
	    }
	}
    }

   if(strstr(etiquetteMilieu,"etiquette")==NULL)
     {
       /*allocation memoire*/
       noeudMilieu=(noeud_t *)malloc(sizeof(noeud_t));
       /*fin allocation memoire*/
       
       noeudMilieu->numero=compteurIndividus;compteurIndividus++;
       strcpy(noeudMilieu->etiquette,etiquetteMilieu);
       noeudMilieu->copains[1]=NULL;
       noeudMilieu->copains[2]=NULL;
       noeudMilieu->copain1=NULL;
       noeudMilieu->copain2=NULL;
     }
   else
     {
       for(j=0;j<*nbIndividus;j++)
	 {
	   if(strcmp((noeudsCourants[j])->etiquette,etiquetteMilieu)==0)
	     {
	       noeudMilieu=noeudsCourants[j];
	       break;
	     }
	 }
     }
   
   if(strstr(etiquetteDroite,"etiquette")==NULL)
     {
       /*allocation memoire*/
       noeudDroite=(noeud_t *)malloc(sizeof(noeud_t));
       /*fin allocation memoire*/
       
       noeudDroite->numero=compteurIndividus;compteurIndividus++;
       strcpy(noeudDroite->etiquette,etiquetteDroite);
       noeudDroite->copains[1]=NULL;
       noeudDroite->copains[2]=NULL;
       noeudDroite->copain1=NULL;
       noeudDroite->copain2=NULL;
     }
   else
     {
       for(j=0;j<*nbIndividus;j++)
	 {
	   if(strcmp((noeudsCourants[j])->etiquette,etiquetteDroite)==0)
	     {
	       noeudDroite=noeudsCourants[j];
	       break;
	     }
	 }
     }

   /*allocation memoire*/
   noeudPere=(noeud_t *)malloc(sizeof(noeud_t));
   /*fin allocation memoire*/
   
   sprintf(nouvelleEtiquette,"etiquette%d",*nbIndividus+nombreNoeuds);

   strcpy(noeudPere->etiquette,nouvelleEtiquette);
  
   noeudPere->numero=*nbIndividus+nombreNoeuds;
   nombreNoeuds++;

   noeudGauche->distances[0]=atof(distanceGauche);
   noeudMilieu->distances[0]=atof(distanceMilieu);
   noeudDroite->distances[0]=atof(distanceDroite);
   noeudPere->distances[0]=0;
   

   noeudPere->copains[1]=noeudGauche;
   noeudPere->copains[2]=noeudMilieu;
   noeudPere->distances[1]=atof(distanceGauche);
   noeudPere->distances[2]=atof(distanceMilieu);

   noeudGauche->copains[0]=noeudPere;
   noeudMilieu->copains[0]=noeudPere;

   /*allocation memoire*/
   *racine=(noeud_t *)malloc(sizeof(noeud_t));
   /*fin allocation memoire*/

   sprintf(nouvelleEtiquette,"etiquette%d",*nbIndividus+nombreNoeuds);

   strcpy((*racine)->etiquette,nouvelleEtiquette);
   (*racine)->numero=*nbIndividus+nombreNoeuds;
   (*racine)->copains[1]=noeudPere;
   (*racine)->copains[2]=noeudDroite;
   (*racine)->copains[0]=NULL;
   /*   (*racine)->distances[0]=0;
	(*racine)->distances[1]=0;
	(*racine)->distances[2]=atof(distanceDroite);*/

   noeudDroite->copains[0]=*racine;
   noeudPere->copains[0]=*racine;

   /*desallocation memoire*/
   free(noeudsCourants);
   /*fin desallocation memoire*/
}

/*******************************************************************************/
/*                                                                             */
/*Procedure de copie d'une portion de chaine de caracteres delimitee de maniere*/
/*stricte par deux pointeurs sur char dans une autre chaine de caracteres      */
/*                                                                             */
/*******************************************************************************/

void copiePartieString(char *debut,char *fin,char *destination)
{
  /*declaration des variables*/
  int i;
  /*fin declaration des variables*/

  i=1;
  while(&(debut[i])!=fin)
    {
      destination[i-1]=debut[i];
      i++;
    } 
  destination[i-1]=0;
} 



/******************************************************************************/
/*                                                                            */
/*Procedure de remplacement d'une partie d'une chaine de caracteres, delimitee*/
/*de maniere large par debut et fin, par une autre chaine de caracteres       */
/*                                                                            */
/******************************************************************************/

void remplacePartieString(char **chaine,char *debut,char *fin,char *chaineDeRemplacement)
{
  /*declaration des variables*/
  char *premierePartie,*secondePartie;
  /*fin declaration des variables*/

  /*allocation memoire*/
  premierePartie=(char *)malloc(sizeof(char)*(strlen(*chaine)+1));
  secondePartie=(char *)malloc(sizeof(char)*(strlen(*chaine)+1));
  /*fin allocation memoire*/

  strcpy(premierePartie,*chaine); 
  premierePartie[strlen(*chaine)-strlen(debut)]=0;
  strcpy(secondePartie,&(fin[1]));

  *chaine=realloc(*chaine,sizeof(char)*(strlen(premierePartie)+strlen(chaineDeRemplacement)+
					strlen(secondePartie)+1));
  (*chaine)[0]=0;
  strcat(*chaine,premierePartie);
  strcat(*chaine,chaineDeRemplacement);
  strcat(*chaine,secondePartie);
  
  /*desallocation memoire*/
  free(premierePartie);
  free(secondePartie);
  /*fin desallocation memoire*/
}

