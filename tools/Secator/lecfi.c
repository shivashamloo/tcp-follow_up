#include "main.h"
#include "lecfi.h"

/***********************************************************/
/*                                                         */
/*Procedure de lecture du fichier d'un arbre phylogenetique*/
/*(au format bionj)                                        */
/*                                                         */
/***********************************************************/

void lectureArbrePhylogenetique(char *nomFichier,char **arbre)
{
  /*declaration des variables*/
  FILE *fichierEntree;
  int nbCaracteres;
  char poubelle[TAILLE_POUBELLE];
  /*fin declaration des variables*/

  /*lecture du fichier de l'arbre*/
  fichierEntree=fopen(nomFichier,"r");
    
  if(fichierEntree==NULL)
    {
      printf("file %s does not exist !\n",nomFichier);
      exit(0);
    }
  
  /*allocation de memoire*/
  *arbre=(char *)malloc(sizeof(char)*TAILLE_POUBELLE);
  if(*arbre==NULL)
    {
      printf("echec de l'allocation de memoire\n");
      exit(0);
    }
  /*fin allocation de memoire*/
  
  (*arbre)[0]=0;
  nbCaracteres=0;

  while(fgets(poubelle,TAILLE_POUBELLE,fichierEntree)!=NULL)
    {
      nbCaracteres+=strlen(poubelle);

      /*reallocation de memoire*/
      (*arbre)=(char *)realloc(*arbre,sizeof(char)*(nbCaracteres+1));
      if(*arbre==NULL)
	{
	  printf("echec de l'allocation de memoire\n");
	  exit(0);
	}
      /*fin reallocation de memoire*/
      
      strcat(*arbre,poubelle);
    }
  fclose(fichierEntree);
  /*fin lecture du fichier de l'arbre*/
}

/********************************************/
/*                                          */
/*Procedure de lecture du fichier contenant */
/*les coordonnees des individus a classifier*/
/*Les individus sont stockes dans un tableau*/
/*de pointeurs sur structure                */
/*                                          */
/********************************************/

void lectureFichierCoordonnees(int *nbIndividus,int *nbDimensions,individu_t ***individus,
				 char *nomFichierEntree,int typeDonnees)
{
  /*declaration des variables*/
  int i,j,retourLecture;
  FILE *fichier;
  char nomIndividu[TAILLE_NOM],ligne[TAILLE_MAX_LIGNE],ligne1[TAILLE_MAX_LIGNE],ligne2[TAILLE_MAX_LIGNE];
  double somme;
  /*fin declaration des variables*/

  /*premiere lecture pour compter le nombre d'individus et de colonnes*/
  fichier=fopen(nomFichierEntree,"r");
  retourLecture=(int)fgets(ligne1,TAILLE_MAX_LIGNE,fichier);
  retourLecture=(int)fgets(ligne2,TAILLE_MAX_LIGNE,fichier);
  
  /*on compte le nombre d'entrees dans la premiere et la deuxieme ligne;
    si le nombre d'entrees differe on efface la premiere ligne (c'est une ligne de commentaires)*/
  char * pch;
  int nbEntreesLigne1,nbEntreesLigne2;
  
  nbEntreesLigne1=0;
  pch=strtok(ligne1," \t");
  while (pch != NULL) {
    nbEntreesLigne1++;
    pch = strtok (NULL, " \t");
  }

  nbEntreesLigne2=0;
  pch=strtok(ligne2," \t");
  while (pch != NULL) {
    nbEntreesLigne2++;
    pch = strtok (NULL, " \t");
  }

  if(nbEntreesLigne1==(nbEntreesLigne2-1)) {
  	/*on ignore la premiere ligne*/
	*nbIndividus=1;
  } else if(nbEntreesLigne1==nbEntreesLigne2) {
 	*nbIndividus=2;
  } else {
        printf("format problem !\n");
  	exit(1);
  }

  while(fgets(ligne,TAILLE_MAX_LIGNE,fichier)!=NULL) {
     (*nbIndividus)++;
  }
  *nbDimensions=nbEntreesLigne2-1;
  fclose(fichier);  
 //  printf("%d\n",nbEntreesLigne2);
  //printf("nbIndividus : %d ; nbDimensions : %d\n",*nbIndividus,*nbDimensions);
  /*deuxieme lecture de la ligne de description des colonnes*/
  fichier=fopen(nomFichierEntree,"r");
   
  /*allocation memoire*/
  *individus = (individu_t **)malloc(sizeof(individu_t *)*(*nbIndividus));
  for(i=0;i<*nbIndividus;i++) {
      (*individus)[i]=(individu_t *)malloc(sizeof(individu_t));
      (*individus)[i]->id=i;
      (*individus)[i]->cluster=-1;
      (*individus)[i]->valeurs = (double *)malloc(sizeof(double)*(*nbDimensions));
      (*individus)[i]->cluster=AUCUN_GROUPE;
     }
  /*fin allocation memoire*/
  for(i=0;i<*nbIndividus;i++) {
      retourLecture=(int)fscanf(fichier,"%s",nomIndividu);
      (*individus)[i]->nom = (char *)malloc(sizeof(char)*(strlen(nomIndividu)+1));
      strcpy((*individus)[i]->nom,nomIndividu);
      somme=0;
      for(j=0;j<*nbDimensions;j++) {
	  retourLecture=(int)fscanf(fichier,"%lf",&((*individus)[i]->valeurs[j]));
	}

      retourLecture=(int)fgets(ligne,TAILLE_MAX_LIGNE,fichier);
      (*individus)[i]->description=(char *)malloc(sizeof(char)*(strlen(ligne)+10));
      strcpy((*individus)[i]->description,ligne);
    }
  fclose(fichier);
}

/************************************************/
/*                                              */
/*Procedure de lecture du fichier contenant     */
/*les distances entre les individus a classifier*/
/*Les individus sont stockes dans un tableau    */
/*de pointeurs sur structure                    */
/*                                              */
/************************************************/

void lectureFichierDistances(int *nbIndividus,int nbDimensions,individu_t ***individus,
			       char *nomFichierEntree,double ***distances)
{
  /*declaration des variables*/
  int i,j,retourLecture;
  FILE *fichier;
  char nomIndividu[TAILLE_NOM],ligne[TAILLE_MAX_LIGNE];
  double somme;
  /*fin declaration des variables*/

  fichier=fopen(nomFichierEntree,"r");
  retourLecture=(int)fscanf(fichier,"%d",nbIndividus);
  retourLecture=(int)fgets(ligne,TAILLE_MAX_LIGNE,fichier);

  /*lecture de la ligne de description des colonnes*/
  retourLecture=(int)fgets(ligne,TAILLE_MAX_LIGNE,fichier);
   
  /*allocation memoire*/
  *individus = (individu_t **)malloc(sizeof(individu_t *)*(*nbIndividus));
  *distances=(double **)malloc(sizeof(double *)*(*nbIndividus));
  for(i=0;i<*nbIndividus;i++) {
      (*individus)[i]=(individu_t *)malloc(sizeof(individu_t));
      (*individus)[i]->id=i;
      (*individus)[i]->cluster=-1;
      (*individus)[i]->valeurs = (double *)malloc(sizeof(double)*nbDimensions);
      (*individus)[i]->cluster=AUCUN_GROUPE;
      (*distances)[i]=(double *)malloc(sizeof(double)*(*nbIndividus));
(*individus)[i]->nbIndividusSimilaires=0;
    }
  /*fin allocation memoire*/
  for(i=0;i<*nbIndividus;i++) {
      retourLecture=(int)fscanf(fichier,"%s",nomIndividu);
      (*individus)[i]->nom = (char *)malloc(sizeof(char)*(strlen(nomIndividu)+1));
      strcpy((*individus)[i]->nom,nomIndividu);
      for(j=0;j<*nbIndividus;j++)
	{
	  retourLecture=(int)fscanf(fichier,"%lf",&((*distances)[i][j]));
	}
      retourLecture=(int)fgets(ligne,TAILLE_MAX_LIGNE,fichier);
      (*individus)[i]->description=(char *)malloc(sizeof(char)*(strlen(ligne)+10));
      strcpy((*individus)[i]->description,ligne);
    }
  fclose(fichier);
}

/**********************************************************/
/*                                                        */
/*Procedure de lecture du fichier d'un alignement multiple*/
/*                                                        */
/**********************************************************/

void lectureFichierAlignement(char *nomAlignement,int *longueurAlignement,int *nbIndividus,
				individu_t ***individus,char ***sequences)
{  
  /*declaration des variables*/
  int i,typeFichier,retourLecture;
  FILE *fichier;
  char poubelle[TAILLE_POUBELLE];
  /*fin declaration des variables*/
  
  
  /*on regarde si on a affaire a un fichier.msf ou a un fichier.tfa*/
  fichier=fopen(nomAlignement,"r");
  typeFichier=FICHIER_TFA;
  retourLecture=(int)fgets(poubelle,TAILLE_POUBELLE,fichier);
  if(poubelle[0]=='>')
    {
      typeFichier=FICHIER_TFA;
    }
  else
    {
      typeFichier=FICHIER_INCONNU;
      if((strstr(poubelle,"MSF:")!=NULL)&&(strstr(poubelle,"..")!=NULL))
	{
	  typeFichier=FICHIER_MSF;
	}
      while(fgets(poubelle,TAILLE_POUBELLE,fichier)!=NULL)
	{
	  if((strstr(poubelle,"MSF:")!=NULL)&&(strstr(poubelle,"..")!=NULL))
	    {
	      typeFichier=FICHIER_MSF;
	      break;
	    }
	}
    }
  fclose(fichier);
  if(typeFichier==FICHIER_INCONNU)
    {
      printf("File type must be either tfa or msf !\n");
      exit(0);
    }

  if(typeFichier==FICHIER_TFA)
    {       
      /*lecture du fichier.tfa*/
      lectureFichierTfa(nomAlignement,longueurAlignement,nbIndividus,sequences,
			  individus);
    }
  else
    {
      /*lecture du fichier.msf*/
      lectureFichierMsf(nomAlignement,longueurAlignement,nbIndividus,sequences,
			  individus);      
    }
}

/***************************************/
/*                                     */
/*Procedure de lecture d'un fichier.tfa*/
/*                                     */
/***************************************/

void lectureFichierTfa(char *nomAlignement,int *longueurAlignement,
		       int *nbIndividus,char ***sequences,individu_t ***individus)
{
  /*declaration des variables*/
  FILE *fichierEntree;
  int i,j,compteur,longueurSequenceCourante,retourLecture;
  char nomFichierTfa[200],poubelle[TAILLE_POUBELLE];
  /*fin declaration des variables*/

  sprintf(nomFichierTfa,"%s",nomAlignement);

  /*premiere lecture pour compter les sequences et calculer la longueur de l'alignement*/
  fichierEntree=fopen(nomFichierTfa,"r");
  
  if(fichierEntree==NULL)
    {
      printf("file %s does not exist !\n",nomFichierTfa);
      exit(0);
    }

  *nbIndividus=0;
  longueurSequenceCourante=0;
  *longueurAlignement=0;
  while(fgets(poubelle,TAILLE_POUBELLE,fichierEntree)!=NULL)
    {
      if(poubelle[0]=='>')
	{
	  (*nbIndividus)++;
	  if(longueurSequenceCourante>*longueurAlignement)
	    {
	      *longueurAlignement=longueurSequenceCourante;
	    }
	  longueurSequenceCourante=0;
	}
      else
	{
	  longueurSequenceCourante+=strlen(poubelle)-1;
	}
    }

  if(longueurSequenceCourante>*longueurAlignement)
    {
      *longueurAlignement=longueurSequenceCourante;
    }
  fclose(fichierEntree);
  /*fin premiere lecture pour compter les sequences et calculer la longueur de l'alignement*/

  /*allocation memoire*/
  *sequences=(char **)malloc(sizeof(char *)*(*nbIndividus+1));
  for(i=0;i<*nbIndividus+1;i++)
    {
      (*sequences)[i]=(char *)malloc(sizeof(char)*(*longueurAlignement));
    }
  *individus = (individu_t **)malloc(sizeof(individu_t *)*(*nbIndividus));
  for(i=0;i<*nbIndividus;i++)
    {
      (*individus)[i]=(individu_t *)malloc(sizeof(individu_t));
      (*individus)[i]->id = i;
      (*individus)[i]->cluster=-1;
      (*individus)[i]->valeurs = (double *)malloc(sizeof(double)*(*nbIndividus));
      (*individus)[i]->nom = (char *)malloc(sizeof(char)*TAILLE_NOM);
      (*individus)[i]->description=NULL;
(*individus)[i]->nbIndividusSimilaires=0;
    }
  /*fin allocation memoire*/

  /*seconde lecture pour lire l'alignement proprement dit*/
  fichierEntree=fopen(nomFichierTfa,"r");
  
  if(fichierEntree==NULL)
    {
      printf("file %s does not exist !\n",nomFichierTfa);
      exit(0);
    }
  
  retourLecture=(int)fgets(poubelle,TAILLE_POUBELLE,fichierEntree);
  for(i=0;i<*nbIndividus;i++)
    {
      sscanf(&(poubelle[1]),"%s",(*individus)[i]->nom);

      compteur=0;
      retourLecture=(int)fgets(poubelle,TAILLE_POUBELLE,fichierEntree);
      while(poubelle[0]!='>')
	{
	  for(j=0;j<strlen(poubelle)-1;j++)
	    {
	      if((poubelle[j]-'A'<0)||(poubelle[j]-'A'>24))
		{
		  (*sequences)[i][compteur]='O';
		}
	      else if((poubelle[j]=='B')||(poubelle[j]=='J')||(poubelle[j]=='O')||
		      (poubelle[j]=='U'))
		{
		  (*sequences)[i][compteur]='O';
		}
	      else
		{
		  (*sequences)[i][compteur]=poubelle[j];
		}
	      compteur++;
	    }
	  if(fgets(poubelle,TAILLE_POUBELLE,fichierEntree)==NULL)
	    {
	      break;
	    }
	}
      for(j=compteur;j<*longueurAlignement;j++)
	{
	  (*sequences)[i][j]='O';
	}
    }

  fclose(fichierEntree);
  /*fin seconde lecture pour lire l'alignement proprement dit*/  

  for(i=0;i<*longueurAlignement;i++)
    {
      (*sequences)[*nbIndividus][i]='P';
    }
}

/***************************************/
/*                                     */
/*Procedure de lecture d'un fichier.msf*/
/*                                     */
/***************************************/

void lectureFichierMsf(char *nomAlignement,int *longueurAlignement,
			 int *nbIndividus,char ***sequences,individu_t ***individus)
{
  /*declaration des variables*/
  FILE *fichierEntree;
  int i,j,position,compteurSequences,positionAGauche,nbColonnes,retourLecture;
  char residu,nomFichierMsf[200],poubelle[TAILLE_POUBELLE];
  char nomPremiereSequence[200],motLu[200];
  /*fin declarations des variables*/

  /*premiere lecture du fichier pour compter le nombre de sequences*/
  sprintf(nomFichierMsf,"%s",nomAlignement);

  fichierEntree=fopen(nomFichierMsf,"r");

  if(fichierEntree==NULL)
    {
      printf("file %s does not exist !\n",nomFichierMsf);
      exit(0);
    }

  retourLecture=(int)fgets(poubelle,TAILLE_POUBELLE,fichierEntree);  
  while(strstr(poubelle,"Name")==NULL)
    {
      retourLecture=(int)fgets(poubelle,TAILLE_POUBELLE,fichierEntree);  
    }
  sscanf(strstr(poubelle,"Name:")+5,"%s",nomPremiereSequence);
  sscanf(strstr(poubelle,"Len:")+4,"%d",longueurAlignement);

  *nbIndividus=0;
  while(strstr(poubelle,"Name")!=NULL)
    {
      (*nbIndividus)++;
      retourLecture=(int)fgets(poubelle,TAILLE_POUBELLE,fichierEntree);
    }
  fclose(fichierEntree);
  /*fin premiere lecture du fichier pour compter le nombre de sequences*/

  /*debut allocation de memoire*/
  *sequences=(char **)malloc(sizeof(char *)*(*nbIndividus+1));
  if(*sequences==NULL)
    {
      printf("probleme d'allocation memoire\n");
      exit(0);
    }
  for(i=0;i<*nbIndividus+1;i++)
    {
      (*sequences)[i]=(char *)malloc(sizeof(char)*(*longueurAlignement));
      if((*sequences)[i]==NULL)
	{
	  printf("probleme d'allocation memoire\n");
	  exit(0);
	}
    }
  *individus = (individu_t **)malloc(sizeof(individu_t *)*(*nbIndividus));
  for(i=0;i<*nbIndividus;i++)
    {
      (*individus)[i]=(individu_t *)malloc(sizeof(individu_t));
      (*individus)[i]->id = i;
      (*individus)[i]->cluster = -1; 
      (*individus)[i]->valeurs = (double *)malloc(sizeof(double)*(*nbIndividus));
      (*individus)[i]->nom = (char *)malloc(sizeof(char)*TAILLE_NOM);
      (*individus)[i]->description=NULL;
(*individus)[i]->nbIndividusSimilaires=0;     
}
  /*fin allocation de memoire*/

  /*deuxieme lecture*/
  fichierEntree=fopen(nomFichierMsf,"r");
  
  retourLecture=(int)fgets(poubelle,TAILLE_POUBELLE,fichierEntree);  
  while(strstr(poubelle,"Name")==NULL)
    {
      retourLecture=(int)fgets(poubelle,TAILLE_POUBELLE,fichierEntree);  
    }

  i=0;
  while(strstr(poubelle,"Name")!=NULL)
    {
      sscanf(strstr(poubelle,"Name")+5,"%s",(*individus)[i]->nom);
      i++;
      retourLecture=(int)fgets(poubelle,TAILLE_POUBELLE,fichierEntree);
    }

  positionAGauche=0;
  while(positionAGauche<*longueurAlignement)
    {
      retourLecture=(int)fscanf(fichierEntree,"%s",motLu);
      while(strcmp(motLu,nomPremiereSequence)!=0)
	{
	  retourLecture=(int)fscanf(fichierEntree,"%s",motLu);
	}
      position=positionAGauche;
     
      /*  for(i=0;(i<50)&&(i+position_a_gauche<*longueurAlignement);i++)*/
      while(position<*longueurAlignement)
	{
	  retourLecture=(int)fscanf(fichierEntree,"%c",&residu);
	  if(residu=='\n')
	    {
	      break;
	    }
	  else if(residu==' ')
	    {
	    }
	  else if((residu-'A'<0)||(residu-'A'>24))
	    {
	      /*O signifie gap*/
	      (*sequences)[0][position]='O';
	      position++; 
	    }
	  else if((residu=='B')||(residu=='J')||(residu=='O')||
		  (residu=='U'))
	    {
	      /*O signifie gap*/
	      (*sequences)[0][position]='O';
	      position++; 
	    }
	  else
	    {
	      (*sequences)[0][position]=residu;
	      position++;
	    }
	}
      for(compteurSequences=1;compteurSequences<*nbIndividus;
	  compteurSequences++)
	{
	  position=positionAGauche;
	  retourLecture=(int)fscanf(fichierEntree,"%s",poubelle);
	  /*	  for(i=0;(i<50)&&(i+position_a_gauche<*longueurAlignement);i++)*/
	  while(position<*longueurAlignement)
	    {
	      retourLecture=(int)fscanf(fichierEntree,"%c",&residu);
	      if(residu=='\n')
		{
		  break;
		}
	      else if(residu==' ')
		{
		}
	      else if((residu-'A'<0)||(residu-'A'>24))
		{
		  /*O signifie gap*/
		  (*sequences)[compteurSequences][position]='O';
		  position++; 
		}
	      else if((residu=='B')||(residu=='J')||(residu=='O')||
		      (residu=='U'))
		{
		  /*O signifie gap*/
		  (*sequences)[compteurSequences][position]='O';
		  position++; 
		}
	      else
		{
		  (*sequences)[compteurSequences][position]=residu;
		  position++;
		}
	    }
	}
       positionAGauche=position;
    }

  fclose(fichierEntree);
  for(i=0;i<*longueurAlignement;i++)
    {
      (*sequences)[*nbIndividus][i]='P';
    }
  /*fin deuxieme lecture*/
}

/*lecture du fichier des etiquettes designant les clusters*/
void lectureEtiquettes(int nbIndividus,individu_t **individus,char *nomFichier) {
  /*declaration des variables*/
  int i,tempi,retourLecture;
  FILE *fichier;
  /*fin declaration des variables*/

  fichier=fopen(nomFichier,"r"); 
  for(i=0;i<nbIndividus;i++) {
    retourLecture=(int)fscanf(fichier,"%d %d",&tempi,&(individus[i]->cluster));
  }
  fclose(fichier);
}
