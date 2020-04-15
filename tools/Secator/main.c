#include "main.h"
#include "lecfi.h"
#include "ecrfi.h"
#include "divide.h"
#include "alignment.h"
#include "treebionj.h"
#include "secatorbionj.h"
#include "treeward.h"
#include "treecut.h"

/**********************************************/
/*                                            */
/*Programme principal du package de clustering*/
/*                                            */
/**********************************************/

int main(int argc,char **argv)
{
  /*declaration de variables*/
  int i,j,typeDonnees,clusteringMethod,weighting,longueurAlignement,nbIndividus;
  int dimension,nbClusters,nbIterationsMax,nbVoisins,individusTousIdentiques,indiceIndividu;
  int nbClustersMax=-1;
  char *buf,*fichierEntree,**sequences,nomFichier[TAILLE_NOM];
  char *outputFileClustering,*outputFile_tfa;
  char *arbreTexte;
  double distanceMin,seuilDissimilarite,**distances;
  individu_t **individus;
  noeud_t *arbre;
  /*fin declaration de variables*/


  /************************************/
  /*Traitement des parametres d'entree*/ 
  /************************************/


  /*parametres par defaut*/
  weighting=NON;
  
  outputFileClustering=NULL;
  outputFile_tfa=NULL;
  nbIterationsMax=-1;
  clusteringMethod=-1;
  /*fin parametres par defaut*/

  if(argc==1) {
     helpMessageSecator();
  }

  fichierEntree=strdup(argv[1]);

  if((strcmp(fichierEntree,"-help")==0)||(strcmp(fichierEntree,"-h")==0)) {
      helpMessageSecator();
    }
  if(fopen(fichierEntree,"r")==NULL) {
  	printf("File %s does not exist !\n",fichierEntree);
  	exit(0);
  }
  for(i=2;i<argc;i++) {
      if(strstr(argv[i],"-dt=")!=NULL)
	{
	  buf=strstr(argv[i],"=")+1;
	  if(strcmp(buf,"alignment")==0)
	    {
	      typeDonnees=ALIGNEMENT;
	    }
	  else if(strcmp(buf,"coordinates")==0)
	    {
	      typeDonnees=COORDONNEES;
	  }
	  else if(strcmp(buf,"distances")==0)
	    {
	      typeDonnees=DISTANCES;
	    }
	  else if(strcmp(buf,"newick")==0) 
	    {
	      typeDonnees=ARBRE;
	    }
	  else 
	    {
	      printf("bad argument for -data_types !\n");
	      exit(1);
	    }
	}
      else if(strstr(argv[i],"-cm=")!=NULL)
	{
	  buf=strstr(argv[i],"=")+1;
	  if(strcmp(buf,"hierar")==0)
	    {
	      clusteringMethod=WARD;
	    }
	  else if(strcmp(buf,"bionj")==0)
	    {
	      clusteringMethod=BIONJ;
	    }
	  else 
	    {
	      printf("bad argument for -clusteringMethod !\n");
	      exit(1);
	    }
	}
      else if(strstr(argv[i],"-weighting")!=NULL)
	{
	  weighting=OUI;
	}
       else if(strstr(argv[i],"-otfa=")!=NULL)
	{
	  outputFile_tfa=strdup(strstr(argv[i],"=")+1);
	} 
       else if(strstr(argv[i],"-oclu=")!=NULL)
	{
	  outputFileClustering=strdup(strstr(argv[i],"=")+1);
	} 
    }
  

  /****************************************/
  /*Fin traitement des parametres d'entree*/
  /****************************************/

  /*********************/
  /*                   */
  /*LECTURE DES DONNEES*/
  /*                   */
  /*********************/

  /*lecture des donnees*/
  if(typeDonnees==ALIGNEMENT)
    {
      /*lecture du fichier de l'alignement multiple*/
      lectureFichierAlignement(fichierEntree,&longueurAlignement,&nbIndividus,
			       &individus,&sequences);
     
     
	  /*nbDimensions=nbIndividus;*/
	  /*calcul des distances entre les sequences*/
	  calculDistancesSequences(nbIndividus,individus,longueurAlignement,sequences);

    }
  else if(typeDonnees==ARBRE)
    {
      /*lecture de l'arbre phylogenetique*/
      lectureArbrePhylogenetique(fichierEntree,&arbreTexte);
     
      /*(re)construction de l'arbre phylogenetique*/
       constructionArbre(&nbIndividus,arbreTexte,&arbre);

       nbIndividus--;
       individus=(individu_t **)malloc(sizeof(individu_t *)*nbIndividus);
       for(i=0;i<nbIndividus;i++) {
	 individus[i]=(individu_t *)malloc(sizeof(individu_t));
	 individus[i]->valeurs=(double *)malloc(sizeof(double)*nbIndividus);
	 individus[i]->cluster=-1;
	 individus[i]->nom = (char *)malloc(sizeof(char)*TAILLE_NOM);
	 individus[i]->description=NULL;
       }
       
       /*enracinement de l'arbre*/ 
       enracinementArbre(nbIndividus,&arbre); 

       /*association des individus aux noeuds*/
       indiceIndividu=0;
       associationIndividusNoeuds2(&indiceIndividu,individus,arbre);

       /*calcul des distances*/
       calculDistances(nbIndividus,individus,arbre);

      /*the number of clusters is found by Secator*/
      secatorBionj(nbIndividus,arbre,&nbClusters,weighting);
    } else if(typeDonnees==COORDONNEES) {
   /*lecture du fichier de coordonnees*/
      lectureFichierCoordonnees(&nbIndividus,&dimension,&individus,fichierEntree,typeDonnees);
  }
  else if(typeDonnees==DISTANCES) {
    /*lecture du fichier de distances*/
    lectureFichierDistances(&nbIndividus,dimension,&individus,fichierEntree,&distances);
  }


  /************/
  /*          */
  /*CLUSTERING*/
  /*          */
  /************/
  
  /*clustering method : KMEANS */
  if(clusteringMethod==WARD) {
      classificationWard(nbIndividus,dimension,individus,typeDonnees,&arbre);
      
      /*the number of clusters is found by Secator*/
	secatorTree(nbIndividus,arbre,&seuilDissimilarite,&nbClusters);

	  /*creation des groupes*/
	  decoupageArbreSeuilDissimilarite(nbIndividus,arbre,seuilDissimilarite,&nbClusters);	
    }  /*clustering method : bionj*/
  else if(clusteringMethod==BIONJ)
    {
      /*on verifie que les individus ne sont pas tous identiques*/
      individusTousIdentiques=OUI;
      for(i=0;i<nbIndividus;i++)
	{
	  for(j=i+1;j<nbIndividus;j++)
	    {
	      if(individus[i]->valeurs[j]>DBL_MIN)
		{
		  individusTousIdentiques=NON;
		}
	    }
	}

      if(individusTousIdentiques==OUI)
	{
	  nbClusters=1;
	  for(i=0;i<nbIndividus;i++)
	    {
	      individus[i]->cluster=0;
	    }
	}
      else
	{
	  if((typeDonnees==ALIGNEMENT)||(typeDonnees=DISTANCES))
	    {
	      creationArbreBionj(nbIndividus,individus,fichierEntree,&arbre);
	    }
	  /*the number of clusters is found by Secator*/
	  
	      secatorBionj(nbIndividus,arbre,&nbClusters,weighting);

	  }
    }
 
  /************************/
  /*                      */
  /*ECRITURE DES RESULTATS*/
  /*                      */
  /************************/

  /*ecriture du fichier resultat presentant les clusters*/
  if(outputFileClustering==NULL)
    {
      strcpy(nomFichier,fichierEntree);
      if(strstr(nomFichier,".")!=NULL)
	{
	  sprintf(strstr(nomFichier,"."),".clu");
	}
      else
	{
	  strcat(nomFichier,".clu");
	}
    }
  else
    {
      strcpy(nomFichier,outputFileClustering);
    }


//  printf("avant ecriture fichier clusters\n");
    ecritureFichierClusters(nbIndividus,dimension,individus,nbClusters,nomFichier);
  

  if(typeDonnees==ALIGNEMENT)
    {
      if(outputFile_tfa==NULL)
	{
	  strcpy(nomFichier,fichierEntree);
	  if(strstr(nomFichier,".")!=NULL)
	  {  
	    sprintf(strstr(nomFichier,"."),"2.tfa");
	   }
          else
	  {
	    strcat(nomFichier,"2.tfa");
	  }
	}
      else
	{
	  strcpy(nomFichier,outputFile_tfa);
	}
     
    
      /*ecriture du fichier de l'alignement multiple avec les sequences*/
      /*reordonnees suivant leurs groupes*/
      ecritureFichierClustersAlignement(nomFichier,longueurAlignement,nbIndividus,
					   individus,sequences,nbClusters);
    }

  /*desallocation memoire*/
  for(i=0;i<nbIndividus;i++)
    {
      free(individus[i]->nom);
      free(individus[i]->valeurs);
      free(individus[i]);
    }
  free(individus);
  /*fin desallocation memoire*/
 // printf("main exit\n");
    exit(0);
}

void helpMessageSecator() {
      printf("Program usage : secator file options\n\n");
      printf("********************OPTIONS BELOW********************\n\n");
      printf("-dt=[newick|alignment|distances] (dt stands for data_type)\n");
      printf("-cm=[hierar|bionj] (cm stands for clusteringMethod))\n");
      printf("[-otfa=outputFile for alignment]\n");
      printf("[-oclu=outputFile for clustering]\n");
      exit(0);
}


