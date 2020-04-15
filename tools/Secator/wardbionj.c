#include "main.h"
#include "wardbionj.h"

/******************************************/
/*                                        */
/*Procedure de classification hierarchique*/
/*                                        */
/******************************************/

void wardSurUnArbreBionj(int nbSequences,noeud_t **adressesNoeuds,
			 noeud_t *noeudFictif,int *compteurIndicesDissimilarite,
			 double *historiqueIndicesDissimilarite,
			 noeud_t **historiqueNoeuds)
{
  /*declaration des variables*/
  int i,j,k,l,groupe1,groupe2,tempGroupe,copain1,copain2,tempCopain,unNoeudEnCommun;
  int *numerosGroupesCourants,nbGroupesCourants,kcourant,lcourant;
  double **indicesDissimilariteAnciens,**indicesDissimilariteNouveaux;
  double plusPetitIndiceDissimilarite;
  /*fin declaration des variables*/

  /*allocation memoire*/
  numerosGroupesCourants=(int *)malloc(sizeof(int)*nbSequences);
  
  indicesDissimilariteAnciens=(double **)malloc(sizeof(double *)*nbSequences);
  for(i=0;i<nbSequences;i++)
    {
      indicesDissimilariteAnciens[i]=(double *)malloc(sizeof(double)*nbSequences);
    }
  
  indicesDissimilariteNouveaux=(double **)malloc(sizeof(double *)*nbSequences);
  for(i=0;i<nbSequences;i++)
    {
      indicesDissimilariteNouveaux[i]=(double *)malloc(sizeof(double)*nbSequences);
    }
  /*fin allocation memoire*/

  /*mise a zero des pertes d'inertie inter-classe*/
  for(i=0;i<nbSequences;i++)
    { 
      adressesNoeuds[i]->dissimilarite=0;
    } 

  nbGroupesCourants=nbSequences;
  
  /*initialisation des indices de dissimilarite*/
  initialisationIndicesDissimilarite(nbSequences,indicesDissimilariteAnciens,adressesNoeuds);

  /*initialisation des numeros des groupes ;*/
  /*les numeros de groupe etant les numeros des derniers noeuds qui ont rejoint leur groupe*/
  for(i=0;i<nbSequences;i++)
    {
      numerosGroupesCourants[i]=i;
    }

  *compteurIndicesDissimilarite=0;
  while(nbGroupesCourants>2)
    {
      /*on recherche le couple de groupes (groupe1,groupe2)*/
      /*ayant le plus petit indice de dissimilarite*/
      plusPetitIndiceDissimilarite=-1.0;
      for(i=0;i<nbGroupesCourants;i++)
	{
	  for(j=i+1;j<nbGroupesCourants;j++)
	    {
	      unNoeudEnCommun=NON;
	      for(k=0;k<3;k++)
		{
		  for(l=0;l<3;l++)
		    {
		      if((adressesNoeuds[numerosGroupesCourants[i]]->copains[k]!=NULL)&&
			 (adressesNoeuds[numerosGroupesCourants[j]]->copains[l]!=NULL)&&
			 (adressesNoeuds[numerosGroupesCourants[i]]->copains[k]==
			  adressesNoeuds[numerosGroupesCourants[j]]->copains[l]))
			{
			  unNoeudEnCommun=OUI;
			  kcourant=k;
			  lcourant=l;
			  break;
			}
		    }
		  if(unNoeudEnCommun==OUI)
		    {
		      break;
		    }
		}

	      if(((indicesDissimilariteAnciens[i][j]<=plusPetitIndiceDissimilarite)||
		  (plusPetitIndiceDissimilarite<0))&&(unNoeudEnCommun==OUI))
	    {
		  plusPetitIndiceDissimilarite=indicesDissimilariteAnciens[i][j];
		  groupe1=i;
		  groupe2=j;
		  copain1=kcourant;
		  copain2=lcourant;
		}
	    }
	}
      /*fin recherche des groupes a reunir*/

      nbGroupesCourants--;
      
      adressesNoeuds[numerosGroupesCourants[groupe1]]->copains[copain1]->dissimilarite=
	plusPetitIndiceDissimilarite;

      /*      adressesNoeuds[numerosGroupesCourants[groupe1]]->distance=
	      adressesNoeuds[numerosGroupesCourants[groupe1]]->distances[copain1];
	      adressesNoeuds[numerosGroupesCourants[groupe2]]->distance=
	      adressesNoeuds[numerosGroupesCourants[groupe2]]->distances[copain2];*/

      adressesNoeuds[numerosGroupesCourants[groupe1]]->copains[copain1]->poids=
	adressesNoeuds[numerosGroupesCourants[groupe1]]->poids+
	adressesNoeuds[numerosGroupesCourants[groupe2]]->poids;

      adressesNoeuds[numerosGroupesCourants[groupe1]]->copains[copain1]->copain1=
	adressesNoeuds[numerosGroupesCourants[groupe1]];
      adressesNoeuds[numerosGroupesCourants[groupe1]]->copains[copain1]->copain2=
	adressesNoeuds[numerosGroupesCourants[groupe2]];
       
      historiqueIndicesDissimilarite[*compteurIndicesDissimilarite]=
	plusPetitIndiceDissimilarite;

      historiqueNoeuds[*compteurIndicesDissimilarite]=
	adressesNoeuds[numerosGroupesCourants[groupe1]]->copains[copain1];

      if(groupe2<groupe1)
	{
	  tempGroupe=groupe1;
	  groupe1=groupe2;
	  groupe2=tempGroupe;

	  tempCopain=copain1;
	  copain1=copain2;
	  copain2=tempCopain;
	}

      /*mise a jour des indices de dissimilarite*/
      miseAJourIndicesWard(groupe1,groupe2,adressesNoeuds,nbGroupesCourants,
			   numerosGroupesCourants,indicesDissimilariteAnciens,
			   indicesDissimilariteNouveaux);

      /*le groupe1 happe le groupe2 qui disparait*/
      for(i=groupe2;i<nbGroupesCourants;i++)
	{
	  numerosGroupesCourants[i]=numerosGroupesCourants[i+1];
	} 
      numerosGroupesCourants[groupe1]=
	adressesNoeuds[numerosGroupesCourants[groupe1]]->copains[copain1]->numero; 

      (*compteurIndicesDissimilarite)++;
    }
  
  for(i=0;i<3;i++)
    {
      if(adressesNoeuds[numerosGroupesCourants[0]]->copains[i]==
	 adressesNoeuds[numerosGroupesCourants[1]])
	{
	  copain1=i;
	  break;
	}
    }
   for(i=0;i<3;i++)
    {
      if(adressesNoeuds[numerosGroupesCourants[1]]->copains[i]==
	 adressesNoeuds[numerosGroupesCourants[0]])
	{
	  copain2=i;
	  break;
	}
    }

   strcpy(noeudFictif->etiquette,"noeudFictif");
   noeudFictif->copains[0]=NULL;
   noeudFictif->copains[1]=adressesNoeuds[numerosGroupesCourants[0]];
   noeudFictif->copains[2]=adressesNoeuds[numerosGroupesCourants[1]];
   noeudFictif->copain1=adressesNoeuds[numerosGroupesCourants[0]];
   noeudFictif->copain2=adressesNoeuds[numerosGroupesCourants[1]];
   noeudFictif->dissimilarite=indicesDissimilariteAnciens[0][1];
   /*   noeudFictif->distances[1]=adressesNoeuds[numerosGroupesCourants[0]]->distance;
	noeudFictif->distances[2]=adressesNoeuds[numerosGroupesCourants[1]]->distance;*/

   adressesNoeuds[numerosGroupesCourants[0]]->copains[copain1]=noeudFictif;
   adressesNoeuds[numerosGroupesCourants[1]]->copains[copain2]=noeudFictif;
   
   /*   adressesNoeuds[numerosGroupesCourants[0]]->distance=
	adressesNoeuds[numerosGroupesCourants[0]]->distances[copain1]/2.0;
	adressesNoeuds[numerosGroupesCourants[1]]->distance=
	adressesNoeuds[numerosGroupesCourants[1]]->distances[copain2]/2.0;*/
   
   historiqueNoeuds[*compteurIndicesDissimilarite]=noeudFictif;
   historiqueIndicesDissimilarite[*compteurIndicesDissimilarite]=noeudFictif->dissimilarite;
   (*compteurIndicesDissimilarite)++;

   /*desallocation memoire*/
   free(numerosGroupesCourants);
   
   for(i=0;i<nbSequences;i++)
     {
       free(indicesDissimilariteAnciens[i]);
     }
   free(indicesDissimilariteAnciens);
   
   for(i=0;i<nbSequences;i++)
     {
       free(indicesDissimilariteNouveaux[i]);
     }
   free(indicesDissimilariteNouveaux);
   /*fin desallocation memoire*/
}

/*********************************************************/
/*                                                       */
/*Procedure d'initialisation des indices de dissimilarite*/
/*                                                       */
/*********************************************************/

void initialisationIndicesDissimilarite(int nbIndividus,double **indicesDissimilariteAnciens,
					noeud_t **adressesNoeuds)
{
  /*declaration des variables*/
  int i,j;
  /*fin declaration des variables*/
 
  /*indices de Ward*/
  for(i=0;i<nbIndividus-1;i++)
    {
      for(j=i+1;j<nbIndividus;j++)
	{
	  indicesDissimilariteAnciens[i][j]=
	    (adressesNoeuds[i]->poids*adressesNoeuds[j]->poids)*
	    pow(adressesNoeuds[i]->individu->valeurs[adressesNoeuds[j]->individu->id],
		2.0)/(adressesNoeuds[i]->poids+adressesNoeuds[j]->poids);
	  indicesDissimilariteAnciens[j][i]=indicesDissimilariteAnciens[i][j];
	}
    } 
}


/***************************************************/
/*                                                 */
/*Procedure de mise a jour des indices de Ward     */
/*lors de la classification hierarchique ascendante*/
/*                                                 */
/***************************************************/

void miseAJourIndicesWard(int groupe1,int groupe2,noeud_t **adressesNoeuds,
			      int nbGroupesCourants,int *numerosGroupesCourants,
			      double **indicesWardAnciens,double **indicesWardNouveaux)
{
  /*declaration des variables*/
  int i,j;
  /*fin declaration des variables*/
  
  /******************************************/
  /*mise a jour des indices de Ward nouveaux*/
  /******************************************/

  /*mise a jour des indices de groupe1*/
  for(i=0;i<groupe2;i++)
    {
      if(i!=groupe1)
	{
	  indicesWardNouveaux[groupe1][i]=
	    (indicesWardAnciens[groupe1][i]*
	     (adressesNoeuds[numerosGroupesCourants[groupe1]]->poids+
	      adressesNoeuds[numerosGroupesCourants[i]]->poids)+
	     indicesWardAnciens[groupe2][i]*
	     (adressesNoeuds[numerosGroupesCourants[groupe2]]->poids+
	      adressesNoeuds[numerosGroupesCourants[i]]->poids)-
	     indicesWardAnciens[groupe1][groupe2]*
	     adressesNoeuds[numerosGroupesCourants[i]]->poids)/
	    (adressesNoeuds[numerosGroupesCourants[groupe1]]->poids+
	     adressesNoeuds[numerosGroupesCourants[groupe2]]->poids+
	     adressesNoeuds[numerosGroupesCourants[i]]->poids);
      
	  indicesWardNouveaux[i][groupe1]=indicesWardNouveaux[groupe1][i];
	}
    }
  indicesWardNouveaux[groupe1][groupe1]=0;

  for(i=groupe2;i<nbGroupesCourants;i++)
    {
      indicesWardNouveaux[groupe1][i]=
	(indicesWardAnciens[groupe1][i+1]*
	 (adressesNoeuds[numerosGroupesCourants[groupe1]]->poids+
	  adressesNoeuds[numerosGroupesCourants[i+1]]->poids)+
	 indicesWardAnciens[groupe2][i+1]*
	 (adressesNoeuds[numerosGroupesCourants[groupe2]]->poids+
	  adressesNoeuds[numerosGroupesCourants[i+1]]->poids)-
	 indicesWardAnciens[groupe1][groupe2]*
	 adressesNoeuds[numerosGroupesCourants[i+1]]->poids)/
	(adressesNoeuds[numerosGroupesCourants[groupe1]]->poids+
	 adressesNoeuds[numerosGroupesCourants[groupe2]]->poids+
	 adressesNoeuds[numerosGroupesCourants[i+1]]->poids);
      
      indicesWardNouveaux[i][groupe1]=indicesWardNouveaux[groupe1][i];
    }
  /*fin mise a jour des indices de groupe1*/

  for(i=0;i<groupe2;i++)
    {
      if(i!=groupe1)
	{
	  for(j=i+1;j<groupe2;j++)
	    {
	      if(j!=groupe1)
		{
		  indicesWardNouveaux[i][j]=indicesWardAnciens[i][j];
		  indicesWardNouveaux[j][i]=indicesWardNouveaux[i][j];
		}
	    }
	}
    }
  
  for(i=0;i<groupe2;i++)
    {
      if(i!=groupe1)
	{
	  for(j=groupe2;j<nbGroupesCourants;j++)
	    {
	      indicesWardNouveaux[i][j]=indicesWardAnciens[i][j+1];
	      indicesWardNouveaux[j][i]=indicesWardNouveaux[i][j];
	    }
	}
    }
  
  for(i=groupe2;i<nbGroupesCourants;i++)
    {
      for(j=i+1;j<nbGroupesCourants;j++)
	{
	  indicesWardNouveaux[i][j]=indicesWardAnciens[i+1][j+1];
	  indicesWardNouveaux[j][i]=indicesWardNouveaux[i][j];
	}
    }

  /**********************************************/
  /*fin mise a jour des indices de Ward nouveaux*/
  /**********************************************/

  /*****************************************/
  /*mise a jour des indices de Ward anciens*/
  /*****************************************/

  for(i=0;i<nbGroupesCourants;i++)
    {
      for(j=0;j<nbGroupesCourants;j++)
	{
	  if(i!=j)
	    {
	      indicesWardAnciens[i][j]=indicesWardNouveaux[i][j];
	    }
	}
    }

  /*****************************************/
  /*mise a jour des indices de Ward anciens*/
  /*****************************************/
}

