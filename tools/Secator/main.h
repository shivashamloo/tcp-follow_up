#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <float.h>

/*constantes booleennes*/
#define OUI 1
#define NON 0
#define NOT_DEFINED -1
#define OK 1
#define PAS_OK 0

/*reussite d'une fonction*/
#define SUCCES 1
#define ECHEC 0
 
/*constantes du type de cluster*/
#define NORMAL 0
#define BRUIT 1
#define VIDE 2

/*constantes des types de donnees*/
#define COORDONNEES 0
#define ALIGNEMENT 1
#define DISTANCES 2
#define ARBRE 3

/*constantes des programmes de clustering*/
#define BIONJ 1
#define WARD 2

/*constantes de types de fichiers*/
#define FICHIER_TFA 0
#define FICHIER_MSF 1
#define FICHIER_INCONNU 2

/*constantes liees a la gestion de la memoire*/
#define TAILLE_NOM 200
#define TAILLE_POUBELLE 40000
#define TAILLE_MAX_LIGNE 40000


/*constantes diverses*/
#define AUCUN_GROUPE -1



/*types de prediction*/
#define CLASSIQUE 0
#define CICRV 1

typedef struct individu_t {
  int id;
  char *nom;
  double *valeurs;
  int cluster;
  char *description;
  int marque; /*0 si pas marque et 1 si marque*/
  int nbIndividusSimilaires;
  double valeurTemp;
  double *distancesPivots;
  struct noeud_t *noeud;}individu_t;

/*type d'un noeud de l'arbre phylogenetique*/
typedef struct noeud_t{
  int numero;
  int groupe;
  int qualite;
  int marquage;
  char etiquette[TAILLE_NOM];
  double poids;
  double dissimilarite;
  double distance;
  double distances[3]; /*en particulier, distances[0] est la distance au noeud superieur*/
  /*                     quand l'arbre est roote*/
  struct noeud_t *copain1;
  struct noeud_t *copain2;
  struct noeud_t *copains[3]; /*en particulier, copains[0] est le pere du noeud*/
  /*                            quand l'arbre est roote*/
  int secable;
  individu_t *individu;
}noeud_t;


typedef struct cluster_t {
  int nbIndividus;
  individu_t **individus;
  int insecable;
  struct cluster_t *precedent;
  struct cluster_t *suivant;
  struct cluster_t *fils1;
  struct cluster_t *fils2;
  struct cluster_t *pere;
  struct mixtureModel_t *mixtureModel;}cluster_t;

typedef struct couple {
  int i;
  int j;}couple_t;


/*****************************/
/*                           */
/*Procedure de message d'aide*/
/*                           */
/*****************************/

void helpMessageSecator();










