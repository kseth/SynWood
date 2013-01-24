/*
 * =====================================================================================
 *
 *       Filename:  migration_model.c
 *
 *    Description:  routines for migration_model.r
 *    			R CMD SHLIB migration_model.c
 *
 *        Version:  1.0
 *        Created:  06/09/2012 05:53:47 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  C. Barbu (CB), corentin.barbu@gmail.com
 *        Company:  CCEB, University of Pennsylvania
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

unsigned long jcong=123124312312;
/***** déclaration CONG********	*/
// UNICONG is a draw of 0<= x < 1
// équivalent à ran() de C mais en 4 fois plus rapide
// sur le code ici, on gagne 26% de temps de calcul
#define TAUS88_MASK   0xffffffffUL /* required on 64 bit machines */
#define CONGMAX 4294967295.
#define CONG  (jcong=69069*jcong+1234567)
#define UNICONG ((CONG & TAUS88_MASK)*1/((double)CONGMAX+1.)) 
// le plus +1 permet d'éviter d'avoir 1 ce qui serait embétant dans le pg
/***** fin déclaration CONG****	*/

int fulldichot (double *coord,double rand,int first,int last){
// la sortie est le numéro de la case contenant 
// le plus petit réel plus grand ou égal à la référence rand
// dans un vecteur comprennant des réels croissants.
// first et last permettent de donner éventuellement les numéros des 
// cases entre lesquelles on doit choisir
//
// utile notamment lorsqu'on a un choix aléatoire, dans ce cas :
// coord est un vecteur de probabilités cumulées ou plus simplement des poids cumulés
// rand la valeur aléatoire tirée ou plus simplement (rapidement) la valeur fois la dernière case des poidscumulés
// first = 0
// last = nb_cases-1
// test de fulldichot :
	int i=last, lasti=0,val_sortie;
	int end=0;
	while(end==0){
		lasti=i;
		i=((last-first)/2)+first;
		// printf("rand=%f coord[%d]=%f\n",rand,i,coord[i]);

		if(coord[i]>=rand){
			if(i==0 || coord[i-1]<rand){
				end=1;
				val_sortie = i;
			}
			else last=i;
		}
		else if(coord[i]<rand){
			if(i==last-1 || coord[i+1]>=rand){
				end=1;
				val_sortie =  i+1;
			}		
			else first=i;
		}
		fprintf(stderr,"lasti %i i %i",lasti,i);
		if(lasti==i){
			fprintf(stderr,"erreur in fulldichot\n");
			return(-1);
		}

	}  
	return val_sortie;
}


// dispatch vector from moving to arrival according
// to lines of probMat
// L: nb of houses (size of moving)
// probMat: L by L matrix, lines of cummulative probability, last being 1
// arrival: vector to be incremented
void basicMigration(int *moving, int *L,double *probMat,int *arrival){
	int orig=0;
	int nbitem=0;
	for ( orig = 0; orig < *L; orig += 1 ) { // for all initial orig
		// // check that probMat line end up to 1
		// printf("i %i endLine %f \n",i,probMat[(i+1)* *L-1]);

		// successively move each item
		for ( nbitem = 0; nbitem < moving[orig]; nbitem += 1 ) { 
			// draw a number between 0 and 1
			double x = UNICONG;
			int dest = fulldichot(probMat+orig* *L,x,0,*L-1);
			// printf("x %f (%i)->(%i) old %i",x,orig,dest,arrival[dest]);
			arrival[dest]+=1;
			// printf("new %i \n",arrival[dest]);
		}
		
	}
}

// dispatch vector from moving to arrival according
// to lines of probMat
// L: nb of houses (size of moving)
// probMat: L by L spam matrix, lines of cummulative probability, last being 1
// arrival: vector to be incremented
// kind of broken, infinite loop if enough individuals to migrate
// left there as no more a priority
void basicMigrationSpam(int *moving, int *L,double *PMentries,int *PMcolindices,int *PMrowpoint,int *arrival){
	int orig=0;
	int nbitem=0;
	for ( orig = 0; orig < *L; orig += 1 ) { // for all initial orig
		// // check that probMat line end up to 1
		// double sum=0;
		// int entry=0;
		// for ( entry = PMrowpoint[orig]; entry < PMrowpoint[orig+1]; entry += 1 ) { 
		// 	sum +=PMentries[entry-1]; 
		// }
		// // if(abs(sum-1)>10e-3){
		// printf("orig %i sum %f (%i,%i)\n",orig,sum,PMrowpoint[orig],PMrowpoint[orig+1]);
		// // }

		// successively move each item
		for ( nbitem = 0; nbitem < moving[orig]; nbitem += 1 ) { 
			// draw a number between 0 and 1
			double x = UNICONG;

			// get the corresponding in the entries of ProbMat
			int PMdest = fulldichot(PMentries,x,PMrowpoint[orig]-1,PMrowpoint[orig+1]-2);
			// get the corresponding in arrival
			int dest = PMcolindices[PMdest]-1;
			printf("n %i x %f (%i)->(%i)->(%i) old %i",nbitem,x,orig,PMdest,dest,arrival[dest]);
			arrival[dest]+=1;
			printf(" new %i p %.3f \n",arrival[dest],*(PMentries+PMdest));
		}
		printf("L %i",*L);
		
	}
}


