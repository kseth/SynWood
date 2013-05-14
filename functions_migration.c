/*
 * =====================================================================================
 *
 *       Filename:  functions_migration.c
 *
 *    Description:  routines for functions_migration.R
 *
 *        Version:  1.0
 *        Created:  09/28/2012
 *       Revision:  none
 *       Compiler:  R CMD SHLIB -lgsl -lgslcblas -lm functions_migration.c
 *       Nota: need gsl (ubuntu package libgsl0-dev)
 *
 * =====================================================================================
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_multifit.h>
// allows polynomial fitting
#include <R_ext/Utils.h>
// allows the use of R_CheckUserInterrupt()

#include "functions_migration.h"
/********************************/
//    Random number generator
/********************************/
unsigned long jcong=123124312;

/***** declare CONG*********/
// UNICONG is a draw of 0<=x<1
// equivalent to rand() in C but 4x as fast (save 26% comp. time)
#define TAUS88_MASK   0xffffffffUL /* required on 64 bit machines */
#define CONGMAX 4294967295.
#define CONG  (jcong=69069*jcong+1234567)
#define UNICONG ((CONG & TAUS88_MASK)*1/((double)CONGMAX+1.)) 
// +1 in the denominator guarantees <1
/***** end declaration CONG*****/

/********************************/
//    Binary search
/********************************/
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
	int i=0, maxit=last-first,val_sortie;
	int end=0;
	while(end==0){
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
		// to avoid infinite loop if wrong arguments
		maxit--;
		if(maxit<0){
			fprintf(stderr,"Infinite loop in fulldichot\n");
			exit(2);
		}

	}  
	return val_sortie;
}


// CB: this is a reimplementation of fulldichot from migration_model.c
//     should be unified but not a priority
// does a binary search to find the correct distance class for distance
// breaks must be in sorted order
// returns the dist class number or -1 if no corresponding distance class
int findIndex(double distance, int nbins, double* breaks, double maxdist){
	if(distance >= maxdist || distance < 0.0)
		return -1;

	int lo = 0;
	int hi = nbins - 1;
	int cur = (lo+hi)/2;
	
	if(distance < *(breaks + cur) && distance >= *(breaks + cur - 1))
		return (cur - 1);
	
	if(distance >= *(breaks + cur))
		lo = cur+1;
	else
		hi = cur;
	
	while(lo < hi)
	{
		cur = (lo+hi)/2;
		if(distance < *(breaks + cur) && distance >= *(breaks + cur - 1))
		return (cur - 1);
	
		if(distance >= *(breaks + cur))
			lo = cur+1;
		else
			hi = cur;	
	}

	return (lo-1);
}

// generate distances
void makeDistMat(double *xc // x of objects
		,int *L		// length of xc/yc
		,double *yc	// y of objects
		,double *dists  // eucl. distance for for each pair ij
		){

	double distance = 0.0;
	for(int i=0; i<*L; i++){
		for(int j=i+1; j<*L; j++){
			distance = hypot(xc[i] - xc[j], yc[i] - yc[j]);
			*(dists + i* *L+j) = distance;
			*(dists + j* *L+i) = distance;
			// printf("i %i j %i dist %.2f \n",i,j,distance);
		}
	}
}


// initial determination of the distances 
// and class indices for each pair
void makeDistClasses(double *xc // x of objects
		,int *L		// length of xc/yc
		,double *yc	// y of objects
		,int *cbin	// number of pairs per class
		,int *indices	// indice of class for each pair ij
		,double *dists  // eucl. distance for for each pair ij
		,int *nbbreaks	// nb breaks between classes
		,double *breaks // breaks between classes
		,double *maxdist// max distance considered for classes
				// should be max(breaks)
		){

	double distance = 0.0;
	int index = 0;	
	for(int i=0; i<*L; i++){
		for(int j=i+1; j<*L; j++){
			distance = hypot(xc[i] - xc[j], yc[i] - yc[j]);
			index = findIndex(distance, *nbbreaks, breaks, *maxdist);
			*(dists + i* *L+j) = distance;
			*(dists + j* *L+i) = distance;
			*(indices + i* *L+j) = index;
			
			if(index != -1)
				cbin[index]++;
			// printf("i %i j %i dist %.2f index %i cbin %i\n",i,j,distance,index,cbin[index]);
		} 
	}
}

// initial determination of the distances 
// and class indices for each pair
// implemented to return count of pairs within block and across streets (sb, as)
void makeDistClassesWithStreets(double *xc // x of objects
		,int *L		// length of xc/yc
		,double *yc	// y of objects
		,int *cbin	// number of pairs per class
		,int *cbinas	// number of pairs across streets per class
		,int *cbinsb	// number of pairs same block per class
		,int *indices	// indice of class for each pair ij
		,double *dists  // eucl. distance for for each pair ij
		,int *nbbreaks	// nb breaks between classes
		,double *breaks // breaks between classes
		,double *maxdist// max distance considered for classes
		,int *blockIndex// block of each house
		){

	double distance = 0.0;
	int index = 0;	
	for(int i=0; i<*L; i++){
		for(int j=i+1; j<*L; j++){
			distance = hypot(xc[i] - xc[j], yc[i] - yc[j]);
			index = findIndex(distance, *nbbreaks, breaks, *maxdist);
			*(dists + i* *L+j) = distance;
			*(dists + j* *L+i) = distance;
			*(indices + i* *L+j) = index;

			if(index != -1){
				cbin[index]++;
				if(*(blockIndex + i) == *(blockIndex + j))	
					cbinsb[index]++;
				else
					cbinas[index]++;
			}
			// printf("i %i j %i dist %.2f index %i cbin %i\n",i,j,distance,index,cbin[index]);
		} 
	}
}

// return prob_mat based on euclideant distances given in dist_mat
// cumul == 1, make cumulative probability matrix
// useDelta == 1, use the delta variable
// blockIndex - an array with block number for each house
void generateProbMat(double* halfDistJ,
			double* halfDistH, 
			int* useDelta, 
			double* delta, 
			double* rateHopInMove, 
			double* rateSkipInMove, 
			double* rateJumpInMove, 
			double* dist_mat,
			double* prob_mat,
			int* blockIndex,
			int* cumul,
			int* L){

	// Generating hop, skip, jump matrices
	double whop[*L];
	double wskip[*L];
	double wjump[*L];
	double weightHop = 0.0;
	double weightSkip = 0.0;
	double weightJump = 0.0;
	double totalWeight = 0.0;
	double distance = 0.0;
	double sumHop = 0.0;
	double sumSkip = 0.0;
	double sumJump = 0.0;
	
	// decreasing spatial link hop/skip
	for(int i=0; i<*L; i++){
		for(int j=0; j<*L; j++){	
			if(i==j){ // reaching end for this house
				//insects can't go from house to same house
				if(*cumul != 1){
					*(whop + j) = 0.0;
					*(wskip + j) = 0.0;
					*(wjump + j) = 0.0;
				}else{
					*(whop + j) = sumHop;
					*(wskip + j) = sumSkip;
					*(wjump + j) = sumJump;
				}

				continue;
			}
			
			distance = *(dist_mat + i* *L+j);

			if(*(blockIndex+i) == *(blockIndex + j)){
				weightHop = exp(-distance/ *halfDistH);
				weightSkip = 0.0;
			}else{
				weightSkip = exp(-distance/ *halfDistH);
				// in fact 
				if(*useDelta == 1){
					weightHop = weightSkip * *delta;
				}else{
					weightHop = 0.0;
				}
			}

			weightJump = exp(-distance/ *halfDistJ);
			
			sumHop += weightHop;
			sumSkip += weightSkip;
			sumJump += weightJump;
		
			// if not cumulative weight, each is just individual weight
			// if computing cumulative weight, each is sum weights	
			if(*cumul != 1){
				*(whop + j) = weightHop;
				*(wskip + j) = weightSkip;
				*(wjump + j) = weightJump;
			}else{
				*(whop + j) = sumHop;
				*(wskip + j) = sumSkip;
				*(wjump + j) = sumJump;
			}
		}
		
		//normalize each of the weights
		//add together to give overall prob_mat
		for(int j=0; j<*L; j++){
			weightHop = *(whop + j) / sumHop;
			weightSkip = *(wskip + j) / sumSkip;
			weightJump = *(wjump + j) / sumJump;
	
			totalWeight = weightHop * *rateHopInMove + weightSkip * *rateSkipInMove + weightJump * *rateJumpInMove;
			*(prob_mat + i* *L+j) = totalWeight;
		}
		
		sumHop = 0.0;
		sumSkip = 0.0;
		sumJump = 0.0;
	}
}

void modBinIt(int* n, int* dist_index, double* inf_data, double* start_inf_data, int* cbin, double* stats, int* nbins, int* endIndex, int* haveBlocks, int* blockIndex, int* cbinsb, int* cbinas){  
  // Calculate all the pair-wise statistics

	int ind=0;
	double prevalence = *endIndex + 1;
	prevalence /= *n;
	double sq_residual_prevalence = 0;
	double v=0.;
	double v_on=0.;
	double v_moran=0.;
	double v_geary=0.;
  	double *vbin = stats; //new new global semivar
  	double *sdbin = stats+ *nbins -1; //new new global sd
	double *vbin_on = sdbin + *nbins -1; //new old global semivar 
	double *sdbin_on = vbin_on + *nbins -1; //new old global sd 
	double *vbin_moran = sdbin_on + *nbins-1; //moran's I
	double *vbin_geary = vbin_moran + *nbins-1; //geary's C
	double *vbin_ripleyl = vbin_geary + *nbins-1; //ripley's L 
	double *vbin_sb_as; //sameblock - acrossstreets semivar
	double *sdbin_sb_as; //sameblock - acrossstreets sd
	double *vbinsb;
	double *vbinas;
	double *sdbinsb;
	double *sdbinas;

	if(*haveBlocks == 1){ //if have blocks, init sb-as machinery
		vbin_sb_as = vbin_ripleyl + *nbins-1;
		sdbin_sb_as = vbin_sb_as + *nbins-1; //sameblock - acrossstreets
		vbinsb = (double *) malloc(sizeof(double)* (*nbins-1));
		vbinas = (double *) malloc(sizeof(double)* (*nbins-1));
		sdbinsb = (double *) malloc(sizeof(double)* (*nbins-1));
		sdbinas = (double *) malloc(sizeof(double)* (*nbins-1));

		for (int class=0; class < (*nbins-1); class++){
			vbinsb[class] = 0;
			vbinas[class] = 0;
			sdbinsb[class] = 0;
			sdbinas[class] = 0;
		}
	}

	// this loop covers only unique i,j pairs
	for (int i=0; i< *n; i++){  // loop on all points
		for (int j=i+1; j<*n; j++){ //loop on half of the points

			// general variogram
			ind = *(dist_index + (i* *n) + j);

			if (ind != -1){
			        //(new-new) semivariance	
				v = inf_data[i] - inf_data[j];
				v = v*v;
				vbin[ind]+= v; 
				sdbin[ind] += v*v;
				
				if(*haveBlocks == 1){ //finding sb-as semivariance
					if(*(blockIndex + i) == *(blockIndex + j)){//same blocks
						vbinsb[ind] += v;
						sdbinsb[ind] += v*v;	
					}else{//across streets
						vbinas[ind] += v;
						sdbinas[ind] += v*v;	
					}
				}

				//(old-new) semivariance
				//need to go both ways (since i->j is not symmetric anymore)
				v_on = start_inf_data[i] - inf_data[j];
				v_on = v_on*v_on;
				vbin_on[ind]+=v_on;
				sdbin_on[ind]+=v_on*v_on;

				v_on = start_inf_data[j] - inf_data[i];
				v_on = v_on*v_on;
				vbin_on[ind]+=v_on;
				sdbin_on[ind]+=v_on*v_on;
			
				// Moran's I
				v_moran = (inf_data[i] - prevalence) * (inf_data[j] - prevalence);
				vbin_moran[ind] += v_moran;

				// Geary's C
				v_geary = inf_data[i] - inf_data[j];
				v_geary = v_geary*v_geary;
				vbin_geary[ind] += v_geary;

				// Ripley's K and L functions
				if(inf_data[i] == 1 && inf_data[j] == 1)
					vbin_ripleyl[ind] += 1;		
			}
		}

		sq_residual_prevalence += (inf_data[i] - prevalence) * (inf_data[i] - prevalence); //normalizing factor for geary + moran
	}

	// for each dist class finalize the computation of semi-variance
	for (int class=0; class < (*nbins-1); class++) 
	{

		if (cbin[class]>0){
			//remember cbin[class]=half number pairs in distance class 
			sdbin[class] = sqrt((sdbin[class] - ((vbin[class] * vbin[class])/cbin[class]))/(4*(cbin[class] - 1)));
			vbin[class] = vbin[class]/(2*cbin[class]); 
			sdbin_on[class] = sqrt((sdbin_on[class] - ((vbin_on[class] * vbin_on[class])/(2*cbin[class])))/(16*(cbin[class] - 1)));
			vbin_on[class] = vbin_on[class]/(4*cbin[class]);
			vbin_moran[class] = (vbin_moran[class] * *n)/(cbin[class] * sq_residual_prevalence);
			vbin_geary[class] = (vbin_geary[class] * (*n - 1))/(2*cbin[class] * sq_residual_prevalence);
			vbin_ripleyl[class] = sqrt(vbin_ripleyl[class]);
		} else {
			sdbin[class]=NAN;
			vbin[class]=NAN;
			sdbin_on[class]=NAN;
			vbin_on[class]=NAN;
			vbin_moran[class]=NAN;
			vbin_geary[class]=NAN;
			vbin_ripleyl[class]=NAN;
		}

		if(*haveBlocks == 1){ //calculate final sb-as statistic
			if(cbinas[class] > 0 && cbinsb[class] > 0){
				sdbinas[class] = sqrt((sdbinas[class] - ((vbinas[class] * vbinas[class])/cbinas[class]))/(4*(cbinas[class] - 1)));
				vbinas[class] = vbinas[class]/(2*cbinas[class]);

				sdbinsb[class] = sqrt((sdbinsb[class] - ((vbinsb[class] * vbinsb[class])/cbinsb[class]))/(4*(cbinsb[class] - 1)));
				vbinsb[class] = vbinsb[class]/(2*cbinsb[class]);

				vbin_sb_as[class] = vbinsb[class] - vbinas[class];
				sdbin_sb_as[class] = sqrt(sdbinsb[class]*sdbinsb[class] + sdbinas[class]*sdbinas[class]);
			}else{
				vbin_sb_as[class] = NAN;
				sdbin_sb_as[class] = NAN;
			}
		}


	}

	if(*haveBlocks == 1){ //free allocated vectors used only for blocks
		free(vbinsb);
		free(vbinas);
		free(sdbinsb);
		free(sdbinas);		
	}

}

// implemented to include streets and calculate semivariance same block and across streets
// has not been updated to handle old-new semivariance
void modBinItWithStreets(int* n, int* dist_index, double* inf_data, double* start_inf_data, int* cbin, int* cbinsb, int* cbinas, double* stats, int* nbins, int* blockIndex){  

	int ind=0;
	double v = 0.;
  	double *vbin = stats; // global semi-variance
  	double *sdbin = stats+ *nbins -1; // sd of the global semi-variance
	double *vbinsb = sdbin + *nbins - 1; // same block semi-variance
	double *sdbinsb = vbinsb + *nbins - 1; // sd same block 
	double *vbinas = sdbinsb + *nbins - 1; // accross streets semi
	double *sdbinas = vbinas + *nbins - 1; // sd a s
	double *vbin_s_d = sdbinas + *nbins - 1; // diff semi-var inter/intra block
	double *sdbin_s_d = vbin_s_d + *nbins - 1;

	// this loop covers only unique i,j pairs
	for (int i=0; i< *n; i++){  // loop on all points
		for (int j=i+1; j<*n; j++){ //loop on half of the points

			ind = *(dist_index + (i* *n) + j);

			if (ind != -1)
			{
				//general variogram   
				v = inf_data[i] - inf_data[j];
				v = v*v;
				vbin[ind]+= v; 
				sdbin[ind] += v*v;
				
				if(*(blockIndex + i) == *(blockIndex + j))
				{//same blocks variogram

					vbinsb[ind] += v;
					sdbinsb[ind] += v*v;	
				}
				else
				{//across streets variogram
					vbinas[ind] += v;
					sdbinas[ind] += v*v;	
				}
				
			}
		}
	}

	// for each dist class finalize the computation of semi-variance
	for (int class=0; class < (*nbins-1); class++) 
	{
		if (cbin[class]>0)
		{
			sdbin[class] = sqrt((sdbin[class] - ((vbin[class] * vbin[class])/cbin[class]))/(4*(cbin[class] - 1)));
			vbin[class] = vbin[class]/(2*cbin[class]);
			
		}
		else
		{
			sdbin[class]=NAN;
			vbin[class]=NAN;
		}

	    	// if as or sb pairs doesn't exist for a distance class
		// set the diff to NAN	
		if(cbinas[class] > 0 && cbinsb[class] > 0)
		{
			sdbinas[class] = sqrt((sdbinas[class] - ((vbinas[class] * vbinas[class])/cbinas[class]))/(4*(cbinas[class] - 1)));
			vbinas[class] = vbinas[class]/(2*cbinas[class]);

			sdbinsb[class] = sqrt((sdbinsb[class] - ((vbinsb[class] * vbinsb[class])/cbinsb[class]))/(4*(cbinsb[class] - 1)));
			vbinsb[class] = vbinsb[class]/(2*cbinsb[class]);


			vbin_s_d[class] = vbinsb[class] - vbinas[class];
			sdbin_s_d[class] = sqrt(sdbinsb[class]*sdbinsb[class] + sdbinas[class]*sdbinas[class]);
		}else{
			sdbinas[class]=NAN;
			vbinas[class]=NAN;
			sdbinsb[class]=NAN;
			vbinsb[class]=NAN;
			vbin_s_d[class]=NAN;
			sdbin_s_d[class]=NAN;
		}
	}
}

void gillespie(int *infested, int *endIndex, int *L, double *probMat, double *endTime, int *indexInfest, double *age, double *movePerTunit, int *seed){

	jcong = (unsigned long)*seed;
	// printf("seed: %li \n",jcong);

	double rand = UNICONG;
	
	//nextEvent - the time to the nextEvent
	double nextEvent = log(1-rand)/(-*movePerTunit * (*endIndex+1));

	//set starting time to zero
	double currentTime = 0;

	//the gillespie loop
	// printf("entering gillespie loop (endtime: %.4f)",*endTime);
	while(currentTime + nextEvent < *endTime && *endIndex+1 < *L){
		// printf("time %f Ninf %i ", currentTime, *endIndex+1);
		fflush(stdout);
		
		currentTime+=nextEvent;

		//pick a location to be infesting house
		rand = UNICONG;
		int index = (int)(rand* (*endIndex+1));
		int house = *(indexInfest + index);

		// printf("infesting: %i; ", house);

		//pick a new house to become infested from the infesting house
		rand = UNICONG;
		double *pinit=probMat+house* *L;
		// printf("h:%i, rand :%.4f in [%.4f;%.4f]; ",
		// 		house,rand,*pinit,*(pinit+*L-1));
		// fflush(stdout);
		
		int dest = fulldichot(pinit, rand, 0, *L-1);
		// ///printf("new infested: %i\n", dest);
		// ///fflush(stdout);
		
		if((*(infested+dest)) != 1){
			*endIndex+=1;
			*(infested+dest) = 1;
			*(indexInfest + *endIndex) = dest;
			*(age + *endIndex) = currentTime;
		}

		//calculate time to next event again
		rand = UNICONG;
		nextEvent = log(1-rand)/(-*movePerTunit * (*endIndex+1));
	}
	*seed=(int)jcong;
	// printf("final seed:%i",*seed);
}

void stratGillespie(int* infested,int* endIndex, int* L, double* rateHopInMove, double* rateSkipInMove, double* rateJumpInMove, int* hopColIndex, int* hopRowPointer, int* skipColIndex, int* skipRowPointer, int* jumpColIndex, int* jumpRowPointer, double* endTime, int* indexInfest, double* age, double* movePerTunit, double* introPerTunit, int* seed){
	
	jcong = (unsigned long)*seed;
	// printf("seed: %li \n",jcong);

	double rand = UNICONG;
	
	//nextEvent - the time to the nextEvent
	double nextEvent = log(1-rand)/(-*movePerTunit * (*endIndex+1) - *introPerTunit);
	double percentIntro = *movePerTunit * (*endIndex+1) + *introPerTunit;
	percentIntro = *introPerTunit/percentIntro;

	//set starting time to zero
	double currentTime = 0;
	
	//variables used in loop
	int index, house, dest, numHouses;

	//the gillespie loop
	// printf("entering gillespie loop (endtime: %.4f)",*endTime);
	while(currentTime + nextEvent < *endTime && *endIndex+1 < *L){
		// printf("time %f Ninf %i ", currentTime, *endIndex+1);
		// fflush(stdout);
		
		currentTime+=nextEvent;
		
		//pick whether new move or new introduction
		rand = UNICONG;
		if(rand < percentIntro){ // new introduction

			rand = UNICONG;
			dest = (int)(rand* *L);	//the introduction location
		}
		else{ // new move

			//pick a location to be infesting house
			rand = UNICONG;
			index = (int)(rand* (*endIndex+1));
			house = *(indexInfest + index);

			// printf("infesting: %i; ", house);
			
	
			//pick whether next move is hop/skip/jump
			
			dest = -1; //the move location
			rand = UNICONG;
			
			if(rand < *rateHopInMove){
				// next move is hop
				// pick new house to become infested
				// printf("hop ");
				rand = UNICONG;
				numHouses = hopRowPointer[house+1] - hopRowPointer[house];
				dest = hopColIndex[hopRowPointer[house] + (int)(rand*numHouses)]; 
			}
			else if(*rateHopInMove < rand && rand < (*rateHopInMove + *rateSkipInMove)){
				// next move is skip
				// pick new house to become infested
				// printf("skip ");
				rand = UNICONG;
				numHouses = skipRowPointer[house+1] - skipRowPointer[house];
				dest = skipColIndex[skipRowPointer[house] + (int)(rand*numHouses)]; 		
			}
			else{
				// next move is jump
				// pick new house to become infested
				// printf("jump ");
				rand = UNICONG;
				numHouses = jumpRowPointer[house+1] - jumpRowPointer[house];
				dest = jumpColIndex[jumpRowPointer[house] + (int)(rand*numHouses)]; 
			}
	
	
			// printf("new infested: %i\n", dest);
			// ///fflush(stdout);
		}			

		if((*(infested+dest)) != 1){ // if the destination is not already infested
			*endIndex+=1;
			*(infested+dest) = 1;
			*(indexInfest + *endIndex) = dest;
			*(age + *endIndex) = currentTime;
		}
	
		//calculate time to next event again
		rand = UNICONG;
		nextEvent = log(1-rand)/(-*movePerTunit * (*endIndex+1) - *introPerTunit);
		percentIntro = *movePerTunit * (*endIndex+1) + *introPerTunit;
		percentIntro = *introPerTunit/percentIntro;

	}

	*seed=(int)jcong;
	// printf("final seed:%i",*seed);
}

//taken from rosettacode.org
//obs, number of points to fit
//degree, degree of polynomial
//dx, array of x coords
//dy, array of y coords
//store, array of length degree to hold constants
void polynomialfit(int obs, int degree, double *dx, double *dy, double *store){
 
  gsl_multifit_linear_workspace *ws;
  gsl_matrix *cov, *X;
  gsl_vector *y, *c;
  double chisq;
 
  int i, j;
 
  X = gsl_matrix_alloc(obs, degree);
  y = gsl_vector_alloc(obs);
  c = gsl_vector_alloc(degree);
  cov = gsl_matrix_alloc(degree, degree);
 
  for(i=0; i < obs; i++) {
    gsl_matrix_set(X, i, 0, 1.0);
    for(j=0; j < degree; j++) {
      gsl_matrix_set(X, i, j, pow(dx[i], j));
    }
    gsl_vector_set(y, i, dy[i]);
  }
 
  ws = gsl_multifit_linear_alloc(obs, degree);
  gsl_multifit_linear(X, y, c, cov, &chisq, ws);
 
  /* store result ... */
  for(i=0; i < degree; i++)
  {
    store[i] = gsl_vector_get(c, i);
  }
 
  gsl_multifit_linear_free(ws);
  gsl_matrix_free(X);
  gsl_matrix_free(cov);
  gsl_vector_free(y);
  gsl_vector_free(c);
  return;
}

//comparator for qsort used in get_stats_grid
int double_compare(const void *a, const void *b){
	if(*(double*)a < *(double*)b) return -1;
	if(*(double*)a == *(double*)b) return 0;
	return 1;	
}

//have not implemented get_stats_grid for blocks
void get_stats_grid(int* rep, int* L, int* endInfest, int* endIndex, int* gridnbStats, int* numDiffGrids, int* gridIndexes, int* gridNumCells, int* gridEmptyCells, int* gridCountCells, int* numCoeffs, double* gridstats){

	//printf("%d \n", *endIndex);

	double* stats = gridstats + (*rep * *gridnbStats);
	int count = 0;
	//initialize all the storage cells (to store how many positive) in gridEmptyCells to 0
	for(int grid=0; grid<*numDiffGrids; grid++){

		for(int cell=0; cell<*(gridNumCells+grid);cell++){
			gridEmptyCells[count++] = 0;
		}
	}

	int currentIndexStartingPoint = 0;
	int currentCellStartingPoint = 0;
	int infestedCell = 0;

	//traverse through all the different grid systems
	for(int grid=0; grid<*numDiffGrids; grid++){

		currentIndexStartingPoint = *L * grid;

		//for each grid system		
		//traverse through all infested houses and populate
		//note that endIndex delineates the last spot that is occupied 
		for(int house=0; house<=*endIndex; house++){
			infestedCell = gridIndexes[currentIndexStartingPoint + endInfest[house]];
			gridEmptyCells[currentCellStartingPoint + infestedCell]++;

			//printf("%03d %03d ", infestedCell, gridEmptyCells[currentCellStartingPoint + infestedCell]);
		}

		//printf("\n");
		currentCellStartingPoint += *(gridNumCells+grid);
	}
  
	//the first stat inserted will be num positive cells
	//the second stat inserted will be variance of %positive 
	//the third-sixth stats will be regression coefficients
	// need to compute variance + store in gridstats
	// the percent positive is the mean %positive over all the cells (this is scale invariant)
	double meanPP = (*endIndex + 1) / *L;
	count = 0;

	// keeps track of number of positive cells in grid
	int positivecount = 0;

	//keeps track of the variance of the percent positive
	double varPP = 0;

	//num positive, num total per cell
	double numPositive = 0;
	double numTotal = 0;
	double cellPP= 0;

	//create *dx and *dy and *coeff (used in regression)
	double *dx, *dy; // coordinats of the initial points in the stats
	double coeff[*numCoeffs]; // the coefficients estimated

	for(int grid=0; grid<*numDiffGrids; grid++){

		//allocate *dx and *dy
		dx = (double *) malloc(sizeof(double)* *(gridNumCells+grid));
		dy = (double *) malloc(sizeof(double)* *(gridNumCells+grid));

		for(int cell=0; cell<*(gridNumCells+grid); cell++){

			numPositive = gridEmptyCells[count];
			numTotal = gridCountCells[count];
			cellPP = numPositive/numTotal;
			dx[cell] = ((double)(cell+1))/((double)(*(gridNumCells+grid)+1));
			dy[cell] = cellPP;

			if(numPositive > 0)
				positivecount++;
			
			varPP += cellPP*cellPP;
			count++;
		}


		//divide variance by number of cells and then subtract (mean percent positive)^2
		varPP = varPP/ *(gridNumCells+grid) - meanPP*meanPP;

		//calculate the regression coefficients
		//sort so that we have a "quantile like distribution"
		// actually making a "weibull plot"
		// dx: the number of the cell in the partition  
		// dy: the number of positive per cell
		qsort(dy, *(gridNumCells+grid), sizeof(double), double_compare);		
		polynomialfit(*(gridNumCells+grid), *numCoeffs, dx, dy, coeff);

		//store positive count in gridstats
		stats[grid* *gridnbStats] = positivecount; 
		//store the variance of the percent positive
		stats[grid* *gridnbStats+1] = varPP;
		//store the regression coefficients
		for(int c=0; c<*numCoeffs; c++){
			stats[grid* *gridnbStats + 2+c] = coeff[c];
		}
	
		//printf("cells: %d coeffs: %f %f %f %f\n", *(gridNumCells+grid), coeff[0], coeff[1], coeff[2], coeff[3]);	
		positivecount = 0;
		varPP = 0;
		free(dx);
		free(dy);	
	}
}

void get_stats_circle(int* rep, int* L, int* endInfest, int* endIndex, int* circlenbStats, int* numDiffCircles, int* numDiffCenters, int* circleIndexes, int* circleCounts, double* circlestats){

	printf("%i\n", *endIndex);

	//store the stats in the right place
	double* stats = circlestats + (*rep * *circlenbStats);

	//store number of positives per circle data
	double numPP[*numDiffCenters][*numDiffCircles];

	//instantiate the doubles to 0
	for(int center=0; center<*numDiffCenters; center++)
		for(int circle=0; circle<*numDiffCircles; circle++)
			numPP[center][circle] = 0;

	//for every house in endInfest (infested house)
	//put in appropriate numPP per center
	int wherePut = 0;
	int whichHouseInfested = 0;
	for(int house = 0; house <= *endIndex; house++)
		for(int center=0; center<*numDiffCenters; center++)
		{
			whichHouseInfested = *(endInfest+house);
			wherePut = circleIndexes[(center* *L)+whichHouseInfested];

			printf("%04d %04d\n", whichHouseInfested, wherePut); 
			if(wherePut != -1)
				numPP[center][wherePut] = numPP[center][wherePut]+1;
		}

	double count = 0;
	double meanPP = 0;
	double varPP = 0;
	for(int circle = 0; circle<*numDiffCircles; circle++)
	{	
	
		for(int center = 0; center<*numDiffCenters; center++)
		{
			count = circleCounts[(center* *numDiffCircles) + circle];
			meanPP = meanPP + numPP[center][circle]/count;
			varPP = varPP + (numPP[center][circle]/count) * (numPP[center][circle]/count);
		}

		meanPP = meanPP / *numDiffCenters;
		varPP = varPP/ *numDiffCenters - meanPP*meanPP;
		
		stats[circle*2] = meanPP;
		stats[circle*2+1] = varPP;
		printf("%f %f \n", meanPP, varPP);
		meanPP = 0;
		varPP = 0;
	}

	for(int num = 0; num < *circlenbStats; num++)
		printf("%f ", stats[num]);	
	
} 

void get_stats_semivar(int *rep, int *nbStats, int* L, int* dist_index, int* infestedInit, int* startInfested, int* cbin, int* cbinas, int* cbinsb, double* stats, int* nbins, int* blockIndex, int* haveBlocks, int* endIndex){  
	
	// cast infestedInit, startInfested from integer to double
  	double semivarianceData[*L];
	double startinfestData[*L];
	for(int h=0;h<*L;h++){
		semivarianceData[h] = infestedInit[h];   		
		startinfestData[h] = startInfested[h];
	}

	int startGVar=*rep* *nbStats;
	modBinIt(L, dist_index, semivarianceData, startinfestData, cbin, (stats+startGVar), nbins, endIndex, haveBlocks, blockIndex, cbinsb, cbinas);

	// if block data is passed (old blocks treatment)
	//	if(*haveBlocks == 1){
	//		// calculate semi-variance stats
	//		modBinItWithStreets(L, dist_index, semivarianceData, startinfestData, cbin, cbinsb, cbinas, (stats+startGVar), nbins, blockIndex); 
	//	
	//	}
}

void get_stats_num_inf(int *rep, int *infnbstats, double* infstats, int* L, int* infestedInit, int* endIndex, int* blockIndex, int* haveBlocks){

	//store the stats in the right place
	double* stats = infstats + (*rep * *infnbstats);	
	*(stats + 0) = *endIndex + 1;

	if(*haveBlocks == 1){	
		//now calculate 2 more stats:
		//number infested blocks, number infested houses/number infested blocks
			
		int currentBlock = 0;
		int maxBlock = -1;
		int minBlock = *L + 1;
	
	    	int infBlockCount = 0;
		for(int spot = 0; spot < *L; spot++){
			if(infestedInit[spot] == 1){ 
	           	 // to count the infested blocks
	            	 // find the maximum number of infested blocks
				currentBlock = blockIndex[spot];
				if(currentBlock > maxBlock)
					maxBlock = currentBlock;
				if(currentBlock < minBlock)
					minBlock = currentBlock;
			}
	    	}
		
		int length = (maxBlock - minBlock) + 1;
		int infBlocks[length];
			
		// initialize infestation of blocks
		for(int spot = 0; spot < length; spot++)
			infBlocks[spot] = 0;
			
	    	// detect infestation for each block
		for(int spot = 0; spot < *L; spot++)
			if(infestedInit[spot] == 1)
				infBlocks[blockIndex[spot] - minBlock]++;
			
		// count number of infested blocks
		for(int spot = 0; spot < length; spot++)
			if(infBlocks[spot] > 0)
				infBlockCount++;
		
		// save the corresponding stats
		*(stats + 1) = (double)infBlockCount;
		*(stats + 2) = ((double)(*endIndex + 1))/((double)infBlockCount);
	}	
}

// use at_risk_stat and polynomial fitting to return the at_risk stats
void get_stats_at_risk(int* rep, int* L, int* pos, int* npos, double* dists, double* trs, int* ntrs, double* at_risk_stat,int* ncoefs){
	int sizeLine = *ncoefs + *ntrs; // line in at_risk_stat
	double* at_risk_current = at_risk_stat+ sizeLine * *rep;
	// compute raw stats
	get_at_risk_stat(at_risk_current,L,pos,npos,dists,trs,ntrs);

	// fit of the stats with polynomial regression
	double coefs[*ncoefs];
	polynomialfit(*ntrs, *ncoefs, trs, at_risk_current, at_risk_current+ *ntrs);
}

//=======================================
// At risk stat
// - according to dist give the number of points 
//   within dist of infested
//=======================================

void get_at_risk_stat(double* at_risk_stats,int *L,int *posnodes, int *nPosnodes, double*dists,double *trs,int *nTr){

	// get matrix of indicator of at risk for each household
	// and each threshold
	int* at_risk_ind;
	at_risk_ind = (int *) malloc(sizeof(int)* (*nTr * *L));
	get_at_risk_indicator(at_risk_ind,L,posnodes,nPosnodes,dists,trs,nTr);

	// transform matrix of indicator in raw stat
	for ( int ih = 0; ih < *L; ih += 1 ) { 
		for(int itr=0;itr< *nTr;itr++){
			at_risk_stats[itr] += at_risk_ind[itr + ih * *nTr ];
		}
	}
	free(at_risk_ind);
}

//=======================================
// in at_risk (n*nTr) return for each tr and each node if at risk
//=======================================
void get_at_risk_indicator(int *at_risk,int *n,int *posnodes, int *nPosnodes, double*dists,double *trs,int *nTr){
  // the distances must be increasing in trs
  
  int summary[*nTr];
  for(int itr=0;itr< *nTr;itr++) summary[itr]=0;
  for(int node =0;node<*n;node++){ 
    int lnode = node * *nTr; // line for the node with != tr
    int atRisk = 0; // indicator that at risk with at least one tr
    for(int itr=0;itr< *nTr;itr++){ 
      int ipos = 0;
      while(ipos<*nPosnodes && atRisk ==0){ 
	int neigh = posnodes[ipos];
	if(dists[node* *n + neigh]<trs[itr]){ 
	  atRisk = 1; // if at risk with smal tr, at risk with bigger tr
	}
	ipos++;
      }
      at_risk[lnode+itr]=atRisk;
      summary[itr]+=atRisk;
    }
  }
  // printf("summary in indicator\n");
  // for(int itr=0;itr< *nTr;itr++) printf("%d ",summary[itr]);
  // printf("\n");
}


//=======================================
// Functions for percolation stat
//=======================================
// int remonte_liens2(int *Grille,int num_case){
int follow_link(int *Grille,int sizeGrille, int num_case){
	// sizeGrille, for security of the while loop: stop if more than sizeGrille iterations
	// remonte les liens de Grille
	// et retourne le numéro de la case initiale 
	// printf("num_case_init %i Grille[num_case] %i \n",num_case,Grille[num_case]);
	int nbit=0;
	while (Grille[num_case]!=num_case || nbit>sizeGrille){
		num_case=Grille[num_case];// passe a la case suivante
		// printf("num_case %i\n",num_case);
		nbit++;
	}
	// printf("num_case_final %i\n",num_case);

	return num_case;
}

void percolation_circle(int *Nodes,int *n,double*dists,double *tr){
	// Nodes has n nodes 0 (should select only the "1" before
	// dists: distance between the nodes
	// group nodes of Nodes==1 within tr (threshold)
	// and corresponding 0 within tr
	// return Nodes with 0: not in a group
	// other integers correspond to cliques
	// printf("entering_perc_circles\n");
	
	// different integer in all non 0 node
	int toUpdate[*n]; // vector with neighbors to update
	int minFlag=0;
	int itu=0;
	for ( int z = 0; z < *n; z += 1 ) { 
	  *(Nodes+z)=z;
	}
	// printf("debut_percolation for %d items\n",*n);
	// create the links
	for (int i = 0; i < *n; i += 1 ){ 
	  // printf("item: %d\n",i);
	  minFlag = Nodes[i];
	  itu = 0;
	  for(int di = 0; di < i; di+=1){
	    // printf("i %d,n %d;",i,di);
	    int neigh = i* *n+di;
	    if(dists[neigh]<*tr){
	      int flag = follow_link(Nodes,*n,di);
	      // printf("flag found: %d\n",flag);
	      if(flag<minFlag){
	         minFlag= flag;
	      }
	      toUpdate[itu] = flag; // add neigh to correct
	      itu++;
	    }
	  }

	  // printf("found: %d neighbors minimum flag %d\n",itu,minFlag);
	  toUpdate[itu] = -1; // set stop flag
	  for(int k = 0; k<itu;k++){
	    Nodes[toUpdate[k]] = minFlag;
	  }
	  Nodes[i] = minFlag;
	}
	// printf("calcul perco ok\n");
	// affiche_tableau(Nodes,h,l);
	// affectation de tous les liens
	for (int k = 0; k < *n; k += 1 ) { 
	    Nodes[k]=follow_link(Nodes,*n,k);
	}
	// printf("percolateur filled\n");

	// affiche_tableau(Percolateur,h,l);
}

void simulObserved(int* L, int* infestedInit, int* endIndex, int* indexInfestInit, double* detectRate, int* seed){

	if(*detectRate<1){ //if we have some error

		// printf("starting with %d inf ", *endIndex + 1);

		double randForHouse = 0;
		int tempHold[*endIndex+1];
		int tempCount = 0;

		srand(*seed); //initialize random number generator in stdlib

		for(int house=0; house<*L; house++){
		randForHouse = rand();
		randForHouse /= RAND_MAX; //get random value 0<=x<1

			if(randForHouse > *detectRate){ //remove house from infestedInit 
				infestedInit[house] = 0; //set house to 0
			}	
		}

		//now need to reorganize indexInfestInit
		for(int house=0; house<=*endIndex; house++){
			if(infestedInit[indexInfestInit[house]] == 1)
				tempHold[tempCount++] = indexInfestInit[house];
		}

		tempCount--;

		for(int house=0; house<=*endIndex; house++){
			if(house<=tempCount)
				indexInfestInit[house] = tempHold[house];
			else
				indexInfestInit[house] = -1;
			
		}

		*endIndex = tempCount;

		// printf("ending with %d inf\n", *endIndex+1);
	}

}

void multiGilStat(double* probMat, int* useProbMat, double* distMat, double* halfDistJ, double* halfDistH, int* useDelta, double* delta, double* rateHopInMove, double* rateSkipInMove, double* rateJumpInMove, int* blockIndex, int *simul, int *infested, double *infestedDens, int *endIndex, int *L, double *endTime, int *indexInfest, double *age, double *scale, int *seed, int *Nrep, int* getStats, int *nbins, int *cbin, int* cbinas, int* cbinsb, int* indices, double* stats, int *nbStats){

	// to pass to get_stats_semivar
	int haveBlocks = 1;

	// if not useProbMat, then generate probMat
	if(*useProbMat != 1){
		int cumul = 1;
		if(*distMat==-1 || *halfDistJ==-1 || *halfDistH ==-1 || *useDelta==-1 || *delta==-1 ||  *rateHopInMove==-1 || *rateSkipInMove==-1 ||  *rateJumpInMove==-1){
			printf("ERROR IN THETA VARIABLES, not all values provided, multiGilStat not executed\n");
			return;
		}
  
		generateProbMat(halfDistJ, halfDistH, useDelta, delta, rateHopInMove, rateSkipInMove, rateJumpInMove, distMat, probMat, blockIndex, &cumul, L);
	}
	
	int valEndIndex = *endIndex;	

	int infestedInit[*L];
  	int indexInfestInit[*L];
 
	for(int rep=0; rep< *Nrep; rep++){ // loop simul/stat
		R_CheckUserInterrupt(); // allow killing from R with Ctrl+c

	 	// initialisation simul
	 	for(int h=0;h<*L;h++){
	 		infestedInit[h]=*(infested+h);
	 		indexInfestInit[h]=*(indexInfest+h);	
	 	}
	 	
	 	*endIndex=valEndIndex; 

	 	if(*simul==1){ // simulation normal 
	 		
	 		gillespie(infestedInit,endIndex,L,probMat,endTime,indexInfestInit,age,scale,seed);

	 		for(int h=0;h<*L;h++){
	 			infestedDens[h]+=infestedInit[h];
	 		}
	 		
	 	}

	 	if(*getStats==1){
	 		get_stats_semivar(&rep, nbStats, L, indices, infestedInit, infested, cbin, cbinas, cbinsb, stats, nbins, blockIndex, &haveBlocks, endIndex);
	 	}

	 	if(*simul==0){ // no simulations, just stats 
	 		break; // to exit loop even if Nrep!=1
	 	}

	}
	
}

void noKernelMultiGilStat(
	int* hopColIndex, int* hopRowPointer, 
	int* skipColIndex, int* skipRowPointer, 
	int* jumpColIndex, int* jumpRowPointer, 
	double* rateHopInMove, double* rateSkipInMove, double* rateJumpInMove, 
	int* blockIndex, 
	int *simul, 
	int *infested, 
	double *infestedDens, 
	int *endIndex, 
	int *L, 
	double *endTime, 
	int *indexInfest, 
	double *age, 
	double *rateMove, double* rateIntro, 
	int *seed, int *Nrep, 
	int* getStats, 
	int* matchStats, int* lengthStats, 
	int *nbins, int *cbin, int* cbinas, int* cbinsb, int* indices, double* semivarstats, int *nbStats,
	int* haveBlocks, 
	int* numDiffGrids, int* gridIndexes, int* gridNumCells, int* gridEmptyCells, int* gridCountCells, int* gridnbStats, int* numCoeffs, double* gridstats, 
	int* numDiffCircles, int* numDiffCenters, int* circleIndexes, int* circleCounts, int* circlenbStats, double* circlestats,
	int* infnbStats, double* infstats,	
	double* atRisk_trs, int* atRisk_ntrs, // thresholds area at Risk stat
	double* at_riskStats, // results area at Risk stat, size (*atRisk_ntrs+*ncoefsAtRisk) * (*Nrep)
	int* ncoefsAtRisk, // number of coefs in poly fit of at_risk
    	double* xs, // Xs of the nodes
    	double* ys, // Ys of the nodes
    	double* detectRate){

	// if no blocks but still pass a rate skip
	// passing rateskip = 0 will prevent gillespie from skipping
	if(*skipColIndex == 0 && *skipRowPointer==0 && *rateSkipInMove != 0){
		printf("no skips given but rateSkipInMove!=0");
		return;
	}

	
	int valEndIndex = *endIndex;	
	int infestedInit[*L];
  	int indexInfestInit[*L];

	// make the distances matrix
	double* dists = (double *) calloc(*L * *L, sizeof(double));  //calloc, or 0 allocate, dists
	if(dists == NULL){
		printf("cannot (c)allocate memory");
		return;
	}
	
	// makeDistMat(xs,L,ys,dists);
 	printf("3");

	for(int rep=0; rep< *Nrep; rep++){ // loop simul/stat
		R_CheckUserInterrupt(); // allow killing from R with Ctrl+c

	 	// initialisation simul
	 	for(int h=0;h<*L;h++){
	 		infestedInit[h]=*(infested+h);
	 		indexInfestInit[h]=*(indexInfest+h);	
	 	}
	 	
	 	*endIndex=valEndIndex; 

	 	if(*simul==1){ // run a normal simulation
	 		
	 		stratGillespie(infestedInit,endIndex,L,rateHopInMove,rateSkipInMove,rateJumpInMove,hopColIndex,hopRowPointer,skipColIndex,skipRowPointer,jumpColIndex,jumpRowPointer,endTime,indexInfestInit,age,rateMove, rateIntro, seed);

			simulObserved(L, infestedInit, endIndex, indexInfestInit, detectRate, seed); // withhold data post gillespie

	 		for(int h=0;h<*L;h++){
	 			infestedDens[h]+=infestedInit[h];
	 		}
	 		
	 	}

	 	if(*getStats==1){ // calculate the statistics

			if(*simul!=1){ // if no simulation, need to withhold data and calculated infestedDens here instead 

				simulObserved(L, infestedInit, endIndex, indexInfestInit, detectRate, seed);

		 		for(int h=0;h<*L;h++){
	 				infestedDens[h]+=infestedInit[h];
	 			}
			}

			// always calculate number infested stats (quick + easy)
			get_stats_num_inf(&rep, infnbStats, infstats, L, infestedInit, endIndex, blockIndex, haveBlocks);

			for(int stat=0;stat<*lengthStats; stat++){

				// for every stat that we want, switch case

				int npos[1];
				npos[0] = *endIndex + 1;

				switch(matchStats[stat]){
	 				case 1:	get_stats_semivar(&rep, nbStats, L, indices, infestedInit, infested, cbin, cbinas, cbinsb, semivarstats, nbins, blockIndex, haveBlocks, endIndex); break;
					case 2: get_stats_grid(&rep, L, indexInfestInit, endIndex, gridnbStats, numDiffGrids, gridIndexes, gridNumCells, gridEmptyCells, gridCountCells, numCoeffs, gridstats); break;
					case 3: get_stats_circle(&rep, L, indexInfestInit, endIndex, circlenbStats, numDiffCircles, numDiffCenters, circleIndexes, circleCounts, circlestats); break; 
					case 4: get_stats_at_risk(&rep, L, indexInfestInit, npos, dists, atRisk_trs, atRisk_ntrs,at_riskStats, ncoefsAtRisk); break; 
					default: printf("stat that isn't supported yet\n"); break;
				}
			}
	 	}

	 	if(*simul==0){ // no simulations, just stats 
	 		break; // to exit loop even if Nrep!=1
	 	}

	}

	free(dists); //free malloc'ed dists
	
}
