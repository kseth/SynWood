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
 *       Compiler:  R CMD SHLIB functions_migration.c
 *
 * =====================================================================================
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R_ext/Utils.h>
// allows the use of R_CheckUserInterrupt()


/********************************/
//    Random number generator
/********************************/
unsigned long jcong=123124312;
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


void binit(int *n, 	   // number of points
		double *xc,  // coordinates
		double *yc, 
		double *sim, // the input data (values for each point) 
	   int *nbins, // nb breaks 
	   double *lims, 	   // breaks vector
	   double *maxdist,  // threshold of considered distances (usually max(lims))
	   // results
	   int *cbin,   // number of house pairs per distance class   
	   double *vbin,	// total semivariance per distance class
			    // sum of square diffs
	   double *sdbin // total standard deviation per distance class      
	   ){
  
  int i, j, ind=0;
  double v=0.0;
  double dist=0.0, dx=0.0, dy=0.0;
 
  // this loop covers only unique i,j pairs
  for (j=0; j < *n; j++){  // loop on all points
      for (i=j+1; i<*n; i++){ //loop on half of the points
	  // euclidian dist
	  dx = xc[i] - xc[j];
	  dy = yc[i] - yc[j];
	  dist = hypot(dx, dy);
	  
	  // if distance i,j < max allowed distance
	  // value (v) = (data at i - data at j)^2/2:
	  if(dist <= *maxdist){
	      v = sim[i] - sim[j];
		v = (v*v)/2.0;
		
	//start at first index
	      ind = 0;
	//loop through indices
	      while (dist >= lims[ind] && ind <= *nbins ) ind++ ;
	//find the minimum distance class where dist is acceptable
	//should be: ind < *nbins (and remove following if)
   
		//printf("index: %i dist: %.2f \n", ind, lims[ind]);	

		// vbin: the semi-variance (sum of the square diffs)
		// cbin: count by distance classes 
		// sdbin: will become the standard dev of dist class	
		if (dist < lims[ind])
		{
		  vbin[(ind-1)]+= v; 
		  cbin[(ind-1)]++;
		  sdbin[(ind-1)] += v*v;
		}
	     }
	}
    }
  
  // for each dist class  finalize the computation of semi-variance
  for (j=0; j < *nbins; j++){
      if (cbin[j]){
	sdbin[j] = sqrt((sdbin[j] - ((vbin[j] * vbin[j])/cbin[j]))/(cbin[j] - 1));
	vbin[j] = vbin[j]/cbin[j];
      }
    }

  // for(int i=0; i<*nbins; i++)
  //    printf("dist: %.2f pairs in bin: %i semivar: %.2f stdev: %.2f \n", lims[i], cbin[i], vbin[i], sdbin[i]);	

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
			// *(dists + j* *L+i) = distance;
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
			// *(dists + j* *L+i) = distance;
			*(indices + i* *L+j) = index;			

			if(index != -1)
			{	
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
	int sameblock = 0;

	
	// decreasing spatial link hop/skip
	for(int i=0; i<*L; i++)
	{
		for(int j=0; j<*L; j++)
		{	
			//insects can't go from house to same house
			if(i==j)
			{
				if(*cumul != 1)
				{
					*(whop + j) = 0.0;
					*(wskip + j) = 0.0;
					*(wjump + j) = 0.0;
				}
				else
				{
					*(whop + j) = sumHop;
					*(wskip + j) = sumSkip;
					*(wjump + j) = sumJump;
				}

				continue;
			}
			
			distance = *(dist_mat + i* *L+j);
			
			if(distance == 0)
				distance = *(dist_mat + j* *L+i);
	
			sameblock = 0;
			if(*(blockIndex+i) == *(blockIndex + j))
				sameblock = 1;
			
			if(sameblock == 1)
			{
				weightHop = exp(-distance/ *halfDistH);
				weightSkip = 0.0;
			}
			else
			{
				weightSkip = exp(-distance/ *halfDistH);
				weightHop = 0.0;
			}
			
			//if we want to use delta, than weighthop+=weightSkip * delta
			if(*useDelta == 1)
			{
				weightHop += weightSkip * *delta;
			}	
			
			weightJump = exp(-distance/ *halfDistJ);
			
			sumHop += weightHop;
			sumSkip += weightSkip;
			sumJump += weightJump;
		
			// if not cumulative weight, each is just individual weight
			// if computing cumulative weight, each is sum weights	
			if(*cumul != 1)
			{
				*(whop + j) = weightHop;
				*(wskip + j) = weightSkip;
				*(wjump + j) = weightJump;
			}
			else
			{
				*(whop + j) = sumHop;
				*(wskip + j) = sumSkip;
				*(wjump + j) = sumJump;
			}
		}
		
		//normalize each of the weights
		//add together to give overall prob_mat
		for(int j=0; j<*L; j++)
		{
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

void modBinIt(int* n, int* dist_index, double* inf_data, int* cbin, double* stats, int* nbins){  

	int ind=0;
	double v=0.;
  	double *vbin = stats;
  	double *sdbin = stats+ *nbins -1;

	// printf("vbin: %p, sdbin %p\n",vbin,sdbin);
	// printf("vbin[0]: %f, sdbin[0] %f\n",vbin[0],sdbin[0]);

	// this loop covers only unique i,j pairs
	for (int i=0; i< *n; i++){  // loop on all points
		for (int j=i+1; j<*n; j++){ //loop on half of the points

			// general variogram
			ind = *(dist_index + (i* *n) + j);

			if (ind != -1){   
				v = inf_data[i] - inf_data[j];
				v = v*v;
				vbin[ind]+= v; 
				sdbin[ind] += v*v;
			}
		}
	}
	// printf("2 vbin[0]: %f, sdbin[0] %f\n",vbin[0],sdbin[0]);

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
	}
	// printf("3 vbin[0]: %f, sdbin[0] %f\n",vbin[0],sdbin[0]);

	// // display results by class there are (*nbins-1) class
	// for(int i=0; i<(*nbins-1); i++){
	// 	printf("index: %i pairs in bin: %i semivar: %.4f stdev: %.4f \n", i, cbin[i], vbin[i], sdbin[i]);	
	// }
}

// implemented to include streets and calculate semivariance same block and across streets
void modBinItWithStreets(int* n, int* dist_index, double* inf_data, int* cbin, int* cbinsb, int* cbinas, double* stats, int* nbins, int* blockIndex){  

	int ind=0;
	double v = 0.;
  	double *vbin = stats; // global semi-variance
  	double *sdbin = stats+ *nbins -1; // sd of the global semi-variance
	double *vbinsb = sdbin + *nbins - 1; // same block semi-variance
	double *sdbinsb = vbinsb + *nbins - 1; // sd same block 
	double *vbinas = sdbinsb + *nbins - 1; // accross streets semi
	double *sdbinas = vbinas + *nbins - 1; // sd a s
	double *vbin_s_d = sdbinas + *nbins - 1; // diff semi-var inter/intra block

	// printf("vbin: %p, sdbin %p\n",vbin,sdbin);
	// printf("vbin[0]: %f, sdbin[0] %f\n",vbin[0],sdbin[0]);

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
	// printf("2 vbin[0]: %f, sdbin[0] %f\n",vbin[0],sdbin[0]);

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
	
		if (cbinas[class]>0)
		{
			sdbinas[class] = sqrt((sdbinas[class] - ((vbinas[class] * vbinas[class])/cbinas[class]))/(4*(cbinas[class] - 1)));
			vbinas[class] = vbinas[class]/(2*cbinas[class]);
			
		}
		else
		{
			sdbinas[class]=NAN;
			vbinas[class]=NAN;
		}

		if (cbinsb[class]>0)
		{
			sdbinsb[class] = sqrt((sdbinsb[class] - ((vbinsb[class] * vbinsb[class])/cbinsb[class]))/(4*(cbinsb[class] - 1)));
			vbinsb[class] = vbinsb[class]/(2*cbinsb[class]);
			
		}
		else
		{
			sdbinsb[class]=NAN;
			vbinsb[class]=NAN;
		}
		
		if(cbinas[class] > 0 && cbinsb[class] > 0)
			vbin_s_d[class] = vbinsb[class] - vbinas[class];
		else
			vbin_s_d[class] = NAN;
	
	}

	// printf("3 vbin[0]: %f, sdbin[0] %f\n",vbin[0],sdbin[0]);

	// // display results by class there are (*nbins-1) class
	// for(int i=0; i<(*nbins-1); i++){
	// 	printf("index: %i pairs in bin: %i semivar: %.4f stdev: %.4f \n", i, cbin[i], vbin[i], sdbin[i]);	
	// }
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
	while(currentTime + nextEvent < *endTime){
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

void get_stats(int *rep, int *nbStats, int* L, int* dist_index, int* infestedInit, int* cbin, int* cbinas, int* cbinsb, int* sizeVvar, double* stats, int* nbins, int* blockIndex){  
	
	// cast infestedInit from integer to double
  	double semivarianceData[*L];
	double nbInfestedHouses=0;
	for(int h=0;h<*L;h++){
		semivarianceData[h] = infestedInit[h];   		
		nbInfestedHouses += infestedInit[h];
	}

	// calculate semi-variance stats
	int startGVar=*rep* *nbStats;
	// modBinIt(L, dist_index, semivarianceData, cbin, (stats+startGVar), nbins);
	modBinItWithStreets(L, dist_index, semivarianceData, cbin, cbinsb, cbinas, (stats+startGVar), nbins, blockIndex); 

	// add the number of houses infested
	*(stats+startGVar+*sizeVvar)=nbInfestedHouses;
}

void multiGilStat(double* probMat,double* distMat, int* blockIndex, int *simul, int *infested, double *infestedDens, int *endIndex, int *L, double *endTime, int *indexInfest, double *age, double *scale, int* sizeScale, int *seed, int *Nrep, int* getStats, int *nbins, int *cbin, int* cbinas, int* cbinsb, int* indices,int* sizeVvar, double* stats, int *nbStats){

	int valEndIndex = *endIndex;	

	int infestedInit[*L];
  	int indexInfestInit[*L];
 
  	double vbinInit[*nbins];
  	double sdbinInit[*nbins];

	for(int numTheta=0; numTheta<*sizeScale;numTheta++){
		for(int rep=0; rep< *Nrep; rep++){ // loop simul/stat
			R_CheckUserInterrupt(); // allow killing from R with Ctrl+c

			// initialisation simul
			for(int h=0;h<*L;h++)
			{
				infestedInit[h]=*(infested+h);
				indexInfestInit[h]=*(indexInfest+h);	
			}

			*endIndex=valEndIndex; 

			// printf("Init simul(%i)\n",*simul);
			if(*simul==1){ // simulation normal 

				gillespie(infestedInit,endIndex,L,probMat,endTime,indexInfestInit,age,scale,seed);

				for(int h=0;h<*L;h++){
					infestedDens[h]+=infestedInit[h];
				}

				// printf("rep: %i endIndex:%i seed:%i \n",rep,*endIndex,*seed);
			}else{ // do stats on initial data
			}
			// printf("simul OK, *getStats: %i\n",*getStats);

			if(*getStats==1){
				//printf("getting stats");
				get_stats(&rep, nbStats, L, indices, infestedInit, cbin, cbinas, cbinsb, sizeVvar, stats, nbins, blockIndex);
			}

			if(*simul==0){ // simulation normal 
				break; // to exit loop even if Nrep!=1
			}

		}
		scale++;
	}
	
}


void multiGilStat_C_ProbMat(double* dist_mat, double* probMat, double* halfDistJ, double* halfDistH, int* useDelta, double* delta, double* rateHopInMove, double* rateSkipInMove, double* rateJumpInMove, int* blockIndex, int *simul, int *infested, double *infestedDens, int *endIndex, int *L, double *endTime, int *indexInfest, double *age, double *scale, int *seed, int *Nrep, int *nbins, int *cbin, int* cbinas, int* cbinsb,int* sizeVvar, int* indices, double* stats,int *nbStats){

	int valEndIndex = *endIndex;	

	int infestedInit[*L];
  	int indexInfestInit[*L];
 
  	double vbinInit[*nbins];
  	double sdbinInit[*nbins];

	int cumul = 1;
	generateProbMat(halfDistJ, halfDistH, useDelta, delta, rateHopInMove, rateSkipInMove, rateJumpInMove, dist_mat, probMat, blockIndex, &cumul, L);

	for(int rep=0; rep< *Nrep; rep++){ // loop simul/stat
		R_CheckUserInterrupt(); // allow killing from R with Ctrl+c

		// initialisation simul
		for(int h=0;h<*L;h++)
		{
			infestedInit[h]=*(infested+h);
			indexInfestInit[h]=*(indexInfest+h);	
		}
		
		*endIndex=valEndIndex; 

		// printf("Init simul(%i)\n",*simul);
		if(*simul==1)
		{ // simulation normal 
			
			gillespie(infestedInit,endIndex,L,probMat,endTime,indexInfestInit,age,scale,seed);

			for(int h=0;h<*L;h++)
			{
				infestedDens[h]+=infestedInit[h];
			}
			// printf("rep: %i endIndex:%i seed:%i \n",rep,*endIndex,*seed);
		}
		else
		{ // do stats on initial data
		}
		// printf("simul OK\n");

		get_stats(&rep, nbStats, L, indices, infestedInit, cbin, cbinas, cbinsb, sizeVvar, stats, nbins, blockIndex);

		if(*simul==0)
		{ // simulation normal 
			break; // to exit loop even if Nrep!=1
		}

	}
	
}


