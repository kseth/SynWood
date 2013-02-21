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
// compute only cumulative matrix
// not normalized: will need to draw between 0 and max when drawing
void generateProbMatGrant(
			double* limitHopSkips, 
			// double* weightHop, // the ref, is 1
			double* weightSkip, 
			double* weightJump, 
			double* dist_mat,
			double* prob_mat,
			int* blockIndex,
			int* L,
			int* cumul
			){

	// Generating hop, skip, jump matrices
	double weight = 0.0;
	double baseWeight = 0.0;
	double distance = 0.0;
	
	// decreasing spatial link hop/skip
	
	for(int i=0; i<*L; i++){
		// set self cell to proba=0
		if(i>0 && *cumul){
			*(prob_mat + i* *L+i) = *(prob_mat + i* *L+i-1);
		}else{
			*(prob_mat + i* *L+i) = 0;
		}

		// set other cells
		for(int j=i+1; j<*L; j++){	
			
			distance = *(dist_mat + i* *L+j);

			// baseWeight=exp(distance/ *limitHopSkips);
			baseWeight=1;
			if(distance<*limitHopSkips){
				if(*(blockIndex+i) == *(blockIndex + j)){
					weight = baseWeight;
				}else{
					weight = baseWeight* *weightSkip;
				}
			}else{
				weight = baseWeight * *weightJump;
			}
			if(*cumul){
				*(prob_mat + i* *L+j) = weight+*(prob_mat + i* *L+j-1); // current line is directly filled up
				if(i>0){
					// transpose of current Line can be 
					// summed to left column only if 
					// not first column
					*(prob_mat + j* *L+i) = weight+*(prob_mat + j* *L+i-1); 
				}else{
					*(prob_mat + j* *L+i) = weight;
				}
			}else{
				*(prob_mat + j* *L+i) = weight;
				*(prob_mat + i* *L+j) = weight;
			}
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

void modBinIt(int* n, int* dist_index, double* inf_data, int* cbin, double* stats, int* nbins){  

	int ind=0;
	double v=0.;
  	double *vbin = stats;
  	double *sdbin = stats+ *nbins -1;

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
			vbin_s_d[class] = NAN;
			sdbin_s_d[class] = NAN;
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

void stratGillespie(int* infested,int* endIndex, int* L, double* rateHopInMove, double* rateSkipInMove, double* rateJumpInMove, int* hopColIndex, int* hopRowPointer, int* skipColIndex, int* skipRowPointer, int* jumpColIndex, int* jumpRowPointer, double* endTime, int* indexInfest, double* age, double* movePerTunit, int* seed){
	
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
		// fflush(stdout);
		
		currentTime+=nextEvent;

		//pick a location to be infesting house
		rand = UNICONG;
		int index = (int)(rand* (*endIndex+1));
		int house = *(indexInfest + index);

		// printf("infesting: %i; ", house);
		

		//pick whether next move is hop/skip/jump
		
		int dest = -1;
		rand = UNICONG;
		
		if(rand < *rateHopInMove){
			// next move is hop
			// pick new house to become infested
			// printf("hop ");
			rand = UNICONG;
			int numHouses = hopRowPointer[house+1] - hopRowPointer[house];
			dest = hopColIndex[hopRowPointer[house] + (int)(rand*numHouses)]; 
		}
		else if(*rateHopInMove < rand && rand < (*rateHopInMove + *rateSkipInMove)){
			// next move is skip
			// pick new house to become infested
			// printf("skip ");
			rand = UNICONG;
			int numHouses = skipRowPointer[house+1] - skipRowPointer[house];
			dest = skipColIndex[skipRowPointer[house] + (int)(rand*numHouses)]; 		
		}
		else{
			// next move is jump
			// pick new house to become infested
			// printf("jump ");
			rand = UNICONG;
			int numHouses = jumpRowPointer[house+1] - jumpRowPointer[house];
			dest = jumpColIndex[jumpRowPointer[house] + (int)(rand*numHouses)]; 
		}


		// printf("new infested: %i\n", dest);
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

//have not implemented get_stats_grid for blocks
void get_stats_grid(int* rep, int* L, int* indexInfestInit, int* endIndex, int* gridnbStats, int* numDiffGrids, int* gridIndexes, int* gridNumCells, int* gridEmptyCells, int* gridCountCells, double* gridstats){

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
			infestedCell = gridIndexes[currentIndexStartingPoint + indexInfestInit[house]];
			gridEmptyCells[currentCellStartingPoint + infestedCell]++;

			//printf("%03d %03d ", infestedCell, gridEmptyCells[currentCellStartingPoint + infestedCell]);
		}

		//printf("\n");
		currentCellStartingPoint += *(gridNumCells+grid);
	}
  
	//determine the number of positive cells per grid system
	//the first stat inserted will be num positive cells
	//the second stat inserted will be variance of %positive 

	// need to compute variance + store in gridstats
	// the percent positive is the mean %positive over all the cells (this is scale invariant)
	double meanPP = (*endIndex + 1) / *L;
	count = 0;

	// how many stats per different grid
	// subtract one because blocks not yet implemented
	// int numStatsPerSystem = (*gridnbStats-1)/ *numDiffGrids;
	
	// keeps track of number of positive cells in grid
	int positivecount = 0;

	//keeps track of the variance of the percent positive
	double varPP = 0;

	//num positive, num total per cell
	double numPositive = 0;
	double numTotal = 0;
	double cellPP= 0;

	for(int grid=0; grid<*numDiffGrids; grid++){

		for(int cell=0; cell<*(gridNumCells+grid); cell++){

			numPositive = gridEmptyCells[count];
			numTotal = gridCountCells[count];
			cellPP = numPositive/numTotal;

			if(numPositive > 0)
				positivecount++;
			
			varPP += cellPP*cellPP;
			count++;
		}


		//divide variance by number of cells and then subtract (mean percent positive)^2
		varPP = varPP/ *(gridNumCells+grid) - meanPP*meanPP;

		//printf("numcells %d poscount %d varpp %f\n", *(gridNumCells+grid), positivecount, varPP);

		//store positive count in gridstats
		stats[grid*2] = positivecount; 
		//store the variance of the percent positive
		stats[grid*2+1] = varPP;
		
		positivecount = 0;
		varPP = 0;	
	}


	// the last statistic is number of positive houses
	stats[*numDiffGrids*2] = *endIndex + 1;
}

void get_stats_semivar(int *rep, int *nbStats, int* L, int* dist_index, int* infestedInit, int* cbin, int* cbinas, int* cbinsb, int* sizeVvar, double* stats, int* nbins, int* blockIndex, int* haveBlocks){  
	
	// cast infestedInit from integer to double
  	double semivarianceData[*L];
	for(int h=0;h<*L;h++){
		semivarianceData[h] = infestedInit[h];   		
	}

	// if block data is passed (note, sizeVvar must be handled appropriately in R code)
	if(*haveBlocks == 1){
		// calculate semi-variance stats
		int startGVar=*rep* *nbStats;
		modBinItWithStreets(L, dist_index, semivarianceData, cbin, cbinsb, cbinas, (stats+startGVar), nbins, blockIndex); 
	
		// move position over by that many stats	
		startGVar += *sizeVvar;
	
		//now calculate 3 more stats:
		//number infested houses, number infested blocks, number infested houses/number infested blocks
		
		int infCount = 0;
		int currentBlock = 0;
		int maxBlock = -1;
		int minBlock = *L + 1;
	
	    	int infBlockCount = 0;
		for(int spot = 0; spot < *L; spot++){
			if(infestedInit[spot] == 1){
				infCount++; // infested houses
	            
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
		*(stats + startGVar) = (double)infCount;
		*(stats + startGVar + 1) = (double)infBlockCount;
		*(stats + startGVar + 2) = ((double)infCount)/((double)infBlockCount);
	}else{ // don't have blocks

		// calculate semi-variance stats
		int startGVar=*rep* *nbStats;
		modBinIt(L, dist_index, semivarianceData, cbin, (stats+startGVar), nbins); 
	
		// move position over by that many stats	
		startGVar += *sizeVvar;
	
		//now calculate 1 more stat:
		//number infested houses		
		int infCount = 0;
		for(int spot = 0; spot < *L; spot++){
			if(infestedInit[spot] == 1){
				infCount++; // infested houses
	           	}
		} 

		// save the corresponding stats
		*(stats + startGVar) = (double)infCount;
	}
}

void multiGilStat(double* probMat, int* useProbMat, double* distMat, double* halfDistJ, double* halfDistH, int* useDelta, double* delta, double* rateHopInMove, double* rateSkipInMove, double* rateJumpInMove, int* blockIndex, int *simul, int *infested, double *infestedDens, int *endIndex, int *L, double *endTime, int *indexInfest, double *age, double *scale, int *seed, int *Nrep, int* getStats, int *nbins, int *cbin, int* cbinas, int* cbinsb, int* indices, double* stats, int *nbStats,int *sizeVvar){

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
	 		get_stats_semivar(&rep, nbStats, L, indices, infestedInit, cbin, cbinas, cbinsb, sizeVvar, stats, nbins, blockIndex, &haveBlocks);
	 	}

	 	if(*simul==0){ // no simulations, just stats 
	 		break; // to exit loop even if Nrep!=1
	 	}

	}
	
}

void noKernelMultiGilStat(int* hopColIndex, int* hopRowPointer, int* skipColIndex, int* skipRowPointer, int* jumpColIndex, int* jumpRowPointer, double* rateHopInMove, double* rateSkipInMove, double* rateJumpInMove, int* blockIndex, int *simul, int *infested, double *infestedDens, int *endIndex, int *L, double *endTime, int *indexInfest, double *age, double *scale, int *seed, int *Nrep, int* getStats, int* matchStats, int* lengthStats, int *nbins, int *cbin, int* cbinas, int* cbinsb, int* indices, double* stats, int *nbStats, int *sizeVvar, int* haveBlocks, int* numDiffGrids, int* gridIndexes, int* gridNumCells, int* gridEmptyCells, int* gridCountCells, int* gridnbStats, double* gridstats){

	// if no blocks but still pass a rate skip
	// passing rateskip = 0 will prevent gillespie from skipping 
	if(*haveBlocks == 0 && *rateSkipInMove != 0){
		printf("no blocks but rateSkipInMove!=0");
		return;
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
	 		
	 		stratGillespie(infestedInit,endIndex,L,rateHopInMove,rateSkipInMove,rateJumpInMove,hopColIndex,hopRowPointer,skipColIndex,skipRowPointer,jumpColIndex,jumpRowPointer,endTime,indexInfestInit,age,scale,seed);

	 		for(int h=0;h<*L;h++){
	 			infestedDens[h]+=infestedInit[h];
	 		}
	 		
	 	}

	 	if(*getStats==1){

			for(int stat=0;stat<*lengthStats; stat++){

				// for every stat that we want, switch (if 1, do semivariance stats; if 2, do grid stats)	
				switch(matchStats[stat]){
	 				case 1:	get_stats_semivar(&rep, nbStats, L, indices, infestedInit, cbin, cbinas, cbinsb, sizeVvar, stats, nbins, blockIndex, haveBlocks); break;
					case 2: get_stats_grid(&rep, L, indexInfestInit, endIndex, gridnbStats, numDiffGrids, gridIndexes, gridNumCells, gridEmptyCells, gridCountCells, gridstats); break;
					default: printf("stat that isn't supported yet\n"); break;
				}
			}
	 	}

	 	if(*simul==0){ // no simulations, just stats 
	 		break; // to exit loop even if Nrep!=1
	 	}

	}
	
}
