#ifndef FUNCTIONS_MIGRATION_H   /* Include guard */
#define FUNCTIONS_MIGRATION_H
unsigned long jcong;

int fulldichot (double *coord,double rand,int first,int last);

int findIndex(double distance, int nbins, double* breaks, double maxdist);

void makeDistClasses(double *xc, int *L, double *yc, int *cbin, int *indices, double *dists, int *nbbreaks, double *breaks, double *maxdist);

void makeDistClassesWithStreets(double *xc, int *L, double *yc, int *cbin, int *cbinas, int *cbinsb, int *indices, double *dists, int *nbbreaks, double *breaks, double *maxdist, int *blockIndex);

void generateProbMat(double* halfDistJ, double* halfDistH, int* useDelta, double* delta, double* rateHopInMove,double* rateSkipInMove, double* rateJumpInMove, double* dist_mat, double* prob_mat, int* blockIndex, int* cumul, int* L);

void modBinIt(int* n, int* dist_index, double* inf_data, double* start_inf_data, int* cbin, double* stats, int* nbins, int* endIndex);

void modBinItWithStreets(int* n, int* dist_index, double* inf_data, double* start_inf_data, int* cbin, int* cbinsb, int* cbinas, double* stats, int* nbins, int* blockIndex);

void gillespie(int *infested, int *endIndex, int *L, double *probMat, double *endTime, int *indexInfest, double *age, double *movePerTunit, int *seed);

void stratGillespie(int* infested,int* endIndex, int* L, double* rateHopInMove, double* rateSkipInMove, double* rateJumpInMove, int* hopColIndex, int* hopRowPointer, int* skipColIndex, int* skipRowPointer, int* jumpColIndex, int* jumpRowPointer, double* endTime, int* indexInfest, double* age, double* movePerTunit, double* introPerTunit, int* seed);

void polynomialfit(int obs, int degree, double *dx, double *dy, double *store);

int double_compare(const void *a, const void *b);

void get_stats_grid(int* rep, int* L, int* endInfest, int* endIndex, int* gridnbStats, int* numDiffGrids, int* gridIndexes, int* gridNumCells, int* gridEmptyCells, int* gridCountCells, double* gridstats);

void get_stats_circle(int* rep, int* L, int* endInfest, int* endIndex, int* circlenbStats, int* numDiffCircles, int* numDiffCenters, int* circleIndexes, int* circleCounts, double* circlestats);

void get_stats_semivar(int *rep, int *nbStats, int* L, int* dist_index, int* infestedInit, int* startInfested, int* cbin, int* cbinas, int* cbinsb, double* stats, int* nbins, int* blockIndex, int* haveBlocks, int* endIndex);

void get_stats_num_inf(int *rep, int *infnbstats, double* infstats, int* L, int* infestedInit, int* endIndex, int* blockIndex, int* haveBlocks);

void at_risk_stat(int *at_risk,int *n,int *posnodes, int *nPosnodes, double*dists, double *trs, int *nTr);

int follow_link(int *Grille, int sizeGrille, int num_case);

void percolation_circle(int *Nodes,int *n,double*dists,double *tr);

void simulObserved(int* L, int* infestedInit, int* endIndex, int* indexInfestInit, double* detectRate, int* seed);

void multiGilStat(double* probMat, int* useProbMat, double* distMat, double* halfDistJ, double* halfDistH, int* useDelta, double* delta, double* rateHopInMove, double* rateSkipInMove, double* rateJumpInMove, int* blockIndex, int *simul, int *infested, double *infestedDens, int *endIndex, int *L, double *endTime, int *indexInfest, double *age, double *scale, int *seed, int *Nrep, int* getStats, int *nbins, int *cbin, int* cbinas, int* cbinsb, int* indices, double* stats, int *nbStats);

void noKernelMultiGilStat(int* hopColIndex, int* hopRowPointer, int* skipColIndex, int* skipRowPointer, int* jumpColIndex, int* jumpRowPointer, double* rateHopInMove, double* rateSkipInMove, double* rateJumpInMove, int* blockIndex, int *simul, int *infested, double *infestedDens, int *endIndex, int *L, double *endTime, int *indexInfest, double *age, double *rateMove, double* rateIntro, int *seed, int *Nrep, int* getStats, int* matchStats, int* lengthStats, int *nbins, int *cbin, int* cbinas, int* cbinsb, int* indices, double* semivarstats, int *nbStats, int* haveBlocks, int* numDiffGrids, int* gridIndexes, int* gridNumCells, int* gridEmptyCells, int* gridCountCells, int* gridnbStats, double* gridstats, int* numDiffCircles, int* numDiffCenters, int* circleIndexes, int* circleCounts, int* circlenbStats, double* circlestats, int* infnbstats, double* infstats, double* xs, double* ys, double* detectRate);
#endif
