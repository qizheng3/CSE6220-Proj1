/*
 * CX 4220 / CSE 6220 Introduction to High Performance Computing
 *              Programming Assignment 1
 * 
 *  Utility function implementation
 * 
 */

#include <mpi.h>
#include "utils.h"
#include <math.h>

//------------------- Timer Functions (Do not change) -------------------//
void set_time(double &t, const int rank, MPI_Comm comm){
    if (rank>=0) // Do not call barrier if rank is negative
        MPI_Barrier(comm);
    if (rank <= 0){ //only 1 processor will set the time
        t = MPI_Wtime();
    }
}

double get_duration(double &t_start, double &t_end){
    return (t_end - t_start);
}
//---------------------------------------------------------------------//

/*********************************************************************
 *                 Implement your own functions here                 *
 *********************************************************************/
int getLocalN(int rank, int total_n, int p){
	int localN = 0;
	int maxN, minN;
	int FloorRank, CeilRank;
	maxN = ceil((total_n/(double)p));
	minN = maxN - 1;
	FloorRank = maxN*p - total_n;
	CeilRank = p - FloorRank;
	localN = (rank<CeilRank) ? maxN : minN;
	//std::cout<<"in getLocalN  p= "<<p<<"  maxN =  "<<maxN <<" minN= "<<minN<<std::endl;
	//std::cout<<"in getLocalN  n_local= "<<localN<<"  total_n = "<< total_n<<" rank= "<<rank<<"  ceilRank="<<CeilRank<<std::endl;
	return localN;
}





// ...