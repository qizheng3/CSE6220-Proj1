/*
 * CX 4220 / CSE 6220 Introduction to High Performance Computing
 *              Programming Assignment 1
 * 
 *  MPI polynomial evaluation algorithm function implementations go here
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include "mpi_evaluator.h"
#include "const.h"
#include "utils.h"

void scatter(const int n, double* scatter_values, int &n_local, double* &local_values, int source_rank, const MPI_Comm comm){
    //Implementation
	int p,rank;
	MPI_Comm_size(comm,&p);
	MPI_Comm_rank(comm,&rank);
	double total_n = 0;
	total_n = broadcast(n,source_rank,comm); //let all proc know the real n, in main, only rank0 has n
	n_local = getLocalN(rank,total_n,p);
	local_values = new double[n_local]; //allocate memory for local values
	MPI_Request send_req[p];
	MPI_Request recv_req;
	MPI_Barrier(comm);
	if (rank == source_rank){
		for (int j = 0; j<n_local;j++){
			local_values[j] = scatter_values[j];
		}
		int CurrentIndex =n_local;
		for (int i = 1; i<p; i++){
			int rank_p_localN = getLocalN(i,total_n,p);
			double* temp = new double[rank_p_localN];
			for (int j = 0; j<rank_p_localN;j++){
				temp[j] = scatter_values[CurrentIndex+j];
			}
			MPI_Isend(temp,rank_p_localN,MPI_DOUBLE,i, i+906,comm,&send_req[i]);
			MPI_Status myStatus;
			MPI_Wait(&send_req[i], &myStatus);
			CurrentIndex = CurrentIndex + rank_p_localN;
			delete(temp);
		}
		
	}else {
		MPI_Irecv(&local_values[0],n_local,MPI_DOUBLE,source_rank,rank+906,comm,&recv_req);
	}
	MPI_Barrier(comm);
	
}

double broadcast(double value, int source_rank, const MPI_Comm comm){
    //Implementation
	int p;
	int rank;
	MPI_Comm_size(comm,&p);
	MPI_Comm_rank(comm,&rank);
	
	double xValue = value;
	int gap = 1;
	int hasValue = 0;
	int hasValueTemp =0;
	double valueTemp;
	int pp;
	pp = pow(2,ceil(log2(p))); //get p' if p is not a power of 2
	MPI_Barrier(comm);
	if (rank == source_rank){
		hasValue = 1;
	}
	while(gap <=pp/2){
			if (rank/gap %2 == 0){
				if ((rank+gap)<p){
					MPI_Send(&xValue, 1, MPI_DOUBLE, rank+gap, 0, comm);
					MPI_Send(&hasValue, 1, MPI_INT, rank+gap, 1, comm);
					MPI_Status status;
					MPI_Recv(&valueTemp, 1, MPI_DOUBLE, rank+gap, 2, comm, &status);
					MPI_Recv(&hasValueTemp, 1, MPI_INT, rank+gap, 3, comm, &status);
				
					if (hasValueTemp ==1){
						hasValue = 1;
						xValue = valueTemp;
					}
				}
			}
			else{
					MPI_Send(&xValue, 1, MPI_DOUBLE, rank-gap, 2, comm);
					MPI_Send(&hasValue, 1, MPI_INT, rank-gap, 3, comm);
					MPI_Status status;
					MPI_Recv(&valueTemp, 1, MPI_DOUBLE, rank-gap, 0, comm, &status);
					MPI_Recv(&hasValueTemp, 1, MPI_INT, rank-gap, 1, comm, &status);
				if (hasValueTemp ==1){
					hasValue = 1;
					xValue = valueTemp;
				}
			}
		MPI_Barrier(comm);
		gap = 2* gap;
	}
    return xValue;
}

void parallel_prefix(const int n_local, const double* values, double* pre, const int OP, const MPI_Comm comm){
    //Implementation ===pre is for prefix_results===
	
	int p;
	int rank;
	MPI_Comm_size(comm,&p);
	MPI_Comm_rank(comm,&rank);
	int pp;
	pp = pow(2,ceil(log2((double)p))); //get p' if p is not a power of 2
	int gap = 1;
	int rank_p;
	for (int i=0;i<n_local;i++){
		pre[i] = values[i];
	}
	for (int i=1;i<n_local;i++){
		pre[i] = (OP == PREFIX_OP_PRODUCT)? (pre[i]*pre[i-1]) : (pre[i]+pre[i-1]);
	}
	double p_sum = pre[n_local-1]; // value for parallel prefix
	double t_sum = pre[n_local-1];
	double temp = 0;
	MPI_Barrier(comm);
	while(gap <=pp/2){
		rank_p = rank ^ gap; //use XOR
		if (rank_p<p) {
			MPI_Send(&t_sum, 1, MPI_DOUBLE, rank_p, OP+1+gap, comm);
			MPI_Status status;
			MPI_Recv(&temp, 1, MPI_DOUBLE, rank_p, OP+1+gap, comm, &status);
			t_sum = (OP == PREFIX_OP_PRODUCT)? (t_sum*temp):(t_sum+temp);
			if (rank_p < rank){
				p_sum = (OP == PREFIX_OP_PRODUCT)? (p_sum*temp):(p_sum+temp);
			}
		}
		MPI_Barrier(comm);
		gap = gap << 1; //use left shift
	}
	MPI_Barrier(comm);
	for (int i=0;i<n_local;i++){
		pre[i] = (OP == PREFIX_OP_PRODUCT)?(pre[i]*p_sum/pre[n_local-1]):(pre[i]+p_sum-pre[n_local-1]);
	}
	MPI_Barrier(comm);
}

double mpi_poly_evaluator(const double x, const int local_n, const double* constants, const MPI_Comm comm){
    //Implementation
	int p;
	int rank;
	MPI_Comm_size(comm,&p);
	MPI_Comm_rank(comm,&rank);
	
	double* prefix_product_results = new double[local_n];
	double* prefix_product_localV = new double[local_n];
	double* prefix_sum = new double[local_n];
	double prefix_total_sum = 0;
	if (rank == 0){
		for (int i=1; i<local_n;i++){
			prefix_product_localV[i] = x;
		}
		prefix_product_localV[0] = 1;
		
	} else
	{
		for (int i=0; i<local_n;i++){
			prefix_product_localV[i] = x;
		}
	}
	MPI_Barrier(comm);
	parallel_prefix(local_n,prefix_product_localV,prefix_product_results,PREFIX_OP_PRODUCT,comm);
	MPI_Barrier(comm);
	for(int i=0;i<local_n;i++){
		prefix_product_results[i]=constants[i]*prefix_product_results[i];
	}
	MPI_Barrier(comm);
	parallel_prefix(local_n,prefix_product_results,prefix_sum,PREFIX_OP_SUM,comm);
	
	MPI_Barrier(comm);
	
	if (rank == (p-1)) {
		MPI_Send(&prefix_sum[local_n-1], 1, MPI_DOUBLE, 0, 1111, comm);
		}
	MPI_Status status;
	if (rank == 0) MPI_Recv(&prefix_total_sum, 1, MPI_DOUBLE, p-1, 1111, comm, &status);
    MPI_Barrier(comm);
	delete(prefix_product_results);
	delete(prefix_product_localV);
	delete(prefix_sum);
	
	return prefix_total_sum;
}
