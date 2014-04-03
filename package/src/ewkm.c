// Entropy Weighted K-Means

// Original C code for k-means by Joshua Huang. Qiang Wang implemented
// the entropy weighted algorithm. Xiaojun prepared the code for
// packaging, and Graham Williams finalised the code for release and
// maintains the package.

// Copyright (c) 2011 Shenzhen Institutes of Advanced Technology
// Chinese Academy of Sciences

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <R.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>

#include <Rmath.h>
#include "sum_squares.h"
double unif_rand(void);

// ***** Support functions *****

// ----- Create the initial prototype (k centroids) -----

void initPrototypes(// Inputs ---------------------------------------------
		    double *x,	// Numeric matrix as vector by col (nr*nc) 
		    int *nr, 	// Number of rows
		    int *nc, 	// Number of columns
		    int *k,	// Number of clusters
		    // Output ---------------------------------------------
		    double **o_prototype) // Numeric prototype matrix (k*nc)
{
  int i, j, l;
  int flag = 0;
  int index;

  int *random_obj_num; 	// Array for randomly selected objects (k)

  // Memory for array of randomly selected objects

  random_obj_num = (int *) malloc(sizeof(int) * (*k));
  if (!random_obj_num) 
  {
    error("Can't allocate memory for random_obj_num matrix\n");
  }

  for (l = 0; l < *k; l++) random_obj_num[l] = -1;

  // Randomly select k objects.
  
  for (l = 0; l < *k; l++) 
  {
    flag = 1;

    while (flag) 
    {
      index = (int) (rand() % (*nr));
      flag = 0;
      for (i = 0; i < l; i++)
	if (random_obj_num[i] == index)
	  flag = 1;
    }

    random_obj_num[l] = index;
    for (j = 0; j < (*nc); j++) o_prototype[l][j] = x[j * (*nr) + index];
  }

  free(random_obj_num);
}

// ----- Calculate the cluster dispersion (objective function) -----

float calcCost(double *x, 	// Numeric matrix as vector by col (nr*nc)
	       int *nr, 	// Number of rows
	       int *nc, 	// Number of columns
	       int *k, 		// Number of clusters
	       double *lambda,	// Learning rate
	       int *partition, 	// Partition matrix (nr)
	       double **o_prototype, // Numeric prototype matrix (k*nc)
	       double **subspace_weights) // Weights for variable/cluster (k*nc)
{
  float 
    dispersion = 0.0,  // Dispersion of current cluster
    entropy = 0.0;
  
  int i, j, l;

  for (i = 0; i < *nr; i++) 
    for (j = 0; j < *nc; j++) 
    {
      dispersion += subspace_weights[partition[i]][j] * 
	pow(x[j * (*nr) + i] - o_prototype[partition[i]][j], 2);
    }

  for (l = 0; l < *k; l++) 
    for (j = 0; j < *nc; j++) 
      entropy += subspace_weights[l][j] * log(subspace_weights[l][j]);
  
  dispersion += *lambda * entropy;
  
  return dispersion;
}

// ----- Partition objects into clusters -----

void updPartition(// Inputs
		  double *x, 	// Numeric matrix as vector by col (nr*nc)
		  int *nr, 	// Number of rows
		  int *nc, 	// Number of columns
		  int *k, 	// Number of clusters
		  double **o_prototype, // Numeric prototype matrix (k*nc)
		  double **subspace_weights, // Weights for variable/cluster (k*nc)
		  // Output
		  int *partition)	// Partition matrix (nr)
{
  int i, j, l; 

  // We record the cluster number with the smallest distance to a
  // certain object and store the smallest distence between clusers.

  double 
    o_dist, 
    temp_dist = 0.0, 
    min_dist; 
  
  for (i = 0; i < *nr; i++) 
  {
    min_dist = INT_MAX;
    partition[i] = 0;
    for (l = 0; l < *k; l++) 
    {
      o_dist = 0.0;
      
      for (j = 0; j < *nc; j++) 
      {
	temp_dist = x[j * (*nr) + i] - o_prototype[l][j];
	o_dist += subspace_weights[l][j] * temp_dist * temp_dist;
      }
      
      if (min_dist >= o_dist) 
      {
	min_dist = o_dist;
	partition[i] = l;
      }
    }
  }
}

// --- Update the prototypes -----

int updPrototypes(// Inputs ---------------------------------------------
		    double *x, 	// Numeric matrix as vector by col (nr*nc)
		    int *nr, 		// Number of rows
		    int *nc, 		// Number of columns
		    int *k, 		// Number of clusters
		    int *partition, 	// Partition matrix (nr)
		    // Output ---------------------------------------------
		    double **o_prototype) // Numeric prototype matrix (k*nc)
{
  int i, j, l;
  int *no_clusters;

  no_clusters = (int *) malloc(*k * sizeof(int));

  for (l = 0; l < *k; l++) 
  {
    no_clusters[l] = 0;
    for (j = 0; j < *nc; j++) 
    {
      o_prototype[l][j] = 0;
    }
  }

  for (i = 0; i < *nr; i++) 
  {
    no_clusters[partition[i]]++;
    for (j = 0; j < *nc; j++) 
      o_prototype[partition[i]][j] += x[j * (*nr) + i];
  }

  int flag = 1;
  for (l = 0; l < *k; l++) 
  {
    if (no_clusters[l] == 0) 
    {
      flag = 0;
      break;
    }
    for (j = 0; j < *nc; j++) 
      o_prototype[l][j] /= no_clusters[l];
  }
  free(no_clusters);
  return flag;
}

// ----- Calculate exp^{a}/sum(exp^{a}) -----

void expNormalize(double *a, int length) 
{
  int i;
  double max;
  double sum = 0;

  max = a[0];

  for (i = 0; i < length; i++) 
    if (a[i] > max) 
      max = a[i];

  for (i = 0; i < length; i++) 
  {
    a[i] = exp(a[i] - max);
    sum += a[i];
  }

  for (i = 0; i < length; i++) 
    a[i] /= sum;
}

// ----- Update subspace weights. -----

void updWeights(// Inputs -------------------------------------------------------
		double *x, 	// Numeric matrix as vector by col (nr*nc)
		int *nr, 	// Number of rows
		int *nc, 	// Number of columns
		int *k, 	// Number of clusters
		double *lambda,	// Learning rate
		int *partition,	// Partition matrix (nr)
		double **o_prototype, // Numeric prototype matrix (k*nc)
		// Output -------------------------------------------------------
		double **subspace_weights) // Weights for variable/cluster (k*nc)
{
  float **DJ;
  int i, j, l;
  void expNormalize(double *a, int length);

  DJ = (float **) malloc(*k * sizeof(float*));
  if (!DJ) 
  {
    error("Can't allocate memory for DJ\n");
  }
  for (l = 0; l < *k; l++) 
  {
    DJ[l] = (float *) malloc((*nc) * sizeof(float));
    if (!DJ[l]) 
    {
      error("Can't allocate memory for o_prototype DJ\n");
    }

    for (j = 0; j < *nc; j++) 
    {
      DJ[l][j] = 0;
    }
  }

  for (i = 0; i < *nr; i++) 
  {
    for (j = 0; j < *nc; j++) 
    {
      DJ[partition[i]][j] += subspace_weights[partition[i]][j] * 
	pow((x[j * (*nr) + i] - o_prototype[partition[i]][j]), 2);
    }
  }

  for (l = 0; l < *k; l++) 
  {
    for (j = 0; j < *nc; j++) 
      subspace_weights[l][j] = -((float) DJ[l][j] / (*lambda));//?
    expNormalize(subspace_weights[l], *nc);
  }

  for (l = 0; l < *k; l++) free(DJ[l]);
  free(DJ);
}

// ***** Primary Interface *****

// This is oriented toward interfacing with R, though a
// separate main.c can be used for stand alone testing.

void ewkm(// Inputs ----------------------------------------------------------
	  double *x, 		// Numeric matrix as vector by col (nr*nc)
	  int *nr, 		// Number of rows
	  int *nc, 		// Number of columns
	  int *k, 		// Number of clusters
	  double *lambda, 	// Learning rate
	  int *maxiter, 	// Maximum number of iterations
	  double *delta, 	// Minimum change below which iteration stops
	  int *maxrestart,      // Maximum number of restarts
	  // Outputs ---------------------------------------------------------
	  int *iterations,	// Number of iterations
	  int *cluster, 	// Cluster assignment for each obs (nr)
	  double *centers, 	// Cluster centers (k*nc)
	  double *weights, 	// Variable weights (k*nc)
	  int *restarts,	// Number of restarts
		int *totiters, //
		double *totss, //
		double *withiness) // Number of iterations including restarts
{
  int i, j, l, full;

  int iteration; // Count of iterations.

  float dispersion = INT_MAX, dispersion1 = INT_MAX; // Objective function value.

  double **o_prototype; 	// Numeric prototype matrix  (k*nc)
  int     *partition; 		// Partition matrix (nr)
  double **subspace_weights; 	// (k*nc)
  
  // Allocate all neccessary memory
  
  // -- numeric prototype matrix --

  o_prototype = (double **) malloc(*k * sizeof(double*));
  
  if (!o_prototype) 
  {
    error("Can't allocate memory for o_prototype matrix\n");
  }
  
  for (l = 0; l < *k; l++) 
  {
    o_prototype[l] = (double *) malloc(*nc * sizeof(double));
    if (!o_prototype[l]) 
    {
      error("Can't allocate memory for o_prototype matrix\n");
    }
  }

  // -- partition matrix --

  partition = (int *) malloc(*nr * sizeof(int));

  if (!partition) 
  {
    error("Can't allocate memory for partition matrix\n");
  }

  // -- subspace_weights array --

  subspace_weights = (double **) malloc(*k * sizeof(double*));
  
  if (!subspace_weights) 
  {
    error("Can't allocate memory for subspace_weights matrix\n");
  }

  for (l = 0; l < *k; l++) 
  {
    subspace_weights[l] = (double *) malloc(*nc * sizeof(double));
    
    if (!subspace_weights[l]) 
    {
      error("Can't allocate memory for subspace_weights matrix\n");
    }
  }

  // Read in (or create) .Random.seed, the R random number data, and
  // then initialise the random sequence.

  GetRNGstate();

  // Initialise a rand sequence.

  srand(unif_rand() * RAND_MAX);

  // Initialize the prototypes.

  initPrototypes(x, nr, nc, k, o_prototype);

  // Initialize the feature weights of a cluster.

  for (l = 0; l < *k; l++) 
    for (j = 0; j < *nc; j++)
      subspace_weights[l][j] = 1.0 / *nc;
  
  // Now cluster

  iteration = 0;
  *totiters = 0;
  *restarts = 0;

  while (++iteration <= *maxiter) 
  {
    dispersion = dispersion1;
    
    updPartition(x, nr, nc, k, o_prototype, subspace_weights, partition);

    // Check if any prototypes are empty, and if so we have to
    // initiate a new search if we have restarts left
    
    full = updPrototypes(x, nr, nc, k, partition, o_prototype); 

    if (! full  && *maxrestart != 0)
    {
      *restarts += 1;
      *maxrestart -= 1;
      *totiters += iteration;
      iteration = 0;
  
      // Initialize the prototypes

      initPrototypes(x, nr, nc, k, o_prototype);
  
      // Initialize the feature weights of a cluster.

      for (l = 0; l < *k; l++)
	for (j = 0; j < *nc; j++)
	  subspace_weights[l][j] = 1.0 / *nc;
    }

    // Update weights of attibutes of each cluster

    updWeights(x, nr, nc, k, lambda, partition, o_prototype,
	       subspace_weights);
    
    // Compute objective function value

    dispersion1 = calcCost(x, nr, nc, k, lambda, partition, o_prototype,
			   subspace_weights);
    
    // Check for convergence

    if (fabs(dispersion - dispersion1) / dispersion1 < *delta) break;
  }
  
  // Record results in output variables for passing back to R.
  
  iterations[0] = iteration-1;
  
  for (i = 0; i < *nr; i++) cluster[i] = partition[i];
  free(partition);
  
  for (l = 0; l < *k; l++) 
  {
    for (j = 0; j < *nc; j++) 
    {
      i = j * (*k) + l;
      centers[i] = o_prototype[l][j];
      weights[i] = subspace_weights[l][j];
    }
  }
  
  for (l = 0; l < *k; l++) 
  {
    free(o_prototype[l]);
    free(subspace_weights[l]);
  }
  
  free(o_prototype);
  free(subspace_weights);
  
  *totiters += iteration;
  // If we have reaced the maximum iterations, the count was already
  // increased.
  if (iteration == *maxiter + 1) *totiters = *totiters - 1;

  // Write out the R random number data.

  PutRNGstate();

	/**
	 * calcuate within class sum of squares and total sum of squares
	 */
	sum_squares(x, nr, nc, k, cluster, centers, totss, withiness);
	// Done.

}

