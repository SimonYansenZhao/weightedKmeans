// Copyright (c) 2012 Shenzhen Institutes of Advanced Technology
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

#include <stdlib.h>

void sum_squares(double *x, // Data matrix
		int *nr, // number of rows
		int *nc, // number of columns
		int *k, // number of clusters
		int *cluster, // A vector of integers (from 1:k) indicating the cluster to which each point is allocated.
		double *centers, // Cluster centers of every cluster
		double *totss, // The total sum of squares.
		double *withiness // Vector of within-cluster sum of squares, one component per cluster.
		) {
  //printf("Running sum of squares\n");
	int i, j;

	// vector to save average of all objects ---------------------------
	double* global_average = (double *) malloc(*nc * sizeof(double));
	for (j = 0; j < *nc; ++j) {
		global_average[j] = 0.0;
	}

	for (i = 0; i < *nr; ++i) {
		for (j = 0; j < *nc; ++j) {
			// global_average[j] += x[i][j]
			global_average[j] += x[j * (*nr) + i];
		}
	}
	for (j = 0; j < *nc; ++j) {
		global_average[j] = global_average[j] / (*nr);
	}

	// calculate total sum of square ---------------------------------

	*totss = 0.0;
	double tmp_ss, temp;
	for (i = 0; i < (*nr); ++i) {
		tmp_ss = 0.0;
		for (j = 0; j < (*nc); ++j) {
			temp = global_average[j] - x[j * (*nr) + i];
			tmp_ss += temp * temp;
		}
		*totss += tmp_ss;
	}

	// calculate withiness --------------------------------------------
	int t;
	for (t = 0; t < (*k); ++t) {
		withiness[t] = 0.0;
	}

	for (i = 0; i < (*nr); ++i) {

		// calculate the sum of square of object i and the center of the cluster which object i in.
		tmp_ss = 0;
		for (j = 0; j < (*nc); ++j) {
			temp = x[j * (*nr) + i] - centers[j * (*k) + cluster[i]];
			tmp_ss += temp * temp;
		}

		// sum the with in class sum of square of the cluster which object i in.
		withiness[cluster[i]] += tmp_ss;
	}


	free(global_average);
	return;
}

