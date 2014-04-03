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

/**
 * This function calculate the total sum squares and
 * within class sum squares of every cluster.
 *
 * @author Longfei Xiao <lf.xiao@siat.ac.cn>
 */

#include <stdlib.h>
void sum_squares(double *x, // Data matrix
		int *nr, // number of rows
		int *nc, // number of columns
		int *k, // number of clusters
		int *cluster, // A vector of integers (from 1:k) indicating the cluster to which each point is allocated.
		double *centers, // Cluster centers of every cluster
		double *totss, // The total sum of squares.
		double *withiness // Vector of within-cluster sum of squares, one component per cluster.
		);
