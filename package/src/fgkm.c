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
#include <limits.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include "sum_squares.h"

inline double eu_distance1(const double f1, const double f2) {
	return (f1 - f2) * (f1 - f2);
}

void init_centers(const double *x, const int *nr, const int *nc, const int *k,
		double *centers) {
	int i, j, l, flag, index, *random_obj_num;
	random_obj_num = (int *) malloc(*k * sizeof(int));

	if (!random_obj_num) {
		error("can't allocate random_obj_num\n");
	}

	for (l = 0; l < *k; ++l) {
		random_obj_num[l] = -1;
	}

	for (l = 0; l < *k; ++l) {
		flag = 1;

		while (flag) {
			index = (int) (rand() % (*nr));
			flag = 0;
			for (i = 0; i < l; ++i) {
				if (random_obj_num[i] == index)
					flag = 1;
			}
		}
		random_obj_num[l] = index;

		for (j = 0; j < *nc; ++j)
			centers[j * (*k) + l] = x[j * (*nr) + index];

	} // for l = 1:k
	free(random_obj_num);
}

void init_featureWeight(double *featureWeight, const int *k, const int *nc) {
	int l, j;
	for (l = 0; l < *k; ++l)
		for (j = 0; j < *nc; ++j)
			featureWeight[j * (*k) + l] = 1.0 / (*nc);
}

void init_groupWeight(double *groupWeight, const int *k, const int *numGroups) {
	int l, t;
	for (l = 0; l < *k; ++l)
		for (t = 0; t < *numGroups; ++t)
			groupWeight[t * (*k) + l] = 1.0 / (*numGroups);
}

void update_cluster(const double *x, const int *nr, const int *nc, const int *k,
		const int *numGroups, const int *groupInfo, int *cluster,
		const double *centers, const double *featureWeight,
		const double *groupWeight) {
	int i, j, l;
	double min_dist, o_dist;

	for (i = 0; i < *nr; ++i) {

		min_dist = (double) INT_MAX;
		for (l = 0; l < *k; ++l) {

			o_dist = 0.0;
			for (j = 0; j < *nc; ++j) {
				o_dist += featureWeight[j * (*k) + l]
						* groupWeight[groupInfo[j] * (*k) + l]
						* eu_distance1(centers[j * (*k) + l], x[j * (*nr) + i]);
			}
			if (o_dist <= min_dist) {
				min_dist = o_dist;
				cluster[i] = l;
			}
		}
	}
}

int update_centers(const double *x, const int *nr, const int *nc, const int *k,
		const int *cluster, double *centers) {
	int i, j, l, *no_cluster;
	no_cluster = (int *) malloc((*k) * sizeof(int));

	for (l = 0; l < *k; ++l) {
		no_cluster[l] = 0;
		for (j = 0; j < *nc; ++j) {
			centers[j * (*k) + l] = 0.0;
		}
	}

	for (i = 0; i < *nr; ++i) {
		no_cluster[cluster[i]]++;
		for (j = 0; j < *nc; ++j) {
			centers[j * (*k) + cluster[i]] += x[j * (*nr) + i];
		}
	}

	int flag = 1;
	for (l = 0; l < *k; ++l) {
		if (no_cluster[l] == 0) {
			flag = 0;
			break;
		}
		for (j = 0; j < *nc; ++j) {
			centers[j * (*k) + l] /= no_cluster[l];
		}
	}
	free(no_cluster);
	return flag;
}

void update_featureWeight(const double *x, const int *nr, const int *nc,
		const int *k, const double *eta, const int *numGroups,
		const int *groupInfo, const int *nums, const int *cluster,
		const double *centers, double *featureWeight, const double *groupWeight) {
	int i, j, l, t;

	double **E;
	E = (double **) malloc(*k * sizeof(double *));
	if (!E) {
		error("can not allocate E[][].\n");
	}
	for (l = 0; l < *k; ++l) {
		E[l] = (double *) malloc(*nc * sizeof(double));
		if (!E[l]) {
			error("can not allocate E[][].\n");
		}
	}

	for (l = 0; l < *k; ++l) {
		for (j = 0; j < *nc; ++j) {
			E[l][j] = 0.0;

			for (i = 0; i < *nr; ++i)
				if (cluster[i] == l)
					E[l][j] += groupWeight[groupInfo[j] * (*k) + l] // groupInfo[j]==t
					* eu_distance1(x[j * (*nr) + i], centers[j * (*k) + l]);
		}

	}

	double sum = 0.0, max = 0.0;
	int f, skip;
	for (l = 0; l < *k; ++l)
		for (j = 0; j < *nc; ++j)
			featureWeight[j * (*k) + l] = (-E[l][j] / (*eta));

	for (l = 0; l < *k; ++l) { // every CLUSTER
		skip = 0;
		for (t = 0; t < *numGroups; ++t) { // for every feature group
			sum = 0.0;
			max = featureWeight[(0 + skip) * (*k) + l];

			for (f = 0; f < nums[t]; ++f) {
				if (featureWeight[(f + skip) * (*k) + l] >= max)
					max = featureWeight[(f + skip) * (*k) + l];
			}

			for (f = 0; f < nums[t]; ++f) {
				featureWeight[(f + skip) * (*k) + l] = exp(
						featureWeight[(f + skip) * (*k) + l] - max);
				sum += featureWeight[(f + skip) * (*k) + l];
			}

			for (f = 0; f < nums[t]; ++f) {
				featureWeight[(f + skip) * (*k) + l] /= sum;
			}

			skip += nums[t]; // skip to next feature group
		}
	}

	for (l = 0; l < *k; ++l) {
		free(E[l]);
	}
	free(E);
}

void update_groupWeight(const double *x, const int *nr, const int *nc,
		const int *k, const double *lambda, const int *numGroups,
		const int *groupInfo, const int *cluster, const double *centers,
		const double *featureWeight, double *groupWeight) {
	int i, j, l, t;
	double **D;
	D = (double **) malloc(*k * sizeof(double *));
	if (!D) {
		error("can not allocate group weight!\n");
	}
	for (l = 0; l < *k; ++l) {
		D[l] = (double *) malloc(*numGroups * sizeof(double));
		if (!D[l]) {
			error("can not allocate group weight!\n");
		}
	}


	for (l = 0; l < *k; ++l) {
		for (t = 0; t < *numGroups; ++t) {

			D[l][t] = 0;

			for (i = 0; i < *nr; ++i)
				for (j = 0; j < *nc; ++j)
					if (cluster[i] == l && groupInfo[j] == t)
						D[l][t] += featureWeight[j * (*k) + l]
								* eu_distance1(centers[j * (*k) + l],
										x[j * (*nr) + i]);
		}
	}

	for (l = 0; l < *k; ++l)
		for (t = 0; t < *numGroups; ++t)
			groupWeight[t * (*k) + l] = (-D[l][t]) / (*lambda);

	double sum = 0.0, max = 0.0;

	// implement expNormalize()
	for (l = 0; l < *k; ++l) {
		sum = 0.0;
		max = groupWeight[0 * (*k) + l]; // initially assign gw[l][0] to max
		for (t = 0; t < *numGroups; ++t) {
			if (groupWeight[t * (*k) + l] >= max)
				max = groupWeight[t * (*k) + l];
		}

		for (t = 0; t < *numGroups; ++t) {
			groupWeight[t * (*k) + l] = exp(groupWeight[t * (*k) + l] - max);
			sum += groupWeight[t * (*k) + l];
		}

		for (t = 0; t < *numGroups; ++t) {
			groupWeight[t * (*k) + l] /= sum;
		}
	}

	for (l = 0; l < *k; ++l) {
		free(D[l]);
	}
	free(D);
}

double calculate_cost(const double *x, const int *nr, const int *nc,
		const int *k, const double *lambda, const double *eta,
		const int *numGroups, const int *groupInfo, const int *cluster,
		const double *centers, const double *featureWeight,
		const double *groupWeight) {

	int i, j, l, t;
	double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0, dispersion;

	for (l = 0; l < *k; ++l) {

		for (i = 0; i < *nr; ++i) {
			for (t = 0; t < *numGroups; ++t) {
				for (j = 0; j < *nc; ++j) {
					if (groupInfo[j] == t && cluster[i] == l)
						sum0 += groupWeight[t * (*k) + l]
								* featureWeight[j * (*k) + l]
								* eu_distance1(centers[j * (*k) + l],
										x[j * (*nr) + i]);
				}
			}
		}
		for (t = 0; t < *numGroups; ++t)
			sum1 += groupWeight[t * (*k) + l] * log(groupWeight[t * (*k) + l]);

		for (j = 0; j < *nc; ++j)
			sum2 += featureWeight[j * (*k) + l]
					* log(featureWeight[j * (*k) + l]);

	} // for l = 1:k

	dispersion = sum0 + sum1 * (*lambda) + sum2 * (*eta);
	return dispersion;
}

void parse_group(const char *strGroup, int *numGroups, int *nums,
		int *groupInfo) {

	int length, first, last, i = 0;
	char Buffer[20], *p1;

	length = strlen(strGroup);
	i = 0;
	*numGroups = 0; // important !!!

	while (i < length) {

		p1 = &Buffer[0];

		while (strGroup[i] != '-') {
			*p1 = strGroup[i];
			p1++;
			i++;

		}

		*p1 = '\0';

		first = atoi(Buffer);

		p1 = &Buffer[0];
		i++;
		while (strGroup[i] != ':' && i < length) {
			*p1 = strGroup[i];
			p1++;
			i++;

		}

		*p1 = '\0';
		last = atoi(Buffer);
		nums[*numGroups] = last - first + 1;
		(*numGroups)++;
		i++;
	}

	int g, t;
	i = 0;
	for (t = 0; t < *numGroups; ++t) {
		for (g = 0; g < nums[t]; ++g) {
			groupInfo[i] = t;
			i++;
		}
	}
}

void fgkm(const double *x, const int *nr, const int *nc, const int *k,
		const double *lambda, const double *eta, const char **strGroup,
		const double *delta, const int *maxiter, const int *maxrestart,
		int *cluster, double *centers, double *featureWeight,
		double *groupWeight, int *iterations, int *restarts, int *totiter,
		double *totalCost, //
		double *totss, //    total sum of squares
		double *withiness // vector of sum of square in every cluster
		) {
	int *numGroups, *nums, *groupInfo;
	numGroups = (int *) malloc(sizeof(int));
	nums = (int *) malloc(100 * sizeof(int));
	groupInfo = (int *) malloc(*nc * sizeof(int));

	parse_group((*strGroup), numGroups, nums, groupInfo); // get 'numGroups'

	double dispersion1, dispersion2;
	int flag_not_restart = 1;

	*restarts = 0;
	*iterations = 0;
	*totiter = 0;

	srand((unsigned) time(NULL));

	while ((*restarts) < *maxrestart) {

		init_centers(x, nr, nc, k, centers); // assign randomly
		init_featureWeight(featureWeight, k, nc); // equal value
		init_groupWeight(groupWeight, k, numGroups); // equal value

		*iterations = 0;
		dispersion2 = (double) INT_MAX;

		while ((*iterations) < *maxiter) {
			Rprintf("*");
			(*iterations)++;
			(*totiter)++;
			dispersion1 = dispersion2;

			update_cluster(x, nr, nc, k, numGroups, groupInfo, cluster, centers,
					featureWeight, groupWeight);

			flag_not_restart = update_centers(x, nr, nc, k, cluster, centers);

			if (!flag_not_restart) {
				(*restarts)++;
				break;
			}

			update_featureWeight(x, nr, nc, k, eta, numGroups, groupInfo, nums,
					cluster, centers, featureWeight, groupWeight);

			update_groupWeight(x, nr, nc, k, lambda, numGroups, groupInfo,
					cluster, centers, featureWeight, groupWeight);

			dispersion2 = calculate_cost(x, nr, nc, k, lambda, eta, numGroups,
					groupInfo, cluster, centers, featureWeight, groupWeight);

			// if change of dispersion below delta or iterations exceed max iterations,
			//	then, terminate and return.
			//  but if dispersion < 0, then ???
			if ((fabs((dispersion1 - dispersion2) / dispersion1)) <= (*delta)
					|| (*iterations) == *maxiter) {
				Rprintf("Clustering converged. Terminate!\n");

				*totalCost = dispersion1;
				free(groupInfo);
				free(nums);
				free(numGroups);

				/**
				 * calculate sum of squares
				 */
				sum_squares(x, nr, nc, k, cluster, centers, totss, withiness);
				// Done.

				return;
			}
		}
	}

	free(groupInfo);
	free(nums);
	free(numGroups);

}
