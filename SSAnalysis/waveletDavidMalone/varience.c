/*
 * Estimate the varience of blocks of length SCALE^m taking burst input.
 * usage: varience family [filename]
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SCALE 2

int
main(int argc, char **argv) {
	FILE *input;
	double x;
	double num, sum, ssum, agg;
	double *h;
	int j;
	int len, datalen;

	if (argc != 1 && argc != 2) {
		fprintf(stderr, "usage: %s [file]\n", argv[0]);
		exit(1);
	}

	/*
	 * Set up input.
	 */
	if (argc == 2) {
		if ((input = fopen(argv[2], "r")) == NULL) {
			perror("Couldn't open input");
			exit(1);
		}
	} else {
		input = stdin;
	}

	/*
	 * Spool in data.
	 */
	len = 0;
	datalen = 0;
	h = NULL;
	while (fscanf(input, "%lf", &x) == 1) {
		if (len >= datalen) {
			int n;
			
			n = datalen * 2 + 10;
			if ((h = realloc(h, n*sizeof(*h))) == NULL) {
				perror("Couldn't expland h");
				exit(1);
			}
			datalen = n;
		}
		h[len++] = x;
	}

	/*
	 * Calculate the varience and aggregate.
	 */
	while (len >= 1) {
		num = 0;
		sum = 0;
		agg = 0;
		ssum = 0;
		for (j = 0; j < len; j++) {
			num++;
			sum += h[j];
			agg += h[j];
			ssum += h[j]*h[j];
			if (j % SCALE == SCALE - 1) {
				h[j/SCALE] = agg/SCALE;
				agg = 0;
			}
		}
		printf("%f\n", log(ssum/num - sum*sum/num/num));
		len = j/SCALE;
	}
	exit(0);
}
