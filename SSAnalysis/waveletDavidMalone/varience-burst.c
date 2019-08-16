/*
 * Estimate the varience of blocks of length SCALE^m taking burst input.
 * usage: varience-burst [filename]
 */

/*
 * Copyright (C) 2001, 2002, 2004, 2006 David Malone <David.Malone@nuim.ie>
 * All rights reserved.
 *
 * The development of this software was supported by Science Foundation Ireland
 * under the National Development Plan.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the project nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * $Id: varience-burst.c,v 1.1 2006/09/22 13:36:56 dwmalone Exp $
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SCALE 2

struct values {
	int type;
#define LIST 0
#define BURST 1
#define NONE 2
	double *l_val;
	int l_num;
	int l_len;
	double b_val;
	double b_len;
	struct values *next;
};

int window_pos = 0;
double agg_val = 0;

void push_val(double x, int len_hint);
void push_burst(double x, double len);

/* The bursts and values for the next level. */
static struct values **nprevp, *cur;

int
main(int argc, char **argv) {
	FILE *input;
	struct values *head, **prevp, *this, *next;
	double blen, bval;
	double num, sum, ssum;
	double total_len;
	int k;

	if (argc != 1 && argc != 2) {
		fprintf(stderr, "usage: %s varience-burst [file]\n", argv[0]);
		exit(1);
	}

	/*
	 * Set up input.
	 */
	if (argc == 2) {
		if ((input = fopen(argv[1], "r")) == NULL) {
			perror("Couldn't open input");
			exit(1);
		}
	} else {
		input = stdin;
	}

	/*
	 * Spool in data.
	 */
	head = NULL;
	prevp = &head;
	total_len = 0;
	while (fscanf(input, "%lf %lf", &blen, &bval) == 2) {
		if ((this = malloc(sizeof(*this))) == NULL) {
			perror("Couldn't allocate next burst");
			exit(1);
		}
		this->type = BURST;
		this->b_len = floor(blen);
		this->b_val = bval;
		this->next = NULL;
		*prevp = this;
		prevp = &(this->next);
		total_len += this->b_len;
	}

	/*
	 * Estimate the varience and the double up.
	 */
	for (; total_len >= 1; total_len /= 2) {
		window_pos = 0;
		agg_val = 0;
		num = 0;
		sum = 0;
		ssum = 0;
		this = head;
		head = NULL;
		cur = NULL;
		nprevp = &head;
		for (; this != NULL; this = next) {
			next = this->next;
			switch (this->type) {
			case BURST:
				/* printf("B %f %f\n", this->b_len, this->b_val); */
				num += this->b_len;
				sum += this->b_val * this->b_len;
				ssum += this->b_val * this->b_val * this->b_len;
				push_burst(this->b_val, this->b_len);
				break;
			case LIST:
				/* printf("L[%d] ", this->l_len); */
				num += this->l_len;
				for (k = 0; k < this->l_len; k++) {
					/* printf("%f ", this->l_val[k]); */
					sum += this->l_val[k];
					ssum += this->l_val[k] * this->l_val[k];
					push_val(this->l_val[k], this->l_len-k);
				}
				/* printf("\n"); */
				free(this->l_val);
				break;
			default:
				fprintf(stderr, "Unknown type %d.\n", this->type);
				exit(1);
			}
			free(this);
		}
		printf("%f\n", log(ssum/num - sum*sum/num/num));
	}
	exit(0);
}

void
push_val(double x, int len_hint)
{
	agg_val += x;
	window_pos++;
	if (window_pos < SCALE)
		return;

	agg_val /= SCALE;

	/* Figure out where to put ht. */
	switch (cur != NULL ? cur->type : NONE) {
	case BURST:
		if (cur->b_val == agg_val) {
			cur->b_len++;
			break;
		}
		/* FALLTHROUGH */
	case NONE:
		if ((cur = malloc(sizeof(*cur))) == NULL) {
			perror("Couldn't allocate next level list");
			exit(1);
		}
		cur->type = LIST;
		cur->l_len = 0;
		cur->l_num = 0;
		cur->l_val = NULL;
		cur->next = NULL;
		*nprevp = cur;
		nprevp = &(cur->next);
		/* FALLTHROUGH */
	case LIST:
		if (cur->l_len >= cur->l_num) {
			int n;
			n = cur->l_len + len_hint + 1;
			if ((cur->l_val =
			     realloc(cur->l_val, n*sizeof(double)))
			    == NULL) {
				perror("Couldn't expand value list");
				exit(1);
			}
			cur->l_num = n;
		}
		cur->l_val[cur->l_len++] = agg_val;
	}

	agg_val = 0;
	window_pos -= SCALE;
}

void
push_burst(double x, double len) {
	int pushes, k;
	double newlen, newval;

	/* Try to fill out window with the value. */
	pushes =  SCALE - window_pos;
	if (len < pushes)
		pushes = len;
	for (k = 0; k < pushes; k++)
		push_val(x, pushes-k);
	len -= pushes;

	/* Have we used all of the burst? */
	if (len < 1)
		return;

	/*
	 * The window has been filled by the burst.
	 * We need to add exactly SCALE more values
	 * to get the next output, which will have the same mean.
	 */
	newlen = floor(len/SCALE);
	pushes = rint(fmod(len,SCALE));
	newval = x;

	if (newlen > 0) {
		if (cur && cur->type == BURST && cur->b_val == newval) {
			cur->b_len += newlen;
		} else {
			if ((cur = malloc(sizeof(*cur))) == NULL) {
				perror("Couldn't allocate next level burst");
				exit(1);
			}
			cur->type = BURST;
			cur->b_len = newlen;
			cur->b_val = newval;
			cur->next = NULL;
			*nprevp = cur;
			nprevp = &(cur->next);
		}
	}

	/*
	 * Push remainder of the burst that didn't make it into the new burst.
	 */
	for (k = 0; k < pushes; k++)
		push_val(x, pushes-k);
}
