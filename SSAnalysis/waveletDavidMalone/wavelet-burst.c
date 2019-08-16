/*
 * Wavelet analysis of time series, prints out power spectrum estimator
 * as described by Abry and Veitch in IEEE Trans Info Th. Vol 44 No. 1 Jan 1988.
 * This program specialises in bursts, and reads a list of
 * (burstlen, value) pairs.
 * usage: wavelet-burst family [filename]
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
 * $Id: wavelet-burst.c,v 1.5 2006/01/25 10:23:41 dwmalone Exp $
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SCALE 2

struct wavelet {
	const char *name;
	int len;
	double *c;
};

void list_wavelets(FILE *out);
struct wavelet *find_wavelet(const char *name);

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

void push_val(double x, int len_hint);
void push_burst(double x, double len);

/* Window for wavelet analysis. */
static double *h;
int window_pos;

/* Estimator so far. */
static double nspect, spect;

/* The bursts and values for the next level. */
static struct values **nprevp, *cur;

/* Wavelet details. */
static int wlen;
static double *c, *crm;

int
main(int argc, char **argv) {
	FILE *input;
	struct wavelet *w;
	struct values *head, **prevp, *this, *next;
	double blen, bval;
	double x;
	int k;

	if (argc != 2 && argc != 3) {
		fprintf(stderr, "usage: %s wavelet [file]\n", argv[0]);
		exit(1);
	}

	if ((w = find_wavelet(argv[1])) == NULL) {
		fprintf(stderr,
		    "'%s' isn't a known wavelet family.\nKnown wavelets:\n",
		    argv[1]);
		list_wavelets(stderr);
		exit(1);
	}

	/*
	 * Set up input.
	 */
	if (argc == 3) {
		if ((input = fopen(argv[2], "r")) == NULL) {
			perror("Couldn't open input");
			exit(1);
		}
	} else {
		input = stdin;
	}

	/*
	 * Set up filters.
	 */
	wlen = w->len;
	c = w->c;
	if ((crm = malloc(sizeof(double) * wlen)) == NULL) {
		perror("Couln't malloc reversed coeff.\n");
		exit(1);
	}
	for (k = 0, x = 1; k < wlen; k++, x *= -1)
		crm[wlen-1-k] = c[k]*x;

	/*
	 * Set up the window.
	 */
	if ((h = malloc(sizeof(double) * wlen)) == NULL) {
		perror("Couln't malloc window.\n");
		exit(1);
	}
	window_pos = 0;

	/*
	 * Spool in data.
	 */
	head = NULL;
	prevp = &head;
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
	}

	/*
	 * Do the convolution.
	 */
	do {
		window_pos = 0;
		spect = 0.0;
		nspect = 0;
		this = head;
		head = NULL;
		cur = NULL;
		nprevp = &head;
		/*
		 * Walk the values list, setting up new values and doing
		 * estimation.
		 */
		for (; this != NULL; this = next) {
			next = this->next;
			switch (this->type) {
			case BURST:
				push_burst(this->b_val, this->b_len);
				break;
			case LIST:
				for (k = 0; k < this->l_len; k++)
					push_val(this->l_val[k], this->l_len-k);
				free(this->l_val);
				break;
			default:
				fprintf(stderr, "Unknown type %d.\n", this->type);
				exit(1);
			}
			free(this);
		}
		printf("%f\n", log(spect/nspect)/log(SCALE));
	} while (nspect >= wlen);
	exit(0);
}

void
push_val(double x, int len_hint)
{
	double ht, dt;
	int k;

	h[window_pos++] = x;
	if (window_pos < wlen)
		return;

	ht = 0.0;
	dt = 0.0;
	for (k = 0; k < wlen; k++) {
		ht += h[k] * c[k];
		dt += h[k] * crm[k];
	}

	/* Figure out where to put ht. */
	switch (cur != NULL ? cur->type : NONE) {
	case BURST:
		if (cur->b_val == ht) {
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
		cur->l_val[cur->l_len++] = ht;
	}
	/* Update the estimate. */
	spect += dt*dt;
	nspect++;

	/* Shift the window. */
	for (k = 0; k < wlen - SCALE; k++)
		h[k] = h[k+SCALE];
	window_pos -= SCALE;
}

void
push_burst(double x, double len) {
	int pushes, k;
	double newlen, newval;

	/* Try to fill out window with the value. */
	pushes =  wlen + (window_pos % SCALE) - SCALE;
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
	 * to get the next output, which will be 0.
	 */
	newlen = floor(len/SCALE);
	pushes = rint(fmod(len,SCALE));
	newval = 1.4142135623730950488*x; /* If you change SCALE, check me. */

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
	spect += 0.0; /* Wavelets are mean zero. */
	nspect += newlen;

	/*
	 * Push remainder of the burst that didn't make it into the new burst.
	 */
	for (k = 0; k < pushes; k++)
		push_val(x, pushes-k);
}

/* 
 * Numerical versions of the filters can be found at sites like:
 * http://deneb.kuee.kyoto-u.ac.jp/~yoneyone/wavelet/
 * http://www.isye.gatech.edu/~brani/filters.html
 * These particular ones were taken from Thomas Nelson's tnimage data files.
 */

double daub4[] = {
	 0.4829629131445341,
	 0.8365163037378079,
	 0.2241438680420134,
	-0.1294095225512604};

double daub6[] = {
	 0.3326705529500825,
	 0.8068915093110924,
	 0.4598775021184914,
	-0.1350110200102546,
	-0.0854412738820267,
	 0.0352262918857095};

double daub8[] = {
	 0.2303778133088964,
	 0.7148465705529154,
	 0.6308807679398587,
	-0.0279837694168599,
	-0.1870348117190931,
	 0.0308413818355607,
	 0.0328830116668852,
	-0.0105974017850690};

double daub10[] = {
	 0.16010239797419231755,
	 0.60382926979718765104,
	 0.72430852843777138172,
	 0.13842814590132229702,
	-0.24229488706637955509,
	-0.03224486958463790298,
	 0.07757149384004533021,
	-0.00624149021279823274,
	-0.01258075199908191381,
	 0.00333572528547376171};

double daub12[] = {
	 0.11154074335010408237,
	 0.49462389039843135397,
	 0.75113390802107049549,
	 0.31525035170920534533,
	-0.22626469396540713208,
	-0.12976686756724714611,
	 0.09750160558731878202,
	 0.02752286553030498448,
	-0.03158203931748377463,
	 0.00055384220116140592,
	 0.00477725751094519764,
	-0.00107730108530840587};

double daub14[] = {
	 0.07785205408501189028,
	 0.39653931948193121837,
	 0.72913209084625985046,
	 0.46978228740520439066,
	-0.14390600392858410306,
	-0.22403618499390165475,
	 0.07130921926682094736,
	 0.08061260915108464653,
	-0.03802993693501364320,
	-0.01657454163066688843,
	 0.01255099855610011805,
	 0.00042957797292139123,
	-0.00180164070404753161,
	 0.00035371379997452691};

double daub16[] = {
	 0.0544158422431072,
	 0.3128715909143166,
	 0.6756307362973195,
	 0.5853546836542159,
	-0.0158291052563823,
	-0.2840155429615824,
	 0.0004724845739124,
	 0.1287474266204893,
	-0.0173693010018090,
	-0.0440882539307971,
	 0.0139810279174001,
	 0.0087460940474065,
	-0.0048703529934520,
	-0.0003917403733770,
	 0.0006754494064506,
	-0.0001174767841248};

double daub18[] = {
	 0.0380779473638778,
	 0.2438346746125858,
	 0.6048231236900955,
	 0.6572880780512736,
	 0.1331973858249883,
	-0.2932737832791663,
	-0.0968407832229492,
	 0.1485407493381256,
	 0.0307256814793385,
	-0.0676328290613279,
	 0.0002509471148340,
	 0.0223616621236798,
	-0.0047232047577518,
	-0.0042815036824635,
	 0.0018476468830563,
	 0.0002303857635232,
	-0.0002519631889427,
	 0.0000393473203163};

double daub20[] = {
	 0.0266700579005473,
	 0.1881768000776347,
	 0.5272011889315757,
	 0.6884590394534363,
	 0.2811723436605715,
	-0.2498464243271598,
	-0.1959462743772862,
	 0.1273693403357541,
	 0.0930573646035547,
	-0.0713941471663501,
	-0.0294575368218399,
	 0.0332126740593612,
	 0.0036065535669870,
	-0.0107331754833007,
	 0.0013953517470688,
	 0.0019924052951925,
	-0.0006858566949564,
	-0.0001164668551285,
	 0.0000935886703202,
	-0.0000132642028945};

double daub68[] = {
	 0.000005770510509196372,
	 0.000129947617286060910,
	 0.001364061360857502900,
	 0.008819889215070597600,
	 0.039048840515836347000,
	 0.124152479453546570000,
	 0.287765053073302140000,
	 0.478478736036207890000,
	 0.530555088298456210000,
	 0.290366323291261220000,
	-0.128246839428697380000,
	-0.331525294405338620000,
	-0.103891913282997320000,
	 0.216907215581275150000,
	 0.166601746989076480000,
	-0.127337355243722190000,
	-0.160924923462106410000,
	 0.077991845108642127000,
	 0.134125957419103990000,
	-0.054482963113175950000,
	-0.102947593726802140000,
	 0.043576097709937658000,
	 0.073185237455317007000,
	-0.037012834328747815000,
	-0.047438560026982685000,
	 0.030739754472021749000,
	 0.027228350059302471000,
	-0.023671732546529257000,
	-0.013143972588625805000,
	 0.016409377976107239000,
	 0.004713643252652899500,
	-0.010045501052522497000,
	-0.000619476731810230350,
	 0.005334950720477601800,
	-0.000769215319245805750,
	-0.002399456109851749600,
	 0.000858994417719352350,
	 0.000875199041890761700,
	-0.000552735446532307170,
	-0.000232673197326024440,
	 0.000265077237852853530,
	 0.000026600548344633845,
	-0.000099146977262800837,
	 0.000013531187821126573,
	 0.000028449515523895963,
	-0.000010576574554128021,
	-0.000005710825840354940,
	 0.000004169871888982280,
	 0.000000497971843697713,
	-0.000001116306485312163,
	 0.000000144819571623039,
	 0.000000202599062939615,
	-0.000000075267015756319,
	-0.000000019903464623409,
	 0.000000017404232955792,
	-0.000000000866574406854,
	-0.000000002316501897466,
	 0.000000000644637807231,
	 0.000000000130041029079,
	-0.000000000099047743256,
	 0.000000000010042087140,
	 0.000000000006080125224,
	-0.000000000002107879064,
	 0.000000000000097994509,
	 0.000000000000085791939,
	-0.000000000000023170837,
	 0.000000000000002587338,
	-0.000000000000000114894};


double coif6[] = {
	-0.07273261951285,
	 0.33789766245781,
	 0.85257202021226,
	 0.38486484686420,
	-0.07273261951285,
	-0.01565572813546};


double coif12[] = {
	 0.01638733646360,
	-0.04146493678197,
	-0.06737255472230,
	 0.38611006682309,
	 0.81272363544961,
	 0.41700518442378,
	-0.07648859907869,
	-0.05943441864675,
	 0.02368017194645,
	 0.00561143481942,
	-0.00182320887071,
	-0.00072054944537};

double coif18[] = {
	-0.00379351286449,
	 0.00778259642733,
	 0.02345269614184,
	-0.06577191128186,
	-0.06112339000267,
	 0.40517690240962,
	 0.79377722262562,
	 0.42848347637762,
	-0.07179982161931,
	-0.08230192710689,
	 0.03455502757306,
	 0.01588054486362,
	-0.00900797613666,
	-0.00257451768875,
	 0.00111751877089,
	 0.00046621696011,
	-0.00007098330314,
	-0.00003459977284};

double coif24[] = {
	 0.00089231366858,
	-0.00162949201260,
	-0.00734616632764,
	 0.01606894396478,
	 0.02668230015605,
	-0.08126669968088,
	-0.05607731331675,
	 0.41530840703043,
	 0.78223893092050,
	 0.43438605649147,
	-0.06662747426343,
	-0.09622044203399,
	 0.03933442712334,
	 0.02508226184486,
	-0.01521173152795,
	-0.00565828668661,
	 0.00375143615728,
	 0.00126656192930,
	-0.00058902075624,
	-0.00025997455249,
	 0.00006233903446,
	 0.00003122987587,
	-0.00000325968024,
	-0.00000178498500};

double coif30[] = {
	-0.00021208083983,
	 0.00035858968793,
	 0.00217823635833,
	-0.00415935878180,
	-0.01013111752086,
	 0.02340815678818,
	 0.02816802897375,
	-0.09192001056889,
	-0.05204316318145,
	 0.42156620673301,
	 0.77428960373039,
	 0.43799162621564,
	-0.06203596396911,
	-0.10557420871390,
	 0.04128920875431,
	 0.03268357427038,
	-0.01976177894455,
	-0.00916423116340,
	 0.00676418544873,
	 0.00243337321290,
	-0.00166286370218,
	-0.00063813134311,
	 0.00030225958184,
	 0.00014054114972,
	-0.00004134043228,
	-0.00002131502681,
	 0.00000373465518,
	 0.00000206376185,
	-0.00000016744289,
	-0.00000009517657};

#define LEN(array) ((int)(sizeof(array)/sizeof(array[0])))
#define WAVELET(f) { #f, LEN(f), f }

struct wavelet wavelets[] = {
	WAVELET(daub4),
	WAVELET(daub6),
	WAVELET(daub8),
	WAVELET(daub10),
	WAVELET(daub12),
	WAVELET(daub14),
	WAVELET(daub16),
	WAVELET(daub18),
	WAVELET(daub20),
	WAVELET(daub68),
	WAVELET(coif6),
	WAVELET(coif12),
	WAVELET(coif18),
	WAVELET(coif24),
	WAVELET(coif30)
};

void
list_wavelets(FILE *out) {
	int i;

	for (i = 0; i < LEN(wavelets); i++) {
		fprintf(out, "\t%s\t(n = %d, c0 = %f)\n",
		    wavelets[i].name,
		    wavelets[i].len,
		    wavelets[i].len <= 0 ? 0.0 : wavelets[i].c[0]);
	}
}

struct wavelet *
find_wavelet(const char *name) {
	int i;

	for (i = 0; i < LEN(wavelets); i++) {
		if (strcmp(name, wavelets[i].name) == 0)
			return &wavelets[i];
	}
	return NULL;
}
