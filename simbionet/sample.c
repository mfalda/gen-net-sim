#include "sample.h"

#define g_HL globali.sample.HL
#define g_q globali.sample.q
#define g_x globali.sample.x
#define g_y globali.sample.y

// Algoritmi di R (sample.c)

/* Equal probability sampling; with-replacement case */

static void SampleReplace(int k, int n, int *y)
{
	int i;
	for (i = 0; i < k; i++)
		y[i] = (int) n * unif_rand(); // + 1;
}

/* Equal probability sampling; without-replacement case */

static void SampleNoReplace(int k, int n, int *y, int *x)
{
	int i, j;
	for (i = 0; i < n; i++)
		x[i] = i;
	for (i = 0; i < k; i++) {
		j = (int) n * unif_rand();
		y[i] = x[j]; // + 1;
		x[j] = x[--n];
	}
}

/* Unequal probability sampling; without-replacement case */

static void ProbSampleNoReplace(int n, double *p, int *perm,
                                int nans, int *ans)
{
	double rT, mass, totalmass;
	int i, j, k, n1;

	/* Record element identities */
	for (i = 0; i < n; i++)
		perm[i] = i; // + 1;

	/* Sort probabilities into descending order */
	/* Order element identities in parallel */
	revsort(p, perm, n);

	/* Compute the sample */
	totalmass = 1;
	for (i = 0, n1 = n - 1; i < nans; i++, n1--) {
		rT = totalmass * unif_rand();
		mass = 0;
		for (j = 0; j < n1; j++) {
			mass += p[j];
			if (rT <= mass)
				break;
		}
		ans[i] = perm[j];
		totalmass -= p[j];
		for (k = j; k < n1; k++) {
			p[k] = p[k + 1];
			perm[k] = perm[k + 1];
		}
	}
}

/*
 *  Unequal Probability Sampling.
 *
 *  Modelled after Fortran code provided by:
 *    E. S. Venkatraman <venkat@biosta.mskcc.org>
 *  but with significant modifications in the
 *  "with replacement" case.
 */

/* Unequal probability sampling; with-replacement case */

static void ProbSampleReplace(int n, double *p, int *perm, int nans, int *ans)
{
	double rU;
	int i, j;
	int nm1 = n - 1;

	/* record element identities */
	for (i = 0; i < n; i++)
		perm[i] = i; // + 1;

	/* sort the probabilities into descending order */
	revsort(p, perm, n);

	/* compute cumulative probabilities */
	for (i = 1 ; i < n; i++)
		p[i] += p[i - 1];

	/* compute the sample */
	for (i = 0; i < nans; i++) {
		rU = unif_rand();
		for (j = 0; j < nm1; j++) {
			if (rU <= p[j])
				break;
		}
		ans[i] = perm[j];
	}
}

/* A  version using Walker's alias method, based on Alg 3.13B in
   Ripley (1987).
 */

#define SMALL 10000
static void
walker_ProbSampleReplace(int n, double *p, int *a, int nans, int *ans)
{
	double *q, rU;
	int i, j, k;
	int *H, *L, *HL;

	/* Create the alias tables.
	   The idea is that for HL[0] ... L-1 label the entries with q < 1
	   and L ... H[n-1] label those >= 1.
	   By rounding error we could have q[i] < 1. or > 1. for all entries.
	 */
	//if (n <= SMALL) {
		/* might do this repeatedly, so speed matters */
		// OTTIMIZZATA da me!
	CREAv_i(g_HL, n); // cosi` non ho problemi nel caso sia troppo piccolo
		// OTTIMIZZATA da me!
	CREAv_d(g_q, n); // cosi` non ho problemi nel caso sia troppo piccolo
		//R_CheckStack();
	//~ }
	//~ else {
		/* Slow enough anyway not to risk overflow */
		//~ HL = Calloc(n, int);
		//~ q = Calloc(n, double);
	//~ }
	HL = g_HL->dati;
	q = g_q->dati;
	H = HL - 1;
	L = HL + n;
	for (i = 0; i < n; i++) {
		q[i] = p[i] * n;
		if (q[i] < 1.)
			*++H = i;
		else
			*--L = i;
	}
	if (H >= HL && L < HL + n) { /* So some q[i] are >= 1 and some < 1 */
		for (k = 0; k < n - 1; k++) {
			i = HL[k];
			j = *L;
			a[i] = j;
			q[j] += q[i] - 1;
			if (q[j] < 1.) L++;
			if (L >= HL + n) break; /* now all are >= 1 */
		}
	}
	for (i = 0; i < n; i++)
		q[i] += i;

	/* generate sample */
	for (i = 0; i < nans; i++) {
		rU = unif_rand() * n;
		k = (int) rU;
		ans[i] = (rU < q[k]) ? k + 1 : a[k] + 1;
	}
}

void FixupProb(double *p, int n, int k, Rboolean replace)
{
	double sum;
	int i, npos;
	npos = 0;
	sum = 0.;
	for (i = 0; i < n; i++) {
		if (p[i] > 0) {
			npos++;
			sum += p[i];
		}
	}
	for (i = 0; i < n; i++)
		p[i] /= sum;
}

// modificato il prototipo per aumentare l'efficienza
VETTOREi *sample1(VETTOREi *ris, int n, int k, int replace, double *p)
{
	int *x;
	int i, nc = 0;

	GetRNGstate();
	CREAv_i(ris, k);
	// g_x mi serve comunque (e sempre di dim. n)
	// OTTIMIZZATA da me!
	CREAv_i(g_x, n); // cosi` non ho problemi nel caso sia troppo piccolo
	x = g_x->dati;
	if (p != NULL) {
		FixupProb(p, n, k, replace);
		if (replace) {
			for (i = 0; i < n; i++) {
				if (n * p[i] > 0.1)
					nc++;
			}
			if (nc > 200)
				walker_ProbSampleReplace(n, p, x, k, ris->dati);
			else
				ProbSampleReplace(n, p, x, k, ris->dati);
		}
		else
			ProbSampleNoReplace(n, p, x, k, ris->dati);
	}
	else {
		/* avoid allocation for a single sample */
		if (replace || k < 2)
			SampleReplace(k, n, ris->dati);
		else { // ho gia` allocato g_x
			SampleNoReplace(k, n, ris->dati, x);
		}
	}
	PutRNGstate();
	return ris;
}

#ifdef MDEBUG
	VETTOREi *_sample_p(VETTOREi *ris, const VETTOREi *x, int k, int replace, VETTOREd *p, const char *chi, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi *_sample_p(VETTOREi *ris, const VETTOREi *x, int k, int replace, VETTOREd *p, const char *chi)
#endif
{
#ifndef DET
	int i;
#else
	char buf[MAX_RIGA], nome_mod[256], tmp[256];
	FILE *fp;
#endif

	_Intestazione("\n*** sample_p ***\n");

	assert(x != NULL && p != NULL);
#ifdef FDEBUG
	_StampaVett_i(x);
	fprintf(fp_fdbg, "k = %d, replace = %d\n", k, replace);
	_StampaVett_d(p);
#endif
	CREAv_i(ris, k);
#ifndef DET
	g_y = sample1(g_y, x->dim, k, replace, p->dati);
	for (i = 1; i <= k; i++) {
		ASSEGNAv_i(ris, i, ACCEDIv_i(x, ACCEDIv_i(g_y, i) + 1));
	}
#else
	id++;
	snprintf(nome_mod, 256, "sample_%d.txt", id);
	fp = fopen(nome_mod, "r");
	if (!fp) {
		snprintf(nome_mod, 256, "sample_%d-vuoto.txt", id);
		fp = fopen(nome_mod, "r");
		if (!fp) {
			snprintf(tmp, 256, "Non riesco a leggere il file '%s'", nome_mod);
			error(tmp);
		}
		fprintf(fp_det, "+++Vettore nullo per '%s' da '%s'\n", chi, nome_mod);
		CREAv_i(ris, 0); // imposto a zero la lunghezza
	}
	fgets(buf, MAX_RIGA, fp);
	ris = leggi_seq_i(ris, buf, k);
	fprintf(fp_det, "+++Letti %d valore/i per '%s' da '%s'\n", LENGTHv_i(ris), chi, nome_mod);
	fclose(fp);
#endif
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif

	StrBilanciam();

	return ris;
}

#ifdef MDEBUG
	VETTOREi *_sample(VETTOREi *ris, const VETTOREi *x, int k, int replace, const char *chi, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi *_sample(VETTOREi *ris, const VETTOREi *x, int k, int replace, const char *chi)
#endif
{
#ifndef DET
	int i;
#else
	char buf[MAX_RIGA], nome_mod[256], tmp[256];
	FILE *fp;
#endif

	_Intestazione("\n***sample***\n");

	assert(x != NULL);
#ifdef FDEBUG
	_StampaVett_i(x);
	fprintf(fp_fdbg, "k = %d, replace = %d\n", k, replace);
#endif
	CREAv_i(ris, k);
#ifndef DET
	g_y = sample1(g_y, x->dim, k, replace, NULL);
	for (i = 1; i <= k; i++) {
		ASSEGNAv_i(ris, i, ACCEDIv_i(x, ACCEDIv_i(g_y, i) + 1));
	}
#else
	id++;
	snprintf(nome_mod, 256, "sample_%d.txt", id);
	fp = fopen(nome_mod, "r");
	if (!fp) {
		snprintf(nome_mod, 256, "sample_%d-vuoto.txt", id);
		fp = fopen(nome_mod, "r");
		if (!fp) {
			snprintf(tmp, 256, "Non riesco a leggere il file '%s'", nome_mod);
			error(tmp);
		}
		fprintf(fp_det, "+++Vettore nullo per '%s' da '%s'\n", chi, nome_mod);
		CREAv_i(ris, 0); // imposto a zero la lunghezza
	}
	fgets(buf, MAX_RIGA, fp);
	ris = leggi_seq_i(ris, buf, k);
	fprintf(fp_det, "+++Letti %d valore/i per '%s' da '%s'\n", LENGTHv_i(ris), chi, nome_mod);
	fclose(fp);
#endif
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif

	StrBilanciam();

	return ris;
}
