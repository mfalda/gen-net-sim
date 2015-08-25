#include <stdio.h>

#include "distrib.h"


VETTOREd *leggi_seq_d(VETTOREd *ris, char *buf, int n)
{
	int i = 1;
	char *tmp;
	double ris1;

	CREAv_d(ris, n);
	tmp = strtok(buf, " "); // primo numero
	do {
		sscanf(tmp, "%lg", &ris1);
		tmp = strtok(NULL, " ");
		ASSEGNAv_d(ris, i, ris1);
		i++;
	} while (tmp != NULL && i <= n);
	return ris;
}

VETTOREi *leggi_seq_i(VETTOREi *ris, char *buf, int n)
{
	int i = 1;
	char *tmp;
	int ris1;

	CREAv_i(ris, n);
	tmp = strtok(buf, " "); // primo numero
	do {
		sscanf(tmp, "%d", &ris1);
		tmp = strtok(NULL, " ");
		ASSEGNAv_i(ris, i, ris1);
		i++;
	} while (tmp != NULL && i <= n);
	return ris;
}

VETTOREd *rnorm_s(VETTOREd *ris, int n, double mean, double sd, const char *chi)
{
#ifdef DET
	char buf[MAX_RIGA], nome_mod[256];
	FILE *fp;

	id++;
	snprintf(nome_mod, 256, "normali_%d.txt", id);
	fp = fopen(nome_mod, "r");
	if (!fp) {
		Rprintf("Errore in lettura per il file '%s'!\n", nome_mod);
		error("");
	}
	fgets(buf, MAX_RIGA, fp);
	ris = leggi_seq_d(ris, buf, n);
	fprintf(fp_det, "+++Letti %d valore/i per '%s' da '%s'\n", LENGTHv_d(ris), chi, nome_mod);
	fclose(fp);
#else
	int i;

	for (i = 1; i <= n; i++)
		ASSEGNAv_d(ris, i, rnorm(mean, sd));
#endif

	return ris;
}

VETTOREd *rlnorm_s(VETTOREd *ris, int n, double mean, double sd, const char *chi)
{
#ifdef DET
	char buf[MAX_RIGA], nome_mod[256];
	FILE *fp;

	id++;
	snprintf(nome_mod, 256, "log-normali_%d.txt", id);
	fp = fopen(nome_mod, "r");
	if (!fp) {
		Rprintf("Errore in lettura per il file '%s'!\n", nome_mod);
		error("");
	}
	fgets(buf, MAX_RIGA, fp);
	ris = leggi_seq_d(ris, buf, n);
	fprintf(fp_det, "+++Letti %d valore/i per '%s' da '%s'\n", LENGTHv_d(ris), chi, nome_mod);
	fclose(fp);
#else
	int i;

	for (i = 1; i <= n; i++)
		ASSEGNAv_d(ris, i, rlnorm(mean, sd));
#endif

	return ris;
}

VETTOREd *runif_s(VETTOREd *ris, int n, double a, double b, const char *chi)
{
#ifdef DET
	char buf[MAX_RIGA], nome_mod[256];
	FILE *fp;

	id++;
	snprintf(nome_mod, 256, "unif_%d.txt", id);
	fp = fopen(nome_mod, "r");
	if (!fp) {
		Rprintf("Errore in lettura per il file '%s'!\n", nome_mod);
		error("");
	}
	fgets(buf, MAX_RIGA, fp);
	ris = leggi_seq_d(ris, buf, n);
	fprintf(fp_det, "+++Letti %d valore/i per '%s' da '%s'\n", LENGTHv_d(ris), chi, nome_mod);
	fclose(fp);
#else
	int i;

	for (i = 1; i <= n; i++)
		ASSEGNAv_d(ris, i, runif(a, b));
#endif

	return ris;
}

