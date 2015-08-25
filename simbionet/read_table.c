#include "read_table.h"


LISTA *read_lst_i(LISTA *ris, int n, const char *nomefile)
{
	int i, j, n1, dim_buf;
	FILE *fp;
	char *buf, buf1[2048], *tmp;

	_Intestazione("\n*** read_lst_i ***\n");

	CreaLISTA(ris, NULL, n);
	fp = fopen(nomefile, "r");
	if (!fp) {
		Rprintf("Cannot read file '%s'", nomefile);
		error("");
		return NULL;
	}
	// non so quanto sara` lunga una riga, quindi leggo pezzi da 2048 byte
	dim_buf = 0;
	do {
		fgets(buf1, 2048, fp);
		dim_buf += strlen(buf1);
	} while (strlen(buf1) == 2048 - 1 && buf1[2048 - 2] != '\n');
	// qui mi posso permettere un buffer locale, perche´ la frequenza delle letture non e` elevata
	buf = mia_alloc(dim_buf, char);
	if (dim_buf > 0 && buf == NULL) {
		Rprintf("Not enough memory (read_lst_i # %d, buf)", __LINE__ - 2);
		error("");
	}
	fgets(buf, dim_buf, fp);
	i = 0;
	n1 = 0;
	// calcolo quanti campi ci sono in base all'intestazione, prima di dimensionare la matrice
	while (buf[i] != '\n') {
		if (buf[i] == '\t')
			n1++;
		i++;
	}
	//~ Rprintf("%s: %d campi, lunghezza totale: %d\n", nomefile, n1, dim_buf);
	i = 1;
   while (fgets(buf, dim_buf, fp)) {
		j = 1;
		tmp = strtok(buf, "\t"); // primo numero
		tmp = strtok(NULL, "\t");
		if (tmp == NULL)
			break;
		CtrlSlst(ris, i);
		CREAv_i(ACCEDIlst(ris, i, vi), n1);
      while (tmp != NULL) {
			ASSEGNAv_i(ACCEDIlst(ris, i, vi), j, atoi(tmp));
			tmp = strtok(NULL, "\t");
			j++;
      }
		i++;
	}
	fclose(fp);
	libera(buf);
#ifdef FDEBUG
	_StampaLista(ris);
#endif
	StrBilanciam();

	return ris;
}

// a differenza dell'originale si dimensiona automaticamente!
MATRICEi *read_m_i(MATRICEi *ris, const char *nomefile)
{
	int i, j, l, r, c, n, dim_buf = 0;
	FILE *fp;
	VETTOREi *tmp_i = NULL;
	char *buf, *tmp, *v, *v1;
	const int DIM_BUF = 10240;
	char tmp_b[DIM_BUF];

	_Intestazione("\n*** read_m_i ***\n");

	fp = fopen(nomefile, "r");
	if (!fp) {
		Rprintf("Cannot read file '%s'", nomefile);
		error("");
		return NULL;
	}
	// cerco la dimensione di una riga e il n. di campi
	n = 1;
	do {
		if (fgets(tmp_b, DIM_BUF, fp) == NULL)
			break;
		l = strlen(tmp_b);
		dim_buf += l;
		for (i = 0; i < l; i++) {
			if (tmp_b[i] == '\t')
				n++;
		}
		v = strchr(tmp_b, '\n');
		v1 = strchr(tmp_b, '\r');
	} while (v == NULL && v1 == NULL);

	//~ printf("La dimensione di una riga e` %d\n", dim_buf);
	//~ printf("Ci sono %d campi\n", n);

	// qui mi posso permettere un buffer locale, perche´ la frequenza delle letture non e` elevata
	buf = mia_alloc(dim_buf, char);
	if (dim_buf > 0 && buf == NULL) {
		Rprintf("Not enough memory (read_m_i # %d, buf)", __LINE__ - 2);
		error("");
	}
	CREAv_i(tmp_i, 0);
	i = 1;
	j = 0;
   while (fgets(buf, dim_buf, fp)) {
		tmp = strtok(buf, "\t"); // primo numero
		tmp_i = accoda1_vs_i(tmp_i, atoi(tmp));
		i++;
		tmp = strtok(NULL, "\t");
      while (tmp != NULL) {
			tmp_i = accoda1_vs_i(tmp_i, atoi(tmp));
			tmp = strtok(NULL, "\t");
			i++;
      }
		j++;
	}
	//~ printf("Ci sono %d righe\n", j);
	CREAm_i(ris, j, n);
	i = 1;
	for (r = 1; r <= j; r++) {
		for (c = 1; c <= n; c++) {
			ASSEGNAm_i(ris, r, c, ACCEDIv_i(tmp_i, i));
			i++;
		}
	}
	CANCELLAv_i(tmp_i);
	fclose(fp);
	libera(buf);
#ifdef FDEBUG
	_StampaMatr_i(ris);
#endif
	StrBilanciam();

	return ris;
}

MATRICEd *read_m_d(MATRICEd *ris, int n1, const char *nomefile)
{
	int i, j, n, dim_buf;
	FILE *fp;
	char *buf, *tmp;

	_Intestazione("\n*** read_m_d ***\n");

	fp = fopen(nomefile, "r");
	if (!fp) {
		Rprintf("Cannot read file '%s'", nomefile);
		error("");
		return NULL;
	}
	dim_buf = 20 * n1; // 20 considera la lunghezza di un float piu` il tab (circa)
	// qui mi posso permettere un buffer locale, perche´ la frequenza delle letture non e` elevata
	buf = mia_alloc(dim_buf, char);
	if (dim_buf > 0 && buf == NULL) {
		Rprintf("Not enough memory (read_m_d # %d, buf)", __LINE__ - 2);
		error("");
	}
	i = 0;
	n = 0;
	fgets(buf, dim_buf, fp);
	// calcolo quanti campi ci sono in base all'intestazione, prima di dimensionare la matrice
	while (buf[i] != '\n') {
		if (buf[i] == '\t')
			n++;
		i++;
	}
	CREAm_d(ris, n + 1, n + 1);
	i = 1;
   while (fgets(buf, dim_buf, fp)) {
		j = 1;
		tmp = strtok(buf, "\t"); // primo numero
		tmp = strtok(NULL, "\t");
      while (tmp != NULL) {
			ASSEGNAm_d(ris, i, j, atof(tmp));
			tmp = strtok(NULL, "\t");
			j++;
      }
		i++;
	}
	fclose(fp);
	libera(buf);
#ifdef FDEBUG
	_StampaMatr_d(ris);
#endif
	StrBilanciam();

	return ris;
}

VETTOREd *read_vn_d(VETTOREd *ris, const char *nomefile)
{
	int i;
	float tmp;
	char buf[MAX_BUF];
	FILE *fp;

	_Intestazione("\n*** read_vn_d ***\n");

	fp = fopen(nomefile, "r");
	if (!fp) {
		Rprintf("Cannot read file '%s'", nomefile);
		error("");
		return NULL;
	}
	fgets(buf, MAX_BUF, fp);
	i = 1;
	while (!feof(fp)) {
		fscanf(fp, "%*d\t%e\n", &tmp);
		ASSEGNAv_d(ris, i, tmp);
		i++;
	}
	fclose(fp);

#ifdef FDEBUG
	_StampaVett_d(ris);
#endif
	StrBilanciam();

	return ris;
}
