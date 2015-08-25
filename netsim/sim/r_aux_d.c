#define VERSIONE_d

#include "r_aux_d.h"

#ifndef NDEBUG
void _stampa_lista_d(const char *pref, GList *lista)
{
	GList *li;
	Allocazione *t;

	fprintf(fp_mdbg_d, "%s", pref);
	for(li = lista; li != NULL; li = li->next) {
		t = (Allocazione *) li->data;
		fprintf(fp_mdbg_d, "%s::%s\t", t->file_da->str, t->nome->str);
	}
	fprintf(fp_mdbg_d, "\n\n");
}
#endif

void _InitDbg_d(bool stdout1)
{
#ifdef MDEBUG
	fpm_v_d = fopen("memoria_v_d.csv", "a");
	rewind(fpm_v_d);
	fprintf(fpm_v_d, "\"ID\";\"nome\";\"da\";\"a\";\"dim.massima\";\"riallocazioni\";\"prima linea\";\"ultima linea\"\n");
	allocVett_d = NULL;
	fpm_m_d = fopen("memoria_m_d.csv", "a");
	rewind(fpm_m_d);
	fprintf(fpm_m_d, "\"ID\";\"nome\";\"da\";\"a\";\"dim.max.righe\";\"dim.max.col.\";\"riallocazioni\";\"prima linea\";\"ultima linea\"\n");
	allocMatr_d = NULL;
	if (stdout1)
		fp_mdbg_d = stdout;
	else {
		fp_mdbg_d = fopen("memoria_d.txt", "a");
		rewind(fp_mdbg_d);
	}
#endif
#ifdef FDEBUG
	if (!fp_fdbg) {
		if (stdout1)
			fp_fdbg = stdout;
		else {
			fp_fdbg = fopen("passi.txt", "a");
			rewind(fp_fdbg);
		}
	}
	return;
#endif
}

void _Fflush_d()
{
#ifdef MDEBUG
	fflush(fp_mdbg_d);
#endif
#ifdef FDEBUG
	fflush(fp_fdbg);
	return;
#endif
}

void _StampaRawVett_d(const VETTOREd *v)
{
#ifdef FDEBUG
	int i;

	if (v == NULL) {
		return;
	}
	for (i = 1; i <= LENGTHv_d(v); i++) {
		if (!ISNA(v->dati[i - 1]))
			fprintf(fp_det, " %.16g", v->dati[i - 1]);
		else
			fprintf(fp_det, " NA");
	}
	fprintf(fp_det, "\n");
#endif
}

void _StampaVett_d(const VETTOREd *v)
{
#ifdef FDEBUG
	int i;
	Allocazione *t;

	if (v == NULL) {
		fprintf(fp_fdbg, "vettore nullo\n");
		return;
	}
	t = (Allocazione *) v->mem->data;
	fprintf(fp_fdbg, "%s (%d : %d): [", t->nome->str, v->dim, v->mia_alloc);
	for (i = 1; i <= LENGTHv_d(v); i++)
		fprintf(fp_fdbg, " %.16g", v->dati[i - 1]);
	fprintf(fp_fdbg, " ]\n");
#endif
}

void _StampaRawMatr_d(const MATRICEd *m)
{
#ifdef FDEBUG
	int r, c;

	if (m == NULL) {
		return;
	}
	for (c = 1; c <= LENGTHm2_d(m); c++) {
		for (r = 1; r <= LENGTHm1_d(m); r++) {
			fprintf(fp_det, " %.16g", m->dati[(r - 1) + m->nr * (c - 1)]);
		}
	}
	fprintf(fp_det, "\n");
#endif
}

void _StampaMatr_d(const MATRICEd *m)
{
#ifdef FDEBUG
	int r, c;
	Allocazione *t;

	if (m == NULL) {
		fprintf(fp_fdbg, "matrice nulla\n");
		return;
	}
	t = (Allocazione *) m->mem->data;
	fprintf(fp_fdbg, "%s (%d x %d : %d x %d): [\n", t->nome->str, m->nr, m->nc, m->alloc_r, m->alloc_c);
	for (r = 1; r <= LENGTHm1_d(m); r++) {
		fprintf(fp_fdbg, "\t");
		for (c = 1; c <= LENGTHm2_d(m); c++) {
			fprintf(fp_fdbg, " %.16g", m->dati[(r - 1) + m->nr * (c - 1)]);
		}
		fprintf(fp_fdbg, "\n");
	}
	fprintf(fp_fdbg, " ]\n");
#endif
}

#ifdef MDEBUG
	void  _InitVett_d(const VETTOREd *v, double val, const char *nomefile, int linea)
#else
	void  _InitVett_d(const VETTOREd *v, double val)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: InitVett_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
#endif
	for (i = 1; i <= LENGTHv_d(v); i++)
		_ASSEGNAv_d(v, i, val);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(v);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	void  _InitMatr_d(const MATRICEd *m, double val, const char *nomefile, int linea)
#else
	void  _InitMatr_d(const MATRICEd *m, double val)
#endif
{
	int r, c;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: InitMatr_d\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
#endif
	for (r = 1; r <= LENGTHm1_d(m); r++) {
		for (c = 1; c <= LENGTHm2_d(m); c++)
			_ASSEGNAm_d(m, r, c, val);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	double  _max_s_d(double a, double b, const char *nomefile, int linea)
#else
	double  _max_s_d(double a, double b)
#endif
{
	return (a > b)? a : b;
}

#ifdef MDEBUG
	double  _min_s_d(double a, double b, const char *nomefile, int linea)
#else
	double  _min_s_d(double a, double b)
#endif
{
	return (a < b)? a : b;
}

#ifdef MDEBUG
	VETTOREd * _creaVett_d(VETTOREd *ris, int dim, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _creaVett_d(VETTOREd *ris, int dim)
#endif
{
#ifdef MDEBUG
	Allocazione *mem;
#endif
	VETTOREd *ris1 = NULL;

	if (ris == NULL) {
		ris1 = mia_alloc(1, VETTOREd);
		if (ris1 == NULL) {
			Rprintf("Not enough memory (creaVett_d # %d, ris1)", __LINE__ - 2);
			error("");
		}
		ris1->dim = dim;
		ris1->mia_alloc = dim; // inizialmente mia_alloc == dim
		ris1->dati = NULL;
		ris1->mem = NULL;
		ris1->r = 0;
	}
	else
		ris1 = ris;
#ifdef MDEBUG
	if (ris != NULL && ris->r)
		error("Un vettore proveniente da R non puo` essere alterato! (CREAv_d)");
	if (ris1->mem == NULL) {
		mem = mia_alloc(1, Allocazione);
		if (mem == NULL) {
			Rprintf("Not enough memory (creaVett_d # %d, ris1->mem)", __LINE__ - 2);
			error("");
		}
		fprintf(fp_mdbg_d, "%d_d - allocazione del vettore '%s'[%d] (%p) dalla linea %s # %d\n", g_list_length(allocVett_d) + 1, nome, ris1->dim, ris1, nomefile, linea);
		mem->indir = (size_t) ris1;
		_CREAstr(mem->nome, nome);
		_CREAstr(mem->file_da, nomefile);
		_CREAstr(mem->file_a, "");
		mem->linea_da = linea;
		mem->linea_a = 0;
		mem->max_dim1 = dim;
		mem->rialloc = 0;
		mem->in_uso = 0;
		_CREAstr(mem->file_prima, "");
		mem->prima_linea = 0;
		_CREAstr(mem->file_ultima, "");
		mem->ultima_linea = 0;
		allocVett_d = g_list_append(allocVett_d, mem);
		mem->indx = g_list_length(allocVett_d);
		_stampa_lista_d("allocVett_d+: ", allocVett_d);
		ris1->mem = g_list_last(allocVett_d);
	}
	else {
		mem = (Allocazione *) ris1->mem->data;
		if (mem->linea_da > 0 && mem->linea_a == 0) {
			fprintf(fp_mdbg_d, "%d_d - il vettore '%s' (linea %s # %d) esiste gia`: e` stato allocato alla linea %s # %d\n\n", mem->indx, nome, nomefile, linea, mem->file_da->str, mem->linea_da);
		}
	}
	if (ris1->mia_alloc < dim) {
		fprintf(fp_mdbg_d, "*** il vettore '%s' verra` riallocato passando da %d a %d (%d riallocazione/i, finora)\n\n", nome, ris1->mia_alloc, ris1->mia_alloc * 2 + dim, mem->rialloc + 1);
		mem->rialloc++;
	}
#endif
	CONTROLLA(dim >= 0);
	if (dim > 0) {
		if (ris1->mia_alloc < dim) {
			if (ris1->mia_alloc > 0) {
				libera(ris1->dati);
				ris1->dati = NULL;
			}
			ris1->mia_alloc = ris1->mia_alloc * 2 + dim;
			ris1->dati = mia_alloc((ris1->mia_alloc * 2 + dim), double);
		}
		else if (ris1->dati == NULL)
			ris1->dati = mia_alloc(dim, double);
		if (ris1->dati == NULL) { // non puo` mai essere zero la dimensione richiesta, dato che dim > 0
			Rprintf("Not enough memory (creaVett_d # %d, ris1->dati)", __LINE__ - 2);
			error("");
		}
	}
	ris1->dim = dim;
	return ris1;
}

#ifdef MDEBUG
	MATRICEd * _creaMatr_d(MATRICEd *ris, int nr, int nc, const char *nome, const char *nomefile, int linea)
#else
	MATRICEd * _creaMatr_d(MATRICEd *ris, int nr, int nc)
#endif
{
#ifdef MDEBUG
	Allocazione *mem;
#endif
	MATRICEd *ris1;

	if (ris == NULL) {
		ris1 = mia_alloc(1, MATRICEd);
		if (ris1 == NULL) {
			Rprintf("Not enough memory (creaMatr_d # %d, ris1)", __LINE__ - 2);
			error("");
		}
		ris1->nr = nr;
		ris1->nc = nc;
		ris1->alloc_r = nr; // inizialmente alloc_r == nr
		ris1->alloc_c = nc; // inizialmente alloc_c == nc
		ris1->dati = NULL;
		ris1->mem = NULL;
		ris1->r = 0;
	}
	else {
		ris1 = ris;
	}
#ifdef MDEBUG
	if (ris != NULL && ris->r)
		error("Una matrice proveniente da R non puo` essere alterata! (CREAm_d)");
	if (ris1->mem == NULL) {
		mem = mia_alloc(1, Allocazione);
		if (mem == NULL) {
			Rprintf("Not enough memory (creaMatr_d # %d, ris1->mem)", __LINE__ - 2);
			error("");
		}
		fprintf(fp_mdbg_d, "%d_d - allocazione della matrice '%s'[%d x %d] (%p) dalla linea %s # %d\n", g_list_length(allocMatr_d) + 1, nome, ris1->nr, ris1->nc, ris1, nomefile, linea);
		mem->indir = (size_t) ris1;
		_CREAstr(mem->nome, nome);
		_CREAstr(mem->file_da, nomefile);
		_CREAstr(mem->file_a, "");
		mem->linea_da = linea;
		mem->linea_a = 0;
		mem->max_dim1 = nr;
		mem->max_dim2 = nc;
		mem->rialloc = 0;
		mem->in_uso = 0;
		_CREAstr(mem->file_prima, "");
		mem->prima_linea = 0;
		_CREAstr(mem->file_ultima, "");
		mem->ultima_linea = 0;
		allocMatr_d = g_list_append(allocMatr_d, mem);
		mem->indx = g_list_length(allocMatr_d);
		_stampa_lista_d("allocMatr_d+: ", allocMatr_d);
		ris1->mem = g_list_last(allocMatr_d);
	}
	else {
		mem = (Allocazione *) ris1->mem->data;
		if (mem->linea_da > 0 && mem->linea_a == 0) {
			fprintf(fp_mdbg_d, "%d_d - la matrice '%s' (linea %s # %d) esiste gia`: e` stata allocata alla linea %s # %d\n\n", mem->indx, nome, nomefile, linea, mem->file_da->str, mem->linea_da);
		}
	}
	if (ris1->alloc_r < nr || ris1->alloc_c < nc) {
		fprintf(fp_mdbg_d, "*** la matrice '%s' verra` riallocata passando da %d x %d a %d x %d (%d riallocazione/i, finora)\n\n", nome, ris1->nr, ris1->nc, ris1->nr * 2 + nr, ris1->nc * 2 + nc, mem->rialloc + 1);
		mem->rialloc++;
	}
#endif
	CONTROLLA(nr >= 0 && nc >= 0);
	if (nr > 0 && nc > 0) {
		if (ris1->alloc_r < nr && ris1->alloc_c < nc) {
			if (ris1->alloc_r > 0 && ris1->alloc_c > 0) {
				libera(ris1->dati);
				ris1->dati = NULL;
			}
			ris1->alloc_r = ris1->alloc_r * 2 + nr;
		   ris1->alloc_c = ris1->alloc_c * 2 + nc;
			ris1->dati = mia_alloc(ris1->alloc_r * ris1->alloc_c, double);
		}
		else if (ris1->alloc_r < nr) {
			libera(ris1->dati);
			ris1->dati = NULL;
			ris1->alloc_r = ris1->alloc_r * 2 + nr;
			ris1->dati = mia_alloc(ris1->alloc_r * ris1->alloc_c, double);
		}
		else if (ris1->alloc_c < nc) {
			libera(ris1->dati);
			ris1->dati = NULL;
		   ris1->alloc_c = ris1->alloc_c * 2 + nc;
			ris1->dati = mia_alloc(ris1->alloc_r * ris1->alloc_c, double);
		}
		else if (ris1->dati == NULL)
			ris1->dati = mia_alloc(nr * nc, double);
		if (ris1->dati == NULL) { // non puo` mai essere zero la dimensione richiesta, dato che prima ho controllato la condizione ris1->alloc_r > 0 && ris1->alloc_c > 0
			Rprintf("Not enough memory (creaMatr_d # %d, ris1->dati)", __LINE__ - 2);
			error("");
		}
	}
	ris1->nr = nr;
	ris1->nc = nc;
	return ris1;
}

#ifdef MDEBUG
	int  _lengthVett_d(const VETTOREd *v, const char *nome, const char *nomefile, int linea)
#else
	int  _lengthVett_d(const VETTOREd *v, const char *nome)
#endif
{
#ifdef MDEBUG
	if (v == NULL) {
		Rprintf("il vettore '%s' (linea %s # %d) non e` allocato: non posso ricavarne le dimensioni!\n", nome, nomefile, linea);
		error("");
	}
#endif
	return v->dim;
}

#ifdef MDEBUG
	int  _righeMatr_d(const MATRICEd *m, const char *nome, const char *nomefile, int linea)
#else
	int  _righeMatr_d(const MATRICEd *m, const char *nome)
#endif
{
#ifdef MDEBUG
	if (m == NULL) {
		Rprintf("la matrice '%s' (linea %s # %d) non e` allocata: non posso ricavarne il numero di righe!\n", nome, nomefile, linea);
		error("");
	}
#endif
	return m->nr;
}

#ifdef MDEBUG
	int  _colonneMatr_d(const MATRICEd *m, const char *nome, const char *nomefile, int linea)
#else
	int  _colonneMatr_d(const MATRICEd *m, const char *nome)
#endif
{
#ifdef MDEBUG
	if (m == NULL) {
		Rprintf("la matrice '%s' (linea %s # %d) non e` allocata: non posso ricavarne il numero di righe!\n", nome, nomefile, linea);
		error("");
	}
#endif
	return m->nc;
}

#ifdef MDEBUG
	double  _accediVett_d(const VETTOREd *v, int indx, const char *nome, const char *nomefile, int linea)
#else
	double  _accediVett_d(const VETTOREd *v, int indx, const char *nome)
#endif
{
#ifdef MDEBUG
	Allocazione *mem;

	if (v == NULL || v->mem == NULL) {
		Rprintf("il vettore '%s' (linea %s # %d) non e` allocato: non posso accedere ad un elemento!\n", nome, nomefile, linea);
		error("");
	}
	if (indx < 0 || indx >= v->dim) {
		Rprintf("tentativo di accedere in lettura all'elemento %d del vettore '%s' [%d] (linea %s # %d)!\n", indx + 1, nome, nomefile, linea);
		error("");
	}
	mem = (Allocazione *) v->mem->data;
	if (mem->max_dim1 < indx) {
		mem->max_dim1 = indx + 1;
	}
	if (!mem->in_uso) {
		g_string_assign(mem->file_prima, nomefile);
		mem->prima_linea = linea;
		mem->in_uso = 1;
	}
	else {
		g_string_assign(mem->file_ultima, nomefile);
		mem->ultima_linea = linea;
	}
#endif
	return v->dati[indx];
}

#ifdef MDEBUG
	double  _accediMatr_d(const MATRICEd *m, int r, int c, const char *nome, const char *nomefile, int linea)
#else
	double  _accediMatr_d(const MATRICEd *m, int r, int c, const char *nome)
#endif
{
#ifdef MDEBUG
	Allocazione *mem;

	if (m == NULL || m->mem == NULL) {
		Rprintf("la matrice '%s' (linea %s # %d) non e` allocata: non posso accedere ad un elemento!\n", nome, nomefile, linea);
		error("");
	}
	if (r < 0 || r >= m->nr || c < 0 || c >= m->nc) {
		Rprintf("tentativo di accedere in lettura all'elemento [%d, %d] della matrice '%s' [%d x %d] (linea %s # %d)!\n", r + 1, c + 1, nome, m->nr, m->nc, nomefile, linea);
		error("");
	}
	mem = (Allocazione *) m->mem->data;
	if (mem->max_dim1 < r || mem->max_dim2 < c) {
		mem->max_dim1 = r + 1;
		mem->max_dim2 = c + 1;
	}
	if (!mem->in_uso) {
		g_string_assign(mem->file_prima, nomefile);
		mem->prima_linea = linea;
		mem->in_uso = 1;
	}
	if (!strcmp(mem->file_ultima->str, nomefile)) {
		if (mem->ultima_linea < linea) {
			mem->ultima_linea = linea;
		}
	}
	else {
		g_string_assign(mem->file_ultima, nomefile);
		mem->ultima_linea = linea;
	}
#endif
	return m->dati[r + LENGTHm1_d(m) * c];
}

#ifdef MDEBUG
	double  _accediMVett_d(const MATRICEd *m, int indx, const char *nome, const char *nomefile, int linea)
#else
	double  _accediMVett_d(const MATRICEd *m, int indx, const char *nome)
#endif
{
#ifdef MDEBUG
	Allocazione *mem;

	if (m == NULL || m->mem == NULL) {
		Rprintf("la matrice '%s' (linea %s # %d) non e` allocata: non posso accedere ad un elemento!\n", nome, nomefile, linea);
		error("");
	}
	if (indx < 0 || indx >= m->nr * m->nc) {
		Rprintf("tentativo di accedere in lettura all'elemento %d della matrice '%s' [%d x %d] (linea %s # %d)!\n", indx + 1, nome, m->nr, m->nc, nomefile, linea);
		error("");
	}
	mem = (Allocazione *) m->mem->data;
	if (mem->max_dim1 < indx % m->nr + 1 || mem->max_dim2 < ceil(indx / m->nr) + 1) {
		mem->max_dim1 = indx % m->nr + 1;
		mem->max_dim2 = ceil(indx / m->nr) + 1;
	}
	if (!mem->in_uso) {
		g_string_assign(mem->file_prima, nomefile);
		mem->prima_linea = linea;
		mem->in_uso = 1;
	}
	if (!strcmp(mem->file_ultima->str, nomefile)) {
		if (mem->ultima_linea < linea) {
			mem->ultima_linea = linea;
		}
	}
	else {
		g_string_assign(mem->file_ultima, nomefile);
		mem->ultima_linea = linea;
	}
#endif
	return m->dati[indx];
}

#ifdef MDEBUG
	void  _assegnaVett_d(const VETTOREd *v, int indx, double val, const char *nome, const char *nomefile, int linea)
#else
	void  _assegnaVett_d(const VETTOREd *v, int indx, double val, const char *nome)
#endif
{
#ifdef MDEBUG
	Allocazione *mem;

	if (v == NULL || v->mem == NULL) {
		Rprintf("il vettore '%s' (linea %s # %d) non e` allocato: non posso assegnare un elemento!\n", nome, nomefile, linea);
		error("");
	}
	if (indx < 0 || indx >= v->mia_alloc) {
		Rprintf("tentativo di accedere in scrittura all'elemento %d del vettore '%s' [allocato: %d] (linea %s # %d)!\n", indx + 1, nome, v->mia_alloc, nomefile, linea);
		error("");
	}
	mem = (Allocazione *) v->mem->data;
	if (mem->max_dim1 < indx) {
		mem->max_dim1 = indx + 1;
	}
	if (!mem->in_uso) {
		g_string_assign(mem->file_prima, nomefile);
		mem->prima_linea = linea;
		mem->in_uso = 1;
	}
	if (!strcmp(mem->file_ultima->str, nomefile)) {
		if (mem->ultima_linea < linea) {
			mem->ultima_linea = linea;
		}
	}
	else {
		g_string_assign(mem->file_ultima, nomefile);
		mem->ultima_linea = linea;
	}
#endif
	v->dati[indx] = val;
}

#ifdef MDEBUG
	void  _assegnaMatr_d(const MATRICEd *m, int r, int c, double val, const char *nome, const char *nomefile, int linea)
#else
	void  _assegnaMatr_d(const MATRICEd *m, int r, int c, double val, const char *nome)
#endif
{
#ifdef MDEBUG
	Allocazione *mem;

	if (m == NULL || m->mem == NULL) {
		Rprintf("la matrice '%s' (linea %s # %d) non e` allocata: non posso assegnare un elemento!\n", nome, nomefile, linea);
		error("");
	}
	if (r < 0 || r >= m->alloc_r || c < 0 || c >= m->alloc_c) {
		Rprintf("tentativo di accedere in scrittura all'elemento [%d, %d] della matrice '%s' [allocata: %d x %d] (linea %s # %d)!\n", r + 1, c + 1, nome, m->alloc_r, m->alloc_c, nomefile, linea);
		error("");
	}
	mem = (Allocazione *) m->mem->data;
	if (mem->max_dim1 < r || mem->max_dim2 < c) {
		mem->max_dim1 = r + 1;
		mem->max_dim2 = c + 1;
	}
	if (!mem->in_uso) {
		g_string_assign(mem->file_prima, nomefile);
		mem->prima_linea = linea;
		mem->in_uso = 1;
	}
	if (!strcmp(mem->file_ultima->str, nomefile)) {
		if (mem->ultima_linea < linea) {
			mem->ultima_linea = linea;
		}
	}
	else {
		g_string_assign(mem->file_ultima, nomefile);
		mem->ultima_linea = linea;
	}
#endif
	m->dati[r + LENGTHm1_d(m) * c] = val;
}

#ifdef MDEBUG
	void  _assegnaMVett_d(const MATRICEd *m, int indx, double val, const char *nome, const char *nomefile, int linea)
#else
	void  _assegnaMVett_d(const MATRICEd *m, int indx, double val, const char *nome)
#endif
{
#ifdef MDEBUG
	Allocazione *mem;
	if (m == NULL || m->mem == NULL) {
		Rprintf("la matrice '%s' (linea %s # %d) non e` allocata: non posso assegnare un elemento!\n", nome, nomefile, linea);
		error("");
	}
	if (indx < 0 || indx >= m->alloc_r * m->alloc_c) {
		Rprintf("il vettore degli indici e` piu` lungo degli elementi della matrice: R lo trasformerebbe in vettore, e` questa l'intenzione?\n");
		error("");
	}
	mem = (Allocazione *) m->mem->data;
	if (mem->max_dim1 < indx % m->nr + 1 || mem->max_dim2 < ceil(indx / m->nr) + 1) {
		mem->max_dim1 = indx % m->nr + 1;
		mem->max_dim2 = ceil(indx / m->nr) + 1;
	}
	if (!mem->in_uso) {
		g_string_assign(mem->file_prima, nomefile);
		mem->prima_linea = linea;
		mem->in_uso = 1;
	}
	if (!strcmp(mem->file_ultima->str, nomefile)) {
		if (mem->ultima_linea < linea) {
			mem->ultima_linea = linea;
		}
	}
	else {
		g_string_assign(mem->file_ultima, nomefile);
		mem->ultima_linea = linea;
	}
#endif
	m->dati[indx] = val;
}

#ifdef MDEBUG
	VETTOREd * _cancellaVett_d(VETTOREd *v, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _cancellaVett_d(VETTOREd *v)
#endif
{
#ifdef MDEBUG
	Allocazione *mem;
#endif

	if (v == NULL) {
		return NULL;
	}
#ifdef MDEBUG
	if (v->mem == NULL) {
		Rprintf("\til vettore '%s' esiste, ma la sua struttura per la traccia delle allocazioni no! (linea %s # %d)\n", nome, nomefile, linea);
		error("");
	}
	mem = (Allocazione *) v->mem->data;
	g_string_assign(mem->file_a, nomefile);
	fprintf(fp_mdbg_d, "%d_d - disallocazione del vettore '%s'[%d] (%p) dalla linea %s # %d\n", mem->indx, mem->nome->str, LENGTHv_d(v), v, nomefile, linea);
	mem->linea_a = linea;
	fprintf(fpm_v_d, "\"%d\";\"'%s'\";\"%s # %d\";\"%s # %d\";\"%d\";\"%d\";\"%s # %d\";\"%s # %d\"\n",  mem->indx, mem->nome->str, mem->file_da->str, mem->linea_da, mem->file_a->str, mem->linea_a, mem->max_dim1, mem->rialloc, mem->file_prima->str, mem->prima_linea, mem->file_ultima->str, mem->ultima_linea);
#endif
	if (LENGTHv_d(v) > 0 && !v->r) {
		libera(v->dati);
		v->dati = NULL;
	}
#ifdef MDEBUG
	_CANCELLAstr(mem->nome);
	_CANCELLAstr(mem->file_da);
	_CANCELLAstr(mem->file_a);
	_CANCELLAstr(mem->file_prima);
	_CANCELLAstr(mem->file_ultima);
	allocVett_d = g_list_remove(allocVett_d, mem);
	libera(mem);
	mem = NULL;
	_stampa_lista_d("allocVett_d-: ", allocVett_d);
#endif
	libera(v);
	v = NULL;
	return NULL;
}

#ifdef MDEBUG
	MATRICEd * _cancellaMatr_d(MATRICEd *m, const char *nome, const char *nomefile, int linea)
#else
	MATRICEd * _cancellaMatr_d(MATRICEd *m)
#endif
{
#ifdef MDEBUG
	Allocazione *mem;
#endif

	if (m == NULL) {
		return NULL;
	}
#ifdef MDEBUG
	if (m->mem == NULL) {
		Rprintf("\tla matrice '%s' esiste, ma la sua struttura per la traccia delle allocazioni no! (linea %s # %d)\n", nome, nomefile, linea);
		error("");
	}
	mem = (Allocazione *) m->mem->data;
	g_string_assign(mem->file_a, nomefile);
	fprintf(fp_mdbg_d, "%d_d - disallocazione della matrice '%s'[%d x %d] (%p) dalla linea %s # %d\n", mem->indx, mem->nome->str, LENGTHm1_d(m), LENGTHm2_d(m), m, nomefile, linea);
	mem->linea_a = linea;
	fprintf(fpm_m_d, "\"%d\";\"'%s'\";\"%s # %d\";\"%s # %d\";\"%d\";\"%d\";\"%d\";\"%s # %d\";\"%s # %d\"\n", mem->indx, mem->nome->str, mem->file_da->str, mem->linea_da, mem->file_a->str, mem->linea_a, mem->max_dim1, mem->max_dim2, mem->rialloc, mem->file_prima->str, mem->prima_linea, mem->file_ultima->str, mem->ultima_linea);
#endif
	if (LENGTHm1_d(m) > 0 && LENGTHm2_d(m) > 0 && !m->r) {
		libera(m->dati);
		m->dati = NULL;
	}
#ifdef MDEBUG
	_CANCELLAstr(mem->nome);
	_CANCELLAstr(mem->file_da);
	_CANCELLAstr(mem->file_a);
	_CANCELLAstr(mem->file_prima);
	_CANCELLAstr(mem->file_ultima);
	allocMatr_d = g_list_remove(allocMatr_d, mem);
	libera(mem);
	mem = NULL;
	_stampa_lista_d("allocMatr_d-: ", allocMatr_d);
#endif
	libera(m);
	m = NULL;
	return NULL;
}

void _infoVett_d(const VETTOREd *v)
{
#ifdef MDEBUG
	Allocazione *mem;

	if (v == NULL) {
		Rprintf("vettore non trovato o mai allocato!\n");
		R_FlushConsole();
		return;
	}
	mem = (Allocazione *) v->mem->data;
	Rprintf("Vettore '%s' (%d):\n", mem->nome->str, v);
	Rprintf("\tdimensione: %d\n", v->dim);
	Rprintf("\tdimensione allocata: %d\n", v->mia_alloc);
	Rprintf("\tdato in ingresso?: %d\n", v->r);
	Rprintf("\tallocato da: %s # %d\n", mem->file_da->str, mem->linea_da);
	if (mem->linea_a == 0)
		Rprintf("\tancora allocato\n");
	else
		Rprintf("\tdisallocato da: %s # %d\n", mem->file_a->str, mem->linea_a);
	Rprintf("\tindice progressivo: %d\n", mem->indx);
	Rprintf("\tdimensione massima finora: %d\n", mem->max_dim1);
	Rprintf("\triallocazioni finora: %d\n", mem->rialloc);
	Rprintf("\tprima riga: %s # %d\n", mem->file_prima->str, mem->prima_linea);
	Rprintf("\tultima riga: %s # %d\n\n", mem->file_ultima->str, mem->ultima_linea);
	R_FlushConsole();
#endif
}

void _infoMatr_d(const MATRICEd *m)
{
#ifdef MDEBUG
	Allocazione *mem;

	if (m == NULL) {
		Rprintf("matrice non trovata o mai allocata!\n");
		R_FlushConsole();
		return;
	}
	mem = (Allocazione *) m->mem->data;
	Rprintf("Matrice '%s' (%d):\n", mem->nome->str, m);
	Rprintf("\trighe: %d\n", m->nr);
	Rprintf("\tcolonne: %d\n", m->nc);
	Rprintf("\tdimensione allocata righe: %d\n", m->alloc_r);
	Rprintf("\tdimensione allocata colonne: %d\n", m->alloc_c);
	Rprintf("\tdata in ingresso?: %d\n", m->r);
	Rprintf("\tallocata da: %s # %d\n", mem->file_da->str, mem->linea_da);
	if (mem->linea_a == 0)
		Rprintf("\tancora allocata\n");
	else
		Rprintf("\tdisallocata da: %s # %d\n", mem->file_a->str, mem->linea_a);
	Rprintf("\tindice progressivo: %d\n", mem->indx);
	Rprintf("\tdimensione massima righe finora: %d\n", mem->max_dim1);
	Rprintf("\tdimensione massima colonne finora: %d\n", mem->max_dim2);
	Rprintf("\triallocazioni finora: %d\n", mem->rialloc);
	Rprintf("\tprima riga: %s # %d\n", mem->file_prima->str, mem->prima_linea);
	Rprintf("\tultima riga: %s # %d\n\n", mem->file_ultima->str, mem->ultima_linea);
	R_FlushConsole();
#endif
}

#ifdef MDEBUG
	void  _controllaCanc_d(const char *nomefile, int linea)
#else
	void  _controllaCanc_d()
#endif
{
#ifdef MDEBUG
	int ok = 1;
	Allocazione *t;

	fprintf(fp_mdbg_d, "--------------------------------------\n");
	fprintf(fpm_v_d, "--------------------------------------\n");
	while (allocVett_d != NULL) {
		t = (Allocazione *) allocVett_d->data;
		Rprintf("Il vettore '%s' di tipo 'double' (allocato alla linea %s # %d, indice %d) non e` stato ancora disallocato: lo faccio adesso.\n", t->nome->str, t->file_da->str, t->linea_da, t->indx);
		t->indir = (size_t) (VETTOREd *) _cancellaVett_d((VETTOREd *) t->indir, t->nome->str, nomefile, linea);
		ok = 0;
		_CANCELLAstr(t->nome);
		_CANCELLAstr(t->file_da);
		_CANCELLAstr(t->file_a);
		_CANCELLAstr(t->file_prima);
		_CANCELLAstr(t->file_ultima);
		libera(t);
		t = NULL;
	}
	g_list_free(allocVett_d);
	fprintf(fpm_v_d, "\n");
	fclose(fpm_v_d);
	fprintf(fpm_m_d, "--------------------------------------\n");
	while (allocMatr_d != NULL) {
		t = (Allocazione *) allocMatr_d->data;
		Rprintf("La matrice '%s' di tipo 'double' (allocata alla linea %s # %d, indice %d) non e` stata ancora disallocata: lo faccio adesso\n", t->nome->str, t->file_da->str, t->linea_da, t->indx);
		t->indir = (size_t) (MATRICEd *) _cancellaMatr_d((MATRICEd *) t->indir, t->nome->str, nomefile, linea);
		ok = 0;
		_CANCELLAstr(t->nome);
		_CANCELLAstr(t->file_da);
		_CANCELLAstr(t->file_a);
		_CANCELLAstr(t->file_ultima);
		_CANCELLAstr(t->file_prima);
		libera(t);
		t = NULL;
	}
	g_list_free(allocMatr_d);
	fprintf(fpm_m_d, "\n");
	fclose(fpm_m_d);
	if (!ok)
		fprintf(fp_mdbg_d, "ATTENZIONE: per il bilanciamento delle allocazioni di elementi di tipo 'double' sono state necessarie disallocazioni automatiche.\n");
	fprintf(fp_mdbg_d, "--------------------------------------\n");
	if (fp_mdbg_d != NULL) {
		fclose(fp_mdbg_d);
		fp_mdbg_d = NULL;
	}
#ifdef FDEBUG
	if (fp_fdbg != NULL) {
		fclose(fp_fdbg);
		fp_fdbg = NULL;
	}
#endif
#endif
	return;
}

#ifdef MDEBUG
	SEXP  _daVETTORE_d(VETTOREd *v, int *nProtected, const char *nomefile, int linea)
#else
	SEXP  _daVETTORE_d(VETTOREd *v, int *nProtected)
#endif
{
	int i;
	double *ris1;
	SEXP ris;
#ifdef MDEBUG
	char *nome_vett;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "Trasformo il vettore ");
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
#endif
	CONTROLLA(v != NULL);
	PROTECT(ris = allocVector(REALSXP, LENGTHv_d(v)));
	(*nProtected)++;
	ris1 = NUMERIC_POINTER(ris);
	for (i = 0; i < LENGTHv_d(v); i++)
		ris1[i] = _ACCEDIv_d(v, i + 1);
#ifdef MDEBUG
	nome_vett = ((Allocazione *) v->mem->data)->nome->str;
	v = _cancellaVett_d(v, nome_vett, nomefile, linea);
#else
	CANCELLAv_d(v);
#endif
	return ris;
}

#ifdef MDEBUG
	SEXP  _daMATRICE_d(MATRICEd *m, int *nProtected, const char *nomefile, int linea)
#else
	SEXP  _daMATRICE_d(MATRICEd *m, int *nProtected)
#endif
{
	int i;
	double *ris1;
	SEXP dim, ris;
#ifdef MDEBUG
	char *nome_matr;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "Trasformo la matrice\n");
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
#endif
	PROTECT(ris = allocMatrix(REALSXP, LENGTHm1_d(m), LENGTHm2_d(m)));
	(*nProtected)++;
	ris1 = NUMERIC_POINTER(ris);
	for (i = 0; i < LENGTHm1_d(m) * LENGTHm2_d(m); i++)
		ris1[i] = _ACCEDImv_d(m, i + 1);
	PROTECT(dim  =  allocVector(INTSXP,  2));
	(*nProtected)++;
#ifdef MDEBUG
	nome_matr = ((Allocazione *) m->mem->data)->nome->str;
	m = _cancellaMatr_d(m, nome_matr, nomefile, linea);
#else
	CANCELLAm_d(m);
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _inVETTORE_d(SEXP s, int *nProtected, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _inVETTORE_d(SEXP s, int *nProtected)
#endif
{
	VETTOREd *v = NULL;
#ifdef MDEBUG
	Allocazione *mem;
#endif

	if (isNull(s)) {
#ifdef FDEBUG
	fprintf(fp_fdbg, "Il vettore '%s' e` nullo\n", nome);
#endif
		return NULL;
	}
	PROTECT(s = AS_NUMERIC(s));
	(*nProtected)++;
#ifdef MDEBUG
	v = _creaVett_d(v, 0, nome, nomefile, linea);
#else
	v = _creaVett_d(v, 0);
#endif
	v->dim = length(s);
	v->mia_alloc = length(s);
	v->dati = NUMERIC_POINTER(s);
	v->r = 1;
#ifdef MDEBUG
	mem = v->mem->data;
	g_string_assign(mem->nome, nome);
	g_string_assign(mem->file_da, nomefile);
	mem->linea_da = linea;
#endif
#ifdef FDEBUG
	fprintf(fp_fdbg, "Ho trasformato il vettore ");
	_StampaVett_d(v);
#endif
	return v;
}

#ifdef MDEBUG
	MATRICEd * _inMATRICE_d(SEXP s, int *nProtected, const char *nome, const char *nomefile, int linea)
#else
	MATRICEd * _inMATRICE_d(SEXP s, int *nProtected)
#endif
{
	MATRICEd *m = NULL;
#ifdef MDEBUG
	Allocazione *mem;
#endif

	if (isNull(s)) {
#ifdef FDEBUG
	fprintf(fp_fdbg, "La matrice '%s' e` nulla\n", nome);
#endif
		return NULL;
	}
	PROTECT(s = AS_NUMERIC(s));
	(*nProtected)++;
#ifdef MDEBUG
	m = _creaMatr_d(m, 0, 0, nome, nomefile, linea);
#else
	m = _creaMatr_d(m, 0, 0);
#endif
	m->nr = nrows(s);
	m->nc = ncols(s);
	m->alloc_r = nrows(s);
	m->alloc_c = ncols(s);
	m->dati = NUMERIC_POINTER(s);
	m->r = 1;
 #ifdef MDEBUG
	mem = m->mem->data;
	g_string_assign(mem->nome, nome);
	g_string_assign(mem->file_da, nomefile);
	mem->linea_da = linea;
#endif
#ifdef FDEBUG
	fprintf(fp_fdbg, "Ho trasformato la matrice ");
	_StampaMatr_d(m);
#endif
	return m;
}

#ifdef MDEBUG
	VETTOREi * _which_m_rowindxne_d(VETTOREi *ris, const MATRICEd *m, int riga, double val, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_m_rowindxne_d(VETTOREi *ris, const MATRICEd *m, int riga, double val)
#endif
{
	int i, j;

	CONTROLLA(m != NULL);
#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_m_rowindxne_d\n", linea);
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "riga = %d, val = %.16g\n", riga, val);
#endif
	CONTROLLA(riga > 0 && riga <= m->nr);
	_CREAv_i(ris, LENGTHm1_d(m) * LENGTHm2_d(m));
	for (i = 1, j = 1; i <= LENGTHm2_d(m); i++) {
		if (DIVERSO(_ACCEDIm_d(m, riga, i), val)) {
			_ASSEGNAv_i(ris, j, i);
			j++;
		}
	}
	ris->dim = j - 1;
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _copia_v_d(VETTOREd *ris, const VETTOREd *da, int st, int end, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _copia_v_d(VETTOREd *ris, const VETTOREd *da, int st, int end)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: copia_v_d\n", linea);
#endif
	CONTROLLA(da != NULL);
#ifdef FDEBUG
	_StampaVett_d(da);
	fprintf(fp_fdbg, "start = %d\n", st);
	fprintf(fp_fdbg, "end = %d\n", end);
#endif
	CONTROLLA(st > 0 && end <= LENGTHv_d(da));
	_CREAv_d(ris, LENGTHv_d(da));
	for (i = st; i <= end; i++)
		_ASSEGNAv_d(ris, i, _ACCEDIv_d(da, i));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _op_ss_seqdiv_d(VETTOREd *ris, int da, int a, double div, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _op_ss_seqdiv_d(VETTOREd *ris, int da, int a, double div)
#endif
{
	int i, n;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: op_ss_seqdiv_d\nda = %d, a = %d, div = %3.3f\n", linea, da, a, div);
#endif
	CONTROLLA(da < a);
	_CREAv_d(ris, a - da + 1);
	if (ISNA(div)) {
#ifdef FDEBUG
		_CREAstr(tmp, "");
		g_string_printf(tmp, "ATTENZIONE (op_ss_seqdiv_d, linea %s # %d): divisione per NA!", nomefile, linea);
		warning(tmp->str);
		fprintf(fp_fdbg, tmp->str);
		_CANCELLAstr(tmp);
#endif
		for (i = 1, n = da; n <= a; n++, i++)
			_ASSEGNAv_d(ris, i, NA_REAL);
	}
	else if (Uguale(div, 0.0)) {
#ifdef FDEBUG
		_CREAstr(tmp, "");
		g_string_printf(tmp, "ATTENZIONE (op_ss_seqdiv_d, linea %s # %d): divisione per zero!\n", nomefile, linea);
		warning(tmp->str);
		fprintf(fp_fdbg, tmp->str);
		_CANCELLAstr(tmp);
#endif
		for (i = 1, n = da; n <= a; n++, i++) {
			if (n > 0)
				_ASSEGNAv_d(ris, i, R_PosInf);
			else
				_ASSEGNAv_d(ris, i, R_NegInf);
		}
	}
	else {
		for (i = 1, n = da; n <= a; n++, i++)
			_ASSEGNAv_d(ris, i, (double) n / div);
	}
	ris->dim = i - 1;
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _complementa_d(VETTOREd *ris, const VETTOREd *v, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _complementa_d(VETTOREd *ris, const VETTOREd *v)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: complementa_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
#endif
	_CREAv_d(ris, LENGTHv_d(v));
	for (i = 1; i <= LENGTHv_d(v); i++)
		_ASSEGNAv_d(ris, i, 1 - _ACCEDIv_d(v, i));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _rep_d(VETTOREd *ris, const VETTOREd *v, int ripetizioni, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _rep_d(VETTOREd *ris, const VETTOREd *v, int ripetizioni)
#endif
{
	int i, j;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: rep_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	fprintf(fp_fdbg, "ripetizioni = %d\n", ripetizioni);
	if (ripetizioni == 0) {
		_CREAstr(tmp, "");
		g_string_printf(tmp, "ATTENZIONE (rep_d, %s # %d): zero ripetizioni!\n", nomefile, linea);
		warning(tmp->str);
		fprintf(fp_fdbg, tmp->str);
		_CANCELLAstr(tmp);
	}
#endif
	_CREAv_d(ris, LENGTHv_d(v) * ripetizioni);
	for (j = 1; j <= LENGTHv_d(v); j++) {
		for (i = 0; i < ripetizioni; i++)
			_ASSEGNAv_d(ris, j + i * LENGTHv_d(v), _ACCEDIv_d(v, j));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _seq_d(VETTOREd *ris, double da, double a, double incremento, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _seq_d(VETTOREd *ris, double da, double a, double incremento)
#endif
{
	int i, tot;
	double n;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: seq_d\nda = %.16g, a = %.16g, incremento = %.16g\n", linea, da, a, incremento);
#endif
	CONTROLLA(da <= a && incremento > 0);
	tot = ceil((a - da + 1) / incremento);
	_CREAv_d(ris, tot);

	for (i = 1, n = da; n <= a && i <= tot; n += incremento, i++)
		_ASSEGNAv_d(ris, i, n);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _vettore2s_d(VETTOREd *ris, double el1, double el2, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _vettore2s_d(VETTOREd *ris, double el1, double el2)
#endif
{
#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: vettore2s_d\nel1 = %.16g, el2 = %.16g\n", linea, el1, el2);
#endif
	_CREAv_d(ris, 2);
	_ASSEGNAv_d(ris, 1, el1);
	_ASSEGNAv_d(ris, 2, el2);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _vettore3s_d(VETTOREd *ris, double el1, double el2, double el3, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _vettore3s_d(VETTOREd *ris, double el1, double el2, double el3)
#endif
{
#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: vettore3s_d\nel1 = %.16g, el2 = %.16g, el3 = %.16g\n", linea, el1, el2, el3);
#endif
	_CREAv_d(ris, 3);
	_ASSEGNAv_d(ris, 1, el1);
	_ASSEGNAv_d(ris, 2, el2);
	_ASSEGNAv_d(ris, 3, el3);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _vettore2v_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _vettore2v_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: vettore2v_d\n", linea);
#endif
	CONTROLLA(v1 != NULL && v2 != NULL);
#ifdef FDEBUG
	_StampaVett_d(v1);
	_StampaVett_d(v2);
#endif
	_CREAv_d(ris, LENGTHv_d(v1) + LENGTHv_d(v2));
	for (i = 1, j = 1; j <= LENGTHv_d(v1); i++, j++)
		_ASSEGNAv_d(ris, i, _ACCEDIv_d(v1, j));
	for (j = 1; j <= LENGTHv_d(v2); i++, j++)
		_ASSEGNAv_d(ris, i, _ACCEDIv_d(v2, j));
	ris->dim = i - 1;
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _vettore3v_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2, const VETTOREd *v3, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _vettore3v_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2, const VETTOREd *v3)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: vettore3v_d\n", linea);
#endif
	CONTROLLA(v1 != NULL && v2 != NULL && v3 != NULL);
#ifdef FDEBUG
	_StampaVett_d(v1);
	_StampaVett_d(v2);
	_StampaVett_d(v3);
#endif
	_CREAv_d(ris, LENGTHv_d(v1) + LENGTHv_d(v2) + LENGTHv_d(v3));
	for (i = 1, j = 1; j <= LENGTHv_d(v1); i++, j++)
		_ASSEGNAv_d(ris, i, _ACCEDIv_d(v1, j));
	for (j = 1; j <= LENGTHv_d(v2); i++, j++)
		_ASSEGNAv_d(ris, i, _ACCEDIv_d(v2, j));
	for (j = 1; j <= LENGTHv_d(v3); i++, j++)
		_ASSEGNAv_d(ris, i, _ACCEDIv_d(v3, j));
	ris->dim = i - 1;
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	double  _max_v_d(const VETTOREd *v, const char *nomefile, int linea)
#else
	double  _max_v_d(const VETTOREd *v)
#endif
{
	int i;
	double max;


#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: max_v_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
#endif
	max = _ACCEDIv_d(v, 1);
	for (i = 2; i <= LENGTHv_d(v); i++) {
		if (_ACCEDIv_d(v, i) > max)
			max = _ACCEDIv_d(v, i);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\nmax = %.16g\n", max);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return max;
}

#ifdef MDEBUG
	VETTOREd * _accoda1_vv_d(VETTOREd *v1, const VETTOREd *v2, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _accoda1_vv_d(VETTOREd *v1, const VETTOREd *v2)
#endif
{
	int i = 1, j;
	VETTOREd *ris = NULL;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: accoda1_vv_d\n", linea);
#endif
	CONTROLLA(v1 != NULL && v2 != NULL);
#ifdef FDEBUG
	_StampaVett_d(v1);
	_StampaVett_d(v2);
#endif
	if (v1->mia_alloc >= v1->dim + v2->dim) {
		for (i = v1->dim + 1, j = 1; j <= LENGTHv_d(v2); i++, j++)
			_ASSEGNAv_d(v1, i, _ACCEDIv_d(v2, j));
		v1->dim += LENGTHv_d(v2);
	}
	else {
#ifdef FDEBUG
		_CREAstr(tmp, "");
		g_string_printf(tmp, "ATTENZIONE (accoda_vv_d, linea %s # %d): ingrandito il vettore da %d a %d!\n", nomefile, linea, LENGTHv_d(v1), LENGTHv_d(v1) + LENGTHv_d(v2));
		warning(tmp->str);
		fprintf(fp_fdbg, tmp->str);
		_CANCELLAstr(tmp);
#endif
		_CREAv_d(ris, 2 * (LENGTHv_d(v1) + LENGTHv_d(v2)));
		ris->dim = LENGTHv_d(v1) + LENGTHv_d(v2);
		for (i = 1, j = 1; j <= LENGTHv_d(v1); i++, j++)
			_ASSEGNAv_d(ris, i, _ACCEDIv_d(v1, j));
		for (j = 1; j <= LENGTHv_d(v2); i++, j++)
			_ASSEGNAv_d(ris, i, _ACCEDIv_d(v2, j));
		_CANCELLAv_d(v1);
		v1 = ris;
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(v1);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return v1;
}

#ifdef MDEBUG
	VETTOREd * _somma_vs_d(VETTOREd *ris, const VETTOREd *v, double s, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _somma_vs_d(VETTOREd *ris, const VETTOREd *v, double s)
#endif
{
	int i;


#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: somma_vs_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	fprintf(fp_fdbg, "s = %.16g\n", s);
#endif
	_CREAv_d(ris, LENGTHv_d(v));
	for (i = 1; i <= LENGTHv_d(v); i++)
		_ASSEGNAv_d(ris, i, _ACCEDIv_d(v, i) + s);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	int  _esiste_v_d(double el, const VETTOREd *v, const char *nomefile, int linea)
#else
	int  _esiste_v_d(double el, const VETTOREd *v)
#endif
{
	int i, indx = 0;


#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: esiste_v_d\nel = %.16g\n", linea, el);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
#endif
	for (i = 1; i <= LENGTHv_d(v); i++) {
		if (UGUALE(el, _ACCEDIv_d(v, i))) {
			indx = i;
			break;
		}
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\nindx = %d\n", indx);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return indx;
}

#ifdef MDEBUG
	void  _elimina_indx_d(VETTOREd *v, int indx, const char *nomefile, int linea)
#else
	void  _elimina_indx_d(VETTOREd *v, int indx)
#endif
{
	int i;


#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: elimina_indx_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	fprintf(fp_fdbg, "indx = %d\n", indx);
#endif
	for (i = indx; i <= LENGTHv_d(v) - 1; i++)
		_ASSEGNAv_d(v, i, _ACCEDIv_d(v, i + 1));
	v->dim--;
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(v);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	VETTOREd * _assegna_v_d(VETTOREd *v, int indx, double val, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _assegna_v_d(VETTOREd *v, int indx, double val)
#endif
{
	int i;
	VETTOREd *ris = NULL;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna_v_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	fprintf(fp_fdbg, "indx = %d, val = %.16g\n", indx, val);
#endif
	CONTROLLA(indx > 0);
	if (indx <= LENGTHv_d(v)) {
		_ASSEGNAv_d(v, indx, val);
#ifdef FDEBUG
		fprintf(fp_fdbg, "\n\n");
		_StampaVett_d(v);
		fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	}
	else {
		_CREAv_d(ris, indx);
#ifdef FDEBUG
		_CREAstr(tmp, "");
		g_string_printf(tmp, "ATTENZIONE (assegnav_d, linea %s # %d): e` stato assegnato un elemento al di fuori dei limiti dell'array; ingrandito il vettore da %d a %d!\n", nomefile, linea, LENGTHv_d(v), indx);
		warning(tmp->str);
		fprintf(fp_fdbg, tmp->str);
		_CANCELLAstr(tmp);
#endif
		for (i = 1; i <= LENGTHv_d(v); i++)
			_ASSEGNAv_d(ris, i, _ACCEDIv_d(v, i));
#ifdef FDEBUG
		if (i < indx) {
			_CREAstr(tmp, "");
			g_string_printf(tmp, "ATTENZIONE (assegnav_d, linea %s # %d): e` stato assegnato un elemento non consecutivo e sara`/nno aggiunto/i %d valore/i NA!\n", nomefile, linea, indx - i);
			warning(tmp->str);
			fprintf(fp_fdbg, tmp->str);
			_CANCELLAstr(tmp);
		}
#endif
		for (; i < indx; i++) {
			_ASSEGNAv_d(ris, i, NA_REAL);
		}
		_ASSEGNAv_d(ris, indx, val);
		_CANCELLAv_d(v);
		ris->dim = indx;
		v = ris;
#ifdef FDEBUG
		fprintf(fp_fdbg, "\n\n");
		_StampaVett_d(ris);
		fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	}
	return v;
}

#ifdef MDEBUG
	VETTOREd * _somma_righe_d(VETTOREd *ris, const MATRICEd *m, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _somma_righe_d(VETTOREd *ris, const MATRICEd *m)
#endif
{
	int r, c;


#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: somma_righe_d\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
#endif
	_CREAv_d(ris, LENGTHm1_d(m));
	for (r = 1; r <= LENGTHm1_d(m); r++) {
		_ASSEGNAv_d(ris, r, 0);
		for (c = 1; c <= LENGTHm2_d(m); c++) {
			_ASSEGNAv_d(ris, r, _ACCEDIv_d(ris, r)  + _ACCEDIm_d(m, r, c));
		}
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	void  _rev1_d(const VETTOREd *v, const char *nomefile, int linea)
#else
	void  _rev1_d(const VETTOREd *v)
#endif
{
	int i;
	double tmp;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: rev_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
#endif
	for (i = 1; i <= LENGTHv_d(v) / 2; i++) {
		tmp = _ACCEDIv_d(v, i);
		_ASSEGNAv_d(v, i, _ACCEDIv_d(v, LENGTHv_d(v) - i + 1));
		_ASSEGNAv_d(v, LENGTHv_d(v) - i + 1, tmp);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(v);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	VETTOREd * _ordina_d(VETTOREd *ris, const VETTOREd *v, bool decr, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _ordina_d(VETTOREd *ris, const VETTOREd *v, bool decr)
#endif
{

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: ordina_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	fprintf(fp_fdbg, "decr = %d\n", decr);
#endif
	_CREAv_d(ris, LENGTHv_d(v));
	ris = copia_v_d(ris, v, 1, LENGTHv_d(v));
	R_rsort(ris->dati, LENGTHv_d(ris));
	if (decr)
		rev1_d(ris);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
		return ris;
}

#ifdef MDEBUG
	void  _ordina1_d(const VETTOREd *v, bool decr, const char *nomefile, int linea)
#else
	void  _ordina1_d(const VETTOREd *v, bool decr)
#endif
{
#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: ordina1_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	fprintf(fp_fdbg, "decr = %d\n", decr);
#endif
	R_rsort(v->dati, LENGTHv_d(v));
	if (decr)
		rev1_d(v);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(v);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	VETTOREd * _riga_d(VETTOREd *ris, const MATRICEd *m, int r, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _riga_d(VETTOREd *ris, const MATRICEd *m, int r)
#endif
{
	int c;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: riga_d\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "r = %d", r);
#endif
	CONTROLLA(r >= 0 && r <= m->nr);
	if (r == 0) {
	#ifdef FDEBUG
		_CREAstr(tmp, "");
		g_string_printf(tmp, "ATTENZIONE (riga_d, linea %s # %d): parametro r = 0!\n\n", nomefile, linea);
		warning(tmp->str);
		fprintf(fp_fdbg, tmp->str);
		_CANCELLAstr(tmp);
#endif
		_CREAv_d(ris, 0);
	}
	else {
		_CREAv_d(ris, LENGTHm2_d(m));
		for (c = 1; c <= LENGTHm2_d(m); c++) {
			_ASSEGNAv_d(ris, c, _ACCEDIm_d(m, r, c));
		}
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _which_v_indxle_d(VETTOREi *ris, const VETTOREd *v, double val, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_v_indxle_d(VETTOREi *ris, const VETTOREd *v, double val)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_v_indxle_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	fprintf(fp_fdbg, "val = %.16g\n", val);
#endif
	_CREAv_i(ris, LENGTHv_d(v));
	for (i = 1, j = 1; i <= LENGTHv_d(v); i++) {
		if (_ACCEDIv_d(v, i) <= val) {
			_ASSEGNAv_i(ris, j, i);
			j++;
		}
	}
	ris->dim = j - 1;
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	double  _min_v_d(const VETTOREd *v, const char *nomefile, int linea)
#else
	double  _min_v_d(const VETTOREd *v)
#endif
{
	int i;
	double min;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: min_v_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
#endif
	min = _ACCEDIv_d(v, 1);
	for (i = 2; i <= LENGTHv_d(v); i++) {
		if (_ACCEDIv_d(v, i) < min)
			min = _ACCEDIv_d(v, i);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\nmin = %.16g\n", min);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return min;
}

#ifdef MDEBUG
	void  _somma1_vs_d(VETTOREd *v, double s, const char *nomefile, int linea)
#else
	void  _somma1_vs_d(VETTOREd *v, double s)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: somma1_vs_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	fprintf(fp_fdbg, "s = %.16g\n", s);
#endif
	for (i = 1; i <= LENGTHv_d(v); i++)
		_ASSEGNAv_d(v, i, _ACCEDIv_d(v, i) + s);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(v);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	VETTOREd * _copia_v_indx_d(VETTOREd *ris, const VETTOREd *v, const VETTOREi *indx, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _copia_v_indx_d(VETTOREd *ris, const VETTOREd *v, const VETTOREi *indx)
#endif
{
	int i, j;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: copia_v_indx_d\n", linea);
#endif
	CONTROLLA(v != NULL && indx != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	_StampaVett_i(indx);
#endif
	_CREAv_d(ris, LENGTHv_i(indx));
	for (i = 1, j = 1; i <= LENGTHv_i(indx); i++) {
		if (_ACCEDIv_i(indx, i) == 0) {
#ifdef FDEBUG
			_CREAstr(tmp, "");
			g_string_printf(tmp, "ATTENZIONE (copia_v_indx_d, linea %s # %d): il vettore degli indici contiene uno 0 alla posizione %d: saltato!\n", nomefile, linea, i);
			warning(tmp->str);
			fprintf(fp_fdbg, tmp->str);
			_CANCELLAstr(tmp);
#endif
			continue;
		}
		_ASSEGNAv_d(ris, j, _ACCEDIv_d(v, _ACCEDIv_i(indx, i)));
		j++;
	}
	ris->dim = j - 1;
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _which_v_indxgt_d(VETTOREi *ris, const VETTOREd *v, double val, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_v_indxgt_d(VETTOREi *ris, const VETTOREd *v, double val)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_v_indxgt_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	fprintf(fp_fdbg, "val = %.16g\n", val);
#endif
	_CREAv_i(ris, LENGTHv_d(v));
	for (i = 1, j = 1; i <= LENGTHv_d(v); i++) {
		if (_ACCEDIv_d(v, i) > val) {
			_ASSEGNAv_i(ris, j, i);
			j++;
		}
	}
	ris->dim = j - 1;
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _which_vv_eq_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _which_vv_eq_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2)
#endif
{
	int i, j = 1;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_vv_eq_d\n", linea);
#endif
	CONTROLLA(v1 != NULL && v2 != NULL);
#ifdef FDEBUG
	_StampaVett_d(v1);
	_StampaVett_d(v2);
#endif
	_CREAv_d(ris, LENGTHv_d(v2));
	for (i = 1; i <= LENGTHv_d(v2); i++) {
		if (esiste_v_d(_ACCEDIv_d(v2, i), v1)) {
			_ASSEGNAv_d(ris, j, _ACCEDIv_d(v2, i));
			j++;
		}
	}
	ris->dim = j - 1;
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	void  _setdiff1_d(VETTOREd *v1, const VETTOREd *v2, const char *nomefile, int linea)
#else
	void  _setdiff1_d(VETTOREd *v1, const VETTOREd *v2)
#endif
{
	int indx, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: setdiff1_d\n", linea);
#endif
	CONTROLLA(v1 != NULL && v2 != NULL);
#ifdef FDEBUG
	_StampaVett_d(v1);
	_StampaVett_d(v2);
#endif
	for (j = 1; j <= LENGTHv_d(v2); j++) {
		indx = esiste_v_d(_ACCEDIv_d(v2, j), v1);
		if (indx > 0)
			elimina_indx_d(v1, indx);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(v1);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	VETTOREi * _which_m_rowindxeq_d(VETTOREi *ris, const MATRICEd *m, int riga, double val, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_m_rowindxeq_d(VETTOREi *ris, const MATRICEd *m, int riga, double val)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_m_rowindxeq_d\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "riga = %d, val = %.16g\n", riga, val);
#endif
	CONTROLLA(riga > 0 && riga <= LENGTHm1_d(m));
	_CREAv_i(ris, LENGTHm1_d(m) * LENGTHm2_d(m));
	for (i = 1, j = 1; i <= LENGTHm2_d(m); i++) {
		if (UGUALE(_ACCEDIm_d(m, riga, i), val)) {
			_ASSEGNAv_i(ris, j, i);
			j++;
		}
	}
	ris->dim = j - 1;
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	void  _assegna1_vs_indx_d(const VETTOREd *v, const VETTOREi *indx, double val, const char *nomefile, int linea)
#else
	void  _assegna1_vs_indx_d(const VETTOREd *v, const VETTOREi *indx, double val)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_vs_indx_d\n", linea);
#endif
	CONTROLLA(v != NULL && indx != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	_StampaVett_i(indx);
	fprintf(fp_fdbg, "val = %.16g\n", val);
#endif
	for (i = 1; i <= LENGTHv_i(indx); i++)
		_ASSEGNAv_d(v, _ACCEDIv_i(indx, i), val);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(v);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	VETTOREd * _assegna_v_indx_d(VETTOREd *ris, const VETTOREd *v, const VETTOREi *indx, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _assegna_v_indx_d(VETTOREd *ris, const VETTOREd *v, const VETTOREi *indx)
#endif
{
	int i, k = 1;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna_v_indx_d\n", linea);
#endif
	CONTROLLA(v != NULL && indx != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	_StampaVett_i(indx);
#endif
	_CREAv_d(ris, LENGTHv_i(indx));
	for (i = 1; i <= LENGTHv_i(indx); i++) {
		_ASSEGNAv_d(ris, k, _ACCEDIv_d(v, _ACCEDIv_i(indx, i)));
		k++;
	}
	ris->dim = k - 1;
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _assegna_v_indxNA_d(VETTOREd *ris, const VETTOREd *v, const VETTOREi *indx, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _assegna_v_indxNA_d(VETTOREd *ris, const VETTOREd *v, const VETTOREi *indx)
#endif
{
	int i, k = 1, pos;
	double tmp1;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna_v_indxNA_d\n", linea);
#endif
	CONTROLLA(v != NULL && indx != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	_StampaVett_i(indx);
#endif
	_CREAv_d(ris, LENGTHv_i(indx));
	for (i = 1; i <= LENGTHv_i(indx); i++) {
		pos = _ACCEDIv_i(indx, i);
		if (pos == NA_INTEGER) {
#ifdef FDEBUG
			_CREAstr(tmp, "");
			g_string_printf(tmp, "ATTENZIONE (assegna_vindxNA_d, linea %s # %d): saltato l'elemento specificato alla posizione %d dell'array degli indici perch� NA!\n", nomefile, linea, i);
			warning(tmp->str);
			fprintf(fp_fdbg, tmp->str);
			_CANCELLAstr(tmp);
#endif
			continue;
		}
		else {
			if (pos < 1 || pos > LENGTHv_d(v)) {
#ifdef FDEBUG
				_CREAstr(tmp, "");
				g_string_printf(tmp, "ATTENZIONE (assegna_vindxNA_d, linea %s # %d): assegnato NA per via di un elemento specificato alla posizione %d che e` al di fuori dei limiti dell'array (%d)!\n", nomefile, linea, _ACCEDIv_i(indx, i), LENGTHv_d(v));
				warning(tmp->str);
				fprintf(fp_fdbg, tmp->str);
				_CANCELLAstr(tmp);
#endif
				tmp1 = NA_REAL;
			}
			else
				tmp1 = _ACCEDIv_d(v, _ACCEDIv_i(indx, i));
			_ASSEGNAv_d(ris, k, tmp1);
			k++;
		}
	}
	ris->dim = k - 1;
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	double  _somma_v_indx_d(const VETTOREd *v, const VETTOREi *indx, const char *nomefile, int linea)
#else
	double  _somma_v_indx_d(const VETTOREd *v, const VETTOREi *indx)
#endif
{
	int i;
	double ris = 0;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: somma_v_indx_d\n", linea);
#endif
	CONTROLLA(v != NULL && indx != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	_StampaVett_i(indx);
#endif
	for (i = 1; i <= LENGTHv_i(indx); i++) {
		if (_ACCEDIv_i(indx, i) >= LENGTHv_d(v))
			return NA_REAL;
		ris += _ACCEDIv_d(v, _ACCEDIv_i(indx, i));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\nsomma = %.16g\n", ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	void  _dividi1_vs_d(VETTOREd *v, double div, const char *nomefile, int linea)
#else
	void  _dividi1_vs_d(VETTOREd *v, double div)
#endif
{
	int i;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: dividi1_vs_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	fprintf(fp_fdbg, "div = %3.3f\n", div);
#endif
	if (ISNA(div)) {
#ifdef FDEBUG
	_CREAstr(tmp, "");
	g_string_printf(tmp, "ATTENZIONE (dividi1_vs_d, linea %s # %d): divisione per NA!", nomefile, linea);
	warning(tmp->str);
	fprintf(fp_fdbg, tmp->str);
	_CANCELLAstr(tmp);
#endif
		for (i = 1; i <= LENGTHv_d(v); i++)
			_ASSEGNAv_d(v, i, NA_REAL);
	}
	else if (Uguale(div, 0.0)) {
#ifdef FDEBUG
	_CREAstr(tmp, "");
	g_string_printf(tmp, "ATTENZIONE (dividi1_vs_d, linea %s # %d): divisione per zero!\n", nomefile, linea);
	warning(tmp->str);
	fprintf(fp_fdbg, tmp->str);
	_CANCELLAstr(tmp);
#endif
		for (i = 1; i <= LENGTHv_d(v); i++) {
			if (_ACCEDIv_d(v, i) > 0)
				_ASSEGNAv_d(v, i, R_PosInf);
			else
				_ASSEGNAv_d(v, i, R_NegInf);
		}
	}
	else {
		for (i = 1; i <= LENGTHv_d(v); i++)
			_ASSEGNAv_d(v, i, _ACCEDIv_d(v, i) / div);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(v);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	void  _assegna1_m_vv_d(const MATRICEd *m, const VETTOREi *vr, const VETTOREi *vc, double val, const char *nomefile, int linea)
#else
	void  _assegna1_m_vv_d(const MATRICEd *m, const VETTOREi *vr, const VETTOREi *vc, double val)
#endif
{
	int r, c;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_m_vv_d\n", linea);
#endif
	CONTROLLA(m != NULL && vr != NULL && vc != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	_StampaVett_i(vr);
	_StampaVett_i(vc);
	fprintf(fp_fdbg, "val = %.16g\n", val);
#endif
	for (r = 1; r <= LENGTHv_i(vr); r++) {
		if (_ACCEDIv_i(vr, r) == 0) {
#ifdef FDEBUG
			_CREAstr(tmp, "");
			g_string_printf(tmp, "ATTENZIONE (assegna1_m_vv_d, linea %s # %d): il vettore delle righe contiene uno 0 alla posizione %d: saltato!\n", nomefile, linea, r);
			warning(tmp->str);
			fprintf(fp_fdbg, tmp->str);
			_CANCELLAstr(tmp);
#endif
			continue;
		}
		for (c = 1; c <= LENGTHv_i(vc); c++) {
			if (_ACCEDIv_i(vc, c) == 0) {
#ifdef FDEBUG
				_CREAstr(tmp, "");
				g_string_printf(tmp, "ATTENZIONE (assegna1_m_vv_d, linea %s # %d): il vettore delle colonne contiene uno 0 alla posizione %d: saltato!\n", nomefile, linea, c);
				warning(tmp->str);
				fprintf(fp_fdbg, tmp->str);
				_CANCELLAstr(tmp);
#endif
				continue;
			}
			_ASSEGNAm_d(m, _ACCEDIv_i(vr, r), _ACCEDIv_i(vc, c), val);
		}
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	MATRICEd * _abs_m_d(MATRICEd *ris, const MATRICEd *m, const char *nome, const char *nomefile, int linea)
#else
	MATRICEd * _abs_m_d(MATRICEd *ris, const MATRICEd *m)
#endif
{
	int r, c;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: abs_m_d\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
#endif
	_CREAm_d(ris, LENGTHm1_d(m), LENGTHm2_d(m));

	for (r = 1; r <= LENGTHm1_d(m); r++) {
		for (c = 1; c <= LENGTHm2_d(m); c++)
			_ASSEGNAm_d(ris, r, c, fabs(_ACCEDIm_d(m, r, c)));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _which_v_indxne_d(VETTOREi *ris, const VETTOREd *v, double val, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_v_indxne_d(VETTOREi *ris, const VETTOREd *v, double val)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_v_indxne_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	fprintf(fp_fdbg, "val = %.16g\n", val);
#endif
	_CREAv_i(ris, LENGTHv_d(v));
	for (i = 1, j = 1; i <= LENGTHv_d(v); i++) {
		if (DIVERSO(_ACCEDIv_d(v, i), val)) {
			_ASSEGNAv_i(ris, j, i);
			j++;
		}
	}
	ris->dim = j - 1;
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _which_m_colindxeq_d(VETTOREi *ris, const MATRICEd *m, int colonna, double val, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_m_colindxeq_d(VETTOREi *ris, const MATRICEd *m, int colonna, double val)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_m_colindxeq_d\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "colonna = %d, val = %.16g\n", colonna, val);
#endif
	CONTROLLA(colonna > 0 && colonna <= LENGTHm2_d(m));
	_CREAv_i(ris, LENGTHm1_d(m) * LENGTHm2_d(m));
	for (i = 1, j = 1; i <= LENGTHm1_d(m); i++) {
		if (UGUALE(_ACCEDIm_d(m, i, colonna), val)) {
			_ASSEGNAv_i(ris, j, i);
			j++;
		}
	}
	ris->dim = j - 1;
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _interseca_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _interseca_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2)
#endif
{
	int i, l, j = 1, indx, indx1;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: interseca_d\n", linea);
#endif
	CONTROLLA(v1 != NULL && v2 != NULL);
#ifdef FDEBUG
	_StampaVett_d(v1);
	_StampaVett_d(v2);
#endif
	l = min_s_i(LENGTHv_d(v1), LENGTHv_d(v2));
	_CREAv_d(ris, l);
	ris->dim = 0;
	for (i = 1; i <= LENGTHv_d(v1); i++) {
		indx = esiste_v_d(_ACCEDIv_d(v1, i), v2);
		indx1 = esiste_v_d(_ACCEDIv_d(v1, i), ris);
		if (indx > 0 && indx1 <= 0) {
			_ASSEGNAv_d(ris, j, _ACCEDIv_d(v1, i));
			ris->dim = j;
			j++;
		}
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	void  _copia1_m_riga_d(const MATRICEd *m1, int riga1, const MATRICEd *m2, int riga2, const char *nomefile, int linea)
#else
	void  _copia1_m_riga_d(const MATRICEd *m1, int riga1, const MATRICEd *m2, int riga2)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: copia1_m_riga_d\n", linea);
#endif
	CONTROLLA(m1 != NULL && m2 != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m1);
	_StampaMatr_d(m2);
	fprintf(fp_fdbg, "riga1 = %d, riga2 = %d\n", riga1, riga2);
#endif
	CONTROLLA(riga1 > 0 && riga1 <= m1->nr && riga2 > 0 && riga2 <= m2->nr);
	for (i = 1; i <= LENGTHm2_d(m1); i++)
		_ASSEGNAm_d(m1, riga1, i, _ACCEDIm_d(m2, riga2, i));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(m1);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return;
}

#ifdef MDEBUG
	MATRICEd * _aggiungi_riga_d(MATRICEd *m1, int riga1, const MATRICEd *m2, int riga2, const char *nome, const char *nomefile, int linea)
#else
	MATRICEd * _aggiungi_riga_d(MATRICEd *m1, int riga1, const MATRICEd *m2, int riga2)
#endif
{
	int i, r, c;
	MATRICEd *tmp_m = NULL;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: aggiungi_riga_d\n", linea);
#endif
	CONTROLLA(m1 != NULL && m2 != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m1);
	_StampaMatr_d(m2);
	fprintf(fp_fdbg, "riga1 = %d, riga2 = %d\n", riga1, riga2);
#endif
	CONTROLLA(riga2 > 0 && riga2 <= m2->nr);
	if (riga1 > m2->nr + 1) {
		Rprintf("la riga da aggiungere non e` consecutiva: R riempirebbe quelle in mezzo con NA, e` questa l'intenzione?\n");
		error("");
	}
	if (riga1 > m1->alloc_r) {
#ifdef FDEBUG
		_CREAstr(tmp, "");
		g_string_printf(tmp, "ATTENZIONE (aggiungi_riga_d, linea %s # %d): ingrandite le righe della matrice da %d a %d!\n", nomefile, linea, LENGTHm1_d(m1), riga1);
		warning(tmp->str);
		fprintf(fp_fdbg, tmp->str);
		_CANCELLAstr(tmp);
#endif
		_CREAm_d(tmp_m, 2 * riga1, m1->nc);
		tmp_m->nr = riga1;
		for (r = 1; r <= LENGTHm1_d(m1); r++) {
			for (c = 1; c <= LENGTHm2_d(m1); c++)
				_ASSEGNAm_d(tmp_m, r, c, _ACCEDIm_d(m1, r, c));
		}
		_CANCELLAm_d(m1);
		m1 = tmp_m;
	}
	for (i = 1; i <= LENGTHm2_d(m1); i++)
		_ASSEGNAm_d(m1, riga1, i, _ACCEDIm_d(m2, riga2, i));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(m1);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return m1;
}

#ifdef MDEBUG
	void  _assegna1_ms_rigaindx_d(const MATRICEd *m, int riga, const VETTOREi *indx, double val, const char *nomefile, int linea)
#else
	void  _assegna1_ms_rigaindx_d(const MATRICEd *m, int riga, const VETTOREi *indx, double val)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_ms_rigaindx_d\n", linea);
#endif
	CONTROLLA(m != NULL && indx != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "riga = %d\n", riga);
	_StampaVett_i(indx);
	fprintf(fp_fdbg, "val = %.16g\n", val);
#endif
	CONTROLLA(riga > 0 && indx != NULL);
	for (i = 1; i <= LENGTHv_i(indx); i++) {
		_ASSEGNAm_d(m, riga, _ACCEDIv_i(indx, i), val);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return;
}

#ifdef MDEBUG
	MATRICEd * _aggiungi_ms_rigaindx_d(MATRICEd *m, int riga, const VETTOREi *indx, double val, const char *nome, const char *nomefile, int linea)
#else
	MATRICEd * _aggiungi_ms_rigaindx_d(MATRICEd *m, int riga, const VETTOREi *indx, double val)
#endif
{
	int i, r, c;
	MATRICEd *tmp_m = NULL;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: aggiungi_ms_rigaindx_d\n", linea);
#endif
	CONTROLLA(m != NULL && indx != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "riga = %d\n", riga);
	_StampaVett_i(indx);
	fprintf(fp_fdbg, "val = %.16g\n", val);
#endif
	if (riga > m->nr + 1) {
		Rprintf("la riga da aggiungere non e` minore o uguale al numero di righe attuali + 1: R riempirebbe quelle in mezzo con NA, e` questa l'intenzione?\n");
		error("");
	}
	if (riga > m->alloc_r) {
#ifdef FDEBUG
		_CREAstr(tmp, "");
		g_string_printf(tmp, "ATTENZIONE (aggiungi_riga_indx_d, linea %s # %d): ingrandite le righe della matrice da %d a %d!\n", nomefile, linea, LENGTHm1_d(m), riga);
		warning(tmp->str);
		fprintf(fp_fdbg, tmp->str);
		_CANCELLAstr(tmp);
#endif
		_CREAm_d(tmp_m, 2 * riga, m->nc);
		tmp_m->nr = riga;
		for (r = 1; r <= LENGTHm1_d(m); r++) {
			for (c = 1; c <= LENGTHm2_d(m); c++)
				_ASSEGNAm_d(tmp_m, r, c, _ACCEDIm_d(m, r, c));
		}
		_CANCELLAm_d(m);
		m = tmp_m;
	}
	for (i = 1; i <= LENGTHv_i(indx); i++) {
		_ASSEGNAm_d(m, riga, _ACCEDIv_i(indx, i), val);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return m;
}

#ifdef MDEBUG
	void  _assegna1_mm_rigaindx_d(const MATRICEd *m1, int riga, const MATRICEd *m2, const VETTOREi *indx, const char *nomefile, int linea)
#else
	void  _assegna1_mm_rigaindx_d(const MATRICEd *m1, int riga, const MATRICEd *m2, const VETTOREi *indx)
#endif
{
	int i;
	double tmp;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_mm_rigaindx_d\n", linea);
#endif
	CONTROLLA(m1 != NULL && m2 != NULL && indx != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m1);
	fprintf(fp_fdbg, "riga = %d\n", riga);
	_StampaMatr_d(m2);
	_StampaVett_i(indx);
#endif
	CONTROLLA(riga > 0 && riga <= m1->nr && indx->dim == m1->nc);
	for (i = 1; i <= LENGTHv_i(indx); i++) {
		if (_ACCEDIv_i(indx, i) > LENGTHm1_d(m2) * LENGTHm2_d(m2))
			tmp = NA_REAL;
		else
			tmp = m2->dati[_ACCEDIv_i(indx, i) - 1];
		_ASSEGNAm_d(m1, riga, i, tmp);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(m1);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return;
}

#ifdef MDEBUG
	void  _assegna1_ms_riga_d(const MATRICEd *m, int riga, double val, const char *nomefile, int linea)
#else
	void  _assegna1_ms_riga_d(const MATRICEd *m, int riga, double val)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_ms_riga_d\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "riga = %d\n", riga);
	fprintf(fp_fdbg, "val = %.16g\n", val);
#endif
	CONTROLLA(riga > 0 && riga <= m->nr);
	for (i = 1; i <= LENGTHm2_d(m); i++)
		_ASSEGNAm_d(m, riga, i, val);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return;
}

#ifdef MDEBUG
	void  _assegna1_m_indxlt_d(const MATRICEd *m, double val1, double val2, const char *nomefile, int linea)
#else
	void  _assegna1_m_indxlt_d(const MATRICEd *m, double val1, double val2)
#endif
{
	int r, c;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_mindxlt_d\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "val1 = %.16g\n", val1);
	fprintf(fp_fdbg, "val2 = %.16g\n", val2);
#endif
	for (r = 1; r <= LENGTHm1_d(m); r++) {
		for (c = 1; c <= LENGTHm2_d(m); c++) {
			if (_ACCEDIm_d(m, r, c) < val1)
				_ASSEGNAm_d(m, r, c, val2);
		}
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return;
}

#ifdef MDEBUG
	void  _assegna1_m_indxgt_d(const MATRICEd *m, double val1, double val2, const char *nomefile, int linea)
#else
	void  _assegna1_m_indxgt_d(const MATRICEd *m, double val1, double val2)
#endif
{
	int r, c;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_m_indxgt_d\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "val1 = %.16g\n", val1);
	fprintf(fp_fdbg, "val2 = %.16g\n", val2);
#endif
	for (r = 1; r <= LENGTHm1_d(m); r++) {
		for (c = 1; c <= LENGTHm2_d(m); c++) {
			if (_ACCEDIm_d(m, r, c) > val1)
				_ASSEGNAm_d(m, r, c, val2);
		}
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return;
}

#ifdef MDEBUG
	double  _copia_m_colindx_d(const MATRICEd *m, const VETTOREi *indx, int colonna, const char *nomefile, int linea)
#else
	double  _copia_m_colindx_d(const MATRICEd *m, const VETTOREi *indx, int colonna)
#endif
{
	double ris;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: copia_m_colindx_d\n", linea);
#endif
	CONTROLLA(m != NULL && indx != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	_StampaVett_i(indx);
	fprintf(fp_fdbg, "colonna = %d\n", colonna);
#endif
	CONTROLLA(colonna > 0 && colonna <= m->nc);
	ris = _ACCEDIm_d(m, _ACCEDIv_i(indx, colonna), colonna);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\nris: %.16g\n", ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	void  _assegna1_mv_riga_d(const MATRICEd *m, int riga, const VETTOREd *v, const char *nomefile, int linea)
#else
	void  _assegna1_mv_riga_d(const MATRICEd *m, int riga, const VETTOREd *v)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_mv_riga_d\n", linea);
#endif
	CONTROLLA(m != NULL && v != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "riga = %d\n", riga);
	_StampaVett_d(v);
#endif
	CONTROLLA(riga > 0 && riga <= m->nr);
	for (i = 1; i <= LENGTHv_d(v); i++)
		_ASSEGNAm_d(m, riga, i, _ACCEDIv_d(v, i));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return;
}

#ifdef MDEBUG
	MATRICEd * _aggiungi_mv_riga_d(MATRICEd *m, int riga, const VETTOREd *v, const char *nome, const char *nomefile, int linea)
#else
	MATRICEd * _aggiungi_mv_riga_d(MATRICEd *m, int riga, const VETTOREd *v)
#endif
{
	int i, r, c;
	MATRICEd *tmp_m = NULL;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: aggiungi_mv_riga_d\n", linea);
#endif
	CONTROLLA(m != NULL && v != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "riga = %d\n", riga);
	_StampaVett_d(v);
#endif
	CONTROLLA(m->nc == v->dim);
	if (riga > m->nr + 1) {
		Rprintf("la riga da aggiungere non e` minore o uguale al numero di righe attuali + 1: R riempirebbe quelle in mezzo con NA, e` questa l'intenzione?\n");
		error("");
	}
	if (riga > LENGTHm1_d(m)) {
		_CREAm_d(tmp_m, riga, m->nc);
		tmp_m->nr = riga;
		for (r = 1; r <= LENGTHm1_d(m); r++) {
			for (c = 1; c <= LENGTHm2_d(m); c++)
				_ASSEGNAm_d(tmp_m, r, c, _ACCEDIm_d(m, r, c));
		}
		_CANCELLAm_d(m);
		m = tmp_m;
	}
	for (i = 1; i <= LENGTHm2_d(m); i++)
		_ASSEGNAm_d(m, riga, i, _ACCEDIv_d(v, i));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return m;
}

#ifdef MDEBUG
	MATRICEd * _aggiungi_ms_riga_d(MATRICEd *m, int riga, double val, const char *nome, const char *nomefile, int linea)
#else
	MATRICEd * _aggiungi_ms_riga_d(MATRICEd *m, int riga, double val)
#endif
{
	int i, r, c;
	MATRICEd *tmp_m = NULL;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: aggiungi_ms_riga_d\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "riga = %d, val = %.16g\n", riga, val);
#endif
	if (riga > m->nr + 1) {
		Rprintf("la riga da aggiungere non e` minore o uguale al numero di righe attuali + 1: R riempirebbe quelle in mezzo con NA, e` questa l'intenzione?\n");
		error("");
	}
	if (riga > m->alloc_r) {
#ifdef FDEBUG
		_CREAstr(tmp, "");
		g_string_printf(tmp, "ATTENZIONE (aggiungi_ms_riga_d, linea %s # %d): ingrandite le righe della matrice da %d a %d\n", nomefile, linea, LENGTHm1_d(m), riga);
		warning(tmp->str);
		fprintf(fp_fdbg, tmp->str);
		_CANCELLAstr(tmp);
#endif
		_CREAm_d(tmp_m, 2 * riga, m->nc);
		tmp_m->nr = riga;
		for (r = 1; r <= LENGTHm1_d(m); r++) {
			for (c = 1; c <= LENGTHm2_d(m); c++)
				_ASSEGNAm_d(tmp_m, r, c, _ACCEDIm_d(m, r, c));
		}
		_CANCELLAm_d(m);
		m = tmp_m;
	}
	for (i = 1; i <= LENGTHm2_d(m); i++)
		_ASSEGNAm_d(m, riga, i, val);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return m;
}

#ifdef MDEBUG
	VETTOREd * _somma_colonne_d(VETTOREd *ris, const MATRICEd *m, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _somma_colonne_d(VETTOREd *ris, const MATRICEd *m)
#endif
{
	int r, c;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: somma_colonne_d\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
#endif
	_CREAv_d(ris, LENGTHm2_d(m));
	for (c = 1; c <= LENGTHm2_d(m); c++) {
		_ASSEGNAv_d(ris, c, 0);
		for (r = 1; r <= LENGTHm1_d(m); r++) {
			_ASSEGNAv_d(ris, c, _ACCEDIv_d(ris, c)  + _ACCEDIm_d(m, r, c));
		}
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _which_v_indxlt_d(VETTOREi *ris, const VETTOREd *v, double val, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_v_indxlt_d(VETTOREi *ris, const VETTOREd *v, double val)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_v_indxlt_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	fprintf(fp_fdbg, "val = %.16g", val);
#endif
	_CREAv_i(ris, LENGTHv_d(v));
	for (i = 1, j = 1; i <= LENGTHv_d(v); i++) {
		if (_ACCEDIv_d(v, i) < val) {
			_ASSEGNAv_i(ris, j, i);
			j++;
		}
	}
	ris->dim = j - 1;
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _unione1_d(VETTOREd *v1, const VETTOREd *v2, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _unione1_d(VETTOREd *v1, const VETTOREd *v2)
#endif
{
	int i, j;
	VETTOREd *tmp_ris = NULL, *ris = NULL;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: unione1_d\n", linea);
#endif
	CONTROLLA(v1 != NULL && v2 != NULL);
#ifdef FDEBUG
	_StampaVett_d(v1);
	_StampaVett_d(v2);
#endif
		_CREAv_d(tmp_ris, LENGTHv_d(v1) + LENGTHv_d(v2));
		for (i = 1; i <= LENGTHv_d(v1); i++)
			_ASSEGNAv_d(tmp_ris, i, _ACCEDIv_d(v1, i));
		for (j = 1; j <= LENGTHv_d(v2); j++) {
			_ASSEGNAv_d(tmp_ris, i, _ACCEDIv_d(v2, j));
			i++;
		}
		_CANCELLAv_d(v1);
		tmp_ris->dim = i - 1;
		ris = elimina_doppi_d(ris, tmp_ris);
		_CANCELLAv_d(tmp_ris);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	double  _somma_v_d(const VETTOREd *v, bool canc_NA, const char *nomefile, int linea)
#else
	double  _somma_v_d(const VETTOREd *v, bool canc_NA)
#endif
{
	int i;
	double ris = 0;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: somma_v_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	fprintf(fp_fdbg, "canc_NA = %d\n\n", canc_NA);
#endif
	for (i = 1; i <= LENGTHv_d(v); i++) {
		if (!canc_NA && ISNA(_ACCEDIv_d(v, i))) {
			ris = NA_REAL;
			break;
		}
		else if (!ISNA(_ACCEDIv_d(v, i)))
			ris += _ACCEDIv_d(v, i);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\nris: %.16g\n", ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	double  _somma_m_d(const MATRICEd *m, const char *nomefile, int linea)
#else
	double  _somma_m_d(const MATRICEd *m)
#endif
{
	int r, c;
	double ris = 0;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: somma_m_d\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
#endif
	for (r = 1; r <= LENGTHm1_d(m); r++) {
		for (c = 1; c <= LENGTHm2_d(m); c++)
			ris += _ACCEDIm_d(m, r, c);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\nris: %.16g\n", ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _setdiff_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _setdiff_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2)
#endif
{
	int indx, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: setdiff_d\n", linea);
#endif
	CONTROLLA(v1 != NULL && v2 != NULL);
#ifdef FDEBUG
	_StampaVett_d(v1);
	_StampaVett_d(v2);
#endif
	_CREAv_d(ris, LENGTHv_d(v1));
	ris = copia_v_d(ris, v1, 1, LENGTHv_d(v1));
	for (j = 1; j <= LENGTHv_d(v2); j++) {
		indx = esiste_v_d(_ACCEDIv_d(v2, j), ris);
		if (indx > 0)
			elimina_indx_d(ris, indx);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _diff_vv_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _diff_vv_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: diff_vv_d\n", linea);
#endif
	CONTROLLA(v1 != NULL && v2 != NULL);
#ifdef FDEBUG
	_StampaVett_d(v1);
	_StampaVett_d(v2);
#endif
	_CREAv_d(ris, LENGTHv_d(v1));
	for (i = 1; i <= LENGTHv_d(v1); i++)
		_ASSEGNAv_d(ris, i, _ACCEDIv_d(v1, i) - _ACCEDIv_d(v2, i));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _abs_v_d(VETTOREd *ris, const VETTOREd *v, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _abs_v_d(VETTOREd *ris, const VETTOREd *v)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: abs_v_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
#endif
	_CREAv_d(ris, LENGTHv_d(v));
	for (i = 1; i <= LENGTHv_d(v); i++)
		_ASSEGNAv_d(ris, i, fabs(_ACCEDIv_d(v, i)));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _which_v_indxna_d(VETTOREi *ris, const VETTOREd *v, bool complemento, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_v_indxna_d(VETTOREi *ris, const VETTOREd *v, bool complemento)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_v_indxna_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	fprintf(fp_fdbg, "complemento = %d\n", complemento);
#endif
	_CREAv_i(ris, LENGTHv_d(v));
	for (i = 1, j = 1; i <= LENGTHv_d(v); i++) {
		if (ISNA(_ACCEDIv_d(v, i)) && !complemento) {
			_ASSEGNAv_i(ris, j, i);
			j++;
		}
		else if (!ISNA(_ACCEDIv_d(v, i)) && complemento) {
			_ASSEGNAv_i(ris, j, i);
			j++;
		}
	}
	ris->dim = j - 1;
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _which_v_indxeq_d(VETTOREi *ris, const VETTOREd *v, double val, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_v_indxeq_d(VETTOREi *ris, const VETTOREd *v, double val)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_v_indxeq_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	fprintf(fp_fdbg, "val = %.16g\n", val);
#endif
	_CREAv_i(ris, LENGTHv_d(v));
	for (i = 1, j = 1; i <= LENGTHv_d(v); i++) {
#ifdef VERSIONE_d
		if (val == R_PosInf && _ACCEDIv_d(v, i) == R_PosInf) {
			_ASSEGNAv_i(ris, j, i);
			j++;
		}
		else if (val == R_NegInf && _ACCEDIv_d(v, i) == R_NegInf) {
			_ASSEGNAv_i(ris, j, i);
			j++;
		} else
#endif
		if (UGUALE(_ACCEDIv_d(v, i), val)) {
			_ASSEGNAv_i(ris, j, i);
			j++;
		}
	}
	ris->dim = j - 1;
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	MATRICEd * _cbind2v_d(MATRICEd *ris, const VETTOREd *v1, const VETTOREd *v2, const char *nome, const char *nomefile, int linea)
#else
	MATRICEd * _cbind2v_d(MATRICEd *ris, const VETTOREd *v1, const VETTOREd *v2)
#endif
{
	int lM, r, r1, r2;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: cbind2v_d\n", linea);
#endif
	CONTROLLA(v1 != NULL && v2 != NULL);
#ifdef FDEBUG
	_StampaVett_d(v1);
	_StampaVett_d(v2);
#endif
	lM = max_s_i(LENGTHv_d(v1), LENGTHv_d(v2));
	_CREAm_d(ris, lM, 2);
	for (r = 1, r1 = 1, r2 = 1; r <= lM; r++, r1++, r2++) {
		if (r > LENGTHv_d(v1))
			r1 = 1;
		_ASSEGNAm_d(ris, r, 1, _ACCEDIv_d(v1, r1));
		if (r > LENGTHv_d(v2))
			r2 = 1;
		_ASSEGNAm_d(ris, r, 2, _ACCEDIv_d(v2, r2));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _max_righe_d(VETTOREd *ris, const MATRICEd *m, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _max_righe_d(VETTOREd *ris, const MATRICEd *m)
#endif
{
	int r, c;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: max_righe_d\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
#endif
	_CREAv_d(ris, LENGTHm1_d(m));
	for (r = 1; r <= LENGTHm1_d(m); r++) {
		_ASSEGNAv_d(ris, r, _ACCEDIm_d(m, r, 1));
		for (c = 2; c <= LENGTHm2_d(m); c++) {
			if (_ACCEDIv_d(ris, r) < _ACCEDIm_d(m, r, c))
				_ASSEGNAv_d(ris, r, _ACCEDIm_d(m, r, c));
		}
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _f_aux_d(VETTOREd *ris, const VETTOREd *a, const VETTOREd *b, const VETTOREd *m, const VETTOREd *t, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _f_aux_d(VETTOREd *ris, const VETTOREd *a, const VETTOREd *b, const VETTOREd *m, const VETTOREd *t)
#endif
{
	int i;
	double tmp1;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: f_aux_d\n", linea);
#endif
	CONTROLLA(a != NULL && b != NULL && m != NULL && t != NULL);
#ifdef FDEBUG
	_StampaVett_d(a);
	_StampaVett_d(b);
	_StampaVett_d(m);
	_StampaVett_d(t);
#endif
	_CREAv_d(ris, LENGTHv_d(a));
	for (i = 1; i <= LENGTHv_d(a); i++) {
		if (UGUALE(_ACCEDIv_d(t, i), 0)) {
#ifdef FDEBUG
			_CREAstr(tmp, "");
			g_string_printf(tmp, "ATTENZIONE (f_aux_d, linea %s # %d): l'elemento %d ha provocato una divisione per zero e gli e` stato assegnato un valore al di fuori del dominio!\n", nomefile, linea, i);
			warning(tmp->str);
			fprintf(fp_fdbg, tmp->str);
			_CANCELLAstr(tmp);
#endif
			if (sign((_ACCEDIv_d(a, i) - _ACCEDIv_d(b, i)) * _ACCEDIv_d(m, i)) > 0)
				tmp1 = R_PosInf;
			else if (sign((_ACCEDIv_d(a, i) - _ACCEDIv_d(b, i)) * _ACCEDIv_d(m, i)) < 0)
				tmp1 = R_NegInf;
			else
				tmp1 = R_NaN;
		}
		else
			tmp1 = (double) sign(_ACCEDIv_d(a, i) - _ACCEDIv_d(b, i)) * _ACCEDIv_d(m, i) / (double) _ACCEDIv_d(t, i);
		_ASSEGNAv_d(ris, i, tmp1);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _somma_vv_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _somma_vv_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: somma_vv_d\n", linea);
#endif
	CONTROLLA(v1 != NULL && v2 != NULL);
#ifdef FDEBUG
	_StampaVett_d(v1);
	_StampaVett_d(v2);
#endif
	_CREAv_d(ris, LENGTHv_d(v1));
	for (i = 1; i <= LENGTHv_d(v1); i++)
		_ASSEGNAv_d(ris, i, _ACCEDIv_d(v1, i) + _ACCEDIv_d(v2, i));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	MATRICEd * _trasponi_d(MATRICEd *ris, const MATRICEd *m, const char *nome, const char *nomefile, int linea)
#else
	MATRICEd * _trasponi_d(MATRICEd *ris, const MATRICEd *m)
#endif
{
	int r, c;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: trasponi_d\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
#endif
	_CREAm_d(ris, LENGTHm2_d(m), LENGTHm1_d(m));

	for (r = 1; r <= LENGTHm1_d(m); r++) {
		for (c = 1; c <= LENGTHm2_d(m); c++)
			_ASSEGNAm_d(ris, c, r, _ACCEDIm_d(m, r, c));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	void  _somma1_m_d(MATRICEd *m1, const MATRICEd *m2, const char *nomefile, int linea)
#else
	void  _somma1_m_d(MATRICEd *m1, const MATRICEd *m2)
#endif
{
	int r, c;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: somma1_m_d\n", linea);
#endif
	CONTROLLA(m1 != NULL && m2 != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m1);
	_StampaMatr_d(m2);
#endif
	CONTROLLA(m1->nr == m2->nr && m1->nc == m2->nc);
	for (r = 1; r <= LENGTHm1_d(m1); r++) {
		for (c = 1; c <= LENGTHm2_d(m1); c++)
			_ASSEGNAm_d(m1, r, c, _ACCEDIm_d(m1, r, c) + _ACCEDIm_d(m2, r, c));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(m1);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	void  _assegna1_ms_indx_d(const MATRICEd *m, const VETTOREi *indx, double val, const char *nomefile, int linea)
#else
	void  _assegna1_ms_indx_d(const MATRICEd *m, const VETTOREi *indx, double val)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_ms_indx_d\n", linea);
#endif
	CONTROLLA(m != NULL && indx != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	_StampaVett_i(indx);
	fprintf(fp_fdbg, "val = %.16g\n", val);
#endif
	for (i = 1; i <= LENGTHv_i(indx); i++)
		_ASSEGNAmv_d(m, _ACCEDIv_i(indx, i), val);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	VETTOREi * _which_m_indxne_d(VETTOREi *ris, const MATRICEd *m, double val, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_m_indxne_d(VETTOREi *ris, const MATRICEd *m, double val)
#endif
{
	int r, c, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_m_indxne_d\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "val = %.16g\n", val);
#endif
	_CREAv_i(ris, LENGTHm1_d(m) * LENGTHm2_d(m));
	j = 1;
	for (c = 1; c <= LENGTHm2_d(m); c++) {
		for (r = 1; r <= LENGTHm1_d(m); r++) {
			if (DIVERSO(_ACCEDIm_d(m, r, c), val)) {
				_ASSEGNAv_i(ris, j, r + m->nr * (c - 1));
				j++;
			}
		}
	}
	ris->dim = j - 1;
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	void  _assegna1_s_diag_d(const MATRICEd *m, double val, const char *nomefile, int linea)
#else
	void  _assegna1_s_diag_d(const MATRICEd *m, double val)
#endif
{
	int l, d;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_s_diag_d\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "val = %.16g\n", val);
#endif
	l = min_s_i(LENGTHm1_d(m), LENGTHm2_d(m));
	for (d = 1; d <= l; d++)
		_ASSEGNAm_d(m, d, d, val);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	double  _media_v_d(VETTOREd *v, const char *nomefile, int linea)
#else
	double  _media_v_d(VETTOREd *v)
#endif
{
	int i, j = 0;
	double ris = 0.0;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: media_v_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
#endif
	for (i = 1; i <= LENGTHv_d(v); i++) {
		if (!ISNA(_ACCEDIv_d(v, i))) {
			ris += _ACCEDIv_d(v, i);
			j++;
		}
	}
	if (j > 0)
		ris /= j;
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\nmedia = %3.3f\n", ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	double  _somma_riga_d(const MATRICEd *m, int riga, const char *nomefile, int linea)
#else
	double  _somma_riga_d(const MATRICEd *m, int riga)
#endif
{
	int c;
	double ris = 0;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: somma_riga_d\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "riga = %d\n", riga);
#endif
	CONTROLLA(riga > 0 && riga <= m->nr);
	for (c = 1; c <= LENGTHm2_d(m); c++)
		ris += _ACCEDIm_d(m, riga, c);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\nris: %.16g\n", ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _ordine_d(VETTOREi *ris, const VETTOREd *v, bool decr, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _ordine_d(VETTOREi *ris, const VETTOREd *v, bool decr)
#endif
{
	int i;
	VETTOREd *tmp = NULL;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: ordine_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	fprintf(fp_fdbg, "decr = %d\n", decr);
#endif
	_CREAv_i(ris, LENGTHv_d(v));
	for (i = 1; i <= LENGTHv_d(v); i++)
		_ASSEGNAv_i(ris, i, i);
#ifdef VERSIONE_d
	tmp = copia_v_d(tmp, v, 1, LENGTHv_d(v));
#else
	_CREAv_d(tmp, LENGTHv_i(v));
	for (i = 1; i <= LENGTHv_i(v); i++)
		_ASSEGNAv_d(tmp, i, (double) _ACCEDIv_i(v, i));
#endif
	if (decr)
		revsort(tmp->dati, ris->dati, LENGTHv_d(tmp));
	else
		rsort_with_index(tmp->dati, ris->dati, LENGTHv_d(tmp));
	_CANCELLAv_d(tmp);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
		return ris;
}

#ifdef MDEBUG
	VETTOREd * _diag_d(VETTOREd *ris, const MATRICEd *m, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _diag_d(VETTOREd *ris, const MATRICEd *m)
#endif
{
	int l, d;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: diag_d\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
#endif
	l = min_s_i(LENGTHm1_d(m), LENGTHm2_d(m));
	_CREAv_d(ris, l);
	for (d = 1; d <= l; d++)
		_ASSEGNAv_d(ris, d, _ACCEDIm_d(m, d, d));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _which_m_colindxne_d(VETTOREi *ris, const MATRICEd *m, int colonna, double val, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_m_colindxne_d(VETTOREi *ris, const MATRICEd *m, int colonna, double val)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_m_colindxne_d\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "colonna = %d, val = %.16g\n", colonna, val);
#endif
	CONTROLLA(colonna > 0 && colonna <= LENGTHm2_d(m));
	_CREAv_i(ris, LENGTHm1_d(m) * LENGTHm2_d(m));
	for (i = 1, j = 1; i <= LENGTHm1_d(m); i++) {
		if (DIVERSO(_ACCEDIm_d(m, i, colonna), val)) {
			_ASSEGNAv_i(ris, j, i);
			j++;
		}
	}
	ris->dim = j - 1;
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _which_m_colindxin_d(VETTOREi *ris, const MATRICEd *m, int c, const VETTOREd *v, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_m_colindxin_d(VETTOREi *ris, const MATRICEd *m, int c, const VETTOREd *v)
#endif
{
	int i, j, indx;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_m_colindxin_d\n", linea);
#endif
	CONTROLLA(m != NULL && v != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "c = %d\n", c);
	_StampaVett_d(v);
#endif
	CONTROLLA(c > 0 && c <= LENGTHm2_d(m));
	_CREAv_i(ris, LENGTHm1_d(m));
	for (i = 1, j = 1; i <= LENGTHm1_d(m); i++) {
		indx = esiste_v_d(_ACCEDIm_d(m, i, c), v);
		if (indx > 0) {
			_ASSEGNAv_i(ris, j, i);
			j++;
		}
	}
	ris->dim = j - 1;
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	MATRICEd * _copia_m_ncol_d(MATRICEd *ris, const MATRICEd *m, int colonne, const char *nome, const char *nomefile, int linea)
#else
	MATRICEd * _copia_m_ncol_d(MATRICEd *ris, const MATRICEd *m, int colonne)
#endif
{
	int r, c, l, k;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: copia_m_ncol_d\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
#endif
	l = LENGTHm1_d(m) * LENGTHm2_d(m);
#ifdef FDEBUG
	if ((l > colonne && l % colonne) || (l < colonne && colonne % l)) {
		_CREAstr(tmp, "");
		g_string_printf(tmp, "ATTENZIONE (copia_m_ncol_d, linea %s # %d): il numero di elementi di m (%d) non e` un multiplo o sotto-multiplo del numero di colonne (%d)!\n", nomefile, linea, l, colonne);
		warning(tmp->str);
		fprintf(fp_fdbg, tmp->str);
		_CANCELLAstr(tmp);
	}
#endif
	if (l >= colonne)
		l /= colonne;
	else
		l = 1;
	_CREAm_d(ris, l, colonne);
	k = 0;
	for (c = 1; c <= colonne; c++) {
		for (r = 1; r <= l; r++)
			_ASSEGNAm_d(ris, r, c, _ACCEDImv_d(m, k + 1));
			k = (k % l) + 1;
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _which_m_colneand2_d(VETTOREi *ris, const MATRICEd *m, int c, double val1, double val2, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_m_colneand2_d(VETTOREi *ris, const MATRICEd *m, int c, double val1, double val2)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_m_colneand2_d\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "c = %d, val1 = %.16g, val2 = %.16g\n", c, val1, val2);
#endif
	CONTROLLA(c > 0 && c <= LENGTHm2_d(m));
	_CREAv_i(ris, LENGTHm1_d(m) * LENGTHm2_d(m));
	for (i = 1, j = 1; i <= LENGTHm1_d(m); i++) {
		if (DIVERSO(_ACCEDIm_d(m, i, c), val1) && DIVERSO(_ACCEDIm_d(m, i, c), val2)) {
			_ASSEGNAv_i(ris, j, i);
			j++;
		}
	}
	ris->dim = j - 1;
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	MATRICEd * _righe_d(MATRICEd *ris, const MATRICEd *m, const VETTOREi *indx, const char *nome, const char *nomefile, int linea)
#else
	MATRICEd * _righe_d(MATRICEd *ris, const MATRICEd *m, const VETTOREi *indx)
#endif
{
	int r, c;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: righe_d\n", linea);
#endif
	CONTROLLA(m != NULL && indx != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	_StampaVett_i(indx);
	fprintf(fp_fdbg, "\n");
#endif
	_CREAm_d(ris, LENGTHv_i(indx), LENGTHm2_d(m));
	for (r = 1; r <= LENGTHv_i(indx); r++) {
		for (c = 1; c <= LENGTHm2_d(m); c++)
			_ASSEGNAm_d(ris, r, c, _ACCEDIm_d(m, _ACCEDIv_i(indx, r), c));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	void  _assegna1_v_diag_d(const MATRICEd *m, const VETTOREd *v, const char *nomefile, int linea)
#else
	void  _assegna1_v_diag_d(const MATRICEd *m, const VETTOREd *v)
#endif
{
	int l, d;

	l = min_s_i(LENGTHm1_d(m), LENGTHm2_d(m));
#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_v_diag_d\n", linea);
#endif
	CONTROLLA(m != NULL && v != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	_StampaVett_d(v);
#endif
	CONTROLLA(l == LENGTHv_d(v));
	for (d = 1; d <= l; d++)
		_ASSEGNAm_d(m, d, d, _ACCEDIv_d(v, d));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	void  _assegna1_ms_indx2_d(const MATRICEd *m, const MATRICEi *indxm, double val, const char *nomefile, int linea)
#else
	void  _assegna1_ms_indx2_d(const MATRICEd *m, const MATRICEi *indxm, double val)
#endif
{
	int i, l;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_ms_indx2_d\n", linea);
#endif
	CONTROLLA(m != NULL && indxm != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	_StampaMatr_i(indxm);
	fprintf(fp_fdbg, "val = %.16g\n", val);
#endif
	CONTROLLA(LENGTHm2_i(indxm) == 2);
	l = min_s_i(LENGTHm1_d(m), LENGTHm1_i(indxm));
	for (i = 1; i <= l; i++)
		_ASSEGNAm_d(m, _ACCEDIm_i(indxm, i, 1), _ACCEDIm_i(indxm, i, 2), val);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	void  _incr1_v_d(const VETTOREd *v, double s, const char *nomefile, int linea)
#else
	void  _incr1_v_d(const VETTOREd *v, double s)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: incr1_v_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	fprintf(fp_fdbg, "s = %.16g\n", s);
#endif
	for (i = 1; i <= LENGTHv_d(v); i++)
		_ASSEGNAv_d(v, i, _ACCEDIv_d(v, i) + s);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(v);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	MATRICEd * _copia_m_d(MATRICEd *ris, const MATRICEd *da, const char *nome, const char *nomefile, int linea)
#else
	MATRICEd * _copia_m_d(MATRICEd *ris, const MATRICEd *da)
#endif
{
	int r, c;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: copia_m_d\n", linea);
#endif
	CONTROLLA(da != NULL);
#ifdef FDEBUG
	_StampaMatr_d(da);
#endif
	_CREAm_d(ris, LENGTHm1_d(da), LENGTHm2_d(da));
	for (r = 1; r <= LENGTHm1_d(da); r++) {
		for (c = 1; c <= LENGTHm2_d(da); c++)
			_ASSEGNAm_d(ris, r, c, _ACCEDIm_d(da, r, c));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	MATRICEd * _somma_mm_d(MATRICEd *ris, const MATRICEd *m1, const MATRICEd *m2, const char *nome, const char *nomefile, int linea)
#else
	MATRICEd * _somma_mm_d(MATRICEd *ris, const MATRICEd *m1, const MATRICEd *m2)
#endif
{
	int r, c;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: somma_mm_d\n", linea);
#endif
	CONTROLLA(m1 != NULL && m2 != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m1);
	_StampaMatr_d(m2);
#endif
	CONTROLLA(m1->nr == m2->nr && m1->nc == m2->nc);
	_CREAm_d(ris, LENGTHm1_d(m1), LENGTHm2_d(m1));

	for (r = 1; r <= LENGTHm1_d(m1); r++) {
		for (c = 1; c <= LENGTHm2_d(m1); c++)
			_ASSEGNAm_d(ris, r, c, _ACCEDIm_d(m1, r, c) + _ACCEDIm_d(m2, r, c));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _accoda1_vs_d(VETTOREd *v, double s, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _accoda1_vs_d(VETTOREd *v, double s)
#endif
{
	int i = 1, j;
	VETTOREd *tmp_v = NULL;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: accoda1_vs_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	fprintf(fp_fdbg, "s: %.16g\n", s);
#endif
	if (v->mia_alloc >= v->dim + 1) {
		_ASSEGNAv_d(v, v->dim + 1, s);
		v->dim++;
	}
	else {
#ifdef FDEBUG
	_CREAstr(tmp, "");
	g_string_printf(tmp, "ATTENZIONE (accoda1_vs_d, linea %s # %d): ingrandito il vettore da %d a %d\n", nomefile, linea, LENGTHv_d(v), LENGTHv_d(v) + 1);
	warning(tmp->str);
	fprintf(fp_fdbg, tmp->str);
	_CANCELLAstr(tmp);
#endif
		_CREAv_d(tmp_v, 2 * (LENGTHv_d(v) + 1));
		tmp_v->dim = LENGTHv_d(v) + 1;
		for (i = 1, j = 1; j <= LENGTHv_d(v); i++, j++)
			_ASSEGNAv_d(tmp_v, i, _ACCEDIv_d(v, j));
		_ASSEGNAv_d(tmp_v, i, s);
		_CANCELLAv_d(v);
		v = tmp_v;
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(v);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return v;
}

#ifdef MDEBUG
	void  _segmento1_v_d(const VETTOREd *v, int st, int end, const char *nomefile, int linea)
#else
	void  _segmento1_v_d(const VETTOREd *v, int st, int end)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: segmento1_v_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	fprintf(fp_fdbg, "start = %d\n", st);
	fprintf(fp_fdbg, "end = %d\n", end);
#endif
	CONTROLLA(st > 0 && end <= LENGTHv_d(v));
	for (i = st, j = 1; i <= end; i++, j++)
		_ASSEGNAv_d(v, j, _ACCEDIv_d(v, i));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(v);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	VETTOREd * _segmento_v_d(VETTOREd *ris, const VETTOREd *v, int st, int end, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _segmento_v_d(VETTOREd *ris, const VETTOREd *v, int st, int end)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: segmento_v_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	fprintf(fp_fdbg, "start = %d\n", st);
	fprintf(fp_fdbg, "end = %d\n", end);
#endif
	CONTROLLA(st > 0 && end <= LENGTHv_d(v));
	_CREAv_d(ris, end - st + 1);
	for (i = st, j = 1; i <= end; i++, j++)
		_ASSEGNAv_d(ris, j, _ACCEDIv_d(v, i));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _which_m_indxeq_d(VETTOREi *ris, const MATRICEd *m, double val, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_m_indxeq_d(VETTOREi *ris, const MATRICEd *m, double val)
#endif
{
	int r, c, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_m_indxeq_d\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "val = %.16g\n", val);
#endif
	_CREAv_i(ris, LENGTHm1_d(m) * LENGTHm2_d(m));
	j = 1;
	for (c = 1; c <= LENGTHm2_d(m); c++) {
		for (r = 1; r <= LENGTHm1_d(m); r++) {
			if (UGUALE(_ACCEDIm_d(m, r, c), val)) {
				_ASSEGNAv_i(ris, j, r + m->nr * (c - 1));
				j++;
			}
		}
	}
	ris->dim = j - 1;
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _dividi_vs_d(VETTOREd *ris, const VETTOREd *v, double div, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _dividi_vs_d(VETTOREd *ris, const VETTOREd *v, double div)
#endif
{
	int i;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: dividi_vs_d\n", linea);
#endif
#ifdef FDEBUG
	_StampaVett_d(v);
	fprintf(fp_fdbg, "div = %3.3f\n", div);
#endif
	_CREAv_d(ris, LENGTHv_d(v));
	if (ISNA(div)) {
#ifdef FDEBUG
	_CREAstr(tmp, "");
	g_string_printf(tmp, "ATTENZIONE (dividi_vs_d, linea %s # %d): divisione per NA!", nomefile, linea);
	warning(tmp->str);
	fprintf(fp_fdbg, tmp->str);
	_CANCELLAstr(tmp);
#endif
		for (i = 1; i <= LENGTHv_d(v); i++)
			_ASSEGNAv_d(ris, i, NA_REAL);
	}
else if (Uguale(div, 0.0)) {
#ifdef FDEBUG
	_CREAstr(tmp, "");
	g_string_printf(tmp, "ATTENZIONE (dividi_vs_d, linea %s # %d): divisione per zero!\n", nomefile, linea);
	warning(tmp->str);
	fprintf(fp_fdbg, tmp->str);
	_CANCELLAstr(tmp);
#endif
		for (i = 1; i <= LENGTHv_d(v); i++) {
			if (_ACCEDIv_d(v, i) > 0)
				_ASSEGNAv_d(ris, i, R_PosInf);
			else
				_ASSEGNAv_d(ris, i, R_NegInf);
		}
	}
	else {
		for (i = 1; i <= LENGTHv_d(v); i++)
			_ASSEGNAv_d(ris, i, _ACCEDIv_d(v, i) / div);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	void  _assegna1_v_indxv_d(const VETTOREd *v, const VETTOREi *indx, double val, const char *nomefile, int linea)
#else
	void  _assegna1_v_indxv_d(const VETTOREd *v, const VETTOREi *indx, double val)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_v_indxv_d\n", linea);
#endif
	CONTROLLA(v != NULL && indx != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	_StampaVett_i(indx);
	fprintf(fp_fdbg, "val = %.16g\n", val);
#endif
	for (i = 1; i <= LENGTHv_i(indx); i++) {
		if (_ACCEDIv_i(indx, i) <= 0 || _ACCEDIv_i(indx, i) > LENGTHv_d(v)) {
#ifdef FDEBUG
			Rprintf("il vettore degli indici contiene alla posizione %d un elemento fuori dai limiti del vettore di riferimento (linea %s # %d)!\n", i, nomefile, linea);
#else
			Rprintf("il vettore degli indici contiene alla posizione %d un elemento fuori dai limiti del vettore di riferimento!\n", i);
#endif
			error("");
		}
		_ASSEGNAv_d(v, _ACCEDIv_i(indx, i), val);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(v);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	VETTOREi * _which_v_andglt_d(VETTOREi *ris, const VETTOREd *v, double val1, double val2, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_v_andglt_d(VETTOREi *ris, const VETTOREd *v, double val1, double val2)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_v_andglt_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	fprintf(fp_fdbg, "val1 = %.16g\n", val1);
	fprintf(fp_fdbg, "val2 = %.16g\n", val2);
#endif
	_CREAv_i(ris, LENGTHv_d(v));
	for (i = 1, j = 1; i <= LENGTHv_d(v); i++) {
		if (_ACCEDIv_d(v, i) < val1 && _ACCEDIv_d(v, i) > val2) {
			_ASSEGNAv_i(ris, j, i);
			j++;
		}
	}
	ris->dim = j - 1;
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _distanza_2dvv_d(VETTOREd *ris, const VETTOREd *x1, const VETTOREd *y1, const VETTOREd *x2, const VETTOREd *y2, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _distanza_2dvv_d(VETTOREd *ris, const VETTOREd *x1, const VETTOREd *y1, const VETTOREd *x2, const VETTOREd *y2)
#endif
{
	int i, l;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: distanza_2dvv_d\n", linea);
#endif
	CONTROLLA(x1 != NULL && y1 != NULL && x2 != NULL && y2 != NULL);
#ifdef FDEBUG
	_StampaVett_d(x1);
	_StampaVett_d(y1);
	_StampaVett_d(x2);
	_StampaVett_d(y2);
#endif
	CONTROLLA(LENGTHv_d(x1) == LENGTHv_d(y1) && LENGTHv_d(x2) == LENGTHv_d(y2));
	l = min_s_i(LENGTHv_d(x1), LENGTHv_d(x2));
	_CREAv_d(ris, l);
	for (i = 1; i <= LENGTHv_d(x1); i++) {
		_ASSEGNAv_d(ris, i, sqrt(pow(_ACCEDIv_d(x2, i) - _ACCEDIv_d(x1, i), 2) + pow(_ACCEDIv_d(y2, i) - _ACCEDIv_d(y1, i), 2)));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _distanza_2dvs_d(VETTOREd *ris, const VETTOREd *x1, const VETTOREd *y1, double x2, double y2, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _distanza_2dvs_d(VETTOREd *ris, const VETTOREd *x1, const VETTOREd *y1, double x2, double y2)
#endif
{
	int i, l;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: distanza_2dvs_d\n", linea);
#endif
	CONTROLLA(x1 != NULL && y1 != NULL);
#ifdef FDEBUG
	_StampaVett_d(x1);
	_StampaVett_d(y1);
	fprintf(fp_fdbg, "x2 = %.16g\n", x2);
	fprintf(fp_fdbg, "y2 = %.16g\n", y2);
#endif
	CONTROLLA(LENGTHv_d(x1) == LENGTHv_d(y1));
	l = LENGTHv_d(x1);
	_CREAv_d(ris, l);
	for (i = 1; i <= LENGTHv_d(x1); i++) {
		_ASSEGNAv_d(ris, i, sqrt(pow(x2 - _ACCEDIv_d(x1, i), 2) + pow(y2 - _ACCEDIv_d(y1, i), 2)));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	void  _assegna1_ms_indxcol_d(const MATRICEd *m, const VETTOREi *indx, int col, double val, const char *nomefile, int linea)
#else
	void  _assegna1_ms_indxcol_d(const MATRICEd *m, const VETTOREi *indx, int col, double val)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_ms_indxcol_d\n", linea);
#endif
	CONTROLLA(m != NULL && indx != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	_StampaVett_i(indx);
	fprintf(fp_fdbg, "colonna = %d\n", col);
	fprintf(fp_fdbg, "val = %.16g\n", val);
#endif
	CONTROLLA(col > 0 && indx != NULL);
	for (i = 1; i <= LENGTHv_i(indx); i++) {
		_ASSEGNAm_d(m, _ACCEDIv_i(indx, i), col, val);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return;
}

#ifdef MDEBUG
	void  _assegna1_mv_indx_d(const MATRICEd *m, const VETTOREi *indx, const VETTOREd *v, const char *nomefile, int linea)
#else
	void  _assegna1_mv_indx_d(const MATRICEd *m, const VETTOREi *indx, const VETTOREd *v)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_mv_indx_d\n", linea);
#endif
	CONTROLLA(m != NULL && indx != NULL && v != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	_StampaVett_i(indx);
	_StampaVett_d(v);
#endif
	CONTROLLA(LENGTHv_i(indx) == LENGTHv_d(v));
	for (i = 1; i <= LENGTHv_i(indx); i++)
		_ASSEGNAmv_d(m, _ACCEDIv_i(indx, i), _ACCEDIv_d(v, i));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	VETTOREi * _arrotonda_v_d(VETTOREi *ris, const VETTOREd *v, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _arrotonda_v_d(VETTOREi *ris, const VETTOREd *v)
#endif
{
#ifdef VERSIONE_d
	int i;
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: arrotonda_v_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
#endif
#ifdef VERSIONE_d
	_CREAv_i(ris, LENGTHv_d(v));
	for (i = 1; i <= LENGTHv_d(v); i++)
		_ASSEGNAv_i(ris, i, (int) _ACCEDIv_d(v, i));
#ifdef FDEBUG
		_CREAstr(tmp, "");
		g_string_printf(tmp, "ATTENZIONE (arrotonda_v_d, linea %s # %d): il vettore e` stato arrotondato da double a intero: possibile perdita di precisione nei calcoli!", nomefile, linea);
		warning(tmp->str);
		fprintf(fp_fdbg, tmp->str);
		_CANCELLAstr(tmp);
#endif
#else
	ris = copia_v_i(ris, v, 1, LENGTHv_i(v));
#endif
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _exp_d(VETTOREd *ris, const VETTOREd *v, double val, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _exp_d(VETTOREd *ris, const VETTOREd *v, double val)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: exp_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	fprintf(fp_fdbg, "val = %.16g\n", val);
#endif
	_CREAv_d(ris, LENGTHv_d(v));
	for (i = 1; i <= LENGTHv_d(v); i++)
		_ASSEGNAv_d(ris, i, pow(_ACCEDIv_d(v, i), val));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _moltiplica_vs_d(VETTOREd *ris, const VETTOREd *v, double val, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _moltiplica_vs_d(VETTOREd *ris, const VETTOREd *v, double val)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: moltiplica_vs_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	fprintf(fp_fdbg, "val = %.16g\n", val);
#endif
	_CREAv_d(ris, LENGTHv_d(v));
	for (i = 1; i <= LENGTHv_d(v); i++)
		_ASSEGNAv_d(ris, i, _ACCEDIv_d(v, i) * val);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	void  _assegna1_v_segms_d(const VETTOREd *v, int st, int end, double val, const char *nomefile, int linea)
#else
	void  _assegna1_v_segms_d(const VETTOREd *v, int st, int end, double val)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna_v_segms_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	fprintf(fp_fdbg, "start = %d\n", st);
	fprintf(fp_fdbg, "end = %d\n", end);
	fprintf(fp_fdbg, "val = %.16g\n", val);
#endif
	CONTROLLA(st > 0 && end <= LENGTHv_d(v));
	for (i = st; i <= end; i++)
		_ASSEGNAv_d(v, i, val);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(v);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	void  _assegna1_v_segmv_d(const VETTOREd *v1, int st, int end, const VETTOREd *v2, const char *nomefile, int linea)
#else
	void  _assegna1_v_segmv_d(const VETTOREd *v1, int st, int end, const VETTOREd *v2)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna_v_segmv_d\n", linea);
#endif
	CONTROLLA(v1 != NULL && v2 != NULL);
#ifdef FDEBUG
	_StampaVett_d(v1);
	fprintf(fp_fdbg, "start = %d\n", st);
	fprintf(fp_fdbg, "end = %d\n", end);
	_StampaVett_d(v2);
#endif
	CONTROLLA(st > 0 && end <= LENGTHv_d(v1) && LENGTHv_d(v2) >= end - st + 1);
	for (i = st, j = 1; i <= end; i++, j++)
		_ASSEGNAv_d(v1, i, _ACCEDIv_d(v2, j));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(v1);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	void  _assegna1_v_indxeq_d(const VETTOREd *v, double val1, double val2, const char *nomefile, int linea)
#else
	void  _assegna1_v_indxeq_d(const VETTOREd *v, double val1, double val2)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_v_indxeq_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	fprintf(fp_fdbg, "val1 = %.16g\n", val1);
	fprintf(fp_fdbg, "val2 = %.16g\n", val2);
#endif
	for (i = 1; i <= LENGTHv_d(v); i++) {
		if (UGUALE(_ACCEDIv_d(v, i), val1))
			_ASSEGNAv_d(v, i, val2);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(v);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return;
}

#ifdef MDEBUG
	void  _assegna1_v_indxNA_d(const VETTOREd *v, double val, bool complemento, const char *nomefile, int linea)
#else
	void  _assegna1_v_indxNA_d(const VETTOREd *v, double val, bool complemento)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_v_indxNA_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	fprintf(fp_fdbg, "val = %.16g\n", val);
	fprintf(fp_fdbg, "complemento = %d\n", complemento);
#endif
	for (i = 1; i <= LENGTHv_d(v); i++) {
		if (ISNA(_ACCEDIv_d(v, i)) && !complemento)
			_ASSEGNAv_d(v, i, val);
		else if (!ISNA(_ACCEDIv_d(v, i)) && complemento)
			_ASSEGNAv_d(v, i, val);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(v);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return;
}

#ifdef MDEBUG
	void  _incr1_v_indx_d(const VETTOREd *v, const VETTOREi *indx, double val, const char *nomefile, int linea)
#else
	void  _incr1_v_indx_d(const VETTOREd *v, const VETTOREi *indx, double val)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: incr1_v_indx_d\n", linea);
#endif
	CONTROLLA(v != NULL && indx != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	_StampaVett_i(indx);
	fprintf(fp_fdbg, "val = %.16g\n", val);
#endif
	for (i = 1; i <= LENGTHv_i(indx); i++)
		_ASSEGNAv_d(v, _ACCEDIv_i(indx, i), _ACCEDIv_d(v, _ACCEDIv_i(indx, i)) + val);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(v);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	VETTOREd * _promuovi_d(VETTOREd *ris, const VETTOREd *v, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _promuovi_d(VETTOREd *ris, const VETTOREd *v)
#endif
{
#ifdef VERSIONE_i
	int i;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: promuovi_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
#endif
#ifdef VERSIONE_i
	_CREAv_d(ris, LENGTHv_i(v));
	for (i = 1; i <= LENGTHv_i(v); i++)
		_ASSEGNAv_d(ris, i, (double) _ACCEDIv_i(v, i));
#ifdef FDEBUG
		_CREAstr(tmp, "");
		g_string_printf(tmp, "ATTENZIONE (promuovi_d, linea %s # %d): il vettore e` stato promosso da intero a double: i calcoli saranno rallentati!", nomefile, linea);
		warning(tmp->str);
		fprintf(fp_fdbg, tmp->str);
		_CANCELLAstr(tmp);
#endif
#else
	ris = copia_v_d(ris, v, 1, LENGTHv_d(v));
#endif
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _f_aux1_d(VETTOREd *ris, const VETTOREd *a, const VETTOREd *b, const VETTOREd *c, int da, int a1, int sgn1, int sgn2, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _f_aux1_d(VETTOREd *ris, const VETTOREd *a, const VETTOREd *b, const VETTOREd *c, int da, int a1, int sgn1, int sgn2)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: f_aux1_d\n", linea);
#endif
	CONTROLLA(a != NULL && b != NULL && c != NULL);
#ifdef FDEBUG
	_StampaVett_d(a);
	_StampaVett_d(b);
	_StampaVett_d(c);
#endif
	CONTROLLA(da > 0 && da < LENGTHv_d(a));
	CONTROLLA(a1 >= da && a1 < LENGTHv_d(a));
	CONTROLLA(LENGTHv_d(a) == LENGTHv_d(b) && LENGTHv_d(b) == LENGTHv_d(c));
#ifdef FDEBUG
	fprintf(fp_fdbg, "da = %d\n", da);
	fprintf(fp_fdbg, "a = %d\n", a1);
	fprintf(fp_fdbg, "sgn1 = %d\n", sgn1);
	fprintf(fp_fdbg, "sgn2 = %d\n", sgn2);
#endif
	_CREAv_d(ris, a1 - da + 1);
	for (i = da, j = 1; i <= a1; i++) {
		_ASSEGNAv_d(ris, j, _ACCEDIv_d(a, i) + sgn1 * _ACCEDIv_d(b, i) + sgn2 * _ACCEDIv_d(c, i));
		j++;
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	int  _which_min_v_d(VETTOREd *v, const char *nomefile, int linea)
#else
	int  _which_min_v_d(VETTOREd *v)
#endif
{
	int i, ultimo = 1;
	double min;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_min_v_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
#endif
	min = _ACCEDIv_d(v, 1);
	for (i = 2; i <= LENGTHv_d(v); i++) {
		if (_ACCEDIv_d(v, i) < min) {
			min = _ACCEDIv_d(v, i);
			ultimo = i;
		}
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\nindx min = %d\n", ultimo);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ultimo;
}

#ifdef MDEBUG
	void  _assegna1_ms_colindx_d(const MATRICEd *m, const VETTOREi *indx, int colonna, double val, const char *nomefile, int linea)
#else
	void  _assegna1_ms_colindx_d(const MATRICEd *m, const VETTOREi *indx, int colonna, double val)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_ms_colindx_d\n", linea);
#endif
	CONTROLLA(m != NULL && indx != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	_StampaVett_i(indx);
	fprintf(fp_fdbg, "colonna = %d\n", colonna);
	fprintf(fp_fdbg, "val = %.16g\n", val);
#endif
	CONTROLLA(colonna > 0 && indx != NULL);
	for (i = 1; i <= LENGTHv_i(indx); i++) {
		_ASSEGNAm_d(m, _ACCEDIv_i(indx, i), colonna, val);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return;
}

#ifdef MDEBUG
	VETTOREi * _f_aux2_d(VETTOREi *ris, const VETTOREd *a, const VETTOREd *b, const VETTOREd *c, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _f_aux2_d(VETTOREi *ris, const VETTOREd *a, const VETTOREd *b, const VETTOREd *c)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: f_aux2_d\n", linea);
#endif
	CONTROLLA(a != NULL && b != NULL && c != NULL);
#ifdef FDEBUG
	_StampaVett_d(a);
	_StampaVett_d(b);
	_StampaVett_d(c);
#endif
	CONTROLLA(LENGTHv_d(a) == LENGTHv_d(b) && LENGTHv_d(b) == LENGTHv_d(c));
	_CREAv_i(ris, LENGTHv_d(a));
	for (i = 1, j = 1; i <= LENGTHv_d(a); i++) {
		if ((DIVERSO(sign(_ACCEDIv_d(a, i) - _ACCEDIv_d(c, i)) - sign(_ACCEDIv_d(b, i) - _ACCEDIv_d(c, i)), 0)) && (DIVERSO(sign(_ACCEDIv_d(b, i) - _ACCEDIv_d(c, i)) - sign(_ACCEDIv_d(a, i) - _ACCEDIv_d(c, i)), 0))) {
			_ASSEGNAv_i(ris, j, i);
			j++;
		}
	}
	ris->dim = j - 1;
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _f_aux3_d(VETTOREi *ris, const VETTOREd *a, const VETTOREd *b, const VETTOREd *c, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _f_aux3_d(VETTOREi *ris, const VETTOREd *a, const VETTOREd *b, const VETTOREd *c)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: f_aux3_d\n", linea);
#endif
	CONTROLLA(a != NULL && b != NULL && c != NULL);
#ifdef FDEBUG
	_StampaVett_d(a);
	_StampaVett_d(b);
	_StampaVett_d(c);
#endif
	CONTROLLA(LENGTHv_d(a) == LENGTHv_d(b) && LENGTHv_d(b) == LENGTHv_d(c));
	_CREAv_i(ris, LENGTHv_d(a));
	for (i = 1, j = 1; i <= LENGTHv_d(a); i++) {
		if ((_ACCEDIv_d(a, i) - _ACCEDIv_d(b, i)) > 0 && _ACCEDIv_d(a, i) < _ACCEDIv_d(c, i)) {
			_ASSEGNAv_i(ris, j, i);
			j++;
		}
	}
	ris->dim = j - 1;
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _unione_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _unione_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2)
#endif
{
	int i, j;
	VETTOREd *tmp_ris = NULL;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: unione_d\n", linea);
#endif
	CONTROLLA(v1 != NULL && v2 != NULL);
#ifdef FDEBUG
	_StampaVett_d(v1);
	_StampaVett_d(v2);
#endif
		_CREAv_d(tmp_ris, LENGTHv_d(v1) + LENGTHv_d(v2));
		for (i = 1; i <= LENGTHv_d(v1); i++)
			_ASSEGNAv_d(tmp_ris, i, _ACCEDIv_d(v1, i));
		for (j = 1; j <= LENGTHv_d(v2); j++) {
			_ASSEGNAv_d(tmp_ris, i, _ACCEDIv_d(v2, j));
			i++;
		}
		tmp_ris->dim = i - 1;
		ris = elimina_doppi_d(ris, tmp_ris);
		_CANCELLAv_d(tmp_ris);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	int  _which_min_d(VETTOREd *v, const char *nomefile, int linea)
#else
	int  _which_min_d(VETTOREd *v)
#endif
{
	int i, ris = 1;
	double min;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_min_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
#endif
	min = _ACCEDIv_d(v, 1);
	for (i = 2; i <= LENGTHv_d(v); i++) {
		if (_ACCEDIv_d(v, i) < min) {
			min = _ACCEDIv_d(v, i);
			ris = i;
		}
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\nris = %d\n", ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;

}

#ifdef MDEBUG
	VETTOREi * _which_m_indxrowindxeq_d(VETTOREi *ris, const MATRICEd *m, const VETTOREi *indx, double val, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_m_indxrowindxeq_d(VETTOREi *ris, const MATRICEd *m, const VETTOREi *indx, double val)
#endif
{
	int i, j, k, r;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_m_indxrowindxeq_d\n", linea);
#endif
	CONTROLLA(m != NULL);
	CONTROLLA(indx != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	_StampaVett_i(indx);
	fprintf(fp_fdbg, "val = %.16g\n", val);
#endif
	CONTROLLA(LENGTHv_i(indx) <= LENGTHm2_d(m));
	_CREAv_i(ris, LENGTHv_i(indx) * LENGTHm2_d(m));
	for (i = 1, j = 1, k = 1; i <= LENGTHm2_d(m); i++) {
		for (r = 1; r <= LENGTHv_i(indx); r++) {
			if (UGUALE(_ACCEDIm_d(m, _ACCEDIv_i(indx, r), i), val)) {
				_ASSEGNAv_i(ris, j, k);
				j++;
			}
			k++;
		}
	}
	ris->dim = j - 1;
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _elimina_doppi_d(VETTOREd *ris, const VETTOREd *v, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _elimina_doppi_d(VETTOREd *ris, const VETTOREd *v)
#endif
{
#ifdef VERSIONE_i
#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: elimina_doppi_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
#endif
	GHashTable *table = g_hash_table_new(g_int_hash, g_int_equal);
	int i, j;

	_CREAv_d(ris, LENGTHv_d(v));
	for (i = 1, j = 1; i <= LENGTHv_d(v); i++) {
		if (!g_hash_table_lookup(table, &(v->dati[i - 1]))) {
			ASSEGNAv_d(ris, j, ACCEDIv_d(v, i));
			g_hash_table_insert(table, &(v->dati[i - 1]), &(v->dati[i - 1]));
			j++;
		}
	}
	ris->dim = j - 1;
	g_hash_table_destroy(table);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
#else
	error("Funzione 'elimina_doppi_d' non implementata!\n");
	return NULL;
#endif
}

#ifdef MDEBUG
	double  _f_aux4_d(const VETTOREd *times, double res, const char *nomefile, int linea)
#else
	double  _f_aux4_d(const VETTOREd *times, double res)
#endif
{
	int i, tmp_i;
	double max = 0.0, tmp_d;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: f_aux4_d\n", linea);
#endif
	CONTROLLA(times != NULL);
#ifdef FDEBUG
	_StampaVett_d(times);
	fprintf(fp_fdbg, "res = %3.3f\n\n", res);
#endif
	if (Uguale(res, 0.0)) {
#ifdef FDEBUG
		_CREAstr(tmp, "");
		g_string_printf(tmp, "ATTENZIONE (f_aux4_d, linea %s # %d): divisione per zero!\n", nomefile, linea);
		warning(tmp->str);
		fprintf(fp_fdbg, tmp->str);
		_CANCELLAstr(tmp);
#endif
		return R_PosInf;
	}
	for (i = 1; i <= LENGTHv_d(times); i++) {
		tmp_d = _ACCEDIv_d(times, i) / res;
		tmp_i = rround(tmp_d, 0);
		tmp_d = fabs(tmp_d - tmp_i);
		if (tmp_d > max) {
			max = tmp_d;
		}
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "ris = %3.3f\n\n", max);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return max;
}

#ifdef MDEBUG
	void  _assegna1_mv_colonna_d(const MATRICEd *m, int colonna, const VETTOREd *v, const char *nomefile, int linea)
#else
	void  _assegna1_mv_colonna_d(const MATRICEd *m, int colonna, const VETTOREd *v)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_mv_colonna_d\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "colonna = %d\n", colonna);
	_StampaVett_d(v);
#endif
	CONTROLLA(colonna > 0 && colonna <= m->nc && v != NULL && m->nr == v->dim);
	for (i = 1; i <= LENGTHm1_d(m); i++)
		_ASSEGNAm_d(m, i, colonna, _ACCEDIv_d(v, i));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return;
}

#ifdef MDEBUG
	VETTOREd * _f_aux5_d(VETTOREd *ris, const VETTOREd *alpha, const VETTOREd *targ, const VETTOREd *theta, const VETTOREd *xmin, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _f_aux5_d(VETTOREd *ris, const VETTOREd *alpha, const VETTOREd *targ, const VETTOREd *theta, const VETTOREd *xmin)
#endif
{
	int i;
	double tmp1;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: f_aux5_d\n", linea);
#endif
	CONTROLLA(alpha != NULL && targ != NULL && theta != NULL && xmin != NULL);
#ifdef FDEBUG
	_StampaVett_d(alpha);
	_StampaVett_d(targ);
	_StampaVett_d(theta);
	_StampaVett_d(xmin);
#endif
	_CREAv_d(ris, LENGTHv_d(alpha));
	for (i = 1; i <= LENGTHv_d(alpha); i++) {
		tmp1 = 1 + pow(2.71828182845904523536, -_ACCEDIv_d(alpha, i) * (_ACCEDIv_d(targ, i) - _ACCEDIv_d(theta, i)));
		if (Uguale(tmp1, 0.0)) {
#ifdef FDEBUG
			_CREAstr(tmp, "");
			g_string_printf(tmp, "ATTENZIONE (f_aux5_d, linea %s # %d): l'elemento %d ha provocato una divisione per zero e gli e` stato assegnato un valore al di fuori del dominio!\n", nomefile, linea, i);
			warning(tmp->str);
			fprintf(fp_fdbg, tmp->str);
			_CANCELLAstr(tmp);
#endif
			tmp1 = R_PosInf;
		}
		else {
			_ASSEGNAv_d(ris, i, (1 / tmp1) * (1 - _ACCEDIv_d(xmin, i)) + _ACCEDIv_d(xmin, i));
		}
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _f_aux6_d(VETTOREd *ris, double res, const VETTOREd *k, const VETTOREd *targetT, const VETTOREd *n, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _f_aux6_d(VETTOREd *ris, double res, const VETTOREd *k, const VETTOREd *targetT, const VETTOREd *n)
#endif
{
	int i;
	double tmp1;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: f_aux6_d\n", linea);
#endif
	CONTROLLA(k != NULL && targetT != NULL && n != NULL);
#ifdef FDEBUG
	fprintf(fp_fdbg, "res = %3.3f\n", res);
	_StampaVett_d(k);
	_StampaVett_d(targetT);
	_StampaVett_d(n);
#endif
	_CREAv_d(ris, LENGTHv_d(k));
	for (i = 1; i <= LENGTHv_d(k); i++) {
		tmp1 = res * _ACCEDIv_d(k, i) * (_ACCEDIv_d(targetT, i) - _ACCEDIv_d(n, i));
		_ASSEGNAv_d(ris, i, tmp1);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	void  _somma1_vv_d(VETTOREd *v1, const VETTOREd *v2, const char *nomefile, int linea)
#else
	void  _somma1_vv_d(VETTOREd *v1, const VETTOREd *v2)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: somma1_vv_d\n", linea);
#endif
	CONTROLLA(v1 != NULL && v2 != NULL);
#ifdef FDEBUG
	_StampaVett_d(v1);
	_StampaVett_d(v2);
#endif
	for (i = 1; i <= LENGTHv_d(v1); i++)
		_ASSEGNAv_d(v1, i, _ACCEDIv_d(v1, i) + _ACCEDIv_d(v2, i));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(v1);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return;
}

#ifdef MDEBUG
	MATRICEd * _aggiungi_mv_colonna_d(MATRICEd *m, int colonna, const VETTOREd *v, const char *nome, const char *nomefile, int linea)
#else
	MATRICEd * _aggiungi_mv_colonna_d(MATRICEd *m, int colonna, const VETTOREd *v)
#endif
{
	int i, r, c;
	MATRICEd *tmp_m = NULL;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: aggiungi_mv_colonna_d\n", linea);
#endif
	CONTROLLA(m != NULL && v != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "colonna = %d\n", colonna);
	_StampaVett_d(v);
#endif
	CONTROLLA(m->nr == v->dim);
	if (colonna > m->nc + 1) {
		Rprintf("la colonna da aggiungere non e` minore o uguale al numero di colonne attuali + 1: R riempirebbe quelle in mezzo con NA, e` questa l'intenzione?\n");
		error("");
	}
	if (colonna > m->alloc_c) {
#ifdef FDEBUG
		_CREAstr(tmp, "");
		g_string_printf(tmp, "ATTENZIONE (aggiungi_mv_colonna_d, linea %s # %d): ingrandite le colonne della matrice da %d a %d\n", nomefile, linea, LENGTHm2_d(m), colonna);
		warning(tmp->str);
		fprintf(fp_fdbg, tmp->str);
		_CANCELLAstr(tmp);
#endif
		_CREAm_d(tmp_m, m->nr, 2 * colonna);
		tmp_m->nc = colonna;
		for (r = 1; r <= LENGTHm1_d(m); r++) {
			for (c = 1; c <= LENGTHm2_d(m); c++)
				_ASSEGNAm_d(tmp_m, r, c, _ACCEDIm_d(m, r, c));
		}
		_CANCELLAm_d(m);
		m = tmp_m;
	}
	for (i = 1; i <= LENGTHm1_d(m); i++)
		_ASSEGNAm_d(m, i, colonna, _ACCEDIv_d(v, i));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return m;
}

#ifdef MDEBUG
	int  _f_aux7_d(const VETTOREd *vettore, double scalare, const char *nomefile, int linea)
#else
	int  _f_aux7_d(const VETTOREd *vettore, double scalare)
#endif
{
	int i, indx = 0;
	double min = R_PosInf, tmp;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: f_aux7_d\n", linea);
#endif
	CONTROLLA(vettore != NULL);
#ifdef FDEBUG
	_StampaVett_d(vettore);
	fprintf(fp_fdbg, "scalare = %.16g\n\n", scalare);
#endif
	for (i = 1; i <= LENGTHv_d(vettore); i++) {
		tmp = _ACCEDIv_d(vettore, i) - scalare;
		if (tmp < min) {
			min = tmp;
			indx = i;
		}
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "ris = %d\n\n", indx);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return indx;
}

#ifdef MDEBUG
	VETTOREd * _f_aux8_d(VETTOREd *ris, const VETTOREd *a, const VETTOREd *b, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _f_aux8_d(VETTOREd *ris, const VETTOREd *a, const VETTOREd *b)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: f_aux8_d\n", linea);
#endif
	CONTROLLA(a != NULL && b != NULL);
#ifdef FDEBUG
	_StampaVett_d(a);
	_StampaVett_d(b);
#endif
	_CREAv_d(ris, LENGTHv_d(a));
	for (i = 1; i <= LENGTHv_d(a); i++) {
		_ASSEGNAv_d(ris, i, _ACCEDIv_d(a, i) * (1 - _ACCEDIv_d(b, i)) + _ACCEDIv_d(b, i));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	MATRICEd * _seleziona_colonne_d(MATRICEd *ris, const MATRICEd *m, const VETTOREi *indici, const char *nome, const char *nomefile, int linea)
#else
	MATRICEd * _seleziona_colonne_d(MATRICEd *ris, const MATRICEd *m, const VETTOREi *indici)
#endif
{
	int r, c;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: seleziona_colonne_d\n", linea);
#endif
	CONTROLLA(m != NULL && indici != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	_StampaVett_i(indici);
#endif
	_CREAm_d(ris, LENGTHm1_d(m), LENGTHv_i(indici));
	for (r = 1; r <= LENGTHm1_d(m); r++) {
		for (c = 1; c <= LENGTHv_i(indici); c++)
		_ASSEGNAm_d(ris, r, c, _ACCEDIm_d(m, r, _ACCEDIv_i(indici, c)));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	void  _cambiadim1_d(MATRICEd *m, const int nr, const int nc, const char *nomefile, int linea)
#else
	void  _cambiadim1_d(MATRICEd *m, const int nr, const int nc)
#endif
{
	int dim_r, dim_c;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: cambiadim1_d\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "nr = %d\n", nr);
	fprintf(fp_fdbg, "nc = %d\n\n", nc);
#endif
	CONTROLLA(nr > 0 || nc > 0);
	if (nr == -1 ) {
		dim_r = m->nr * m->nc / nc;
		dim_c = nc;
	}
	else {
		dim_r = nr;
		dim_c = m->nr * m->nc / nr;
	}
	if (dim_r * dim_c != m->nr * m->nc) {
		Rprintf("le nuove dimensioni non sono sottomultipli interi delle precedenti!\n");
		error("");
	}
	m->nr = dim_r;
	m->nc = dim_c;
	m->alloc_r = dim_r; // come se la riallocassi
	m->alloc_c = dim_c;
}

#ifdef MDEBUG
	MATRICEd * _moltiplica_mm_d(MATRICEd *ris, const MATRICEd *m1, const MATRICEd *m2, const char *nome, const char *nomefile, int linea)
#else
	MATRICEd * _moltiplica_mm_d(MATRICEd *ris, const MATRICEd *m1, const MATRICEd *m2)
#endif
{
	int r, c;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: moltiplica_mm_d\n", linea);
#endif
	CONTROLLA(m1 != NULL && m2 != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m1);
	_StampaMatr_d(m2);
#endif
	CONTROLLA(m1->nr == m2->nr && m1->nc == m2->nc);
	_CREAm_d(ris, m1->nr, m1->nc);
	for (r = 1; r <= LENGTHm1_d(m1); r++) {
		for (c = 1; c <= LENGTHm2_d(m1); c++)
			_ASSEGNAm_d(ris, r, c, _ACCEDIm_d(m1, r, c) * _ACCEDIm_d(m2, r, c));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	void  _moltiplica1_mv_d(const MATRICEd *m, const VETTOREd *v, const char *nomefile, int linea)
#else
	void  _moltiplica1_mv_d(const MATRICEd *m, const VETTOREd *v)
#endif
{
	int r, c, i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: moltiplica_mv_d\n", linea);
#endif
	CONTROLLA(m != NULL && v != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	_StampaVett_d(v);
#endif
	if (((m->nr * m->nc) % v->dim) != 0) {
		error("Il numero di elementi della matrice non e` un multiplo del numero di elementi del vettore!\n");
	}
	for (r = 1, i = 1; r <= LENGTHm1_d(m); r++) {
		for (c = 1; c <= LENGTHm2_d(m); c++) {
			_ASSEGNAm_d(m, r, c, _ACCEDIm_d(m, r, c) * _ACCEDIv_d(v, i));
			i++;
			if (i > LENGTHv_d(v))
				i = 1;
		}
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	VETTOREd * _dividi_vv_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _dividi_vv_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2)
#endif
{
	int i;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: dividi_vv_d\n", linea);
#endif
	CONTROLLA(v1 != NULL && v2 != NULL);
#ifdef FDEBUG
	_StampaVett_d(v1);
	_StampaVett_d(v2);
#endif
	_CREAv_d(ris, LENGTHv_d(v1));
	for (i = 1; i <= LENGTHv_d(v1); i++) {
		if (ISNA(_ACCEDIv_d(v2, i))) {
#ifdef FDEBUG
			_CREAstr(tmp, "");
			g_string_printf(tmp, "ATTENZIONE (dividi_vv_d, linea %s # %d): l'elemento %d del vettore divisore ha causato una divisione per NA!", nomefile, linea, i);
			warning(tmp->str);
			fprintf(fp_fdbg, tmp->str);
			_CANCELLAstr(tmp);
#endif
			_ASSEGNAv_d(ris, i, NA_REAL);
		}
		else if (UGUALE(_ACCEDIv_d(v2, i), 0)) {
#ifdef FDEBUG
			_CREAstr(tmp, "");
			g_string_printf(tmp, "ATTENZIONE (dividi_vv_d, linea %s # %d): l'elemento %d del vettore divisore ha causato una divisione per zero!\n", nomefile, linea, i);
			warning(tmp->str);
			fprintf(fp_fdbg, tmp->str);
			_CANCELLAstr(tmp);
#endif
			if (_ACCEDIv_d(v1, i) > 0)
				_ASSEGNAv_d(ris, i, R_PosInf);
			else
				_ASSEGNAv_d(ris, i, R_NegInf);
		}
		else
			_ASSEGNAv_d(ris, i, _ACCEDIv_d(v1, i) / _ACCEDIv_d(v2, i));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	MATRICEd * _sign_m_d(MATRICEd *ris, const MATRICEd *m, const char *nome, const char *nomefile, int linea)
#else
	MATRICEd * _sign_m_d(MATRICEd *ris, const MATRICEd *m)
#endif
{
	int r, c;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: sign_m_d\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
#endif
	_CREAm_d(ris, LENGTHm1_d(m), LENGTHm2_d(m));

	for (r = 1; r <= LENGTHm1_d(m); r++) {
		for (c = 1; c <= LENGTHm2_d(m); c++) {
			if (_ACCEDIm_d(m, r, c) > 0)
				_ASSEGNAm_d(ris, r, c, 1);
			else if (_ACCEDIm_d(m, r, c) < 0)
				_ASSEGNAm_d(ris, r, c, -1);
			else
				_ASSEGNAm_d(ris, r, c, 0);
		}
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	void  _assegna1_ms_colonna_d(const MATRICEd *m, int colonna, double val, const char *nomefile, int linea)
#else
	void  _assegna1_ms_colonna_d(const MATRICEd *m, int colonna, double val)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_ms_colonna_d\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "colonna = %d\n", colonna);
	fprintf(fp_fdbg, "val = %.16g\n", val);
#endif
	CONTROLLA(colonna > 0 && colonna <= m->nc);
	for (i = 1; i <= LENGTHm1_d(m); i++)
		_ASSEGNAm_d(m, i, colonna, val);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return;
}

#ifdef MDEBUG
	void  _dividi1_vv_d(VETTOREd *v1, const VETTOREd *v2, const char *nomefile, int linea)
#else
	void  _dividi1_vv_d(VETTOREd *v1, const VETTOREd *v2)
#endif
{
	int i;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: dividi1_vv_d\n", linea);
#endif
	CONTROLLA(v1 != NULL && v2 != NULL);
#ifdef FDEBUG
	_StampaVett_d(v1);
	_StampaVett_d(v2);
#endif
	for (i = 1; i <= LENGTHv_d(v1); i++) {
		if (ISNA(_ACCEDIv_d(v2, i))) {
#ifdef FDEBUG
			_CREAstr(tmp, "");
			g_string_printf(tmp, "ATTENZIONE (dividi_vv_d, linea %s # %d): l'elemento %d del vettore divisore ha causato una divisione per NA!", nomefile, linea, i);
			warning(tmp->str);
			fprintf(fp_fdbg, tmp->str);
			_CANCELLAstr(tmp);
#endif
			_ASSEGNAv_d(v1, i, NA_REAL);
		}
		else if (UGUALE(_ACCEDIv_d(v2, i), 0)) {
#ifdef FDEBUG
			_CREAstr(tmp, "");
			g_string_printf(tmp, "ATTENZIONE (dividi_vv_d, linea %s # %d): l'elemento %d del vettore divisore ha causato una divisione per zero!\n", nomefile, linea, i);
			warning(tmp->str);
			fprintf(fp_fdbg, tmp->str);
			_CANCELLAstr(tmp);
#endif
			if (_ACCEDIv_d(v1, i) > 0)
				_ASSEGNAv_d(v1, i, R_PosInf);
			else
				_ASSEGNAv_d(v1, i, R_NegInf);
		}
		else
			_ASSEGNAv_d(v1, i, _ACCEDIv_d(v1, i) / _ACCEDIv_d(v2, i));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(v1);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return;
}

#ifdef MDEBUG
	void  _abs1_v_d(VETTOREd *v, const char *nomefile, int linea)
#else
	void  _abs1_v_d(VETTOREd *v)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: abs1_v_d\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
#endif
	for (i = 1; i <= LENGTHv_d(v); i++)
		_ASSEGNAv_d(v, i, fabs(_ACCEDIv_d(v, i)));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(v);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	MATRICEi * _arrotonda_m_d(MATRICEi *ris, const MATRICEd *m, const char *nome, const char *nomefile, int linea)
#else
	MATRICEi * _arrotonda_m_d(MATRICEi *ris, const MATRICEd *m)
#endif
{
#ifdef VERSIONE_d
	int r, c;
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: arrotonda_m_d\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
#endif
#ifdef VERSIONE_d
	_CREAm_i(ris, LENGTHm1_d(m), LENGTHm2_d(m));
	for (r = 1; r <= LENGTHm1_d(m); r++) {
		for (c = 1; c <= LENGTHm2_d(m); c++)
			_ASSEGNAm_i(ris, r, c, (int) _ACCEDIm_d(m, r, c));
	}
#ifdef FDEBUG
		_CREAstr(tmp, "");
		g_string_printf(tmp, "ATTENZIONE (arrotonda_m_d, linea %s # %d): la matrice e` stata arrotondata da double a intero: possibile perdita di precisione nei calcoli!", nomefile, linea);
		warning(tmp->str);
		fprintf(fp_fdbg, tmp->str);
		_CANCELLAstr(tmp);
#endif
#else
	ris = copia_m_i(ris, m);
#endif
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#undef VERSIONE_d
