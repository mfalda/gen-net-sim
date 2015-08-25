#define VERSIONE_i

#include "r_aux_i.h"

#ifndef NDEBUG
void _stampa_lista_i(const char *pref, GList *lista)
{
	GList *li;
	Allocazione *t;

	fprintf(fp_mdbg_i, "%s", pref);
	for(li = lista; li != NULL; li = li->next) {
		t = (Allocazione *) li->data;
		fprintf(fp_mdbg_i, "%s::%s\t", t->file_da->str, t->nome->str);
	}
	fprintf(fp_mdbg_i, "\n\n");
}
#endif

void _InitDbg_i(bool stdout1)
{
#ifdef MDEBUG
	fpm_v_i = fopen("memoria_v_i.csv", "a");
	rewind(fpm_v_i);
	fprintf(fpm_v_i, "\"ID\";\"nome\";\"da\";\"a\";\"dim.massima\";\"riallocazioni\";\"prima linea\";\"ultima linea\"\n");
	allocVett_i = NULL;
	fpm_m_i = fopen("memoria_m_i.csv", "a");
	rewind(fpm_m_i);
	fprintf(fpm_m_i, "\"ID\";\"nome\";\"da\";\"a\";\"dim.max.righe\";\"dim.max.col.\";\"riallocazioni\";\"prima linea\";\"ultima linea\"\n");
	allocMatr_i = NULL;
	if (stdout1)
		fp_mdbg_i = stdout;
	else {
		fp_mdbg_i = fopen("memoria_i.txt", "a");
		rewind(fp_mdbg_i);
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

void _Fflush_i()
{
#ifdef MDEBUG
	fflush(fp_mdbg_i);
#endif
#ifdef FDEBUG
	fflush(fp_fdbg);
	return;
#endif
}

void _StampaRawVett_i(const VETTOREi *v)
{
#ifdef FDEBUG
	int i;

	if (v == NULL) {
		return;
	}
	for (i = 1; i <= LENGTHv_i(v); i++) {
		if (!(v->dati[i - 1] == NA_INTEGER))
			fprintf(fp_det, " %d", v->dati[i - 1]);
		else
			fprintf(fp_det, " NA");
	}
	fprintf(fp_det, "\n");
#endif
}

void _StampaVett_i(const VETTOREi *v)
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
	for (i = 1; i <= LENGTHv_i(v); i++)
		fprintf(fp_fdbg, " %d", v->dati[i - 1]);
	fprintf(fp_fdbg, " ]\n");
#endif
}

void _StampaRawMatr_i(const MATRICEi *m)
{
#ifdef FDEBUG
	int r, c;

	if (m == NULL) {
		return;
	}
	for (c = 1; c <= LENGTHm2_i(m); c++) {
		for (r = 1; r <= LENGTHm1_i(m); r++) {
			fprintf(fp_det, " %d", m->dati[(r - 1) + m->nr * (c - 1)]);
		}
	}
	fprintf(fp_det, "\n");
#endif
}

void _StampaMatr_i(const MATRICEi *m)
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
	for (r = 1; r <= LENGTHm1_i(m); r++) {
		fprintf(fp_fdbg, "\t");
		for (c = 1; c <= LENGTHm2_i(m); c++) {
			fprintf(fp_fdbg, " %d", m->dati[(r - 1) + m->nr * (c - 1)]);
		}
		fprintf(fp_fdbg, "\n");
	}
	fprintf(fp_fdbg, " ]\n");
#endif
}

#ifdef MDEBUG
	void  _InitVett_i(const VETTOREi *v, int val, const char *nomefile, int linea)
#else
	void  _InitVett_i(const VETTOREi *v, int val)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: InitVett_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
#endif
	for (i = 1; i <= LENGTHv_i(v); i++)
		_ASSEGNAv_i(v, i, val);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(v);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	void  _InitMatr_i(const MATRICEi *m, int val, const char *nomefile, int linea)
#else
	void  _InitMatr_i(const MATRICEi *m, int val)
#endif
{
	int r, c;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: InitMatr_i\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
#endif
	for (r = 1; r <= LENGTHm1_i(m); r++) {
		for (c = 1; c <= LENGTHm2_i(m); c++)
			_ASSEGNAm_i(m, r, c, val);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	int  _max_s_i(int a, int b, const char *nomefile, int linea)
#else
	int  _max_s_i(int a, int b)
#endif
{
	return (a > b)? a : b;
}

#ifdef MDEBUG
	int  _min_s_i(int a, int b, const char *nomefile, int linea)
#else
	int  _min_s_i(int a, int b)
#endif
{
	return (a < b)? a : b;
}

#ifdef MDEBUG
	VETTOREi * _creaVett_i(VETTOREi *ris, int dim, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _creaVett_i(VETTOREi *ris, int dim)
#endif
{
#ifdef MDEBUG
	Allocazione *mem;
#endif
	VETTOREi *ris1 = NULL;

	if (ris == NULL) {
		ris1 = mia_alloc(1, VETTOREi);
		if (ris1 == NULL) {
			Rprintf("Not enough memory (creaVett_i # %d, ris1)", __LINE__ - 2);
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
		error("Un vettore proveniente da R non puo` essere alterato! (CREAv_i)");
	if (ris1->mem == NULL) {
		mem = mia_alloc(1, Allocazione);
		if (mem == NULL) {
			Rprintf("Not enough memory (creaVett_i # %d, ris1->mem)", __LINE__ - 2);
			error("");
		}
		fprintf(fp_mdbg_i, "%d_i - allocazione del vettore '%s'[%d] (%p) dalla linea %s # %d\n", g_list_length(allocVett_i) + 1, nome, ris1->dim, ris1, nomefile, linea);
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
		allocVett_i = g_list_append(allocVett_i, mem);
		mem->indx = g_list_length(allocVett_i);
		_stampa_lista_i("allocVett_i+: ", allocVett_i);
		ris1->mem = g_list_last(allocVett_i);
	}
	else {
		mem = (Allocazione *) ris1->mem->data;
		if (mem->linea_da > 0 && mem->linea_a == 0) {
			fprintf(fp_mdbg_i, "%d_i - il vettore '%s' (linea %s # %d) esiste gia`: e` stato allocato alla linea %s # %d\n\n", mem->indx, nome, nomefile, linea, mem->file_da->str, mem->linea_da);
		}
	}
	if (ris1->mia_alloc < dim) {
		fprintf(fp_mdbg_i, "*** il vettore '%s' verra` riallocato passando da %d a %d (%d riallocazione/i, finora)\n\n", nome, ris1->mia_alloc, ris1->mia_alloc * 2 + dim, mem->rialloc + 1);
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
			ris1->dati = mia_alloc((ris1->mia_alloc * 2 + dim), int);
		}
		else if (ris1->dati == NULL)
			ris1->dati = mia_alloc(dim, int);
		if (ris1->dati == NULL) { // non puo` mai essere zero la dimensione richiesta, dato che dim > 0
			Rprintf("Not enough memory (creaVett_i # %d, ris1->dati)", __LINE__ - 2);
			error("");
		}
	}
	ris1->dim = dim;
	return ris1;
}

#ifdef MDEBUG
	MATRICEi * _creaMatr_i(MATRICEi *ris, int nr, int nc, const char *nome, const char *nomefile, int linea)
#else
	MATRICEi * _creaMatr_i(MATRICEi *ris, int nr, int nc)
#endif
{
#ifdef MDEBUG
	Allocazione *mem;
#endif
	MATRICEi *ris1;

	if (ris == NULL) {
		ris1 = mia_alloc(1, MATRICEi);
		if (ris1 == NULL) {
			Rprintf("Not enough memory (creaMatr_i # %d, ris1)", __LINE__ - 2);
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
		error("Una matrice proveniente da R non puo` essere alterata! (CREAm_i)");
	if (ris1->mem == NULL) {
		mem = mia_alloc(1, Allocazione);
		if (mem == NULL) {
			Rprintf("Not enough memory (creaMatr_i # %d, ris1->mem)", __LINE__ - 2);
			error("");
		}
		fprintf(fp_mdbg_i, "%d_i - allocazione della matrice '%s'[%d x %d] (%p) dalla linea %s # %d\n", g_list_length(allocMatr_i) + 1, nome, ris1->nr, ris1->nc, ris1, nomefile, linea);
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
		allocMatr_i = g_list_append(allocMatr_i, mem);
		mem->indx = g_list_length(allocMatr_i);
		_stampa_lista_i("allocMatr_i+: ", allocMatr_i);
		ris1->mem = g_list_last(allocMatr_i);
	}
	else {
		mem = (Allocazione *) ris1->mem->data;
		if (mem->linea_da > 0 && mem->linea_a == 0) {
			fprintf(fp_mdbg_i, "%d_i - la matrice '%s' (linea %s # %d) esiste gia`: e` stata allocata alla linea %s # %d\n\n", mem->indx, nome, nomefile, linea, mem->file_da->str, mem->linea_da);
		}
	}
	if (ris1->alloc_r < nr || ris1->alloc_c < nc) {
		fprintf(fp_mdbg_i, "*** la matrice '%s' verra` riallocata passando da %d x %d a %d x %d (%d riallocazione/i, finora)\n\n", nome, ris1->nr, ris1->nc, ris1->nr * 2 + nr, ris1->nc * 2 + nc, mem->rialloc + 1);
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
			ris1->dati = mia_alloc(ris1->alloc_r * ris1->alloc_c, int);
		}
		else if (ris1->alloc_r < nr) {
			libera(ris1->dati);
			ris1->dati = NULL;
			ris1->alloc_r = ris1->alloc_r * 2 + nr;
			ris1->dati = mia_alloc(ris1->alloc_r * ris1->alloc_c, int);
		}
		else if (ris1->alloc_c < nc) {
			libera(ris1->dati);
			ris1->dati = NULL;
		   ris1->alloc_c = ris1->alloc_c * 2 + nc;
			ris1->dati = mia_alloc(ris1->alloc_r * ris1->alloc_c, int);
		}
		else if (ris1->dati == NULL)
			ris1->dati = mia_alloc(nr * nc, int);
		if (ris1->dati == NULL) { // non puo` mai essere zero la dimensione richiesta, dato che prima ho controllato la condizione ris1->alloc_r > 0 && ris1->alloc_c > 0
			Rprintf("Not enough memory (creaMatr_i # %d, ris1->dati)", __LINE__ - 2);
			error("");
		}
	}
	ris1->nr = nr;
	ris1->nc = nc;
	return ris1;
}

#ifdef MDEBUG
	int  _lengthVett_i(const VETTOREi *v, const char *nome, const char *nomefile, int linea)
#else
	int  _lengthVett_i(const VETTOREi *v, const char *nome)
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
	int  _righeMatr_i(const MATRICEi *m, const char *nome, const char *nomefile, int linea)
#else
	int  _righeMatr_i(const MATRICEi *m, const char *nome)
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
	int  _colonneMatr_i(const MATRICEi *m, const char *nome, const char *nomefile, int linea)
#else
	int  _colonneMatr_i(const MATRICEi *m, const char *nome)
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
	int  _accediVett_i(const VETTOREi *v, int indx, const char *nome, const char *nomefile, int linea)
#else
	int  _accediVett_i(const VETTOREi *v, int indx, const char *nome)
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
	int  _accediMatr_i(const MATRICEi *m, int r, int c, const char *nome, const char *nomefile, int linea)
#else
	int  _accediMatr_i(const MATRICEi *m, int r, int c, const char *nome)
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
	return m->dati[r + LENGTHm1_i(m) * c];
}

#ifdef MDEBUG
	int  _accediMVett_i(const MATRICEi *m, int indx, const char *nome, const char *nomefile, int linea)
#else
	int  _accediMVett_i(const MATRICEi *m, int indx, const char *nome)
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
	void  _assegnaVett_i(const VETTOREi *v, int indx, int val, const char *nome, const char *nomefile, int linea)
#else
	void  _assegnaVett_i(const VETTOREi *v, int indx, int val, const char *nome)
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
	void  _assegnaMatr_i(const MATRICEi *m, int r, int c, int val, const char *nome, const char *nomefile, int linea)
#else
	void  _assegnaMatr_i(const MATRICEi *m, int r, int c, int val, const char *nome)
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
	m->dati[r + LENGTHm1_i(m) * c] = val;
}

#ifdef MDEBUG
	void  _assegnaMVett_i(const MATRICEi *m, int indx, int val, const char *nome, const char *nomefile, int linea)
#else
	void  _assegnaMVett_i(const MATRICEi *m, int indx, int val, const char *nome)
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
	VETTOREi * _cancellaVett_i(VETTOREi *v, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _cancellaVett_i(VETTOREi *v)
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
	fprintf(fp_mdbg_i, "%d_i - disallocazione del vettore '%s'[%d] (%p) dalla linea %s # %d\n", mem->indx, mem->nome->str, LENGTHv_i(v), v, nomefile, linea);
	mem->linea_a = linea;
	fprintf(fpm_v_i, "\"%d\";\"'%s'\";\"%s # %d\";\"%s # %d\";\"%d\";\"%d\";\"%s # %d\";\"%s # %d\"\n",  mem->indx, mem->nome->str, mem->file_da->str, mem->linea_da, mem->file_a->str, mem->linea_a, mem->max_dim1, mem->rialloc, mem->file_prima->str, mem->prima_linea, mem->file_ultima->str, mem->ultima_linea);
#endif
	if (LENGTHv_i(v) > 0 && !v->r) {
		libera(v->dati);
		v->dati = NULL;
	}
#ifdef MDEBUG
	_CANCELLAstr(mem->nome);
	_CANCELLAstr(mem->file_da);
	_CANCELLAstr(mem->file_a);
	_CANCELLAstr(mem->file_prima);
	_CANCELLAstr(mem->file_ultima);
	allocVett_i = g_list_remove(allocVett_i, mem);
	libera(mem);
	mem = NULL;
	_stampa_lista_i("allocVett_i-: ", allocVett_i);
#endif
	libera(v);
	v = NULL;
	return NULL;
}

#ifdef MDEBUG
	MATRICEi * _cancellaMatr_i(MATRICEi *m, const char *nome, const char *nomefile, int linea)
#else
	MATRICEi * _cancellaMatr_i(MATRICEi *m)
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
	fprintf(fp_mdbg_i, "%d_i - disallocazione della matrice '%s'[%d x %d] (%p) dalla linea %s # %d\n", mem->indx, mem->nome->str, LENGTHm1_i(m), LENGTHm2_i(m), m, nomefile, linea);
	mem->linea_a = linea;
	fprintf(fpm_m_i, "\"%d\";\"'%s'\";\"%s # %d\";\"%s # %d\";\"%d\";\"%d\";\"%d\";\"%s # %d\";\"%s # %d\"\n", mem->indx, mem->nome->str, mem->file_da->str, mem->linea_da, mem->file_a->str, mem->linea_a, mem->max_dim1, mem->max_dim2, mem->rialloc, mem->file_prima->str, mem->prima_linea, mem->file_ultima->str, mem->ultima_linea);
#endif
	if (LENGTHm1_i(m) > 0 && LENGTHm2_i(m) > 0 && !m->r) {
		libera(m->dati);
		m->dati = NULL;
	}
#ifdef MDEBUG
	_CANCELLAstr(mem->nome);
	_CANCELLAstr(mem->file_da);
	_CANCELLAstr(mem->file_a);
	_CANCELLAstr(mem->file_prima);
	_CANCELLAstr(mem->file_ultima);
	allocMatr_i = g_list_remove(allocMatr_i, mem);
	libera(mem);
	mem = NULL;
	_stampa_lista_i("allocMatr_i-: ", allocMatr_i);
#endif
	libera(m);
	m = NULL;
	return NULL;
}

void _infoVett_i(const VETTOREi *v)
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

void _infoMatr_i(const MATRICEi *m)
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
	void  _controllaCanc_i(const char *nomefile, int linea)
#else
	void  _controllaCanc_i()
#endif
{
#ifdef MDEBUG
	int ok = 1;
	Allocazione *t;

	fprintf(fp_mdbg_i, "--------------------------------------\n");
	fprintf(fpm_v_i, "--------------------------------------\n");
	while (allocVett_i != NULL) {
		t = (Allocazione *) allocVett_i->data;
		Rprintf("Il vettore '%s' di tipo 'int' (allocato alla linea %s # %d, indice %d) non e` stato ancora disallocato: lo faccio adesso.\n", t->nome->str, t->file_da->str, t->linea_da, t->indx);
		t->indir = (size_t) (VETTOREi *) _cancellaVett_i((VETTOREi *) t->indir, t->nome->str, nomefile, linea);
		ok = 0;
		_CANCELLAstr(t->nome);
		_CANCELLAstr(t->file_da);
		_CANCELLAstr(t->file_a);
		_CANCELLAstr(t->file_prima);
		_CANCELLAstr(t->file_ultima);
		libera(t);
		t = NULL;
	}
	g_list_free(allocVett_i);
	fprintf(fpm_v_i, "\n");
	fclose(fpm_v_i);
	fprintf(fpm_m_i, "--------------------------------------\n");
	while (allocMatr_i != NULL) {
		t = (Allocazione *) allocMatr_i->data;
		Rprintf("La matrice '%s' di tipo 'int' (allocata alla linea %s # %d, indice %d) non e` stata ancora disallocata: lo faccio adesso\n", t->nome->str, t->file_da->str, t->linea_da, t->indx);
		t->indir = (size_t) (MATRICEi *) _cancellaMatr_i((MATRICEi *) t->indir, t->nome->str, nomefile, linea);
		ok = 0;
		_CANCELLAstr(t->nome);
		_CANCELLAstr(t->file_da);
		_CANCELLAstr(t->file_a);
		_CANCELLAstr(t->file_ultima);
		_CANCELLAstr(t->file_prima);
		libera(t);
		t = NULL;
	}
	g_list_free(allocMatr_i);
	fprintf(fpm_m_i, "\n");
	fclose(fpm_m_i);
	if (!ok)
		fprintf(fp_mdbg_i, "ATTENZIONE: per il bilanciamento delle allocazioni di elementi di tipo 'int' sono state necessarie disallocazioni automatiche.\n");
	fprintf(fp_mdbg_i, "--------------------------------------\n");
	if (fp_mdbg_i != NULL) {
		fclose(fp_mdbg_i);
		fp_mdbg_i = NULL;
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
	SEXP  _daVETTORE_i(VETTOREi *v, int *nProtected, const char *nomefile, int linea)
#else
	SEXP  _daVETTORE_i(VETTOREi *v, int *nProtected)
#endif
{
	int i;
	int *ris1;
	SEXP ris;
#ifdef MDEBUG
	char *nome_vett;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "Trasformo il vettore ");
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
#endif
	CONTROLLA(v != NULL);
	PROTECT(ris = allocVector(INTSXP, LENGTHv_i(v)));
	(*nProtected)++;
	ris1 = INTEGER_POINTER(ris);
	for (i = 0; i < LENGTHv_i(v); i++)
		ris1[i] = _ACCEDIv_i(v, i + 1);
#ifdef MDEBUG
	nome_vett = ((Allocazione *) v->mem->data)->nome->str;
	v = _cancellaVett_i(v, nome_vett, nomefile, linea);
#else
	CANCELLAv_i(v);
#endif
	return ris;
}

#ifdef MDEBUG
	SEXP  _daMATRICE_i(MATRICEi *m, int *nProtected, const char *nomefile, int linea)
#else
	SEXP  _daMATRICE_i(MATRICEi *m, int *nProtected)
#endif
{
	int i;
	int *ris1;
	SEXP dim, ris;
#ifdef MDEBUG
	char *nome_matr;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "Trasformo la matrice\n");
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
#endif
	PROTECT(ris = allocMatrix(INTSXP, LENGTHm1_i(m), LENGTHm2_i(m)));
	(*nProtected)++;
	ris1 = INTEGER_POINTER(ris);
	for (i = 0; i < LENGTHm1_i(m) * LENGTHm2_i(m); i++)
		ris1[i] = _ACCEDImv_i(m, i + 1);
	PROTECT(dim  =  allocVector(INTSXP,  2));
	(*nProtected)++;
#ifdef MDEBUG
	nome_matr = ((Allocazione *) m->mem->data)->nome->str;
	m = _cancellaMatr_i(m, nome_matr, nomefile, linea);
#else
	CANCELLAm_i(m);
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _inVETTORE_i(SEXP s, int *nProtected, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _inVETTORE_i(SEXP s, int *nProtected)
#endif
{
	VETTOREi *v = NULL;
#ifdef MDEBUG
	Allocazione *mem;
#endif

	if (isNull(s)) {
#ifdef FDEBUG
	fprintf(fp_fdbg, "Il vettore '%s' e` nullo\n", nome);
#endif
		return NULL;
	}
	PROTECT(s = AS_INTEGER(s));
	(*nProtected)++;
#ifdef MDEBUG
	v = _creaVett_i(v, 0, nome, nomefile, linea);
#else
	v = _creaVett_i(v, 0);
#endif
	v->dim = length(s);
	v->mia_alloc = length(s);
	v->dati = INTEGER_POINTER(s);
	v->r = 1;
#ifdef MDEBUG
	mem = v->mem->data;
	g_string_assign(mem->nome, nome);
	g_string_assign(mem->file_da, nomefile);
	mem->linea_da = linea;
#endif
#ifdef FDEBUG
	fprintf(fp_fdbg, "Ho trasformato il vettore ");
	_StampaVett_i(v);
#endif
	return v;
}

#ifdef MDEBUG
	MATRICEi * _inMATRICE_i(SEXP s, int *nProtected, const char *nome, const char *nomefile, int linea)
#else
	MATRICEi * _inMATRICE_i(SEXP s, int *nProtected)
#endif
{
	MATRICEi *m = NULL;
#ifdef MDEBUG
	Allocazione *mem;
#endif

	if (isNull(s)) {
#ifdef FDEBUG
	fprintf(fp_fdbg, "La matrice '%s' e` nulla\n", nome);
#endif
		return NULL;
	}
	PROTECT(s = AS_INTEGER(s));
	(*nProtected)++;
#ifdef MDEBUG
	m = _creaMatr_i(m, 0, 0, nome, nomefile, linea);
#else
	m = _creaMatr_i(m, 0, 0);
#endif
	m->nr = nrows(s);
	m->nc = ncols(s);
	m->alloc_r = nrows(s);
	m->alloc_c = ncols(s);
	m->dati = INTEGER_POINTER(s);
	m->r = 1;
 #ifdef MDEBUG
	mem = m->mem->data;
	g_string_assign(mem->nome, nome);
	g_string_assign(mem->file_da, nomefile);
	mem->linea_da = linea;
#endif
#ifdef FDEBUG
	fprintf(fp_fdbg, "Ho trasformato la matrice ");
	_StampaMatr_i(m);
#endif
	return m;
}

#ifdef MDEBUG
	VETTOREi * _which_m_rowindxne_i(VETTOREi *ris, const MATRICEi *m, int riga, int val, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_m_rowindxne_i(VETTOREi *ris, const MATRICEi *m, int riga, int val)
#endif
{
	int i, j;

	CONTROLLA(m != NULL);
#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_m_rowindxne_i\n", linea);
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "riga = %d, val = %d\n", riga, val);
#endif
	CONTROLLA(riga > 0 && riga <= m->nr);
	_CREAv_i(ris, LENGTHm1_i(m) * LENGTHm2_i(m));
	for (i = 1, j = 1; i <= LENGTHm2_i(m); i++) {
		if (DIVERSO(_ACCEDIm_i(m, riga, i), val)) {
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
	VETTOREi * _copia_v_i(VETTOREi *ris, const VETTOREi *da, int st, int end, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _copia_v_i(VETTOREi *ris, const VETTOREi *da, int st, int end)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: copia_v_i\n", linea);
#endif
	CONTROLLA(da != NULL);
#ifdef FDEBUG
	_StampaVett_i(da);
	fprintf(fp_fdbg, "start = %d\n", st);
	fprintf(fp_fdbg, "end = %d\n", end);
#endif
	CONTROLLA(st > 0 && end <= LENGTHv_i(da));
	_CREAv_i(ris, LENGTHv_i(da));
	for (i = st; i <= end; i++)
		_ASSEGNAv_i(ris, i, _ACCEDIv_i(da, i));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _op_ss_seqdiv_i(VETTOREd *ris, int da, int a, double div, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _op_ss_seqdiv_i(VETTOREd *ris, int da, int a, double div)
#endif
{
	int i, n;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: op_ss_seqdiv_i\nda = %d, a = %d, div = %3.3f\n", linea, da, a, div);
#endif
	CONTROLLA(da < a);
	_CREAv_d(ris, a - da + 1);
	if (ISNA(div)) {
#ifdef FDEBUG
		_CREAstr(tmp, "");
		g_string_printf(tmp, "ATTENZIONE (op_ss_seqdiv_i, linea %s # %d): divisione per NA!", nomefile, linea);
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
		g_string_printf(tmp, "ATTENZIONE (op_ss_seqdiv_i, linea %s # %d): divisione per zero!\n", nomefile, linea);
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
			_ASSEGNAv_d(ris, i, (int) n / div);
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
	VETTOREi * _complementa_i(VETTOREi *ris, const VETTOREi *v, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _complementa_i(VETTOREi *ris, const VETTOREi *v)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: complementa_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
#endif
	_CREAv_i(ris, LENGTHv_i(v));
	for (i = 1; i <= LENGTHv_i(v); i++)
		_ASSEGNAv_i(ris, i, 1 - _ACCEDIv_i(v, i));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _rep_i(VETTOREi *ris, const VETTOREi *v, int ripetizioni, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _rep_i(VETTOREi *ris, const VETTOREi *v, int ripetizioni)
#endif
{
	int i, j;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: rep_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
	fprintf(fp_fdbg, "ripetizioni = %d\n", ripetizioni);
	if (ripetizioni == 0) {
		_CREAstr(tmp, "");
		g_string_printf(tmp, "ATTENZIONE (rep_i, %s # %d): zero ripetizioni!\n", nomefile, linea);
		warning(tmp->str);
		fprintf(fp_fdbg, tmp->str);
		_CANCELLAstr(tmp);
	}
#endif
	_CREAv_i(ris, LENGTHv_i(v) * ripetizioni);
	for (j = 1; j <= LENGTHv_i(v); j++) {
		for (i = 0; i < ripetizioni; i++)
			_ASSEGNAv_i(ris, j + i * LENGTHv_i(v), _ACCEDIv_i(v, j));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _seq_i(VETTOREi *ris, int da, int a, int incremento, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _seq_i(VETTOREi *ris, int da, int a, int incremento)
#endif
{
	int i, tot;
	int n;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: seq_i\nda = %d, a = %d, incremento = %d\n", linea, da, a, incremento);
#endif
	CONTROLLA(da <= a && incremento > 0);
	tot = ceil((a - da + 1) / incremento);
	_CREAv_i(ris, tot);

	for (i = 1, n = da; n <= a && i <= tot; n += incremento, i++)
		_ASSEGNAv_i(ris, i, n);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _vettore2s_i(VETTOREi *ris, int el1, int el2, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _vettore2s_i(VETTOREi *ris, int el1, int el2)
#endif
{
#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: vettore2s_i\nel1 = %d, el2 = %d\n", linea, el1, el2);
#endif
	_CREAv_i(ris, 2);
	_ASSEGNAv_i(ris, 1, el1);
	_ASSEGNAv_i(ris, 2, el2);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _vettore3s_i(VETTOREi *ris, int el1, int el2, int el3, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _vettore3s_i(VETTOREi *ris, int el1, int el2, int el3)
#endif
{
#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: vettore3s_i\nel1 = %d, el2 = %d, el3 = %d\n", linea, el1, el2, el3);
#endif
	_CREAv_i(ris, 3);
	_ASSEGNAv_i(ris, 1, el1);
	_ASSEGNAv_i(ris, 2, el2);
	_ASSEGNAv_i(ris, 3, el3);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _vettore2v_i(VETTOREi *ris, const VETTOREi *v1, const VETTOREi *v2, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _vettore2v_i(VETTOREi *ris, const VETTOREi *v1, const VETTOREi *v2)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: vettore2v_i\n", linea);
#endif
	CONTROLLA(v1 != NULL && v2 != NULL);
#ifdef FDEBUG
	_StampaVett_i(v1);
	_StampaVett_i(v2);
#endif
	_CREAv_i(ris, LENGTHv_i(v1) + LENGTHv_i(v2));
	for (i = 1, j = 1; j <= LENGTHv_i(v1); i++, j++)
		_ASSEGNAv_i(ris, i, _ACCEDIv_i(v1, j));
	for (j = 1; j <= LENGTHv_i(v2); i++, j++)
		_ASSEGNAv_i(ris, i, _ACCEDIv_i(v2, j));
	ris->dim = i - 1;
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _vettore3v_i(VETTOREi *ris, const VETTOREi *v1, const VETTOREi *v2, const VETTOREi *v3, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _vettore3v_i(VETTOREi *ris, const VETTOREi *v1, const VETTOREi *v2, const VETTOREi *v3)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: vettore3v_i\n", linea);
#endif
	CONTROLLA(v1 != NULL && v2 != NULL && v3 != NULL);
#ifdef FDEBUG
	_StampaVett_i(v1);
	_StampaVett_i(v2);
	_StampaVett_i(v3);
#endif
	_CREAv_i(ris, LENGTHv_i(v1) + LENGTHv_i(v2) + LENGTHv_i(v3));
	for (i = 1, j = 1; j <= LENGTHv_i(v1); i++, j++)
		_ASSEGNAv_i(ris, i, _ACCEDIv_i(v1, j));
	for (j = 1; j <= LENGTHv_i(v2); i++, j++)
		_ASSEGNAv_i(ris, i, _ACCEDIv_i(v2, j));
	for (j = 1; j <= LENGTHv_i(v3); i++, j++)
		_ASSEGNAv_i(ris, i, _ACCEDIv_i(v3, j));
	ris->dim = i - 1;
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	int  _max_v_i(const VETTOREi *v, const char *nomefile, int linea)
#else
	int  _max_v_i(const VETTOREi *v)
#endif
{
	int i;
	int max;


#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: max_v_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
#endif
	max = _ACCEDIv_i(v, 1);
	for (i = 2; i <= LENGTHv_i(v); i++) {
		if (_ACCEDIv_i(v, i) > max)
			max = _ACCEDIv_i(v, i);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\nmax = %d\n", max);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return max;
}

#ifdef MDEBUG
	VETTOREi * _accoda1_vv_i(VETTOREi *v1, const VETTOREi *v2, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _accoda1_vv_i(VETTOREi *v1, const VETTOREi *v2)
#endif
{
	int i = 1, j;
	VETTOREi *ris = NULL;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: accoda1_vv_i\n", linea);
#endif
	CONTROLLA(v1 != NULL && v2 != NULL);
#ifdef FDEBUG
	_StampaVett_i(v1);
	_StampaVett_i(v2);
#endif
	if (v1->mia_alloc >= v1->dim + v2->dim) {
		for (i = v1->dim + 1, j = 1; j <= LENGTHv_i(v2); i++, j++)
			_ASSEGNAv_i(v1, i, _ACCEDIv_i(v2, j));
		v1->dim += LENGTHv_i(v2);
	}
	else {
#ifdef FDEBUG
		_CREAstr(tmp, "");
		g_string_printf(tmp, "ATTENZIONE (accoda_vv_i, linea %s # %d): ingrandito il vettore da %d a %d!\n", nomefile, linea, LENGTHv_i(v1), LENGTHv_i(v1) + LENGTHv_i(v2));
		warning(tmp->str);
		fprintf(fp_fdbg, tmp->str);
		_CANCELLAstr(tmp);
#endif
		_CREAv_i(ris, 2 * (LENGTHv_i(v1) + LENGTHv_i(v2)));
		ris->dim = LENGTHv_i(v1) + LENGTHv_i(v2);
		for (i = 1, j = 1; j <= LENGTHv_i(v1); i++, j++)
			_ASSEGNAv_i(ris, i, _ACCEDIv_i(v1, j));
		for (j = 1; j <= LENGTHv_i(v2); i++, j++)
			_ASSEGNAv_i(ris, i, _ACCEDIv_i(v2, j));
		_CANCELLAv_i(v1);
		v1 = ris;
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(v1);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return v1;
}

#ifdef MDEBUG
	VETTOREi * _somma_vs_i(VETTOREi *ris, const VETTOREi *v, int s, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _somma_vs_i(VETTOREi *ris, const VETTOREi *v, int s)
#endif
{
	int i;


#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: somma_vs_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
	fprintf(fp_fdbg, "s = %d\n", s);
#endif
	_CREAv_i(ris, LENGTHv_i(v));
	for (i = 1; i <= LENGTHv_i(v); i++)
		_ASSEGNAv_i(ris, i, _ACCEDIv_i(v, i) + s);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	int  _esiste_v_i(int el, const VETTOREi *v, const char *nomefile, int linea)
#else
	int  _esiste_v_i(int el, const VETTOREi *v)
#endif
{
	int i, indx = 0;


#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: esiste_v_i\nel = %d\n", linea, el);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
#endif
	for (i = 1; i <= LENGTHv_i(v); i++) {
		if (UGUALE(el, _ACCEDIv_i(v, i))) {
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
	void  _elimina_indx_i(VETTOREi *v, int indx, const char *nomefile, int linea)
#else
	void  _elimina_indx_i(VETTOREi *v, int indx)
#endif
{
	int i;


#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: elimina_indx_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
	fprintf(fp_fdbg, "indx = %d\n", indx);
#endif
	for (i = indx; i <= LENGTHv_i(v) - 1; i++)
		_ASSEGNAv_i(v, i, _ACCEDIv_i(v, i + 1));
	v->dim--;
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(v);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	VETTOREi * _assegna_v_i(VETTOREi *v, int indx, int val, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _assegna_v_i(VETTOREi *v, int indx, int val)
#endif
{
	int i;
	VETTOREi *ris = NULL;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna_v_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
	fprintf(fp_fdbg, "indx = %d, val = %d\n", indx, val);
#endif
	CONTROLLA(indx > 0);
	if (indx <= LENGTHv_i(v)) {
		_ASSEGNAv_i(v, indx, val);
#ifdef FDEBUG
		fprintf(fp_fdbg, "\n\n");
		_StampaVett_i(v);
		fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	}
	else {
		_CREAv_i(ris, indx);
#ifdef FDEBUG
		_CREAstr(tmp, "");
		g_string_printf(tmp, "ATTENZIONE (assegnav_i, linea %s # %d): e` stato assegnato un elemento al di fuori dei limiti dell'array; ingrandito il vettore da %d a %d!\n", nomefile, linea, LENGTHv_i(v), indx);
		warning(tmp->str);
		fprintf(fp_fdbg, tmp->str);
		_CANCELLAstr(tmp);
#endif
		for (i = 1; i <= LENGTHv_i(v); i++)
			_ASSEGNAv_i(ris, i, _ACCEDIv_i(v, i));
#ifdef FDEBUG
		if (i < indx) {
			_CREAstr(tmp, "");
			g_string_printf(tmp, "ATTENZIONE (assegnav_i, linea %s # %d): e` stato assegnato un elemento non consecutivo e sara`/nno aggiunto/i %d valore/i NA!\n", nomefile, linea, indx - i);
			warning(tmp->str);
			fprintf(fp_fdbg, tmp->str);
			_CANCELLAstr(tmp);
		}
#endif
		for (; i < indx; i++) {
			_ASSEGNAv_i(ris, i, NA_INTEGER);
		}
		_ASSEGNAv_i(ris, indx, val);
		_CANCELLAv_i(v);
		ris->dim = indx;
		v = ris;
#ifdef FDEBUG
		fprintf(fp_fdbg, "\n\n");
		_StampaVett_i(ris);
		fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	}
	return v;
}

#ifdef MDEBUG
	VETTOREi * _somma_righe_i(VETTOREi *ris, const MATRICEi *m, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _somma_righe_i(VETTOREi *ris, const MATRICEi *m)
#endif
{
	int r, c;


#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: somma_righe_i\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
#endif
	_CREAv_i(ris, LENGTHm1_i(m));
	for (r = 1; r <= LENGTHm1_i(m); r++) {
		_ASSEGNAv_i(ris, r, 0);
		for (c = 1; c <= LENGTHm2_i(m); c++) {
			_ASSEGNAv_i(ris, r, _ACCEDIv_i(ris, r)  + _ACCEDIm_i(m, r, c));
		}
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	void  _rev1_i(const VETTOREi *v, const char *nomefile, int linea)
#else
	void  _rev1_i(const VETTOREi *v)
#endif
{
	int i;
	int tmp;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: rev_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
#endif
	for (i = 1; i <= LENGTHv_i(v) / 2; i++) {
		tmp = _ACCEDIv_i(v, i);
		_ASSEGNAv_i(v, i, _ACCEDIv_i(v, LENGTHv_i(v) - i + 1));
		_ASSEGNAv_i(v, LENGTHv_i(v) - i + 1, tmp);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(v);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	VETTOREi * _ordina_i(VETTOREi *ris, const VETTOREi *v, bool decr, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _ordina_i(VETTOREi *ris, const VETTOREi *v, bool decr)
#endif
{

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: ordina_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
	fprintf(fp_fdbg, "decr = %d\n", decr);
#endif
	_CREAv_i(ris, LENGTHv_i(v));
	ris = copia_v_i(ris, v, 1, LENGTHv_i(v));
	R_isort(ris->dati, LENGTHv_i(ris));
	if (decr)
		rev1_i(ris);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
		return ris;
}

#ifdef MDEBUG
	void  _ordina1_i(const VETTOREi *v, bool decr, const char *nomefile, int linea)
#else
	void  _ordina1_i(const VETTOREi *v, bool decr)
#endif
{
#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: ordina1_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
	fprintf(fp_fdbg, "decr = %d\n", decr);
#endif
	R_isort(v->dati, LENGTHv_i(v));
	if (decr)
		rev1_i(v);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(v);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	VETTOREi * _riga_i(VETTOREi *ris, const MATRICEi *m, int r, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _riga_i(VETTOREi *ris, const MATRICEi *m, int r)
#endif
{
	int c;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: riga_i\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "r = %d", r);
#endif
	CONTROLLA(r >= 0 && r <= m->nr);
	if (r == 0) {
	#ifdef FDEBUG
		_CREAstr(tmp, "");
		g_string_printf(tmp, "ATTENZIONE (riga_i, linea %s # %d): parametro r = 0!\n\n", nomefile, linea);
		warning(tmp->str);
		fprintf(fp_fdbg, tmp->str);
		_CANCELLAstr(tmp);
#endif
		_CREAv_i(ris, 0);
	}
	else {
		_CREAv_i(ris, LENGTHm2_i(m));
		for (c = 1; c <= LENGTHm2_i(m); c++) {
			_ASSEGNAv_i(ris, c, _ACCEDIm_i(m, r, c));
		}
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _which_v_indxle_i(VETTOREi *ris, const VETTOREi *v, int val, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_v_indxle_i(VETTOREi *ris, const VETTOREi *v, int val)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_v_indxle_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
	fprintf(fp_fdbg, "val = %d\n", val);
#endif
	_CREAv_i(ris, LENGTHv_i(v));
	for (i = 1, j = 1; i <= LENGTHv_i(v); i++) {
		if (_ACCEDIv_i(v, i) <= val) {
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
	int  _min_v_i(const VETTOREi *v, const char *nomefile, int linea)
#else
	int  _min_v_i(const VETTOREi *v)
#endif
{
	int i;
	int min;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: min_v_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
#endif
	min = _ACCEDIv_i(v, 1);
	for (i = 2; i <= LENGTHv_i(v); i++) {
		if (_ACCEDIv_i(v, i) < min)
			min = _ACCEDIv_i(v, i);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\nmin = %d\n", min);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return min;
}

#ifdef MDEBUG
	void  _somma1_vs_i(VETTOREi *v, int s, const char *nomefile, int linea)
#else
	void  _somma1_vs_i(VETTOREi *v, int s)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: somma1_vs_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
	fprintf(fp_fdbg, "s = %d\n", s);
#endif
	for (i = 1; i <= LENGTHv_i(v); i++)
		_ASSEGNAv_i(v, i, _ACCEDIv_i(v, i) + s);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(v);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	VETTOREi * _copia_v_indx_i(VETTOREi *ris, const VETTOREi *v, const VETTOREi *indx, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _copia_v_indx_i(VETTOREi *ris, const VETTOREi *v, const VETTOREi *indx)
#endif
{
	int i, j;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: copia_v_indx_i\n", linea);
#endif
	CONTROLLA(v != NULL && indx != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
	_StampaVett_i(indx);
#endif
	_CREAv_i(ris, LENGTHv_i(indx));
	for (i = 1, j = 1; i <= LENGTHv_i(indx); i++) {
		if (_ACCEDIv_i(indx, i) == 0) {
#ifdef FDEBUG
			_CREAstr(tmp, "");
			g_string_printf(tmp, "ATTENZIONE (copia_v_indx_i, linea %s # %d): il vettore degli indici contiene uno 0 alla posizione %d: saltato!\n", nomefile, linea, i);
			warning(tmp->str);
			fprintf(fp_fdbg, tmp->str);
			_CANCELLAstr(tmp);
#endif
			continue;
		}
		_ASSEGNAv_i(ris, j, _ACCEDIv_i(v, _ACCEDIv_i(indx, i)));
		j++;
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
	VETTOREi * _which_v_indxgt_i(VETTOREi *ris, const VETTOREi *v, int val, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_v_indxgt_i(VETTOREi *ris, const VETTOREi *v, int val)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_v_indxgt_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
	fprintf(fp_fdbg, "val = %d\n", val);
#endif
	_CREAv_i(ris, LENGTHv_i(v));
	for (i = 1, j = 1; i <= LENGTHv_i(v); i++) {
		if (_ACCEDIv_i(v, i) > val) {
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
	VETTOREi * _which_vv_eq_i(VETTOREi *ris, const VETTOREi *v1, const VETTOREi *v2, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_vv_eq_i(VETTOREi *ris, const VETTOREi *v1, const VETTOREi *v2)
#endif
{
	int i, j = 1;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_vv_eq_i\n", linea);
#endif
	CONTROLLA(v1 != NULL && v2 != NULL);
#ifdef FDEBUG
	_StampaVett_i(v1);
	_StampaVett_i(v2);
#endif
	_CREAv_i(ris, LENGTHv_i(v2));
	for (i = 1; i <= LENGTHv_i(v2); i++) {
		if (esiste_v_i(_ACCEDIv_i(v2, i), v1)) {
			_ASSEGNAv_i(ris, j, _ACCEDIv_i(v2, i));
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
	void  _setdiff1_i(VETTOREi *v1, const VETTOREi *v2, const char *nomefile, int linea)
#else
	void  _setdiff1_i(VETTOREi *v1, const VETTOREi *v2)
#endif
{
	int indx, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: setdiff1_i\n", linea);
#endif
	CONTROLLA(v1 != NULL && v2 != NULL);
#ifdef FDEBUG
	_StampaVett_i(v1);
	_StampaVett_i(v2);
#endif
	for (j = 1; j <= LENGTHv_i(v2); j++) {
		indx = esiste_v_i(_ACCEDIv_i(v2, j), v1);
		if (indx > 0)
			elimina_indx_i(v1, indx);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(v1);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	VETTOREi * _which_m_rowindxeq_i(VETTOREi *ris, const MATRICEi *m, int riga, int val, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_m_rowindxeq_i(VETTOREi *ris, const MATRICEi *m, int riga, int val)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_m_rowindxeq_i\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "riga = %d, val = %d\n", riga, val);
#endif
	CONTROLLA(riga > 0 && riga <= LENGTHm1_i(m));
	_CREAv_i(ris, LENGTHm1_i(m) * LENGTHm2_i(m));
	for (i = 1, j = 1; i <= LENGTHm2_i(m); i++) {
		if (UGUALE(_ACCEDIm_i(m, riga, i), val)) {
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
	void  _assegna1_vs_indx_i(const VETTOREi *v, const VETTOREi *indx, int val, const char *nomefile, int linea)
#else
	void  _assegna1_vs_indx_i(const VETTOREi *v, const VETTOREi *indx, int val)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_vs_indx_i\n", linea);
#endif
	CONTROLLA(v != NULL && indx != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
	_StampaVett_i(indx);
	fprintf(fp_fdbg, "val = %d\n", val);
#endif
	for (i = 1; i <= LENGTHv_i(indx); i++)
		_ASSEGNAv_i(v, _ACCEDIv_i(indx, i), val);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(v);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	VETTOREi * _assegna_v_indx_i(VETTOREi *ris, const VETTOREi *v, const VETTOREi *indx, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _assegna_v_indx_i(VETTOREi *ris, const VETTOREi *v, const VETTOREi *indx)
#endif
{
	int i, k = 1;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna_v_indx_i\n", linea);
#endif
	CONTROLLA(v != NULL && indx != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
	_StampaVett_i(indx);
#endif
	_CREAv_i(ris, LENGTHv_i(indx));
	for (i = 1; i <= LENGTHv_i(indx); i++) {
		_ASSEGNAv_i(ris, k, _ACCEDIv_i(v, _ACCEDIv_i(indx, i)));
		k++;
	}
	ris->dim = k - 1;
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _assegna_v_indxNA_i(VETTOREi *ris, const VETTOREi *v, const VETTOREi *indx, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _assegna_v_indxNA_i(VETTOREi *ris, const VETTOREi *v, const VETTOREi *indx)
#endif
{
	int i, k = 1, pos;
	int tmp1;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna_v_indxNA_i\n", linea);
#endif
	CONTROLLA(v != NULL && indx != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
	_StampaVett_i(indx);
#endif
	_CREAv_i(ris, LENGTHv_i(indx));
	for (i = 1; i <= LENGTHv_i(indx); i++) {
		pos = _ACCEDIv_i(indx, i);
		if (pos == NA_INTEGER) {
#ifdef FDEBUG
			_CREAstr(tmp, "");
			g_string_printf(tmp, "ATTENZIONE (assegna_vindxNA_i, linea %s # %d): saltato l'elemento specificato alla posizione %d dell'array degli indici perch� NA!\n", nomefile, linea, i);
			warning(tmp->str);
			fprintf(fp_fdbg, tmp->str);
			_CANCELLAstr(tmp);
#endif
			continue;
		}
		else {
			if (pos < 1 || pos > LENGTHv_i(v)) {
#ifdef FDEBUG
				_CREAstr(tmp, "");
				g_string_printf(tmp, "ATTENZIONE (assegna_vindxNA_i, linea %s # %d): assegnato NA per via di un elemento specificato alla posizione %d che e` al di fuori dei limiti dell'array (%d)!\n", nomefile, linea, _ACCEDIv_i(indx, i), LENGTHv_i(v));
				warning(tmp->str);
				fprintf(fp_fdbg, tmp->str);
				_CANCELLAstr(tmp);
#endif
				tmp1 = NA_INTEGER;
			}
			else
				tmp1 = _ACCEDIv_i(v, _ACCEDIv_i(indx, i));
			_ASSEGNAv_i(ris, k, tmp1);
			k++;
		}
	}
	ris->dim = k - 1;
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	int  _somma_v_indx_i(const VETTOREi *v, const VETTOREi *indx, const char *nomefile, int linea)
#else
	int  _somma_v_indx_i(const VETTOREi *v, const VETTOREi *indx)
#endif
{
	int i;
	int ris = 0;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: somma_v_indx_i\n", linea);
#endif
	CONTROLLA(v != NULL && indx != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
	_StampaVett_i(indx);
#endif
	for (i = 1; i <= LENGTHv_i(indx); i++) {
		if (_ACCEDIv_i(indx, i) >= LENGTHv_i(v))
			return NA_INTEGER;
		ris += _ACCEDIv_i(v, _ACCEDIv_i(indx, i));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\nsomma = %d\n", ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	void  _dividi1_vs_i(VETTOREd *v, double div, const char *nomefile, int linea)
#else
	void  _dividi1_vs_i(VETTOREd *v, double div)
#endif
{
	int i;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: dividi1_vs_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
	fprintf(fp_fdbg, "div = %3.3f\n", div);
#endif
	if (ISNA(div)) {
#ifdef FDEBUG
	_CREAstr(tmp, "");
	g_string_printf(tmp, "ATTENZIONE (dividi1_vs_i, linea %s # %d): divisione per NA!", nomefile, linea);
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
	g_string_printf(tmp, "ATTENZIONE (dividi1_vs_i, linea %s # %d): divisione per zero!\n", nomefile, linea);
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
	void  _assegna1_m_vv_i(const MATRICEi *m, const VETTOREi *vr, const VETTOREi *vc, int val, const char *nomefile, int linea)
#else
	void  _assegna1_m_vv_i(const MATRICEi *m, const VETTOREi *vr, const VETTOREi *vc, int val)
#endif
{
	int r, c;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_m_vv_i\n", linea);
#endif
	CONTROLLA(m != NULL && vr != NULL && vc != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	_StampaVett_i(vr);
	_StampaVett_i(vc);
	fprintf(fp_fdbg, "val = %d\n", val);
#endif
	for (r = 1; r <= LENGTHv_i(vr); r++) {
		if (_ACCEDIv_i(vr, r) == 0) {
#ifdef FDEBUG
			_CREAstr(tmp, "");
			g_string_printf(tmp, "ATTENZIONE (assegna1_m_vv_i, linea %s # %d): il vettore delle righe contiene uno 0 alla posizione %d: saltato!\n", nomefile, linea, r);
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
				g_string_printf(tmp, "ATTENZIONE (assegna1_m_vv_i, linea %s # %d): il vettore delle colonne contiene uno 0 alla posizione %d: saltato!\n", nomefile, linea, c);
				warning(tmp->str);
				fprintf(fp_fdbg, tmp->str);
				_CANCELLAstr(tmp);
#endif
				continue;
			}
			_ASSEGNAm_i(m, _ACCEDIv_i(vr, r), _ACCEDIv_i(vc, c), val);
		}
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	MATRICEi * _abs_m_i(MATRICEi *ris, const MATRICEi *m, const char *nome, const char *nomefile, int linea)
#else
	MATRICEi * _abs_m_i(MATRICEi *ris, const MATRICEi *m)
#endif
{
	int r, c;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: abs_m_i\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
#endif
	_CREAm_i(ris, LENGTHm1_i(m), LENGTHm2_i(m));

	for (r = 1; r <= LENGTHm1_i(m); r++) {
		for (c = 1; c <= LENGTHm2_i(m); c++)
			_ASSEGNAm_i(ris, r, c, abs(_ACCEDIm_i(m, r, c)));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _which_v_indxne_i(VETTOREi *ris, const VETTOREi *v, int val, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_v_indxne_i(VETTOREi *ris, const VETTOREi *v, int val)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_v_indxne_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
	fprintf(fp_fdbg, "val = %d\n", val);
#endif
	_CREAv_i(ris, LENGTHv_i(v));
	for (i = 1, j = 1; i <= LENGTHv_i(v); i++) {
		if (DIVERSO(_ACCEDIv_i(v, i), val)) {
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
	VETTOREi * _which_m_colindxeq_i(VETTOREi *ris, const MATRICEi *m, int colonna, int val, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_m_colindxeq_i(VETTOREi *ris, const MATRICEi *m, int colonna, int val)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_m_colindxeq_i\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "colonna = %d, val = %d\n", colonna, val);
#endif
	CONTROLLA(colonna > 0 && colonna <= LENGTHm2_i(m));
	_CREAv_i(ris, LENGTHm1_i(m) * LENGTHm2_i(m));
	for (i = 1, j = 1; i <= LENGTHm1_i(m); i++) {
		if (UGUALE(_ACCEDIm_i(m, i, colonna), val)) {
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
	VETTOREi * _interseca_i(VETTOREi *ris, const VETTOREi *v1, const VETTOREi *v2, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _interseca_i(VETTOREi *ris, const VETTOREi *v1, const VETTOREi *v2)
#endif
{
	int i, l, j = 1, indx, indx1;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: interseca_i\n", linea);
#endif
	CONTROLLA(v1 != NULL && v2 != NULL);
#ifdef FDEBUG
	_StampaVett_i(v1);
	_StampaVett_i(v2);
#endif
	l = min_s_i(LENGTHv_i(v1), LENGTHv_i(v2));
	_CREAv_i(ris, l);
	ris->dim = 0;
	for (i = 1; i <= LENGTHv_i(v1); i++) {
		indx = esiste_v_i(_ACCEDIv_i(v1, i), v2);
		indx1 = esiste_v_i(_ACCEDIv_i(v1, i), ris);
		if (indx > 0 && indx1 <= 0) {
			_ASSEGNAv_i(ris, j, _ACCEDIv_i(v1, i));
			ris->dim = j;
			j++;
		}
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	void  _copia1_m_riga_i(const MATRICEi *m1, int riga1, const MATRICEi *m2, int riga2, const char *nomefile, int linea)
#else
	void  _copia1_m_riga_i(const MATRICEi *m1, int riga1, const MATRICEi *m2, int riga2)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: copia1_m_riga_i\n", linea);
#endif
	CONTROLLA(m1 != NULL && m2 != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m1);
	_StampaMatr_i(m2);
	fprintf(fp_fdbg, "riga1 = %d, riga2 = %d\n", riga1, riga2);
#endif
	CONTROLLA(riga1 > 0 && riga1 <= m1->nr && riga2 > 0 && riga2 <= m2->nr);
	for (i = 1; i <= LENGTHm2_i(m1); i++)
		_ASSEGNAm_i(m1, riga1, i, _ACCEDIm_i(m2, riga2, i));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(m1);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return;
}

#ifdef MDEBUG
	MATRICEi * _aggiungi_riga_i(MATRICEi *m1, int riga1, const MATRICEi *m2, int riga2, const char *nome, const char *nomefile, int linea)
#else
	MATRICEi * _aggiungi_riga_i(MATRICEi *m1, int riga1, const MATRICEi *m2, int riga2)
#endif
{
	int i, r, c;
	MATRICEi *tmp_m = NULL;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: aggiungi_riga_i\n", linea);
#endif
	CONTROLLA(m1 != NULL && m2 != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m1);
	_StampaMatr_i(m2);
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
		g_string_printf(tmp, "ATTENZIONE (aggiungi_riga_i, linea %s # %d): ingrandite le righe della matrice da %d a %d!\n", nomefile, linea, LENGTHm1_i(m1), riga1);
		warning(tmp->str);
		fprintf(fp_fdbg, tmp->str);
		_CANCELLAstr(tmp);
#endif
		_CREAm_i(tmp_m, 2 * riga1, m1->nc);
		tmp_m->nr = riga1;
		for (r = 1; r <= LENGTHm1_i(m1); r++) {
			for (c = 1; c <= LENGTHm2_i(m1); c++)
				_ASSEGNAm_i(tmp_m, r, c, _ACCEDIm_i(m1, r, c));
		}
		_CANCELLAm_i(m1);
		m1 = tmp_m;
	}
	for (i = 1; i <= LENGTHm2_i(m1); i++)
		_ASSEGNAm_i(m1, riga1, i, _ACCEDIm_i(m2, riga2, i));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(m1);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return m1;
}

#ifdef MDEBUG
	void  _assegna1_ms_rigaindx_i(const MATRICEi *m, int riga, const VETTOREi *indx, int val, const char *nomefile, int linea)
#else
	void  _assegna1_ms_rigaindx_i(const MATRICEi *m, int riga, const VETTOREi *indx, int val)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_ms_rigaindx_i\n", linea);
#endif
	CONTROLLA(m != NULL && indx != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "riga = %d\n", riga);
	_StampaVett_i(indx);
	fprintf(fp_fdbg, "val = %d\n", val);
#endif
	CONTROLLA(riga > 0 && indx != NULL);
	for (i = 1; i <= LENGTHv_i(indx); i++) {
		_ASSEGNAm_i(m, riga, _ACCEDIv_i(indx, i), val);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return;
}

#ifdef MDEBUG
	MATRICEi * _aggiungi_ms_rigaindx_i(MATRICEi *m, int riga, const VETTOREi *indx, int val, const char *nome, const char *nomefile, int linea)
#else
	MATRICEi * _aggiungi_ms_rigaindx_i(MATRICEi *m, int riga, const VETTOREi *indx, int val)
#endif
{
	int i, r, c;
	MATRICEi *tmp_m = NULL;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: aggiungi_ms_rigaindx_i\n", linea);
#endif
	CONTROLLA(m != NULL && indx != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "riga = %d\n", riga);
	_StampaVett_i(indx);
	fprintf(fp_fdbg, "val = %d\n", val);
#endif
	if (riga > m->nr + 1) {
		Rprintf("la riga da aggiungere non e` minore o uguale al numero di righe attuali + 1: R riempirebbe quelle in mezzo con NA, e` questa l'intenzione?\n");
		error("");
	}
	if (riga > m->alloc_r) {
#ifdef FDEBUG
		_CREAstr(tmp, "");
		g_string_printf(tmp, "ATTENZIONE (aggiungi_riga_indx_i, linea %s # %d): ingrandite le righe della matrice da %d a %d!\n", nomefile, linea, LENGTHm1_i(m), riga);
		warning(tmp->str);
		fprintf(fp_fdbg, tmp->str);
		_CANCELLAstr(tmp);
#endif
		_CREAm_i(tmp_m, 2 * riga, m->nc);
		tmp_m->nr = riga;
		for (r = 1; r <= LENGTHm1_i(m); r++) {
			for (c = 1; c <= LENGTHm2_i(m); c++)
				_ASSEGNAm_i(tmp_m, r, c, _ACCEDIm_i(m, r, c));
		}
		_CANCELLAm_i(m);
		m = tmp_m;
	}
	for (i = 1; i <= LENGTHv_i(indx); i++) {
		_ASSEGNAm_i(m, riga, _ACCEDIv_i(indx, i), val);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return m;
}

#ifdef MDEBUG
	void  _assegna1_mm_rigaindx_i(const MATRICEi *m1, int riga, const MATRICEi *m2, const VETTOREi *indx, const char *nomefile, int linea)
#else
	void  _assegna1_mm_rigaindx_i(const MATRICEi *m1, int riga, const MATRICEi *m2, const VETTOREi *indx)
#endif
{
	int i;
	int tmp;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_mm_rigaindx_i\n", linea);
#endif
	CONTROLLA(m1 != NULL && m2 != NULL && indx != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m1);
	fprintf(fp_fdbg, "riga = %d\n", riga);
	_StampaMatr_i(m2);
	_StampaVett_i(indx);
#endif
	CONTROLLA(riga > 0 && riga <= m1->nr && indx->dim == m1->nc);
	for (i = 1; i <= LENGTHv_i(indx); i++) {
		if (_ACCEDIv_i(indx, i) > LENGTHm1_i(m2) * LENGTHm2_i(m2))
			tmp = NA_INTEGER;
		else
			tmp = m2->dati[_ACCEDIv_i(indx, i) - 1];
		_ASSEGNAm_i(m1, riga, i, tmp);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(m1);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return;
}

#ifdef MDEBUG
	void  _assegna1_ms_riga_i(const MATRICEi *m, int riga, int val, const char *nomefile, int linea)
#else
	void  _assegna1_ms_riga_i(const MATRICEi *m, int riga, int val)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_ms_riga_i\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "riga = %d\n", riga);
	fprintf(fp_fdbg, "val = %d\n", val);
#endif
	CONTROLLA(riga > 0 && riga <= m->nr);
	for (i = 1; i <= LENGTHm2_i(m); i++)
		_ASSEGNAm_i(m, riga, i, val);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return;
}

#ifdef MDEBUG
	void  _assegna1_m_indxlt_i(const MATRICEi *m, int val1, int val2, const char *nomefile, int linea)
#else
	void  _assegna1_m_indxlt_i(const MATRICEi *m, int val1, int val2)
#endif
{
	int r, c;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_mindxlt_i\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "val1 = %d\n", val1);
	fprintf(fp_fdbg, "val2 = %d\n", val2);
#endif
	for (r = 1; r <= LENGTHm1_i(m); r++) {
		for (c = 1; c <= LENGTHm2_i(m); c++) {
			if (_ACCEDIm_i(m, r, c) < val1)
				_ASSEGNAm_i(m, r, c, val2);
		}
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return;
}

#ifdef MDEBUG
	void  _assegna1_m_indxgt_i(const MATRICEi *m, int val1, int val2, const char *nomefile, int linea)
#else
	void  _assegna1_m_indxgt_i(const MATRICEi *m, int val1, int val2)
#endif
{
	int r, c;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_m_indxgt_i\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "val1 = %d\n", val1);
	fprintf(fp_fdbg, "val2 = %d\n", val2);
#endif
	for (r = 1; r <= LENGTHm1_i(m); r++) {
		for (c = 1; c <= LENGTHm2_i(m); c++) {
			if (_ACCEDIm_i(m, r, c) > val1)
				_ASSEGNAm_i(m, r, c, val2);
		}
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return;
}

#ifdef MDEBUG
	int  _copia_m_colindx_i(const MATRICEi *m, const VETTOREi *indx, int colonna, const char *nomefile, int linea)
#else
	int  _copia_m_colindx_i(const MATRICEi *m, const VETTOREi *indx, int colonna)
#endif
{
	int ris;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: copia_m_colindx_i\n", linea);
#endif
	CONTROLLA(m != NULL && indx != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	_StampaVett_i(indx);
	fprintf(fp_fdbg, "colonna = %d\n", colonna);
#endif
	CONTROLLA(colonna > 0 && colonna <= m->nc);
	ris = _ACCEDIm_i(m, _ACCEDIv_i(indx, colonna), colonna);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\nris: %d\n", ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	void  _assegna1_mv_riga_i(const MATRICEi *m, int riga, const VETTOREi *v, const char *nomefile, int linea)
#else
	void  _assegna1_mv_riga_i(const MATRICEi *m, int riga, const VETTOREi *v)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_mv_riga_i\n", linea);
#endif
	CONTROLLA(m != NULL && v != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "riga = %d\n", riga);
	_StampaVett_i(v);
#endif
	CONTROLLA(riga > 0 && riga <= m->nr);
	for (i = 1; i <= LENGTHv_i(v); i++)
		_ASSEGNAm_i(m, riga, i, _ACCEDIv_i(v, i));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return;
}

#ifdef MDEBUG
	MATRICEi * _aggiungi_mv_riga_i(MATRICEi *m, int riga, const VETTOREi *v, const char *nome, const char *nomefile, int linea)
#else
	MATRICEi * _aggiungi_mv_riga_i(MATRICEi *m, int riga, const VETTOREi *v)
#endif
{
	int i, r, c;
	MATRICEi *tmp_m = NULL;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: aggiungi_mv_riga_i\n", linea);
#endif
	CONTROLLA(m != NULL && v != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "riga = %d\n", riga);
	_StampaVett_i(v);
#endif
	CONTROLLA(m->nc == v->dim);
	if (riga > m->nr + 1) {
		Rprintf("la riga da aggiungere non e` minore o uguale al numero di righe attuali + 1: R riempirebbe quelle in mezzo con NA, e` questa l'intenzione?\n");
		error("");
	}
	if (riga > LENGTHm1_i(m)) {
		_CREAm_i(tmp_m, riga, m->nc);
		tmp_m->nr = riga;
		for (r = 1; r <= LENGTHm1_i(m); r++) {
			for (c = 1; c <= LENGTHm2_i(m); c++)
				_ASSEGNAm_i(tmp_m, r, c, _ACCEDIm_i(m, r, c));
		}
		_CANCELLAm_i(m);
		m = tmp_m;
	}
	for (i = 1; i <= LENGTHm2_i(m); i++)
		_ASSEGNAm_i(m, riga, i, _ACCEDIv_i(v, i));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return m;
}

#ifdef MDEBUG
	MATRICEi * _aggiungi_ms_riga_i(MATRICEi *m, int riga, int val, const char *nome, const char *nomefile, int linea)
#else
	MATRICEi * _aggiungi_ms_riga_i(MATRICEi *m, int riga, int val)
#endif
{
	int i, r, c;
	MATRICEi *tmp_m = NULL;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: aggiungi_ms_riga_i\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "riga = %d, val = %d\n", riga, val);
#endif
	if (riga > m->nr + 1) {
		Rprintf("la riga da aggiungere non e` minore o uguale al numero di righe attuali + 1: R riempirebbe quelle in mezzo con NA, e` questa l'intenzione?\n");
		error("");
	}
	if (riga > m->alloc_r) {
#ifdef FDEBUG
		_CREAstr(tmp, "");
		g_string_printf(tmp, "ATTENZIONE (aggiungi_ms_riga_i, linea %s # %d): ingrandite le righe della matrice da %d a %d\n", nomefile, linea, LENGTHm1_i(m), riga);
		warning(tmp->str);
		fprintf(fp_fdbg, tmp->str);
		_CANCELLAstr(tmp);
#endif
		_CREAm_i(tmp_m, 2 * riga, m->nc);
		tmp_m->nr = riga;
		for (r = 1; r <= LENGTHm1_i(m); r++) {
			for (c = 1; c <= LENGTHm2_i(m); c++)
				_ASSEGNAm_i(tmp_m, r, c, _ACCEDIm_i(m, r, c));
		}
		_CANCELLAm_i(m);
		m = tmp_m;
	}
	for (i = 1; i <= LENGTHm2_i(m); i++)
		_ASSEGNAm_i(m, riga, i, val);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return m;
}

#ifdef MDEBUG
	VETTOREi * _somma_colonne_i(VETTOREi *ris, const MATRICEi *m, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _somma_colonne_i(VETTOREi *ris, const MATRICEi *m)
#endif
{
	int r, c;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: somma_colonne_i\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
#endif
	_CREAv_i(ris, LENGTHm2_i(m));
	for (c = 1; c <= LENGTHm2_i(m); c++) {
		_ASSEGNAv_i(ris, c, 0);
		for (r = 1; r <= LENGTHm1_i(m); r++) {
			_ASSEGNAv_i(ris, c, _ACCEDIv_i(ris, c)  + _ACCEDIm_i(m, r, c));
		}
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _which_v_indxlt_i(VETTOREi *ris, const VETTOREi *v, int val, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_v_indxlt_i(VETTOREi *ris, const VETTOREi *v, int val)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_v_indxlt_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
	fprintf(fp_fdbg, "val = %d", val);
#endif
	_CREAv_i(ris, LENGTHv_i(v));
	for (i = 1, j = 1; i <= LENGTHv_i(v); i++) {
		if (_ACCEDIv_i(v, i) < val) {
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
	VETTOREi * _unione1_i(VETTOREi *v1, const VETTOREi *v2, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _unione1_i(VETTOREi *v1, const VETTOREi *v2)
#endif
{
	int i, j;
	VETTOREi *tmp_ris = NULL, *ris = NULL;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: unione1_i\n", linea);
#endif
	CONTROLLA(v1 != NULL && v2 != NULL);
#ifdef FDEBUG
	_StampaVett_i(v1);
	_StampaVett_i(v2);
#endif
		_CREAv_i(tmp_ris, LENGTHv_i(v1) + LENGTHv_i(v2));
		for (i = 1; i <= LENGTHv_i(v1); i++)
			_ASSEGNAv_i(tmp_ris, i, _ACCEDIv_i(v1, i));
		for (j = 1; j <= LENGTHv_i(v2); j++) {
			_ASSEGNAv_i(tmp_ris, i, _ACCEDIv_i(v2, j));
			i++;
		}
		_CANCELLAv_i(v1);
		tmp_ris->dim = i - 1;
		ris = elimina_doppi_i(ris, tmp_ris);
		_CANCELLAv_i(tmp_ris);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	int  _somma_v_i(const VETTOREi *v, bool canc_NA, const char *nomefile, int linea)
#else
	int  _somma_v_i(const VETTOREi *v, bool canc_NA)
#endif
{
	int i;
	int ris = 0;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: somma_v_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
	fprintf(fp_fdbg, "canc_NA = %d\n\n", canc_NA);
#endif
	for (i = 1; i <= LENGTHv_i(v); i++) {
		if (!canc_NA && (_ACCEDIv_i(v, i) == NA_INTEGER)) {
			ris = NA_INTEGER;
			break;
		}
		else if (!(_ACCEDIv_i(v, i) == NA_INTEGER))
			ris += _ACCEDIv_i(v, i);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\nris: %d\n", ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	int  _somma_m_i(const MATRICEi *m, const char *nomefile, int linea)
#else
	int  _somma_m_i(const MATRICEi *m)
#endif
{
	int r, c;
	int ris = 0;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: somma_m_i\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
#endif
	for (r = 1; r <= LENGTHm1_i(m); r++) {
		for (c = 1; c <= LENGTHm2_i(m); c++)
			ris += _ACCEDIm_i(m, r, c);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\nris: %d\n", ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _setdiff_i(VETTOREi *ris, const VETTOREi *v1, const VETTOREi *v2, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _setdiff_i(VETTOREi *ris, const VETTOREi *v1, const VETTOREi *v2)
#endif
{
	int indx, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: setdiff_i\n", linea);
#endif
	CONTROLLA(v1 != NULL && v2 != NULL);
#ifdef FDEBUG
	_StampaVett_i(v1);
	_StampaVett_i(v2);
#endif
	_CREAv_i(ris, LENGTHv_i(v1));
	ris = copia_v_i(ris, v1, 1, LENGTHv_i(v1));
	for (j = 1; j <= LENGTHv_i(v2); j++) {
		indx = esiste_v_i(_ACCEDIv_i(v2, j), ris);
		if (indx > 0)
			elimina_indx_i(ris, indx);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _diff_vv_i(VETTOREi *ris, const VETTOREi *v1, const VETTOREi *v2, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _diff_vv_i(VETTOREi *ris, const VETTOREi *v1, const VETTOREi *v2)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: diff_vv_i\n", linea);
#endif
	CONTROLLA(v1 != NULL && v2 != NULL);
#ifdef FDEBUG
	_StampaVett_i(v1);
	_StampaVett_i(v2);
#endif
	_CREAv_i(ris, LENGTHv_i(v1));
	for (i = 1; i <= LENGTHv_i(v1); i++)
		_ASSEGNAv_i(ris, i, _ACCEDIv_i(v1, i) - _ACCEDIv_i(v2, i));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _abs_v_i(VETTOREi *ris, const VETTOREi *v, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _abs_v_i(VETTOREi *ris, const VETTOREi *v)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: abs_v_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
#endif
	_CREAv_i(ris, LENGTHv_i(v));
	for (i = 1; i <= LENGTHv_i(v); i++)
		_ASSEGNAv_i(ris, i, abs(_ACCEDIv_i(v, i)));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _which_v_indxna_i(VETTOREi *ris, const VETTOREi *v, bool complemento, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_v_indxna_i(VETTOREi *ris, const VETTOREi *v, bool complemento)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_v_indxna_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
	fprintf(fp_fdbg, "complemento = %d\n", complemento);
#endif
	_CREAv_i(ris, LENGTHv_i(v));
	for (i = 1, j = 1; i <= LENGTHv_i(v); i++) {
		if ((_ACCEDIv_i(v, i) == NA_INTEGER) && !complemento) {
			_ASSEGNAv_i(ris, j, i);
			j++;
		}
		else if (!(_ACCEDIv_i(v, i) == NA_INTEGER) && complemento) {
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
	VETTOREi * _which_v_indxeq_i(VETTOREi *ris, const VETTOREi *v, int val, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_v_indxeq_i(VETTOREi *ris, const VETTOREi *v, int val)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_v_indxeq_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
	fprintf(fp_fdbg, "val = %d\n", val);
#endif
	_CREAv_i(ris, LENGTHv_i(v));
	for (i = 1, j = 1; i <= LENGTHv_i(v); i++) {
#ifdef VERSIONE_d
		if (val == R_PosInf && _ACCEDIv_i(v, i) == R_PosInf) {
			_ASSEGNAv_i(ris, j, i);
			j++;
		}
		else if (val == R_NegInf && _ACCEDIv_i(v, i) == R_NegInf) {
			_ASSEGNAv_i(ris, j, i);
			j++;
		} else
#endif
		if (UGUALE(_ACCEDIv_i(v, i), val)) {
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
	MATRICEi * _cbind2v_i(MATRICEi *ris, const VETTOREi *v1, const VETTOREi *v2, const char *nome, const char *nomefile, int linea)
#else
	MATRICEi * _cbind2v_i(MATRICEi *ris, const VETTOREi *v1, const VETTOREi *v2)
#endif
{
	int lM, r, r1, r2;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: cbind2v_i\n", linea);
#endif
	CONTROLLA(v1 != NULL && v2 != NULL);
#ifdef FDEBUG
	_StampaVett_i(v1);
	_StampaVett_i(v2);
#endif
	lM = max_s_i(LENGTHv_i(v1), LENGTHv_i(v2));
	_CREAm_i(ris, lM, 2);
	for (r = 1, r1 = 1, r2 = 1; r <= lM; r++, r1++, r2++) {
		if (r > LENGTHv_i(v1))
			r1 = 1;
		_ASSEGNAm_i(ris, r, 1, _ACCEDIv_i(v1, r1));
		if (r > LENGTHv_i(v2))
			r2 = 1;
		_ASSEGNAm_i(ris, r, 2, _ACCEDIv_i(v2, r2));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _max_righe_i(VETTOREi *ris, const MATRICEi *m, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _max_righe_i(VETTOREi *ris, const MATRICEi *m)
#endif
{
	int r, c;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: max_righe_i\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
#endif
	_CREAv_i(ris, LENGTHm1_i(m));
	for (r = 1; r <= LENGTHm1_i(m); r++) {
		_ASSEGNAv_i(ris, r, _ACCEDIm_i(m, r, 1));
		for (c = 2; c <= LENGTHm2_i(m); c++) {
			if (_ACCEDIv_i(ris, r) < _ACCEDIm_i(m, r, c))
				_ASSEGNAv_i(ris, r, _ACCEDIm_i(m, r, c));
		}
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _f_aux_i(VETTOREd *ris, const VETTOREi *a, const VETTOREi *b, const VETTOREi *m, const VETTOREi *t, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _f_aux_i(VETTOREd *ris, const VETTOREi *a, const VETTOREi *b, const VETTOREi *m, const VETTOREi *t)
#endif
{
	int i;
	double tmp1;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: f_aux_i\n", linea);
#endif
	CONTROLLA(a != NULL && b != NULL && m != NULL && t != NULL);
#ifdef FDEBUG
	_StampaVett_i(a);
	_StampaVett_i(b);
	_StampaVett_i(m);
	_StampaVett_i(t);
#endif
	_CREAv_d(ris, LENGTHv_i(a));
	for (i = 1; i <= LENGTHv_i(a); i++) {
		if (UGUALE(_ACCEDIv_i(t, i), 0)) {
#ifdef FDEBUG
			_CREAstr(tmp, "");
			g_string_printf(tmp, "ATTENZIONE (f_aux_i, linea %s # %d): l'elemento %d ha provocato una divisione per zero e gli e` stato assegnato un valore al di fuori del dominio!\n", nomefile, linea, i);
			warning(tmp->str);
			fprintf(fp_fdbg, tmp->str);
			_CANCELLAstr(tmp);
#endif
			if (Segno((_ACCEDIv_i(a, i) - _ACCEDIv_i(b, i)) * _ACCEDIv_i(m, i)) > 0)
				tmp1 = R_PosInf;
			else if (Segno((_ACCEDIv_i(a, i) - _ACCEDIv_i(b, i)) * _ACCEDIv_i(m, i)) < 0)
				tmp1 = R_NegInf;
			else
				tmp1 = R_NaN;
		}
		else
			tmp1 = (double) Segno(_ACCEDIv_i(a, i) - _ACCEDIv_i(b, i)) * _ACCEDIv_i(m, i) / (double) _ACCEDIv_i(t, i);
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
	VETTOREi * _somma_vv_i(VETTOREi *ris, const VETTOREi *v1, const VETTOREi *v2, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _somma_vv_i(VETTOREi *ris, const VETTOREi *v1, const VETTOREi *v2)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: somma_vv_i\n", linea);
#endif
	CONTROLLA(v1 != NULL && v2 != NULL);
#ifdef FDEBUG
	_StampaVett_i(v1);
	_StampaVett_i(v2);
#endif
	_CREAv_i(ris, LENGTHv_i(v1));
	for (i = 1; i <= LENGTHv_i(v1); i++)
		_ASSEGNAv_i(ris, i, _ACCEDIv_i(v1, i) + _ACCEDIv_i(v2, i));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	MATRICEi * _trasponi_i(MATRICEi *ris, const MATRICEi *m, const char *nome, const char *nomefile, int linea)
#else
	MATRICEi * _trasponi_i(MATRICEi *ris, const MATRICEi *m)
#endif
{
	int r, c;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: trasponi_i\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
#endif
	_CREAm_i(ris, LENGTHm2_i(m), LENGTHm1_i(m));

	for (r = 1; r <= LENGTHm1_i(m); r++) {
		for (c = 1; c <= LENGTHm2_i(m); c++)
			_ASSEGNAm_i(ris, c, r, _ACCEDIm_i(m, r, c));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	void  _somma1_m_i(MATRICEi *m1, const MATRICEi *m2, const char *nomefile, int linea)
#else
	void  _somma1_m_i(MATRICEi *m1, const MATRICEi *m2)
#endif
{
	int r, c;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: somma1_m_i\n", linea);
#endif
	CONTROLLA(m1 != NULL && m2 != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m1);
	_StampaMatr_i(m2);
#endif
	CONTROLLA(m1->nr == m2->nr && m1->nc == m2->nc);
	for (r = 1; r <= LENGTHm1_i(m1); r++) {
		for (c = 1; c <= LENGTHm2_i(m1); c++)
			_ASSEGNAm_i(m1, r, c, _ACCEDIm_i(m1, r, c) + _ACCEDIm_i(m2, r, c));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(m1);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	void  _assegna1_ms_indx_i(const MATRICEi *m, const VETTOREi *indx, int val, const char *nomefile, int linea)
#else
	void  _assegna1_ms_indx_i(const MATRICEi *m, const VETTOREi *indx, int val)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_ms_indx_i\n", linea);
#endif
	CONTROLLA(m != NULL && indx != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	_StampaVett_i(indx);
	fprintf(fp_fdbg, "val = %d\n", val);
#endif
	for (i = 1; i <= LENGTHv_i(indx); i++)
		_ASSEGNAmv_i(m, _ACCEDIv_i(indx, i), val);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	VETTOREi * _which_m_indxne_i(VETTOREi *ris, const MATRICEi *m, int val, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_m_indxne_i(VETTOREi *ris, const MATRICEi *m, int val)
#endif
{
	int r, c, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_m_indxne_i\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "val = %d\n", val);
#endif
	_CREAv_i(ris, LENGTHm1_i(m) * LENGTHm2_i(m));
	j = 1;
	for (c = 1; c <= LENGTHm2_i(m); c++) {
		for (r = 1; r <= LENGTHm1_i(m); r++) {
			if (DIVERSO(_ACCEDIm_i(m, r, c), val)) {
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
	void  _assegna1_s_diag_i(const MATRICEi *m, int val, const char *nomefile, int linea)
#else
	void  _assegna1_s_diag_i(const MATRICEi *m, int val)
#endif
{
	int l, d;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_s_diag_i\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "val = %d\n", val);
#endif
	l = min_s_i(LENGTHm1_i(m), LENGTHm2_i(m));
	for (d = 1; d <= l; d++)
		_ASSEGNAm_i(m, d, d, val);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	double  _media_v_i(VETTOREi *v, const char *nomefile, int linea)
#else
	double  _media_v_i(VETTOREi *v)
#endif
{
	int i, j = 0;
	double ris = 0.0;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: media_v_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
#endif
	for (i = 1; i <= LENGTHv_i(v); i++) {
		if (!(_ACCEDIv_i(v, i) == NA_INTEGER)) {
			ris += _ACCEDIv_i(v, i);
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
	int  _somma_riga_i(const MATRICEi *m, int riga, const char *nomefile, int linea)
#else
	int  _somma_riga_i(const MATRICEi *m, int riga)
#endif
{
	int c;
	int ris = 0;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: somma_riga_i\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "riga = %d\n", riga);
#endif
	CONTROLLA(riga > 0 && riga <= m->nr);
	for (c = 1; c <= LENGTHm2_i(m); c++)
		ris += _ACCEDIm_i(m, riga, c);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\nris: %d\n", ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _ordine_i(VETTOREi *ris, const VETTOREi *v, bool decr, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _ordine_i(VETTOREi *ris, const VETTOREi *v, bool decr)
#endif
{
	int i;
	VETTOREd *tmp = NULL;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: ordine_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
	fprintf(fp_fdbg, "decr = %d\n", decr);
#endif
	_CREAv_i(ris, LENGTHv_i(v));
	for (i = 1; i <= LENGTHv_i(v); i++)
		_ASSEGNAv_i(ris, i, i);
#ifdef VERSIONE_d
	tmp = copia_v_d(tmp, v, 1, LENGTHv_i(v));
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
	VETTOREi * _diag_i(VETTOREi *ris, const MATRICEi *m, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _diag_i(VETTOREi *ris, const MATRICEi *m)
#endif
{
	int l, d;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: diag_i\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
#endif
	l = min_s_i(LENGTHm1_i(m), LENGTHm2_i(m));
	_CREAv_i(ris, l);
	for (d = 1; d <= l; d++)
		_ASSEGNAv_i(ris, d, _ACCEDIm_i(m, d, d));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _which_m_colindxne_i(VETTOREi *ris, const MATRICEi *m, int colonna, int val, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_m_colindxne_i(VETTOREi *ris, const MATRICEi *m, int colonna, int val)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_m_colindxne_i\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "colonna = %d, val = %d\n", colonna, val);
#endif
	CONTROLLA(colonna > 0 && colonna <= LENGTHm2_i(m));
	_CREAv_i(ris, LENGTHm1_i(m) * LENGTHm2_i(m));
	for (i = 1, j = 1; i <= LENGTHm1_i(m); i++) {
		if (DIVERSO(_ACCEDIm_i(m, i, colonna), val)) {
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
	VETTOREi * _which_m_colindxin_i(VETTOREi *ris, const MATRICEi *m, int c, const VETTOREi *v, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_m_colindxin_i(VETTOREi *ris, const MATRICEi *m, int c, const VETTOREi *v)
#endif
{
	int i, j, indx;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_m_colindxin_i\n", linea);
#endif
	CONTROLLA(m != NULL && v != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "c = %d\n", c);
	_StampaVett_i(v);
#endif
	CONTROLLA(c > 0 && c <= LENGTHm2_i(m));
	_CREAv_i(ris, LENGTHm1_i(m));
	for (i = 1, j = 1; i <= LENGTHm1_i(m); i++) {
		indx = esiste_v_i(_ACCEDIm_i(m, i, c), v);
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
	MATRICEi * _copia_m_ncol_i(MATRICEi *ris, const MATRICEi *m, int colonne, const char *nome, const char *nomefile, int linea)
#else
	MATRICEi * _copia_m_ncol_i(MATRICEi *ris, const MATRICEi *m, int colonne)
#endif
{
	int r, c, l, k;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: copia_m_ncol_i\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
#endif
	l = LENGTHm1_i(m) * LENGTHm2_i(m);
#ifdef FDEBUG
	if ((l > colonne && l % colonne) || (l < colonne && colonne % l)) {
		_CREAstr(tmp, "");
		g_string_printf(tmp, "ATTENZIONE (copia_m_ncol_i, linea %s # %d): il numero di elementi di m (%d) non e` un multiplo o sotto-multiplo del numero di colonne (%d)!\n", nomefile, linea, l, colonne);
		warning(tmp->str);
		fprintf(fp_fdbg, tmp->str);
		_CANCELLAstr(tmp);
	}
#endif
	if (l >= colonne)
		l /= colonne;
	else
		l = 1;
	_CREAm_i(ris, l, colonne);
	k = 0;
	for (c = 1; c <= colonne; c++) {
		for (r = 1; r <= l; r++)
			_ASSEGNAm_i(ris, r, c, _ACCEDImv_i(m, k + 1));
			k = (k % l) + 1;
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _which_m_colneand2_i(VETTOREi *ris, const MATRICEi *m, int c, int val1, int val2, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_m_colneand2_i(VETTOREi *ris, const MATRICEi *m, int c, int val1, int val2)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_m_colneand2_i\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "c = %d, val1 = %d, val2 = %d\n", c, val1, val2);
#endif
	CONTROLLA(c > 0 && c <= LENGTHm2_i(m));
	_CREAv_i(ris, LENGTHm1_i(m) * LENGTHm2_i(m));
	for (i = 1, j = 1; i <= LENGTHm1_i(m); i++) {
		if (DIVERSO(_ACCEDIm_i(m, i, c), val1) && DIVERSO(_ACCEDIm_i(m, i, c), val2)) {
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
	MATRICEi * _righe_i(MATRICEi *ris, const MATRICEi *m, const VETTOREi *indx, const char *nome, const char *nomefile, int linea)
#else
	MATRICEi * _righe_i(MATRICEi *ris, const MATRICEi *m, const VETTOREi *indx)
#endif
{
	int r, c;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: righe_i\n", linea);
#endif
	CONTROLLA(m != NULL && indx != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	_StampaVett_i(indx);
	fprintf(fp_fdbg, "\n");
#endif
	_CREAm_i(ris, LENGTHv_i(indx), LENGTHm2_i(m));
	for (r = 1; r <= LENGTHv_i(indx); r++) {
		for (c = 1; c <= LENGTHm2_i(m); c++)
			_ASSEGNAm_i(ris, r, c, _ACCEDIm_i(m, _ACCEDIv_i(indx, r), c));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	void  _assegna1_v_diag_i(const MATRICEi *m, const VETTOREi *v, const char *nomefile, int linea)
#else
	void  _assegna1_v_diag_i(const MATRICEi *m, const VETTOREi *v)
#endif
{
	int l, d;

	l = min_s_i(LENGTHm1_i(m), LENGTHm2_i(m));
#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_v_diag_i\n", linea);
#endif
	CONTROLLA(m != NULL && v != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	_StampaVett_i(v);
#endif
	CONTROLLA(l == LENGTHv_i(v));
	for (d = 1; d <= l; d++)
		_ASSEGNAm_i(m, d, d, _ACCEDIv_i(v, d));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	void  _assegna1_ms_indx2_i(const MATRICEi *m, const MATRICEi *indxm, int val, const char *nomefile, int linea)
#else
	void  _assegna1_ms_indx2_i(const MATRICEi *m, const MATRICEi *indxm, int val)
#endif
{
	int i, l;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_ms_indx2_i\n", linea);
#endif
	CONTROLLA(m != NULL && indxm != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	_StampaMatr_i(indxm);
	fprintf(fp_fdbg, "val = %d\n", val);
#endif
	CONTROLLA(LENGTHm2_i(indxm) == 2);
	l = min_s_i(LENGTHm1_i(m), LENGTHm1_i(indxm));
	for (i = 1; i <= l; i++)
		_ASSEGNAm_i(m, _ACCEDIm_i(indxm, i, 1), _ACCEDIm_i(indxm, i, 2), val);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	void  _incr1_v_i(const VETTOREi *v, int s, const char *nomefile, int linea)
#else
	void  _incr1_v_i(const VETTOREi *v, int s)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: incr1_v_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
	fprintf(fp_fdbg, "s = %d\n", s);
#endif
	for (i = 1; i <= LENGTHv_i(v); i++)
		_ASSEGNAv_i(v, i, _ACCEDIv_i(v, i) + s);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(v);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	MATRICEi * _copia_m_i(MATRICEi *ris, const MATRICEi *da, const char *nome, const char *nomefile, int linea)
#else
	MATRICEi * _copia_m_i(MATRICEi *ris, const MATRICEi *da)
#endif
{
	int r, c;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: copia_m_i\n", linea);
#endif
	CONTROLLA(da != NULL);
#ifdef FDEBUG
	_StampaMatr_i(da);
#endif
	_CREAm_i(ris, LENGTHm1_i(da), LENGTHm2_i(da));
	for (r = 1; r <= LENGTHm1_i(da); r++) {
		for (c = 1; c <= LENGTHm2_i(da); c++)
			_ASSEGNAm_i(ris, r, c, _ACCEDIm_i(da, r, c));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	MATRICEi * _somma_mm_i(MATRICEi *ris, const MATRICEi *m1, const MATRICEi *m2, const char *nome, const char *nomefile, int linea)
#else
	MATRICEi * _somma_mm_i(MATRICEi *ris, const MATRICEi *m1, const MATRICEi *m2)
#endif
{
	int r, c;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: somma_mm_i\n", linea);
#endif
	CONTROLLA(m1 != NULL && m2 != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m1);
	_StampaMatr_i(m2);
#endif
	CONTROLLA(m1->nr == m2->nr && m1->nc == m2->nc);
	_CREAm_i(ris, LENGTHm1_i(m1), LENGTHm2_i(m1));

	for (r = 1; r <= LENGTHm1_i(m1); r++) {
		for (c = 1; c <= LENGTHm2_i(m1); c++)
			_ASSEGNAm_i(ris, r, c, _ACCEDIm_i(m1, r, c) + _ACCEDIm_i(m2, r, c));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _accoda1_vs_i(VETTOREi *v, int s, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _accoda1_vs_i(VETTOREi *v, int s)
#endif
{
	int i = 1, j;
	VETTOREi *tmp_v = NULL;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: accoda1_vs_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
	fprintf(fp_fdbg, "s: %d\n", s);
#endif
	if (v->mia_alloc >= v->dim + 1) {
		_ASSEGNAv_i(v, v->dim + 1, s);
		v->dim++;
	}
	else {
#ifdef FDEBUG
	_CREAstr(tmp, "");
	g_string_printf(tmp, "ATTENZIONE (accoda1_vs_i, linea %s # %d): ingrandito il vettore da %d a %d\n", nomefile, linea, LENGTHv_i(v), LENGTHv_i(v) + 1);
	warning(tmp->str);
	fprintf(fp_fdbg, tmp->str);
	_CANCELLAstr(tmp);
#endif
		_CREAv_i(tmp_v, 2 * (LENGTHv_i(v) + 1));
		tmp_v->dim = LENGTHv_i(v) + 1;
		for (i = 1, j = 1; j <= LENGTHv_i(v); i++, j++)
			_ASSEGNAv_i(tmp_v, i, _ACCEDIv_i(v, j));
		_ASSEGNAv_i(tmp_v, i, s);
		_CANCELLAv_i(v);
		v = tmp_v;
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(v);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return v;
}

#ifdef MDEBUG
	void  _segmento1_v_i(const VETTOREi *v, int st, int end, const char *nomefile, int linea)
#else
	void  _segmento1_v_i(const VETTOREi *v, int st, int end)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: segmento1_v_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
	fprintf(fp_fdbg, "start = %d\n", st);
	fprintf(fp_fdbg, "end = %d\n", end);
#endif
	CONTROLLA(st > 0 && end <= LENGTHv_i(v));
	for (i = st, j = 1; i <= end; i++, j++)
		_ASSEGNAv_i(v, j, _ACCEDIv_i(v, i));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(v);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	VETTOREi * _segmento_v_i(VETTOREi *ris, const VETTOREi *v, int st, int end, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _segmento_v_i(VETTOREi *ris, const VETTOREi *v, int st, int end)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: segmento_v_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
	fprintf(fp_fdbg, "start = %d\n", st);
	fprintf(fp_fdbg, "end = %d\n", end);
#endif
	CONTROLLA(st > 0 && end <= LENGTHv_i(v));
	_CREAv_i(ris, end - st + 1);
	for (i = st, j = 1; i <= end; i++, j++)
		_ASSEGNAv_i(ris, j, _ACCEDIv_i(v, i));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _which_m_indxeq_i(VETTOREi *ris, const MATRICEi *m, int val, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_m_indxeq_i(VETTOREi *ris, const MATRICEi *m, int val)
#endif
{
	int r, c, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_m_indxeq_i\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "val = %d\n", val);
#endif
	_CREAv_i(ris, LENGTHm1_i(m) * LENGTHm2_i(m));
	j = 1;
	for (c = 1; c <= LENGTHm2_i(m); c++) {
		for (r = 1; r <= LENGTHm1_i(m); r++) {
			if (UGUALE(_ACCEDIm_i(m, r, c), val)) {
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
	VETTOREd * _dividi_vs_i(VETTOREd *ris, const VETTOREi *v, double div, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _dividi_vs_i(VETTOREd *ris, const VETTOREi *v, double div)
#endif
{
	int i;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: dividi_vs_i\n", linea);
#endif
#ifdef FDEBUG
	_StampaVett_i(v);
	fprintf(fp_fdbg, "div = %3.3f\n", div);
#endif
	_CREAv_d(ris, LENGTHv_i(v));
	if (ISNA(div)) {
#ifdef FDEBUG
	_CREAstr(tmp, "");
	g_string_printf(tmp, "ATTENZIONE (dividi_vs_i, linea %s # %d): divisione per NA!", nomefile, linea);
	warning(tmp->str);
	fprintf(fp_fdbg, tmp->str);
	_CANCELLAstr(tmp);
#endif
		for (i = 1; i <= LENGTHv_i(v); i++)
			_ASSEGNAv_d(ris, i, NA_REAL);
	}
else if (Uguale(div, 0.0)) {
#ifdef FDEBUG
	_CREAstr(tmp, "");
	g_string_printf(tmp, "ATTENZIONE (dividi_vs_i, linea %s # %d): divisione per zero!\n", nomefile, linea);
	warning(tmp->str);
	fprintf(fp_fdbg, tmp->str);
	_CANCELLAstr(tmp);
#endif
		for (i = 1; i <= LENGTHv_i(v); i++) {
			if (_ACCEDIv_i(v, i) > 0)
				_ASSEGNAv_d(ris, i, R_PosInf);
			else
				_ASSEGNAv_d(ris, i, R_NegInf);
		}
	}
	else {
		for (i = 1; i <= LENGTHv_i(v); i++)
			_ASSEGNAv_d(ris, i, _ACCEDIv_i(v, i) / div);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	void  _assegna1_v_indxv_i(const VETTOREi *v, const VETTOREi *indx, int val, const char *nomefile, int linea)
#else
	void  _assegna1_v_indxv_i(const VETTOREi *v, const VETTOREi *indx, int val)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_v_indxv_i\n", linea);
#endif
	CONTROLLA(v != NULL && indx != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
	_StampaVett_i(indx);
	fprintf(fp_fdbg, "val = %d\n", val);
#endif
	for (i = 1; i <= LENGTHv_i(indx); i++) {
		if (_ACCEDIv_i(indx, i) <= 0 || _ACCEDIv_i(indx, i) > LENGTHv_i(v)) {
#ifdef FDEBUG
			Rprintf("il vettore degli indici contiene alla posizione %d un elemento fuori dai limiti del vettore di riferimento (linea %s # %d)!\n", i, nomefile, linea);
#else
			Rprintf("il vettore degli indici contiene alla posizione %d un elemento fuori dai limiti del vettore di riferimento!\n", i);
#endif
			error("");
		}
		_ASSEGNAv_i(v, _ACCEDIv_i(indx, i), val);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(v);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	VETTOREi * _which_v_andglt_i(VETTOREi *ris, const VETTOREi *v, int val1, int val2, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_v_andglt_i(VETTOREi *ris, const VETTOREi *v, int val1, int val2)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_v_andglt_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
	fprintf(fp_fdbg, "val1 = %d\n", val1);
	fprintf(fp_fdbg, "val2 = %d\n", val2);
#endif
	_CREAv_i(ris, LENGTHv_i(v));
	for (i = 1, j = 1; i <= LENGTHv_i(v); i++) {
		if (_ACCEDIv_i(v, i) < val1 && _ACCEDIv_i(v, i) > val2) {
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
	VETTOREd * _distanza_2dvv_i(VETTOREd *ris, const VETTOREi *x1, const VETTOREi *y1, const VETTOREi *x2, const VETTOREi *y2, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _distanza_2dvv_i(VETTOREd *ris, const VETTOREi *x1, const VETTOREi *y1, const VETTOREi *x2, const VETTOREi *y2)
#endif
{
	int i, l;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: distanza_2dvv_i\n", linea);
#endif
	CONTROLLA(x1 != NULL && y1 != NULL && x2 != NULL && y2 != NULL);
#ifdef FDEBUG
	_StampaVett_i(x1);
	_StampaVett_i(y1);
	_StampaVett_i(x2);
	_StampaVett_i(y2);
#endif
	CONTROLLA(LENGTHv_i(x1) == LENGTHv_i(y1) && LENGTHv_i(x2) == LENGTHv_i(y2));
	l = min_s_i(LENGTHv_i(x1), LENGTHv_i(x2));
	_CREAv_d(ris, l);
	for (i = 1; i <= LENGTHv_i(x1); i++) {
		_ASSEGNAv_d(ris, i, sqrt(pow(_ACCEDIv_i(x2, i) - _ACCEDIv_i(x1, i), 2) + pow(_ACCEDIv_i(y2, i) - _ACCEDIv_i(y1, i), 2)));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREd * _distanza_2dvs_i(VETTOREd *ris, const VETTOREi *x1, const VETTOREi *y1, int x2, int y2, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _distanza_2dvs_i(VETTOREd *ris, const VETTOREi *x1, const VETTOREi *y1, int x2, int y2)
#endif
{
	int i, l;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: distanza_2dvs_i\n", linea);
#endif
	CONTROLLA(x1 != NULL && y1 != NULL);
#ifdef FDEBUG
	_StampaVett_i(x1);
	_StampaVett_i(y1);
	fprintf(fp_fdbg, "x2 = %d\n", x2);
	fprintf(fp_fdbg, "y2 = %d\n", y2);
#endif
	CONTROLLA(LENGTHv_i(x1) == LENGTHv_i(y1));
	l = LENGTHv_i(x1);
	_CREAv_d(ris, l);
	for (i = 1; i <= LENGTHv_i(x1); i++) {
		_ASSEGNAv_d(ris, i, sqrt(pow(x2 - _ACCEDIv_i(x1, i), 2) + pow(y2 - _ACCEDIv_i(y1, i), 2)));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	void  _assegna1_ms_indxcol_i(const MATRICEi *m, const VETTOREi *indx, int col, int val, const char *nomefile, int linea)
#else
	void  _assegna1_ms_indxcol_i(const MATRICEi *m, const VETTOREi *indx, int col, int val)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_ms_indxcol_i\n", linea);
#endif
	CONTROLLA(m != NULL && indx != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	_StampaVett_i(indx);
	fprintf(fp_fdbg, "colonna = %d\n", col);
	fprintf(fp_fdbg, "val = %d\n", val);
#endif
	CONTROLLA(col > 0 && indx != NULL);
	for (i = 1; i <= LENGTHv_i(indx); i++) {
		_ASSEGNAm_i(m, _ACCEDIv_i(indx, i), col, val);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return;
}

#ifdef MDEBUG
	void  _assegna1_mv_indx_i(const MATRICEi *m, const VETTOREi *indx, const VETTOREi *v, const char *nomefile, int linea)
#else
	void  _assegna1_mv_indx_i(const MATRICEi *m, const VETTOREi *indx, const VETTOREi *v)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_mv_indx_i\n", linea);
#endif
	CONTROLLA(m != NULL && indx != NULL && v != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	_StampaVett_i(indx);
	_StampaVett_i(v);
#endif
	CONTROLLA(LENGTHv_i(indx) == LENGTHv_i(v));
	for (i = 1; i <= LENGTHv_i(indx); i++)
		_ASSEGNAmv_i(m, _ACCEDIv_i(indx, i), _ACCEDIv_i(v, i));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	VETTOREi * _arrotonda_v_i(VETTOREi *ris, const VETTOREi *v, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _arrotonda_v_i(VETTOREi *ris, const VETTOREi *v)
#endif
{
#ifdef VERSIONE_d
	int i;
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: arrotonda_v_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
#endif
#ifdef VERSIONE_d
	_CREAv_i(ris, LENGTHv_i(v));
	for (i = 1; i <= LENGTHv_i(v); i++)
		_ASSEGNAv_i(ris, i, (int) _ACCEDIv_i(v, i));
#ifdef FDEBUG
		_CREAstr(tmp, "");
		g_string_printf(tmp, "ATTENZIONE (arrotonda_v_i, linea %s # %d): il vettore e` stato arrotondato da double a intero: possibile perdita di precisione nei calcoli!", nomefile, linea);
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
	VETTOREi * _exp_i(VETTOREi *ris, const VETTOREi *v, int val, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _exp_i(VETTOREi *ris, const VETTOREi *v, int val)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: exp_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
	fprintf(fp_fdbg, "val = %d\n", val);
#endif
	_CREAv_i(ris, LENGTHv_i(v));
	for (i = 1; i <= LENGTHv_i(v); i++)
		_ASSEGNAv_i(ris, i, pow(_ACCEDIv_i(v, i), val));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	VETTOREi * _moltiplica_vs_i(VETTOREi *ris, const VETTOREi *v, int val, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _moltiplica_vs_i(VETTOREi *ris, const VETTOREi *v, int val)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: moltiplica_vs_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
	fprintf(fp_fdbg, "val = %d\n", val);
#endif
	_CREAv_i(ris, LENGTHv_i(v));
	for (i = 1; i <= LENGTHv_i(v); i++)
		_ASSEGNAv_i(ris, i, _ACCEDIv_i(v, i) * val);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	void  _assegna1_v_segms_i(const VETTOREi *v, int st, int end, int val, const char *nomefile, int linea)
#else
	void  _assegna1_v_segms_i(const VETTOREi *v, int st, int end, int val)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna_v_segms_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
	fprintf(fp_fdbg, "start = %d\n", st);
	fprintf(fp_fdbg, "end = %d\n", end);
	fprintf(fp_fdbg, "val = %d\n", val);
#endif
	CONTROLLA(st > 0 && end <= LENGTHv_i(v));
	for (i = st; i <= end; i++)
		_ASSEGNAv_i(v, i, val);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(v);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	void  _assegna1_v_segmv_i(const VETTOREi *v1, int st, int end, const VETTOREi *v2, const char *nomefile, int linea)
#else
	void  _assegna1_v_segmv_i(const VETTOREi *v1, int st, int end, const VETTOREi *v2)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna_v_segmv_i\n", linea);
#endif
	CONTROLLA(v1 != NULL && v2 != NULL);
#ifdef FDEBUG
	_StampaVett_i(v1);
	fprintf(fp_fdbg, "start = %d\n", st);
	fprintf(fp_fdbg, "end = %d\n", end);
	_StampaVett_i(v2);
#endif
	CONTROLLA(st > 0 && end <= LENGTHv_i(v1) && LENGTHv_i(v2) >= end - st + 1);
	for (i = st, j = 1; i <= end; i++, j++)
		_ASSEGNAv_i(v1, i, _ACCEDIv_i(v2, j));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(v1);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	void  _assegna1_v_indxeq_i(const VETTOREi *v, int val1, int val2, const char *nomefile, int linea)
#else
	void  _assegna1_v_indxeq_i(const VETTOREi *v, int val1, int val2)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_v_indxeq_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
	fprintf(fp_fdbg, "val1 = %d\n", val1);
	fprintf(fp_fdbg, "val2 = %d\n", val2);
#endif
	for (i = 1; i <= LENGTHv_i(v); i++) {
		if (UGUALE(_ACCEDIv_i(v, i), val1))
			_ASSEGNAv_i(v, i, val2);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(v);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return;
}

#ifdef MDEBUG
	void  _assegna1_v_indxNA_i(const VETTOREi *v, int val, bool complemento, const char *nomefile, int linea)
#else
	void  _assegna1_v_indxNA_i(const VETTOREi *v, int val, bool complemento)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_v_indxNA_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
	fprintf(fp_fdbg, "val = %d\n", val);
	fprintf(fp_fdbg, "complemento = %d\n", complemento);
#endif
	for (i = 1; i <= LENGTHv_i(v); i++) {
		if ((_ACCEDIv_i(v, i) == NA_INTEGER) && !complemento)
			_ASSEGNAv_i(v, i, val);
		else if (!(_ACCEDIv_i(v, i) == NA_INTEGER) && complemento)
			_ASSEGNAv_i(v, i, val);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(v);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return;
}

#ifdef MDEBUG
	void  _incr1_v_indx_i(const VETTOREi *v, const VETTOREi *indx, int val, const char *nomefile, int linea)
#else
	void  _incr1_v_indx_i(const VETTOREi *v, const VETTOREi *indx, int val)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: incr1_v_indx_i\n", linea);
#endif
	CONTROLLA(v != NULL && indx != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
	_StampaVett_i(indx);
	fprintf(fp_fdbg, "val = %d\n", val);
#endif
	for (i = 1; i <= LENGTHv_i(indx); i++)
		_ASSEGNAv_i(v, _ACCEDIv_i(indx, i), _ACCEDIv_i(v, _ACCEDIv_i(indx, i)) + val);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(v);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	VETTOREd * _promuovi_i(VETTOREd *ris, const VETTOREi *v, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _promuovi_i(VETTOREd *ris, const VETTOREi *v)
#endif
{
#ifdef VERSIONE_i
	int i;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: promuovi_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
#endif
#ifdef VERSIONE_i
	_CREAv_d(ris, LENGTHv_i(v));
	for (i = 1; i <= LENGTHv_i(v); i++)
		_ASSEGNAv_d(ris, i, (double) _ACCEDIv_i(v, i));
#ifdef FDEBUG
		_CREAstr(tmp, "");
		g_string_printf(tmp, "ATTENZIONE (promuovi_i, linea %s # %d): il vettore e` stato promosso da intero a double: i calcoli saranno rallentati!", nomefile, linea);
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
	VETTOREi * _f_aux1_i(VETTOREi *ris, const VETTOREi *a, const VETTOREi *b, const VETTOREi *c, int da, int a1, int sgn1, int sgn2, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _f_aux1_i(VETTOREi *ris, const VETTOREi *a, const VETTOREi *b, const VETTOREi *c, int da, int a1, int sgn1, int sgn2)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: f_aux1_i\n", linea);
#endif
	CONTROLLA(a != NULL && b != NULL && c != NULL);
#ifdef FDEBUG
	_StampaVett_i(a);
	_StampaVett_i(b);
	_StampaVett_i(c);
#endif
	CONTROLLA(da > 0 && da < LENGTHv_i(a));
	CONTROLLA(a1 >= da && a1 < LENGTHv_i(a));
	CONTROLLA(LENGTHv_i(a) == LENGTHv_i(b) && LENGTHv_i(b) == LENGTHv_i(c));
#ifdef FDEBUG
	fprintf(fp_fdbg, "da = %d\n", da);
	fprintf(fp_fdbg, "a = %d\n", a1);
	fprintf(fp_fdbg, "sgn1 = %d\n", sgn1);
	fprintf(fp_fdbg, "sgn2 = %d\n", sgn2);
#endif
	_CREAv_i(ris, a1 - da + 1);
	for (i = da, j = 1; i <= a1; i++) {
		_ASSEGNAv_i(ris, j, _ACCEDIv_i(a, i) + sgn1 * _ACCEDIv_i(b, i) + sgn2 * _ACCEDIv_i(c, i));
		j++;
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	int  _which_min_v_i(VETTOREi *v, const char *nomefile, int linea)
#else
	int  _which_min_v_i(VETTOREi *v)
#endif
{
	int i, ultimo = 1;
	int min;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_min_v_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
#endif
	min = _ACCEDIv_i(v, 1);
	for (i = 2; i <= LENGTHv_i(v); i++) {
		if (_ACCEDIv_i(v, i) < min) {
			min = _ACCEDIv_i(v, i);
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
	void  _assegna1_ms_colindx_i(const MATRICEi *m, const VETTOREi *indx, int colonna, int val, const char *nomefile, int linea)
#else
	void  _assegna1_ms_colindx_i(const MATRICEi *m, const VETTOREi *indx, int colonna, int val)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_ms_colindx_i\n", linea);
#endif
	CONTROLLA(m != NULL && indx != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	_StampaVett_i(indx);
	fprintf(fp_fdbg, "colonna = %d\n", colonna);
	fprintf(fp_fdbg, "val = %d\n", val);
#endif
	CONTROLLA(colonna > 0 && indx != NULL);
	for (i = 1; i <= LENGTHv_i(indx); i++) {
		_ASSEGNAm_i(m, _ACCEDIv_i(indx, i), colonna, val);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return;
}

#ifdef MDEBUG
	VETTOREi * _f_aux2_i(VETTOREi *ris, const VETTOREi *a, const VETTOREi *b, const VETTOREi *c, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _f_aux2_i(VETTOREi *ris, const VETTOREi *a, const VETTOREi *b, const VETTOREi *c)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: f_aux2_i\n", linea);
#endif
	CONTROLLA(a != NULL && b != NULL && c != NULL);
#ifdef FDEBUG
	_StampaVett_i(a);
	_StampaVett_i(b);
	_StampaVett_i(c);
#endif
	CONTROLLA(LENGTHv_i(a) == LENGTHv_i(b) && LENGTHv_i(b) == LENGTHv_i(c));
	_CREAv_i(ris, LENGTHv_i(a));
	for (i = 1, j = 1; i <= LENGTHv_i(a); i++) {
		if ((DIVERSO(Segno(_ACCEDIv_i(a, i) - _ACCEDIv_i(c, i)) - Segno(_ACCEDIv_i(b, i) - _ACCEDIv_i(c, i)), 0)) && (DIVERSO(Segno(_ACCEDIv_i(b, i) - _ACCEDIv_i(c, i)) - Segno(_ACCEDIv_i(a, i) - _ACCEDIv_i(c, i)), 0))) {
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
	VETTOREi * _f_aux3_i(VETTOREi *ris, const VETTOREi *a, const VETTOREi *b, const VETTOREi *c, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _f_aux3_i(VETTOREi *ris, const VETTOREi *a, const VETTOREi *b, const VETTOREi *c)
#endif
{
	int i, j;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: f_aux3_i\n", linea);
#endif
	CONTROLLA(a != NULL && b != NULL && c != NULL);
#ifdef FDEBUG
	_StampaVett_i(a);
	_StampaVett_i(b);
	_StampaVett_i(c);
#endif
	CONTROLLA(LENGTHv_i(a) == LENGTHv_i(b) && LENGTHv_i(b) == LENGTHv_i(c));
	_CREAv_i(ris, LENGTHv_i(a));
	for (i = 1, j = 1; i <= LENGTHv_i(a); i++) {
		if ((_ACCEDIv_i(a, i) - _ACCEDIv_i(b, i)) > 0 && _ACCEDIv_i(a, i) < _ACCEDIv_i(c, i)) {
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
	VETTOREi * _unione_i(VETTOREi *ris, const VETTOREi *v1, const VETTOREi *v2, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _unione_i(VETTOREi *ris, const VETTOREi *v1, const VETTOREi *v2)
#endif
{
	int i, j;
	VETTOREi *tmp_ris = NULL;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: unione_i\n", linea);
#endif
	CONTROLLA(v1 != NULL && v2 != NULL);
#ifdef FDEBUG
	_StampaVett_i(v1);
	_StampaVett_i(v2);
#endif
		_CREAv_i(tmp_ris, LENGTHv_i(v1) + LENGTHv_i(v2));
		for (i = 1; i <= LENGTHv_i(v1); i++)
			_ASSEGNAv_i(tmp_ris, i, _ACCEDIv_i(v1, i));
		for (j = 1; j <= LENGTHv_i(v2); j++) {
			_ASSEGNAv_i(tmp_ris, i, _ACCEDIv_i(v2, j));
			i++;
		}
		tmp_ris->dim = i - 1;
		ris = elimina_doppi_i(ris, tmp_ris);
		_CANCELLAv_i(tmp_ris);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	int  _which_min_i(VETTOREi *v, const char *nomefile, int linea)
#else
	int  _which_min_i(VETTOREi *v)
#endif
{
	int i, ris = 1;
	int min;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_min_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
#endif
	min = _ACCEDIv_i(v, 1);
	for (i = 2; i <= LENGTHv_i(v); i++) {
		if (_ACCEDIv_i(v, i) < min) {
			min = _ACCEDIv_i(v, i);
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
	VETTOREi * _which_m_indxrowindxeq_i(VETTOREi *ris, const MATRICEi *m, const VETTOREi *indx, int val, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _which_m_indxrowindxeq_i(VETTOREi *ris, const MATRICEi *m, const VETTOREi *indx, int val)
#endif
{
	int i, j, k, r;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: which_m_indxrowindxeq_i\n", linea);
#endif
	CONTROLLA(m != NULL);
	CONTROLLA(indx != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	_StampaVett_i(indx);
	fprintf(fp_fdbg, "val = %d\n", val);
#endif
	CONTROLLA(LENGTHv_i(indx) <= LENGTHm2_i(m));
	_CREAv_i(ris, LENGTHv_i(indx) * LENGTHm2_i(m));
	for (i = 1, j = 1, k = 1; i <= LENGTHm2_i(m); i++) {
		for (r = 1; r <= LENGTHv_i(indx); r++) {
			if (UGUALE(_ACCEDIm_i(m, _ACCEDIv_i(indx, r), i), val)) {
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
	VETTOREi * _elimina_doppi_i(VETTOREi *ris, const VETTOREi *v, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _elimina_doppi_i(VETTOREi *ris, const VETTOREi *v)
#endif
{
#ifdef VERSIONE_i
#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: elimina_doppi_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
#endif
	GHashTable *table = g_hash_table_new(g_int_hash, g_int_equal);
	int i, j;

	_CREAv_i(ris, LENGTHv_i(v));
	for (i = 1, j = 1; i <= LENGTHv_i(v); i++) {
		if (!g_hash_table_lookup(table, &(v->dati[i - 1]))) {
			ASSEGNAv_i(ris, j, ACCEDIv_i(v, i));
			g_hash_table_insert(table, &(v->dati[i - 1]), &(v->dati[i - 1]));
			j++;
		}
	}
	ris->dim = j - 1;
	g_hash_table_destroy(table);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
#else
	error("Funzione 'elimina_doppi_d' non implementata!\n");
	return NULL;
#endif
}

#ifdef MDEBUG
	double  _f_aux4_i(const VETTOREi *times, double res, const char *nomefile, int linea)
#else
	double  _f_aux4_i(const VETTOREi *times, double res)
#endif
{
	int i, tmp_i;
	double max = 0.0, tmp_d;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: f_aux4_i\n", linea);
#endif
	CONTROLLA(times != NULL);
#ifdef FDEBUG
	_StampaVett_i(times);
	fprintf(fp_fdbg, "res = %3.3f\n\n", res);
#endif
	if (Uguale(res, 0.0)) {
#ifdef FDEBUG
		_CREAstr(tmp, "");
		g_string_printf(tmp, "ATTENZIONE (f_aux4_i, linea %s # %d): divisione per zero!\n", nomefile, linea);
		warning(tmp->str);
		fprintf(fp_fdbg, tmp->str);
		_CANCELLAstr(tmp);
#endif
		return R_PosInf;
	}
	for (i = 1; i <= LENGTHv_i(times); i++) {
		tmp_d = _ACCEDIv_i(times, i) / res;
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
	void  _assegna1_mv_colonna_i(const MATRICEi *m, int colonna, const VETTOREi *v, const char *nomefile, int linea)
#else
	void  _assegna1_mv_colonna_i(const MATRICEi *m, int colonna, const VETTOREi *v)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_mv_colonna_i\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "colonna = %d\n", colonna);
	_StampaVett_i(v);
#endif
	CONTROLLA(colonna > 0 && colonna <= m->nc && v != NULL && m->nr == v->dim);
	for (i = 1; i <= LENGTHm1_i(m); i++)
		_ASSEGNAm_i(m, i, colonna, _ACCEDIv_i(v, i));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return;
}

#ifdef MDEBUG
	VETTOREd * _f_aux5_i(VETTOREd *ris, const VETTOREi *alpha, const VETTOREi *targ, const VETTOREi *theta, const VETTOREi *xmin, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _f_aux5_i(VETTOREd *ris, const VETTOREi *alpha, const VETTOREi *targ, const VETTOREi *theta, const VETTOREi *xmin)
#endif
{
	int i;
	double tmp1;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: f_aux5_i\n", linea);
#endif
	CONTROLLA(alpha != NULL && targ != NULL && theta != NULL && xmin != NULL);
#ifdef FDEBUG
	_StampaVett_i(alpha);
	_StampaVett_i(targ);
	_StampaVett_i(theta);
	_StampaVett_i(xmin);
#endif
	_CREAv_d(ris, LENGTHv_i(alpha));
	for (i = 1; i <= LENGTHv_i(alpha); i++) {
		tmp1 = 1 + pow(2.71828182845904523536, -_ACCEDIv_i(alpha, i) * (_ACCEDIv_i(targ, i) - _ACCEDIv_i(theta, i)));
		if (Uguale(tmp1, 0.0)) {
#ifdef FDEBUG
			_CREAstr(tmp, "");
			g_string_printf(tmp, "ATTENZIONE (f_aux5_i, linea %s # %d): l'elemento %d ha provocato una divisione per zero e gli e` stato assegnato un valore al di fuori del dominio!\n", nomefile, linea, i);
			warning(tmp->str);
			fprintf(fp_fdbg, tmp->str);
			_CANCELLAstr(tmp);
#endif
			tmp1 = R_PosInf;
		}
		else {
			_ASSEGNAv_d(ris, i, (1 / tmp1) * (1 - _ACCEDIv_i(xmin, i)) + _ACCEDIv_i(xmin, i));
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
	VETTOREd * _f_aux6_i(VETTOREd *ris, double res, const VETTOREi *k, const VETTOREi *targetT, const VETTOREi *n, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _f_aux6_i(VETTOREd *ris, double res, const VETTOREi *k, const VETTOREi *targetT, const VETTOREi *n)
#endif
{
	int i;
	double tmp1;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: f_aux6_i\n", linea);
#endif
	CONTROLLA(k != NULL && targetT != NULL && n != NULL);
#ifdef FDEBUG
	fprintf(fp_fdbg, "res = %3.3f\n", res);
	_StampaVett_i(k);
	_StampaVett_i(targetT);
	_StampaVett_i(n);
#endif
	_CREAv_d(ris, LENGTHv_i(k));
	for (i = 1; i <= LENGTHv_i(k); i++) {
		tmp1 = res * _ACCEDIv_i(k, i) * (_ACCEDIv_i(targetT, i) - _ACCEDIv_i(n, i));
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
	void  _somma1_vv_i(VETTOREi *v1, const VETTOREi *v2, const char *nomefile, int linea)
#else
	void  _somma1_vv_i(VETTOREi *v1, const VETTOREi *v2)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: somma1_vv_i\n", linea);
#endif
	CONTROLLA(v1 != NULL && v2 != NULL);
#ifdef FDEBUG
	_StampaVett_i(v1);
	_StampaVett_i(v2);
#endif
	for (i = 1; i <= LENGTHv_i(v1); i++)
		_ASSEGNAv_i(v1, i, _ACCEDIv_i(v1, i) + _ACCEDIv_i(v2, i));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(v1);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return;
}

#ifdef MDEBUG
	MATRICEi * _aggiungi_mv_colonna_i(MATRICEi *m, int colonna, const VETTOREi *v, const char *nome, const char *nomefile, int linea)
#else
	MATRICEi * _aggiungi_mv_colonna_i(MATRICEi *m, int colonna, const VETTOREi *v)
#endif
{
	int i, r, c;
	MATRICEi *tmp_m = NULL;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: aggiungi_mv_colonna_i\n", linea);
#endif
	CONTROLLA(m != NULL && v != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "colonna = %d\n", colonna);
	_StampaVett_i(v);
#endif
	CONTROLLA(m->nr == v->dim);
	if (colonna > m->nc + 1) {
		Rprintf("la colonna da aggiungere non e` minore o uguale al numero di colonne attuali + 1: R riempirebbe quelle in mezzo con NA, e` questa l'intenzione?\n");
		error("");
	}
	if (colonna > m->alloc_c) {
#ifdef FDEBUG
		_CREAstr(tmp, "");
		g_string_printf(tmp, "ATTENZIONE (aggiungi_mv_colonna_i, linea %s # %d): ingrandite le colonne della matrice da %d a %d\n", nomefile, linea, LENGTHm2_i(m), colonna);
		warning(tmp->str);
		fprintf(fp_fdbg, tmp->str);
		_CANCELLAstr(tmp);
#endif
		_CREAm_i(tmp_m, m->nr, 2 * colonna);
		tmp_m->nc = colonna;
		for (r = 1; r <= LENGTHm1_i(m); r++) {
			for (c = 1; c <= LENGTHm2_i(m); c++)
				_ASSEGNAm_i(tmp_m, r, c, _ACCEDIm_i(m, r, c));
		}
		_CANCELLAm_i(m);
		m = tmp_m;
	}
	for (i = 1; i <= LENGTHm1_i(m); i++)
		_ASSEGNAm_i(m, i, colonna, _ACCEDIv_i(v, i));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return m;
}

#ifdef MDEBUG
	int  _f_aux7_i(const VETTOREi *vettore, int scalare, const char *nomefile, int linea)
#else
	int  _f_aux7_i(const VETTOREi *vettore, int scalare)
#endif
{
	int i, indx = 0;
	double min = R_PosInf, tmp;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: f_aux7_i\n", linea);
#endif
	CONTROLLA(vettore != NULL);
#ifdef FDEBUG
	_StampaVett_i(vettore);
	fprintf(fp_fdbg, "scalare = %d\n\n", scalare);
#endif
	for (i = 1; i <= LENGTHv_i(vettore); i++) {
		tmp = _ACCEDIv_i(vettore, i) - scalare;
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
	VETTOREi * _f_aux8_i(VETTOREi *ris, const VETTOREi *a, const VETTOREi *b, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi * _f_aux8_i(VETTOREi *ris, const VETTOREi *a, const VETTOREi *b)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: f_aux8_i\n", linea);
#endif
	CONTROLLA(a != NULL && b != NULL);
#ifdef FDEBUG
	_StampaVett_i(a);
	_StampaVett_i(b);
#endif
	_CREAv_i(ris, LENGTHv_i(a));
	for (i = 1; i <= LENGTHv_i(a); i++) {
		_ASSEGNAv_i(ris, i, _ACCEDIv_i(a, i) * (1 - _ACCEDIv_i(b, i)) + _ACCEDIv_i(b, i));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	MATRICEi * _seleziona_colonne_i(MATRICEi *ris, const MATRICEi *m, const VETTOREi *indici, const char *nome, const char *nomefile, int linea)
#else
	MATRICEi * _seleziona_colonne_i(MATRICEi *ris, const MATRICEi *m, const VETTOREi *indici)
#endif
{
	int r, c;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: seleziona_colonne_i\n", linea);
#endif
	CONTROLLA(m != NULL && indici != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	_StampaVett_i(indici);
#endif
	_CREAm_i(ris, LENGTHm1_i(m), LENGTHv_i(indici));
	for (r = 1; r <= LENGTHm1_i(m); r++) {
		for (c = 1; c <= LENGTHv_i(indici); c++)
		_ASSEGNAm_i(ris, r, c, _ACCEDIm_i(m, r, _ACCEDIv_i(indici, c)));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	void  _cambiadim1_i(MATRICEi *m, const int nr, const int nc, const char *nomefile, int linea)
#else
	void  _cambiadim1_i(MATRICEi *m, const int nr, const int nc)
#endif
{
	int dim_r, dim_c;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: cambiadim1_i\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
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
	MATRICEi * _moltiplica_mm_i(MATRICEi *ris, const MATRICEi *m1, const MATRICEi *m2, const char *nome, const char *nomefile, int linea)
#else
	MATRICEi * _moltiplica_mm_i(MATRICEi *ris, const MATRICEi *m1, const MATRICEi *m2)
#endif
{
	int r, c;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: moltiplica_mm_i\n", linea);
#endif
	CONTROLLA(m1 != NULL && m2 != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m1);
	_StampaMatr_i(m2);
#endif
	CONTROLLA(m1->nr == m2->nr && m1->nc == m2->nc);
	_CREAm_i(ris, m1->nr, m1->nc);
	for (r = 1; r <= LENGTHm1_i(m1); r++) {
		for (c = 1; c <= LENGTHm2_i(m1); c++)
			_ASSEGNAm_i(ris, r, c, _ACCEDIm_i(m1, r, c) * _ACCEDIm_i(m2, r, c));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	void  _moltiplica1_mv_i(const MATRICEi *m, const VETTOREi *v, const char *nomefile, int linea)
#else
	void  _moltiplica1_mv_i(const MATRICEi *m, const VETTOREi *v)
#endif
{
	int r, c, i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: moltiplica_mv_i\n", linea);
#endif
	CONTROLLA(m != NULL && v != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	_StampaVett_i(v);
#endif
	if (((m->nr * m->nc) % v->dim) != 0) {
		error("Il numero di elementi della matrice non e` un multiplo del numero di elementi del vettore!\n");
	}
	for (r = 1, i = 1; r <= LENGTHm1_i(m); r++) {
		for (c = 1; c <= LENGTHm2_i(m); c++) {
			_ASSEGNAm_i(m, r, c, _ACCEDIm_i(m, r, c) * _ACCEDIv_i(v, i));
			i++;
			if (i > LENGTHv_i(v))
				i = 1;
		}
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	VETTOREd * _dividi_vv_i(VETTOREd *ris, const VETTOREi *v1, const VETTOREi *v2, const char *nome, const char *nomefile, int linea)
#else
	VETTOREd * _dividi_vv_i(VETTOREd *ris, const VETTOREi *v1, const VETTOREi *v2)
#endif
{
	int i;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: dividi_vv_i\n", linea);
#endif
	CONTROLLA(v1 != NULL && v2 != NULL);
#ifdef FDEBUG
	_StampaVett_i(v1);
	_StampaVett_i(v2);
#endif
	_CREAv_d(ris, LENGTHv_i(v1));
	for (i = 1; i <= LENGTHv_i(v1); i++) {
		if ((_ACCEDIv_i(v2, i) == NA_INTEGER)) {
#ifdef FDEBUG
			_CREAstr(tmp, "");
			g_string_printf(tmp, "ATTENZIONE (dividi_vv_i, linea %s # %d): l'elemento %d del vettore divisore ha causato una divisione per NA!", nomefile, linea, i);
			warning(tmp->str);
			fprintf(fp_fdbg, tmp->str);
			_CANCELLAstr(tmp);
#endif
			_ASSEGNAv_d(ris, i, NA_REAL);
		}
		else if (UGUALE(_ACCEDIv_i(v2, i), 0)) {
#ifdef FDEBUG
			_CREAstr(tmp, "");
			g_string_printf(tmp, "ATTENZIONE (dividi_vv_i, linea %s # %d): l'elemento %d del vettore divisore ha causato una divisione per zero!\n", nomefile, linea, i);
			warning(tmp->str);
			fprintf(fp_fdbg, tmp->str);
			_CANCELLAstr(tmp);
#endif
			if (_ACCEDIv_i(v1, i) > 0)
				_ASSEGNAv_d(ris, i, R_PosInf);
			else
				_ASSEGNAv_d(ris, i, R_NegInf);
		}
		else
			_ASSEGNAv_d(ris, i, _ACCEDIv_i(v1, i) / _ACCEDIv_i(v2, i));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	MATRICEi * _sign_m_i(MATRICEi *ris, const MATRICEi *m, const char *nome, const char *nomefile, int linea)
#else
	MATRICEi * _sign_m_i(MATRICEi *ris, const MATRICEi *m)
#endif
{
	int r, c;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: sign_m_i\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
#endif
	_CREAm_i(ris, LENGTHm1_i(m), LENGTHm2_i(m));

	for (r = 1; r <= LENGTHm1_i(m); r++) {
		for (c = 1; c <= LENGTHm2_i(m); c++) {
			if (_ACCEDIm_i(m, r, c) > 0)
				_ASSEGNAm_i(ris, r, c, 1);
			else if (_ACCEDIm_i(m, r, c) < 0)
				_ASSEGNAm_i(ris, r, c, -1);
			else
				_ASSEGNAm_i(ris, r, c, 0);
		}
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}

#ifdef MDEBUG
	void  _assegna1_ms_colonna_i(const MATRICEi *m, int colonna, int val, const char *nomefile, int linea)
#else
	void  _assegna1_ms_colonna_i(const MATRICEi *m, int colonna, int val)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: assegna1_ms_colonna_i\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "colonna = %d\n", colonna);
	fprintf(fp_fdbg, "val = %d\n", val);
#endif
	CONTROLLA(colonna > 0 && colonna <= m->nc);
	for (i = 1; i <= LENGTHm1_i(m); i++)
		_ASSEGNAm_i(m, i, colonna, val);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_i(m);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return;
}

#ifdef MDEBUG
	void  _dividi1_vv_i(VETTOREd *v1, const VETTOREi *v2, const char *nomefile, int linea)
#else
	void  _dividi1_vv_i(VETTOREd *v1, const VETTOREi *v2)
#endif
{
	int i;
#ifdef FDEBUG
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: dividi1_vv_i\n", linea);
#endif
	CONTROLLA(v1 != NULL && v2 != NULL);
#ifdef FDEBUG
	_StampaVett_d(v1);
	_StampaVett_i(v2);
#endif
	for (i = 1; i <= LENGTHv_d(v1); i++) {
		if ((_ACCEDIv_i(v2, i) == NA_INTEGER)) {
#ifdef FDEBUG
			_CREAstr(tmp, "");
			g_string_printf(tmp, "ATTENZIONE (dividi_vv_i, linea %s # %d): l'elemento %d del vettore divisore ha causato una divisione per NA!", nomefile, linea, i);
			warning(tmp->str);
			fprintf(fp_fdbg, tmp->str);
			_CANCELLAstr(tmp);
#endif
			_ASSEGNAv_d(v1, i, NA_REAL);
		}
		else if (UGUALE(_ACCEDIv_i(v2, i), 0)) {
#ifdef FDEBUG
			_CREAstr(tmp, "");
			g_string_printf(tmp, "ATTENZIONE (dividi_vv_i, linea %s # %d): l'elemento %d del vettore divisore ha causato una divisione per zero!\n", nomefile, linea, i);
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
			_ASSEGNAv_d(v1, i, _ACCEDIv_d(v1, i) / _ACCEDIv_i(v2, i));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_d(v1);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return;
}

#ifdef MDEBUG
	void  _abs1_v_i(VETTOREi *v, const char *nomefile, int linea)
#else
	void  _abs1_v_i(VETTOREi *v)
#endif
{
	int i;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: abs1_v_i\n", linea);
#endif
	CONTROLLA(v != NULL);
#ifdef FDEBUG
	_StampaVett_i(v);
#endif
	for (i = 1; i <= LENGTHv_i(v); i++)
		_ASSEGNAv_i(v, i, abs(_ACCEDIv_i(v, i)));
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(v);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
}

#ifdef MDEBUG
	MATRICEi * _arrotonda_m_i(MATRICEi *ris, const MATRICEi *m, const char *nome, const char *nomefile, int linea)
#else
	MATRICEi * _arrotonda_m_i(MATRICEi *ris, const MATRICEi *m)
#endif
{
#ifdef VERSIONE_d
	int r, c;
	GString *tmp = NULL;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: arrotonda_m_i\n", linea);
#endif
	CONTROLLA(m != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
#endif
#ifdef VERSIONE_d
	_CREAm_i(ris, LENGTHm1_i(m), LENGTHm2_i(m));
	for (r = 1; r <= LENGTHm1_i(m); r++) {
		for (c = 1; c <= LENGTHm2_i(m); c++)
			_ASSEGNAm_i(ris, r, c, (int) _ACCEDIm_i(m, r, c));
	}
#ifdef FDEBUG
		_CREAstr(tmp, "");
		g_string_printf(tmp, "ATTENZIONE (arrotonda_m_i, linea %s # %d): la matrice e` stata arrotondata da double a intero: possibile perdita di precisione nei calcoli!", nomefile, linea);
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

#undef VERSIONE_i
