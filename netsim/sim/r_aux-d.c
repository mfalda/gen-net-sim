#include "r_aux.h"


int Segno(int x)
{
	if (x > 0)
		return 1;
	else if (x < 0)
		return -1;
	else
		return 0;
}

// interna
bool Uguale(double x, double y)
{
	if (x == R_PosInf) {
		if (y == R_PosInf)
			return true;
		else
			return false;
	}
	else if (x == R_NegInf) {
		if (y == R_NegInf)
			return true;
		else
			return false;
	}
	// Knuth section 4.2.2 pages 217-218
	return fabs(x - y) <= DBL_EPSILON * fabs(x);
}

// alternativa, � Audric Thevenet
/*
inline int GetExpoBase2(double d)
{
	int i = 0;
	((short *)(&i))[0] = (((short *)(&d))[3] & (short)32752); // _123456789ab____ & 0111111111110000
	return (i >> 4) - 1023;
}

bool Uguale(double x, double y)
{
	if (d1 == d2)
		return true;
	int e1 = GetExpoBase2(d1);
	int e2 = GetExpoBase2(d2);
	int e3 = GetExpoBase2(d1 - d2);
	if ((e3 - e2 < -48) && (e3 - e1 < -48))
		return true;
	return false;
}
*/

// interna
void _errore(int val, const char *msg, const char *nomefile, int linea)
{
	char err[64];

	if (!val) {
		snprintf(err, 64, "ASSERZIONE FALLITA (linea %s # %d): %s !\n", nomefile, linea, msg);
		error(err);
	}
}

void _StampaLista(const LISTA *l)
{
#ifdef FDEBUG
	int i;
	Allocazione *t;

	if (l == NULL) {
		fprintf(fp_fdbg, "lista nulla\n");
		return;
	}
	assert(l != NULL);
	t = (Allocazione *) l->mem->data;
	fprintf(fp_fdbg, "%s (%d : %d): {\n", t->nome->str, l->dim, l->mia_alloc);
	for (i = 0; i < l->dim; i++) {
		fprintf(fp_fdbg, "\t");
		switch (l->tipi[i]) {
			case INTERO:
			case VETTi:
				_StampaVett_i(l->dati[i].vi);
				break;
			case REALE:
			case VETTd:
				_StampaVett_d(l->dati[i].vd);
				break;
			case MATRi:
				_StampaMatr_i(l->dati[i].mi);
				break;
			case MATRd:
				_StampaMatr_d(l->dati[i].md);
				break;
			case STRINGA:
				fprintf(fp_fdbg, "'%s'\n", l->dati[i].str->str);
				break;
		}
	}
	fprintf(fp_fdbg, "}\n\n");
#endif
}

void _StampaRawLista(const LISTA *l)
{
#ifdef DET
	int i;
	Allocazione *t;

	if (l == NULL) {
		return;
	}
	assert(l != NULL);
	t = (Allocazione *) l->mem->data;
	for (i = 0; i < l->dim; i++) {
		fprintf(fp_det, " ");
		switch (l->tipi[i]) {
			case INTERO:
			case VETTi:
				_StampaRawVett_i(l->dati[i].vi);
				break;
			case REALE:
			case VETTd:
				_StampaRawVett_d(l->dati[i].vd);
				break;
			case MATRi:
				_StampaRawMatr_i(l->dati[i].mi);
				break;
			case MATRd:
				_StampaRawMatr_d(l->dati[i].md);
				break;
			case STRINGA:
				fprintf(fp_det, "'%s'\n", l->dati[i].str->str);
				break;
		}
	}
	fprintf(fp_det, "\n");
#endif
}

// interna
#ifdef MDEBUG
void _stampa_lista_lst(const char *pref, GList *lista)
{
	GList *li;
	Allocazione *t;


	fprintf(fp_mdbg_lst, "%s", pref);
	for(li = lista; li != NULL; li = li->next) {
		t = (Allocazione *) li->data;
		fprintf(fp_mdbg_lst, "%s::%s\t", t->file_da->str, t->nome->str);
	}
	fprintf(fp_mdbg_lst, "\n\n");
}
#endif

#ifdef MDEBUG
	int  _lengthList(const LISTA *l, const char *nome, const char *nomefile, int linea)
#else
	int  _lengthList(const LISTA *l, const char *nome)
#endif
{
#ifdef MDEBUG
	char err[128];
	if (l == NULL) {
		snprintf(err, 64, "la lista '%s' (linea %s # %d) non e` allocata: non posso ricavarne le dimensioni!\n", nome, nomefile, linea);
		error(err);
	}
#endif
	return l->dim;
}

void _Intestazione(const char *msg)
{
#ifdef FDEBUG
	fprintf(fp_fdbg, msg);
#endif
#ifdef MDEBUG
	fprintf(fp_mdbg_i, msg);
	fprintf(fp_mdbg_d, msg);
#endif
#ifdef DET
	fprintf(fp_det, msg);
#endif
}

// interna
void _InitDbg_lst(bool stdout1)
{
#ifdef MDEBUG
	tot_str_alloc = 0;
	tot_str_dealloc = 0;
	tot_str_intalloc = 0;
	tot_str_intdealloc = 0;
	fpm_lst = fopen("memoria_lst.csv", "a");
	rewind(fpm_lst);
	fprintf(fpm_lst, "\"ID\";\"nome\";\"da\";\"a\";\"dim.\";\"riallocazioni\";\"prima linea\";\"ultima linea\"\n");
	allocLst = NULL;
	if (stdout1)
		fp_mdbg_lst = stdout;
	else {
		fp_mdbg_lst = fopen("memoria_lst.txt", "a");
		rewind(fp_mdbg_lst);
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
#endif
#ifdef DET
	if (!fp_det) {
		if (stdout1)
			fp_det = stdout;
		else {
			fp_det = fopen("moduli.txt", "a");
			rewind(fp_det);
		}
	}
#endif
	return;
}

void _InitDbg(bool stdout1_v, bool stdout1_m, bool stdout1_lst)
{
	_InitDbg_i(stdout1_v);
	_InitDbg_d(stdout1_m);
	_InitDbg_lst(stdout1_lst);
}

void _Fflush_ls()
{
#ifdef MDEBUG
	fflush(fp_mdbg_lst);
#endif
#ifdef FDEBUG
	fflush(fp_fdbg);
#endif
#ifdef DET
	fflush(fp_det);
#endif
	return;
}

#ifdef MDEBUG
	LISTA *_inLISTA(SEXP s, int *nProtected, int len, const enum TIPO *t, GString **nomi, const char *nome, const char *nomefile, int linea)
#else
	LISTA *_inLISTA(SEXP s, int *nProtected, int len, const enum TIPO *t)
#endif
{
	int i, ll;
	LISTA *l = NULL;
	VETTOREi *tmpv_i = NULL;
	VETTOREd *tmpv_d = NULL;
	MATRICEi *tmpm_i = NULL;
	MATRICEd *tmpm_d = NULL;
	SEXP tmp;
#ifdef MDEBUG
	Allocazione *mem;
#endif

	if (isNull(s)) {
#ifdef FDEBUG
	fprintf(fp_fdbg, "La lista '%s' e` nulla\n", nome);
#endif
		return NULL;
	}
	ll = length(s);
#ifdef MDEBUG
	l = _CreaLista(l, t, ll, nome, nomefile, linea);
#else
	l = _CreaLista(l, t, len);
#endif
	for (i = 0; i < ll; i++) {
		tmp = VECTOR_ELT(s, i);
		if (t == NULL) {
#ifdef MDEBUG
			tmpv_i = _inVETTORE_i(tmp, nProtected, nomi[i]->str, nomefile, linea);
#else
			tmpv_i = _inVETTORE_i(tmp, nProtected);
#endif
			l->dati[i].vi = tmpv_i;
		}
		else {
			switch (t[i]) {
				case INTERO:
				case VETTi:
	#ifdef MDEBUG
					tmpv_i = _inVETTORE_i(tmp, nProtected, nomi[i]->str, nomefile, linea);
	#else
					tmpv_i = _inVETTORE_i(tmp, nProtected);
	#endif
					l->dati[i].vi = tmpv_i;
					break;
				case REALE:
				case VETTd:
	#ifdef MDEBUG
					tmpv_d = _inVETTORE_d(tmp, nProtected, nomi[i]->str, nomefile, linea);
	#else
					tmpv_d = _inVETTORE_d(tmp, nProtected);
	#endif
					l->dati[i].vd = tmpv_d;
					break;
				case MATRi:
	#ifdef MDEBUG
					tmpm_i = _inMATRICE_i(tmp, nProtected, nomi[i]->str, nomefile, linea);
	#else
					tmpm_i = _inMATRICE_i(tmp, nProtected);
	#endif
					l->dati[i].mi = tmpm_i;
					break;
				case MATRd:
	#ifdef MDEBUG
					tmpm_d = _inMATRICE_d(tmp, nProtected, nomi[i]->str, nomefile, linea);
	#else
					tmpm_d = _inMATRICE_d(tmp, nProtected);
	#endif
					l->dati[i].md = tmpm_d;
					break;
				case STRINGA:
					l->dati[i].str = inSTRINGA(tmp, nProtected, nomi[i]->str);
					break;
			}
		}
	}
	l->r = true;
#ifdef MDEBUG
	mem = l->mem->data;
	g_string_assign(mem->nome, nome);
	g_string_assign(mem->file_da, nomefile);
	mem->linea_da = linea;
#endif
#ifdef FDEBUG
	fprintf(fp_fdbg, "Ho trasformato la lista ");
	_StampaLista(l);
#endif
	return l;
}

#ifdef MDEBUG
	LISTA *_CreaLista(LISTA *ris, const enum TIPO *t, int len, const char *nome, const char *nomefile, int linea)
#else
	LISTA *_CreaLista(LISTA *ris, const enum TIPO *t, int len)
#endif
{
	int i;
	LISTA *ris1;
	bool nuova = false;
	char err[64];
#ifdef MDEBUG
	Allocazione *mem;
#endif

	// se non esiste la creo
	if (ris == NULL) {
		nuova = true;
		ris1 = mia_alloc(1, LISTA);
		if (ris1 == NULL) {
			snprintf(err, 64, "Not enough memory (CreaLISTA # %d, ris1)", __LINE__ - 2);
		error(err);
		}
		ris1->dim = len;
		ris1->mia_alloc = len; // inizialmente mia_alloc == len
		ris1->dati = NULL;
		ris1->tipi = NULL;
		ris1->mem = NULL;
		ris1->r = false;
	}
	// altrimenti punto all'esistente
	else
		ris1 = ris;
#ifdef MDEBUG
	// se e` da R e` un errore, perche� non la posso alterare!
	if (ris != NULL && ris->r)
		error("Una lista proveniente da R non puo` essere alterata! (CreaLISTA)");
	// se non l'ho trovata
	if (ris1->mem == NULL) {
		mem = mia_alloc(1, Allocazione);
		if (mem == NULL) {
			snprintf(err, 64, "Not enough memory (CreaLISTA # %d, ris1->mem)", __LINE__ - 2);
			error(err);
		}
		fprintf(fp_mdbg_lst, "%d_lst - allocazione della lista '%s'[%d] (%p) dalla linea %s # %d\n", g_list_length(allocLst) + 1, nome, ris1->dim, ris1, nomefile, linea);
		mem->indir = (size_t) ris1;
		_CREAstr(mem->nome, nome);
		_CREAstr(mem->file_da, nomefile);
		_CREAstr(mem->file_a, "");
		mem->linea_da = linea;
		mem->linea_a = 0;
		mem->max_dim1 = len;
		mem->rialloc = 0;
		mem->in_uso = 0;
		_CREAstr(mem->file_prima, "");
		mem->prima_linea = 0;
		_CREAstr(mem->file_ultima, "");
		mem->ultima_linea = 0;
		allocLst = g_list_append(allocLst, mem);
		mem->indx = g_list_length(allocLst);
		_stampa_lista_lst("allocLst+: ", allocLst);
		ris1->mem = g_list_last(allocLst);
	}
	else {
		mem = (Allocazione *) ris1->mem->data;
		if (mem->linea_da > 0 && mem->linea_a == 0) {
			fprintf(fp_mdbg_lst, "%d_lst - la lista '%s' (linea %s # %d) esiste gia`: e` stata allocata alla linea %s # %d\n\n", mem->indx, nome, nomefile, linea, mem->file_da->str, mem->linea_da);
		}
	}
	// se non e` abbastanza grande segnalo il fatto (la ridimensionero` effettivamente solo alla fine, nel codice comune)
	if (ris1->mia_alloc < len) {
		fprintf(fp_mdbg_lst, "*** la lista '%s' verra` riallocata passando da %d a %d (%d riallocazione/i, finora)\n\n", nome, ris1->mia_alloc, ris1->mia_alloc * 2 + len, mem->rialloc + 1);
		mem->rialloc++;
	}
#endif
	CONTROLLA(len >= 0);
	if (len > 0) {
		// se non e` abbastanza grande la rialloco raddoppiando mia_alloc (puo` capitare solo se la sto riutilizzando)
		if (ris1->mia_alloc < len) {
			if (ris1->mia_alloc > 0) {
				warning("La lista sara` riallocata e tutti gli elementi saranno cancellati!\n");
				CancellaLISTA(ris1, true);
				libera(ris1->dati);
				libera(ris1->tipi);
				ris1->dati = NULL;
				ris1->tipi = NULL;
			}
			ris1->mia_alloc = ris1->mia_alloc * 2 + len;
			ris1->dati = mia_alloc((ris1->mia_alloc), union Dati);
			ris1->tipi = mia_alloc((ris1->mia_alloc), enum TIPO);
		}
		// se non esiste la alloco
		else if (ris1->dati == NULL || ris1->tipi == NULL) {
			ris1->dati = mia_alloc(len, union Dati);
			ris1->tipi = mia_alloc(len, enum TIPO);
		}
		// in ogni caso a questo punto se e` ancora NULL e` un errore!
		if (ris1->dati == NULL || ris1->tipi == NULL) { // non puo` mai essere zero la dimensione richiesta, dato che dim > 0
			snprintf(err, 64, "Not enough memory (CreaLISTA # %d, ris1->dati)", __LINE__ - 2);
			error(err);
		}
	}
	ris1->dim = len;
	// alla fine assegno i tipi
	for (i = 0; i < len; i++) {
		if (nuova)
			ris1->dati[i].vi = NULL;
		if (t == NULL)
			ris1->tipi[i] = VETTi;
		else
			ris1->tipi[i] = t[i];
	}
	return ris1;
}

#ifdef MDEBUG
	LISTA *_CancellaLista(LISTA *l, const bool tutto, const char *nome, const char *nomefile, int linea)
#else
	LISTA *_CancellaLista(LISTA *l, const bool tutto)
#endif
{
	int i;
#ifdef MDEBUG
	char *nome_elem;
	Allocazione *mem;
	char err[128];
#endif

	if (l == NULL)
		return NULL;
#ifdef MDEBUG
	// se non e` stato trovato o da == 0 e` un errore!
	if (l->mem == NULL) {
		snprintf(err, 128, "\tla lista '%s' non e` allocata! (linea %s # %d)\n", nome, nomefile, linea);
		error(err);
	}
	mem = (Allocazione *) l->mem->data;
	// qui si` il nome puo` avere piu` senso
	g_string_assign(mem->file_a, nomefile);
	fprintf(fp_mdbg_lst, "%d_lst - disallocazione della lista '%s'[%d] (%p) dalla linea %s # %d\n", mem->indx, mem->nome->str, l->dim, l, nomefile, linea);
	mem->linea_a = linea;
	fprintf(fpm_lst, "\"%d\";\"'%s'\";\"%s # %d\";\"%s # %d\";\"%d\";\"%d\";\"%s # %d\";\"%s # %d\"\n",  mem->indx, mem->nome->str, mem->file_da->str, mem->linea_da, mem->file_a->str, mem->linea_a, mem->max_dim1, mem->rialloc, mem->file_prima->str, mem->prima_linea, mem->file_ultima->str, mem->ultima_linea);
#endif
	// il contenuto della lista l'ho comunque creato io
	if (tutto && LENGTHlst(l) > 0) {
		for (i = 0; i < l->dim; i++) {
			switch (l->tipi[i]) {
				case INTERO:
				case VETTi:
#ifdef MDEBUG
					if (l->dati[i].vi != NULL) {
						nome_elem = ((Allocazione *) l->dati[i].vi->mem->data)->nome->str;
						l->dati[i].vi = _cancellaVett_i(l->dati[i].vi, nome_elem, nomefile, linea);
					}
#else
					CANCELLAv_i(l->dati[i].vi);
#endif
					break;
				case REALE:
				case VETTd:
#ifdef MDEBUG
					if (l->dati[i].vd != NULL) {
						nome_elem = ((Allocazione *) l->dati[i].vd->mem->data)->nome->str;
						l->dati[i].vd = _cancellaVett_d(l->dati[i].vd, nome_elem, nomefile, linea);
					}
#else
					CANCELLAv_d(l->dati[i].vd);
#endif
					break;
				case MATRi:
#ifdef MDEBUG
					if (l->dati[i].mi != NULL) {
						nome_elem = ((Allocazione *) l->dati[i].mi->mem->data)->nome->str;
						l->dati[i].mi = _cancellaMatr_i(l->dati[i].mi, nome_elem, nomefile, linea);
					}
#else
					CANCELLAm_i(l->dati[i].mi);
#endif
					break;
				case MATRd:
#ifdef MDEBUG
					if (l->dati[i].md != NULL) {
						nome_elem = ((Allocazione *) l->dati[i].md->mem->data)->nome->str;
						l->dati[i].md = _cancellaMatr_d(l->dati[i].md, nome_elem, nomefile, linea);
					}
#else
					CANCELLAm_d(l->dati[i].md);
#endif
					break;
				case STRINGA:
					CANCELLAstr(l->dati[i].str);
					break;
			}
		}
	}
#ifdef MDEBUG
	_CANCELLAstr(mem->nome);
	_CANCELLAstr(mem->file_da);
	_CANCELLAstr(mem->file_a);
	_CANCELLAstr(mem->file_prima);
	_CANCELLAstr(mem->file_ultima);
	allocLst = g_list_remove(allocLst, mem);
	libera(mem);
	mem = NULL;
	_stampa_lista_lst("allocLst-: ", allocLst);
#endif
	libera(l->dati);
	l->dati = NULL;
	libera(l->tipi);
	l->tipi = NULL;
	libera(l);
	return NULL;
}

#ifdef MDEBUG
	SEXP _daLISTA(LISTA *l, int *nProtected, const char *nomefile, int linea)
#else
	SEXP _daLISTA(LISTA *l, int *nProtected)
#endif
{
	int i;
	SEXP tmp = NULL, ris;
#ifdef MDEBUG
	char *nome_lista;
#endif

#ifdef FDEBUG
	fprintf(fp_fdbg, "Trasformo la lista\n");
#endif
	CONTROLLA(l != NULL);
#ifdef FDEBUG
	_StampaLista(l);
#endif
	PROTECT(ris = allocVector(VECSXP,  l->dim));
	(*nProtected)++;
	for (i = 0; i < l->dim; i++) {
		switch (l->tipi[i]) {
			case INTERO:
			case VETTi:
#ifdef MDEBUG
				tmp = _daVETTORE_i(l->dati[i].vi, nProtected, nomefile, linea);
#else
				tmp = daVETTORE_i(l->dati[i].vi, nProtected);
#endif
				break;
			case REALE:
			case VETTd:
#ifdef MDEBUG
				tmp = _daVETTORE_d(l->dati[i].vd, nProtected, nomefile, linea);
#else
				tmp = daVETTORE_d(l->dati[i].vd, nProtected);
#endif
				break;
			case MATRi:
#ifdef MDEBUG
				tmp = _daMATRICE_i(l->dati[i].mi, nProtected, nomefile, linea);
#else
				tmp = daMATRICE_i(l->dati[i].mi, nProtected);
#endif
				break;
			case MATRd:
#ifdef MDEBUG
				tmp = _daMATRICE_d(l->dati[i].md, nProtected, nomefile, linea);
#else
				tmp = daMATRICE_d(l->dati[i].md, nProtected);
#endif
				break;
			case STRINGA:
				tmp = daSTRINGA(l->dati[i].str, nProtected);
				break;
		}
		SET_VECTOR_ELT(ris, i, tmp);
	}
#ifdef MDEBUG
	nome_lista = ((Allocazione *) l->mem->data)->nome->str;
	_CancellaLista(l, false, nome_lista, nomefile, linea);
#else
	CancellaLISTA(l, false);
#endif
	return ris;
}

void _infoLista(const LISTA *l)
{
#ifdef MDEBUG
	Allocazione *mem;

	if (l == NULL) {
		Rprintf("lista non trovata o mai allocata!\n");
		R_FlushConsole();
		return;
	}
	mem = (Allocazione *) l->mem->data;
	Rprintf("Lista '%s' (%d):\n", mem->nome->str, l);
	Rprintf("\tnumero di elementi: %d\n", l->dim);
	Rprintf("\tdimensione allocata elementi: %d\n", l->mia_alloc);
	Rprintf("\tdata in ingresso?: %d\n", l->r);
	Rprintf("\tallocata da: %s # %d\n", mem->file_da->str, mem->linea_da);
	if (mem->linea_a == 0)
		Rprintf("\tancora allocata\n");
	else
		Rprintf("\tdisallocata da: %s # %d\n", mem->file_a->str, mem->linea_a);
	Rprintf("\tindice progressivo: %d\n", mem->indx);
	Rprintf("\tdimensione massima finora: %d\n", mem->max_dim1);
	Rprintf("\triallocazioni finora: %d\n", mem->rialloc);
	Rprintf("\tprima riga: %s # %d\n", mem->file_prima->str, mem->prima_linea);
	Rprintf("\tultima riga: %s # %d\n\n", mem->file_ultima->str, mem->ultima_linea);
	R_FlushConsole();
#endif
}

#ifdef MDEBUG
	void  _controllaCanc(const char *nomefile, int linea)
#else
	void  _controllaCanc()
#endif
{
#ifdef MDEBUG
	_controllaCanc_i(nomefile, linea);
	_controllaCanc_d(nomefile, linea);
	_controllaCanc_lst(nomefile, linea);
#endif
}

// interna
#ifdef MDEBUG
	void  _controllaCanc_lst(const char *nomefile, int linea)
#else
	void  _controllaCanc_lst()
#endif
{
#ifdef MDEBUG
	int ok = 1;
	Allocazione *t;

	fprintf(fp_mdbg_lst, "--------------------------------------\n");
	fprintf(fpm_lst, "--------------------------------------\n");
	while (allocLst != NULL) {
		t = (Allocazione *) allocLst->data;
		Rprintf("La lista '%s' (allocata alla linea %s # %d, indice %d) non e` stata ancora disallocata: lo faccio adesso.\n", t->nome->str, t->file_da->str, t->linea_da, t->indx);
		t->indir = (size_t) (LISTA *) _CancellaLista((LISTA *) t->indir, false, t->nome->str, nomefile, linea);
		ok = 0;
		_CANCELLAstr(t->nome);
		_CANCELLAstr(t->file_da);
		_CANCELLAstr(t->file_a);
		_CANCELLAstr(t->file_prima);
		_CANCELLAstr(t->file_ultima);
		libera(t);
		t = NULL;
	}
	g_list_free(allocLst);
	fprintf(fpm_lst, "\n");
	fclose(fpm_lst);
	if (!ok)
		fprintf(fp_mdbg_lst, "ATTENZIONE: per il bilanciamento delle allocazioni di elementi di tipo 'TIPO' sono state necessarie disallocazioni automatiche.\n");
	fprintf(fp_mdbg_lst, "--------------------------------------\n");
	if (fp_mdbg_lst != NULL) {
		fclose(fp_mdbg_lst);
		fp_mdbg_lst = NULL;
	}
#ifdef FDEBUG
	if (fp_fdbg != NULL) {
		fclose(fp_fdbg);
		fp_fdbg = NULL;
	}
#endif
#ifdef DET
	if (fp_det != NULL) {
		fclose(fp_det);
		fp_det = NULL;
	}
#endif
#endif
	return;
}

// qui in realta` non faccio niente!
#ifdef MDEBUG
	void _controllaLLista(const LISTA *l, int e, const char *nome, const char *nomefile, int linea)
#else
	void _controllaLLista(const LISTA *l, int e, const char *nome)
#endif
{
#ifdef MDEBUG
	Allocazione *mem;
	char err[128];

	if (l == NULL || l->mem == NULL) {
		snprintf(err, 128, "la lista '%s' (linea %s # %d) non e` allocata!\n", nome, nomefile, linea);
		error(err);
	}
	if (e < 0 || e >= l->dim) {
		snprintf(err, 128, "tentativo di accedere in lettura all'elemento [%d] della lista '%s' [%d] (linea %s # %d)!\n", e + 1, nome, l->dim, nomefile, linea);
		error(err);
	}
	mem = (Allocazione *) l->mem->data;
	if (mem->max_dim1 < e) {
		mem->max_dim1 = e + 1;
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
	return;
#endif
}

// qui in realta` non faccio niente!
#ifdef MDEBUG
	void  _controllaSLista(const LISTA *l, int e, const char *nome, const char *nomefile, int linea)
#else
	void  _controllaSLista(const LISTA *l, int e, const char *nome)
#endif
{
#ifdef MDEBUG
	Allocazione *mem;
	char err[128];

	if (l == NULL || l->mem == NULL) {
		snprintf(err, 128, "la lista '%s' (linea %s # %d) non e` allocata!\n", nome, nomefile, linea);
		error(err);
	}
	if (e < 0 || e >= l->mia_alloc) {
		snprintf(err, 128, "tentativo di accedere in scrittura all'elemento [%d] della lista '%s' [allocata: %d] (linea %s # %d)!\n", e + 1, nome, l->mia_alloc, nomefile, linea);
		error(err);
	}
	mem = (Allocazione *) l->mem->data;
	if (mem->max_dim1 < e) {
		mem->max_dim1 = e + 1;
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
}

void StrBilanciam()
{
#ifdef FDEBUG
	fprintf(fp_fdbg, "\nStringhe allocate: %d (%d interne); stringhe deallocate: %d (%d interne)\n\n", tot_str_alloc, tot_str_intalloc, tot_str_dealloc, tot_str_intdealloc);
#endif
}

#ifdef MDEBUG
	SEXP _daSTRINGA(GString *str, int *nProtected, const char *nome)
#else
	SEXP _daSTRINGA(GString *str, int *nProtected)
#endif
{
	SEXP ris;

#ifdef FDEBUG
	fprintf(fp_fdbg, "Trasformo la stringa '%s' che vale '%s'\n", nome, str->str);
#endif

	PROTECT(ris = allocVector(STRSXP,  1));
	(*nProtected)++;
	SET_STRING_ELT(ris,  0,  mkChar(str->str));
	CANCELLAstr(str);
	return ris;
}

#ifdef MDEBUG
	GString *_inSTRINGA(SEXP s, int *nProtected, const char *nome)
#else
	GString *_inSTRINGA(SEXP s, int *nProtected)
#endif
{
	GString *str = NULL;

	if (isNull(s))
		return NULL;
	PROTECT(s = AS_CHARACTER(s));
	(*nProtected)++;
	// sarebbe interna ma va cancellata all'esterno!
	CREAstr(str, CHAR(STRING_ELT(s, 0)));

#ifdef FDEBUG
	fprintf(fp_fdbg, "Ho trasformato la stringa '%s' che vale '%s'\n", nome, str->str);
#endif
	return str;
}

#ifdef MDEBUG
	MATRICEd * _moltiplica_mm_di(MATRICEd *ris, const MATRICEd *m1, const MATRICEi *m2, const char *nome, const char *nomefile, int linea)
#else
	MATRICEd * _moltiplica_mm_di(MATRICEd *ris, const MATRICEd *m1, const MATRICEi *m2)
#endif
{
	int r, c;

#ifdef FDEBUG
	fprintf(fp_fdbg, "%d: moltiplica_mm_di\n", linea);
#endif
	CONTROLLA(m1 != NULL && m2 != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m1);
	_StampaMatr_i(m2);
#endif
	CONTROLLA(m1->nr == m2->nr && m1->nc == m2->nc);
	_CREAm_d(ris, m1->nr, m1->nc);
	for (r = 1; r <= LENGTHm1_d(m1); r++) {
		for (c = 1; c <= LENGTHm2_d(m1); c++)
			_ASSEGNAm_d(ris, r, c, _ACCEDIm_d(m1, r, c) * _ACCEDIm_i(m2, r, c));
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaMatr_d(ris);
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif
	return ris;
}
