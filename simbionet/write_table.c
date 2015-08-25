#include "write_table.h"


void write_m_i(const char *nomefile, MATRICEi *m)
{
	int i = 0, j;
	FILE *fp;
	char err[256];

	_Intestazione("\n*** write_m_i ***\n");

	assert(m != NULL);
#ifdef FDEBUG
	_StampaMatr_i(m);
#endif
	fp = fopen(nomefile, "w");
	if (!fp) {
		snprintf(err, 256, "Cannot write file '%s'", nomefile);
		error(err);
		return;
	}
	for (j = 1; j <= LENGTHm2_i(m); j++)
		fprintf(fp, "\tV%d", j);
	fprintf(fp, "\n");
	for (i = 1; i <= LENGTHm1_i(m); i++) {
		fprintf(fp, "%d", i);
		for (j = 1; j <= LENGTHm2_i(m); j++)
			fprintf(fp, "\t%d", ACCEDIm_i(m, i, j));
		fprintf(fp, "\n");
	}
	fclose(fp);

	StrBilanciam();
}

SEXP write_m_int(SEXP nome, SEXP m)
{
	int nProtected = 0;
	GString *nome1;
	MATRICEi *m1;

	_InitDbg(false, false, false);

	_Intestazione("\n*** write_m_int ***\n");

	nome1 = inSTRINGA(nome, &nProtected, "nomefile");
	m1 = inMATRICE_i(m, &nProtected);

	write_m_i(nome1->str, m1);

	CANCELLAstr(nome1);
	CANCELLAm_i(m1);

	StrBilanciam();
	ControllaCanc();

	UNPROTECT(nProtected);

	return m;
}

void write_m_d(const char *nomefile, MATRICEd *m)
{
	int i = 0, j;
	FILE *fp;
	char err[256];

	_Intestazione("\n*** write_m_d ***\n");

	assert(m != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
#endif
	fp = fopen(nomefile, "w");
	if (!fp) {
		snprintf(err, 256, "Cannot write file '%s'", nomefile);
		error(err);
		return;
	}
	for (j = 1; j <= LENGTHm2_d(m); j++)
		fprintf(fp, "\tV%d", j);
	fprintf(fp, "\n");
	for (i = 1; i <= LENGTHm1_d(m); i++) {
		fprintf(fp, "%d", i);
		for (j = 1; j <= LENGTHm2_d(m); j++)
			fprintf(fp, "\t%e", ACCEDIm_d(m, i, j));
		fprintf(fp, "\n");
	}
	fclose(fp);

	StrBilanciam();
}

SEXP write_m_double(SEXP nome, SEXP m)
{
	int nProtected = 0;
	GString *nome1;
	MATRICEd *m1;

	_InitDbg(false, false, false);

	_Intestazione("\n*** write_m_double ***\n");

	nome1 = inSTRINGA(nome, &nProtected, "nomefile");
	m1 = inMATRICE_d(m, &nProtected);

	write_m_d(nome1->str, m1);

	CANCELLAstr(nome1);
	CANCELLAm_d(m1);

	StrBilanciam();
	ControllaCanc();

	UNPROTECT(nProtected);

	return m;
}

//~ #ifdef MDEBUG
	//~ void _write_mn_d(const GString *nome, MATRICEd *m, const LISTA *nomi)
//~ #else
	//~ void _write_mn_d(const GString *nome, MATRICEd *m, const LISTA *nomi)
//~ #endif
//~ {
	//~ int i = 0, j;
	//~ FILE *fp;

//~ #ifdef FDEBUG
	//~ fprintf(fp_fdbg, "*** _write_mn_d ***\n");
	//~ fprintf(fp_mdbg_i, "*** _write_mn_d ***\n");
	//~ fprintf(fp_mdbg_d, "*** _write_mn_d ***\n");
//~ #endif
	//~ assert(m != NULL);
//~ #ifdef FDEBUG
	//~ _StampaMatr_d(m);
	//~ _stampa_lista_str("nomi: ", nomi);
//~ #endif
	//~ fp = fopen(nome->str, "w");
	//~ for (j = 0; j < nomi->dim; j++)
		//~ fprintf(fp, "\t%s", nomi->dati[j].str->str);
	//~ fprintf(fp, "\n");
	//~ for (i = 1; i <= LENGTHm1_d(m); i++) {
		//~ fprintf(fp, "%d", i);
		//~ for (j = 1; j < LENGTHm2_d(m); j++)
			//~ fprintf(fp, "\t%f", ACCEDIm_d(m, i, j));
		//~ fprintf(fp, "\n");
	//~ }
	//~ fclose(fp);
//~ }

//~ SEXP write_param_double(SEXP nome, SEXP m, SEXP nomi)
//~ {
	//~ int nProtected = 0, ll, i;
	//~ GString *nome1;
	//~ GString **nomi_elem;
	//~ LISTA *nomi1;
	//~ MATRICEd *m1;
	//~ enum TIPO *tipi;

	//~ _InitDbg_i(0);
	//~ _InitDbg_d(0);

//~ #ifdef FDEBUG
	//~ fprintf(fp_fdbg, "*** write_param_double ***\n");
	//~ fprintf(fp_mdbg_i, "*** write_param_double ***\n");
	//~ fprintf(fp_mdbg_d, "*** write_param_double ***\n");
//~ #endif

	//~ nome1 = inSTRINGA(nome, &nProtected, "nomefile");
	//~ m1 = inMATRICE_d(m, &nProtected);
	//~ ll = LENGTH(nomi);
	//~ tipi = (enum TIPO *) g_mallocll);
	//~ nomi_elem = (GString **) g_mallocll);
	//~ for (i = 0; i < ll; i++) {
		//~ tipi[i] = STRINGA;
		//~ nomi_elem[i] = g_string_new("nome");
	//~ }
	//~ nomi1 = inLISTA(nomi, &nProtected, ll, tipi, nomi_elem);
	//~ _write_mn_d(nome1, m1, nomi1);

	//~ CANCELLAm_d(m1);
	//~ CancellaLISTA(nomi1);

	//~ UNPROTECT(nProtected);
	//~ controllaCanc_i();
	//~ controllaCanc_d();
	//~ return m;
//~ }

void write_vn_d(const char *nomefile, VETTOREd *v, const char *nome)
{
	int i = 0;
	FILE *fp;
	char err[256];

	_Intestazione("\n*** write_vn_d ***\n");

	assert(v != NULL);
#ifdef FDEBUG
	_StampaVett_d(v);
#endif
	fp = fopen(nomefile, "w");
	if (!fp) {
		snprintf(err, 256, "Cannot write file '%s'", nomefile);
		error(err);
		return;
	}
	fprintf(fp, "\t%s", nome);
	fprintf(fp, "\n");
	for (i = 1; i <= LENGTHv_d(v); i++)
		fprintf(fp, "%d\t%e\n", i, ACCEDIv_d(v, i));
	fclose(fp);

	StrBilanciam();
}

SEXP write_vparam_double(SEXP nome_file, SEXP v, SEXP nome)
{
	int nProtected = 0;
	GString *nome_file1, *nome1;
	VETTOREd *v1;

	_InitDbg(false, false, false);

	_Intestazione("\n*** write_vparam_double ***\n");

	nome_file1 = inSTRINGA(nome_file, &nProtected, "nomefile");
	v1 = inVETTORE_d(v, &nProtected);
	nome1 = inSTRINGA(nome, &nProtected, "nome");

	write_vn_d(nome_file1->str, v1, nome1->str);

	CANCELLAstr(nome_file1);
	CANCELLAstr(nome1);
	CANCELLAv_d(v1);

	StrBilanciam();
	ControllaCanc();

	UNPROTECT(nProtected);

	return NULL;
}

void write_mn_d(const char *nomefile, MATRICEd *m, const char *nome)
{
	int i = 0;
	FILE *fp;
	char err[256];

	_Intestazione("\n*** write_mn_d ***\n");

	assert(m != NULL);
#ifdef FDEBUG
	_StampaMatr_d(m);
#endif
	fp = fopen(nomefile, "w");
	if (!fp) {
		snprintf(err, 256, "Cannot write file '%s'", nome);
		error(err);
		return;
	}
	fprintf(fp, "\t%s", nome);
	fprintf(fp, "\n");
	for (i = 1; i <= LENGTHm1_d(m); i++)
		fprintf(fp, "%d\t%e\n", i, ACCEDIm_d(m, i, 1));
	fclose(fp);

	StrBilanciam();
}

SEXP write_mparam_double(SEXP nome_file, SEXP m, SEXP nome)
{
	int nProtected = 0;
	GString *nome_file1, *nome1;
	MATRICEd *m1;

	_InitDbg(false, false, false);

	_Intestazione("\n*** write_mparam_double ***\n");

	nome_file1 = inSTRINGA(nome_file, &nProtected, "nomefile");
	m1 = inMATRICE_d(m, &nProtected);
	nome1 = inSTRINGA(nome, &nProtected, "nome");

	write_mn_d(nome_file1->str, m1, nome1->str);

	CANCELLAstr(nome_file1);
	CANCELLAstr(nome1);
	CANCELLAm_d(m1);

	StrBilanciam();
	ControllaCanc();

	UNPROTECT(nProtected);

	return NULL;
}
