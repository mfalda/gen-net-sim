#define MEM_PROFILE

#include "interfacce.h"

#define g_v2 globali.ret_vettore.v2
#define g_v3 globali.ret_vettore.v3

#define g_v21 globali.ret_lista.v2
#define g_m1 globali.ret_lista.m1

VETTOREd *ret_vettore1(VETTOREd *ris, const VETTOREi *v1, double num)
{
	//~ VETTOREd *v2 = globali.interfacce.v2, *v3 = globali.interfacce.v3;

#ifdef FDEBUG
	// salvo v1 e lo ricontrollo poi per sicurezza
	VETTOREi *v1_org = NULL;

	// scrivo le intestazioni per le tracce (tanto e` solo in debug, quindi e` lo stesso se la metto qui dentro)
	_Intestazione("\n*** ret_vettore1 ***\n");
	v1_org = copia_v_i(v1_org, v1, 1, LENGTHv_i(v1));
#endif

	// informazioni sul vettore g_v3 a video (la prima volta non esiste, perche´ sara` allocato alla riga 34)
	_infoVett_d(g_v3);
	// errore? si`, mi da` un warning!!
	//ASSEGNAv_i(v1, 1.0, 10);
	// creo un vettore g_v2 della stessa lunghezza di v1
	CREAv_d(g_v2, LENGTHv_i(v1));
	// leggo (eventualmente) i valori salvati da una versione R precedente se e` stata definita DEF
	g_v2 = runif_s(g_v2, LENGTHv_d(g_v2), 0, 10, "ret_vettore1");
	// stampo il vettore g_v2 nel file "passi.txt"
	_StampaVett_d(g_v2);
	CREAv_d(g_v3, 0); // se metto zero NON lo cancello!
	// stampo il vettore g_v2 nel file "passi.txt"
	_StampaVett_d(g_v3);
	g_v3 = promuovi_i(g_v3, v1);
	ris = somma_vs_d(ris, g_v3, num);

	CREAv_d(g_v3, 10);
	InitVett_d(g_v3, 10.0);
	_StampaVett_d(g_v3);
	//~ CANCELLAv_d(v3);
	//~ CANCELLAv_d(v2);

	// stampo le statistiche sulle stringhe allocate
	StrBilanciam();

#ifdef FDEBUG
	// prima di uscire mi assicuro di non aver alterato v1
	assert(Uguali_v_i(v1, v1_org) == true);
	CANCELLAv_i(v1_org);
#endif
	return ris;
}

SEXP ret_vettore(SEXP vett, SEXP num)
{
	VETTOREi *v1 = NULL;
	VETTOREd *ris1 = NULL;
	int i, nProtected = 0;
	double num1;
	SEXP ris;

	_InitDbg(false, false, false);

	_Intestazione("\n*** ret_vettore ***\n");

	v1 = inVETTORE_i(vett, &nProtected);
	num1 = NUMERIC_VALUE(num);
	InitGlobali();
	for (i = 0; i < 3; i++)
		ris1 = ret_vettore1(ris1, v1, num1);
	CancGlobali();

	ris = daVETTORE_d(ris1, &nProtected);

	CANCELLAv_i(v1);

	StrBilanciam();
	ControllaCanc();

	UNPROTECT(nProtected);

	return ris;
}

MATRICEd *ret_matrice1(MATRICEd *ris, const MATRICEd *m1)
{
	MATRICEd *m3 = NULL;

#ifdef FDEBUG
	// salvo m1 e la ricontrollo poi per sicurezza
	MATRICEd *m1_org = NULL;

	_Intestazione("\n*** ret_matrice1 ***\n");
	m1_org = copia_m_d(m1_org, m1);
#endif

	ris = abs_m_d(ris, m1);

	CREAm_d(m3, 10, 4);
	InitMatr_d(m3, 10.5);
	_StampaMatr_d(m3);
	CANCELLAm_d(m3);

	StrBilanciam();

#ifdef FDEBUG
	// prima di uscire mi assicuro di non aver alterato m1
	assert(Uguali_m_d(m1, m1_org) == true);
	CANCELLAm_d(m1_org);
#endif

	return ris;
}


SEXP ret_matrice(SEXP matr)
{
	MATRICEd *m1 = NULL, *ris1 =NULL;
	int nProtected = 0;
	SEXP ris;

	_InitDbg(false, false, false);

	_Intestazione("\n*** ret_matrice ***\n");

	InitGlobali();
	m1 = inMATRICE_d(matr, &nProtected);
	ris1 = ret_matrice1(ris1, m1);
	ris = daMATRICE_d(ris1, &nProtected);
	CancGlobali();

	CANCELLAm_d(m1);

	StrBilanciam();
	ControllaCanc();

	UNPROTECT(nProtected);

	return ris;
}

LISTA *ret_lista1(LISTA *ris, const LISTA *l1)
{
	enum TIPO tipi[3];
	VETTOREd *v1 = NULL;
	GString *str1;
#ifdef FDEBUG

	// salvo l1 e la ricontrollo poi per sicurezza
	LISTA *l1_org = NULL;

	_Intestazione("\n*** ret_lista1 ***\n");
	l1_org = CopiaLISTA(l1_org, l1);
#endif

	// controllo di poter accedere in lettura al primo elemento
	CtrlLlst(l1, 1);
	// divido il primo elemento per 2 e pongo il risultato nel vettore v1
	v1 = dividi_vs_d(v1, ACCEDIlst(l1, 1, vd), 2.0);
	g_v21 = arrotonda_v_d(g_v21, v1);
	// controllo di poter accedere in lettura al secondo elemento
	CtrlLlst(l1, 2);
	g_m1 = abs_m_i(g_m1, ACCEDIlst(l1, 2, mi));
	// controllo di poter accedere in lettura al terzo elemento
	CtrlLlst(l1, 3);
	// creo una stringa pari al terzo elemento
	CREAstr(str1, ACCEDIlst(l1, 3, str)->str);
	g_string_ascii_up(str1);

	// creo la stringa di uscita
	tipi[0] = VETTi;
	tipi[1] = MATRi;
	tipi[2] = STRINGA;
	CreaLISTA(ris, tipi, 3);
	// controllo di poter accedere in scrittura al primo elemento
	CtrlSlst(ris, 1);
	// stampo il primo elemento (ancora nullo)
	_StampaVett_i(ACCEDIlst(ris, 1, vi));
	// creo un nuovo vettore per il primo elemento del risultato
	CREAv_i(ACCEDIlst(ris, 1, vi), LENGTHv_i(g_v21));
	// COPIO il vettore g_v21 nel primo elemento appena creato (NB non potrei usare copia_v_i)
	ACCEDIlst(ris, 1, vi) = copia_v_i(ACCEDIlst(ris, 1, vi), g_v21, 1, LENGTHv_i(g_v21));
	//~ for (i = 1; i <= LENGTHv_i(g_v21); i++)
		//~ ASSEGNAv_i(ACCEDIlst(ris, 1, vi), i, ACCEDIv_i(g_v21, i));
	// controllo di poter accedere in scrittura al secondo elemento
	CtrlSlst(ris, 2);
	// vi assegno la matrice g_m1 (solo il puntatore!!)
	ASSEGNAlst(ris, 2, mi, g_m1);
	// controllo di poter accedere in scrittura al terzo elemento
	CtrlSlst(ris, 3);
	// vi assegno la stringa str1 (solo il puntatore!, e` locale ma in memoria dinamica)
	ASSEGNAlst(ris, 3, str, str1);

	_StampaLista(ris);

	// v1 e` temporaneo e lo cancello (invece g_v21 e` globale e verra` cancellato alla fine)
	CANCELLAv_d(v1);

	StrBilanciam();

#ifdef FDEBUG
	// prima di uscire mi assicuro di non aver alterato l1
	assert(UGUALIlst(l1, l1_org) == true);
	CancellaLISTA(l1_org, true);
#endif

	return ris;
}

SEXP ret_lista(SEXP lista)
{
	enum TIPO *tipi; // qui saprei che sono solo tre, ma e` per mostrare l'allocazione dinamica (NB: se sono tutti VETTi basta passare NULL)
	LISTA *l1 = NULL, *ris1 = NULL;
#ifndef NDEBUG
	GString **nomi;
#endif
	int i;
	int nProtected = 0;
	SEXP ris;

	_InitDbg(false, false, false);

	_Intestazione("\n*** ret_lista ***\n");

	tipi = mia_alloc(3, enum TIPO);
	if (3 > 0 && tipi == NULL) { // sempre vero che 3 > 0, ma e` solo per mostrare il modello generale
		Rprintf("Not enough memory (interfacce # %d, tipi)", __LINE__ - 2);
		error("");
	}
	tipi[0] = VETTd;
	tipi[1] = MATRi;
	tipi[2] = STRINGA;
#ifndef NDEBUG
	nomi = mia_alloc(3, GString *);
	if (3 > 0 && nomi == NULL) {
		Rprintf("Not enough memory (interfacce # %d, nomi)", __LINE__ - 2);
		error("");
	}
	CREAstr(nomi[0], "elemento 1");
	CREAstr(nomi[1], "elemento 2");
	CREAstr(nomi[2], "elemento 3");
#endif
	l1 = inLISTA(lista, &nProtected, 3, tipi, nomi);
	libera(tipi);
#ifndef NDEBUG
	CANCELLAstr(nomi[0]);
	CANCELLAstr(nomi[1]);
	CANCELLAstr(nomi[2]);
	libera(nomi);
#endif
	InitGlobali();
	for (i = 0; i < 3; i++)
		ris1 = ret_lista1(ris1, l1);

	ris = daLISTA(ris1, &nProtected);

	CancellaLISTA(l1, true);
	CancGlobali();

	// qui risultano due allocazioni di stringhe in piu` (in effetti le perdo perche´ non le ricopio!)
	StrBilanciam();
	ControllaCanc();

	UNPROTECT(nProtected);

	return ris;
}
