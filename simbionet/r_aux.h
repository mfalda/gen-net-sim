/*! \file r_aux.h
	\brief Header \em file principale.

	include automaticamente gli altri \em header \em file principali.

	\defgroup oggetti Gestione degli oggetti
	\defgroup convvett Conversioni per vettori
	\defgroup convmatr Conversioni per matrici
	\defgroup convlst Conversioni per liste
	\defgroup convstr Conversioni per stringhe
	\defgroup debug Debugging
	\defgroup wvett Which per vettori
	\defgroup wmatr Which per matrici
	\defgroup setvett Assegnamenti per vettori
	\defgroup setmatr Assegnamenti per matrici
	\defgroup opnum Operazioni su numeri scalari
	\defgroup opvett Operazioni su vettori
	\defgroup opmatr Operazioni su matrici
	\defgroup oplst Operazioni su liste
	\defgroup sequenze Sequenze su vettori
	\defgroup aus Funzioni ausiliarie per Netsim
	\defgroup aus1 Funzioni ausiliarie per HMM

*/

// convenzione per i nomi: [_]nome[1]_inp_funzetipo

#ifndef R_AUX_H
#define R_AUX_H

//~ #define DET
//~ #define NDEBUG
//~ #define MDEBUG
//~ #define FDEBUG

#ifdef NDEBUG
#warning "versione di release"
#else
#warning "versione di debug"
#define FDEBUG
#define MDEBUG
#endif

#ifdef DET
#warning "versione deterministica"
#endif

#ifdef __WIN32__
#pragma message "definisco bool"
typedef unsigned int bool;
#define false   0
#define true    1
#else
#include <stdbool.h>
#endif

#include <glib.h>

#include <stdio.h>
#include <assert.h>
#include <stddef.h>
#include <string.h>

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Random.h>
#include <Rmath.h>
#include <R_ext/Arith.h> // nan
#include <R_ext/Utils.h> // sort
#include <errno.h>

#define MAX_NOME 64
#define MAX_STR 5

#define ISNAi(x) (x == NA_INTEGER)
#define ISNANi(x) (x == NA_INTEGER)

// malloc sarebbe piu` veloce, ma e` meno sicuro
#define mia_alloc(A, B) (B*) malloc(A * sizeof(B))
//~ #define alloc(A, B) Calloc(A, B)
#define libera(A) if (A != NULL) free(A)

#ifdef FDEBUG
	#define CONTROLLA(E) _errore((E), #E, nomefile, linea)
#else
	#define CONTROLLA(E) assert(E)
#endif

int id;
int tot_str_alloc, tot_str_dealloc;
int tot_str_intalloc, tot_str_intdealloc;
FILE *fp_fdbg;
FILE *fp_mdbg_lst;
FILE *fpm_lst;
FILE *fp_det;

typedef struct {
	GString *nome;
	size_t indir;
	GString *file_da;
	int linea_da;
	GString *file_a;
	int linea_a;
	int indx;
	int max_dim1;
	int max_dim2;
	int rialloc;
	GString *file_ultima;
	int ultima_linea;
	GString *file_prima;
	int prima_linea;
	bool in_uso;
} Allocazione;

//! Struttura per i vettori di tipo 'int'.
typedef struct {
	int dim; /*!< dimensioni del vettore */
	int alloc; /*!< uso interno */
	bool r; /*!< uso interno */
	int *dati; /*!< elementi del vettore */
	GList *mem; /*!< uso interno */
} VETTOREi;

//! Struttura per i vettori di tipo 'double'.
typedef struct {
	int dim; /*!< dimensioni del vettore */
	int alloc; /*!< uso interno */
	bool r; /*!< uso interno */
	double *dati; /*!< elementi del vettore */
	GList *mem; /*!< uso interno */
} VETTOREd;

GList *allocLst;

#include "r_aux_i.h"
#include "r_aux_d.h"

//! Tipi di elementi per le liste.
enum TIPO {
	INTERO, /*!< un numero intero */
	REALE, /*!< un numero reale */
	VETTi, /*!< un vettore di interi */
	VETTd, /*!< un vettore di reali */
	MATRi, /*!< una matrice di interi */
	MATRd, /*!< una matrice di reali */
	STRINGA /*!< una stringa */
};

//! Unione di supporto
union Dati {
	VETTOREi *vi;  /*!< dati di tipo VETTOREi */
	VETTOREd *vd;  /*!< dati di tipo VETTOREd */
	MATRICEi *mi;  /*!< dati di tipo MATRICEi */
	MATRICEd *md;  /*!< dati di tipo MATRICEd */
	GString *str;  /*!< dati di tipo stringa */
};

//! Struttura per le liste.
typedef struct {
	enum TIPO *tipi;  /*!< tipo di elemento */
	union Dati *dati; /*!< dati */
	bool r; /*!< uso interno */
	int dim; /*!< dimensioni della lista */
	int alloc; /*!< dimensioni allocate */
	GList *mem; /*!< uso interno */
} LISTA;

//! Controlla il segno di un numero intero.
/*!
	\ingroup opnum

	\param x il numero di riferimento

	\return il segno di x.
*/
int Segno(int x);

//! Controlla se due numeri sono uguali tenendo conto di approssimazioni
/*!
	\ingroup opnum

   \param x il primo numero
   \param y il secondo numero

	\return vero se uguali, falso se diversi
*/
bool Uguale(double x, double y);

//interna
void _errore(int val, const char *msg, const char *nomefile, int linea);

//! Stampa una lista.
/*!
	\ingroup debug

	\param l la lista da stampare
*/
void _StampaLista(const LISTA *l); // NB: "TIPO const*" e` lo stesso che const TIPO* inoltre, nelle funzioni, "const TIPO* const" e` ridondante, dato che il puntatore viene copiato!

//! Stampa una lista come R.
/*!
	\ingroup debug

	\param l la lista da stampare
*/
void _StampaRawLista(const LISTA *l);

// interna
void _stampa_lista_lst(const char *pref, GList *lista);

//! Scrive un messaggio di intestazione in tutti i file di debug
/*!
	\ingroup debug

	\param msg il messaggio da scrivere

*/
void _Intestazione(const char *msg);

//! Inizializza il \em debugging per gli oggetti di tipo vettore, matrice e lista
/*!
	\ingroup debug

	\param stdout1_v se vero invia a STDOUT le informazioni sui vettori, altrimenti su file
	\param stdout1_m se vero invia a STDOUT le informazioni sulle matrici, altrimenti su file
	\param stdout1_lst se vero invia a STDOUT le informazioni sulle liste, altrimenti su file

*/
void _InitDbg(bool stdout1_v, bool stdout1_m, bool stdout1_lst);

void _infoLista(const LISTA *l);

//! Stampa un riepilogo delle stringhe usate (solo in MDEBUG).
/*!
	\ingroup debug

	\note possono risultare piu` deallocazioni che allocazioni: capita quando le stringhe non sono copie ma puntatori (tuttavia non e` pericoloso). Usare prima di ControllaCanc, percheï¿½ le serve un handle chiuso da quella
*/
void StrBilanciam();

//! Stampa un riepilogo della memoria usata (solo in MDEBUG).
/*!
	\ingroup debug
*/
#ifdef MDEBUG
	#define ControllaCanc() _controllaCanc(__FILE__, __LINE__)
#else
	#define ControllaCanc() _controllaCanc()
#endif
#ifdef MDEBUG
	void  _controllaCanc(const char *nomefile, int linea);
#else
	void  _controllaCanc();
#endif

#ifdef MDEBUG
	#define controllaCanc_lst() _controllaCanc_lst(__FILE__, __LINE__)
#else
	#define controllaCanc_lst() _controllaCanc_lst()
#endif
#ifdef MDEBUG
	void  _controllaCanc_lst(const char *nomefile, int linea);
#else
	void  _controllaCanc_lst();
#endif

//! Lunghezza di una lista.
/*!
	\ingroup oggetti

   \param L la lista di riferimento
*/
#ifdef MDEBUG
	#define LENGTHlst(L) _lengthList(L, #L, __FILE__, __LINE__)
#else
	#define LENGTHlst(L) L->dim
#endif

//! predispone una LISTA da una espressione SEXP.
/*!
	\ingroup convlst

	\param s l'espressione SEXP da convertire
	\param nProtected un contatore da passare tra tutti gli oggetti costruiti nella stessa funzione esportata
	\param t un vettore con il tipo di elementi (NULL per creare una lista di vettori di interi)
	\param nomi un vettore con il nome degli elementi
	\param len la lunghezza del vettore t

	\return la lista LISTA da restituire
*/

#ifdef MDEBUG
	#define inLISTA(s, nProtected, len, t, nomi) _inLISTA(s, nProtected, len, t, nomi, #s, __FILE__, __LINE__)
#else
	#define inLISTA(s, nProtected, len, t, nomi) _inLISTA(s, nProtected, len, t)
#endif
#ifdef MDEBUG
	LISTA *_inLISTA(SEXP s, int *nProtected, int len, const enum TIPO *t, GString **nomi, const char *nome, const char *nomefile, int linea);
#else
	LISTA *_inLISTA(SEXP s, int *nProtected, int len, const enum TIPO *t);
#endif
//! crea una LISTA.
/*!
	\ingroup oggetti

	\param L un puntatore alla nuova lista
	\param T un vettore con il tipo di elementi (NULL per creare una lista di vettori di interi)
	\param len la lunghezza del vettore t

	\return il puntatore alla nuova lista
*/
#ifdef MDEBUG
	#define CreaLISTA(L, T, len) L = _CreaLista(L, T, len, #L, __FILE__, __LINE__)
	#define _CreaLISTA(L, T, len) L = _CreaLista(L, T, len, #L, nomefile, linea)
#else
	#define CreaLISTA(L, T, len) L = _CreaLista(L, T, len)
	#define _CreaLISTA(L, T, len) L = _CreaLista(L, T, len)
#endif
#ifdef MDEBUG
	LISTA *_CreaLista(LISTA *ris, const enum TIPO *t, int len, const char *nome, const char *nomefile, int linea);
#else
	LISTA *_CreaLista(LISTA *ris, const enum TIPO *t, int len);
#endif

//! cancella una LISTA (per inLISTA basta daLISTA).
/*!
	\ingroup oggetti

	\param L la lista da cancellare
	\param T cancellare tutto il contenuto?
*/
#ifdef MDEBUG
	#define CancellaLISTA(L, T) L = _CancellaLista(L, T, #L, __FILE__, __LINE__)
	#define _CancellaLISTA(L, T) L = _CancellaLista(L, T, nome, nomefile, linea)
#else
	#define CancellaLISTA(L, T) L = _CancellaLista(L, T)
	#define _CancellaLISTA(L, T) L = _CancellaLista(L, T)
#endif
#ifdef MDEBUG
	LISTA *_CancellaLista(LISTA *l, bool tutto, const char *nome,const char *nomefile, int linea);
#else
	LISTA *_CancellaLista(LISTA *l, bool tutto);
#endif

//! converte una LISTA in una espressione SEXP.
/*!
	\ingroup convlst

	\param v la lista da convertire
	\param nProtected un contatore da passare tra tutti gli oggetti costruiti nella stessa funzione esportata

	\return l'espressione convertita
*/
#ifdef MDEBUG
	#define daLISTA(v, nProtected) _daLISTA(v, nProtected, __FILE__, __LINE__)
#else
	#define daLISTA(v, nProtected) _daLISTA(v, nProtected)
#endif
#ifdef MDEBUG
	SEXP _daLISTA(LISTA *l, int *nProtected, const char *nomefile, int linea);
#else
	SEXP _daLISTA(LISTA *l, int *nProtected);
#endif

//! Controlla (in debug) che la posizione di lettura di un elemento di una lista sia corretto.
/*!
	\ingroup debug

   \param L la lista di riferimento
   \param I l'indice
*/
#ifdef MDEBUG
	#define CtrlLlst(L, I) _controllaLLista(L, I - 1, #L, __FILE__, __LINE__)
	#define _CtrlLlst(L, I) _controllaLLista(L, I - 1, #L, nomefile, linea)
#else
	#define CtrlLlst(L, I)
	#define _CtrlLlst(L, I)
#endif
#ifdef MDEBUG
	void _controllaLLista(const LISTA *l, int e, const char *nome, const char *nomefile, int linea);
#else
	void _controllaLLista(const LISTA *l, int e, const char *nome);
#endif

//! Accede ad un elemento di una lista.
/*!
	\ingroup oggetti

   \param L la lista di riferimento
   \param I l'indice
   \param T il tipo di elemento (vi, vd, mi, md, str)

	\note fare attenzione, perche´ il tipo non viene controllato!
*/
#define ACCEDIlst(L, I, T) L->dati[I - 1].T
#define _ACCEDIlst(L, I, T) L->dati[I - 1].T

//! Controlla (in debug) che la posizione di scrittura di un elemento di una lista sia corretto.
/*!
	\ingroup debug

   \param L la lista di riferimento
   \param I l'indice

	\note fare attenzione, percheï¿½ il tipo non viene controllato!
*/
#ifdef MDEBUG
	#define CtrlSlst(L, I) _controllaSLista(L, I - 1, #L, __FILE__, __LINE__)
	#define _CtrlSlst(L, I) _controllaSLista(L, I - 1, #L, nomefile, linea)
#else
	#define CtrlSlst(L, I)
	#define _CtrlSlst(L, I)
#endif
#ifdef MDEBUG
	void  _controllaSLista(const LISTA *const l, int e, const char *nome, const char *nomefile, int linea);
#else
	void  _controllaSLista(const LISTA *l, int e, const char *nome);
#endif

//! Assegna un elemento di una lista.
/*!
	\ingroup oggetti

   \param L la lista di riferimento
   \param I l'indice
   \param T il tipo di elemento (vi, vd, mi, md, str)
   \param Val il valore da assegnare

	\note fare attenzione, perche´ il tipo non viene controllato!
*/
#define ASSEGNAlst(L, I, T, Val) L->dati[I - 1].T = Val
#define _ASSEGNAlst(L, I, T, Val) L->dati[I - 1].T = Val

// privata, usata da LENGTHlst
#ifdef MDEBUG
	int  _lengthList(const LISTA *l, const char *nome, const char *nomefile, int linea);
#else
	int  _lengthList(const LISTA *l, const char *nome);
#endif

//! copia una LISTA
/*!
	\ingroup oplst

	\param L1 un puntatore per la nuova lista
	\param L2 la lista da cui copiare

	\return il puntatore alla nuova lista
*/
#ifdef MDEBUG
	#define CopiaLISTA(L1, L2) L1 = _CopiaLista(L1, L2, #L1, __FILE__, __LINE__)
	#define _CopiaLISTA(L1, L2) L1 = _CopiaLista(L1, L2, nome, nomefile, linea)
#else
	#define CopiaLISTA(L1, L2) L1 = _CopiaLista(L1, L2)
	#define _CopiaLISTA(L1, L2) L1 = _CopiaLista(L1, L2)
#endif
#ifdef MDEBUG
	LISTA *_CopiaLista(LISTA *ris, const LISTA *l2, const char *nome,const char *nomefile, int linea);
#else
	LISTA *_CopiaLista(LISTA *ris, const LISTA *l2);
#endif

//! Controlla se due liste sono uguali tenendo conto di approssimazioni nel caso di 'double'
/*!
	\ingroup oplst

	\param L1 la prima lista
	\param L2 la seconda lista

	\return 1 se uguali, 0 se diverse
*/
#ifdef MDEBUG
	#define UGUALIlst(L1, L2) L1 = _UgualiLst(L1, L2, #L1, __FILE__, __LINE__)
	#define _UGUALIlst(L1, L2) L1 = _UgualiLst(L1, L2, nome, nomefile, linea)
#else
	#define UGUALIlst(L1, L2) L1 = _UgualiLst(L1, L2)
	#define _UGUALIlst(L1, L2) L1 = _UgualiLst(L1, L2)
#endif
#ifdef MDEBUG
	bool _UgualiLst(const LISTA *l1, const LISTA *l2, const char *nome, const char *nomefile, int linea);
#else
	bool _UgualiLst(const LISTA *l1, const LISTA *l2);
#endif

//! Alloca una stringa e tiene traccia del numero di allocazioni (in debug)
/*!
	\ingroup oggetti

   \param N il nome della stringa
   \param V il valore iniziale
*/
#ifdef MDEBUG
	#define CREAstr(N, V) N = g_string_new(V); tot_str_alloc++;
	#define _CREAstr(N, V) N = g_string_new(V); tot_str_intalloc++;
#else
	#define CREAstr(N, V) N = g_string_new(V);
#endif

//! Disalloca una stringa e tiene traccia del numero di disallocazioni (in debug)
/*!
	\ingroup oggetti

   \param N il nome della stringa
*/
#ifdef MDEBUG
	#define CANCELLAstr(N) g_string_free(N, TRUE); N = NULL; tot_str_dealloc++;
	#define _CANCELLAstr(N) g_string_free(N, TRUE); N = NULL; tot_str_intdealloc++;
#else
	#define CANCELLAstr(N) g_string_free(N, TRUE); N = NULL;
	#define _CANCELLAstr(N) g_string_free(N, TRUE); N = NULL;
#endif

//! converte un puntatore ad una stringa in una espressione SEXP.
/*!
	\ingroup convstr

	\param str il puntatore da convertire
	\param nProtected un contatore da passare tra tutti gli oggetti costruiti nella stessa funzione esportata

	\return l'espressione convertita
*/
#ifdef MDEBUG
	#define daSTRINGA(str, nProtected) _daSTRINGA(str, nProtected, #str)
#else
	#define daSTRINGA(str, nProtected) _daSTRINGA(str, nProtected)
#endif
#ifdef MDEBUG
	SEXP _daSTRINGA(GString *str, int *nProtected, const char *nome);
#else
	SEXP _daSTRINGA(GString *str, int *nProtected);
#endif

//! predispone un puntatore a stringa da una espressione SEXP.
/*!
	\ingroup convstr

	\param s l'espressione SEXP da convertire
	\param nProtected un contatore da passare tra tutti gli oggetti costruiti nella stessa funzione esportata
	\param nome un nome per la stringa (ignorato in release)

	\return il puntatore alla stringa
*/
#ifdef MDEBUG
	#define inSTRINGA(s, nProtected, nome) _inSTRINGA(s, nProtected, nome)
#else
	#define inSTRINGA(s, nProtected, nome) _inSTRINGA(s, nProtected)
#endif
#ifdef MDEBUG
	GString *_inSTRINGA(SEXP s, int *nProtected, const char *nome);
#else
	GString *_inSTRINGA(SEXP s, int *nProtected);
#endif

//! moltiplica gli elementi corrispondenti di una matrice di tipo 'double' con una di tipo 'int'.
/*!
	\ingroup opmatr

	\param ris una matrice per il risultato
	\param m1 la prima matrice
	\param m2 la seconda matrice

	\test ris <- m1 * m2 # con m1 di double e m2 di interi
*/
#ifdef MDEBUG
	#define moltiplica_mm_di(ris, m1, m2) _moltiplica_mm_di(ris, m1, m2, #ris, __FILE__, __LINE__)
#else
	#define moltiplica_mm_di(ris, m1, m2) _moltiplica_mm_di(ris, m1, m2)
#endif
#ifdef MDEBUG
	MATRICEd * _moltiplica_mm_di(MATRICEd *ris, const MATRICEd *m1, const MATRICEi *m2, const char *nome, const char *nomefile, int linea);
#else
	MATRICEd * _moltiplica_mm_di(MATRICEd *ris, const MATRICEd *m1, const MATRICEi *m2);
#endif

//! moltiplica gli elementi corrispondenti di una matrice di tipo 'int' con una di tipo 'double'.
/*!
	\ingroup opmatr

	\param ris una matrice per il risultato
	\param m1 la prima matrice
	\param m2 la seconda matrice

	\test ris <- m1 * m2 # con m1 di interi e m2 di double
*/
#ifdef MDEBUG
	#define moltiplica_mm_id(ris, m1, m2) _moltiplica_mm_id(ris, m1, m2, #ris, __FILE__, __LINE__)
#else
	#define moltiplica_mm_id(ris, m1, m2) _moltiplica_mm_id(ris, m1, m2)
#endif
#ifdef MDEBUG
	MATRICEi * _moltiplica_mm_id(MATRICEi *ris, const MATRICEi *m1, const MATRICEd *m2, const char *nome, const char *nomefile, int linea);
#else
	MATRICEi * _moltiplica_mm_id(MATRICEi *ris, const MATRICEi *m1, const MATRICEd *m2);
#endif

#endif
