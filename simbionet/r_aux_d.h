/*! \file r_aux_d.h
	\brief Auxiliary functions for 'double'.
*/
#include "r_aux.h"

#ifndef R_AUX_d_H
#define R_AUX_d_H

#define MIN_MERGESORT_LIST_SIZE_d    32

//! Controlla il segno di un numero.
/*!
	\ingroup rc

	\param x il numero di riferimento

	\return il segno di x.
*/
#ifdef MDEBUG
	#define segno_s_d(x) _segno_s_d(x, __FILE__, __LINE__)
#else
	#define segno_s_d(x) _segno_s_d(x)
#endif
#ifdef MDEBUG
	int  _segno_s_d(double x, const char *nomefile, int linea);
#else
	int  _segno_s_d(double x);
#endif

//! Alloca un vettore
/*!
	\ingroup oggetti

   \param N il nome del vettore
   \param D il numero di elementi
*/
#ifdef MDEBUG
	#define CREAv_d(N, D) N = _creaVett_d(N, (D), #N, __FILE__, __LINE__)
	#define _CREAv_d(N, D) N = _creaVett_d(N, (D), nome, nomefile, linea)
#else
	#define CREAv_d(N, D) N = _creaVett_d(N, (D))
	#define _CREAv_d(N, D) N = _creaVett_d(N, (D))
#endif

//! Alloca una matrice
/*!
	\ingroup oggetti

   \param M il nome della matrice
   \param R il numero di righe
   \param C il numero di colonne
*/
#ifdef MDEBUG
	#define CREAm_d(M, R, C) M = _creaMatr_d(M, (R), (C), #M, __FILE__, __LINE__)
	#define _CREAm_d(M, R, C) M = _creaMatr_d(M, (R), (C), nome, nomefile, linea)
#else
	#define CREAm_d(M, R, C) M = _creaMatr_d(M, (R), (C))
	#define _CREAm_d(M, R, C) M = _creaMatr_d(M, (R), (C))
#endif

//! Disalloca un vettore
/*!
	\ingroup oggetti

   \param V il nome del vettore
*/
#ifdef MDEBUG
	#define CANCELLAv_d(V) V = _cancellaVett_d(V, #V, __FILE__, __LINE__)
	#define _CANCELLAv_d(V) V = _cancellaVett_d(V, nome, nomefile, linea)
#else
	#define CANCELLAv_d(V) V = _cancellaVett_d(V)
	#define _CANCELLAv_d(V) V = _cancellaVett_d(V)
#endif

//! Disalloca una matrice
/*!
	\ingroup oggetti

   \param M il nome della matrice
*/
#ifdef MDEBUG
	#define CANCELLAm_d(M) M = _cancellaMatr_d(M, #M, __FILE__, __LINE__)
	#define _CANCELLAm_d(M) M = _cancellaMatr_d(M, nome, nomefile, linea)
#else
	#define CANCELLAm_d(M) M = _cancellaMatr_d(M)
	#define _CANCELLAm_d(M) M = _cancellaMatr_d(M)
#endif

//! Dimensione di un vettore di tipo 'double'.
/*!
	\ingroup oggetti

   \param V il vettore di riferimento
*/
#ifdef MDEBUG
	#define LENGTHv_d(V) _lengthVett_d(V, #V, __FILE__, __LINE__)
#else
	#define LENGTHv_d(V) V->dim
#endif

//! Numero di righe di una matrice di tipo 'double'.
/*!
	\ingroup oggetti

   \param M la matrice di riferimento
*/
#ifdef MDEBUG
	#define LENGTHm1_d(M) _righeMatr_d(M, #M, __FILE__, __LINE__)
#else
	#define LENGTHm1_d(M) M->nr
#endif

//! Numero di colonne di una matrice di tipo 'double'.
/*!
	\ingroup oggetti

   \param M la matrice di riferimento
*/
#ifdef MDEBUG
	#define LENGTHm2_d(M) _colonneMatr_d(M, #M, __FILE__, __LINE__)
#else
	#define LENGTHm2_d(M) M->nc
#endif

//! Accede ad un elemento di un vettore di tipo 'double'.
/*!
   \param V il vettore di riferimento
   \param I l'indice
*/
#ifdef MDEBUG
	#define ACCEDIv_d(V, I) _accediVett_d(V, I - 1, #V, __FILE__, __LINE__)
	#define _ACCEDIv_d(V, I) _accediVett_d(V, I - 1, #V, nomefile, linea)
#else
	#define ACCEDIv_d(V, I) V->dati[I - 1]
	#define _ACCEDIv_d(V, I) V->dati[I - 1]
#endif

//! Accede ad un elemento di una matrice di tipo 'double'.
/*!
	\ingroup oggetti

   \param M la matrice di riferimento
   \param R l'indice della riga
   \param C l'indice della colonna
*/
#ifdef MDEBUG
	#define ACCEDIm_d(M, R, C) _accediMatr_d(M, R - 1, C - 1, #M, __FILE__, __LINE__)
	#define _ACCEDIm_d(M, R, C) _accediMatr_d(M, R - 1, C - 1, #M, nomefile, linea)
#else
	#define ACCEDIm_d(M, R, C) M->dati[(R - 1) + M->nr * (C - 1)]
	#define _ACCEDIm_d(M, R, C) M->dati[(R - 1) + M->nr * (C - 1)]
#endif

//! Accede ad un elemento di una matrice vettore di tipo 'double' come se fosse un vettore.
/*!
	\ingroup oggetti

   \param M la matrice di riferimento
   \param I l'indice
*/
#ifdef MDEBUG
	#define ACCEDImv_d(M, I) _accediMVett_d(M, I - 1, #M, __FILE__, __LINE__)
	#define _ACCEDImv_d(M, I) _accediMVett_d(M, I - 1, #M, nomefile, linea)
#else
	#define ACCEDImv_d(M, I) M->dati[I - 1]
	#define _ACCEDImv_d(M, I) M->dati[I - 1]
#endif

//! Assegna un elemento di un vettore di tipo 'double'.
/*!
	\ingroup oggetti

   \param V il vettore di riferimento
   \param I l'indice
   \param Val il valore da assegnare
*/
#ifdef MDEBUG
	#define ASSEGNAv_d(V, I, Val) _assegnaVett_d(V, I - 1, Val, #V, __FILE__, __LINE__)
	#define _ASSEGNAv_d(V, I, Val) _assegnaVett_d(V, I - 1, Val, #V, nomefile, linea)
#else
	#define ASSEGNAv_d(V, I, Val) V->dati[I - 1] = Val
	#define _ASSEGNAv_d(V, I, Val) V->dati[I - 1] = Val
#endif

//! Assegna un elemento di una matrice di tipo 'double'.
/*!
	\ingroup oggetti

   \param M la matrice di riferimento
   \param R l'indice della riga
   \param C l'indice della colonna
   \param V il valore da assegnare
*/
#ifdef MDEBUG
	#define ASSEGNAm_d(M, R, C, V) _assegnaMatr_d(M, R - 1, C - 1, V, #M, __FILE__, __LINE__)
	#define _ASSEGNAm_d(M, R, C, V) _assegnaMatr_d(M, R - 1, C - 1, V, #M, nomefile, linea)
#else
	#define ASSEGNAm_d(M, R, C, V) M->dati[(R - 1) + M->nr * (C - 1)] = V
	#define _ASSEGNAm_d(M, R, C, V) M->dati[(R - 1) + M->nr * (C - 1)] = V
#endif

//! Assegna un elemento di una matrice vettore di tipo 'double' indicizzandola come se fosse un vettore.
/*!
	\ingroup oggetti

   \param M la matrice di riferimento
   \param I l'indice
   \param V il valore da assegnare
*/
#ifdef MDEBUG
	#define ASSEGNAmv_d(M, I, V) _assegnaMVett_d(M, I - 1, V, #M, __FILE__, __LINE__)
	#define _ASSEGNAmv_d(M, I, V) _assegnaMVett_d(M, I - 1, V, #M, nomefile, linea)
#else
	#define ASSEGNAmv_d(M, I, V) M->dati[I - 1] = V
	#define _ASSEGNAmv_d(M, I, V) M->dati[I - 1] = V
#endif

// struttura interna per mergesort
typedef struct {
	double val;
	int pos;
} Elem_d;

//! Struttura per le matrici di tipo 'double'.
typedef struct  {
	int nr; /*!< numero di righe */
	int nc; /*!< numero di colonne */
	int alloc_r; /*!< uso interno */
	int alloc_c; /*!< uso interno */
	bool r; /*!< uso interno */
	double *dati; /*!< elementi della matrice */
	GList *mem; /*!< uso interno */
} MATRICEd;

GList *allocVett_d;
GList *allocMatr_d;

FILE *fp_mdbg_d;
FILE *fpm_v_d;
FILE *fpm_m_d;

// funzione interna per ordine_d
void _mergesort_array_d(Elem_d a[], int size, Elem_d temp[], bool decr);

// interna
void _stampa_lista_d(const char *pref, GList *lista);

//! Controlla se due numeri sono uguali tenendo conto di approssimazioni nel caso di 'double'
/*!
	\ingroup opnum

   \param X il primo numero
   \param Y il secondo numero

	\return 1 se uguali, 0 se diversi
*/
#ifdef VERSIONE_i
	#define UGUALE(X, Y) (X == Y)
#else
	#define UGUALE(X, Y) Uguale(X, Y)
#endif

//! Controlla se due numeri sono uguali diversi conto di approssimazioni nel caso di 'double'
/*!
	\ingroup opnum

   \param X il primo numero
   \param Y il secondo numero

	\return 1 se diversi, 0 se uguali
*/
#ifdef VERSIONE_i
	#define DIVERSO(X, Y) (X != Y)
#else
	#define DIVERSO(X, Y) !Uguale(X, Y)
#endif

//! Controlla se due vettori sono uguali tenendo conto di approssimazioni nel caso di 'double'
/*!
	\ingroup opvett

   \param v1 il primo vettore
   \param v2 il secondo vettore

	\return 1 se uguali, 0 se diversi
*/
#ifdef MDEBUG
	#define Uguali_v_d(v1, v2) _Uguali_v_d(v1, v2, __FILE__, __LINE__)
#else
	#define Uguali_v_d(v1, v2) _Uguali_v_d(v1, v2)
#endif
#ifdef MDEBUG
	bool  _Uguali_v_d(const VETTOREd *v1, const VETTOREd *v2, const char *nomefile, int linea);
#else
	bool  _Uguali_v_d(const VETTOREd *v1, const VETTOREd *v2);
#endif

//! Controlla se due matrici sono uguali tenendo conto di approssimazioni nel caso di 'double'
/*!
   \param m1 la prima matrice
   \param m2 la seconda matrice

	\return 1 se uguali, 0 se diverse
*/
#ifdef MDEBUG
	#define Uguali_m_d(m1, m2) _Uguali_m_d(m1, m2, __FILE__, __LINE__)
#else
	#define Uguali_m_d(m1, m2) _Uguali_m_d(m1, m2)
#endif
#ifdef MDEBUG
	bool  _Uguali_m_d(const MATRICEd *m1, const MATRICEd *m2, const char *nomefile, int linea);
#else
	bool  _Uguali_m_d(const MATRICEd *m1, const MATRICEd *m2);
#endif

// interna
void _InitDbg_d(bool std_out);

//! Stampa un vettore come farebbe la funzione "write" di R
/*!
	\ingroup oggetti

	\param v il vettore da stampare
*/
void _StampaRawVett_d(const VETTOREd *v);

//! Stampa un vettore.
/*!
	\ingroup debug

	\param v il vettore da stampare
*/
void _StampaVett_d(const VETTOREd *v);

//! Stampa una matrice come farebbe la funzione "write" di R.
/*!
	\ingroup debug

	\param m la matrice da stampare
*/
void _StampaMatr_d(const MATRICEd *m);

//! Stampa una matrice.
/*!
	\ingroup debug

	\param m la matrice da stampare
*/
void _StampaRawMatr_d(const MATRICEd *m);

//! Inizializza un vettore.
/*!
	\ingroup opvett

	\param v il vettore da inizializzare
	\param val il valore da assegnare

*/
#ifdef MDEBUG
	#define InitVett_d(v, val) _InitVett_d(v, val, __FILE__, __LINE__)
#else
	#define InitVett_d(v, val) _InitVett_d(v, val)
#endif
#ifdef MDEBUG
	void  _InitVett_d(VETTOREd *v, double val, const char *nomefile, int linea);
#else
	void  _InitVett_d(VETTOREd *v, double val);
#endif

//! Inizializza una matrice.
/*!
	\ingroup opmatr

	\param m la matrice da inizializzare
	\param val il valore da assegnare

*/
#ifdef MDEBUG
	#define InitMatr_d(m, val) _InitMatr_d(m, val, __FILE__, __LINE__)
#else
	#define InitMatr_d(m, val) _InitMatr_d(m, val)
#endif
#ifdef MDEBUG
	void  _InitMatr_d(MATRICEd *m, double val, const char *nomefile, int linea);
#else
	void  _InitMatr_d(MATRICEd *m, double val);
#endif

//! Scrive il buffer di debug
/*!
	\ingroup debug
*/
void _Fflush_d();

//! Calcola il massimo.
/*!
	\ingroup opnum

	\param a primo valore
	\param b secondo valore

	\return il massimo tra i due valori
*/
#ifdef MDEBUG
	#define max_s_d(a, b) _max_s_d(a, b, __FILE__, __LINE__)
#else
	#define max_s_d(a, b) _max_s_d(a, b)
#endif
#ifdef MDEBUG
	double  _max_s_d(double a, double b, const char *nomefile, int linea);
#else
	double  _max_s_d(double a, double b);
#endif

//! Calcola il minimo.
/*!
	\ingroup opnum

	\param a primo valore
	\param b secondo valore

	\return il minimo tra i due valori
*/
#ifdef MDEBUG
	#define min_s_d(a, b) _min_s_d(a, b, __FILE__, __LINE__)
#else
	#define min_s_d(a, b) _min_s_d(a, b)
#endif
#ifdef MDEBUG
	double  _min_s_d(double a, double b, const char *nomefile, int linea);
#else
	double  _min_s_d(double a, double b);
#endif

#ifdef MDEBUG
	#define creaVett_d(ris, dim) _creaVett_d(ris, dim, #ris, __FILE__, __LINE__)
#else
	#define creaVett_d(ris, dim) _creaVett_d(ris, dim)
#endif
#ifdef MDEBUG
	VETTOREd * _creaVett_d(VETTOREd *ris, int dim, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _creaVett_d(VETTOREd *ris, int dim);
#endif

#ifdef MDEBUG
	#define creaMatr_d(ris, nr, nc) _creaMatr_d(ris, nr, nc, #ris, __FILE__, __LINE__)
#else
	#define creaMatr_d(ris, nr, nc) _creaMatr_d(ris, nr, nc)
#endif
#ifdef MDEBUG
	MATRICEd * _creaMatr_d(MATRICEd *ris, int nr, int nc, const char *nome, const char *nomefile, int linea);
#else
	MATRICEd * _creaMatr_d(MATRICEd *ris, int nr, int nc);
#endif

// privata, usata da LENGTHv_d
#ifdef MDEBUG
	#define lengthVett_d(v, nome) _lengthVett_d(v, nome, __FILE__, __LINE__)
#else
	#define lengthVett_d(v, nome) _lengthVett_d(v, nome)
#endif
#ifdef MDEBUG
	int  _lengthVett_d(const VETTOREd *v, const char *nome, const char *nomefile, int linea);
#else
	int  _lengthVett_d(const VETTOREd *v, const char *nome);
#endif
// privata, usata da LENGTHm1_d
#ifdef MDEBUG
	#define righeMatr_d(m, nome) _righeMatr_d(m, nome, __FILE__, __LINE__)
#else
	#define righeMatr_d(m, nome) _righeMatr_d(m, nome)
#endif
#ifdef MDEBUG
	int  _righeMatr_d(const MATRICEd *m, const char *nome, const char *nomefile, int linea);
#else
	int  _righeMatr_d(const MATRICEd *m, const char *nome);
#endif
// privata, usata da LENGTHm2_d
#ifdef MDEBUG
	#define colonneMatr_d(m, nome) _colonneMatr_d(m, nome, __FILE__, __LINE__)
#else
	#define colonneMatr_d(m, nome) _colonneMatr_d(m, nome)
#endif
#ifdef MDEBUG
	int  _colonneMatr_d(const MATRICEd *m, const char *nome, const char *nomefile, int linea);
#else
	int  _colonneMatr_d(const MATRICEd *m, const char *nome);
#endif

// privata, usata da ACCEDIv_d
#ifdef MDEBUG
	#define accediVett_d(v, indx, nome) _accediVett_d(v, indx, nome, __FILE__, __LINE__)
#else
	#define accediVett_d(v, indx, nome) _accediVett_d(v, indx, nome)
#endif
#ifdef MDEBUG
	double  _accediVett_d(const VETTOREd *v, int indx, const char *nome, const char *nomefile, int linea);
#else
	double  _accediVett_d(const VETTOREd *v, int indx, const char *nome);
#endif
// privata, usata da ACCEDIm_d
#ifdef MDEBUG
	#define accediMatr_d(m, r, c, nome) _accediMatr_d(m, r, c, nome, __FILE__, __LINE__)
#else
	#define accediMatr_d(m, r, c, nome) _accediMatr_d(m, r, c, nome)
#endif
#ifdef MDEBUG
	double  _accediMatr_d(const MATRICEd *m, int r, int c, const char *nome, const char *nomefile, int linea);
#else
	double  _accediMatr_d(const MATRICEd *m, int r, int c, const char *nome);
#endif
// privata, usata da ACCEDIvm_d
#ifdef MDEBUG
	#define accediMVett_d(m, indx, nome) _accediMVett_d(m, indx, nome, __FILE__, __LINE__)
#else
	#define accediMVett_d(m, indx, nome) _accediMVett_d(m, indx, nome)
#endif
#ifdef MDEBUG
	double  _accediMVett_d(const MATRICEd *m, int indx, const char *nome, const char *nomefile, int linea);
#else
	double  _accediMVett_d(const MATRICEd *m, int indx, const char *nome);
#endif

// privata, usata da ASSEGNAv_d
#ifdef MDEBUG
	#define assegnaVett_d(v, indx, val, nome) _assegnaVett_d(v, indx, val, nome, __FILE__, __LINE__)
#else
	#define assegnaVett_d(v, indx, val, nome) _assegnaVett_d(v, indx, val, nome)
#endif
#ifdef MDEBUG
	void  _assegnaVett_d(VETTOREd *v, int indx, double val, const char *nome, const char *nomefile, int linea);
#else
	void  _assegnaVett_d(VETTOREd *v, int indx, double val, const char *nome);
#endif
// privata, usata da ASSEGNAm_d
#ifdef MDEBUG
	#define assegnaMatr_d(m, r, c, val, nome) _assegnaMatr_d(m, r, c, val, nome, __FILE__, __LINE__)
#else
	#define assegnaMatr_d(m, r, c, val, nome) _assegnaMatr_d(m, r, c, val, nome)
#endif
#ifdef MDEBUG
	void  _assegnaMatr_d(MATRICEd *m, int r, int c, double val, const char *nome, const char *nomefile, int linea);
#else
	void  _assegnaMatr_d(MATRICEd *m, int r, int c, double val, const char *nome);
#endif
// privata, usata da ASSEGNAvm_d
#ifdef MDEBUG
	#define assegnaMVett_d(m, indx, val, nome) _assegnaMVett_d(m, indx, val, nome, __FILE__, __LINE__)
#else
	#define assegnaMVett_d(m, indx, val, nome) _assegnaMVett_d(m, indx, val, nome)
#endif
#ifdef MDEBUG
	void  _assegnaMVett_d(MATRICEd *m, int indx, double val, const char *nome, const char *nomefile, int linea);
#else
	void  _assegnaMVett_d(MATRICEd *m, int indx, double val, const char *nome);
#endif

// privata, usata da CANCELLAv_d
#ifdef MDEBUG
	#define cancellaVett_d(ris) _cancellaVett_d(ris, #ris, __FILE__, __LINE__)
#else
	#define cancellaVett_d(ris) _cancellaVett_d(ris)
#endif
#ifdef MDEBUG
	VETTOREd * _cancellaVett_d(VETTOREd *ris, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _cancellaVett_d(VETTOREd *ris);
#endif
// privata, usata da CANCELLAm_d
#ifdef MDEBUG
	#define cancellaMatr_d(ris) _cancellaMatr_d(ris, #ris, __FILE__, __LINE__)
#else
	#define cancellaMatr_d(ris) _cancellaMatr_d(ris)
#endif
#ifdef MDEBUG
	MATRICEd * _cancellaMatr_d(MATRICEd *ris, const char *nome, const char *nomefile, int linea);
#else
	MATRICEd * _cancellaMatr_d(MATRICEd *ris);
#endif

//! Stampa le informazioni su un vettore di tipo 'double'.
/*!
	\ingroup debug

	\param v il vettore di interesse
*/
void _infoVett_d(const VETTOREd *v);

//! Stampa le informazioni su una matrice di tipo 'double'.
/*!
	\ingroup debug

	\param m la matrice di interesse
*/
void _infoMatr_d(const MATRICEd *m);

#ifdef MDEBUG
	#define controllaCanc_d() _controllaCanc_d(__FILE__, __LINE__)
#else
	#define controllaCanc_d() _controllaCanc_d()
#endif
#ifdef MDEBUG
	void  _controllaCanc_d(const char *nomefile, int linea);
#else
	void  _controllaCanc_d();
#endif

//! converte un VETTOREd in una espressione SEXP.
/*!
	\ingroup convvett

	\param v il vettore da convertire
	\param nProtected un contatore da passare tra tutti gli oggetti costruiti nella stessa funzione esportata

	\return l'espressione convertita
*/
#ifdef MDEBUG
	#define daVETTORE_d(v, nProtected) _daVETTORE_d(v, nProtected, __FILE__, __LINE__)
#else
	#define daVETTORE_d(v, nProtected) _daVETTORE_d(v, nProtected)
#endif
#ifdef MDEBUG
	SEXP  _daVETTORE_d(VETTOREd *v, int *nProtected, const char *nomefile, int linea);
#else
	SEXP  _daVETTORE_d(VETTOREd *v, int *nProtected);
#endif

//! converte una MATRICEd in una espressione SEXP.
/*!
	\ingroup convmatr

	\param m la matrice da convertire
	\param nProtected un contatore da passare tra tutti gli oggetti costruiti nella stessa funzione esportata

	\return l'espressione convertita
*/
#ifdef MDEBUG
	#define daMATRICE_d(m, nProtected) _daMATRICE_d(m, nProtected, __FILE__, __LINE__)
#else
	#define daMATRICE_d(m, nProtected) _daMATRICE_d(m, nProtected)
#endif
#ifdef MDEBUG
	SEXP  _daMATRICE_d(MATRICEd *m, int *nProtected, const char *nomefile, int linea);
#else
	SEXP  _daMATRICE_d(MATRICEd *m, int *nProtected);
#endif

//! predispone un VETTOREd da una espressione SEXP.
/*!
	\ingroup convvett

	\param ris l'espressione SEXP da convertire
	\param nProtected un contatore da passare tra tutti gli oggetti costruiti nella stessa funzione esportata

	\return il vettore VETTOREd da restituire
*/
#ifdef MDEBUG
	#define inVETTORE_d(ris, nProtected) _inVETTORE_d(ris, nProtected, #ris, __FILE__, __LINE__)
#else
	#define inVETTORE_d(ris, nProtected) _inVETTORE_d(ris, nProtected)
#endif
#ifdef MDEBUG
	VETTOREd * _inVETTORE_d(SEXP ris, int *nProtected, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _inVETTORE_d(SEXP ris, int *nProtected);
#endif

//! predispone una MATRICEd da una espressione SEXP.
/*!
	\ingroup convmatr

	\param ris l'espressione SEXP da convertire
	\param nProtected un contatore da passare tra tutti gli oggetti costruiti nella stessa funzione esportata

	\return la matrice MATRICEd da restituire
*/
#ifdef MDEBUG
	#define inMATRICE_d(ris, nProtected) _inMATRICE_d(ris, nProtected, #ris, __FILE__, __LINE__)
#else
	#define inMATRICE_d(ris, nProtected) _inMATRICE_d(ris, nProtected)
#endif
#ifdef MDEBUG
	MATRICEd * _inMATRICE_d(SEXP ris, int *nProtected, const char *nome, const char *nomefile, int linea);
#else
	MATRICEd * _inMATRICE_d(SEXP ris, int *nProtected);
#endif

//! restituisce un vettore di tipo 'double' corrispondente agli indici degli elementi diversi da \c val nella riga \c r della matrice m.
/*!
	\ingroup wmatr

	\param ris un vettore di tipo 'int' per il risultato
	\param m la matrice di riferimento
	\param r la riga di riferimento
	\param val il valore di riferimento

	\test ris <- which(m[r, ] != val)
*/
#ifdef MDEBUG
	#define which_m_rowindxne_d(ris, m, r, val) _which_m_rowindxne_d(ris, m, r, val, #ris, __FILE__, __LINE__)
#else
	#define which_m_rowindxne_d(ris, m, r, val) _which_m_rowindxne_d(ris, m, r, val)
#endif
#ifdef MDEBUG
	VETTOREi * _which_m_rowindxne_d(VETTOREi *ris, const MATRICEd *m, int r, double val, const char *nome, const char *nomefile, int linea);
#else
	VETTOREi * _which_m_rowindxne_d(VETTOREi *ris, const MATRICEd *m, int r, double val);
#endif

//! copia un vettore di tipo 'double' in un altro vettore dello stesso tipo.
/*!
	\ingroup segmenti

	\param ris il vettore di destinazione
	\param da il vettore sorgente
	\param st l'indice di partenza
	\param end l'indice di arrivo

	\return il vettore copiato

	\test ris <- da[st:end]
*/
#ifdef MDEBUG
	#define copia_v_d(ris, da, st, end) _copia_v_d(ris, da, st, end, #ris, __FILE__, __LINE__)
#else
	#define copia_v_d(ris, da, st, end) _copia_v_d(ris, da, st, end)
#endif
#ifdef MDEBUG
	VETTOREd * _copia_v_d(VETTOREd *ris, const VETTOREd *da, int st, int end, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _copia_v_d(VETTOREd *ris, const VETTOREd *da, int st, int end);
#endif

//! divide una sequenza di tipo 'double' per un numero.
/*!
	\ingroup segmenti

	\param ris un vettore per il risultato
	\param da il numero da cui partire
	\param a il numero a cui arrivare (incremento di 1)
	\param div il numero per cui dividere (diverso da zero)

	\test ris <- (da:a) / div
*/
#ifdef MDEBUG
	#define op_ss_seqdiv_d(ris, da, a, div) _op_ss_seqdiv_d(ris, da, a, div, #ris, __FILE__, __LINE__)
#else
	#define op_ss_seqdiv_d(ris, da, a, div) _op_ss_seqdiv_d(ris, da, a, div)
#endif
#ifdef MDEBUG
	VETTOREd * _op_ss_seqdiv_d(VETTOREd *ris, int da, int a, double div, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _op_ss_seqdiv_d(VETTOREd *ris, int da, int a, double div);
#endif

//! complementa un vettore di tipo 'double'.
/*!
	\ingroup opvett

	\param ris un vettore per il risultato
	\param v il vettore da complementare

	\test ris <- 1 - v
*/
#ifdef MDEBUG
	#define complementa_d(ris, v) _complementa_d(ris, v, #ris, __FILE__, __LINE__)
#else
	#define complementa_d(ris, v) _complementa_d(ris, v)
#endif
#ifdef MDEBUG
	VETTOREd * _complementa_d(VETTOREd *ris, const VETTOREd *v, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _complementa_d(VETTOREd *ris, const VETTOREd *v);
#endif

//! crea un vettore ripetendo un vettore di tipo 'double'.
/*!
	\ingroup opvett

	\param ris un vettore per il risultato
	\param v il vettore da ripetere
	\param ripetizioni il numero di ripetizioni

	\test ris <- rep(v, ripetizioni)
*/
#ifdef MDEBUG
	#define rep_v_d(ris, v, ripetizioni) _rep_v_d(ris, v, ripetizioni, #ris, __FILE__, __LINE__)
#else
	#define rep_v_d(ris, v, ripetizioni) _rep_v_d(ris, v, ripetizioni)
#endif
#ifdef MDEBUG
	VETTOREd * _rep_v_d(VETTOREd *ris, const VETTOREd *v, int ripetizioni, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _rep_v_d(VETTOREd *ris, const VETTOREd *v, int ripetizioni);
#endif

//! crea un vettore ripetendo un numero di tipo 'double'.
/*!
	\ingroup opvett

	\param ris un vettore per il risultato
	\param num il numero da ripetere
	\param ripetizioni il numero di ripetizioni

	\test ris <- rep(num, ripetizioni)
*/
#ifdef MDEBUG
	#define rep_s_d(ris, num, ripetizioni) _rep_s_d(ris, num, ripetizioni, #ris, __FILE__, __LINE__)
#else
	#define rep_s_d(ris, num, ripetizioni) _rep_s_d(ris, num, ripetizioni)
#endif
#ifdef MDEBUG
	VETTOREd * _rep_s_d(VETTOREd *ris, double num, int ripetizioni, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _rep_s_d(VETTOREd *ris, double num, int ripetizioni);
#endif

//! crea una sequenza di tipo 'double'.
/*!
	\ingroup segmenti

	\param ris un vettore di tipo 'int' per il risultato
	\param da il valore di partenza
	\param a il valore d'arrivo
	\param incremento la quantita' da addizionare

	\test ris <- seq(da, a, incremento)
*/
#ifdef MDEBUG
	#define seq_d(ris, da, a, incremento) _seq_d(ris, da, a, incremento, #ris, __FILE__, __LINE__)
#else
	#define seq_d(ris, da, a, incremento) _seq_d(ris, da, a, incremento)
#endif
#ifdef MDEBUG
	VETTOREd * _seq_d(VETTOREd *ris, double da, double a, double incremento, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _seq_d(VETTOREd *ris, double da, double a, double incremento);
#endif

//! crea un vettore di tipo 'double' da due elementi dello stesso tipo.
/*!
	\ingroup oggetti

	\param ris un vettore per il risultato
	\param el1 il primo elemento
	\param el2 il secondo elemento

	\test ris <- c(el1, el2)
*/
#ifdef MDEBUG
	#define vettore2s_d(ris, el1, el2) _vettore2s_d(ris, el1, el2, #ris, __FILE__, __LINE__)
#else
	#define vettore2s_d(ris, el1, el2) _vettore2s_d(ris, el1, el2)
#endif
#ifdef MDEBUG
	VETTOREd * _vettore2s_d(VETTOREd *ris, double el1, double el2, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _vettore2s_d(VETTOREd *ris, double el1, double el2);
#endif

//! crea un vettore di tipo 'double' da tre elementi dello stesso tipo.
/*!
	\ingroup oggetti

	\param ris un vettore per il risultato
	\param el1 il primo elemento
	\param el2 il secondo elemento
	\param el3 terzo elemento

	\test ris <- c(el1, el2, el3)
*/
#ifdef MDEBUG
	#define vettore3s_d(ris, el1, el2, el3) _vettore3s_d(ris, el1, el2, el3, #ris, __FILE__, __LINE__)
#else
	#define vettore3s_d(ris, el1, el2, el3) _vettore3s_d(ris, el1, el2, el3)
#endif
#ifdef MDEBUG
	VETTOREd * _vettore3s_d(VETTOREd *ris, double el1, double el2, double el3, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _vettore3s_d(VETTOREd *ris, double el1, double el2, double el3);
#endif

//! crea un vettore di tipo 'double' da due vettori dello stesso tipo.
/*!
	\ingroup oggetti

	\param ris un vettore  per il risultato
	\param v1 il primo vettore
	\param v2 il secondo vettore

	\test ris <- c(v1, v2)
*/
#ifdef MDEBUG
	#define vettore2v_d(ris, v1, v2) _vettore2v_d(ris, v1, v2, #ris, __FILE__, __LINE__)
#else
	#define vettore2v_d(ris, v1, v2) _vettore2v_d(ris, v1, v2)
#endif
#ifdef MDEBUG
	VETTOREd * _vettore2v_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _vettore2v_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2);
#endif

//! crea un vettore di tipo 'double' da tre vettori dello stesso tipo.
/*!
	\ingroup oggetti

	\param ris un vettore per il risultato
	\param v1 il primo vettore
	\param v2 il secondo vettore
	\param v3 il terzo vettore

	\test ris <- c(v1, v2, v3)
*/
#ifdef MDEBUG
	#define vettore3v_d(ris, v1, v2, v3) _vettore3v_d(ris, v1, v2, v3, #ris, __FILE__, __LINE__)
#else
	#define vettore3v_d(ris, v1, v2, v3) _vettore3v_d(ris, v1, v2, v3)
#endif
#ifdef MDEBUG
	VETTOREd * _vettore3v_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2, const VETTOREd *v3, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _vettore3v_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2, const VETTOREd *v3);
#endif

//! trova il massimo in un vettore di tipo 'double'.
/*!
	\ingroup opvett

	\param v il vettore in cui cercare

	\return l'elemento massimo nel vettore

	\test ris <- max(v)
*/
#ifdef MDEBUG
	#define max_v_d(v) _max_v_d(v, __FILE__, __LINE__)
#else
	#define max_v_d(v) _max_v_d(v)
#endif
#ifdef MDEBUG
	double  _max_v_d(const VETTOREd *v, const char *nomefile, int linea);
#else
	double  _max_v_d(const VETTOREd *v);
#endif

//! accoda un vettore di tipo 'double' ad un altro dello stesso tipo.
/*!
	\ingroup opvett

	\param ris il vettore di riferimento
	\param v2 il vettore da accodare

	\test ris <- c(ris, v2)
*/
#ifdef MDEBUG
	#define accoda1_vv_d(ris, v2) _accoda1_vv_d(ris, v2, #ris, __FILE__, __LINE__)
#else
	#define accoda1_vv_d(ris, v2) _accoda1_vv_d(ris, v2)
#endif
#ifdef MDEBUG
	VETTOREd * _accoda1_vv_d(VETTOREd *ris, const VETTOREd *v2, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _accoda1_vv_d(VETTOREd *ris, const VETTOREd *v2);
#endif

//! somma uno scalare di tipo 'double' ad un vettore di tipo 'double'.
/*!
	\ingroup opvett

	\param ris un vettore per il risultato
	\param v il vettore di riferimento
	\param s il valore da addizionare

	\test ris <- v + s
*/
#ifdef MDEBUG
	#define somma_vs_d(ris, v, s) _somma_vs_d(ris, v, s, #ris, __FILE__, __LINE__)
#else
	#define somma_vs_d(ris, v, s) _somma_vs_d(ris, v, s)
#endif
#ifdef MDEBUG
	VETTOREd * _somma_vs_d(VETTOREd *ris, const VETTOREd *v, double s, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _somma_vs_d(VETTOREd *ris, const VETTOREd *v, double s);
#endif

//! cerca un elemento in un vettore di tipo 'double'.
/*!
	\ingroup opvett

	\param el l'elemento da cercare
	\param v il vettore in cui cercare

	\return l'indice dell'elemento cercato (da uno) o zero se l'elemento non esiste

	\test el %IN% v
*/
#ifdef MDEBUG
	#define esiste_v_d(el, v) _esiste_v_d(el, v, __FILE__, __LINE__)
#else
	#define esiste_v_d(el, v) _esiste_v_d(el, v)
#endif
#ifdef MDEBUG
	int  _esiste_v_d(double el, const VETTOREd *v, const char *nomefile, int linea);
#else
	int  _esiste_v_d(double el, const VETTOREd *v);
#endif

//! elimina un elemento da un vettore di tipo 'double' tramite un indice.
/*!
	\ingroup opvett

	\param v il vettore che contiene l'elemento da eliminare (non e' costante perche' si potrebbe restringere)
	\param indx l'indice

	\test v <- v[-indx]
*/
#ifdef MDEBUG
	#define elimina1_indx_d(v, indx) _elimina1_indx_d(v, indx, __FILE__, __LINE__)
#else
	#define elimina1_indx_d(v, indx) _elimina1_indx_d(v, indx)
#endif
#ifdef MDEBUG
	void  _elimina1_indx_d(VETTOREd *v, int indx, const char *nomefile, int linea);
#else
	void  _elimina1_indx_d(VETTOREd *v, int indx);
#endif

//! assegna un un elemento di tipo 'double' ad un vettore di tipo 'double' in base ad un indice e si assicura che ci stia.
/*!
	\ingroup setvett

	\param ris il vettore di riferimento (non e' costante perche' potrebbe essere ricreato)
	\param indx l'indice da assegnare
	\param val il valore da assegnare

	\note modifica le dimensioni e da' un warning se il vettore viene ingrandito, ai valori intermedi viene assegnato il valore speciale NA

	\test ris[indx] <- val
*/
#ifdef MDEBUG
	#define assegna_v_d(ris, indx, val) _assegna_v_d(ris, indx, val, #ris, __FILE__, __LINE__)
#else
	#define assegna_v_d(ris, indx, val) _assegna_v_d(ris, indx, val)
#endif
#ifdef MDEBUG
	VETTOREd * _assegna_v_d(VETTOREd *ris, int indx, double val, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _assegna_v_d(VETTOREd *ris, int indx, double val);
#endif

//! somma le righe di una matrice di tipo 'double'.
/*!
	\ingroup opmatr

	\param ris un vettore per il risultato
	\param m la matrice di riferimento

	\test  ris <- apply(m, 1, sum)
*/
#ifdef MDEBUG
	#define somma_righe_d(ris, m) _somma_righe_d(ris, m, #ris, __FILE__, __LINE__)
#else
	#define somma_righe_d(ris, m) _somma_righe_d(ris, m)
#endif
#ifdef MDEBUG
	VETTOREd * _somma_righe_d(VETTOREd *ris, const MATRICEd *m, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _somma_righe_d(VETTOREd *ris, const MATRICEd *m);
#endif

//! rovescia un vettore di tipo 'double' "sul posto".
/*!
	\ingroup opvett

	\param v il vettore da rovesciare

	\test v <- rev(v)
*/
#ifdef MDEBUG
	#define rev1_d(v) _rev1_d(v, __FILE__, __LINE__)
#else
	#define rev1_d(v) _rev1_d(v)
#endif
#ifdef MDEBUG
	void  _rev1_d(VETTOREd *v, const char *nomefile, int linea);
#else
	void  _rev1_d(VETTOREd *v);
#endif

//! ordina gli elementi di un vettore di tipo 'double'.
/*!
	\ingroup opvett

	\param ris un vettore per il risultato
	\param v il vettore da ordinare
	\param decr se vero ordina in maniera decrescente

	\test ris <- sort(v, decr)
*/
#ifdef MDEBUG
	#define ordina_d(ris, v, decr) _ordina_d(ris, v, decr, #ris, __FILE__, __LINE__)
#else
	#define ordina_d(ris, v, decr) _ordina_d(ris, v, decr)
#endif
#ifdef MDEBUG
	VETTOREd * _ordina_d(VETTOREd *ris, const VETTOREd *v, bool decr, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _ordina_d(VETTOREd *ris, const VETTOREd *v, bool decr);
#endif

//! ordina gli elementi di un vettore di tipo 'double' "sul posto".
/*!
	\ingroup opvett

	\param v il vettore da ordinare
	\param decr se vero ordina in maniera decrescente

	\test v <- sort(v, decr)
*/
#ifdef MDEBUG
	#define ordina1_d(v, decr) _ordina1_d(v, decr, __FILE__, __LINE__)
#else
	#define ordina1_d(v, decr) _ordina1_d(v, decr)
#endif
#ifdef MDEBUG
	void  _ordina1_d(VETTOREd *v, bool decr, const char *nomefile, int linea);
#else
	void  _ordina1_d(VETTOREd *v, bool decr);
#endif

//! estrae una riga da una matrice di tipo 'double'.
/*!
	\ingroup opmatr

	\param ris un vettore per il risultato
	\param m la matrice di riferimento
	\param riga la riga di riferimento

	\test ris <- M[riga,]
*/
#ifdef MDEBUG
	#define riga_d(ris, m, riga) _riga_d(ris, m, riga, #ris, __FILE__, __LINE__)
#else
	#define riga_d(ris, m, riga) _riga_d(ris, m, riga)
#endif
#ifdef MDEBUG
	VETTOREd * _riga_d(VETTOREd *ris, const MATRICEd *m, int riga, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _riga_d(VETTOREd *ris, const MATRICEd *m, int riga);
#endif

//! restituisce un vettore corrispondente agli indici degli elementi minori o uguali di un determinato valore in un vettore di tipo 'double'.
/*!
	\ingroup wvett

	\param ris un vettore di tipo 'int' per il risultato
	\param v il vettore in cui cercare
	\param val il valore di riferimento

	\test ris <- which(v <= val)
*/
#ifdef MDEBUG
	#define which_v_indxle_d(ris, v, val) _which_v_indxle_d(ris, v, val, #ris, __FILE__, __LINE__)
#else
	#define which_v_indxle_d(ris, v, val) _which_v_indxle_d(ris, v, val)
#endif
#ifdef MDEBUG
	VETTOREi * _which_v_indxle_d(VETTOREi *ris, const VETTOREd *v, double val, const char *nome, const char *nomefile, int linea);
#else
	VETTOREi * _which_v_indxle_d(VETTOREi *ris, const VETTOREd *v, double val);
#endif

//! calcola il minimo in un vettore di tipo 'double'.
/*!
	\ingroup opvett

	\param v il vettore di riferimento

	\return il minimo del vettore

	\test ris <- min(v)
*/
#ifdef MDEBUG
	#define min_v_d(v) _min_v_d(v, __FILE__, __LINE__)
#else
	#define min_v_d(v) _min_v_d(v)
#endif
#ifdef MDEBUG
	double  _min_v_d(const VETTOREd *v, const char *nomefile, int linea);
#else
	double  _min_v_d(const VETTOREd *v);
#endif

//! somma uno scalare ad un vettore di tipo 'double' "sul posto".
/*!
	\ingroup opvett

	\param v il vettore di riferimento
	\param s il valore da addizionare

	\test v <- v + s
*/
#ifdef MDEBUG
	#define somma1_vs_d(v, s) _somma1_vs_d(v, s, __FILE__, __LINE__)
#else
	#define somma1_vs_d(v, s) _somma1_vs_d(v, s)
#endif
#ifdef MDEBUG
	void  _somma1_vs_d(VETTOREd *v, double s, const char *nomefile, int linea);
#else
	void  _somma1_vs_d(VETTOREd *v, double s);
#endif

//! crea un vettore di tipo 'double' a partire dagli elementi di un vettore di tipo 'double' corrispondenti ad un vettore di indici.
/*!
	\ingroup wvett

	\param ris un vettore per il risultato
	\param v il vettore da cui trarre i valori
	\param indx il vettore degli indici

	\test  v1 <- v2[indx]
*/
#ifdef MDEBUG
	#define copia_v_indx_d(ris, v, indx) _copia_v_indx_d(ris, v, indx, #ris, __FILE__, __LINE__)
#else
	#define copia_v_indx_d(ris, v, indx) _copia_v_indx_d(ris, v, indx)
#endif
#ifdef MDEBUG
	VETTOREd * _copia_v_indx_d(VETTOREd *ris, const VETTOREd *v, const VETTOREi *indx, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _copia_v_indx_d(VETTOREd *ris, const VETTOREd *v, const VETTOREi *indx);
#endif

//! restituisce un vettore corrispondente agli indici degli elementi maggiori di un determinato valore in un vettore di tipo 'double'.
/*!
	\ingroup wvett

	\param ris un vettore di tipo 'int' per il risultato
	\param v il vettore in cui cercare
	\param val il valore di riferimento

	\test ris <- which(v > val)
*/
#ifdef MDEBUG
	#define which_v_indxgt_d(ris, v, val) _which_v_indxgt_d(ris, v, val, #ris, __FILE__, __LINE__)
#else
	#define which_v_indxgt_d(ris, v, val) _which_v_indxgt_d(ris, v, val)
#endif
#ifdef MDEBUG
	VETTOREi * _which_v_indxgt_d(VETTOREi *ris, const VETTOREd *v, double val, const char *nome, const char *nomefile, int linea);
#else
	VETTOREi * _which_v_indxgt_d(VETTOREi *ris, const VETTOREd *v, double val);
#endif

//! restituisce un vettore di indici corrispondenti agli elementi del secondo vettore presenti nel primo.
/*!
	\ingroup wvett

	\param ris un vettore di tipo 'double' per il risultato
	\param v1 il primo vettore in cui cercare
	\param v2 il secondo vettore in cui cercare

	\test ris <- which(v1 == v2)
*/
#ifdef MDEBUG
	#define which_indx_vv_eq_d(ris, v1, v2) _which_indx_vv_eq_d(ris, v1, v2, #ris, __FILE__, __LINE__)
#else
	#define which_indx_vv_eq_d(ris, v1, v2) _which_indx_vv_eq_d(ris, v1, v2)
#endif
#ifdef MDEBUG
	VETTOREd * _which_indx_vv_eq_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _which_indx_vv_eq_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2);
#endif

//! calcola la differenza insiemistica tra un vettore di tipo 'double' e un altro vettore dello stesso tipo.
/*!
	\ingroup opvett

	\param v1 il vettore da cui sottrarre (non e' costante perche' si potrebbe restringere)
	\param v2 il vettore da sottrarre

	\test v1 <- setdiff(v1, v2)
*/
#ifdef MDEBUG
	#define setdiff1_d(v1, v2) _setdiff1_d(v1, v2, __FILE__, __LINE__)
#else
	#define setdiff1_d(v1, v2) _setdiff1_d(v1, v2)
#endif
#ifdef MDEBUG
	void  _setdiff1_d(VETTOREd *v1, const VETTOREd *v2, const char *nomefile, int linea);
#else
	void  _setdiff1_d(VETTOREd *v1, const VETTOREd *v2);
#endif

//! restituisce un vettore corrispondente agli indici degli elementi uguali ad un determinato valore in una certa riga di una matrice di tipo 'double'.
/*!
	\ingroup wmatr

	\param ris un vettore di tipo 'int' per il risultato
	\param m la matrice in cui cercare
	\param riga la riga di riferimento
	\param val il valore di riferimento

	\test ris <- which(m[riga, ] == val)
*/
#ifdef MDEBUG
	#define which_m_rowindxeq_d(ris, m, riga, val) _which_m_rowindxeq_d(ris, m, riga, val, #ris, __FILE__, __LINE__)
#else
	#define which_m_rowindxeq_d(ris, m, riga, val) _which_m_rowindxeq_d(ris, m, riga, val)
#endif
#ifdef MDEBUG
	VETTOREi * _which_m_rowindxeq_d(VETTOREi *ris, const MATRICEd *m, int riga, double val, const char *nome, const char *nomefile, int linea);
#else
	VETTOREi * _which_m_rowindxeq_d(VETTOREi *ris, const MATRICEd *m, int riga, double val);
#endif

//! assegna un valore agli elementi di un vettore di tipo 'double' corrispondenti agli indici contenuti in un altro vettore.
/*!
	\ingroup setvett

	\param v1 il vettore da assegnare
	\param v2 il vettore degli indici
	\param val il valore da assegnare

	\test v1[v2] <- val
*/
#ifdef MDEBUG
	#define assegna1_vs_indx_d(v1, v2, val) _assegna1_vs_indx_d(v1, v2, val, __FILE__, __LINE__)
#else
	#define assegna1_vs_indx_d(v1, v2, val) _assegna1_vs_indx_d(v1, v2, val)
#endif
#ifdef MDEBUG
	void  _assegna1_vs_indx_d(VETTOREd *v1, const VETTOREi *v2, double val, const char *nomefile, int linea);
#else
	void  _assegna1_vs_indx_d(VETTOREd *v1, const VETTOREi *v2, double val);
#endif

//! crea un vettore da un altro vettore di tipo 'double' in base ad un vettore di indici.
/*!
	\ingroup wvett

	\param ris il vettore da restutuire
	\param v il vettore di riferimento
	\param indx il vettore degli indici

	\test ris <- v[indx] # senza sforamento dei limiti
*/
#ifdef MDEBUG
	#define assegna_v_indx_d(ris, v, indx) _assegna_v_indx_d(ris, v, indx, #ris, __FILE__, __LINE__)
#else
	#define assegna_v_indx_d(ris, v, indx) _assegna_v_indx_d(ris, v, indx)
#endif
#ifdef MDEBUG
	VETTOREd * _assegna_v_indx_d(VETTOREd *ris, const VETTOREd *v, const VETTOREi *indx, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _assegna_v_indx_d(VETTOREd *ris, const VETTOREd *v, const VETTOREi *indx);
#endif

//! crea un vettore da un altro vettore di tipo 'double' in base ad un vettore di indici; nel caso in cui si esca dal limite superiore del vettore assegna NA (comportamento di R).
/*!
	\ingroup  wvett

	\param ris il vettore da restutuire
	\param v il vettore di riferimento
	\param indx il vettore degli indici

	\test ris <- v[indx] # con sforamento dei limiti
*/
#ifdef MDEBUG
	#define assegna_v_indxNA_d(ris, v, indx) _assegna_v_indxNA_d(ris, v, indx, #ris, __FILE__, __LINE__)
#else
	#define assegna_v_indxNA_d(ris, v, indx) _assegna_v_indxNA_d(ris, v, indx)
#endif
#ifdef MDEBUG
	VETTOREd * _assegna_v_indxNA_d(VETTOREd *ris, const VETTOREd *v, const VETTOREi *indx, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _assegna_v_indxNA_d(VETTOREd *ris, const VETTOREd *v, const VETTOREi *indx);
#endif

//! somma gli elementi di un vettore di tipo 'double' corrispondenti agli indici di un vettore.
/*!
	\ingroup wvett

	\param v il vettore in cui cercare
	\param indx il vettore degli indici

	\return la somma risultante

	\test ris <- sum(v[indx])
*/
#ifdef MDEBUG
	#define somma_v_indx_d(v, indx) _somma_v_indx_d(v, indx, __FILE__, __LINE__)
#else
	#define somma_v_indx_d(v, indx) _somma_v_indx_d(v, indx)
#endif
#ifdef MDEBUG
	double  _somma_v_indx_d(const VETTOREd *v, const VETTOREi *indx, const char *nomefile, int linea);
#else
	double  _somma_v_indx_d(const VETTOREd *v, const VETTOREi *indx);
#endif

//! divide un vettore di tipo 'double' "sul posto".
/*!
	\ingroup opvett

	\param n il vettore da dividere
	\param div il valore per cui dividere

	\test v <- v / div
*/
#ifdef MDEBUG
	#define dividi1_vs_d(n, div) _dividi1_vs_d(n, div, __FILE__, __LINE__)
#else
	#define dividi1_vs_d(n, div) _dividi1_vs_d(n, div)
#endif
#ifdef MDEBUG
	void  _dividi1_vs_d(VETTOREd *n, double div, const char *nomefile, int linea);
#else
	void  _dividi1_vs_d(VETTOREd *n, double div);
#endif

//! assegna un valore agli elementi di una matrice di tipo 'double' in base a due vettori di indici per le righe e per le colonne.
/*!
	\ingroup setmatr

	\param m la matrice cui assegnare i valori
	\param vr il vettore contenente gli indici per le righe
	\param vc il vettore contenente gli indici per le colonne
	\param val il valore da assegnare

	\test m[vr, vc] <- val
*/
#ifdef MDEBUG
	#define assegna1_m_vv_d(m, vr, vc, val) _assegna1_m_vv_d(m, vr, vc, val, __FILE__, __LINE__)
#else
	#define assegna1_m_vv_d(m, vr, vc, val) _assegna1_m_vv_d(m, vr, vc, val)
#endif
#ifdef MDEBUG
	void  _assegna1_m_vv_d(MATRICEd *m, const VETTOREi *vr, const VETTOREi *vc, double val, const char *nomefile, int linea);
#else
	void  _assegna1_m_vv_d(MATRICEd *m, const VETTOREi *vr, const VETTOREi *vc, double val);
#endif

//! calcola il valore assoluto degli elementi di una matrice di tipo 'double'.
/*!
	\ingroup opmatr

	\param ris la matrice risultante
	\param m la matrice d'ingresso

	\test ris <- abs(m)
*/
#ifdef MDEBUG
	#define abs_m_d(ris, m) _abs_m_d(ris, m, #ris, __FILE__, __LINE__)
#else
	#define abs_m_d(ris, m) _abs_m_d(ris, m)
#endif
#ifdef MDEBUG
	MATRICEd * _abs_m_d(MATRICEd *ris, const MATRICEd *m, const char *nome, const char *nomefile, int linea);
#else
	MATRICEd * _abs_m_d(MATRICEd *ris, const MATRICEd *m);
#endif

//! restituisce un vettore corrispondente agli indici degli elementi diversi da un determinato valore in un vettore di tipo 'double'.
/*!
	\ingroup wvett

	\param ris un vettore di tipo 'int' per il risultato
	\param v il vettore in cui cercare
	\param val il valore di riferimento

	\test ris <- which(v != val)
*/
#ifdef MDEBUG
	#define which_v_indxne_d(ris, v, val) _which_v_indxne_d(ris, v, val, #ris, __FILE__, __LINE__)
#else
	#define which_v_indxne_d(ris, v, val) _which_v_indxne_d(ris, v, val)
#endif
#ifdef MDEBUG
	VETTOREi * _which_v_indxne_d(VETTOREi *ris, const VETTOREd *v, double val, const char *nome, const char *nomefile, int linea);
#else
	VETTOREi * _which_v_indxne_d(VETTOREi *ris, const VETTOREd *v, double val);
#endif

//! restituisce un vettore corrispondente agli indici degli elementi uguali ad un determinato valore nella colonna \c c di una matrice di tipo 'double'.
/*!
	\ingroup wvett

	\param ris un vettore di tipo 'int' per il risultato
	\param m la matrice in cui cercare
	\param c la colonna di riferimento
	\param val il valore di riferimento

	\test ris <- which(m[, c] == val)
*/
#ifdef MDEBUG
	#define which_m_colindxeq_d(ris, m, c, val) _which_m_colindxeq_d(ris, m, c, val, #ris, __FILE__, __LINE__)
#else
	#define which_m_colindxeq_d(ris, m, c, val) _which_m_colindxeq_d(ris, m, c, val)
#endif
#ifdef MDEBUG
	VETTOREi * _which_m_colindxeq_d(VETTOREi *ris, const MATRICEd *m, int c, double val, const char *nome, const char *nomefile, int linea);
#else
	VETTOREi * _which_m_colindxeq_d(VETTOREi *ris, const MATRICEd *m, int c, double val);
#endif

//! interseca due vettori di tipo 'double'.
/*!
	\ingroup opvett

	\param ris un vettore per il risultato
	\param v1 il primo vettore
	\param v2 il secondo vettore

	\test  ris <- intersect(v1, v2)
*/
#ifdef MDEBUG
	#define interseca_d(ris, v1, v2) _interseca_d(ris, v1, v2, #ris, __FILE__, __LINE__)
#else
	#define interseca_d(ris, v1, v2) _interseca_d(ris, v1, v2)
#endif
#ifdef MDEBUG
	VETTOREd * _interseca_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _interseca_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2);
#endif

//! copia una riga da una matrice ad un'altra.
/*!
	\ingroup setmatr

	\param m1 la prima matrice
	\param riga1 il numero della riga di destinazione
	\param m2 la seconda matrice
	\param riga2 il numero della riga di origine

	\test m1[riga1,] <- m2[riga2,]
*/
#ifdef MDEBUG
	#define copia1_m_riga_d(m1, riga1, m2, riga2) _copia1_m_riga_d(m1, riga1, m2, riga2, __FILE__, __LINE__)
#else
	#define copia1_m_riga_d(m1, riga1, m2, riga2) _copia1_m_riga_d(m1, riga1, m2, riga2)
#endif
#ifdef MDEBUG
	void  _copia1_m_riga_d(MATRICEd *m1, int riga1, const MATRICEd *m2, int riga2, const char *nomefile, int linea);
#else
	void  _copia1_m_riga_d(MATRICEd *m1, int riga1, const MATRICEd *m2, int riga2);
#endif

//! aggiunge una riga da un'altra matrice di tipo 'double'.
/*!
	\ingroup setmatr

	\param ris la prima matrice
	\param riga1 il numero della riga di destinazione
	\param m2 la seconda matrice
	\param riga2 il numero della riga di origine

	\return la matrice modificata

	\note il parametro \c  riga1 dev'essere minore o uguale ad 1 + il numero di righe della prima matrice, mentre il parametro \c riga2 dev'essere compreso nel numero di righe della seconda matrice

	\test ris[riga1,] <- m2[riga2,] # con riga1 <= dim(m2)[0] + 1
*/
#ifdef MDEBUG
	#define aggiungi_riga_d(ris, riga1, m2, riga2) _aggiungi_riga_d(ris, riga1, m2, riga2, #ris, __FILE__, __LINE__)
#else
	#define aggiungi_riga_d(ris, riga1, m2, riga2) _aggiungi_riga_d(ris, riga1, m2, riga2)
#endif
#ifdef MDEBUG
	MATRICEd * _aggiungi_riga_d(MATRICEd *ris, int riga1, const MATRICEd *m2, int riga2, const char *nome, const char *nomefile, int linea);
#else
	MATRICEd * _aggiungi_riga_d(MATRICEd *ris, int riga1, const MATRICEd *m2, int riga2);
#endif

//! restituisce l'elemento di una matrice all'intersezione tra una colonna e un vettore di indici.
/*!
	\ingroup wmatr

	\param m la matrice di riferimento
	\param indx un vettore di indici
	\param colonna la colonna di riferimento

	\return l'elemento all'intersezione

	\test ris <- m[indx, colonna]
*/
#ifdef MDEBUG
	#define copia_m_colindx_d(m, indx, colonna) _copia_m_colindx_d(m, indx, colonna, __FILE__, __LINE__)
#else
	#define copia_m_colindx_d(m, indx, colonna) _copia_m_colindx_d(m, indx, colonna)
#endif
#ifdef MDEBUG
	double  _copia_m_colindx_d(const MATRICEd *m, const VETTOREi *indx, int colonna, const char *nomefile, int linea);
#else
	double  _copia_m_colindx_d(const MATRICEd *m, const VETTOREi *indx, int colonna);
#endif

//! assegna un valore agli elementi di una riga di una matrice di tipo 'double' corrispondenti agli indici contenuti in un vettore.
/*!
	\ingroup setmatr

	\param m la matrice da assegnare
	\param riga la riga da assegnare
	\param indx il vettore degli indici
	\param val il valore da assegnare

	\test m[riga, indx] <- val
*/
#ifdef MDEBUG
	#define assegna1_ms_rigaindx_d(m, riga, indx, val) _assegna1_ms_rigaindx_d(m, riga, indx, val, __FILE__, __LINE__)
#else
	#define assegna1_ms_rigaindx_d(m, riga, indx, val) _assegna1_ms_rigaindx_d(m, riga, indx, val)
#endif
#ifdef MDEBUG
	void  _assegna1_ms_rigaindx_d(MATRICEd *m, int riga, const VETTOREi *indx, double val, const char *nomefile, int linea);
#else
	void  _assegna1_ms_rigaindx_d(MATRICEd *m, int riga, const VETTOREi *indx, double val);
#endif

//! aggiunge una riga ad una matrice di tipo 'double' assegnando gli elementi corrispondenti agli indici contenuti in un vettore.
/*!
	\ingroup setmatr

	\param ris la matrice da assegnare
	\param riga la riga da assegnare
	\param indx il vettore degli indici
	\param val il valore da assegnare

	\note il parametro \c riga dev'essere minore o uguale ad 1 + il numero di righe della matrice

	\test ris[riga, indx] <- val # con riga <= dim(ris)[0] + 1
*/
#ifdef MDEBUG
	#define aggiungi_ms_rigaindx_d(ris, riga, indx, val) _aggiungi_ms_rigaindx_d(ris, riga, indx, val, #ris, __FILE__, __LINE__)
#else
	#define aggiungi_ms_rigaindx_d(ris, riga, indx, val) _aggiungi_ms_rigaindx_d(ris, riga, indx, val)
#endif
#ifdef MDEBUG
	MATRICEd * _aggiungi_ms_rigaindx_d(MATRICEd *ris, int riga, const VETTOREi *indx, double val, const char *nome, const char *nomefile, int linea);
#else
	MATRICEd * _aggiungi_ms_rigaindx_d(MATRICEd *ris, int riga, const VETTOREi *indx, double val);
#endif

//! assegna agli elementi di una riga di una matrice di tipo 'double' i valori di un'altra matrice, considerata come un vettore , in base a un vettore di indici.
/*!
	\ingroup setmatr

	\param m1 la matrice da assegnare
	\param riga la riga da assegnare
	\param m2 la matrice da cui copiare
	\param indx il vettore degli indici

	\test m1[riga, ] <- m2[indx]
*/
#ifdef MDEBUG
	#define assegna1_mm_rigaindx_d(m1, riga, m2, indx) _assegna1_mm_rigaindx_d(m1, riga, m2, indx, __FILE__, __LINE__)
#else
	#define assegna1_mm_rigaindx_d(m1, riga, m2, indx) _assegna1_mm_rigaindx_d(m1, riga, m2, indx)
#endif
#ifdef MDEBUG
	void  _assegna1_mm_rigaindx_d(MATRICEd *m1, int riga, const MATRICEd *m2, const VETTOREi *indx, const char *nomefile, int linea);
#else
	void  _assegna1_mm_rigaindx_d(MATRICEd *m1, int riga, const MATRICEd *m2, const VETTOREi *indx);
#endif

//! assegna un valore agli elementi di una riga di una matrice di tipo 'double'.
/*!
	\ingroup setmatr

	\param m la matrice da assegnare
	\param riga la riga da assegnare
	\param val il valore da assegnare

	\test m[riga, ] <- val
*/
#ifdef MDEBUG
	#define assegna1_ms_riga_d(m, riga, val) _assegna1_ms_riga_d(m, riga, val, __FILE__, __LINE__)
#else
	#define assegna1_ms_riga_d(m, riga, val) _assegna1_ms_riga_d(m, riga, val)
#endif
#ifdef MDEBUG
	void  _assegna1_ms_riga_d(MATRICEd *m, int riga, double val, const char *nomefile, int linea);
#else
	void  _assegna1_ms_riga_d(MATRICEd *m, int riga, double val);
#endif

//! somma le colonne di una matrice di tipo 'double'.
/*!
	\ingroup opmatr

	\param ris un vettore per il risultato
	\param m la matrice di riferimento

	\test  ris <- apply(m, 2, sum)
*/
#ifdef MDEBUG
	#define somma_colonne_d(ris, m) _somma_colonne_d(ris, m, #ris, __FILE__, __LINE__)
#else
	#define somma_colonne_d(ris, m) _somma_colonne_d(ris, m)
#endif
#ifdef MDEBUG
	VETTOREd * _somma_colonne_d(VETTOREd *ris, const MATRICEd *m, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _somma_colonne_d(VETTOREd *ris, const MATRICEd *m);
#endif

//! restituisce un vettore corrispondente agli indici degli elementi minori di un determinato valore in un vettore di tipo 'double'.
/*!
	\ingroup wvett

	\param ris un vettore di tipo 'int' per il risultato
	\param v il vettore in cui cercare
	\param val il valore di riferimento

	\test ris <- which(v < val)
*/
#ifdef MDEBUG
	#define which_v_indxlt_d(ris, v, val) _which_v_indxlt_d(ris, v, val, #ris, __FILE__, __LINE__)
#else
	#define which_v_indxlt_d(ris, v, val) _which_v_indxlt_d(ris, v, val)
#endif
#ifdef MDEBUG
	VETTOREi * _which_v_indxlt_d(VETTOREi *ris, const VETTOREd *v, double val, const char *nome, const char *nomefile, int linea);
#else
	VETTOREi * _which_v_indxlt_d(VETTOREi *ris, const VETTOREd *v, double val);
#endif

//! unisce due vettori di tipo 'double' ed elimina gli elementi ripetuti
/*!
	\ingroup opvett

	\param ris il primo vettore di partenza (non e' costante perche' potrebbe ingrandirsi)
	\param v2 il secondo vettore di partenza

	\test ris <- union(ris, v2)

	\note devo restituirlo perche' altrimenti non funziona (si potrebbe ingrandire)
*/
#ifdef MDEBUG
	#define unione1_d(ris, v2) _unione1_d(ris, v2, #ris, __FILE__, __LINE__)
#else
	#define unione1_d(ris, v2) _unione1_d(ris, v2)
#endif
#ifdef MDEBUG
	VETTOREd * _unione1_d(VETTOREd *ris, const VETTOREd *v2, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _unione1_d(VETTOREd *ris, const VETTOREd *v2);
#endif

//! somma gli elementi di un vettore di tipo 'double'.
/*!
	\ingroup opvett

	\param v il vettore di riferimento
	\param canc_NA se vero ignora gli elementi NA, altrimenti rida` NA se presenti (falso per default)

	\return la somma risultante

	\test ris <- sum(v)
*/
#ifdef MDEBUG
	#define somma_v_d(v, canc_NA) _somma_v_d(v, canc_NA, __FILE__, __LINE__)
#else
	#define somma_v_d(v, canc_NA) _somma_v_d(v, canc_NA)
#endif
#ifdef MDEBUG
	double  _somma_v_d(const VETTOREd *v, bool canc_NA, const char *nomefile, int linea);
#else
	double  _somma_v_d(const VETTOREd *v, bool canc_NA);
#endif

//! somma gli elementi di una matrice di tipo 'double'.
/*!
	\ingroup opmatr

	\param m la matrice di riferimento

	\return la somma risultante

	\test ris <- sum(m)
*/
#ifdef MDEBUG
	#define somma_m_d(m) _somma_m_d(m, __FILE__, __LINE__)
#else
	#define somma_m_d(m) _somma_m_d(m)
#endif
#ifdef MDEBUG
	double  _somma_m_d(const MATRICEd *m, const char *nomefile, int linea);
#else
	double  _somma_m_d(const MATRICEd *m);
#endif

//! calcola la differenza insiemistica tra due vettori di tipo 'double'.
/*!
	\ingroup opvett

	\param ris il vettore per il risultato
	\param v1 il vettore da cui sottrarre
	\param v2 il vettore da sottrarre

	\test ris <- setdiff(v1, v2)
*/
#ifdef MDEBUG
	#define setdiff_d(ris, v1, v2) _setdiff_d(ris, v1, v2, #ris, __FILE__, __LINE__)
#else
	#define setdiff_d(ris, v1, v2) _setdiff_d(ris, v1, v2)
#endif
#ifdef MDEBUG
	VETTOREd * _setdiff_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _setdiff_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2);
#endif

//! assegna ad una matrice una riga copiandola da un vettore
/*!
	\ingroup setmatr

	\param m la matrice di riferimento
	\param riga la riga
	\param v il vettore da assegnare

	\test m[riga,] <- v
*/
#ifdef MDEBUG
	#define assegna1_mv_riga_d(m, riga, v) _assegna1_mv_riga_d(m, riga, v, __FILE__, __LINE__)
#else
	#define assegna1_mv_riga_d(m, riga, v) _assegna1_mv_riga_d(m, riga, v)
#endif
#ifdef MDEBUG
	void  _assegna1_mv_riga_d(MATRICEd *m, int riga, const VETTOREd *v, const char *nomefile, int linea);
#else
	void  _assegna1_mv_riga_d(MATRICEd *m, int riga, const VETTOREd *v);
#endif

//! aggiunge ad una matrice una riga copiandola da un vettore
/*!
	\ingroup setmatr

	\param ris la matrice di riferimento
	\param riga la riga di riferimento
	\param v il vettore da assegnare

	\return la matrice con la nuova riga

	\note il parametro \c riga dev'essere minore o uguale ad 1 + il numero di righe della matrice

	\test ris[riga,] <- v # con riga <= dim(ris)[0] + 1
*/
#ifdef MDEBUG
	#define aggiungi_mv_riga_d(ris, riga, v) _aggiungi_mv_riga_d(ris, riga, v, #ris, __FILE__, __LINE__)
#else
	#define aggiungi_mv_riga_d(ris, riga, v) _aggiungi_mv_riga_d(ris, riga, v)
#endif
#ifdef MDEBUG
	MATRICEd * _aggiungi_mv_riga_d(MATRICEd *ris, int riga, const VETTOREd *v, const char *nome, const char *nomefile, int linea);
#else
	MATRICEd * _aggiungi_mv_riga_d(MATRICEd *ris, int riga, const VETTOREd *v);
#endif

//! aggiunge ad una matrice una riga assegnando agli elementi un valore costante
/*!
	\ingroup setmatr

	\param ris la matrice di riferimento
	\param riga la riga
	\param val il valore da assegnare

	\return la matrice con la nuova riga

	\note il parametro \c riga1 dev'essere minore o uguale ad 1 + il numero di righe della matrice

	\test ris[riga,] <- val # con riga <= dim(ris)[0] + 1
*/
#ifdef MDEBUG
	#define aggiungi_ms_riga_d(ris, riga, val) _aggiungi_ms_riga_d(ris, riga, val, #ris, __FILE__, __LINE__)
#else
	#define aggiungi_ms_riga_d(ris, riga, val) _aggiungi_ms_riga_d(ris, riga, val)
#endif
#ifdef MDEBUG
	MATRICEd * _aggiungi_ms_riga_d(MATRICEd *ris, int riga, double val, const char *nome, const char *nomefile, int linea);
#else
	MATRICEd * _aggiungi_ms_riga_d(MATRICEd *ris, int riga, double val);
#endif

//! assegna agli elementi di una matrice minori di un determinato valore un altro valore.
/*!
	\ingroup setmatr

	\param m1 la matrice di riferimento
	\param m2 la matrice da cui copiare i valori
	\param val1 il valore da confrontare
	\param val2 il valore da assegnare

	\test m1[which(m2<val1)] <- val2
*/
#ifdef MDEBUG
	#define assegna1_m_indxlt_d(m1, m2, val1, val2) _assegna1_m_indxlt_d(m1, m2, val1, val2, __FILE__, __LINE__)
#else
	#define assegna1_m_indxlt_d(m1, m2, val1, val2) _assegna1_m_indxlt_d(m1, m2, val1, val2)
#endif
#ifdef MDEBUG
	void  _assegna1_m_indxlt_d(MATRICEd *m1, MATRICEd *m2, double val1, double val2, const char *nomefile, int linea);
#else
	void  _assegna1_m_indxlt_d(MATRICEd *m1, MATRICEd *m2, double val1, double val2);
#endif

//! assegna agli elementi di una matrice maggiori di un determinato valore un altro valore.
/*!
	\ingroup setmatr

	\param m1 la matrice di riferimento
	\param m2 la matrice da cui copiare i valori
	\param val1 il valore da confrontare
	\param val2 il valore da assegnare

	\test m1[which(m2>val1)] <- val2
*/
#ifdef MDEBUG
	#define assegna1_m_indxgt_d(m1, m2, val1, val2) _assegna1_m_indxgt_d(m1, m2, val1, val2, __FILE__, __LINE__)
#else
	#define assegna1_m_indxgt_d(m1, m2, val1, val2) _assegna1_m_indxgt_d(m1, m2, val1, val2)
#endif
#ifdef MDEBUG
	void  _assegna1_m_indxgt_d(MATRICEd *m1, MATRICEd *m2, double val1, double val2, const char *nomefile, int linea);
#else
	void  _assegna1_m_indxgt_d(MATRICEd *m1, MATRICEd *m2, double val1, double val2);
#endif

//! ordina gli elementi di un vettore di tipo 'double' e ne restituisce le posizioni in un vettore
/*!
	\ingroup opvett

	\param ris un vettore per il risultato
	\param v il vettore da ordinare
	\param decr se vero ordina in maniera decrescente

	\return un vettore con gli indici degli elementi ordinati

	\test ris <- order(v, decr)
*/
#ifdef MDEBUG
	#define ordine_d(ris, v, decr) _ordine_d(ris, v, decr, #ris, __FILE__, __LINE__)
#else
	#define ordine_d(ris, v, decr) _ordine_d(ris, v, decr)
#endif
#ifdef MDEBUG
	VETTOREi * _ordine_d(VETTOREi *ris, const VETTOREd *v, bool decr, const char *nome, const char *nomefile, int linea);
#else
	VETTOREi * _ordine_d(VETTOREi *ris, const VETTOREd *v, bool decr);
#endif

//! sottrae due vettori di tipo 'double' elemento per elemento.
/*!
	\ingroup opvett

	\param ris un vettore per il risultato
	\param v1 il primo vettore
	\param v2 il secondo vettore

	\test  ris <- v1 - v2
*/
#ifdef MDEBUG
	#define diff_vv_d(ris, v1, v2) _diff_vv_d(ris, v1, v2, #ris, __FILE__, __LINE__)
#else
	#define diff_vv_d(ris, v1, v2) _diff_vv_d(ris, v1, v2)
#endif
#ifdef MDEBUG
	VETTOREd * _diff_vv_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _diff_vv_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2);
#endif

//! calcola il valore assoluto degli elementi di un vettore di tipo 'double'.
/*!
	\ingroup opvett

	\param ris il vettore risultante
	\param m il vettore d'ingresso

	\test ris <- abs(v)
*/
#ifdef MDEBUG
	#define abs_v_d(ris, m) _abs_v_d(ris, m, #ris, __FILE__, __LINE__)
#else
	#define abs_v_d(ris, m) _abs_v_d(ris, m)
#endif
#ifdef MDEBUG
	VETTOREd * _abs_v_d(VETTOREd *ris, const VETTOREd *m, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _abs_v_d(VETTOREd *ris, const VETTOREd *m);
#endif

//! restituisce un vettore corrispondente agli indici degli elementi uguali o diversi da NA.
/*!
	\ingroup wvett

	\param ris un vettore di tipo 'int' per il risultato
	\param v il vettore in cui cercare
	\param complemento vero se servono i non NA

	\test ris <- which(is_na(v)) oppure which(!is_na(v))
*/
#ifdef MDEBUG
	#define which_v_indxNA_d(ris, v, complemento) _which_v_indxNA_d(ris, v, complemento, #ris, __FILE__, __LINE__)
#else
	#define which_v_indxNA_d(ris, v, complemento) _which_v_indxNA_d(ris, v, complemento)
#endif
#ifdef MDEBUG
	VETTOREi * _which_v_indxNA_d(VETTOREi *ris, const VETTOREd *v, bool complemento, const char *nome, const char *nomefile, int linea);
#else
	VETTOREi * _which_v_indxNA_d(VETTOREi *ris, const VETTOREd *v, bool complemento);
#endif

//! restituisce un vettore corrispondente agli indici degli elementi  uguali ad un determinato valore in un vettore di tipo 'double'.
/*!
	\ingroup wvett

	\param ris un vettore di tipo 'int' per il risultato
	\param v il vettore in cui cercare
	\param val il valore di riferimento

	\test ris <- which(v == val)
*/
#ifdef MDEBUG
	#define which_v_indxeq_d(ris, v, val) _which_v_indxeq_d(ris, v, val, #ris, __FILE__, __LINE__)
#else
	#define which_v_indxeq_d(ris, v, val) _which_v_indxeq_d(ris, v, val)
#endif
#ifdef MDEBUG
	VETTOREi * _which_v_indxeq_d(VETTOREi *ris, const VETTOREd *v, double val, const char *nome, const char *nomefile, int linea);
#else
	VETTOREi * _which_v_indxeq_d(VETTOREi *ris, const VETTOREd *v, double val);
#endif

//! affianca due vettori colonna creando una matrice.
/*!
	\ingroup wvett

	\param ris una matrice per il risultato
	\param v1 il primo vettore
	\param v2 il secondo vettore

	\test  ris <- cbind(v1, v2)
*/
#ifdef MDEBUG
	#define cbind2v_d(ris, v1, v2) _cbind2v_d(ris, v1, v2, #ris, __FILE__, __LINE__)
#else
	#define cbind2v_d(ris, v1, v2) _cbind2v_d(ris, v1, v2)
#endif
#ifdef MDEBUG
	MATRICEd * _cbind2v_d(MATRICEd *ris, const VETTOREd *v1, const VETTOREd *v2, const char *nome, const char *nomefile, int linea);
#else
	MATRICEd * _cbind2v_d(MATRICEd *ris, const VETTOREd *v1, const VETTOREd *v2);
#endif

//! calcola il massimo per ciascuna riga di una matrice di tipo 'double'.
/*!
	\ingroup opmatr

	\param ris un vettore per il risultato
	\param m la matrice di riferimento

	\test  ris <- apply(m, 1, max)
*/
#ifdef MDEBUG
	#define max_righe_d(ris, m) _max_righe_d(ris, m, #ris, __FILE__, __LINE__)
#else
	#define max_righe_d(ris, m) _max_righe_d(ris, m)
#endif
#ifdef MDEBUG
	VETTOREd * _max_righe_d(VETTOREd *ris, const MATRICEd *m, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _max_righe_d(VETTOREd *ris, const MATRICEd *m);
#endif

//! calcola una particolare funzione vettoriale usata in scoremodular.
/*!
	\ingroup aus

	\param ris un vettore per il risultato
	\param a  parametro della funzione
	\param b  parametro della funzione
	\param m  parametro della funzione
	\param t parametro della funzione

	\test  ris <- (sign(a - b)) * m / t
*/
#ifdef MDEBUG
	#define f_aux_d(ris, a, b, m, t) _f_aux_d(ris, a, b, m, t, #ris, __FILE__, __LINE__)
#else
	#define f_aux_d(ris, a, b, m, t) _f_aux_d(ris, a, b, m, t)
#endif
#ifdef MDEBUG
	VETTOREd * _f_aux_d(VETTOREd *ris, const VETTOREd *a, const VETTOREd *b, const VETTOREd *m, const VETTOREd *t, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _f_aux_d(VETTOREd *ris, const VETTOREd *a, const VETTOREd *b, const VETTOREd *m, const VETTOREd *t);
#endif

//! restituisce un vettore corrispondente agli indici degli elementi maggiori di un determinato valore in un vettore di tipo 'double'.
/*!
	\ingroup wvett

	\param ris un vettore di tipo 'int' per il risultato
	\param v il vettore in cui cercare
	\param val il valore di riferimento

	\test ris <- which(v > val)
*/
#ifdef MDEBUG
	#define which_v_indxgt_d(ris, v, val) _which_v_indxgt_d(ris, v, val, #ris, __FILE__, __LINE__)
#else
	#define which_v_indxgt_d(ris, v, val) _which_v_indxgt_d(ris, v, val)
#endif
#ifdef MDEBUG
	VETTOREi * _which_v_indxgt_d(VETTOREi *ris, const VETTOREd *v, double val, const char *nome, const char *nomefile, int linea);
#else
	VETTOREi * _which_v_indxgt_d(VETTOREi *ris, const VETTOREd *v, double val);
#endif

//! somma due vettori di tipo 'double'.
/*!
	\ingroup opvett

	\param ris un vettore per il risultato
	\param v1 il primo vettore
	\param v2 il secondo vettore

	\test  ris <- v1 + v2
*/
#ifdef MDEBUG
	#define somma_vv_d(ris, v1, v2) _somma_vv_d(ris, v1, v2, #ris, __FILE__, __LINE__)
#else
	#define somma_vv_d(ris, v1, v2) _somma_vv_d(ris, v1, v2)
#endif
#ifdef MDEBUG
	VETTOREd * _somma_vv_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _somma_vv_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2);
#endif

//! traspone la matrice in ingresso.
/*!
	\ingroup opmatr

	\param ris la matrice risultante
	\param m la matrice d'ingresso

	\test ris <- t(m)
*/
#ifdef MDEBUG
	#define trasponi_d(ris, m) _trasponi_d(ris, m, #ris, __FILE__, __LINE__)
#else
	#define trasponi_d(ris, m) _trasponi_d(ris, m)
#endif
#ifdef MDEBUG
	MATRICEd * _trasponi_d(MATRICEd *ris, const MATRICEd *m, const char *nome, const char *nomefile, int linea);
#else
	MATRICEd * _trasponi_d(MATRICEd *ris, const MATRICEd *m);
#endif

//! somma una matrice di tipo 'double' ad un'altra.
/*!
	\ingroup opmatr

	\param m1 la prima matrice
	\param m2 la seconda matrice

	\test  m1 <- m1 + m2
*/
#ifdef MDEBUG
	#define somma1_m_d(m1, m2) _somma1_m_d(m1, m2, __FILE__, __LINE__)
#else
	#define somma1_m_d(m1, m2) _somma1_m_d(m1, m2)
#endif
#ifdef MDEBUG
	void  _somma1_m_d(MATRICEd *m1, const MATRICEd *m2, const char *nomefile, int linea);
#else
	void  _somma1_m_d(MATRICEd *m1, const MATRICEd *m2);
#endif

//! assegna un valore agli elementi di una matrice di tipo 'double' corrispondenti agli indici contenuti in un vettore.
/*!
	\ingroup setmatr

	\param m la matrice da assegnare
	\param indx il vettore degli indici
	\param val il valore da assegnare

	\test m[indx] <- val
*/
#ifdef MDEBUG
	#define assegna1_ms_indx_d(m, indx, val) _assegna1_ms_indx_d(m, indx, val, __FILE__, __LINE__)
#else
	#define assegna1_ms_indx_d(m, indx, val) _assegna1_ms_indx_d(m, indx, val)
#endif
#ifdef MDEBUG
	void  _assegna1_ms_indx_d(MATRICEd *m, const VETTOREi *indx, double val, const char *nomefile, int linea);
#else
	void  _assegna1_ms_indx_d(MATRICEd *m, const VETTOREi *indx, double val);
#endif

//! restituisce un vettore di tipo 'double' corrispondente agli indici degli elementi diversi da val nella matrice m.
/*!
	\ingroup wmatr

	\param ris un vettore di tipo 'int' per il risultato
	\param m la matrice di riferimento
	\param val il valore di riferimento

	\test ris <- which(m != val)
*/
#ifdef MDEBUG
	#define which_m_indxne_d(ris, m, val) _which_m_indxne_d(ris, m, val, #ris, __FILE__, __LINE__)
#else
	#define which_m_indxne_d(ris, m, val) _which_m_indxne_d(ris, m, val)
#endif
#ifdef MDEBUG
	VETTOREi * _which_m_indxne_d(VETTOREi *ris, const MATRICEd *m, double val, const char *nome, const char *nomefile, int linea);
#else
	VETTOREi * _which_m_indxne_d(VETTOREi *ris, const MATRICEd *m, double val);
#endif

//! assegna un valore alla diagonale una matrice di tipo 'double'.
/*!
	\ingroup setmatr

	\param m la matrice di riferimento
	\param val il valore di riferimento

	\test diag(m) <- val
*/
#ifdef MDEBUG
	#define assegna1_s_diag_d(m, val) _assegna1_s_diag_d(m, val, __FILE__, __LINE__)
#else
	#define assegna1_s_diag_d(m, val) _assegna1_s_diag_d(m, val)
#endif
#ifdef MDEBUG
	void  _assegna1_s_diag_d(MATRICEd *m, double val, const char *nomefile, int linea);
#else
	void  _assegna1_s_diag_d(MATRICEd *m, double val);
#endif

//! calcola la media degli elementi di un vettore di tipo 'double' saltando i NA
/*!
	\ingroup opvett

	\param v il vettore in cui cercare

	\return la media risultante

	\test ris <- mean(v)
*/
#ifdef MDEBUG
	#define media_v_d(v) _media_v_d(v, __FILE__, __LINE__)
#else
	#define media_v_d(v) _media_v_d(v)
#endif
#ifdef MDEBUG
	double  _media_v_d(VETTOREd *v, const char *nomefile, int linea);
#else
	double  _media_v_d(VETTOREd *v);
#endif

//! calcola la somma di una riga di una matrice di tipo 'double'.
/*!
	\ingroup wmatr

	\param m la matrice di riferimento
	\param riga la riga di riferimento

	\return la somma della riga

	\test  ris <- sum(m[riga, ])
*/
#ifdef MDEBUG
	#define somma_riga_d(m, riga) _somma_riga_d(m, riga, __FILE__, __LINE__)
#else
	#define somma_riga_d(m, riga) _somma_riga_d(m, riga)
#endif
#ifdef MDEBUG
	double  _somma_riga_d(const MATRICEd *m, int riga, const char *nomefile, int linea);
#else
	double  _somma_riga_d(const MATRICEd *m, int riga);
#endif

//! restituisce un vettore di tipo 'double' con i valori della diagonale.
/*!
	\ingroup wmatr

	\param ris un vettore di tipo 'double' per il risultato
	\param m la matrice di riferimento

	\return il vettore contenente la diagonale

	\test ris <- diag(m)
*/
#ifdef MDEBUG
	#define diag_d(ris, m) _diag_d(ris, m, #ris, __FILE__, __LINE__)
#else
	#define diag_d(ris, m) _diag_d(ris, m)
#endif
#ifdef MDEBUG
	VETTOREd * _diag_d(VETTOREd *ris, const MATRICEd *m, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _diag_d(VETTOREd *ris, const MATRICEd *m);
#endif

//! restituisce un vettore di tipo 'double' corrispondente agli indici degli elementi diversi da val nella colonna \c c della matrice m.
/*!
	\ingroup wmatr

	\param ris un vettore di tipo 'int' per il risultato
	\param m la matrice di riferimento
	\param c la colonna di riferimento
	\param val il valore di riferimento

	\test ris <- which(m[, c] != val)
*/
#ifdef MDEBUG
	#define which_m_colindxne_d(ris, m, c, val) _which_m_colindxne_d(ris, m, c, val, #ris, __FILE__, __LINE__)
#else
	#define which_m_colindxne_d(ris, m, c, val) _which_m_colindxne_d(ris, m, c, val)
#endif
#ifdef MDEBUG
	VETTOREi * _which_m_colindxne_d(VETTOREi *ris, const MATRICEd *m, int c, double val, const char *nome, const char *nomefile, int linea);
#else
	VETTOREi * _which_m_colindxne_d(VETTOREi *ris, const MATRICEd *m, int c, double val);
#endif

//! restituisce un vettore corrispondente agli indici degli elementi presenti in un altro vettore nella colonna \c c di una matrice di tipo 'double'.
/*!
	\ingroup wmatr

	\param ris un vettore di tipo 'int' per il risultato
	\param m la matrice in cui cercare
	\param c la colonna di riferimento
	\param v il vettore in cui cercare

	\test ris <- which(m[, c] %IN% v)
*/
#ifdef MDEBUG
	#define which_m_colindxin_d(ris, m, c, v) _which_m_colindxin_d(ris, m, c, v, #ris, __FILE__, __LINE__)
#else
	#define which_m_colindxin_d(ris, m, c, v) _which_m_colindxin_d(ris, m, c, v)
#endif
#ifdef MDEBUG
	VETTOREi * _which_m_colindxin_d(VETTOREi *ris, const MATRICEd *m, int c, const VETTOREd *v, const char *nome, const char *nomefile, int linea);
#else
	VETTOREi * _which_m_colindxin_d(VETTOREi *ris, const MATRICEd *m, int c, const VETTOREd *v);
#endif

//! crea una matrice di tipo 'double' costituita da \c c colonne a partire da una matrice di tipo 'double' letta come un vettore.
/*!
	\ingroup oggetti

	\param ris la matrice per il risultato
	\param m la matrice di riferimento
	\param c il numero di colonne della nuova matrice

	\return la nuova matrice

	\test ris <- matrix(m, ncol=c)
*/
#ifdef MDEBUG
	#define copia_m_ncol_d(ris, m, c) _copia_m_ncol_d(ris, m, c, #ris, __FILE__, __LINE__)
#else
	#define copia_m_ncol_d(ris, m, c) _copia_m_ncol_d(ris, m, c)
#endif
#ifdef MDEBUG
	MATRICEd * _copia_m_ncol_d(MATRICEd *ris, const MATRICEd *m, int c, const char *nome, const char *nomefile, int linea);
#else
	MATRICEd * _copia_m_ncol_d(MATRICEd *ris, const MATRICEd *m, int c);
#endif

//! restituisce un vettore corrispondente agli indici degli elementi per cui vale la formula sotto riportata per gli elementi della colonna \c c di una matrice di tipo 'double'.
/*!
	\ingroup wmatr

	\param ris un vettore di tipo 'int' per il risultato
	\param m la matrice in cui cercare
	\param c la colonna di riferimento
	\param val1 il primo valore di confronto
	\param val2 il secondo valore di confronto

	\test ris <- which((m[, c] != val1) & (m[, c] != val2))
*/
#ifdef MDEBUG
	#define which_m_colneand2_d(ris, m, c, val1, val2) _which_m_colneand2_d(ris, m, c, val1, val2, #ris, __FILE__, __LINE__)
#else
	#define which_m_colneand2_d(ris, m, c, val1, val2) _which_m_colneand2_d(ris, m, c, val1, val2)
#endif
#ifdef MDEBUG
	VETTOREi * _which_m_colneand2_d(VETTOREi *ris, const MATRICEd *m, int c, double val1, double val2, const char *nome, const char *nomefile, int linea);
#else
	VETTOREi * _which_m_colneand2_d(VETTOREi *ris, const MATRICEd *m, int c, double val1, double val2);
#endif

//! estrae le righe da una matrice di tipo 'double' in base ad un vettore di indici
/*!
	\ingroup wmatr

	\param ris un vettore per il risultato
	\param m la matrice di riferimento
	\param indx il vettore di indici

	\test ris <- M[righe,]
*/
#ifdef MDEBUG
	#define righe_d(ris, m, indx) _righe_d(ris, m, indx, #ris, __FILE__, __LINE__)
#else
	#define righe_d(ris, m, indx) _righe_d(ris, m, indx)
#endif
#ifdef MDEBUG
	MATRICEd * _righe_d(MATRICEd *ris, const MATRICEd *m, const VETTOREi *indx, const char *nome, const char *nomefile, int linea);
#else
	MATRICEd * _righe_d(MATRICEd *ris, const MATRICEd *m, const VETTOREi *indx);
#endif

//! assegna un vettore alla diagonale una matrice di tipo 'double'.
/*!
	\ingroup setmatr

	\param m la matrice di riferimento
	\param v il vettore di riferimento

	\test diag(m) <- v
*/
#ifdef MDEBUG
	#define assegna1_v_diag_d(m, v) _assegna1_v_diag_d(m, v, __FILE__, __LINE__)
#else
	#define assegna1_v_diag_d(m, v) _assegna1_v_diag_d(m, v)
#endif
#ifdef MDEBUG
	void  _assegna1_v_diag_d(MATRICEd *m, const VETTOREd *v, const char *nomefile, int linea);
#else
	void  _assegna1_v_diag_d(MATRICEd *m, const VETTOREd *v);
#endif

//! assegna un valore agli elementi di una matrice di tipo 'double' corrispondenti agli indici di riga e colonna contenuti in una matrice di interi r x 2.
/*!
	\ingroup setmatr

	\param m la matrice da assegnare
	\param indxm la matrice di coppie di indici
	\param val il valore da assegnare

	\test m[indx] <- val
*/
#ifdef MDEBUG
	#define assegna1_ms_indx2_d(m, indxm, val) _assegna1_ms_indx2_d(m, indxm, val, __FILE__, __LINE__)
#else
	#define assegna1_ms_indx2_d(m, indxm, val) _assegna1_ms_indx2_d(m, indxm, val)
#endif
#ifdef MDEBUG
	void  _assegna1_ms_indx2_d(MATRICEd *m, const MATRICEi *indxm, double val, const char *nomefile, int linea);
#else
	void  _assegna1_ms_indx2_d(MATRICEd *m, const MATRICEi *indxm, double val);
#endif

//! incrementa (o decrementa) un vettore di tipo 'double' di \c n.
/*!
	\ingroup opvett

	\param v il vettore di riferimento
	\param s il valore da addizionare

	\test ris <- v + s
*/
#ifdef MDEBUG
	#define incr1_v_d(v, s) _incr1_v_d(v, s, __FILE__, __LINE__)
#else
	#define incr1_v_d(v, s) _incr1_v_d(v, s)
#endif
#ifdef MDEBUG
	void  _incr1_v_d(VETTOREd *v, double s, const char *nomefile, int linea);
#else
	void  _incr1_v_d(VETTOREd *v, double s);
#endif

//! copia una matrice di tipo 'double' in un altra già allocata.
/*!
	\ingroup opmatr

	\param ris la matrice di destinazione
	\param da la matrice sorgente

	\return la matrice copiata

	\test ris <- da
*/
#ifdef MDEBUG
	#define copia_m_d(ris, da) _copia_m_d(ris, da, #ris, __FILE__, __LINE__)
#else
	#define copia_m_d(ris, da) _copia_m_d(ris, da)
#endif
#ifdef MDEBUG
	MATRICEd * _copia_m_d(MATRICEd *ris, const MATRICEd *da, const char *nome, const char *nomefile, int linea);
#else
	MATRICEd * _copia_m_d(MATRICEd *ris, const MATRICEd *da);
#endif

//! somma due matrici di tipo 'double'.
/*!
	\ingroup opmatr

	\param ris la matrice risultante
	\param m1 la prima matrice
	\param m2 la seconda matrice

	\return la matrice risultante

	\test  ris <- m1 + m2
*/
#ifdef MDEBUG
	#define somma_mm_d(ris, m1, m2) _somma_mm_d(ris, m1, m2, #ris, __FILE__, __LINE__)
#else
	#define somma_mm_d(ris, m1, m2) _somma_mm_d(ris, m1, m2)
#endif
#ifdef MDEBUG
	MATRICEd * _somma_mm_d(MATRICEd *ris, const MATRICEd *m1, const MATRICEd *m2, const char *nome, const char *nomefile, int linea);
#else
	MATRICEd * _somma_mm_d(MATRICEd *ris, const MATRICEd *m1, const MATRICEd *m2);
#endif

//! accoda uno scalare di tipo 'double' ad un vettore.
/*!
	\ingroup opvett

	\param ris il vettore di riferimento (non e' costante perche' potrebbe ingrandirsi)
	\param s lo scalare da accodare

	\test ris <- c(ris, s)
*/
#ifdef MDEBUG
	#define accoda1_vs_d(ris, s) _accoda1_vs_d(ris, s, #ris, __FILE__, __LINE__)
#else
	#define accoda1_vs_d(ris, s) _accoda1_vs_d(ris, s)
#endif
#ifdef MDEBUG
	VETTOREd * _accoda1_vs_d(VETTOREd *ris, double s, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _accoda1_vs_d(VETTOREd *ris, double s);
#endif

//! sostituisce un vettore con un suo segmento.
/*!
	\ingroup sequenze

	\param v il vettore di riferimento
	\param st l'indice di partenza
	\param end l'indice di arrivo

	\test v <- v[st:end]
*/
#ifdef MDEBUG
	#define segmento1_v_d(v, st, end) _segmento1_v_d(v, st, end, __FILE__, __LINE__)
#else
	#define segmento1_v_d(v, st, end) _segmento1_v_d(v, st, end)
#endif
#ifdef MDEBUG
	void  _segmento1_v_d(VETTOREd *v, int st, int end, const char *nomefile, int linea);
#else
	void  _segmento1_v_d(VETTOREd *v, int st, int end);
#endif

//! restituisce un segmento di un vettore di tipo 'double'.
/*!
	\ingroup sequenze

	\param ris un vettore per il risultato
	\param v il vettore di riferimento
	\param st l'indice di partenza
	\param end l'indice di arrivo

	\return il nuovo vettore

	\test ris <- v[st:end]
*/
#ifdef MDEBUG
	#define segmento_v_d(ris, v, st, end) _segmento_v_d(ris, v, st, end, #ris, __FILE__, __LINE__)
#else
	#define segmento_v_d(ris, v, st, end) _segmento_v_d(ris, v, st, end)
#endif
#ifdef MDEBUG
	VETTOREd * _segmento_v_d(VETTOREd *ris, const VETTOREd *v, int st, int end, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _segmento_v_d(VETTOREd *ris, const VETTOREd *v, int st, int end);
#endif

//! divide un vettore di tipo 'double' per una costante di tipo 'double'.
/*!
	\ingroup opvett

	\param ris il vettore per il risultato
	\param n il vettore da dividere
	\param div il valore per cui dividere

	\return il vettore risultante

	\test ris <- v / div
*/
#ifdef MDEBUG
	#define dividi_vs_d(ris, n, div) _dividi_vs_d(ris, n, div, #ris, __FILE__, __LINE__)
#else
	#define dividi_vs_d(ris, n, div) _dividi_vs_d(ris, n, div)
#endif
#ifdef MDEBUG
	VETTOREd * _dividi_vs_d(VETTOREd *ris, const VETTOREd *n, double div, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _dividi_vs_d(VETTOREd *ris, const VETTOREd *n, double div);
#endif

//! assegna un valore ad un vettore di tipo 'double' in base ad un indice.
/*!
	\ingroup setvett

	\param v il vettore di riferimento
	\param indx l'indice da assegnare
	\param val il valore da assegnare

	\test v[indx] <- val
*/
#ifdef MDEBUG
	#define assegna1_v_indx_d(v, indx, val) _assegna1_v_indx_d(v, indx, val, __FILE__, __LINE__)
#else
	#define assegna1_v_indx_d(v, indx, val) _assegna1_v_indx_d(v, indx, val)
#endif
#ifdef MDEBUG
	void  _assegna1_v_indx_d(VETTOREd *v, const VETTOREi *indx, double val, const char *nomefile, int linea);
#else
	void  _assegna1_v_indx_d(VETTOREd *v, const VETTOREi *indx, double val);
#endif

//! restituisce un vettore di tipo 'double' corrispondente agli indici degli elementi uguali a val nella matrice m.
/*!
	\ingroup wmatr

	\param ris un vettore di tipo 'int' per il risultato
	\param m la matrice di riferimento
	\param val il valore di riferimento

	\test ris <- which(m == val)
*/
#ifdef MDEBUG
	#define which_m_indxeq_d(ris, m, val) _which_m_indxeq_d(ris, m, val, #ris, __FILE__, __LINE__)
#else
	#define which_m_indxeq_d(ris, m, val) _which_m_indxeq_d(ris, m, val)
#endif
#ifdef MDEBUG
	VETTOREi * _which_m_indxeq_d(VETTOREi *ris, const MATRICEd *m, double val, const char *nome, const char *nomefile, int linea);
#else
	VETTOREi * _which_m_indxeq_d(VETTOREi *ris, const MATRICEd *m, double val);
#endif

//! restituisce un vettore corrispondente agli indici degli elementi per cui vale la formula sotto riportata per un vettore di tipo 'double'.
/*!
	\ingroup wvett

	\param ris un vettore di tipo 'int' per il risultato
	\param v il vettore in cui cercare
	\param val1 il primo valore di confronto
	\param val2 il secondo valore di confronto

	\return il vettore degli indici

	\test ris <- which((v < val1) & (v > val2))
*/
#ifdef MDEBUG
	#define which_v_andglt_d(ris, v, val1, val2) _which_v_andglt_d(ris, v, val1, val2, #ris, __FILE__, __LINE__)
#else
	#define which_v_andglt_d(ris, v, val1, val2) _which_v_andglt_d(ris, v, val1, val2)
#endif
#ifdef MDEBUG
	VETTOREi * _which_v_andglt_d(VETTOREi *ris, const VETTOREd *v, double val1, double val2, const char *nome, const char *nomefile, int linea);
#else
	VETTOREi * _which_v_andglt_d(VETTOREi *ris, const VETTOREd *v, double val1, double val2);
#endif

//! restituisce un vettore delle distanze euclidee tra due vettori di punti di tipo 'double'.
/*!
	\ingroup opvett

	\param ris un vettore per il risultato
	\param x1 le x del primo vettore
	\param y1 le y del primo vettore
	\param x2 le x del secondo vettore
	\param y2 le y del secondo vettore

	\return il vettore con le distanze

	\test ris <- sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2)
*/
#ifdef MDEBUG
	#define distanza_2dvv_d(ris, x1, y1, x2, y2) _distanza_2dvv_d(ris, x1, y1, x2, y2, #ris, __FILE__, __LINE__)
#else
	#define distanza_2dvv_d(ris, x1, y1, x2, y2) _distanza_2dvv_d(ris, x1, y1, x2, y2)
#endif
#ifdef MDEBUG
	VETTOREd * _distanza_2dvv_d(VETTOREd *ris, const VETTOREd *x1, const VETTOREd *y1, const VETTOREd *x2, const VETTOREd *y2, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _distanza_2dvv_d(VETTOREd *ris, const VETTOREd *x1, const VETTOREd *y1, const VETTOREd *x2, const VETTOREd *y2);
#endif

//! restituisce un vettore delle distanze euclidee tra un vettore di punti di tipo 'double' ed un punto.
/*!
	\ingroup opvett

	\param ris un vettore per il risultato
	\param x1 le x del primo vettore
	\param y1 le y del primo vettore
	\param x2 la x del secondo punto
	\param y2 la y del secondo punto

	\return il vettore con le distanze

	\test ris <- sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2)
*/
#ifdef MDEBUG
	#define distanza_2dvs_d(ris, x1, y1, x2, y2) _distanza_2dvs_d(ris, x1, y1, x2, y2, #ris, __FILE__, __LINE__)
#else
	#define distanza_2dvs_d(ris, x1, y1, x2, y2) _distanza_2dvs_d(ris, x1, y1, x2, y2)
#endif
#ifdef MDEBUG
	VETTOREd * _distanza_2dvs_d(VETTOREd *ris, const VETTOREd *x1, const VETTOREd *y1, double x2, double y2, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _distanza_2dvs_d(VETTOREd *ris, const VETTOREd *x1, const VETTOREd *y1, double x2, double y2);
#endif

//! assegna un valore agli elementi di una colonna di una matrice di tipo 'double' corrispondenti agli indici contenuti in un vettore.
/*!
	\ingroup setmatr

	\param m la matrice da assegnare
	\param indx il vettore degli indici
	\param col la colonna da assegnare
	\param val il valore da assegnare

	\test m[indx, col] <- val
*/
#ifdef MDEBUG
	#define assegna1_ms_indxcol_d(m, indx, col, val) _assegna1_ms_indxcol_d(m, indx, col, val, __FILE__, __LINE__)
#else
	#define assegna1_ms_indxcol_d(m, indx, col, val) _assegna1_ms_indxcol_d(m, indx, col, val)
#endif
#ifdef MDEBUG
	void  _assegna1_ms_indxcol_d(MATRICEd *m, const VETTOREi *indx, int col, double val, const char *nomefile, int linea);
#else
	void  _assegna1_ms_indxcol_d(MATRICEd *m, const VETTOREi *indx, int col, double val);
#endif

//! assegna un vettore agli elementi di una matrice di tipo 'double' corrispondenti agli indici contenuti in un vettore.
/*!
	\ingroup setmatr

	\param m la matrice da assegnare
	\param indx il vettore degli indici
	\param v il vettore da assegnare

	\test m[indx] <- v
*/
#ifdef MDEBUG
	#define assegna1_mv_indx_d(m, indx, v) _assegna1_mv_indx_d(m, indx, v, __FILE__, __LINE__)
#else
	#define assegna1_mv_indx_d(m, indx, v) _assegna1_mv_indx_d(m, indx, v)
#endif
#ifdef MDEBUG
	void  _assegna1_mv_indx_d(MATRICEd *m, const VETTOREi *indx, const VETTOREd *v, const char *nomefile, int linea);
#else
	void  _assegna1_mv_indx_d(MATRICEd *m, const VETTOREi *indx, const VETTOREd *v);
#endif

//! converte un vettore di tipo 'double' in un vettore di tipo 'int' (la versione intera non fa nulla).
/*!
	\ingroup opvett

	\param ris il vettore da restituire
	\param v il vettore di riferimento

	\return il vettore arrotondato

	\test implicita
*/
#ifdef MDEBUG
	#define arrotonda_v_d(ris, v) _arrotonda_v_d(ris, v, #ris, __FILE__, __LINE__)
#else
	#define arrotonda_v_d(ris, v) _arrotonda_v_d(ris, v)
#endif
#ifdef MDEBUG
	VETTOREi * _arrotonda_v_d(VETTOREi *ris, const VETTOREd *v, const char *nome, const char *nomefile, int linea);
#else
	VETTOREi * _arrotonda_v_d(VETTOREi *ris, const VETTOREd *v);
#endif

//! eleva gli elementi di un vettore di tipo 'double' ad un valore.
/*!
	\ingroup opvett

	\param ris il vettore da restituire
	\param v il vettore di riferimento
	\param val l'esponente

	\return il vettore risultante

	\test ris <- v ^ val
*/
#ifdef MDEBUG
	#define exp_d(ris, v, val) _exp_d(ris, v, val, #ris, __FILE__, __LINE__)
#else
	#define exp_d(ris, v, val) _exp_d(ris, v, val)
#endif
#ifdef MDEBUG
	VETTOREd * _exp_d(VETTOREd *ris, const VETTOREd *v, double val, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _exp_d(VETTOREd *ris, const VETTOREd *v, double val);
#endif

//! moltiplica un vettore di tipo 'double' per un valore di tipo 'double'.
/*!
	\ingroup opvett

	\param ris un vettore per il risultato
	\param v il vettore da moltiplicare
	\param val il fattore (scalare)

	\test ris <- v * val
*/
#ifdef MDEBUG
	#define moltiplica_vs_d(ris, v, val) _moltiplica_vs_d(ris, v, val, #ris, __FILE__, __LINE__)
#else
	#define moltiplica_vs_d(ris, v, val) _moltiplica_vs_d(ris, v, val)
#endif
#ifdef MDEBUG
	VETTOREd * _moltiplica_vs_d(VETTOREd *ris, const VETTOREd *v, double val, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _moltiplica_vs_d(VETTOREd *ris, const VETTOREd *v, double val);
#endif

//! assegna ad un segmento di un vettore di tipo 'double' un valore 'double' costante.
/*!
	\ingroup sequenze

	\param v il vettore di riferimento
	\param st l'indice di partenza
	\param end l'indice di arrivo
	\param val la costante da assegnare

	\test v[st:end] <- val
*/
#ifdef MDEBUG
	#define assegna1_v_segm_d(v, st, end, val) _assegna1_v_segm_d(v, st, end, val, __FILE__, __LINE__)
#else
	#define assegna1_v_segm_d(v, st, end, val) _assegna1_v_segm_d(v, st, end, val)
#endif
#ifdef MDEBUG
	void  _assegna1_v_segm_d(VETTOREd *v, int st, int end, double val, const char *nomefile, int linea);
#else
	void  _assegna1_v_segm_d(VETTOREd *v, int st, int end, double val);
#endif

//! assegna ad un segmento di un vettore di tipo 'double' un altro vettore.
/*!
	\ingroup sequenze

	\param v1 il vettore di riferimento
	\param st l'indice di partenza
	\param end l'indice di arrivo
	\param v2 il vettore da assegnare

	\test v1[st:end] <- v2
*/
#ifdef MDEBUG
	#define assegna1_v_segmv_d(v1, st, end, v2) _assegna1_v_segmv_d(v1, st, end, v2, __FILE__, __LINE__)
#else
	#define assegna1_v_segmv_d(v1, st, end, v2) _assegna1_v_segmv_d(v1, st, end, v2)
#endif
#ifdef MDEBUG
	void  _assegna1_v_segmv_d(VETTOREd *v1, int st, int end, const VETTOREd *v2, const char *nomefile, int linea);
#else
	void  _assegna1_v_segmv_d(VETTOREd *v1, int st, int end, const VETTOREd *v2);
#endif

//! assegna agli elementi di un vettore uguali ad un determinato valore un altro valore.
/*!
	\ingroup setvett

	\param v1 il vettore di riferimento
	\param v2 il vettore da cui copiare i valori
	\param val1 il valore da confrontare
	\param val2 il valore da assegnare

	\test v1[which(v2==val1)] <- val2
*/
#ifdef MDEBUG
	#define assegna1_v_indxeq_d(v1, v2, val1, val2) _assegna1_v_indxeq_d(v1, v2, val1, val2, __FILE__, __LINE__)
#else
	#define assegna1_v_indxeq_d(v1, v2, val1, val2) _assegna1_v_indxeq_d(v1, v2, val1, val2)
#endif
#ifdef MDEBUG
	void  _assegna1_v_indxeq_d(VETTOREd *v1, VETTOREd *v2, double val1, double val2, const char *nomefile, int linea);
#else
	void  _assegna1_v_indxeq_d(VETTOREd *v1, VETTOREd *v2, double val1, double val2);
#endif

//! assegna agli elementi di un vettore uguali o diversi da NA un altro valore.
/*!
	\ingroup setvett

	\param v1 il vettore di riferimento
	\param v2 il vettore da cui copiare i valori
	\param val il valore da assegnare
	\param complemento se vero assegna agli elementi diversi da NA

	\test v1[which(is.na(v2))] <- val oppure v1[which(!is.na(v2))] <- val
*/
#ifdef MDEBUG
	#define assegna1_v_indxNA_d(v1, v2, val, complemento) _assegna1_v_indxNA_d(v1, v2, val, complemento, __FILE__, __LINE__)
#else
	#define assegna1_v_indxNA_d(v1, v2, val, complemento) _assegna1_v_indxNA_d(v1, v2, val, complemento)
#endif
#ifdef MDEBUG
	void  _assegna1_v_indxNA_d(VETTOREd *v1, VETTOREd *v2, double val, bool complemento, const char *nomefile, int linea);
#else
	void  _assegna1_v_indxNA_d(VETTOREd *v1, VETTOREd *v2, double val, bool complemento);
#endif

//! incrementa (o decrementa) gli elementi di un vettore di tipo 'double' di \c n in base ad un vettore di indici
/*!
	\ingroup setvett

	\param v il vettore di riferimento
	\param indx il vettore degli indici
	\param s il valore da addizionare

	\test v[indx] <- v[indx] + s
*/
#ifdef MDEBUG
	#define incr1_v_indx_d(v, indx, s) _incr1_v_indx_d(v, indx, s, __FILE__, __LINE__)
#else
	#define incr1_v_indx_d(v, indx, s) _incr1_v_indx_d(v, indx, s)
#endif
#ifdef MDEBUG
	void  _incr1_v_indx_d(VETTOREd *v, const VETTOREi * indx, double s, const char *nomefile, int linea);
#else
	void  _incr1_v_indx_d(VETTOREd *v, const VETTOREi * indx, double s);
#endif

//! promuove un vettore a 'double' (non accade nulla se già 'double').
/*!
	\ingroup opvett

	\param ris il vettore per il risultato
	\param v il vettore da promuovere

	\return il nuovo vettore

	\test implicita
*/
#ifdef MDEBUG
	#define promuovi_d(ris, v) _promuovi_d(ris, v, #ris, __FILE__, __LINE__)
#else
	#define promuovi_d(ris, v) _promuovi_d(ris, v)
#endif
#ifdef MDEBUG
	VETTOREd * _promuovi_d(VETTOREd *ris, const VETTOREd *v, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _promuovi_d(VETTOREd *ris, const VETTOREd *v);
#endif

//! calcola una particolare funzione vettoriale usata in connectivity_scalefree.
/*!
	\ingroup aus

	\param ris un vettore per il risultato
	\param a  parametro della funzione
	\param b  parametro della funzione
	\param c  parametro della funzione
	\param da indice di partenza
	\param a1 indice di arrivo
	\param sgn1 segno della prima operazione di arrivo
	\param sgn2 segno della seconda operazione

	\test  ris <- a[da:a1] + b[da:a1] * sgn1 + c[da:a1] * sgn2
*/
#ifdef MDEBUG
	#define f_aux1_d(ris, a, b, c, da, a1, sgn1, sgn2) _f_aux1_d(ris, a, b, c, da, a1, sgn1, sgn2, #ris, __FILE__, __LINE__)
#else
	#define f_aux1_d(ris, a, b, c, da, a1, sgn1, sgn2) _f_aux1_d(ris, a, b, c, da, a1, sgn1, sgn2)
#endif
#ifdef MDEBUG
	VETTOREd * _f_aux1_d(VETTOREd *ris, const VETTOREd *a, const VETTOREd *b, const VETTOREd *c, int da, int a1, int sgn1, int sgn2, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _f_aux1_d(VETTOREd *ris, const VETTOREd *a, const VETTOREd *b, const VETTOREd *c, int da, int a1, int sgn1, int sgn2);
#endif

//! restituisce l'indice del primo elemento minimo saltando gli NA.
/*!
	\ingroup wvett

	\param v il vettore in cui cercare

	\return il'indice dell'elemento minimo

	\test ris <- which.min(v)
*/
#ifdef MDEBUG
	#define which_v_indxmin_d(v) _which_v_indxmin_d(v, __FILE__, __LINE__)
#else
	#define which_v_indxmin_d(v) _which_v_indxmin_d(v)
#endif
#ifdef MDEBUG
	int  _which_v_indxmin_d(VETTOREd *v, const char *nomefile, int linea);
#else
	int  _which_v_indxmin_d(VETTOREd *v);
#endif

//! assegna un valore agli elementi di una colonna di una matrice di tipo 'double' corrispondenti agli indici contenuti in un vettore.
/*!
	\ingroup wmatr

	\param m la matrice da assegnare
	\param indx il vettore degli indici
	\param colonna la colonna da assegnare
	\param val il valore da assegnare

	\test m[indx, col] <- val
*/
#ifdef MDEBUG
	#define assegna1_ms_colindx_d(m, indx, colonna, val) _assegna1_ms_colindx_d(m, indx, colonna, val, __FILE__, __LINE__)
#else
	#define assegna1_ms_colindx_d(m, indx, colonna, val) _assegna1_ms_colindx_d(m, indx, colonna, val)
#endif
#ifdef MDEBUG
	void  _assegna1_ms_colindx_d(MATRICEd *m, const VETTOREi *indx, int colonna, double val, const char *nomefile, int linea);
#else
	void  _assegna1_ms_colindx_d(MATRICEd *m, const VETTOREi *indx, int colonna, double val);
#endif

//! calcola una particolare funzione vettoriale usata in score_sf.
/*!
	\ingroup aus

	\param ris un vettore per il risultato
	\param a  parametro della funzione
	\param b  parametro della funzione
	\param c  parametro della funzione

	\return il vettore degli indici

	\test  ris <- which( ((sign(a-c)-sign(b-c))!=0)&(sign(b-c)!=0)&(sign(a-c)!=0) )
*/
#ifdef MDEBUG
	#define f_aux2_d(ris, a, b, c) _f_aux2_d(ris, a, b, c, #ris, __FILE__, __LINE__)
#else
	#define f_aux2_d(ris, a, b, c) _f_aux2_d(ris, a, b, c)
#endif
#ifdef MDEBUG
	VETTOREi * _f_aux2_d(VETTOREi *ris, const VETTOREd *a, const VETTOREd *b, const VETTOREd *c, const char *nome, const char *nomefile, int linea);
#else
	VETTOREi * _f_aux2_d(VETTOREi *ris, const VETTOREd *a, const VETTOREd *b, const VETTOREd *c);
#endif

//! calcola una particolare funzione vettoriale usata in score_sf.
/*!
	\ingroup aus

	\param ris un vettore per il risultato
	\param a  parametro della funzione
	\param b  parametro della funzione
	\param c  parametro della funzione

	\return il vettore degli indici

	\test  ris <- which( ((a - b) > 0) & (a < c) )
*/
#ifdef MDEBUG
	#define f_aux3_d(ris, a, b, c) _f_aux3_d(ris, a, b, c, #ris, __FILE__, __LINE__)
#else
	#define f_aux3_d(ris, a, b, c) _f_aux3_d(ris, a, b, c)
#endif
#ifdef MDEBUG
	VETTOREi * _f_aux3_d(VETTOREi *ris, const VETTOREd *a, const VETTOREd *b, const VETTOREd *c, const char *nome, const char *nomefile, int linea);
#else
	VETTOREi * _f_aux3_d(VETTOREi *ris, const VETTOREd *a, const VETTOREd *b, const VETTOREd *c);
#endif

//! unisce due vettori di tipo 'double' ed elimina gli elementi ripetuti
/*!
	\ingroup opvett

	\param ris il vettore per il risultato
	\param v1 il primo vettore
	\param v2 il secondo vettore

	\return il vettore risultante

	\test ris <- union(v1, v2)
*/
#ifdef MDEBUG
	#define unione_d(ris, v1, v2) _unione_d(ris, v1, v2, #ris, __FILE__, __LINE__)
#else
	#define unione_d(ris, v1, v2) _unione_d(ris, v1, v2)
#endif
#ifdef MDEBUG
	VETTOREd * _unione_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _unione_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2);
#endif

//! restituisce un vettore di tipo 'double' corrispondente agli indici degli elementi diversi da val nella matrice m alla riga specificata da un vettore di indici.
/*!
	\ingroup wvett

	\param ris un vettore di tipo 'int' per il risultato
	\param m la matrice di riferimento
	\param v il vettore per gli indici delle righe
	\param val il valore di riferimento

	\test ris <- which(m[v, ] == val)
*/
#ifdef MDEBUG
	#define which_m_indxrowindxeq_d(ris, m, v, val) _which_m_indxrowindxeq_d(ris, m, v, val, #ris, __FILE__, __LINE__)
#else
	#define which_m_indxrowindxeq_d(ris, m, v, val) _which_m_indxrowindxeq_d(ris, m, v, val)
#endif
#ifdef MDEBUG
	VETTOREi * _which_m_indxrowindxeq_d(VETTOREi *ris, const MATRICEd *m, const VETTOREi *v, double val, const char *nome, const char *nomefile, int linea);
#else
	VETTOREi * _which_m_indxrowindxeq_d(VETTOREi *ris, const MATRICEd *m, const VETTOREi *v, double val);
#endif

//! elimina gli elementi doppi da un vettore di tipo 'double'.
/*!
	\ingroup opvett

	\param ris un vettore di tipo 'double' per il risultato
	\param v il vettore di riferimento

	\test ris <- union(v, v)
*/
#ifdef MDEBUG
	#define elimina_doppi_d(ris, v) _elimina_doppi_d(ris, v, #ris, __FILE__, __LINE__)
#else
	#define elimina_doppi_d(ris, v) _elimina_doppi_d(ris, v)
#endif
#ifdef MDEBUG
	VETTOREd * _elimina_doppi_d(VETTOREd *ris, const VETTOREd *v, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _elimina_doppi_d(VETTOREd *ris, const VETTOREd *v);
#endif

//! calcola una particolare funzione scalare usata in dinamica.
/*!
	\ingroup aus

	\param times  parametro della funzione
	\param res  parametro della funzione

	\test  ris <- max(abs(times / res - round(times / res)))
*/
#ifdef MDEBUG
	#define f_aux4_d(times, res) _f_aux4_d(times, res, __FILE__, __LINE__)
#else
	#define f_aux4_d(times, res) _f_aux4_d(times, res)
#endif
#ifdef MDEBUG
	double  _f_aux4_d(const VETTOREd *times, double res, const char *nomefile, int linea);
#else
	double  _f_aux4_d(const VETTOREd *times, double res);
#endif

//! assegna ad una matrice una colonna copiandola da un vettore
/*!
	\ingroup setmatr

	\param m la matrice di riferimento
	\param colonna la colonna
	\param v il vettore da assegnare

	\test m[,colonna] <- v
*/
#ifdef MDEBUG
	#define assegna1_mv_colonna_d(m, colonna, v) _assegna1_mv_colonna_d(m, colonna, v, __FILE__, __LINE__)
#else
	#define assegna1_mv_colonna_d(m, colonna, v) _assegna1_mv_colonna_d(m, colonna, v)
#endif
#ifdef MDEBUG
	void  _assegna1_mv_colonna_d(MATRICEd *m, int colonna, const VETTOREd *v, const char *nomefile, int linea);
#else
	void  _assegna1_mv_colonna_d(MATRICEd *m, int colonna, const VETTOREd *v);
#endif

//! calcola una particolare funzione vettoriale usata in dinamica.
/*!
	\ingroup aus

	\param ris un vettore per il risultato
	\param alpha  parametro della funzione
	\param targ  parametro della funzione
	\param theta  parametro della funzione
	\param xmin  parametro della funzione

	\test  ris <- 1 / (1 + exp(-alpha * (targ - theta))) * (1 - xmin) + xmin
*/
#ifdef MDEBUG
	#define f_aux5_d(ris, alpha, targ, theta, xmin) _f_aux5_d(ris, alpha, targ, theta, xmin, #ris, __FILE__, __LINE__)
#else
	#define f_aux5_d(ris, alpha, targ, theta, xmin) _f_aux5_d(ris, alpha, targ, theta, xmin)
#endif
#ifdef MDEBUG
	VETTOREd * _f_aux5_d(VETTOREd *ris, const VETTOREd *alpha, const VETTOREd *targ, const VETTOREd *theta, const VETTOREd *xmin, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _f_aux5_d(VETTOREd *ris, const VETTOREd *alpha, const VETTOREd *targ, const VETTOREd *theta, const VETTOREd *xmin);
#endif

//! calcola una particolare funzione vettoriale usata in dinamica.
/*!
	\ingroup aus

	\param ris un vettore per il risultato
	\param res  parametro della funzione
	\param k parametro della funzione
	\param targetT  parametro della funzione
	\param n parametro della funzione

	\test  ris <- matrix(res * k * (targetT - n), ncol=1)
*/
#ifdef MDEBUG
	#define f_aux6_d(ris, res, k, targetT, n) _f_aux6_d(ris, res, k, targetT, n, #ris, __FILE__, __LINE__)
#else
	#define f_aux6_d(ris, res, k, targetT, n) _f_aux6_d(ris, res, k, targetT, n)
#endif
#ifdef MDEBUG
	VETTOREd * _f_aux6_d(VETTOREd *ris, double res, const VETTOREd *k, const VETTOREd *targetT, const VETTOREd *n, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _f_aux6_d(VETTOREd *ris, double res, const VETTOREd *k, const VETTOREd *targetT, const VETTOREd *n);
#endif

//! somma ad un vettore di tipo 'double' un altro vettore di tipo 'double' elemento per elemento.
/*!
	\ingroup opvett

	\param v1 il primo vettore
	\param v2 il secondo vettore

	\test  v1 <- v1 + v2
*/
#ifdef MDEBUG
	#define somma1_vv_d(v1, v2) _somma1_vv_d(v1, v2, __FILE__, __LINE__)
#else
	#define somma1_vv_d(v1, v2) _somma1_vv_d(v1, v2)
#endif
#ifdef MDEBUG
	void  _somma1_vv_d(VETTOREd *v1, const VETTOREd *v2, const char *nomefile, int linea);
#else
	void  _somma1_vv_d(VETTOREd *v1, const VETTOREd *v2);
#endif

//! aggiunge ad una matrice una colonna copiandola da un vettore
/*!
	\ingroup setmatr

	\param ris la matrice di riferimento
	\param colonna la colonna di riferimento
	\param v il vettore da assegnare

	\return la matrice con la nuova riga

	\note il parametro \c colonna dev'essere minore o uguale ad 1 + il numero di colonna della matrice

	\test ris[,colonna] <- v # con colonna <= dim(ris)[1] + 1
*/
#ifdef MDEBUG
	#define aggiungi_mv_colonna_d(ris, colonna, v) _aggiungi_mv_colonna_d(ris, colonna, v, #ris, __FILE__, __LINE__)
#else
	#define aggiungi_mv_colonna_d(ris, colonna, v) _aggiungi_mv_colonna_d(ris, colonna, v)
#endif
#ifdef MDEBUG
	MATRICEd * _aggiungi_mv_colonna_d(MATRICEd *ris, int colonna, const VETTOREd *v, const char *nome, const char *nomefile, int linea);
#else
	MATRICEd * _aggiungi_mv_colonna_d(MATRICEd *ris, int colonna, const VETTOREd *v);
#endif

//! calcola una particolare funzione scalare usata in dinamica.
/*!
	\ingroup aus

	\param vettore  parametro della funzione
	\param scalare  parametro della funzione

	\test  ris <- which.min(abs(vettore - scalare))
*/
#ifdef MDEBUG
	#define f_aux7_d(vettore, scalare) _f_aux7_d(vettore, scalare, __FILE__, __LINE__)
#else
	#define f_aux7_d(vettore, scalare) _f_aux7_d(vettore, scalare)
#endif
#ifdef MDEBUG
	int  _f_aux7_d(const VETTOREd *vettore, double scalare, const char *nomefile, int linea);
#else
	int  _f_aux7_d(const VETTOREd *vettore, double scalare);
#endif

//! calcola una particolare funzione vettoriale usata in dinamica.
/*!
	\ingroup aus

	\param ris un vettore per il risultato
	\param a  il primo vettore
	\param b  il secondo vettore

	\test  ris <- a * (1 - b) + b
*/
#ifdef MDEBUG
	#define f_aux8_d(ris, a, b) _f_aux8_d(ris, a, b, #ris, __FILE__, __LINE__)
#else
	#define f_aux8_d(ris, a, b) _f_aux8_d(ris, a, b)
#endif
#ifdef MDEBUG
	VETTOREd * _f_aux8_d(VETTOREd *ris, const VETTOREd *a, const VETTOREd *b, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _f_aux8_d(VETTOREd *ris, const VETTOREd *a, const VETTOREd *b);
#endif

//! crea una nuova matrice selezionandone gli indici da un vettore.
/*!
	\ingroup wmatr

	\param ris un matrice di tipo 'double' per il risultato
	\param m la matrice di riferimento
	\param indici il vettore degli indici

	\test ris <- matrix(m, indici)
*/
#ifdef MDEBUG
	#define seleziona_colonne_d(ris, m, indici) _seleziona_colonne_d(ris, m, indici, #ris, __FILE__, __LINE__)
#else
	#define seleziona_colonne_d(ris, m, indici) _seleziona_colonne_d(ris, m, indici)
#endif
#ifdef MDEBUG
	MATRICEd * _seleziona_colonne_d(MATRICEd *ris, const MATRICEd *m, const VETTOREi *indici, const char *nome, const char *nomefile, int linea);
#else
	MATRICEd * _seleziona_colonne_d(MATRICEd *ris, const MATRICEd *m, const VETTOREi *indici);
#endif

//! cambia le dimensioni di una matrice  di tipo 'double' mantenendone il prodotto.
/*!
	\ingroup oggetti

	\param m la matrice di riferimento
	\param nr il nuovo numero di righe (-1 per farlo derivare da \c nc)
	\param nc il nuovo numero di colonne (-1 per farlo derivare da \c nr)

	\test m <- matrix(m, nr=nr)
*/
#ifdef MDEBUG
	#define cambiadim1_d(m, nr, nc) _cambiadim1_d(m, nr, nc, __FILE__, __LINE__)
#else
	#define cambiadim1_d(m, nr, nc) _cambiadim1_d(m, nr, nc)
#endif
#ifdef MDEBUG
	void  _cambiadim1_d(MATRICEd *m, const int nr, const int nc, const char *nomefile, int linea);
#else
	void  _cambiadim1_d(MATRICEd *m, const int nr, const int nc);
#endif

//! moltiplica gli elementi corrispondenti di due matrici di tipo 'double'.
/*!
	\ingroup opmatr

	\param ris una matrice per il risultato
	\param m1 la prima matrice
	\param m2 la seconda matrice

	\test ris <- m1 * m2
*/
#ifdef MDEBUG
	#define moltiplica_mm_d(ris, m1, m2) _moltiplica_mm_d(ris, m1, m2, #ris, __FILE__, __LINE__)
#else
	#define moltiplica_mm_d(ris, m1, m2) _moltiplica_mm_d(ris, m1, m2)
#endif
#ifdef MDEBUG
	MATRICEd * _moltiplica_mm_d(MATRICEd *ris, const MATRICEd *m1, const MATRICEd *m2, const char *nome, const char *nomefile, int linea);
#else
	MATRICEd * _moltiplica_mm_d(MATRICEd *ris, const MATRICEd *m1, const MATRICEd *m2);
#endif

//! moltiplica riga per riga gli elementi corrispondenti di una matrice di tipo 'double' per un vettore.
/*!
	\ingroup opmatr

	\param m la matrice
	\param v il vettore

	\test ris <- m * v
*/
#ifdef MDEBUG
	#define moltiplica1_mv_d(m, v) _moltiplica1_mv_d(m, v, __FILE__, __LINE__)
#else
	#define moltiplica1_mv_d(m, v) _moltiplica1_mv_d(m, v)
#endif
#ifdef MDEBUG
	void  _moltiplica1_mv_d(MATRICEd *m, const VETTOREd *v, const char *nomefile, int linea);
#else
	void  _moltiplica1_mv_d(MATRICEd *m, const VETTOREd *v);
#endif

//! divide un vettore di tipo 'double' per un altro vettore di tipo 'double' elemento per elemento.
/*!
	\ingroup opvett

	\param ris un vettore per il risultato
	\param v1 il primo vettore
	\param v2 il secondo vettore

	\test  ris <- v1 / v2
*/
#ifdef MDEBUG
	#define dividi_vv_d(ris, v1, v2) _dividi_vv_d(ris, v1, v2, #ris, __FILE__, __LINE__)
#else
	#define dividi_vv_d(ris, v1, v2) _dividi_vv_d(ris, v1, v2)
#endif
#ifdef MDEBUG
	VETTOREd * _dividi_vv_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _dividi_vv_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2);
#endif

//! calcola il segno degli elementi di una matrice di tipo 'double'.
/*!
	\ingroup opmatr

	\param ris la matrice risultante
	\param m la matrice d'ingresso

	\test ris <- sign(m)
*/
#ifdef MDEBUG
	#define segno_m_d(ris, m) _segno_m_d(ris, m, #ris, __FILE__, __LINE__)
#else
	#define segno_m_d(ris, m) _segno_m_d(ris, m)
#endif
#ifdef MDEBUG
	MATRICEi * _segno_m_d(MATRICEi *ris, const MATRICEd *m, const char *nome, const char *nomefile, int linea);
#else
	MATRICEi * _segno_m_d(MATRICEi *ris, const MATRICEd *m);
#endif

//! assegna un valore agli elementi di una colonna di una matrice di tipo 'double'.
/*!
	\ingroup setmatr

	\param m la matrice da assegnare
	\param colonna la colonna da assegnare
	\param val il valore da assegnare

	\test m[, colonna] <- val
*/
#ifdef MDEBUG
	#define assegna1_ms_colonna_d(m, colonna, val) _assegna1_ms_colonna_d(m, colonna, val, __FILE__, __LINE__)
#else
	#define assegna1_ms_colonna_d(m, colonna, val) _assegna1_ms_colonna_d(m, colonna, val)
#endif
#ifdef MDEBUG
	void  _assegna1_ms_colonna_d(MATRICEd *m, int colonna, double val, const char *nomefile, int linea);
#else
	void  _assegna1_ms_colonna_d(MATRICEd *m, int colonna, double val);
#endif

//! divide un vettore di tipo 'double' per un altro vettore di tipo 'double' elemento per elemento.
/*!
	\ingroup opvett

	\param v1 il primo vettore
	\param v2 il secondo vettore

	\test  v1 <- v1 / v2
*/
#ifdef MDEBUG
	#define dividi1_vv_d(v1, v2) _dividi1_vv_d(v1, v2, __FILE__, __LINE__)
#else
	#define dividi1_vv_d(v1, v2) _dividi1_vv_d(v1, v2)
#endif
#ifdef MDEBUG
	void  _dividi1_vv_d(VETTOREd *v1, const VETTOREd *v2, const char *nomefile, int linea);
#else
	void  _dividi1_vv_d(VETTOREd *v1, const VETTOREd *v2);
#endif

//! calcola il valore assoluto degli elementi di un vettore di tipo 'double' "sul posto".
/*!
	\ingroup opvett

	\param v il vettore di riferimento

	\test v <- abs(v)
*/
#ifdef MDEBUG
	#define abs1_v_d(v) _abs1_v_d(v, __FILE__, __LINE__)
#else
	#define abs1_v_d(v) _abs1_v_d(v)
#endif
#ifdef MDEBUG
	void  _abs1_v_d(VETTOREd *v, const char *nomefile, int linea);
#else
	void  _abs1_v_d(VETTOREd *v);
#endif

//! converte una matrice di tipo 'double' in una matrice di tipo 'int' (la versione intera non fa nulla).
/*!
	\ingroup opmatr

	\param ris la matrice da restituire
	\param m la matrice di riferimento

	\return la matrice arrotondata

	\test implicita
*/
#ifdef MDEBUG
	#define arrotonda_m_d(ris, m) _arrotonda_m_d(ris, m, #ris, __FILE__, __LINE__)
#else
	#define arrotonda_m_d(ris, m) _arrotonda_m_d(ris, m)
#endif
#ifdef MDEBUG
	MATRICEi * _arrotonda_m_d(MATRICEi *ris, const MATRICEd *m, const char *nome, const char *nomefile, int linea);
#else
	MATRICEi * _arrotonda_m_d(MATRICEi *ris, const MATRICEd *m);
#endif

//! elimina un elemento da una matrice di tipo 'double' tramite un indice.
/*!
	\ingroup opmatr

	\param m la matrice che contiene l'elemento da eliminare (non e' costante perche' si potrebbe restringere)
	\param riga l'indice

	\test m <- m[-riga,]
*/
#ifdef MDEBUG
	#define elimina1_riga_d(m, riga) _elimina1_riga_d(m, riga, __FILE__, __LINE__)
#else
	#define elimina1_riga_d(m, riga) _elimina1_riga_d(m, riga)
#endif
#ifdef MDEBUG
	void  _elimina1_riga_d(MATRICEd *m, int riga, const char *nomefile, int linea);
#else
	void  _elimina1_riga_d(MATRICEd *m, int riga);
#endif

//! estrae una colonna da una matrice di tipo 'double'.
/*!
	\ingroup opmatr

	\param ris un vettore per il risultato
	\param m la matrice di riferimento
	\param colonna la colonna di riferimento

	\test ris <- M[,colonna]
*/
#ifdef MDEBUG
	#define colonna_d(ris, m, colonna) _colonna_d(ris, m, colonna, #ris, __FILE__, __LINE__)
#else
	#define colonna_d(ris, m, colonna) _colonna_d(ris, m, colonna)
#endif
#ifdef MDEBUG
	VETTOREd * _colonna_d(VETTOREd *ris, const MATRICEd *m, int colonna, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _colonna_d(VETTOREd *ris, const MATRICEd *m, int colonna);
#endif

//! restituisce un vettore di tipo 'double' corrispondente agli indici degli elementi non simmetrici nella matrice m.
/*!
	\ingroup wmatr

	\param ris un vettore di tipo 'int' per il risultato
	\param m la matrice di riferimento

	\test ris <- which((m == t(m)) == 0) # cioe` which(m != t(m))
*/
#ifdef MDEBUG
	#define which_m_indxnsimm_d(ris, m) _which_m_indxnsimm_d(ris, m, #ris, __FILE__, __LINE__)
#else
	#define which_m_indxnsimm_d(ris, m) _which_m_indxnsimm_d(ris, m)
#endif
#ifdef MDEBUG
	VETTOREi * _which_m_indxnsimm_d(VETTOREi *ris, const MATRICEd *m, const char *nome, const char *nomefile, int linea);
#else
	VETTOREi * _which_m_indxnsimm_d(VETTOREi *ris, const MATRICEd *m);
#endif

//! assegna agli elementi di una matrice uguali o diversi da NA un altro valore.
/*!
	\ingroup setmatr

	\param m1 la matrice di riferimento
	\param m2 la matrice da cui copiare
	\param val il valore da assegnare
	\param complemento se vero assegna agli elementi diversi da NA

	\test m1[which(is.na(m2))] <- val oppure m1[which(!is.na(m2))] <- val
*/
#ifdef MDEBUG
	#define assegna1_m_indxNA_d(m1, m2, val, complemento) _assegna1_m_indxNA_d(m1, m2, val, complemento, __FILE__, __LINE__)
#else
	#define assegna1_m_indxNA_d(m1, m2, val, complemento) _assegna1_m_indxNA_d(m1, m2, val, complemento)
#endif
#ifdef MDEBUG
	void  _assegna1_m_indxNA_d(MATRICEd *m1, MATRICEd *m2, double val, bool complemento, const char *nomefile, int linea);
#else
	void  _assegna1_m_indxNA_d(MATRICEd *m1, MATRICEd *m2, double val, bool complemento);
#endif

//! assegna un vettore agli elementi di una matrice di tipo 'double'.
/*!
	\ingroup opmatr

	\param m la matrice da assegnare
	\param v il vettore da assegnare

	\test m <- v
*/
#ifdef MDEBUG
	#define assegna1_mv_d(m, v) _assegna1_mv_d(m, v, __FILE__, __LINE__)
#else
	#define assegna1_mv_d(m, v) _assegna1_mv_d(m, v)
#endif
#ifdef MDEBUG
	void  _assegna1_mv_d(MATRICEd *m, const VETTOREd *v, const char *nomefile, int linea);
#else
	void  _assegna1_mv_d(MATRICEd *m, const VETTOREd *v);
#endif

//! calcola una particolare funzione vettoriale usata in HMM.
/*!
	\ingroup aus1

	\param ris il primo vettore
	\param val1 il primo
	\param v2 il secondo vettore
	\param indx il vettore degli indici
	\param val2 il secondo numero

	\test  v1[indx] <- (val1 - abs(val2 - v2[indx])) + v1[indx]
*/
#ifdef MDEBUG
	#define f_aux9_d(ris, val1, v2, indx, val2) _f_aux9_d(ris, val1, v2, indx, val2, #ris, __FILE__, __LINE__)
#else
	#define f_aux9_d(ris, val1, v2, indx, val2) _f_aux9_d(ris, val1, v2, indx, val2)
#endif
#ifdef MDEBUG
	VETTOREd * _f_aux9_d(VETTOREd *ris, double val1, VETTOREi *v2, VETTOREi *indx, double val2, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _f_aux9_d(VETTOREd *ris, double val1, VETTOREi *v2, VETTOREi *indx, double val2);
#endif

//! calcola il segno degli elementi di un vettore di tipo 'double'.
/*!
	\ingroup opvett

	\param ris il vettore risultante
	\param v il vettore d'ingresso

	\test ris <- sign(v)
*/
#ifdef MDEBUG
	#define segno_v_d(ris, v) _segno_v_d(ris, v, #ris, __FILE__, __LINE__)
#else
	#define segno_v_d(ris, v) _segno_v_d(ris, v)
#endif
#ifdef MDEBUG
	VETTOREi * _segno_v_d(VETTOREi *ris, const VETTOREd *v, const char *nome, const char *nomefile, int linea);
#else
	VETTOREi * _segno_v_d(VETTOREi *ris, const VETTOREd *v);
#endif

//! assegna agli elementi di una matrice minori di un determinato valore un altro valore.
/*!
	\ingroup setmatr

	\param m1 la matrice di riferimento
	\param m2 la matrice da cui copiare i valori
	\param val1 il valore da confrontare
	\param val2 il valore da assegnare

	\test m1[which(m2<val1)] <- val2
*/
#ifdef MDEBUG
	#define assegna1_m_indxle_d(m1, m2, val1, val2) _assegna1_m_indxle_d(m1, m2, val1, val2, __FILE__, __LINE__)
#else
	#define assegna1_m_indxle_d(m1, m2, val1, val2) _assegna1_m_indxle_d(m1, m2, val1, val2)
#endif
#ifdef MDEBUG
	void  _assegna1_m_indxle_d(MATRICEd *m1, MATRICEd *m2, double val1, double val2, const char *nomefile, int linea);
#else
	void  _assegna1_m_indxle_d(MATRICEd *m1, MATRICEd *m2, double val1, double val2);
#endif

//! assegna agli elementi di un vettore maggiori di un determinato valore un altro valore.
/*!
	\ingroup setvett

	\param v1 il vettore di riferimento
	\param v2 il vettore da cui copiare i valori
	\param val1 il valore da confrontare
	\param val2 il valore da assegnare

	\test v1[which(v2>val1)] <- val2
*/
#ifdef MDEBUG
	#define assegna1_v_indxgt_d(v1, v2, val1, val2) _assegna1_v_indxgt_d(v1, v2, val1, val2, __FILE__, __LINE__)
#else
	#define assegna1_v_indxgt_d(v1, v2, val1, val2) _assegna1_v_indxgt_d(v1, v2, val1, val2)
#endif
#ifdef MDEBUG
	void  _assegna1_v_indxgt_d(VETTOREd *v1, VETTOREd *v2, double val1, double val2, const char *nomefile, int linea);
#else
	void  _assegna1_v_indxgt_d(VETTOREd *v1, VETTOREd *v2, double val1, double val2);
#endif

//! moltiplica due vettori di tipo 'double' elemento per elemento.
/*!
	\ingroup opvett

	\param ris un vettore per il risultato
	\param v1 il primo vettore da moltiplicare
	\param v2 il secondo vettore da moltiplicare

	\test ris <- v1 * v2
*/
#ifdef MDEBUG
	#define moltiplica_vv_d(ris, v1, v2) _moltiplica_vv_d(ris, v1, v2, #ris, __FILE__, __LINE__)
#else
	#define moltiplica_vv_d(ris, v1, v2) _moltiplica_vv_d(ris, v1, v2)
#endif
#ifdef MDEBUG
	VETTOREd * _moltiplica_vv_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _moltiplica_vv_d(VETTOREd *ris, const VETTOREd *v1, const VETTOREd *v2);
#endif

//! assegna agli elementi di un vettore minori o uguali di un determinato valore un altro valore.
/*!
	\ingroup setvett

	\param v1 il vettore di riferimento
	\param v2 il vettore da cui copiare i valori
	\param val1 il valore da confrontare
	\param val2 il valore da assegnare

	\test v1[which(v2<=val1)] <- val2
*/
#ifdef MDEBUG
	#define assegna1_v_indxle_d(v1, v2, val1, val2) _assegna1_v_indxle_d(v1, v2, val1, val2, __FILE__, __LINE__)
#else
	#define assegna1_v_indxle_d(v1, v2, val1, val2) _assegna1_v_indxle_d(v1, v2, val1, val2)
#endif
#ifdef MDEBUG
	void  _assegna1_v_indxle_d(VETTOREd *v1, VETTOREd *v2, double val1, double val2, const char *nomefile, int linea);
#else
	void  _assegna1_v_indxle_d(VETTOREd *v1, VETTOREd *v2, double val1, double val2);
#endif

//! calcola una particolare funzione vettoriale usata in scoremodular.
/*!
	\ingroup aus1

	\param ris un vettore per il risultato
	\param indx  parametro della funzione
	\param a  parametro della funzione
	\param b  parametro della funzione
	\param m  parametro della funzione

	\test  ris[indx] <- (sign(a[indx]-b[indx]))*m[indx]/1
*/
#ifdef MDEBUG
	#define f_aux10_d(ris, indx, a, b, m) _f_aux10_d(ris, indx, a, b, m, #ris, __FILE__, __LINE__)
#else
	#define f_aux10_d(ris, indx, a, b, m) _f_aux10_d(ris, indx, a, b, m)
#endif
#ifdef MDEBUG
	VETTOREd * _f_aux10_d(VETTOREd *ris, const VETTOREi *indx, const VETTOREd *a, const VETTOREd *b, const VETTOREd *m, const char *nome, const char *nomefile, int linea);
#else
	VETTOREd * _f_aux10_d(VETTOREd *ris, const VETTOREi *indx, const VETTOREd *a, const VETTOREd *b, const VETTOREd *m);
#endif

#endif
