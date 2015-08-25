#ifndef PARS_H
#define PARS_H

#define _MATH_DEFINES_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <glib.h>

/* malloc sarebbe pi? veloce, ma ? meno sicuro */
#define alloc(A, B) (B *) malloc(A * sizeof(B))
/* #define alloc(A, B) Calloc(A, B) */
#define libera(A) if (A != NULL) free(A)

#define N_F 21
#define N_V 3
#define N_C 2
#define MN_F N_F + 20
#define MN_V N_V + 7
#define MN_C N_C + 8
#define M_PI 3.14159265358979323846

enum TipoElem { ERRORE, COSTANTE, VARIABILE, FUNZ };

enum TipoOF { PREFISSO, INFISSO, POSTFISSO, FUNZIONE };

enum Stato {
   INIZIO,
   OPERANDO,
   OP_INFISSO,
   OP_POSTFISSO,
   PAR_A,
   PAR_C,
   OP_PREFISSO,
   OP_FUNZIONE,
   SEP
};

enum OpFunzione {
   OP_NEG, OP_SOMMA, OP_SOTTR,
   OP_MOLT, OP_DIV, OP_POT,
   OP_ABS, OP_SEN, OP_COS,
   OP_LOG, OP_ARCTAN, OP_TAN,
   OP_MOD, OP_SE, OP_ASEN,
   OP_ACOS, OP_FATT, OP_UGUALE,
   OP_NOT, OP_MINORE, OP_MAGG,
   OP_MINORE_U, OP_MAGG_U, OP_DIVERSO,
   OP_SENH, OP_COSH, OP_TANH,
   OP_ASENH, OP_ACOSH, OP_ATANH,
   OP_RADQ, OP_ANGOLO, OP_LN,
   OP_SEGNO, OP_APPROX, OP_TRONCA
};

struct Elem {
   GString *nome;
   double val; /* deve rimanere double perche´ e` usato dagli operandi */
   int pr;
   int paren;
   enum TipoOF tipo;
   int n_arg;
   struct Elem *def;
   unsigned len;
};

GQueue *p;

GHashTable *variabili;
GHashTable *costanti;
GHashTable *funzioni;

struct Elem f[MN_F];
struct Elem v[MN_V];
struct Elem c[MN_C];

int nf, nv, nc;
int stato[9];

void Ripristina();
void Cancella();
struct Elem *FCV(const char *nome, int tipo);
double ValV(const char *nome);
double ValC(const char *nome);
void DefC(const char *nome, double val);
void DefV(const char *nome);
GString *ModificaPF(const char *nome, int lev);
GString *Informa();
GString *DefF(const char *nome, const char *espr, int args, int lev);
double Fatt(int n);
GString *Calcola(const char *nome, double *a, double *ris);

#endif
