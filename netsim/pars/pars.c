#include "pars.h"


void Ripristina()
{
	struct Elem *tmp_c;
	struct Elem *tmp_v;
	struct Elem *tmp_f;

	nf = 33;
	nv = 3;
	nc = 2;

#define INIZIO1 (1 << (INIZIO + 1))
#define OPERANDO1 (1 << (OPERANDO + 1))
#define OP_INFISSO1 (1 << (OP_INFISSO + 1))
#define OP_POSTFISSO1 (1 << (OP_POSTFISSO + 1))
#define PAR_A1 (1 << (PAR_A + 1))
#define PAR_C1 (1 << (PAR_C + 1))
#define OP_PREFISSO1 (1 << (OP_PREFISSO + 1))
#define OP_FUNZIONE1 (1 << (OP_FUNZIONE + 1))
#define SEP1 (1 << (SEP + 1))

	stato[INIZIO] = INIZIO1 | OPERANDO1 | PAR_A1 | OP_PREFISSO1 | OP_FUNZIONE1;
	stato[OPERANDO] = OP_INFISSO1 | OP_POSTFISSO1 | PAR_C1 | SEP1;
	stato[OP_INFISSO] = OPERANDO1 | PAR_A1 | OP_PREFISSO1 | OP_FUNZIONE1;
	stato[OP_POSTFISSO] = OP_INFISSO1 | OP_POSTFISSO1 | PAR_C1 | SEP1;
	stato[PAR_A] = OPERANDO1 | PAR_A1 | OP_PREFISSO1 | OP_FUNZIONE1;
	stato[PAR_C] = OP_INFISSO1 | OP_POSTFISSO1 | PAR_C1 | SEP1;
	stato[OP_PREFISSO] = OPERANDO1 | PAR_A1 | OP_PREFISSO1 | OP_FUNZIONE1;
	stato[OP_FUNZIONE] = PAR_A1;
	stato[SEP] = OPERANDO1 | PAR_A1 | OP_PREFISSO1 | OP_FUNZIONE1;

	costanti = g_hash_table_new(g_str_hash, g_str_equal);
	tmp_c = alloc(1, struct Elem);
	tmp_c->nome = g_string_new("Pi");
	tmp_c->val = M_PI;
	tmp_c->n_arg = 0; /* le costanti hanno zero argomenti ed un nome */
	g_hash_table_insert(costanti, "Pi", (gpointer) tmp_c);
	/* lo sovrascrivo ma tanto l'ho memorizzato nella ht */
	tmp_c = alloc(1, struct Elem);
	tmp_c->nome = g_string_new("e");
	tmp_c->val = exp(1);
	tmp_c->n_arg = 0; /* le costanti hanno zero argomenti ed un nome */
	g_hash_table_insert(costanti, "e", (gpointer) tmp_c);

	variabili = g_hash_table_new(g_str_hash, g_str_equal);
	tmp_v = alloc(1, struct Elem);
	tmp_v->nome = g_string_new("x");
	tmp_v->paren = 1; /* questo indicizza le variabili */
	tmp_v->n_arg = 0; /* le variabili hanno zero argomenti e saranno senza nome */
	g_hash_table_insert(variabili, "x", (gpointer) tmp_v);
	tmp_v = alloc(1, struct Elem);
	tmp_v->nome = g_string_new("y");
	tmp_v->paren = 2; /* questo indicizza le variabili */
	tmp_v->n_arg = 0; /* le variabili hanno zero argomenti e saranno senza nome */
	g_hash_table_insert(variabili, "y", (gpointer) tmp_v);
	tmp_v = alloc(1, struct Elem);
	tmp_v->nome = g_string_new("z");
	tmp_v->paren = 3; /* questo indicizza le veriabili */
	tmp_v->n_arg = 0; /* le variabili hanno zero argomenti e saranno senza nome */
	g_hash_table_insert(variabili, "z", (gpointer) tmp_v);

	/* le funzioni hanno almeno un argomento e hanno un nome */
	funzioni = g_hash_table_new(g_str_hash, g_str_equal);
	/* questo e` il segno meno relativo, che e` ambiguo; lo sistemero` dopo */
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("?");
	tmp_f->n_arg = 1;
	tmp_f->val = OP_NEG;
	tmp_f->tipo = PREFISSO;
	tmp_f->pr = 1;
	g_hash_table_insert(funzioni, "§", (gpointer) tmp_f);
	tmp_f = malloc(1 * sizeof(struct Elem));
	tmp_f->nome = g_string_new("+");
	tmp_f->n_arg = 2;
	tmp_f->val = OP_SOMMA;
	tmp_f->tipo = INFISSO;
	tmp_f->pr = 3;
	g_hash_table_insert(funzioni, "+", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("-");
	tmp_f->n_arg = 2;
	tmp_f->val = OP_SOTTR;
	tmp_f->tipo = INFISSO;
	tmp_f ->pr = 3;
	g_hash_table_insert(funzioni, "-", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("*");
	tmp_f->n_arg = 2;
	tmp_f->val = OP_MOLT;
	tmp_f->tipo = INFISSO;
	tmp_f->pr = 2;
	g_hash_table_insert(funzioni, "*", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("/");
	tmp_f->n_arg = 2;
	tmp_f->val = OP_DIV;
	tmp_f->tipo = INFISSO;
	tmp_f->pr = 2;
	g_hash_table_insert(funzioni, "/", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("^");
	tmp_f->n_arg = 2;
	tmp_f->val = OP_POT;
	tmp_f->tipo = INFISSO;
	tmp_f->pr = 1;
	g_hash_table_insert(funzioni, "^", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("abs");
	tmp_f->n_arg = 1;
	tmp_f->val = OP_ABS;
	tmp_f->tipo = FUNZIONE;
	tmp_f->pr = 1;
	g_hash_table_insert(funzioni, "abs", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("sin");
	tmp_f->n_arg = 1;
	tmp_f->val = OP_SEN;
	tmp_f->tipo = FUNZIONE;
	tmp_f->pr = 1;
	g_hash_table_insert(funzioni, "sin", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("cos");
	tmp_f->n_arg = 1;
	tmp_f->val = OP_COS;
	tmp_f->tipo = FUNZIONE;
	tmp_f->pr = 1;
	g_hash_table_insert(funzioni, "cos", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("log");
	tmp_f->n_arg = 2;
	tmp_f->val = OP_LOG;
	tmp_f->tipo = FUNZIONE;
	tmp_f->pr = 1;
	g_hash_table_insert(funzioni, "log", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("atan");
	tmp_f->n_arg = 1;
	tmp_f->val = OP_ARCTAN;
	tmp_f->tipo = FUNZIONE;
	tmp_f->pr = 1;
	g_hash_table_insert(funzioni, "atan", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("tan");
	tmp_f->n_arg = 1;
	tmp_f->val = OP_TAN;
	tmp_f->tipo = FUNZIONE;
	tmp_f->pr = 1;
	g_hash_table_insert(funzioni, "tan", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("%");
	tmp_f->n_arg = 2;
	tmp_f->val = OP_MOD;
	tmp_f->tipo = INFISSO;
	tmp_f->pr = 2;
	g_hash_table_insert(funzioni, "%", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("if");
	tmp_f->n_arg = 3;
	tmp_f->val = OP_SE;
	tmp_f->tipo = FUNZIONE;
	tmp_f->pr = 1;
	g_hash_table_insert(funzioni, "if", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("asin");
	tmp_f->n_arg = 1;
	tmp_f->val = OP_ASEN;
	tmp_f->tipo = FUNZIONE;
	tmp_f->pr = 1;
	g_hash_table_insert(funzioni, "asin", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("acos");
	tmp_f->n_arg = 1;
	tmp_f->val = OP_ACOS;
	tmp_f->tipo = FUNZIONE;
	tmp_f->pr = 1;
	g_hash_table_insert(funzioni, "acos", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("!");
	tmp_f->n_arg = 1;
	tmp_f->val = OP_FATT;
	tmp_f->tipo = POSTFISSO;
	tmp_f->pr = 1;
	g_hash_table_insert(funzioni, "!", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("=");
	tmp_f->n_arg = 2;
	tmp_f->val = OP_UGUALE;
	tmp_f->tipo = INFISSO;
	tmp_f->pr = 1;
	g_hash_table_insert(funzioni, "=", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("~");
	tmp_f->n_arg = 1;
	tmp_f->val = OP_NOT;
	tmp_f->tipo = PREFISSO;
	tmp_f->pr = 1;
	g_hash_table_insert(funzioni, "~", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("<");
	tmp_f->n_arg = 2;
	tmp_f->val = OP_MINORE;
	tmp_f->tipo = INFISSO;
	tmp_f->pr = 1;
	g_hash_table_insert(funzioni, "<", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new(">");
	tmp_f->n_arg = 2;
	tmp_f->val = OP_MAGG;
	tmp_f->tipo = INFISSO;
	tmp_f->pr = 1;
	g_hash_table_insert(funzioni, ">", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("<=");
	tmp_f->n_arg = 2;
	tmp_f->val = OP_MINORE_U;
	tmp_f->tipo = INFISSO;
	tmp_f->pr = 1;
	g_hash_table_insert(funzioni, "<=", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new(">=");
	tmp_f->n_arg = 2;
	tmp_f->val = OP_MAGG_U;
	tmp_f->tipo = INFISSO;
	tmp_f->pr = 1;
	g_hash_table_insert(funzioni, ">=", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("!=");
	tmp_f->n_arg = 2;
	tmp_f->val = OP_DIVERSO;
	tmp_f->tipo = INFISSO;
	tmp_f->pr = 1;
	g_hash_table_insert(funzioni, "!=", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("sinh");
	tmp_f->n_arg = 1;
	tmp_f->val = OP_SENH;
	tmp_f->tipo = FUNZIONE;
	tmp_f->pr = 1;
	g_hash_table_insert(funzioni, "sinh", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("cosh");
	tmp_f->n_arg = 1;
	tmp_f->val = OP_COSH;
	tmp_f->tipo = FUNZIONE;
	tmp_f->pr = 1;
	g_hash_table_insert(funzioni, "cosh", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("tanh");
	tmp_f->n_arg = 1;
	tmp_f->val = OP_TANH;
	tmp_f->tipo = FUNZIONE;
	tmp_f->pr = 1;
	g_hash_table_insert(funzioni, "tanh", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("asinh");
	tmp_f->n_arg = 1;
	tmp_f->val = OP_ASENH;
	tmp_f->tipo = FUNZIONE;
	tmp_f->pr = 1;
	g_hash_table_insert(funzioni, "asinh", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("acosh");
	tmp_f->n_arg = 1;
	tmp_f->val = OP_ACOSH;
	tmp_f->tipo = FUNZIONE;
	tmp_f->pr = 1;
	g_hash_table_insert(funzioni, "acosh", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("atanh");
	tmp_f->n_arg = 1;
	tmp_f->val = OP_ATANH;
	tmp_f->tipo = FUNZIONE;
	tmp_f->pr = 1;
	g_hash_table_insert(funzioni, "atanh", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("sqrt");
	tmp_f->n_arg = 1;
	tmp_f->val = OP_RADQ;
	tmp_f->tipo = FUNZIONE;
	tmp_f->pr = 1;
	g_hash_table_insert(funzioni, "sqrt", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("angle");
	tmp_f->n_arg = 2;
	tmp_f->val = OP_ANGOLO;
	tmp_f->tipo = FUNZIONE;
	tmp_f->pr = 1;
	g_hash_table_insert(funzioni, "angle", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("ln");
	tmp_f->n_arg = 1;
	tmp_f->val = OP_LN;
	tmp_f->tipo = FUNZIONE;
	tmp_f->pr = 1;
	g_hash_table_insert(funzioni, "ln", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("sgn");
	tmp_f->n_arg = 1;
	tmp_f->val = OP_SEGNO;
	tmp_f->tipo = FUNZIONE;
	tmp_f->pr = 1;
	g_hash_table_insert(funzioni, "sgn", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("round");
	tmp_f->n_arg = 1;
	tmp_f->val = OP_APPROX;
	tmp_f->tipo = FUNZIONE;
	tmp_f->pr = 1;
	g_hash_table_insert(funzioni, "round", (gpointer) tmp_f);
	tmp_f = alloc(1, struct Elem);
	tmp_f->nome = g_string_new("trunc");
	tmp_f->n_arg = 1;
	tmp_f->val = OP_TRONCA;
	tmp_f->tipo = FUNZIONE;
	tmp_f->pr = 1;
	g_hash_table_insert(funzioni, "trunc", (gpointer) tmp_f);
}

static void CancellaNome(gpointer key, gpointer value, gpointer user_data)
{
	if (((struct Elem *) value)->nome != NULL)
		g_string_free(((struct Elem *) value)->nome, TRUE);
}

void Cancella()
{
	g_hash_table_foreach(costanti, CancellaNome, NULL);
	g_hash_table_foreach(variabili, CancellaNome, NULL);
	g_hash_table_foreach(funzioni, CancellaNome, NULL);
}

struct Elem *FCV(const char *nome, int tipo)
{
	struct Elem *ris = NULL;

	if (tipo | 1)
		ris = (struct Elem *) g_hash_table_lookup(funzioni, nome);
	if (tipo | 2)
		ris = (struct Elem *) g_hash_table_lookup(costanti, nome);
	if (tipo | 4)
		ris = (struct Elem *) g_hash_table_lookup(variabili, nome);
	return ris;
}

double ValV(const char *nome)
{
	struct Elem *ris = (struct Elem *) g_hash_table_lookup(variabili, nome);

	if (ris != NULL)
		return ris->val;
	else
		return 0.0;
}

double ValC(const char *nome)
{
	struct Elem *ris = (struct Elem *) g_hash_table_lookup(costanti, nome);

	if (ris != NULL)
		return ris->val;
	else
		return 0.0;
}

void DefC(const char *nome, double val)
{
	struct Elem *new_c;

	new_c = g_hash_table_lookup(costanti, nome);
	if (new_c == NULL) {
		new_c = alloc(1, struct Elem);
		new_c->nome = g_string_new(nome);
		new_c->n_arg = 0; /* le costanti hanno zero argomenti ed un nome */
		g_hash_table_insert(costanti, nome, (gpointer) new_c);
	}
	new_c->val = val;
	return;
}

void DefV(const char *nome)
{
	struct Elem *new_v;

	new_v = g_hash_table_lookup(variabili, nome);
	if (new_v == NULL) {
		new_v = alloc(1, struct Elem);
		new_v->nome = g_string_new(nome);
		new_v->paren = nv + 1; /* questo indicizza le variabili */
		new_v->n_arg = 0; /* le variabili hanno zero argomenti e saranno senza nome */
		g_hash_table_insert(variabili, nome, (gpointer) new_v);
	}
	return;
}

GString *ModificaPF(const char *nome, int lev)
{
	struct Elem *tmp_f = FCV(nome, 1);
	GString *errore = g_string_new("");

	if (tmp_f == NULL) {
		g_string_printf(errore, "function '%s' does not exists!", nome);
		return errore;
	}
	tmp_f->pr = lev;
	g_string_free(errore, TRUE);
	return NULL;
}

void append_key(void *k, void *v1, GList **keys)
{
	*keys = g_list_append(*keys, v1);
}

GString *Informa()
{
	GString *msg = g_string_new("");
	GList *lista = NULL, *li = NULL;
	struct Elem *el = NULL;

	g_string_printf(msg, "Predefined constants: Pi, e\n");
	g_string_append_printf(msg, "Predefined variables: x, y, z\n");
	g_string_append_printf(msg, "Predefined functions: name/args (priority)\n");
	g_hash_table_foreach(funzioni, append_key, &lista);
	for (li = lista; li != NULL && li->data != NULL; li = li->next) {
		el = (struct Elem *) li->data;
		g_string_append_printf(msg, "\t%s/%d (%d)\n", el->nome->str, el->n_arg, el->pr);
	}
	g_list_free(lista);
	return msg;
}

/* args = 0, lev = 0 */
GString *DefF(const char *nome, const char *expr, int args, int lev)
{
	int t = 0, l = 0, len;
	int e, arg = 0;
	int vir, par = 0;
	GString *tmp = g_string_new("");
	GString *errore = g_string_new("");
	struct Elem *tmpe, *tmpe1, *new_f;
	enum TipoElem vfc = ERRORE;
	enum Stato s = INIZIO;

	new_f = g_hash_table_lookup(funzioni, nome);
	if (new_f == NULL) {
		new_f = alloc(1, struct Elem);
		g_hash_table_insert(funzioni, nome, (gpointer) new_f);
		new_f->def = NULL;
	}
	if (new_f->def != NULL) {
		libera(new_f->def);
		new_f->def = NULL;
	}
	len = strlen(expr);
	new_f->def = alloc(len, struct Elem);
	new_f->len = len;
	new_f->nome = g_string_new(nome);
	new_f->n_arg = args;
	new_f->pr = lev;
	p = g_queue_new();
	/* g_queue_push_tail(p, (gpointer) new_f); */
	while (t < len) {
		while (expr[t] == ' ')
			t++;
		if (expr[t] == '\0')
			break;
		else if (g_ascii_isdigit(expr[t])) {
			if ((stato[s] & OPERANDO1) == 0) {
				g_string_printf(errore, "syntax error in the expression '%s': unexpected number at position %d!", expr, t + 1);
				Cancella();
				return errore;
			}
			s = OPERANDO;
			new_f->def[l].n_arg = 0;
			new_f->def[l].nome = g_string_new("");
			new_f->def[l].val = 0;
			vir = 0;
			e = 10;
			do {
				if (expr[t] == '.') {
					t++;
					if (vir == 1) {
						g_string_printf(errore, "syntax error: double decimal at position %d!", t);
						Cancella();
						return errore;
					}
					vir = 1;
				}
				if (vir == 1) {
					new_f->def[l].val += (double)(expr[t] - 48) / e;
					e *= 10;
				} else {
					new_f->def[l].val = new_f->def[l].val * 10 + (expr[t] - 48);
				}
				if (t > len)
					break;
				t++;
			} while (g_ascii_isdigit(expr[t]) || expr[t] == '.');
			l++;
			arg++;
		}
		else if (expr[t] == '(') {
			if ((stato[s] & PAR_A1) == 0) {
				g_string_printf(errore, "syntax error in the expression '%s': unexpected open parenthesis at position %d!", expr, t + 1);
				Cancella();
				return errore;
			}
			s = PAR_A;
			par++;
			t++;
		}
		else if (expr[t] == ')' || expr[t] == ',') {
			if (expr[t] == ')') {
				if ((stato[s] & PAR_C1) == 0) {
					g_string_printf(errore, "syntax error in the expression '%s': unexpected closed parenthesis at position %d!", expr, t + 1);
					Cancella();
					return errore;
				}
				s = PAR_C;
			}
			else {
				if ((stato[s] & SEP) == 0) {
					g_string_printf(errore, "syntax error in the expression '%s': unexpected separator at position %d!", expr, t + 1);
					Cancella();
					return errore;
				}
				s = SEP;
			}
			tmpe = (struct Elem *) g_queue_peek_tail(p);
			while (tmpe != NULL && tmpe->paren >= par && tmpe->n_arg > 0) {
				arg -= tmpe->n_arg;
				arg++;
				new_f->def[l] = *(struct Elem *) g_queue_pop_tail(p);
				l++;
				tmpe = (struct Elem *) g_queue_peek_tail(p);
			}
			if (expr[t] == ')')
				par--;
			t++;
		}
		else if (!g_ascii_isdigit(expr[t]) && expr[t] != '(' && expr[t] != ')' && expr[t] != ',') {
			g_string_assign(tmp, "");
			/* gli operatori hanno comunque un carattere, quindi lo memorizzo */
			do {
				g_string_append_c(tmp, expr[t]);
				t++;
			} while (t < len && g_ascii_isalpha(expr[t - 1]) && g_ascii_isalpha(expr[t]));
			/* a questo punto se sono uscito ho solo un carattere: ok */
			/* gestisco il "-" unario (lo e` se all'inizio o dopo un operatore prefisso o '(' o un operatore infisso o un argomento) */
			if ((s == INIZIO || s == OP_PREFISSO || s == PAR_A || s == OP_INFISSO || s == SEP) && expr[t - 1] == '-')
				g_string_assign(tmp, "§");
			/* gestisco <= */
			if (expr[t - 1] == '<' && expr[t] == '=') {
				g_string_assign(tmp, "<=");
				t++;
			}
			/* gestisco >= */
			if (expr[t - 1] == '>' && expr[t] == '=') {
				g_string_assign(tmp, ">=");
				t++;
			}
			/* gestisco != */
			if (expr[t - 1] == '!' && expr[t] == '=') {
				g_string_assign(tmp, "!=");
				t++;
			}
			if ((tmpe1 = g_hash_table_lookup(funzioni, tmp->str)) != NULL)
				vfc = FUNZ;
			else if ((tmpe1 = g_hash_table_lookup(costanti, tmp->str)) != NULL)
				vfc = COSTANTE;
			else if ((tmpe1 = g_hash_table_lookup(variabili, tmp->str)) != NULL)
				vfc = VARIABILE;
			if (tmpe1 == NULL) {
				g_string_printf(errore, "syntax error in the expression '%s': symbol '%s' does not exists", expr, tmp->str);
				Cancella();
				return errore;
			}
			if (vfc == COSTANTE) {
				if ((stato[s] & OPERANDO1) == 0) {
					g_string_printf(errore, "syntax error in the expression '%s': unexpected constant at position %d!", expr, t);
					Cancella();
					return errore;
				}
				s = OPERANDO;
				new_f->def[l].nome = g_string_new(tmpe1->nome->str);
				new_f->def[l].n_arg = 0;
				new_f->def[l].val = tmpe1->val;
				l++;
				 arg++;
			}
			else if (vfc == VARIABILE) {
				if ((stato[s] & OPERANDO1) == 0) {
					g_string_printf(errore, "syntax error in the expression '%s': unexpected variable at position %d!", expr, t);
					Cancella();
					return errore;
				}
				s = OPERANDO;
				new_f->def[l].nome = NULL;
				new_f->def[l].n_arg = 0;
				new_f->def[l].paren = tmpe1->paren;
				l++;
				arg++;
			}
			else if (vfc == FUNZ) {
				switch (tmpe1->tipo) {
					case PREFISSO:
						if ((stato[s] & OP_PREFISSO1) == 0) {
							g_string_printf(errore, "syntax error: unexpected prefix operator at position %d!", t);
							Cancella();
							return errore;
						}
						s = OP_PREFISSO;
						break;
					case FUNZIONE:
						if ((stato[s] & OP_FUNZIONE1) == 0) {
							g_string_printf(errore, "syntax error: unexpected function at position %d!", t);
							Cancella();
							return errore;
						}
						s = OP_FUNZIONE;
						break;
					case INFISSO:
						if ((stato[s] & OP_INFISSO1) == 0) {
							g_string_printf(errore, "syntax error: unexpected infix operator at position %d!", t);
							Cancella();
							return errore;
						}
						s = OP_INFISSO;
						break;
					case POSTFISSO:
						if ((stato[s] & OP_POSTFISSO1) == 0) {
							g_string_printf(errore, "syntax error: unexpected postfix operator at position %d!", t);
							Cancella();
							return errore;
						}
						s = OP_POSTFISSO;
						break;
				}
				tmpe = (struct Elem *) g_queue_peek_tail(p);
				tmpe1->paren = par;
				if (tmpe1 != NULL && (tmpe == NULL || (tmpe1->pr > tmpe->pr && par <= tmpe->paren))) {
					while (tmpe != NULL && tmpe->pr <= tmpe1->pr) {
						tmpe = (struct Elem *) g_queue_pop_tail(p);
						arg -= tmpe->n_arg;
						new_f->def[l] = *tmpe;
						l++;
						arg++;
						tmpe = (struct Elem *) g_queue_peek_tail(p);
					}
				}
				g_queue_push_tail(p, tmpe1);
			}
		}
	}
	tmpe = (struct Elem *) g_queue_pop_tail(p);
	while (tmpe != NULL && tmpe->pr > 0) {
		arg -= tmpe->n_arg;
		new_f->def[l] = *tmpe;
		l++;
		arg++;
		tmpe = (struct Elem *) g_queue_pop_tail(p);
	}
	if (arg < 1) {
		g_string_assign(errore, "missing arguments!");
		Cancella();
		return errore;
	}
	if (arg > 1) {
		g_string_assign(errore, "too much arguments!");
		Cancella();
		return errore;
	}
	if (par != 0) {
		g_string_assign(errore, "parentheses do not match!");
		Cancella();
		return errore;
	}
	new_f->len = l;
	g_queue_free(p);
	g_string_free(tmp, TRUE);
	g_string_free(errore, TRUE);
	return NULL;
}

double Fatt(int n)
{
	int k;

	if (n == 0)
		return 1;
	for (k = n - 1; k > 0; k--)
		n *= k;
	return n;
}

/* a = NULL */
GString *Calcola(const char *nome, double *a, double *ris)
{
	double args[10];
	GString *op_nome = g_string_new("");
	GString *errore = g_string_new("");
	int op, pos = 0;
	double r, *p;
	int k;
	unsigned int t = 0;
	struct Elem *f;

	f = (struct Elem *) g_hash_table_lookup(funzioni, nome);
	if (f == NULL) {
		g_string_printf(errore, "semantic error: function '%s' does not exists!", nome);
		return errore;
	}
	p = alloc(f->len, double);
	while (t < f->len) {
		if (f->def[t].n_arg == 0) {
			if (f->def[t].nome != NULL) {
				p[pos++] = f->def[t].val;
			} else {
				p[pos++] = a[f->def[t].paren - 1];
			}
			t++;
		}
		else {
			for (k = f->def[t].n_arg - 1; k > -1; k--)
				args[k] = p[--pos];
			op = (int) f->def[t].val;
			g_string_assign(op_nome, f->def[t].nome->str);
			switch (op) {
				case OP_NEG:
					r = -args[0];
					break;
				case OP_SOMMA:
					r = args[0] + args[1];
					break;
				case OP_SOTTR:
					r = args[0] - args[1];
					break;
				case OP_MOLT:
					r = args[0] * args[1];
					break;
				case OP_DIV:
					r = args[0] / args[1];
					break;
				case OP_POT:
					r = pow(args[0], args[1]);
					break;
				case OP_ABS:
					r = fabs(args[0]);
					break;
				case OP_SEN:
					r = sin(args[0]);
					break;
				case OP_COS:
					r = cos(args[0]);
					break;
				case OP_LOG:
					r = log(args[0]) / log(args[1]);
					break;
				case OP_ARCTAN:
					r = atan(args[0]);
					break;
				case OP_TAN:
					r = tan(args[0]);
					break;
				case OP_MOD:
					r = fmod(args[0], args[1]);
					break;
				case OP_SE:
					r = (args[0] > 0.0) ? args[1] : args[2];
					break;
				case OP_ASEN:
					r = asin(args[0]);
					break;
				case OP_ACOS:
					r = acos(args[0]);
					break;
				case OP_FATT:
					r = Fatt((int) args[0]);
					break;
				case OP_UGUALE:
					r = fabs(args[0] - args[1]) <= DBL_EPSILON * fabs(args[0]);
					break;
				case OP_NOT:
					r = !(args[0] > 0.0);
					break;
				case OP_MAGG:
					r = args[0] > args[1];
					break;
				case OP_MINORE:
					r = args[0] < args[1];
					break;
				case OP_MAGG_U:
					r = args[0] >= args[1];
					break;
				case OP_MINORE_U:
					r = args[0] <= args[1];
					break;
				case OP_DIVERSO:
					r = args[0] <= args[1];
					break;
				case OP_SENH:
					r = sinh(args[0]);
					break;
				case OP_COSH:
					r = cosh(args[0]);
					break;
				case OP_TANH:
					r = tanh(args[0]);
					break;
				case OP_ASENH:
					r = log(args[0] + sqrt(args[0] * args[0] + 1));
					break;
				case OP_ACOSH:
					r = log(args[0] + sqrt(args[0] * args[0] - 1));
					break;
				case OP_ATANH:
					r = log((1 + args[0]) / (1 - args[0])) / 2;
					break;
				case OP_RADQ:
					r = sqrt(args[0]);
					break;
				case OP_ANGOLO:
					r = atan2(args[0], args[1]);
					break;
				case OP_LN:
					r = log(args[0]);
					break;
				case OP_SEGNO:
					r = (args[0] >= 0)? ((args[0] > 0)? 1 : 0) : -1;
					break;
				case OP_APPROX:
					r = floor(args[0] + 0.5);
					break;
				case OP_TRONCA:
					r = (int) args[0];
					break;
				default:
					errore = Calcola(op_nome->str, args, &r);
					if (errore != NULL) {
						g_string_printf(errore, "error: %s!", errore->str);
						return errore;
					}
			}
			p[pos++] = r;
			if (t == f->len)
				break;
			t++;
		}
	}
	*ris = p[--pos];
	libera(p);
	g_string_free(op_nome, TRUE);
	g_string_free(errore, TRUE);
	return NULL;
}
