#include "pars.h"


int main(int argc, char *argv[])
{
	char *formula;
	GString *errore, *inf;
	double args[4], ris;

	Ripristina();
	if (argc < 2) {
		printf("Sintassi: prova \"funzione in x\"\n");
		inf = Informa();
		printf("%s\n", inf->str);
		g_string_free(inf, TRUE);
		Cancella();
		return 0;
	}
	DefC("T", 1);
	DefV("t");
	formula = argv[1];
	printf("Per x = 4.5 e t = -4.5 => %s = ", formula);
	errore = DefF("main", formula, 2, 0);
	if (errore == NULL) {
		args[0] = 4.5;
		args[3] = -4.5; /* le prime tre sono x, y e z */
		errore = Calcola("main", args, &ris);
		if (errore == NULL) {
			printf("%f\n", ris);
			Cancella(); /* cancella solo se OK, altrimenti lo fa automaticamente! */
		}
		else
			printf("%s\n", errore->str);
	}
	else
		printf("%s\n", errore->str);
	return 0;
}

