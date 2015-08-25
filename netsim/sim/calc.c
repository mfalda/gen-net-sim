#include "calc.h"


SEXP calc(SEXP str, SEXP args)
{
	int nProtected = 0;
	double *args1 = NULL, *ris1 = NULL;
	SEXP ris = NULL;

	muChar_t szLine[100];
	muFloat_t fVal = 0, afVarVal[1]; // values of the parser variables

	muParserHandle_t hParser = InitCalc();

	PROTECT(str = AS_CHARACTER(str));
	nProtected++;
	strncpy(szLine, CHAR(STRING_ELT(str, 0)), 100);
	PROTECT(args = AS_NUMERIC(args));
	nProtected++;
	args1 = NUMERIC_POINTER(args);
	afVarVal[0] = args1[0];
	// Define parser variables and bind them to C++ variables [optional]
	mupDefineVar(hParser, "x", &afVarVal[0]);
	PROTECT(ris = NEW_NUMERIC(1));
	nProtected++;
	ris1 = NUMERIC_POINTER(ris);

	mupSetExpr(hParser, szLine);

   fVal = mupEval(hParser);
	if (!mupError(hParser))
		ris1[0] = fVal;

	DeInitCalc(hParser);

	UNPROTECT(nProtected);

	return ris;
}

SEXP calcv(SEXP str, SEXP args)
{
	int nProtected = 0, i;
	R_len_t n_args;
	double *args1 = NULL, *ris1 = NULL;
	SEXP ris = NULL;

	muChar_t szLine[100];
	muFloat_t fVal = 0, afVarVal[1]; // values of the parser variables

	muParserHandle_t hParser = InitCalc();              // initialize the parser

	PROTECT(str = AS_CHARACTER(str));
	nProtected++;
	strncpy(szLine, CHAR(STRING_ELT(str, 0)), 100);
	PROTECT(args = AS_NUMERIC(args));
	nProtected++;
	args1 = NUMERIC_POINTER(args);
	n_args = LENGTH(args);
	PROTECT(ris = NEW_NUMERIC(n_args));
	nProtected++;
	ris1 = NUMERIC_POINTER(ris);
	// Define parser variables and bind them to C++ variables [optional]
	mupDefineVar(hParser, "x", &afVarVal[0]);

	mupSetExpr(hParser, szLine);

	for (i = 0; i < n_args; i++) {
		afVarVal[0] = args1[i];
		fVal = mupEval(hParser);
		if (!mupError(hParser))
			ris1[i] = fVal;
	}

	// finally free the parser resources
	DeInitCalc(hParser);

	UNPROTECT(nProtected);

	return ris;
}
