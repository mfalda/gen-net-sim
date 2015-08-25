#include "ext_parser.h"

// Callback function for parser errors
void OnError(muParserHandle_t hParser)
{
	char err[256];

	snprintf(err, 256, "Error:\n");
	snprintf(err, 256, "%s------\n", err);
	snprintf(err, 256, "%sMessage:  \"%s\"\n", err, mupGetErrorMsg(hParser));
	snprintf(err, 256, "%sToken:    \"%s\"\n", err, mupGetErrorToken(hParser));
	snprintf(err, 256, "%sPosition: %d\n", err, mupGetErrorPos(hParser));
	snprintf(err, 256, "%sErrc:     %d\n", err, mupGetErrorCode(hParser));

	error(err);
}

// Function callbacks
muFloat_t Not(muFloat_t x)
{
	return (ceil(x) == 0.0);
}

muFloat_t Mod(muFloat_t x, muFloat_t y)
{
	return fmod(x, y);
}

muFloat_t Angle(muFloat_t x, muFloat_t y)
{
	return atan2(x, y);
}

muFloat_t f_aux(muFloat_t x, muFloat_t n)
{
	if (x < 0.5)
		return pow(x * 2.0, n) / 2.0;
	else
		return 1.0 - pow((1.0 - x) * 2.0, n) / 2.0;
}

muFloat_t PRamp(muFloat_t x, muFloat_t n, muFloat_t pos, muFloat_t neg)
{
	if (x - floor(x) < 0.5)
		return pos * f_aux((x - floor(x)) * 2.0, n);
	else
		return neg * f_aux((x - 0.5 - floor(x - 0.5)) * 2.0, n);
}

muFloat_t Periodic1(muFloat_t x, muFloat_t n, muFloat_t a, muFloat_t b)
{
	if (x - floor(x) < a)
		return f_aux((x - floor(x)) / a, n);
	else if (x - floor(x) > b)
		return -f_aux((1 - x - floor(1 - x)) / (1 - b), n);
	else
		return 1;
}

muFloat_t Periodic2(muFloat_t x, muFloat_t n, muFloat_t a, muFloat_t b)
{
	if (x - floor(x) < a)
		return f_aux((x - floor(x)) / a, n);
	else if (x - floor(x) > b)
		return f_aux((1 - x - floor(1 - x)) / (1 - b), n);
	else
		return 1;
}

muFloat_t Pulse(muFloat_t x, muFloat_t p, muFloat_t r)
{
	return (x - floor(x) >= p - p / r && x - floor(x) <= p + p / r)? 1 : 0;
}

// Callback for the postfix "!" operator
muFloat_t Fact(muFloat_t n)
{
	int k;

	if (n == 0)
		return 1;
	for (k = n - 1; k > 0; k--)
		n *= k;
	return n;
}

muParserHandle_t InitCalc()
{
	muParserHandle_t hParser = mupCreate(muBASETYPE_FLOAT);              // initialize the parser

	// Set an error handler [optional]
	mupSetErrorHandler(hParser, OnError);

	// Define postfix operators [optional]
	mupDefinePostfixOprt(hParser, "!", Fact, 0);

	// Define infix operator [optional]
	mupDefineInfixOprt(hParser, "~", Not, 0);

	// Define binary operators [optional]
	mupDefineOprt(hParser, "%", Mod, 0, muOPRT_ASCT_LEFT, 0);

	// Define functions [optional]
	mupDefineFun2(hParser, "angle", Angle, 2);
	mupDefineFun4(hParser, "pramp", PRamp, 4);
	mupDefineFun4(hParser, "periodic1", Periodic1, 4);
	mupDefineFun4(hParser, "periodic2", Periodic2, 4);
	mupDefineFun3(hParser, "pulse", Pulse, 3);

	return hParser;
}

void DeInitCalc(muParserHandle_t hParser)
{
	// finally free the parser resources
	mupRelease(hParser);
}
