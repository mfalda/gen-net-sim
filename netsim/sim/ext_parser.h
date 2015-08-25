#ifndef EXT_PARS_H
#define EXT_PARS_H

#include <R.h>
#include <Rdefines.h>

#include "muParserDLL.h"

#define PARSER_CONST_PI  3.141592653589793238462643
#define PARSER_CONST_E   2.718281828459045235360287

void OnError(muParserHandle_t hParser);

muFloat_t Not(muFloat_t x);
muFloat_t Mod(muFloat_t v1, muFloat_t v2);
muFloat_t Angle(muFloat_t x, muFloat_t y);
muFloat_t PRamp(muFloat_t x, muFloat_t n, muFloat_t pos, muFloat_t neg);
muFloat_t Periodic1(muFloat_t x, muFloat_t n, muFloat_t a, muFloat_t b);
muFloat_t Periodic2(muFloat_t x, muFloat_t n, muFloat_t a, muFloat_t b);
muFloat_t Fact(muFloat_t n);
muFloat_t Pulse(muFloat_t x, muFloat_t p, muFloat_t r);

muParserHandle_t InitCalc();
void DeInitCalc(muParserHandle_t hParser);

#endif
