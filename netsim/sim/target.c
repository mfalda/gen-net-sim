#include "target.h"

#define g_tipi globali.tipi
#define g_ntipi globali.n_tipi

#define g_aus globali.target.aus
#define g_tmp_d globali.target.tmp_d
#define g_tmp1_d globali.target.tmp1_d
#define g_ind globali.target.ind
#define g_aus globali.target.aus
#define g_ind_reg globali.target.ind_reg

/*double calcola_f(int n, double t)
{
	switch (n) {
		case 1:
			return 1;
		case 2:
			return t;
		case 3:
			return sin(t);
		case 4:
			return cos(t);
		default:
			error("funzione per target non definita!\n");
			return -1;
	}
}*/

// target<-function(n,R,M,N,EXT.IN,EXT.FUN,t)
VETTOREd *target1(VETTOREd *ris, const VETTOREd *n, const LISTA *R, const MATRICEd *m, const MATRICEi *N, const MATRICEd *ext_in, const LISTA *ext_fun, double t, muParserHandle_t *hParsers, double sd_noise)
{
	int Ln, i, j, L;
	double expon, val, valr, signr, args[1], ris_calc, fVal;
	char nome[9] = "ext_fun0";
	GString *errore = NULL;
#ifdef MDEBUG
	GString *tmp = NULL;
#endif

	_Intestazione("\n***target***\n");
#ifdef DET
	fprintf(fp_det, "input:\n");
	fprintf(fp_det, "\tn = ");
	_StampaRawVett_d(n);
	fprintf(fp_det, "\tR = ");
	_StampaRawLista(R);
	fprintf(fp_det, "\tM = ");
	_StampaRawMatr_d(m);
	fprintf(fp_det, "\tN = ");
	_StampaRawMatr_i(N);
	fprintf(fp_det, "\tt =  %.16g\n", t);
#endif

	// Ln <- length(n)
	Ln = LENGTHv_d(n);
	//  T<-numeric(Ln)  #== rep(0,Ln)
	CREAv_d(ris, Ln);
	//  if (length(EXT.FUN)==0)
	if (ext_fun->dim == 0) {
		//   {for(i in 1:Ln)
		for (i = 1; i <= Ln; i++) {
			//     {rule<-R[[i]]
			CtrlLlst(R, i);
			// 	   g_aus<-boole.result(rule,n,M,N,i)
			g_aus = boole_result(g_aus, ACCEDIlst(R, i, vi), n, m, N, i);
			//     if (g_aus[1]==(-1))
			if (Uguale(ACCEDIv_d(g_aus, 1), -1.0))
				// T[i]<-1-g_aus[2]
				ASSEGNAv_d(ris, i, 1 - ACCEDIv_d(g_aus, 2));
			//      else    T[i]<-g_aus[2]
			else
				ASSEGNAv_d(ris, i, ACCEDIv_d(g_aus, 2));
			//     }
		}
		//    }
	}
	//  else
	else {
		//   {g_ind.reg<-which(apply(abs(EXT.IN),1,sum)!=0)
		g_tmp_d = abs_m_d(g_tmp_d, ext_in);
		g_tmp1_d = somma_righe_d(g_tmp1_d, g_tmp_d);
		g_ind_reg = which_v_indxne_d(g_ind_reg, g_tmp1_d, 0.0);

		for (i = 0; i < LENGTHm2_d(ext_in); i++) {
			nome[7] = (char) (49 + i);
			CtrlLlst(ext_fun, i + 1);
#ifdef FDEBUG
			fprintf(fp_fdbg, "Compilo la funzione %s definita come '%s' ... ",
				nome, (ACCEDIlst(ext_fun, i + 1, str))->str);
#endif
			mupSetExpr(hParsers[i], (ACCEDIlst(ext_fun, i + 1, str))->str);
			mupDefineVar(hParsers[i], "t", &args[0]);
#ifdef FDEBUG
			fprintf(fp_fdbg, "ok\n");
#endif
		}
		//    for(i in 1:Ln)
		for (i = 1; i <= Ln; i++) {
			//      {rule<-R[[i]]
			CtrlLlst(R, i);
			// 	   g_aus<-boole.result(rule,n,M,N,i)
			g_aus = boole_result(g_aus, ACCEDIlst(R, i, vi), n, m, N, i);
			//       g_aus<-boole.result(rule,n,M,N,i)
			//       if (i %in% g_ind.reg)
			if (esiste_v_i(i, g_ind_reg) > 0) {
				//         {g_ind<-which(EXT.IN[i,]!=0)
				g_ind = which_m_rowindxne_d(g_ind, ext_in, i, 0.0);
				//          L<-length(g_ind)
				L = LENGTHv_i(g_ind);
				//          val<-0
				val = 0;
				//          for (j in (1:L))
				for (j = 1; j <= L; j++) {
					//           {expon<-1/abs(EXT.IN[i,g_ind[j]])
#ifdef MDEBUG
					if (Uguale(ACCEDIm_d(ext_in, i, ACCEDIv_i(g_ind, j)), 0.0)) {
						CREAstr(tmp, "");
						g_string_printf(tmp, "ATTENZIONE (target.c, linea 110): divisione per zero!\n");
						warning(tmp->str);
						fprintf(fp_fdbg, tmp->str);
						CANCELLAstr(tmp);
					}
#endif
					expon = (double) 1.0 / fabs(ACCEDIm_d(ext_in, i, ACCEDIv_i(g_ind, j)));
					//            valr<-EXT.FUN[[g_ind[j]]](t)^expon
					args[0] = t;
#ifdef FDEBUG
					nome[7] = (char) (48 + j);
	fprintf(fp_fdbg, "Calcolo la funzione %s in %.5e: ", nome, args[0]);
#endif
					fVal = mupEval(hParsers[j - 1]);
					if (!mupError(hParsers[j - 1]))
						ris_calc = fVal;
#ifdef FDEBUG
	fprintf(fp_fdbg, "%.5e\n", ris_calc);
#endif
					valr = pow(ris_calc, expon);
					//            signr<-sign(EXT.IN[i,g_ind[j]])
					signr = sign(ACCEDIm_d(ext_in, i, ACCEDIv_i(g_ind, j)));
					//            val<-val+valr*signr
					val += valr * signr;
					//           }
				}
				//          val<-val+g_aus[1]*g_aus[2]
				val += ACCEDIv_d(g_aus, 1) * ACCEDIv_d(g_aus, 2);
				//          if  (val<0) T[i]<-max(0,(1+val))
				if (val < 0.0)
					ASSEGNAv_d(ris, i, max_s_d(0.0, (1.0 + val)));
				//           else T[i]<-min(1,val)
				else
					ASSEGNAv_d(ris, i, min_s_d(1.0, val));
				//         }
			}
			//        else {if (g_aus[1]==(-1)) T[i]<-1-g_aus[2]
			else {
				if (Uguale(ACCEDIv_d(g_aus, 1), -1.0))
					ASSEGNAv_d(ris, i, 1 - ACCEDIv_d(g_aus, 2));
				//               else    T[i]<-g_aus[2]
				else
					ASSEGNAv_d(ris, i, ACCEDIv_d(g_aus, 2));
				//             }
			}
			// if (!is.null(sd_noise)) T[i]<-min(1,max((T[i]+rnorm(1,0,sd_noise)),0))
			if (sd_noise > 0.0)
				ASSEGNAv_d(ris, i, min_s_d(1.0, max_s_d(ACCEDIv_d(ris, i) + rnorm(1.0, sd_noise), 0.0)));
			//       }
		}
		//~ CANCELLAv_i(g_tmp1_d);
		//~ CANCELLAv_i(g_tmp_d);
		//~ CANCELLAv_i(g_ind);
		//~ CANCELLAv_d(g_aus);
		//~ CANCELLAv_i(g_ind_reg);
		//   }
	}

	StrBilanciam();

#ifdef DET
	fprintf(fp_det, "target output:\n");
	fprintf(fp_det, "\tT = ");
	_StampaRawVett_d(ris);
#endif

	// return(T)
	return ris;
	// }
}

SEXP target(SEXP n, SEXP r, SEXP m, SEXP N, SEXP ext_in, SEXP ext_fun, SEXP t)
{
	int i, ll, llr, llf, nProtected = 0;
	MATRICEi *n2 = NULL;
	MATRICEd *m1 = NULL;
	VETTOREd *n1 = NULL, *ris2 = NULL;
	MATRICEd *ext_in1 = NULL;
	LISTA *r1 = NULL, *ext_fun1 = NULL;
	double t1;
#ifdef MDEBUG
	GString **nomi;
	char tmp[10];
#endif
	SEXP ris;

	_InitDbg(false, false, false);

	_Intestazione("\n*** target ***\n");

	llr = length(r);
	llf = length(ext_fun);
	ll = max_s_i(llr, llf);
#ifdef MDEBUG
	nomi = mia_alloc(ll, GString *);
	if (ll > 0 && nomi == NULL) {
		Rprintf("Not enough memory (target # %d, nomi)", __LINE__ - 2);
		error("");
	}
#endif
	n1 = inVETTORE_d(n, &nProtected);
	m1 = inMATRICE_d(m, &nProtected);
	n2 = inMATRICE_i(N, &nProtected);
	ext_in1 = inMATRICE_d(ext_in, &nProtected);
	t1 = NUMERIC_VALUE(t);
	for (i = 0; i < ll; i++) {
#ifdef MDEBUG
		if (i < llr)
			snprintf(tmp, 10, "Lista r %d", i + 1);
		else
			tmp[0] = '\0';
		CREAstr(nomi[i], tmp);
#endif
	}
	r1 = inLISTA(r, &nProtected, llr, NULL, nomi);
	if (g_tipi == NULL) {
		g_ntipi = llf;
		g_tipi = mia_alloc(llf, enum TIPO);
		if (llf > 0 && g_tipi == NULL) {
			Rprintf("Not enough memory (target # %d, tipi)", __LINE__ - 2);
		error("");
		}
	}
	assert(g_ntipi == llf);
	for (i = 0; i < llf; i++) {
		g_tipi[i] = STRINGA;
#ifdef MDEBUG
		g_string_printf(nomi[i], "Lista ext_fun %d", i + 1);
#endif
	}
	ext_fun1 = inLISTA(ext_fun, &nProtected, llf, g_tipi, nomi);
#ifdef MDEBUG
	for (i = 0; i < ll; i++)
		CANCELLAstr(nomi[i]);
	libera(nomi);
#endif

	ris2 = target1(ris2, n1, r1, m1, n2, ext_in1, ext_fun1, t1, NULL, 0);
	ris = daVETTORE_d(ris2, &nProtected);

	// qui va bene perché non la ripasso a daLISTA
	CancellaLISTA(r1, true);
	CancellaLISTA(ext_fun1, true);

	CANCELLAm_d(m1);
	CANCELLAm_i(n2);
	CANCELLAv_d(n1);
	CANCELLAm_d(ext_in1);

	StrBilanciam();
	ControllaCanc();

	UNPROTECT(nProtected);

	return ris;
}
