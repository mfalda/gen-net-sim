#include "dinamica.h"

#define g_tipi globali.tipi
#define g_ntipi globali.n_tipi

// arrotonda a d cifre significative
double signif(const double n, const unsigned int d)
{
	double scala, ris;

	if (Uguale(n, 0.0))
		ris = n;
	else {
		// trovo un fattore di dieci che sposti tutte le cifre significative dopo la virgola
		scala = pow(10.0, (double)((int) d - 1 - (int) floor(log10(fabs(n)))));
		// arrotonda
      ris = floor(n * scala + 0.5) / scala;
   }
	//~ Rprintf("signif (%d): %f -> %f\g_n", d, g_n, ris);
	return ris;
}

#define g_D globali.dinamica.D
#define g_ind globali.dinamica.ind
#define g_y_prec globali.dinamica.y_prec
#define g_targetT globali.dinamica.targetT
#define g_targ globali.dinamica.targ
#define g_incr globali.dinamica.incr
#define g_n globali.dinamica.n
#define g_aus globali.dinamica.aus
#define g_tmp_ris globali.dinamica.tmp_ris

// dinamica<-function(parms,n0,times){
MATRICEd *dinamica1(MATRICEd *ris, LISTA *parms, LISTA *regole, LISTA *ext_fun, VETTOREd *n0, VETTOREd *times, muParserHandle_t *hParsers, double stat_thr, double stat_width, double sd_noise)
{
	double res, s, media;
	int i, j, passi, l, indx;
	int st = 0, st_all;
	double y1 = 0.0;
	double max_stat;
	GString *tmp = NULL;
	GString* act_fun = NULL;
	// le prossime sono solo puntatori
	MATRICEd *M = NULL, *ext_in = NULL;
	MATRICEi *N = NULL;
	VETTOREd *k = NULL, *alpha = NULL, *theta = NULL, *xmin = NULL;

	_Intestazione("\n*** dinamica1 ***\n");

	// M<-parms[[1]]
	CtrlLlst(parms, 1);
	M = ACCEDIlst(parms, 1, md);
	// R<-parms[[2]]
	// da parametro
	// N<-parms[[3]]
	CtrlLlst(parms, 2);
	N = ACCEDIlst(parms, 2, mi);
	// k<-parms[[4]]
	CtrlLlst(parms, 3);
	k = ACCEDIlst(parms, 3, vd);
	// alpha<-parms[[5]]
	CtrlLlst(parms, 4);
	alpha = ACCEDIlst(parms, 4, vd);
	// theta<-parms[[6]]
	CtrlLlst(parms, 5);
	theta = ACCEDIlst(parms, 5, vd);
	// act.fun<-parms[[7]]
	CtrlLlst(parms, 6);
	act_fun = ACCEDIlst(parms, 6, str);
	// EXT.IN<-parms[[8]]
	CtrlLlst(parms, 7);
	ext_in	 = ACCEDIlst(parms, 7, md);
	// EXT.FUN<-parms[[9]]
	// da parametro
	// l'ottavo non mi serve
	CtrlLlst(parms, 8);
	xmin = ACCEDIlst(parms, 8, vd);

	// res<-signif((log(2)/mean(k))/length(times),1)
	media = media_v_d(k);
#ifdef MDEBUG
	if (Uguale(media, 0.0)) {
		CREAstr(tmp, "");
		g_string_printf(tmp, "ATTENZIONE (dinamica.c, linea %d): divisione per zero!\n", __LINE__);
		warning(tmp->str);
		fprintf(fp_fdbg, tmp->str);
		CANCELLAstr(tmp);
	}
#endif
	res = signif(0.693147180559945309417 / media / LENGTHv_d(times), 1);
	// S<-max(abs(times/res-round(times/res)))
	s = f_aux4_d(times, res);
	// while (S>1e-8)
	while (s > 1e-8 && s != R_PosInf) { // aggiunto un ulteriore controllo
	//    {res<-res/10; S<-max(abs(times/res-round(times/res)))}
		res /= 10.0;
		s = f_aux4_d(times, res);
	}

	// passi<-max(times)/res
	passi = (int) (max_v_d(times) / res);
	if (res < 1e-8) {
		passi = 1000;
		res = max_v_d(times) / passi;
		CREAstr(tmp, "");
		g_string_printf(tmp, "resolution for Euler algorithm has been automatically set to %.16g!\n", res);
		warning(tmp->str);
#ifdef FDEBUG
		fprintf(fp_fdbg, tmp->str);
#endif
		CANCELLAstr(tmp);
	}

	// g_D<-matrix(n0,length(n0))
	CREAm_d(g_tmp_ris, LENGTHv_d(n0), passi);
	CREAv_d(g_y_prec, LENGTHv_d(n0));

	if (stat_width == 0.0)
		max_stat = LENGTHv_d(times);
	else
		max_stat = (double) stat_width * LENGTHv_d(times);
	g_n = copia_v_d(g_n, n0, 1, LENGTHv_d(n0));
	// for(i in 1:passi){
	for (i = 1; (!stat_width || st < max_stat) && i <= passi; i++) {
	// 	g_n<-g_D[,i]
	// non lo rileggo da D
	// 	g_targ<-target(g_n,R,M,N,EXT.IN,EXT.FUN,res*(i-1))
		g_targ = target1(g_targ, g_n, regole, M, N, ext_in, ext_fun, res * (i - 1), hParsers, sd_noise);
	//         if (act.fun=="linear")   g_targetT<-targ
		if (!strncmp(act_fun->str, "linear", 6))
			g_targetT = f_aux8_d(g_targetT, g_targ, xmin);
	//          else  {#norm.f<-(1-theta)^alpha/((1-theta)^alpha+0.5^alpha)
		else  {//norm_f=(1-theta)^alpha/((1-theta)^alpha+0_5^alpha);
	// 		#g_aus<-(g_targ-theta)
			//g_aus=(g_targ-theta);
	// 		#g_aus[which(g_aus<0)]<-0
			//g_aus[which(g_aus<0)]=0;
	// 		#g_targetT<-(g_aus^alpha/(g_aus^alpha+0.5^alpha))/norm.f
			//g_targetT=(g_aus^alpha/(g_aus^alpha+0_5^alpha))/norm_f;
	// 		g_targetT<-1/(1+exp(-alpha*(g_targ-theta))) * (1 - xmin) + xmin
			g_targetT = f_aux5_d(g_targetT, alpha, g_targ, theta, xmin);
	// 		#g_aus<-apply(M,1,sum)
			//g_aus=apply(M,1,sum);
	// 		#g_targetT[which(g_aus==0)]<-0
			//g_targetT[which(g_aus==0)]=0;
// 	       }
		}
	// 	# Euler Algorithm
	// 	g_incr<-matrix(res*k*(g_targetT-g_n),ncol=1)
		g_incr = f_aux6_d(g_incr, res, k ,g_targetT, g_n);
	// 	g_n<-g_n+incr
		somma1_vv_d(g_n, g_incr);
	// 	g_D<-cbind(g_D,g_n)
		st_all = 1;
		for (j = 1; j <= LENGTHv_d(n0); j++) {
			ASSEGNAm_d(g_tmp_ris, j, i, ACCEDIv_d(g_n, j));
			st_all &= (fabs(ACCEDIv_d(g_n, j) - ACCEDIv_d(g_y_prec, j)) < stat_thr);
		}
		if (st_all)
			st++;
		g_y_prec = copia_v_d(g_y_prec, g_n, 1, LENGTHv_d(g_n));
	// }
	}

	CREAm_d(g_D, LENGTHv_d(n0), i - 1);
	assegna1_mv_colonna_d(g_D, 1, n0);
	for (j = 1; j < i; j++)
		copia1_m_colonna_d(g_D, j, g_tmp_ris, j);

	// L<-length(times)
	l = LENGTHv_d(times);
	// g_ind<-rep(0,L)
	CREAv_i(g_ind, l);
	InitVett_i(g_ind, 0);
	// g_aus<-seq(0,max(times),res)
	g_aus = seq_d(g_aus, 0.0, max_v_d(times), res);
	// for (i in (1:L))
	for (i = 1; i <= l; i++) {
	//   {g_ind[i]<-which.min(abs(g_aus-times[i]))
		indx = f_aux7_d(g_aus, ACCEDIv_d(times, i));
		ASSEGNAv_i(g_ind, i, indx);
	//    if ( (g_aus[g_ind[i]]-times[i])>1e-8) stop("error")
		if (ACCEDIv_d(g_aus, ACCEDIv_i(g_ind, i)) - ACCEDIv_d(times, i) > 1e-8)
			error("error");
	//   }
	}
	// g_D<-g_D[,g_ind]
	ris = seleziona_colonne_d(ris, g_D, g_ind);
	// if (length(n0)==1) g_D<-matrix(g_D,nrow=1)
	if (LENGTHv_d(n0) == 1) {
		cambiadim1_d(ris, 1, -1);
	}
	//~ CANCELLAm_d(g_D);
	//~ CANCELLAv_i(g_ind);
	//~ CANCELLAv_d(g_targetT);
	//~ CANCELLAv_d(g_targ);
	//~ CANCELLAv_d(g_incr);
	//~ CANCELLAv_d(g_n);
	//~ CANCELLAv_d(g_aus);

	StrBilanciam();

	// return(g_D)
	return ris;
}

SEXP dinamica(SEXP parms, SEXP regole, SEXP ext_fun, SEXP n0, SEXP times)
{
	int i, ll, ll1, ll2, ll3, nProtected = 0;
	VETTOREd *n01 = NULL, *times1 = NULL;
	MATRICEd *ris1 = NULL;
	LISTA *parms1 = NULL, *regole1 = NULL, *ext_fun1 = NULL;
#ifdef MDEBUG
	GString **nomi;
	char tmp[10];
#endif
	enum TIPO tipi[7];
	SEXP ris;

	_InitDbg(false, false, false);

	_Intestazione("\n*** dinamica ***\n");

	ll1 = length(parms);
	ll2 = length(regole);
	ll3 = length(ext_fun);
	ll = max_s_i(ll1, max_s_i(ll2, ll3));
#ifdef MDEBUG
	nomi = mia_alloc(ll, GString *);
	if (nomi == NULL) {
		Rprintf("Not enough memory (dinamica1 # %d, nomi)", __LINE__ - 2);
		error("");
	}
#endif
	n01 = inVETTORE_d(n0, &nProtected);
	times1 = inVETTORE_d(times, &nProtected);
	tipi[0] = MATRi;
	tipi[1] = MATRi;
	tipi[2] = VETTd;
	tipi[3] = VETTd;
	tipi[4] = VETTd;
	tipi[5] = STRINGA;
	tipi[6] = MATRd;
	for (i = 0; i < ll; i++) {
#ifdef MDEBUG
		if (i < ll1)
			snprintf(tmp, 14, "Lista parms %d", i + 1);
		else
			tmp[0] = '\0';
		CREAstr(nomi[i], tmp);
#endif
	}
	parms1 = inLISTA(parms, &nProtected, ll1, tipi, nomi);
	for (i = 0; i < ll2; i++) {
#ifdef MDEBUG
		g_string_printf(nomi[i], "Lista regole %d", i + 1);
#endif
	}
	regole1 = inLISTA(regole, &nProtected, ll2, NULL, nomi);
	if (g_tipi == NULL) {
		g_ntipi = ll3;
		g_tipi = mia_alloc(ll3, enum TIPO);
		if (ll3 > 0 && g_tipi == NULL) {
			Rprintf("Not enough memory (dinamica # %d, tipi)", __LINE__ - 2);
		error("");
		}
	}
	assert(g_ntipi == ll3);
	for (i = 0; i < ll3; i++) {
		g_tipi[i] = STRINGA;
#ifdef MDEBUG
		g_string_printf(nomi[i], "Lista ext_fun %d", i + 1);
#endif
	}
	ext_fun1 = inLISTA(ext_fun, &nProtected, ll3, g_tipi, nomi);
#ifdef MDEBUG
	for (i = 0; i < ll; i++)
		CANCELLAstr(nomi[i]);
	libera(nomi);
#endif
	ris1 = dinamica1(ris1, parms1, regole1, ext_fun1, n01, times1, NULL, 0, 0, 0);
	ris = daMATRICE_d(ris1, &nProtected);

	// qui va bene perche´ non la ripasso a daLISTA
	CancellaLISTA(parms1, true);
	CancellaLISTA(regole1, true);
	CancellaLISTA(ext_fun1, true);

	CANCELLAv_d(n01);
	CANCELLAv_d(times1);

	StrBilanciam();
	ControllaCanc();

	UNPROTECT(nProtected);

	return ris;
}
