#include "connectivity_random.h"

int Fattoriale(int n)
{
	int ris = 1;

	while (n > 0) {
		ris *= n;
		n--;
	}
	return ris;
}

#define g_M globali.connectivity_random.M
#define g_Mdiscr globali.connectivity_random.Mdiscr
#define g_aus globali.connectivity_random.aus
#define g_tmp_d globali.connectivity_random.tmp_d
#define g_scalare_i globali.connectivity_random.scalare_i
#define g_Pnum globali.connectivity_random.Pnum
#define g_ind globali.connectivity_random.ind
#define g_num globali.connectivity_random.num
#define g_tmp_i globali.connectivity_random.tmp_i

// connectivityrandom<-function(N=50,max_con=12,k=3,weight_mean=1, weight_sd=0.1)
LISTA *connectivity_random1(LISTA *ris, int N, int max_con, double k, double weight_mean, double weight_sd)
{
	double p;
	int i, j, a, L, F;
	enum TIPO tipi[2];

	_Intestazione("\n***connectivity_random***\n");
#ifdef DET
	fprintf(fp_det, "input:\n");
	fprintf(fp_det, "\tN =  %d\n", N);
	fprintf(fp_det, "\tmax_con =  %d\n", max_con);
	fprintf(fp_det, "\tk =  %.16g\n", k);
	fprintf(fp_det, "\tweight_mean =  %.16g\n", weight_mean);
	fprintf(fp_det, "\tweight_sd =  %.16g\n", weight_sd);
#endif

	GetRNGstate();
	//  if (max_con>N) max_con<-N
	if (max_con > N)
		max_con = N;
	//  g_Mdiscr<-g_M<-matrix(0,ncol=N,nrow=N)
	CREAm_i(g_Mdiscr, N, N);
	InitMatr_i(g_Mdiscr, 0);
	CREAm_d(g_M, N, N);
	InitMatr_d(g_M, 0.0);
	CREAv_i(g_scalare_i, 1);
	//  p<-k/N
	p = (double) k / N;
	//  g_Pnum<-rep(0,max_con)
	CREAv_d(g_Pnum, max_con);
	InitVett_d(g_Pnum, 0.0);
	//  for (j in (1:max_con))
	for (j = 1; j <= max_con; j++) {
		//    {a<-N-j+1
		a = N - j + 1;
		//     L<-N-a
		L = N - a;
		//     F<-1
		F = 1;
		//     if (L>0) {for(i in (0:L)) F<-F*(a+i)}
		if (L > 0) {
			for (i = 0; i <= L; i++)
				F *= a + i;
		}
		//     g_Pnum[j]<-(F/factorial(j))*(p^j)*((1-p)^(N-j))
		ASSEGNAv_d(g_Pnum, j, ((double) F / Fattoriale(j)) * (pow(p, j)) * (pow((1 - p), (N - j))));
		//    }
	}
	//  for (i in (1:N))
	for (i = 1; i <= N; i++) {
		//   {g_num<-sampleB(seq(1,max_con,1),1,prob=g_Pnum)
		ASSEGNAv_i(g_scalare_i, 1, 1);
		g_tmp_i = seq_i(g_tmp_i, 1, max_con, 1);
		// senza ripetizioni
		g_num = sampleB_p(g_num, g_tmp_i, 1, 0, g_Pnum);
		//    g_ind<-sampleB(seq(1,N,1),g_num)
		g_tmp_i = seq_i(g_tmp_i, 1, N, 1);
		// senza ripetizioni
		g_ind = sampleB(g_ind, g_tmp_i, ACCEDIv_i(g_num, 1), 0);
		//    g_Mdiscr[i,g_ind]<-1
		assegna1_ms_rigaindx_i(g_Mdiscr, i, g_ind, 1);
		//   }
	}
	// #######assign weights to the connectivity matrix##########
	//######assign weights to the connectivity matrix##########;
	// g_ind<-which(g_Mdiscr==1,arr.g_ind=TRUE)
	g_ind = which_m_indxeq_i(g_ind, g_Mdiscr, 1);
	// L<-dim(g_ind)[1] // ???
	L = LENGTHv_i(g_ind);
	// g_aus<-abs(rnorm_s(L,weight_mean, weight_sd))
	CREAv_d(g_tmp_d, L);
	g_tmp_d = rnorm_s(g_tmp_d, L, weight_mean, weight_sd, "cr");
	g_aus = abs_v_d(g_aus, g_tmp_d);
	// g_M[g_ind]<-aus
	assegna1_mv_indx_d(g_M, g_ind, g_aus);
	// return(list(g_M,g_Mdiscr))
	PutRNGstate();
	// return(list(g_M,g_Mdiscr))
	tipi[0] = MATRi;
	tipi[1] = MATRi;
	CreaLISTA(ris, tipi, 2);
	ris->dati[0].md = g_M;
	ris->dati[1].mi = g_Mdiscr;
	//~ CANCELLAv_d(g_aus);
	//~ CANCELLAv_d(g_tmp_d);
	//~ CANCELLAv_i(g_scalare_i);
	//~ CANCELLAv_d(g_Pnum);
	//~ CANCELLAv_i(g_ind);
	//~ CANCELLAv_i(g_num);
	//~ CANCELLAv_i(g_tmp_i);

	StrBilanciam();

#ifdef DET
	fprintf(fp_det, "cr output:\n");
	fprintf(fp_det, "\tlist(M,Mdiscr) = ");
	_StampaRawMatr_d(g_M);
	_StampaRawMatr_i(g_Mdiscr);
#endif

	return ris;
}

SEXP connectivity_random(SEXP N, SEXP max_con, SEXP k, SEXP weight_mean, SEXP weight_sd)
{
	int nProtected = 0;
	int N1, max_con1;
	double k1, weight_mean1, weight_sd1;
	LISTA *l = NULL;
	SEXP ris;

	_InitDbg(false, false, false);

	_Intestazione("\n*** connectivity_random ***\n");

	N1 = INTEGER_VALUE(N);
	max_con1 = INTEGER_VALUE(max_con);
	k1 = NUMERIC_VALUE(k);
	weight_mean1 = NUMERIC_VALUE(weight_mean);
	weight_sd1 = NUMERIC_VALUE(weight_sd);

	l = connectivity_random1(l, N1, max_con1, k1,weight_mean1, weight_sd1);
	ris = daLISTA(l, &nProtected);

	StrBilanciam();
	ControllaCanc();

	UNPROTECT(nProtected);

	return ris;
}
