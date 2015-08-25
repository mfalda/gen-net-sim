#include "module3.h"

#define g_mm globali.module3.mm
#define g_tmp1_i globali.module3.tmp1_i
#define g_Ng globali.module3.Ng
#define g_Ng_UP globali.module3.Ng_UP
#define g_conn_matr3 globali.module3.conn_matr
#define g_score_matr3 globali.probmod.score_matr3
#define g_indices3 globali.module3.indices


void module32(int num3, int Nrim, const MATRICEi *Mdiscr, const VETTOREi *h, const VETTOREi *h_new, const VETTOREi *Sin, const VETTOREi *Sout, const VETTOREd *STin, const VETTOREd *STout, const VETTOREd *Freq_in, const VETTOREd *Freq_out, int max_con, double Cf_c, const VETTOREd *toll, struct RisProbM *ris)
{
	int controllo = 0;
	double Pm3;
	int max_g3, max_g3_UP;
	int Ng_DOWN = 0;

	_Intestazione("\n***module3***\n");
#ifdef DET
	fprintf(fp_det, "input:\n");
	fprintf(fp_det, "\tnum3 =  %d\n", num3);
	fprintf(fp_det, "\tNrim =  %d\n", Nrim);
	fprintf(fp_det, "\tMdiscr = ");
	_StampaRawMatr_i(Mdiscr);
	fprintf(fp_det, "\th = ");
	_StampaRawVett_i(h);
	fprintf(fp_det, "\th_new = ");
	_StampaRawVett_i(h_new);
	fprintf(fp_det, "\tSin = ");
	_StampaRawVett_i(Sin);
	fprintf(fp_det, "\tSout = ");
	_StampaRawVett_i(Sout);
	fprintf(fp_det, "\tSTin = ");
	_StampaRawVett_d(STin);
	fprintf(fp_det, "\tSTout = ");
	_StampaRawVett_d(STout);
	fprintf(fp_det, "\tFreq_in = ");
	_StampaRawVett_d(Freq_in);
	fprintf(fp_det, "\tFreq_out = ");
	_StampaRawVett_d(Freq_out);
	fprintf(fp_det, "\tmax_con =  %d\n", max_con);
	fprintf(fp_det, "\tCf_c =  %.16g\n", Cf_c);
	fprintf(fp_det, "\ttoll = ");
	_StampaRawVett_d(toll);
#endif

	// controllo<-0
	// while (controllo==0)
	while (controllo == 0) {
		// max.g3<-min(num3,Nrim)
		max_g3 = min_s_i(num3, Nrim);
		// max.g3.UP<-min(max.con,max.g3-1)
		max_g3_UP = min_s_i(max_con, max_g3 - 1);
		// g_Ng<-sampleB(seq(2,max.g3,1),1)
		g_tmp1_i = seq_i(g_tmp1_i, 2, max_g3, 1);
		g_Ng = sampleB(g_Ng, g_tmp1_i, 1, 0);
		// g_Ng.UP<-sampleB(seq(1,min(max.g3.UP,(g_Ng-1)),1),1)
		g_tmp1_i = seq_i(g_tmp1_i, 1, min_s_i(max_g3_UP, ACCEDIv_i(g_Ng, 1) - 1), 1);
		// replace == FALSE
		g_Ng_UP = sampleB(g_Ng_UP, g_tmp1_i, 1, 0);
		// g_Ng.DOWN<-g_Ng-g_Ng.UP
		Ng_DOWN = ACCEDIv_i(g_Ng, 1) - ACCEDIv_i(g_Ng_UP, 1);
		// g_mm<-MOD3(g_Ng.UP,g_Ng.DOWN,STin=STin,STout=STout,max.con,Cf.c)
		g_conn_matr3 = mod3(g_conn_matr3, ACCEDIv_i(g_Ng_UP, 1), Ng_DOWN, STin, STout, max_con, Cf_c);
		// aus.p<-probmod(M=g_mm,h=h,Sin=Sin,Sout=Sout,STin=STin,STout=STout,Freq.in=Freq.in,Freq.out=Freq.out,toll=toll)
		g_score_matr3 = probmod2(g_score_matr3, g_conn_matr3, h, Sin, Sout, STin, STout, Freq_in, Freq_out, toll, ris);
		// Pm3<-aus.p[[1]]
		Pm3 = ris->score;
		// qui dovrei ricopiarla, in realta` Sc e` passata per riferimento, per cui non serve farlo (anzi, sembra non essere mai usata!!)
		//Sc<-aus.p[[2]]
		// Sc = aus_p->dati[2].md; // [[2]] ?
		// if (!is.na(Pm3))   controllo<-1
		if (!ISNA(Pm3))
			controllo = 1;
		// else {
		else {
			// num3<-num3-1
			num3--;
			// if (num3<2)
			if (num3 < 2) {
				// {controllo<-1
				controllo = 1;
				// aus.p[[1]] <- Pm3 <- NA
				Pm3 = NA_REAL;
				ris->score = Pm3;
			}
		}
	}
	// il 2' elemento non e` creato in probmod, ma in module*
	// il 5' elemento non e` allocato da probmod e potrebbe gia` essere stato allocato!
	//hubs = g_indices;
	// aus.p[[5]]<-seq(g_Ng.DOWN+1,g_Ng.UP+g_Ng.DOWN,1)
	g_indices3 = seq_i(g_indices3, Ng_DOWN + 1, ACCEDIv_i(g_Ng_UP, 1) + Ng_DOWN, 1);
	// non devo cancellare g_mm qui, lo farò in seguito
	//~ // CANCELLAm_i(g_mm);
	//~ CANCELLAv_i(g_tmp1_i);
	//~ CANCELLAv_i(g_Ng);
	//~ CANCELLAv_i(g_Ng_UP);

	StrBilanciam();

#ifdef DET
	fprintf(fp_det, "module3 output:\n");
	fprintf(fp_det, "\taus.p = %.16g\n", ris->score);
	_StampaRawMatr_i(g_conn_matr3);
	_StampaRawMatr_d(g_score_matr3);
	fprintf(fp_det, " %s\n", ris->label);
	_StampaRawVett_i(g_indices3);
#endif

	return;
}

LISTA *module31(LISTA* ris, int num3, int Nrim, MATRICEi *Mdiscr, VETTOREi *h, VETTOREi *h_new, VETTOREi *Sin, VETTOREi *Sout, VETTOREd *STin, VETTOREd *STout, VETTOREd *Freq_in, VETTOREd *Freq_out, int max_con, double Cf_c, VETTOREd *toll)
{
	int controllo = 0;
	double Pm3;
	int max_g3, max_g3_UP;
	int NDOWN = 0;
	VETTOREi *NUP = NULL, *Ng = NULL, *tmp1_i = NULL;
	MATRICEi *mm = NULL;
	//MATRICEd *Sc = NULL;

	_Intestazione("\n*** module31 ***\n");

	while (controllo == 0) {
		max_g3 = min_s_i(num3, Nrim);
		max_g3_UP = min_s_i(max_con, max_g3 - 1);
		tmp1_i = seq_i(tmp1_i, 2, max_g3, 1);
		Ng = sampleB(Ng, tmp1_i, 1, 0);
		tmp1_i = seq_i(tmp1_i, 1, min_s_i(max_g3_UP, ACCEDIv_i(Ng, 1) - 1), 1);
		NUP = sampleB(NUP, tmp1_i, 1, 0);
		NDOWN = ACCEDIv_i(Ng, 1) - ACCEDIv_i(NUP, 1);

		mm = mod3(mm, ACCEDIv_i(NUP, 1), NDOWN, STin, STout, max_con, Cf_c);
		ris = probmod1(ris, mm, h, Sin, Sout, STin, STout, Freq_in, Freq_out, toll);
		CtrlLlst(ris, 1);
		Pm3 = ACCEDIv_d(ACCEDIlst(ris, 1, vd), 1);
		// qui dovrei ricopiarla, in realta` Sc e` passata per riferimento, per cui non serve farlo (anzi, sembra non essere mai usata!!)
		// Sc = aus_p->dati[2].md; // [[2]] ?
		if (!ISNA(Pm3))
			controllo = 1;
		else {
			num3--;
			if (num3 < 2) {
				controllo = 1;
				// aus.p[[1]] <- Pm3 <- NA
				Pm3 = NA_REAL;
				CtrlLlst(ris, 1);
				ASSEGNAv_d(ACCEDIlst(ris, 1, vd), 1, Pm3);
			}
		}

	}
	tmp1_i = seq_i(tmp1_i, NDOWN + 1, ACCEDIv_i(NUP, 1) + NDOWN, 1); //hubs
	// non devo cancellare mm, perche´ sarebbe la M di aus_p->dati[1]
	CANCELLAv_i(Ng);
	CANCELLAv_i(NUP);
	CtrlSlst(ris, 5);
	ASSEGNAlst(ris, 5, vi, tmp1_i);

	StrBilanciam();

	return ris;
}

SEXP module3(SEXP num1, SEXP Nrim, SEXP Mdiscr, SEXP h, SEXP h_new, SEXP Sin, SEXP Sout, SEXP STin, SEXP STout, SEXP Freq_in, SEXP Freq_out, SEXP max_con, SEXP Cf_c, SEXP toll)
{
	int nProtected = 0;
	int num11, Nrim1, max_con1;
	double Cf_c1;
	MATRICEi *Mdiscr1;
	VETTOREi *h1, *h_new1, *Sin1, *Sout1;
	VETTOREd *STin1, *STout1, *Freq_in1, *Freq_out1, *toll1;
	LISTA *l = NULL;
	SEXP ris;

	_InitDbg(false, false, false);

	_Intestazione("\n*** module3 ***\n");

	num11 = INTEGER_VALUE(num1);
	Nrim1 = INTEGER_VALUE(Nrim);
	Mdiscr1 = inMATRICE_i(Mdiscr, &nProtected);
	h1 = inVETTORE_i(h, &nProtected);
	h_new1 = inVETTORE_i(h_new, &nProtected);
	Sin1 = inVETTORE_i(Sin, &nProtected);
	Sout1 = inVETTORE_i(Sout, &nProtected);
	STin1 = inVETTORE_d(STin, &nProtected);
	STout1 = inVETTORE_d(STout, &nProtected);
	Freq_in1 = inVETTORE_d(Freq_in, &nProtected);
	Freq_out1 = inVETTORE_d(Freq_out, &nProtected);
	max_con1 = INTEGER_VALUE(max_con);
	Cf_c1 = NUMERIC_VALUE(Cf_c);
	toll1 = inVETTORE_d(toll, &nProtected);

	l = module31(l, num11, Nrim1, Mdiscr1, h1, h_new1, Sin1, Sout1, STin1, STout1, Freq_in1, Freq_out1, max_con1, Cf_c1, toll1);
	ris = daLISTA(l, &nProtected);

	CANCELLAv_i(h1);
	CANCELLAv_i(h_new1);
	CANCELLAv_i(Sin1);
	CANCELLAv_i(Sout1);
	CANCELLAm_i(Mdiscr1);
	CANCELLAv_d(STin1);
	CANCELLAv_d(Freq_in1);
	CANCELLAv_d(Freq_out1);
	CANCELLAv_d(STout1);
	CANCELLAv_d(toll1);

	StrBilanciam();
	ControllaCanc();

	UNPROTECT(nProtected);

	return ris;
}




