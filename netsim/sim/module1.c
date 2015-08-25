#include "module1.h"

#define g_scalare_i globali.module1.scalare_i
#define g_Ng globali.module1.Ng
#define g_mm globali.module1.mm
#define g_tmp1_i globali.module1.tmp1_i
#define g_conn_matr1 globali.module1.conn_matr
#define g_score_matr1 globali.probmod.score_matr1
#define g_indices1 globali.module1.indices

void module12(int num1, int Nrim, const MATRICEi *Mdiscr, const VETTOREi *h, const VETTOREi *h_new, const VETTOREi *Sin, const VETTOREi *Sout, const VETTOREd *STin, const VETTOREd *STout, const VETTOREd *Freq_in, const VETTOREd *Freq_out, int max_con, double Cf_c, const VETTOREd *toll, struct RisProbM *ris)
{
	int controllo = 0, max_g1;
	double Pm1;

	_Intestazione("\n***module1***\n");
#ifdef DET
	fprintf(fp_det, "input:\n");
	fprintf(fp_det, "\tnum1 =  %d\n", num1);
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

	CREAv_i(g_scalare_i, 1);
	while (controllo == 0) {
		max_g1 = min_s_i(num1, Nrim);
		g_tmp1_i = seq_i(g_tmp1_i, 1, max_g1, 1);
		g_Ng = sampleB(g_Ng, g_tmp1_i, 1, 0);
		g_conn_matr1 = mod1(g_conn_matr1, ACCEDIv_i(g_Ng, 1), Cf_c, max_con);
		g_score_matr1 = probmod2(g_score_matr1, g_conn_matr1, h, Sin, Sout, STin, STout, Freq_in, Freq_out, toll, ris);
		// Pm1<-aus.p[[1]]
		Pm1 = ris->score;
		// qui dovrei ricopiarla, in realta` Sc e` passata per riferimento, per cui non serve farlo (anzi, sembra non essere mai usata!!)
		// Sc<-aus.p[[2]]
		// Sc = aus_p->dati[2].md; // [[2]] ?
		if (!ISNA(Pm1))
			controllo = 1;
		else {
			num1--;
			if (num1 < 1) {
				controllo = 1;
				// aus.p[[1]]<-Pm1<-NA
				Pm1 = NA_REAL;
				// già allocato in probmod1
				ris->score = Pm1;
			}
		}
	}
	//~ CANCELLAv_i(g_scalare_i);
	//~ CANCELLAv_i(g_Ng);
	//~ // non devo cancellare g_mm qui, lo farò in seguito
	//~ // CANCELLAv_i(g_mm);
	//~ CANCELLAv_i(g_tmp1_i);
	// aus_p[[5]] = c(1) #hubs
	// il 2' elemento non e` creato in probmod, ma in module*
	// il 5' elemento non e' allocato da probmod e potrebbe gia` essere stato allocato!
	// hubs = g_indices;
	CREAv_i(g_indices1, 1);
	ASSEGNAv_i(g_indices1, 1, 1);

	StrBilanciam();

#ifdef DET
	fprintf(fp_det, "module1 output:\n");
	fprintf(fp_det, "\taus.p = %.16g\n", ris->score);
	_StampaRawMatr_i(g_conn_matr1);
	_StampaRawMatr_d(g_score_matr1);
	fprintf(fp_det, " %s\n", ris->label);
	_StampaRawVett_i(g_indices1);
#endif

	return;
}

LISTA *module11(LISTA *ris, int num1, int Nrim, MATRICEi *Mdiscr, VETTOREi *h, VETTOREi *h_new, VETTOREi *Sin, VETTOREi *Sout, VETTOREd *STin, VETTOREd *STout, VETTOREd *Freq_in, VETTOREd *Freq_out, int max_con, double Cf_c, VETTOREd *toll)
{
	int controllo = 0, max_g1;
	double Pm1;
	VETTOREi *Ng = NULL, *tmp1_i = NULL, *scalare_i = NULL;
	MATRICEi *mm = NULL;
	//MATRICEd *Sc = NULL;

	_Intestazione("\n*** module11 ***\n");

	CREAv_i(scalare_i, 1);
	while (controllo == 0) {
		max_g1 = min_s_i(num1, Nrim);
		tmp1_i = seq_i(tmp1_i, 1, max_g1, 1);
		Ng = sampleB(Ng, tmp1_i, 1, 0);
		mm = mod1(mm, ACCEDIv_i(Ng, 1), Cf_c, max_con);
		ris = probmod1(ris, mm, h, Sin, Sout, STin, STout, Freq_in, Freq_out, toll);
		// Pm1<-aus.p[[1]]
		CtrlLlst(ris, 1);
		Pm1 = ACCEDIv_d(ACCEDIlst(ris, 1, vd), 1);
		// qui dovrei ricopiarla, in realta` Sc e` passata per riferimento, per cui non serve farlo (anzi, sembra non essere mai usata!!)
		// Sc<-aus.p[[2]]
		// Sc = aus_p->dati[2].md; // [[2]] ?
		if (!ISNA(Pm1))
			controllo = 1;
		else {
			num1--;
			if (num1 < 1) {
				controllo = 1;
				// aus.p[[1]]<-Pm1<-NA
				Pm1 = NA_REAL;
				// già allocato in probmod1
				CtrlLlst(ris, 1);
				ASSEGNAv_d(ACCEDIlst(ris, 1, vd), 1, Pm1);
			}
		}
	}
	CANCELLAv_i(scalare_i);
	CANCELLAv_i(Ng);
	// non devo cancellare mm, perche´ sarebbe la M di aus_p->dati[1]
	CANCELLAv_i(tmp1_i);
	// aus_p[[5]] = c(1) #hubs
	// il 5' elemento non e' allocato da probmod
	CtrlLlst(ris, 5);
	CREAv_i(ACCEDIlst(ris, 5, vi), 1);
	ASSEGNAv_i(ACCEDIlst(ris, 5, vi), 1, 1);

	StrBilanciam();

	return ris;
}

SEXP module1(SEXP num1, SEXP Nrim, SEXP Mdiscr, SEXP h, SEXP h_new, SEXP Sin, SEXP Sout, SEXP STin, SEXP STout, SEXP Freq_in, SEXP Freq_out, SEXP max_con, SEXP Cf_c, SEXP toll)
{
	int nProtected = 0;
	int num11, Nrim1, max_con1;
	double Cf_c1;
	MATRICEi *Mdiscr1;
	VETTOREi *h1, *h_new1, *Sin1, *Sout1;
	VETTOREd *STin1, *STout1, *Freq_in1, *Freq_out1, *toll1;
	LISTA *l = NULL;
	SEXP ris;

	_Intestazione("\n*** module1 ***\n");

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

	l = module11(l, num11, Nrim1, Mdiscr1, h1, h_new1, Sin1, Sout1, STin1, STout1, Freq_in1, Freq_out1, max_con1, Cf_c1, toll1);
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
