#include "boole_result.h"

#define g_valr globali.boole_result.valr
#define g_signr globali.boole_result.signr
#define g_tmp_r globali.boole_result.tmp_r

VETTOREd *boole_result(VETTOREd *ris, const VETTOREi *r, const VETTOREd *n, const MATRICEd *M, const MATRICEi *N, int regulated)
{
	int i, l = 0;
	double expon;
#ifdef MDEBUG
	GString *tmp = NULL;
#endif

	_Intestazione("\n***boole_result***\n");
#ifdef DET
	fprintf(fp_det, "input:\n");
	fprintf(fp_det, "\tr = ");
	_StampaRawVett_i(r);
	fprintf(fp_det, "\tn = ");
	_StampaRawVett_d(n);
	fprintf(fp_det, "\tM = ");
	_StampaRawMatr_d(M);
	fprintf(fp_det, "\tN = ");
	_StampaRawMatr_i(N);
	fprintf(fp_det, "\tregulated =  %d\n", regulated);
#endif

	// if (is.nan(r[1])) {signr[1]<-1; valr[1]<-0 }
	if (r == NULL) {
		ris = vettore2s_d(ris, 1.0, 0.0);
		StrBilanciam();
		return ris;
	}
	else {
		// l<-length(r)
		l = LENGTHv_i(r);
		// devo COPIARE r, dato che lo modifichero`!
		g_tmp_r = copia_v_i(g_tmp_r, r, 1, l);
		CREAv_d(g_valr, l);
		CREAv_i(g_signr, l);
		// signr<-valr<-rep(0,l)
		InitVett_d(g_valr, 0.0);
		InitVett_i(g_signr, 0);
		//  {for(i in 1:l)
		for (i = 1; i <= l; i++) {
			//    {if(r[i]>0) {expon<-1/abs(M[regulated,r[i]]); valr[i]<-n[r[i]]^expon; signr[i]<-sign(N[regulated,r[i]])}}
			if (ACCEDIv_i(g_tmp_r, i) > 0) {
#ifdef MDEBUG
				if (Uguale(ACCEDIm_d(M, regulated, ACCEDIv_i(g_tmp_r, i)), 0.0)) {
					CREAstr(tmp, "");
					g_string_printf(tmp, "ATTENZIONE (boole_result.c, linea 42): divisione per zero!\n");
					warning(tmp->str);
					fprintf(fp_fdbg, tmp->str);
					CANCELLAstr(tmp);
				}
#endif
				expon = (double) 1.0 / fabs(ACCEDIm_d(M, regulated, ACCEDIv_i(g_tmp_r, i)));
				ASSEGNAv_d(g_valr, i, pow(ACCEDIv_d(n, ACCEDIv_i(g_tmp_r, i)), expon));
				ASSEGNAv_i(g_signr, i, Segno(ACCEDIm_i(N, regulated, ACCEDIv_i(g_tmp_r, i))));
			}
		}
		//   while (next.op(r)!=-1){
		while ((i = next_op(g_tmp_r)) && i != -1) {
			// 	i<-next.op(r)
			// 	if(r[i]==-3){	#OR
			if (ACCEDIv_i(g_tmp_r, i) == -3) {	// OR
				// 		if (signr[2*i]==signr[2*i+1]) {valr[i]<-min(1,sum(valr[2*i],valr[2*i+1])); signr[i]<-signr[2*i]}
				if (ACCEDIv_i(g_signr, 2 * i) == ACCEDIv_i(g_signr, 2 * i + 1)) {
					ASSEGNAv_d(g_valr, i, min_s_d(1.0, ACCEDIv_d(g_valr, 2 * i) + ACCEDIv_d(g_valr, 2 * i + 1)));
					ASSEGNAv_i(g_signr, i, ACCEDIv_i(g_signr, 2 * i));
				}
				// 		 else  {valr[i]<-min(1,1+signr[2*i]*valr[2*i]+signr[2*i+1]*valr[2*i+1]); signr[i]<-1}
				else {
					ASSEGNAv_d(g_valr, i, min_s_d(1.0, 1 + ACCEDIv_i(g_signr, 2 * i) * ACCEDIv_d(g_valr, 2 * i) + ACCEDIv_i(g_signr, 2 * i + 1) * ACCEDIv_d(g_valr, 2 * i + 1)));
					ASSEGNAv_i(g_signr, i, 1);
				}
				//   	r->dati[i]<-0
				ASSEGNAv_i(g_tmp_r, i, 0);
				// 		r->dati[2*i]<-(-1)
				ASSEGNAv_i(g_tmp_r, 2 * i, -1);
				// 		r->dati[2*i+1]<-(-1)
				ASSEGNAv_i(g_tmp_r, 2 * i + 1, -1);
				// 	}
			}
			// 	if(r->dati[i]==-2){	#AND
			if (ACCEDIv_i(g_tmp_r, i) == -2) {	// AND
				// 		if (signr->dati[2*i]==signr->dati[2*i+1]) {valr->dati[i]<-min(valr->dati[2*i],valr->dati[2*i+1]); signr->dati[i]<-signr->dati[2*i]}
				if (ACCEDIv_i(g_signr, 2 * i) == ACCEDIv_i(g_signr, 2 * i + 1)) {
					ASSEGNAv_d(g_valr, i, min_s_d(ACCEDIv_d(g_valr, 2 * i), ACCEDIv_d(g_valr, 2 * i + 1)));
					ASSEGNAv_i(g_signr, i, ACCEDIv_i(g_signr, 2 * i));
				}
				// 		 else  {valr->dati[i]<-max(0,signr->dati[2*i]*valr->dati[2*i]+signr->dati[2*i+1]*valr->dati[2*i+1]); signr->dati[i]<-1}
				else {
					ASSEGNAv_d(g_valr, i, max_s_d(0.0, ACCEDIv_i(g_signr, 2 * i) * ACCEDIv_d(g_valr, 2 * i) + ACCEDIv_i(g_signr, 2 * i + 1) * ACCEDIv_d(g_valr, 2 * i + 1)));
					ASSEGNAv_i(g_signr, i, 1);
				}
				// 		r->dati[i]<-0
				ASSEGNAv_i(g_tmp_r, i, 0);
				// 		r->dati[2*i]<-(-1)
				ASSEGNAv_i(g_tmp_r, 2 * i, -1);
				// 		r->dati[2*i+1]<-(-1)
				ASSEGNAv_i(g_tmp_r, 2 * i + 1, -1);
				// 	}
			}
			//    }#end while
		} // end while
		//   }
	}
	// return(c(signr->dati[1],valr->dati[1]))
	ris = vettore2s_d(ris, (double) ACCEDIv_i(g_signr, 1), ACCEDIv_d(g_valr, 1));

	//~ CANCELLAv_d(valr);
	//~ CANCELLAv_i(signr);

	StrBilanciam();

#ifdef DET
	fprintf(fp_det, "boole_result output:\n");
	fprintf(fp_det, "\tc(signr[1],valr[1]) = ");
	_StampaRawVett_d(ris);
#endif

	return ris;
	// }
}
//
