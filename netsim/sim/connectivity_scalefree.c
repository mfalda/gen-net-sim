#include "connectivity_scalefree.h"

#define g_M globali.connectivity_scalefree.M
#define g_Mdiscr globali.connectivity_scalefree.Mdiscr
#define g_x globali.connectivity_scalefree.x
#define g_y globali.connectivity_scalefree.y
#define g_d globali.connectivity_scalefree.d
#define g_s globali.connectivity_scalefree.s
#define g_o globali.connectivity_scalefree.o
#define g_regulatedind globali.connectivity_scalefree.regulatedind
#define g_indL globali.connectivity_scalefree.indL
#define g_Sr globali.connectivity_scalefree.Sr
#define g_inthenet globali.connectivity_scalefree.inthenet
#define g_not_inthenet globali.connectivity_scalefree.not_inthenet
#define g_not_regulated globali.connectivity_scalefree.not_regulated
#define g_numposs globali.connectivity_scalefree.numposs
#define g_give_outlink globali.connectivity_scalefree.give_outlink
#define g_aus_give_outlink globali.connectivity_scalefree.aus_give_outlink
#define g_indInf globali.connectivity_scalefree.indInf
#define g_ind_Sc globali.connectivity_scalefree.ind_Sc
#define g_a1 globali.connectivity_scalefree.a1
#define g_a2 globali.connectivity_scalefree.a2
#define g_primi globali.connectivity_scalefree.primi
#define g_indici globali.connectivity_scalefree.indici
#define g_num_v globali.connectivity_scalefree.num_v
#define g_mem_o globali.connectivity_scalefree.mem_o
#define g_available globali.connectivity_scalefree.available
#define g_campione globali.connectivity_scalefree.campione
#define g_linked globali.connectivity_scalefree.linked
#define g_ind_s globali.connectivity_scalefree.ind_s
#define g_Sout globali.connectivity_scalefree.Sout
#define g_Sin globali.connectivity_scalefree.Sin
#define g_ind globali.connectivity_scalefree.ind
#define g_indok globali.connectivity_scalefree.indok
#define g_ind1 globali.connectivity_scalefree.ind1
#define g_ind0 globali.connectivity_scalefree.ind0
#define g_tmp1_i globali.connectivity_scalefree.tmp1_i
#define g_tmp2_i globali.connectivity_scalefree.tmp2_i
#define g_tmp1_d globali.connectivity_scalefree.tmp1_d
#define g_tmp2_d globali.connectivity_scalefree.tmp2_d
#define g_scalare_i globali.connectivity_scalefree.scalare_i
#define g_scalare_d globali.connectivity_scalefree.scalare_d
#define g_Prob globali.connectivity_scalefree.Prob
#define g_Freq_in globali.connectivity_scalefree.Freq_in
#define g_Freq_out globali.connectivity_scalefree.Freq_out
#define g_STin globali.connectivity_scalefree.STin
#define g_STout globali.connectivity_scalefree.STout
#define g_p globali.connectivity_scalefree.p
#define g_toll1 globali.connectivity_scalefree.toll1
#define g_toll_in globali.connectivity_scalefree.toll_in
#define g_toll_out globali.connectivity_scalefree.toll_out
#define g_Sc globali.connectivity_scalefree.Sc
#define g_p_ind globali.connectivity_scalefree.p_ind
#define g_p_out globali.connectivity_scalefree.p_out
#define g_aus globali.connectivity_scalefree.aus

// connectivityscalefree<-function(N=50,max.con=12,gamma=2.2,r.tol=0.1,a.tol=1,  weight.mean=1, weight.sd=0.1)
LISTA *connectivity_scalefree1(LISTA *ris, int N, int max_con, double gamma, double r_tol, double a_tol, double weight_mean, double weight_sd)
{
	int i, j, h, L, ri, n, NN, m, mx, Lp, m_c, num = 0, min_out_i, L_m, L_go, mj;
	double somma = 0.0, tmp_d = 0.0;
	enum TIPO tipi[2];

	_Intestazione("\n***connectivity_scalefree***\n");
#ifdef DET
	fprintf(fp_det, "input:\n");
	fprintf(fp_det, "\tN =  %d\n", N);
	fprintf(fp_det, "\tmax_con =  %d\n", max_con);
	fprintf(fp_det, "\tgamma =  %.16g\n", gamma);
	fprintf(fp_det, "\tr_tol =  %.16g\n", r_tol);
	fprintf(fp_det, "\ta_tol =  %.16g\n", a_tol);
	fprintf(fp_det, "\tweight_mean =  %.16g\n", weight_mean);
	fprintf(fp_det, "\tweight_sd =  %.16g\n", weight_sd);
#endif

	//  g_Mdiscr<-matrix(0,ncol=N,nrow=N)
	CREAm_i(g_Mdiscr, N, N);
	InitMatr_i(g_Mdiscr, 0);
	//  if (max.con>N) max.con<-N
	if (max_con > N)
		max_con = N;
	//  g_Prob<-c(seq(1,N,1)^(-gamma),0)
	g_tmp1_d = seq_d(g_tmp1_d, 1.0, (double) N, 1.0);
	g_tmp2_d = exp_d(g_tmp2_d, g_tmp1_d, -gamma);
	CREAv_i(g_scalare_i, 1);
	// g_linked ha almeno un elemento
	CREAv_i(g_linked, 1);
	g_linked->dim = 0;
	// g_not_regulated ha almeno un elemento
	CREAv_i(g_not_regulated, 1);
	g_not_regulated->dim = 0;
	CREAv_d(g_scalare_d, 1);
	CREAv_i(g_o, 1);
	ASSEGNAv_d(g_scalare_d, 1, 0.0);
	g_Prob = vettore2v_d(g_Prob, g_tmp2_d, g_scalare_d);
	//  g_Prob<-g_Prob/(sum(g_Prob)) // cancella NA = 0
	dividi1_vs_d(g_Prob, somma_v_d(g_Prob, false));
	//  Freq.out<-Freq.in<-c(N,rep(0,N+1))
	CREAv_d(g_tmp2_d, N + 1);
	InitVett_d(g_tmp2_d, 0.0);
	ASSEGNAv_d(g_scalare_d, 1, (double) N);
	g_Freq_out = vettore2v_d(g_Freq_out, g_scalare_d, g_tmp2_d);
	g_Freq_in = copia_v_d(g_Freq_in, g_Freq_out, 1, LENGTHv_d(g_Freq_out)); // copia, come nell'instruzione di R
	//  g_STout<-g_STin<-c(NA,g_Prob*N)
	g_tmp1_d = moltiplica_vs_d(g_tmp1_d, g_Prob, (double) N);
	ASSEGNAv_d(g_scalare_d, 1, NA_REAL);
	g_STout = vettore2v_d(g_STout, g_scalare_d, g_tmp1_d);
	g_STin = copia_v_d(g_STin, g_STout, 1, LENGTHv_d(g_STout)); // copia, come nell'instruzione di R
	//  g_p<-g_Prob[1:max.con]
	g_p = segmento_v_d(g_p, g_Prob, 1, max_con);
	//  g_STin[(max.con+2):(N+2)]<-0
	assegna1_v_segm_d(g_STin, max_con + 2, N + 2, 0.0);
	//  g_STin[2:(max.con+1)]<-N*g_p/sum(g_p)
	g_tmp2_d = moltiplica_vs_d(g_tmp2_d, g_p, (double) N);
	// in somma cancNA = 0
	dividi1_vs_d(g_tmp2_d, somma_v_d(g_p, false));
	assegna1_v_segmv_d(g_STin, 2, max_con + 1, g_tmp2_d);
	//  g_STin[1]<-0
	ASSEGNAv_d(g_STin, 1, 0.0);
	//  g_aus<-cbind(g_STout*r.tol,rep(a.tol,length(g_STout)))
	CREAv_d(g_tmp1_d, LENGTHv_d(g_STout));
	InitVett_d(g_tmp1_d, a_tol);
	g_tmp2_d = moltiplica_vs_d(g_tmp2_d, g_STout, r_tol);
	g_aus = cbind2v_d(g_aus, g_tmp2_d, g_tmp1_d);
	//  toll.out<-apply(g_aus,1,max)
	g_toll_out = max_righe_d(g_toll_out, g_aus);
	//  toll.out[which(g_STout==0)]<-0		# When g_STin==0 toll is set to 0
	g_tmp1_i = which_v_indxeq_d(g_tmp1_i, g_STout, 0.0);
	assegna1_v_indx_d(g_toll_out, g_tmp1_i, 0.0);
	//  g_aus<-cbind(g_STin*r.tol,rep(a.tol,length(g_STin)))
	CREAv_d(g_tmp1_d, LENGTHv_d(g_STout));
	InitVett_d(g_tmp1_d, a_tol);
	g_tmp2_d = moltiplica_vs_d(g_tmp2_d, g_STin, r_tol);
	g_aus = cbind2v_d(g_aus, g_tmp2_d, g_tmp1_d);
	//  toll.in<-apply(g_aus,1,max)
	g_toll_in = max_righe_d(g_toll_in, g_aus);
	//  toll.in[which(g_STin==0)]<-0		# When g_STin==0 toll is set to 0
	g_tmp1_i = which_v_indxeq_d(g_tmp1_i, g_STin, 0.0);
	assegna1_v_indx_d(g_toll_in, g_tmp1_i, 0.0);
	//  g_Mdiscr[1,2]<-1
	ASSEGNAm_i(g_Mdiscr, 1, 2, 1);
	//  g_inthenet<-1
	CREAv_i(g_inthenet, 1);
	ASSEGNAv_i(g_inthenet, 1, 1);
	//  not.g_inthenet<-seq(2,N,1)
	g_not_inthenet = seq_i(g_not_inthenet, 2, N, 1);
	//  NN<-length(not.g_inthenet)
	NN = LENGTHv_i(g_not_inthenet);
	//  while (NN>0)
	while (NN > 0) {
		//    { not.regulated<-give.outlink<-vector()
		CREAv_i(g_not_regulated, 0);
		CREAv_i(g_give_outlink, 0);
		//      for (j in (1:NN))
		for (j = 1; j <= NN; j++) {
			//        {i<-not.g_inthenet[j]
			i = ACCEDIv_i(g_not_inthenet, j);
			//  	g_Sin<-apply(g_Mdiscr,1,sum)
			g_Sin = somma_righe_i(g_Sin, g_Mdiscr);
			//    	m<-max(g_Sin)+1
			m = max_v_i(g_Sin) + 1;
			//    	Freq.in[1:m]<-hist(g_Sin,breaks=seq(0,m,1),right=FALSE,plot=FALSE)$counts
			g_tmp2_i = seq_i(g_tmp2_i, 0, m, 1);
			g_tmp1_i = hist1(g_tmp1_i, g_Sin, g_tmp2_i, 0, 1, 0);
			g_tmp1_d = promuovi_i(g_tmp1_d, g_tmp1_i);
			assegna1_v_segmv_d(g_Freq_in, 1, m, g_tmp1_d);
			//    	g_Sout<-apply(g_Mdiscr,2,sum)
			g_Sout = somma_colonne_i(g_Sout, g_Mdiscr);
			//    	m<-max(g_Sout)+1
			m = max_v_i(g_Sout) + 1;
			//   	Freq.out[1:m]<-hist(g_Sout,breaks=seq(0,m,1),right=FALSE,plot=FALSE)$counts
			g_tmp2_i = seq_i(g_tmp2_i, 0, m, 1);
			g_tmp1_i = hist1(g_tmp1_i, g_Sout, g_tmp2_i, 0, 1, 0);
			g_tmp1_d = promuovi_i(g_tmp1_d, g_tmp1_i);
			assegna1_v_segmv_d(g_Freq_out, 1, m, g_tmp1_d);
			//    	mx<-length(g_inthenet)
			mx = LENGTHv_i(g_inthenet);
			//    	g_aus.in<-Freq.in[2:(mx+1)]-g_STin[2:(mx+1)]-toll.in[2:(mx+1)]
			g_tmp1_d = f_aux1_d(g_tmp1_d, g_Freq_in, g_STin, g_toll_in, 2, mx + 1, -1, -1);
			//    	g_aus.in[which(is.na(g_aus.in))]<-1
			assegna1_v_indxNA_d(g_tmp1_d, g_tmp1_d, 1.0, false);
			//    	g_numposs<-which(g_aus.in<0)
			g_numposs = which_v_indxlt_d(g_numposs, g_tmp1_d, 0.0);
			//    	Lp<-min(length(g_numposs),max.con)
			Lp = min_s_i(LENGTHv_i(g_numposs), max_con);
			// 	if (Lp>0)
			if (Lp > 0) {
				// 	  {g_aus.in<-Freq.in[2:(mx+1)]-g_STin[2:(mx+1)]+toll.in[2:(mx+1)]
				g_tmp1_d = f_aux1_d(g_tmp1_d, g_Freq_in, g_STin, g_toll_in, 2, mx + 1, -1, 1);
				//    	   g_aus.in[which(is.na(g_aus.in))]<-1
				assegna1_v_indxNA_d(g_tmp1_d, g_tmp1_d, 1.0, false);
				// 	   aus2<-g_aus.in
				g_tmp2_d = copia_v_d(g_tmp2_d, g_tmp1_d, 1, LENGTHv_d(g_tmp1_d)); // copia, come nell'instruzione di R
				// 	   g_primi<-which(aus2<0)
				// cambiato in > per tener conto della tolleranza
				g_primi = which_v_indxlt_d(g_primi, g_tmp2_d, 0.0);
				//       	   if (length(g_primi)>0)
				if (LENGTHv_i(g_primi) > 0) {
					// 	     {num<-which.min(aus2)
					num = which_v_indxmin_d(g_tmp2_d);
					// 	      if (num>max.con) num<-0
					if (num > max_con)
						num = 0;
					// 	     }
				}
				// 	    else
				else {
					// 	     {g_p<-rep(0,Lp)
					CREAv_d(g_p, Lp);
					InitVett_d(g_p, 0.0);
					//       	      for (n in (1:Lp)) g_p[n]<-Score.sf(S=0,ST=g_STin,Freq=Freq.in,n=n,toll=toll.in)
					for (n = 1; n <= Lp ; n++) {
						CREAv_i(g_tmp1_i, n);
						InitVett_i(g_tmp1_i, 0);
						g_tmp1_d = score_sf1(g_tmp1_d, g_tmp1_i, g_STin, g_Freq_in, n, g_toll_in);
						ASSEGNAv_d(g_p, n, ACCEDIv_d(g_tmp1_d, 1));
					}
					//               g_indInf<-which(g_p==(-Inf))
					g_indInf = which_v_indxeq_d(g_indInf, g_p, R_NegInf);
					//               g_p[g_indInf]<-0
					assegna1_v_indx_d(g_p, g_indInf, 0.0);
					//               if (min(g_p)<=0) g_p<-g_p-min(g_p)*(1.01)
					tmp_d = min_v_d(g_p);
					if (tmp_d <= 0.0)
						incr1_v_d(g_p, - tmp_d * 1.01);
					// g_p[g_indInf]<-0
					assegna1_v_indx_d(g_p, g_indInf, 0.0);
					// if (sum(g_p)>0) num<-sampleB(seq(1,Lp,1),1,prob=g_p/sum(g_p))
						// NA = 0
					somma = somma_v_d(g_p, false);
					if (somma > 0.0) {
						g_tmp1_i = seq_i(g_tmp1_i, 1, Lp, 1);
						g_tmp1_d = dividi_vs_d(g_tmp1_d, g_p, somma);
						// replace = 0
						g_num_v = sampleB_p(g_num_v, g_tmp1_i, 1, 0, g_tmp1_d);
						num = ACCEDIv_i(g_num_v, 1);
					}
					// 	      else num<-0
					else
						num = 0;
					//    	     }
				}
				// 	  }
			}
			// 	if (num==0) not.regulated<-c(not.regulated,i)
			if (num == 0)
				g_not_regulated = accoda1_vs_i(g_not_regulated, i);
			// 	g_linked<-vector()
			g_linked->dim = 0;
			//    	g_aus.give.outlink<-union(give.outlink,give.outlink)
			g_aus_give_outlink = elimina_doppi_i(g_aus_give_outlink, g_give_outlink);
			//    	mem.g_o<-vector()
			CREAv_i(g_mem_o, 0);
			//    	while (num>0)
			while (num > 0) {
				//      	  {if (length(g_aus.give.outlink)>0)
				if (LENGTHv_i(g_aus_give_outlink) > 0) {
					// 	     {g_o<-g_aus.give.outlink[1]
					ASSEGNAv_i(g_o, 1, ACCEDIv_i(g_aus_give_outlink, 1));
					// 	      g_aus.give.outlink<-setdiff(g_aus.give.outlink,g_o)
					setdiff1_i(g_aus_give_outlink, g_o);
					// 	      mem.g_o<-c(mem.g_o,g_o)
					g_mem_o = accoda1_vv_i(g_mem_o, g_o);
					// 	      g_Mdiscr[i,g_o]<-1
					assegna1_ms_rigaindx_i(g_Mdiscr, i, g_o, 1);
					// 	      num<-num-1
					num--;
					// 	      g_linked<-c(g_linked,g_o)
					g_linked = accoda1_vv_i(g_linked, g_o);
					// 	     }
				}
				//           else
				else {
					// 	     {g_Sout<-apply(g_Mdiscr,2,sum)
					g_Sout = somma_colonne_i(g_Sout, g_Mdiscr);
					//    	      m<-max(g_Sout)+1
					m = max_v_i(g_Sout) + 1;
					//    	      Freq.out[1:m]<-hist(g_Sout,breaks=seq(0,m,1),right=FALSE,plot=FALSE)$counts
					g_tmp2_i = seq_i(g_tmp2_i, 0, m, 1);
					g_tmp1_i = hist1(g_tmp1_i, g_Sout, g_tmp2_i, 0, 1, 0);
					g_tmp1_d = promuovi_i(g_tmp1_d, g_tmp1_i);
					assegna1_v_segmv_d(g_Freq_out, 1, m, g_tmp1_d);
					// 	      g_Sc<-Score.sf(S=g_Sout[g_inthenet],ST=g_STout,Freq=Freq.out,n=1,toll=toll.out)
					g_tmp1_i = copia_v_indx_i(g_tmp1_i, g_Sout, g_inthenet);
					g_Sc = score_sf1(g_Sc, g_tmp1_i, g_STout, g_Freq_out, 1, g_toll_out);
					//       	      g_indInf<-which(g_Sc==(-Inf))
					g_indInf = which_v_indxeq_d(g_indInf, g_Sc, R_NegInf);
					//               g_Sc[g_indInf]<-0
					assegna1_v_indx_d(g_Sc, g_indInf, 0.0);
					//               if (min(g_Sc)<=0) g_Sc<-g_Sc-min(g_Sc)*(1.01)
					tmp_d = min_v_d(g_Sc);
					if (tmp_d <= 0.0)
						incr1_v_d(g_Sc, -tmp_d * 1.01);
					//               g_Sc[g_indInf]<-0
					assegna1_v_indx_d(g_Sc, g_indInf, 0.0);
					//               g_ind<-setdiff(which(g_Sc>0),g_linked)
					g_tmp1_i = which_v_indxgt_d(g_tmp1_i, g_Sc, 0.0);
					g_ind = setdiff_i(g_ind, g_tmp1_i, g_linked);
					//               m.c<-length(g_ind)
					m_c = LENGTHv_i(g_ind);
					//               if (m.c>0)
					if (m_c > 0) {
						// 	   	{g_p.g_ind<-g_Sc[g_ind]/sum(g_Sc[g_ind]) # NA = 0
						g_tmp1_d = copia_v_indx_d(g_tmp1_d, g_Sc, g_ind);
						g_p_ind = dividi_vs_d(g_p_ind, g_tmp1_d, somma_v_d(g_tmp1_d, false));
						//             	 g_s<-sampleB(g_ind,1,prob=g_p.g_ind)
						g_s = sampleB_p(g_s, g_ind, 1, 0, g_p_ind);
						// 	    	 g_Mdiscr[i,g_s]<-1
						assegna1_ms_rigaindx_i(g_Mdiscr, i, g_s, 1);

						//             	 num<-num-1
						num--;
						// 	    	 g_linked<-c(g_linked,g_s)
						g_linked = accoda1_vv_i(g_linked, g_s);
						//            	} #end if (m.c>0)
					} //end if (m_c>0)
					//               else #(m.c==0)
					else {
						// 	        {g_aus<-Freq.out[2:(mx+1)]-g_STout[2:(mx+1)]+toll.out[2:(mx+1)]
						g_tmp1_d = f_aux1_d(g_tmp1_d, g_Freq_out, g_STout, g_toll_out, 2, mx + 1, -1, 1);
						//    	         g_aus[which(is.na(g_aus))]<-1
						assegna1_v_indxNA_d(g_tmp1_d, g_tmp1_d, 1.0, false);
						// 	         g_available<-which(g_aus<0)
						g_available = which_v_indxlt_d(g_available, g_tmp1_d, 0.0);
						// 		 if (length(g_available)>0) min.out.i<-min(g_available)
						if (LENGTHv_i(g_available) > 0)
							min_out_i = min_v_i(g_available);
						// 		 else
						else {
							// 		   {g_aus<-Freq.out[2:(mx+1)]-g_STout[2:(mx+1)]-toll.out[2:(mx+1)]
							g_tmp1_d = f_aux1_d(g_tmp1_d, g_Freq_out, g_STout, g_toll_out, 2, mx + 1, -1, -1);
							//    	            g_aus[which(is.na(g_aus))]<-1
							assegna1_v_indxNA_d(g_tmp1_d, g_tmp1_d, 1.0, false);
							// 	            g_available<-which(g_aus<0).
							g_available = which_v_indxlt_d(g_available, g_tmp1_d, 0.0);
							// 		    min.out.i<-min(g_available)
							min_out_i = min_v_i(g_available);
							// 		   }
						}
						// 		 g_campione<-vector()
						CREAv_i(g_campione, 0);
						// 	         g_indici<-inthenet
						g_indici = copia_v_i(g_indici, g_inthenet, 1, LENGTHv_i(g_inthenet)); // copia, come nell'instruzione di R
						//                  while ((length(g_campione)==0)&(length(g_indici)>0))
						while ((LENGTHv_i(g_campione) == 0) && (LENGTHv_i(g_indici) > 0)) {

							// 	     	   {g_aus.regulated<-sampleB(g_indici,1)
							g_tmp1_i = sampleB(g_tmp1_i, g_indici, 1, 0);
							//             	    g_campione<-setdiff(which(g_Mdiscr[g_aus.regulated,]==1),g_linked)
							g_tmp2_i = which_m_indxrowindxeq_i(g_tmp2_i, g_Mdiscr, g_tmp1_i, 1);
							g_campione = setdiff_i(g_campione, g_tmp2_i, g_linked);
							// 	    	    if (length(g_campione)==0) g_indici<-setdiff(g_indici,g_aus.regulated)
							if (LENGTHv_i(g_campione) == 0)
								setdiff1_i(g_indici, g_tmp1_i);
							// 	    	   }
						}
						// 		 if (length(g_campione)>0)
						if (LENGTHv_i(g_campione) > 0) {
							// 		   {g_aus.regulator<-sampleB(g_campione,1)
							g_tmp2_i = sampleB(g_tmp2_i, g_campione, 1, 0);
							//             	    g_Mdiscr[i,g_aus.regulator]<-1
							assegna1_ms_rigaindx_i(g_Mdiscr, i, g_tmp2_i, 1);
							// 	    	    g_Mdiscr[g_aus.regulated,g_aus.regulator]<-0
							assegna1_m_vv_i(g_Mdiscr, g_tmp1_i, g_tmp2_i, 0);
							// 	    	    g_Mdiscr[g_aus.regulated,i]<-1
							assegna1_ms_colindx_i(g_Mdiscr, g_tmp1_i, i, 1);
							// 		    num<-num-1
							num--;
							// 	 	    g_linked<-c(g_linked,g_aus.regulator)
							g_linked = accoda1_vv_i(g_linked, g_tmp2_i);
							// 		    if (min.out.i>1) give.outlink<-c(give.outlink,rep(i,(min.out.i-1)))
							if (min_out_i > 1) {
								CREAv_i(g_tmp1_i, min_out_i - 1);
								InitVett_i(g_tmp1_i, i);
								g_give_outlink = accoda1_vv_i(g_give_outlink, g_tmp1_i);
								// 		   }
							}
						}
						// 		  else
						else {

							// 		    {num<-0
							num = 0;
							// 		     g_Mdiscr[i,]<-0
							assegna1_ms_riga_i(g_Mdiscr, i, 0);
							// 		     not.regulated<-c(not.regulated,i)
							g_not_regulated = accoda1_vs_i(g_not_regulated, i);
							// 		    }
						}
						// 	   	} #end else #(m.c==0)
					} //end else #(m_c==0);
					// 	     }
				}
				// 	  } #end while (num>0)
			}
			// 	if (sum(g_Mdiscr[i,])>0) g_inthenet<-c(g_inthenet,i)
			if (somma_riga_i(g_Mdiscr, i) > 0)
				g_inthenet = accoda1_vs_i(g_inthenet, i);
			// 	L.m<-length(mem.g_o)
			L_m = LENGTHv_i(g_mem_o);
			//         if (L.m>0)
			if (L_m > 0) {
				// 	  {mj<-1
				mj = 1;
				//            L.go<-length(give.outlink)
				L_go = LENGTHv_i(g_give_outlink);
				//            for (h in (1:L.m))
				for (h = 1; h <= L_m; h++) {
					// 	     {g_ind<-which(give.outlink==mem.g_o[h])[1]
					g_tmp1_i = which_v_indxeq_i(g_tmp1_i, g_give_outlink, ACCEDIv_i(g_mem_o, h));
					ASSEGNAv_i(g_ind, 1, ACCEDIv_i(g_tmp1_i, 1));
					//       	      g_indok<-setdiff(seq(1,L.go,1),g_ind)
					g_tmp1_i = seq_i(g_tmp1_i, 1, L_go, 1);
					g_indok = setdiff_i(g_indok, g_tmp1_i, g_ind);
					//               give.outlink<-give.outlink[g_indok]
					g_tmp1_i = copia_v_indx_i(g_tmp1_i, g_give_outlink, g_indok);
					g_give_outlink = copia_v_i(g_give_outlink, g_tmp1_i, 1, LENGTHv_i(g_tmp1_i)); // copia, come nell'instruzione di R
					// 	      L.go<-L.go-1
					L_go--;
					// 	     }
				}
				//           }
			}
			//        } #end for (j in (1:NN))
		} //end for (j in (1:NN));
		//     if (length(not.regulated)>0)
		if (LENGTHv_i(g_not_regulated) > 0) {
			//      {if (length(intersect(not.g_inthenet,not.regulated))==NN)
			g_tmp1_i = interseca_i(g_tmp1_i, g_not_inthenet, g_not_regulated);
			if (LENGTHv_i(g_tmp1_i) == NN) {
				// 	{#cat("\n The algorithm converged to a topology in which in-degree and out-degree constraint are not compatible. \n tollerance constraints are relaxed.\n")
				//cat("\n The algorithm converged to a topology in which in-degree and out-degree constraint are not compatible_ \n tollerance constraints are relaxed_\n");
				// 	 g_toll1<-which(toll.in!=0)
				g_toll1 = which_v_indxne_d(g_toll1, g_toll_in, 0.0);
				//          toll.in[g_toll1]<-toll.in[g_toll1]+1
				incr1_v_indx_d(g_toll_in, g_toll1, 1.0);
				// 	 g_toll1<-which(toll.out!=0)
				g_toll1 = which_v_indxne_d(g_toll1, g_toll_out, 0.0);
				//          toll.out[g_toll1]<-toll.out[g_toll1]+1
				incr1_v_indx_d(g_toll_out, g_toll1, 1.0);
				// 	}
			}
			//      }
		}
		// not.g_inthenet<-not.regulated
		g_not_inthenet = copia_v_i(g_not_inthenet, g_not_regulated, 1, LENGTHv_i(g_not_regulated));
		//     NN<-length(not.g_inthenet)
		NN = LENGTHv_i(g_not_inthenet);
		//    } #END while (NN>0)
	} //END while (NN>0);
	// ###check that every gene has at least 1 regulator
	//##check that every gene has at least 1 regulator;
	// g_Sr<-apply(g_Mdiscr,1,sum)
	g_Sr = somma_righe_i(g_Sr, g_Mdiscr);
	// g_ind<-which(g_Sr==0)
	g_ind = which_v_indxeq_i(g_ind, g_Sr, 0);
	// L<-length(g_ind)
	L = LENGTHv_i(g_ind);
	// if (L>0)
	if (L > 0) {
		//  {for (i in (1:L))
		for (i = 1; i <= L ; i++) {
			//    {ri<-g_ind[i]
			ri = ACCEDIv_i(g_ind, i);
			//     num<-1
			num = 1;
			//     g_Sout<-apply(g_Mdiscr,2,sum)
			g_Sout = somma_colonne_i(g_Sout, g_Mdiscr);
			// g_Sc<-Score(S=g_Sout,ST=g_STout,Freq=Freq.out,n=1,toll)
			g_Sc = score1(g_Sc, g_Sout, g_STout, g_Freq_out, 1, g_toll_out);
			//     g_ind.g_Sc<-which(g_Sc==(-Inf))
			g_ind_Sc = which_v_indxeq_d(g_ind_Sc, g_Sc, R_NegInf);
			//     g_Sc[g_ind.g_Sc]<-0
			assegna1_v_indx_d(g_Sc, g_ind_Sc, 0.0);
			//     if (length(which(g_Sc<=0))>0) g_Sc<-g_Sc-min(g_Sc)+1/(N^2)
			g_tmp1_i = which_v_indxle_d(g_tmp1_i, g_Sc, 0.0);
			if (LENGTHv_i(g_tmp1_i) > 0)
				incr1_v_d(g_Sc, -min_v_d(g_Sc) + (double) 1.0 / (N * N));
			//     g_Sc[g_ind.g_Sc]<-0
			assegna1_v_indx_d(g_Sc, g_ind_Sc, 0.0);
			//     if (sum(g_Sc)==0) # NA = 0
			if (Uguale(somma_v_d(g_Sc, false), 0.0)) {
				// 	{#cat("\n WARNING: tolerance parameters were relaxed to allow each node to have at least 1 regulator")
				//cat("\n WARNING: tolerance parameters were relaxed to allow each node to have at least 1 regulator");
				// 	 g_Sc<-Score(S=g_Sout,ST=g_STout,Freq=Freq.out,n=1,toll=rep(Inf,N))
				//CREAv_d(g_tmp1_d, N);
				InitVett_d(g_tmp1_d, R_PosInf);
				g_Sc = score1(g_Sc, g_Sout, g_STout, g_Freq_out, 1, g_tmp1_d);
				// 	 g_ind.g_Sc<-which(g_Sc==(-Inf))
				g_ind_Sc = which_v_indxeq_d(g_ind_Sc, g_Sc, R_NegInf);
				//   	 g_Sc[g_ind.g_Sc]<-0
				assegna1_v_indx_d(g_Sc, g_ind_Sc, 0.0);
				//   	 if (length(which(g_Sc<=0))>0) g_Sc<-g_Sc-min(g_Sc)+1/(N^2)
				g_tmp1_i = which_v_indxle_d(g_tmp1_i, g_Sc, 0.0);
				if (LENGTHv_i(g_tmp1_i) > 0)
					incr1_v_d(g_Sc, -min_v_d(g_Sc) + (double) 1.0 / (N * N));
				//   	 g_Sc[g_ind.g_Sc]<-0
				assegna1_v_indx_d(g_Sc, g_ind_Sc, 0.0);
				// 	}
			}
			//     g_p.out<-g_Sc/sum(g_Sc) #NA = 0
			g_p_out = dividi_vs_d(g_p_out, g_Sc, somma_v_d(g_Sc, false));
			//     g_ind.g_s<-sampleB(seq(1,N,1),num,prob=g_p.out)
			g_tmp1_i = seq_i(g_tmp1_i, 1, N, 1);
			g_ind_s = sampleB_p(g_ind_s, g_tmp1_i, num, 0, g_p_out);
			//     g_Mdiscr[ri,g_ind.g_s]<-1
			assegna1_ms_rigaindx_i(g_Mdiscr, ri, g_ind_s, 1);
			//     a1<-g_Sout[g_ind.g_s]
			g_a1 = copia_v_indx_i(g_a1, g_Sout, g_ind_s);
			//     a2<-a1+1
			g_a2 = somma_vs_i(g_a2, g_a1, 1);
			incr1_v_i(g_a1, 1);
			incr1_v_indx_d(g_Freq_out, g_a1, -1.0);
			//     Freq.out[a2+1]<-Freq.out[a2+1]+1
			incr1_v_i(g_a2, 1);
			incr1_v_indx_d(g_Freq_out, g_a2, 1.0);
			//    }
		}
		//  }
	}
	// #######assign weights to the connectivity matrix##########
	//######assign weights to the connectivity matrix##########;
	// g_ind<-which(g_Mdiscr==1,arr.g_ind=TRUE)
	g_ind = which_m_indxeq_i(g_ind, g_Mdiscr, 1);
	// L<-dim(g_ind)[1]
	L = LENGTHv_i(g_ind);
	// g_aus<-abs(rnorm_s(L,weight.mean, weight.sd))
	CREAv_d(g_tmp1_d, L);
	g_tmp1_d = rnorm_s(g_tmp1_d, L, weight_mean, weight_sd, "csf");
	g_tmp2_d = abs_v_d(g_tmp2_d, g_tmp1_d);
	// g_M<-matrix(0,ncol=N,nrow=N)
	CREAm_d(g_M, N, N);
	InitMatr_d(g_M, 0.0);
	// g_M[g_ind]<-aus
	assegna1_mv_indx_d(g_M, g_ind, g_tmp2_d);
	PutRNGstate();
	// return(list(g_M,g_Mdiscr))
	tipi[0] = MATRd;
	tipi[1] = MATRi;
	CreaLISTA(ris, tipi, 2);
	ris->dati[0].md = g_M;
	ris->dati[1].mi = g_Mdiscr;
	//~ CANCELLAv_d(g_x);
	//~ CANCELLAv_d(g_y);
	//~ CANCELLAv_d(g_d);
	//~ CANCELLAv_i(g_s);
	//~ CANCELLAv_i(g_o);
	//~ CANCELLAv_i(g_regulatedind);
	//~ CANCELLAv_i(g_indL);
	//~ CANCELLAv_i(g_Sr);
	//~ CANCELLAv_i(g_inthenet);
	//~ CANCELLAv_i(g_not_inthenet);
	//~ CANCELLAv_i(g_not_regulated);
	//~ CANCELLAv_i(g_numposs);
	//~ CANCELLAv_i(g_give_outlink);
	//~ CANCELLAv_i(g_aus_give_outlink);
	//~ CANCELLAv_i(g_indInf);
	//~ CANCELLAv_i(g_primi);
	//~ CANCELLAv_i(g_indici);
	//~ CANCELLAv_i(g_num_v);
	//~ CANCELLAv_i(g_mem_o);
	//~ CANCELLAv_i(g_available);
	//~ CANCELLAv_i(g_campione);
	//~ CANCELLAv_i(g_linked);
	//~ CANCELLAv_i(g_ind_s);
	//~ CANCELLAv_i(g_Sout);
	//~ CANCELLAv_i(g_Sin);
	//~ CANCELLAv_i(g_ind);
	//~ CANCELLAv_i(g_indok);
	//~ CANCELLAv_i(g_ind1);
	//~ CANCELLAv_i(g_ind0);
	//~ CANCELLAv_i(g_tmp1_i);
	//~ CANCELLAv_i(g_tmp2_i);
	//~ CANCELLAv_d(g_tmp1_d);
	//~ CANCELLAv_d(g_tmp2_d);
	//~ CANCELLAv_i(g_scalare_i);
	//~ CANCELLAv_d(g_scalare_d);
	//~ CANCELLAv_d(g_Prob);
	//~ CANCELLAv_d(g_Freq_in);
	//~ CANCELLAv_d(g_Freq_out);
	//~ CANCELLAv_d(g_STin);
	//~ CANCELLAv_d(g_STout);
	//~ CANCELLAv_d(g_p);
	//~ CANCELLAv_i(g_toll1);
	//~ CANCELLAv_d(g_toll_in);
	//~ CANCELLAv_d(g_toll_out);
	//~ CANCELLAv_d(g_Sc);
	//~ CANCELLAv_d(g_p_ind);
	//~ CANCELLAm_d(g_aus);

	StrBilanciam();

#ifdef DET
	fprintf(fp_det, "csf output:\n");
	fprintf(fp_det, "\tlist(M,Mdiscr) = ");
	_StampaRawMatr_d(g_M);
	_StampaRawMatr_i(g_Mdiscr);
#endif

	return ris;
}

SEXP connectivity_scalefree(SEXP N, SEXP max_con, SEXP gamma, SEXP r_tol, SEXP a_tol, SEXP weight_mean, SEXP weight_sd)
{
	int nProtected = 0;
	int N1, max_con1;
	double gamma1, r_tol1, a_tol1, weight_mean1, weight_sd1;
	LISTA *l = NULL;
	SEXP ris;

	_InitDbg(false, false, false);

	_Intestazione("\n*** connectivity_scalefree ***\n");

	N1 = INTEGER_VALUE(N);
	max_con1 = INTEGER_VALUE(max_con);
	gamma1 = NUMERIC_VALUE(gamma);
	r_tol1 = NUMERIC_VALUE(r_tol);
	a_tol1 = NUMERIC_VALUE(a_tol);
	weight_mean1 = NUMERIC_VALUE(weight_mean);
	weight_sd1 = NUMERIC_VALUE(weight_sd);

	l = connectivity_scalefree1(l, N1, max_con1, gamma1, r_tol1, a_tol1, weight_mean1, weight_sd1);
	ris = daLISTA(l, &nProtected);

	StrBilanciam();
	ControllaCanc();

	UNPROTECT(nProtected);

	return ris;
}
