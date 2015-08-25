#include "connectivity_modular.h"

#define g_M globali.connectivity_modular.M
#define g_Mdiscr globali.connectivity_modular.Mdiscr
#define g_scalare_i globali.connectivity_modular.scalare_i
#define g_tmp1_i globali.connectivity_modular.tmp1_i
#define g_tmp2_i globali.connectivity_modular.tmp2_i
#define g_tmp1_d globali.connectivity_modular.tmp1_d
#define g_tmp2_d globali.connectivity_modular.tmp2_d
#define g_Prob globali.connectivity_modular.Prob
#define g_Freq_out globali.connectivity_modular.Freq_out
#define g_Freq_in globali.connectivity_modular.Freq_in
#define g_STout globali.connectivity_modular.STout
#define g_STin globali.connectivity_modular.STin
#define g_toll globali.connectivity_modular.toll
#define g_p_out globali.connectivity_modular.p_out
#define g_scalare_d globali.connectivity_modular.scalare_d
#define g_Sr globali.connectivity_modular.Sr
#define g_p globali.connectivity_modular.p
#define g_h globali.connectivity_modular.h
#define g_prob_mod globali.connectivity_modular.prob_mod
#define g_aus globali.connectivity_modular.aus
#define g_Sc_v globali.connectivity_modular.Sc_v
#define g_Cg globali.connectivity_modular.Cg
#define g_Sin globali.connectivity_modular.Sin
#define g_Sout globali.connectivity_modular.Sout
#define g_a1 globali.connectivity_modular.a1
#define g_a2 globali.connectivity_modular.a2
#define g_ind globali.connectivity_modular.ind
#define g_ind_Sc globali.connectivity_modular.ind_Sc
#define g_h_new globali.connectivity_modular.h_new
#define g_mod_type globali.connectivity_modular.mod_type
#define g_ind_s globali.connectivity_modular.ind_s
#define g_aus0 globali.connectivity_modular.aus0
#define g_conn_matr1 globali.module1.conn_matr
#define g_score_matr1 globali.probmod.score_matr1
#define g_indices1 globali.module1.indices
#define g_conn_matr2 globali.module2.conn_matr
#define g_score_matr2 globali.probmod.score_matr2
#define g_indices2 globali.module2.indices
#define g_conn_matr3 globali.module3.conn_matr
#define g_score_matr3 globali.probmod.score_matr3
#define g_indices3 globali.module3.indices

// connectivitymodular<-function(N=50,max.con=12,gamma=2.2,INdegree=c("free","out"), Cf.cl=0.4, num.subnet=c(5,5,10),r.tol=0.1,a.tol=1, weight.mean=1, weight.sd=0.1)
LISTA *connectivity_modular1(LISTA *ris, int N, int max_con, double gamma, GString *INdegree, double Cf_cl, VETTOREi *num_subnet, double r_tol, double a_tol, double weight_mean, double weight_sd)
{
	int num1, num2, num3, it, Nrim, controllo, m, Ng, i, L, ri, num;
	double Cf_c, Pm1, Pm2, Pm3, CC, CCi, mn;
	const char IN[3] = "in";
	const char OUT[4] = "out";
	enum TIPO tipi[2];
	// queste servono da puntatori
	MATRICEi *M1 = NULL;
	MATRICEd *Sc_m = NULL;
	VETTOREi *hubs = NULL;
	struct RisProbM pm_ris1, pm_ris2, pm_ris3;

	_Intestazione("\n***connectivity_modular***\n");
#ifdef DET
	fprintf(fp_det, "input:\n");
	fprintf(fp_det, "\tN =  %d\n", N);
	fprintf(fp_det, "\tmax_con =  %d\n", max_con);
	fprintf(fp_det, "\tgamma =  %.16g\n", gamma);
	fprintf(fp_det, "\tINdegree =  %s\n", INdegree->str);
	fprintf(fp_det, "\tCf_cl =  %.16g\n", Cf_cl);
	fprintf(fp_det, "\tnum_subnet = ");
	_StampaRawVett_i(num_subnet);
	fprintf(fp_det, "\tr_tol =  %.16g\n", r_tol);
	fprintf(fp_det, "\ta_tol =  %.16g\n", a_tol);
	fprintf(fp_det, "\tweight_mean =  %.16g\n", weight_mean);
	fprintf(fp_det, "\tweight_sd =  %.16g\n", weight_sd);
#endif

	//~ g_ris1.conn_matr = NULL;
	//~ g_ris1.score_matr = NULL;
	//~ g_ris1.indices = NULL;
	//~ CREAv_i(g_indices, 1);
	//~ g_ris2.conn_matr = NULL;
	//~ g_ris2.score_matr = NULL;
	//~ g_ris2.indices = NULL;
	//~ CREAv_i(g_ris2.indices, 1);
	//~ g_ris3.conn_matr = NULL;
	//~ g_ris3.score_matr = NULL;
	//~ g_ris3.indices = NULL;
	//~ CREAv_i(g_ris3.indices, 1);
	CREAm_i(g_Mdiscr, N, N);
	CREAv_i(g_scalare_i, 1);
	CREAv_d(g_scalare_d, 1);
	InitMatr_i(g_Mdiscr, 0);
	//  if (max.con>N) max.con<-N
	if (max_con > N)
		max_con = N;
	//  g_Prob<-c(seq(1,N,1)^(-gamma),0)
	g_tmp1_d = seq_d(g_tmp1_d, 1.0, (double) N, 1.0);
	g_tmp2_d = exp_d(g_tmp2_d, g_tmp1_d, -gamma);
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
	//  g_aus<-cbind(g_STout*r.tol,rep(a.tol,length(g_STout)))
	CREAv_d(g_tmp1_d, LENGTHv_d(g_STout));
	InitVett_d(g_tmp1_d, a_tol);
	g_tmp2_d = moltiplica_vs_d(g_tmp2_d, g_STout, r_tol);
	g_aus = cbind2v_d(g_aus, g_tmp2_d, g_tmp1_d);
	//  g_toll<-apply(g_aus,1,max)
	g_toll = max_righe_d(g_toll, g_aus);
	//  g_STin[(max.con+2):(N+2)]<-0
	assegna1_v_segm_d(g_STin, max_con + 2, N + 2, 0.0);
	//  if (INdegree=="out") g_STin[2:(max.con+1)]<-N*g_Prob[1:max.con]/sum(g_Prob[1:max.con])
	if (!strncmp(INdegree->str, OUT, 3)) {
		g_tmp1_d = segmento_v_d(g_tmp1_d, g_Prob, 1, max_con);
		g_tmp2_d = moltiplica_vs_d(g_tmp2_d, g_tmp1_d, (double) N);
		// in somma cancNA = 0
		dividi1_vs_d(g_tmp2_d, somma_v_d(g_tmp1_d, false));
		assegna1_v_segmv_d(g_STin, 2, max_con + 1, g_tmp2_d);
		//   else g_STin[2:(max.con+1)]<-NA
	} else
		assegna1_v_segm_d(g_STin, 2, max_con + 1, NA_REAL);
	//  num1<-num.subnet[1]
	num1 = ACCEDIv_i(num_subnet, 1);
	//  num2<-num.subnet[2]
	num2 = ACCEDIv_i(num_subnet, 2);
	//  num3<-num.subnet[3]
	num3 = ACCEDIv_i(num_subnet, 3);
	//  g_h<-seq(1,N,1)
	g_h = seq_i(g_h, 1, N, 1);
	//  g_h.new<-vector()
	CREAv_i(g_h_new, 0);
	//  Cf.c<-Cf.cl
	Cf_c = Cf_cl;
	//  it<-0
	it = 0;
	//  Pm1<-Pm2<-Pm3<-1/3
	Pm1 = Pm2 = Pm3 = (double) 1.0 / 3.0;
	//  while (length(g_h)>1)
	while (LENGTHv_i(g_h) > 1) {
		//    {Nrim<-length(g_h)
		Nrim = LENGTHv_i(g_h);
		//     g_Sin<-apply(Mdiscr,1,sum)
		g_Sin = somma_righe_i(g_Sin, g_Mdiscr);
		//     g_Sout<-apply(Mdiscr,2,sum)
		g_Sout = somma_colonne_i(g_Sout, g_Mdiscr);
		//     if (Cf.c>0) CC<-cluster.coeff(Mdiscr)[[1]]    else CC<-0
		if (Cf_c > 0.0) {
			g_Cg = cluster_coeff2(g_Cg, g_Mdiscr, &CC);
		} else
			CC = 0.0;
		//     if (CC>=Cf.c) CCi<-0 else CCi<-Cf.c
		if (CC >= Cf_c)
			CCi = 0.0;
		else
			CCi = Cf_c;
		//     controllo<-1
		controllo = 1;
		//     while (controllo==1)
		while (controllo == 1) {
			//      {if (!is.na(Pm1))  {aus1<-MODULE1(num1,Nrim,Mdiscr,g_h,g_h.new,g_Sin,g_Sout,g_STin, g_STout,Freq.in, Freq.out, max.con, CCi,g_toll);Pm1<-aus1[[1]]}
			if (Pm1 != NA_REAL) {
				module12(num1, Nrim, g_Mdiscr, g_h, g_h_new, g_Sin, g_Sout, g_STin, g_STout, g_Freq_in, g_Freq_out, max_con, CCi, g_toll, &pm_ris1);
				Pm1 = pm_ris1.score;
			}
			//       if (!is.na(Pm2))  {aus2<-MODULE2(num2,Nrim,Mdiscr,g_h,g_h.new,g_Sin,g_Sout,g_STin, g_STout,Freq.in, Freq.out, max.con, CCi,g_toll);Pm2<-aus2[[1]]}
			if (Pm2 != NA_REAL) {
				module22(num2, Nrim, g_Mdiscr, g_h, g_h_new, g_Sin, g_Sout, g_STin, g_STout, g_Freq_in, g_Freq_out, max_con, CCi, g_toll, &pm_ris2);
				Pm2 = pm_ris2.score;
			}
			//       if (!is.na(Pm3))  {aus3<-MODULE3(num3,Nrim,Mdiscr,g_h,g_h.new,g_Sin,g_Sout,g_STin, g_STout,Freq.in, Freq.out, max.con, CCi,g_toll);Pm3<-aus3[[1]]}
			if (Pm3 != NA_REAL) {
				module32(num3, Nrim, g_Mdiscr, g_h, g_h_new, g_Sin, g_Sout, g_STin, g_STout, g_Freq_in, g_Freq_out, max_con, CCi, g_toll, &pm_ris3);
				Pm3 = pm_ris3.score;
			}
			//       prob.mod<-c(Pm1,Pm2,Pm3)
			g_prob_mod = vettore3s_d(g_prob_mod, Pm1, Pm2, Pm3);
			//       prob.mod[which(prob.mod==0)]<-NA
			assegna1_v_indxeq_d(g_prob_mod, g_prob_mod, 0.0, NA_REAL);
			//       g_ind<-which(!is.na(prob.mod))
			g_ind = which_v_indxNA_d(g_ind, g_prob_mod, true);
			//       if (length(g_ind>0))
			if (LENGTHv_i(g_ind) > 0) {
				//         {controllo<-0
				controllo = 0;
				// 	 mn<-min(prob.mod[g_ind])
				g_tmp1_d = copia_v_indx_d(g_tmp1_d, g_prob_mod, g_ind);
				mn = min_v_d(g_tmp1_d);
				//          if (mn<=0) prob.mod<-prob.mod-mn+1/9
				if (mn <= 0) {
					somma1_vs_d(g_prob_mod, -mn + (double) 1.0 / 9.0);
				}
				//          prob.mod[which(is.na(prob.mod))]<-0
				assegna1_v_indxNA_d(g_prob_mod, g_prob_mod, 0.0, false);
				// 	 prob.mod<-prob.mod/sum(prob.mod) # cancNA = 0
				dividi1_vs_d(g_prob_mod, somma_v_d(g_prob_mod, false));
				//          mod.type<-sampleB(c(1,2,3),1,prob=prob.mod)
				g_tmp1_i = vettore3s_i(g_tmp1_i, 1, 2, 3);
				// ripet = 0
				g_mod_type = sampleB_p(g_mod_type, g_tmp1_i, 1, 0, g_prob_mod);
				//          if (mod.type==1) {M<-aus1[[2]]; Sc<-aus1[[3]];  hubs<-aus1[[5]]}
				// prima di riassegnare i puntatori cancello i dati precedenti
				// CANCELLAm_d(Sc_m); questa non devo cancellarla, perche´ la riuso
				// CANCELLAv_i(hubs);
				// questi "if" potrei anche toglierli...
				if (ACCEDIv_i(g_mod_type, 1) == 1) {
					M1 = g_conn_matr1;
					Sc_m = g_score_matr1;
					hubs = g_indices1;
				}
				// 	 if (mod.type==2) {M<-aus2[[2]]; Sc<-aus2[[3]];  hubs<-aus2[[5]]}
				else if (ACCEDIv_i(g_mod_type, 1) == 2) {
					M1 = g_conn_matr2;
					Sc_m = g_score_matr2;
					hubs = g_indices2;
				}
				// 	 if (mod.type==3) {M<-aus3[[2]]; Sc<-aus3[[3]];  hubs<-aus3[[5]]}
				else if (ACCEDIv_i(g_mod_type, 1) == 3) {
					M1 = g_conn_matr3;
					Sc_m = g_score_matr3;
					hubs = g_indices3;
				}
				// 	}
			}
			//       else
			else {

				// 	{Pm1<-Pm2<-Pm3<-1/3
				Pm1 = Pm2 = Pm3 = (double) 1.0 / 3.0;
				// 	 if (CCi>0)  {CCi<-max(0,CCi-0.1); Cf.c<-max(0,Cf.c-0.1)}
				if (CCi > 0.0) {
					CCi = max_s_d(0.0, CCi - 0.1);
					Cf_c = max_s_d(0.0, Cf_c - 0.1);
				}
				//           else
				else {
					// 	    {if ( (aus1[[4]]=="in")&(aus2[[4]]=="in")&(aus3[[4]]=="in")&(max.con<N))
					if (!strncmp(pm_ris1.label, IN, 2) && !strncmp(pm_ris2.label, IN, 2) && !strncmp(pm_ris3.label, IN, 2) && max_con < N) {
						// 		{max.con<-max.con+1
						max_con++;
						// 		 if (INdegree=="out") g_STin[2:(max.con+1)]<-N*g_Prob[1:max.con]/sum(g_Prob[1:max.con])
						if (!strncmp(INdegree->str, OUT, 3)) {
							g_tmp1_d = segmento_v_d(g_tmp1_d, g_Prob, 1, max_con);
							g_tmp2_d = moltiplica_vs_d(g_tmp2_d, g_tmp1_d, (double) N);
							// in somma cancNA = 0
							dividi1_vs_d(g_tmp2_d, somma_v_d(g_tmp1_d, false));
							assegna1_v_segmv_d(g_STin, 2, max_con + 1, g_tmp2_d);
						}
						//   		  else g_STin[2:(max.con+1)]<-NA
						else
							assegna1_v_segm_d(g_STin, 2, max_con + 1, NA_REAL);
						// 		 Cf.c<-Cf.cl
						Cf_c = Cf_cl;
						// 		 if (CC>=Cf.c) CCi<-0 else CCi<-Cf.c
						if (CC >= Cf_c)
							CCi = 0;
						else
							CCi = Cf_c;
						// 		 #cat("\n WARNING: maximum in.degree set to",max.con,"in iteration",it+1,"connecting the remaining",length(g_h)+length(g_h.new),"hubs, because lower in.degree is not compatible with module structure and other parameters setting\n")
						//cat("\n WARNING: maximum in_degree set to",max_con,"in iteration",it+1,"connecting the remaining",length(g_h)+length(g_h_new),"hubs, because lower in_degree is not compatible with module structure and other parameters setting\n");
						// 		}
					}
					// 	     else
					else {

						// 		{if (gamma==0) stop("computation failed! Try with different parameter settings!")
						if (Uguale(gamma, 0.0))
							error("computation failed! Try with different parameter settings!");
						// 		 gamma<-max(gamma-0.2,0)
						gamma = max_s_d(gamma - 0.2, 0.0);
						//                  if (gamma!=0)
						if (DIVERSO(gamma, 0.0)) {
							// 		   {g_Prob<-c(seq(1,N,1)^(-gamma),0)
							g_tmp1_d = seq_d(g_tmp1_d, 1.0, (double) N, 1.0);
							g_tmp2_d = exp_d(g_tmp2_d, g_tmp1_d, -gamma);
							ASSEGNAv_d(g_scalare_d, 1, 0.0);
							g_Prob = vettore2v_d(g_Prob, g_tmp1_d, g_scalare_d);
							//  		    g_Prob<-g_Prob/(sum(g_Prob)) // cancella NA = 0
							dividi1_vs_d(g_Prob, somma_v_d(g_Prob, false));
							//  		    #g_p<-g_Prob[seq(1,max.con,1)]/sum(g_Prob[seq(1,max.con,1)])
							g_tmp1_i = seq_i(g_tmp1_i, 1, max_con, 1);
							g_tmp2_d = copia_v_indx_d(g_tmp2_d, g_Prob, g_tmp1_i);
							// cancNA in somma = 0
							g_p = dividi_vs_d(g_p, g_tmp2_d, somma_v_d(g_tmp2_d, false));
							//  		    g_STout<-c(NA,g_Prob*N)
							g_tmp1_d = moltiplica_vs_d(g_tmp1_d, g_Prob, (double) N);
							ASSEGNAv_d(g_scalare_d, 1, NA_REAL);
							g_STout = vettore2v_d(g_STout, g_scalare_d, g_tmp1_d);
							// 		    if (INdegree=="out") g_STin[2:(max.con+1)]<-N*g_Prob[1:max.con]/sum(g_Prob[1:max.con])
							if (!strncmp(INdegree->str, OUT, 3)) {
								g_tmp1_d = segmento_v_d(g_tmp1_d, g_Prob, 1, max_con);
								g_tmp2_d = moltiplica_vs_d(g_tmp2_d, g_tmp1_d, (double) N);
								// in somma cancNA = 0
								dividi1_vs_d(g_tmp2_d, somma_v_d(g_tmp1_d, false));
								assegna1_v_segmv_d(g_STin, 2, max_con + 1, g_tmp2_d);
								//   else g_STin[2:(max.con+1)]<-NA
							} else
								assegna1_v_segm_d(g_STin, 2, max_con + 1, NA_REAL);
							// 		    Cf.c<-Cf.cl
							Cf_c = Cf_cl;
							// 		    if (CC>=Cf.c) CCi<-0 else CCi<-Cf.c
							if (CC >= Cf_c)
								CCi = 0;
							else
								CCi = Cf_c;
							// 		    #cat("\n WARNING: gamma set to",gamma,"in iteration",it+1,"connecting the remaining",length(g_h)+length(g_h.new),"hubs, because higher gamma is not compatible with module structure and other parameters setting\n")
							//cat("\n WARNING: gamma set to",gamma,"in iteration",it+1,"connecting the remaining",length(g_h)+length(g_h_new),"hubs, because higher gamma is not compatible with module structure and other parameters setting\n");
							// 		   }
						}
						// 		  else
						else {
							// 		   {g_Prob<-c(rep(1/N,N),0)
							CREAv_d(g_tmp1_d, N);
							InitVett_d(g_tmp1_d, (double) 1.0 / N);
							ASSEGNAv_d(g_scalare_d, 1, 0.0);
							g_Prob = vettore2v_d(g_Prob, g_tmp1_d, g_scalare_d);
							// g_p<-g_Prob[seq(1,max.con,1)]/sum(g_Prob[seq(1,max.con,1)])
							g_tmp1_i = seq_i(g_tmp1_i, 1, max_con, 1);
							g_tmp2_d = copia_v_indx_d(g_tmp2_d, g_Prob, g_tmp1_i);
							// cancNA = 0 in somma
							g_p = dividi_vs_d(g_p, g_tmp2_d, somma_v_d(g_tmp2_d, false));
							//  		    g_STout<-c(NA,g_Prob*N)
							g_tmp1_d = moltiplica_vs_d(g_tmp1_d, g_Prob, (double) N);
							ASSEGNAv_d(g_scalare_d, 1, NA_REAL);
							g_STout = vettore2v_d(g_STout, g_scalare_d, g_tmp1_d);

							// 		    if (INdegree=="out") g_STin[2:(max.con+1)]<-N*g_Prob[1:max.con]/sum(g_Prob[1:max.con])
							if (!strncmp(INdegree->str, OUT, 3)) {
								g_tmp1_d = segmento_v_d(g_tmp1_d, g_Prob, 1, max_con);
								g_tmp2_d = moltiplica_vs_d(g_tmp2_d, g_tmp1_d, (double) N);
								// in somma cancNA = 0
								dividi1_vs_d(g_tmp2_d, somma_v_d(g_tmp1_d, false));
								assegna1_v_segmv_d(g_STin, 2, max_con + 1, g_tmp2_d);
								//   else g_STin[2:(max.con+1)]<-NA
							} else
								assegna1_v_segm_d(g_STin, 2, max_con + 1, NA_REAL);
							// 		    Cf.c<-Cf.cl
							Cf_c = Cf_cl;
							// 		    if (CC>=Cf.c) CCi<-0 else CCi<-Cf.c
							if (CC >= Cf_c)
								CCi = 0;
							else
								CCi = Cf_c;
							// 		    #cat("\n WARNING: power law distribution replaced by flat distribution in iteration",it+1,"connecting the remaining",length(g_h)+length(g_h.new),"hubs, because power-law distribution is not compatible with module structure and other parameters setting\n")
							//cat("\n WARNING: power law distribution replaced by flat distribution in iteration",it+1,"connecting the remaining",length(g_h)+length(g_h_new),"hubs, because power-law distribution is not compatible with module structure and other parameters setting\n");
							// 		   }
						}
						// 		 }
					}
					// 	     }
				}
				//         }
			}
			//      }#end while (controllo==1)
		} //end while (controllo==1);
		//     Ng<-dim(M)[1]
		Ng = LENGTHm1_i(M1);
		//     g_aus<-assign.nodes(M,Mdiscr,g_h,hubs=hubs,Sc=Sc,g_Sin,max.con)
		g_aus0 = assign_nodes2(g_aus0, M1, g_Mdiscr, g_h, hubs, Sc_m, g_Sin, max_con);
		// Mdiscr e g_h vengono modificate, ma sono passate per riferimento, per cui non serve aggiornarle
		//     g_Mdiscr<-g_aus[[1]]
		// Mdiscr = copia_m_d(Mdiscr, g_aus0->dati[0].md);
		//     g_h<-g_aus[[2]]
		// g_h = copia_v_i(g_h, g_aus0->dati[1].vi, 1, LENGTHv_i(g_h));
		//     g_h.new<-c(g_h.new,g_aus[[3]])
		g_h_new = accoda1_vv_i(g_h_new, g_aus0);
		//     g_Sout<-apply(Mdiscr,2,sum)
		g_Sout = somma_colonne_i(g_Sout, g_Mdiscr);
		//     m<-max(g_Sout)+1
		m = max_v_i(g_Sout) + 1;
		//     Freq.out[1:m]<-hist(g_Sout,breaks=seq(0,m,1),right=FALSE,plot=FALSE)$counts
		g_tmp2_i = seq_i(g_tmp2_i, 0, m, 1);
		g_tmp1_i = hist1(g_tmp1_i, g_Sout, g_tmp2_i, 0, 1, 0);
		g_tmp1_d = promuovi_i(g_tmp1_d, g_tmp1_i);
		assegna1_v_segmv_d(g_Freq_out, 1, m, g_tmp1_d);
		//     g_Sin<-apply(Mdiscr,1,sum)
		g_Sin = somma_righe_i(g_Sin, g_Mdiscr);
		//     m<-max(g_Sin)+1
		m = max_v_i(g_Sin) + 1;
		//     Freq.in[1:m]<-hist(g_Sin,breaks=seq(0,m,1),right=FALSE,plot=FALSE)$counts
		g_tmp2_i = seq_i(g_tmp2_i, 0, m, 1);
		g_tmp1_i = hist1(g_tmp1_i, g_Sin, g_tmp2_i, 0, 1, 0);
		g_tmp1_d = promuovi_i(g_tmp1_d, g_tmp1_i);
		assegna1_v_segmv_d(g_Freq_in, 1, m, g_tmp1_d);
		//     Nrim<-Nrim-Ng
		Nrim -= Ng;
		//     if (Nrim<=1)
		if (Nrim <= 1) {
			//      {it<-it+1
			it++;
			//       if ((it==1)&(Nrim==1)) g_h<-c(g_h,g_h.new)
			if (it == 1 && Nrim == 1)
				g_h = accoda1_vv_i(g_h, g_h_new);
			//        else g_h<-g_h.new
			else
				g_h = copia_v_i(g_h, g_h_new, 1, LENGTHv_i(g_h_new)); // copia, come nell'instruzione di R
			//       g_h.new<-vector()
			CREAv_i(g_h_new, 0);
			//       if (length(g_h)>0)
			if (LENGTHv_i(g_h) > 0) {
				// g_h<-g_h[which(g_Sin[g_h]!=max.con)]
				g_tmp1_i = copia_v_indx_i(g_tmp1_i, g_Sin, g_h);
				g_tmp2_i = which_v_indxne_i(g_tmp2_i, g_tmp1_i, max_con);
				// gli indici sono sempre in ordine, quindi non c'e` problema
				g_h = copia_v_indx_i(g_h, g_h, g_tmp2_i);
				//       Pm1<-Pm2<-Pm3<-1/3
				Pm1 = Pm2 = Pm3 = (double) 1.0 / 3.0;
			}
			//      }
		}
		//    }# END WHILE  (length(g_h)>1)
	} // END WHILE  (length(g_h)>1);
	// ###check that every gene has at least 1 regulator
	//##check that every gene has at least 1 regulator;
	// g_Sr<-apply(Mdiscr,1,sum)
	g_Sr = somma_righe_i(g_Sr, g_Mdiscr);
	// g_ind<-which(g_Sr==0)
	g_ind = which_v_indxeq_i(g_ind, g_Sr, 0);
	// L<-length(g_ind)
	L = LENGTHv_i(g_ind);
	// if (L>0)
	if (L > 0) {
		//  {for (i in (1:L))
		for (i = 1; i <= L; i++) {
			//    {ri<-g_ind[i]
			ri = ACCEDIv_i(g_ind, i);
			//     num<-1
			num = 1;
			//     g_Sout<-apply(Mdiscr,2,sum)
			g_Sout = somma_colonne_i(g_Sout, g_Mdiscr);
			//     Sc<-Score(S=g_Sout,ST=g_STout,Freq=Freq.out,n=1,g_toll)
			g_Sc_v = score1(g_Sc_v, g_Sout, g_STout, g_Freq_out, 1, g_toll);
			//     g_ind.Sc<-which(Sc==(-Inf))
			g_ind_Sc = which_v_indxeq_d(g_ind_Sc, g_Sc_v, R_NegInf);
			//     Sc[g_ind.Sc]<-0
			assegna1_v_indx_d(g_Sc_v, g_ind_Sc, 0.0);
			//     if (length(which(Sc<=0))>0) Sc<-Sc-min(Sc)+1/(N^2)
			g_tmp1_i = which_v_indxle_d(g_tmp1_i, g_Sc_v, 0.0);
			if (LENGTHv_i(g_tmp1_i) > 0)
				somma1_vs_d(g_Sc_v,  -min_v_d(g_Sc_v) + (double) 1.0 / (N * N));
			//     Sc[g_ind.Sc]<-0
			assegna1_v_indx_d(g_Sc_v, g_ind_Sc, 0.0);
			//     if (sum(Sc)==0) # cancNA = 0
			if (Uguale(somma_v_d(g_Sc_v, false), 0.0)) {
				// 	{#cat("\n WARNING: tolerance parameters were relaxed to allow each node to have at least 1 regulator")
				//cat("\n WARNING: tolerance parameters were relaxed to allow each node to have at least 1 regulator");
				// 	 Sc<-Score(S=g_Sout,ST=g_STout,Freq=Freq.out,n=1,g_toll=rep(Inf,N))
				CREAv_d(g_tmp1_d, N);
				InitVett_d(g_tmp1_d, R_PosInf);
				g_Sc_v = score1(g_Sc_v, g_Sout, g_STout, g_Freq_out, 1, g_tmp1_d);
				// 	 g_ind.Sc<-which(Sc==(-Inf))
				g_ind_Sc = which_v_indxeq_d(g_ind_Sc, g_Sc_v, R_NegInf);
				//   	 Sc[g_ind.Sc]<-0
				assegna1_v_indx_d(g_Sc_v, g_ind_Sc, 0.0);
				//   	 if (length(which(Sc<=0))>0) Sc<-Sc-min(Sc)+1/(N^2)
				g_tmp1_i = which_v_indxle_d(g_tmp1_i, g_Sc_v, 0.0);
				if (LENGTHv_i(g_tmp1_i) > 0)
					somma1_vs_d(g_Sc_v, min_v_d(g_Sc_v) + (double) 1.0 / (N * N));
				//   	 Sc[g_ind.Sc]<-0
				assegna1_v_indx_d(g_Sc_v, g_ind_Sc, 0.0);
				// 	}
			}
			//     g_p.out<-Sc/sum(Sc) # cancNA = 0
			g_p_out = dividi_vs_d(g_p_out, g_Sc_v, somma_v_d(g_Sc_v, false));
			//     g_ind.s<-sampleB(seq(1,N,1),num,prob=g_p.out)
			g_tmp1_i = seq_i(g_tmp1_i, 1, N, 1);
			// senza ripetizioni
			g_ind_s = sampleB_p(g_ind_s, g_tmp1_i, num, 0, g_p_out);
			//     Mdiscr[ri,g_ind.s]<-1
			assegna1_ms_rigaindx_i(g_Mdiscr, ri, g_ind_s, 1);
			//     g_a1<-g_Sout[g_ind.s]
			g_a1 = copia_v_indx_i(g_a1, g_Sout, g_ind_s);
			//     g_a2<-g_a1+1
			g_a2 = somma_vs_i(g_a2, g_a1, 1);
			//     Freq.out[g_a1+1]<-Freq.out[g_a1+1]-1
			incr1_v_i(g_a1, 1);
			incr1_v_indx_d(g_Freq_out, g_a1, -1.0);
			//     Freq.out[g_a2+1]<-Freq.out[g_a2+1]+1
			incr1_v_i(g_a2, 1);
			incr1_v_indx_d(g_Freq_out, g_a2, 1.0);
			//    }
		}
		//  }
	}
	// #######assign weights to the connectivity matrix##########
	//######assign weights to the connectivity matrix##########;
	// g_ind<-which(Mdiscr==1,arr.g_ind=TRUE)
	g_ind = which_m_indxeq_i(g_ind, g_Mdiscr, 1);
	// L<-dim(g_ind)[1]
	L = LENGTHv_i(g_ind);
	// g_aus<-abs(rnorm_s(L,weight.mean, weight.sd))
	CREAv_d(g_tmp1_d, L);
	g_tmp1_d = rnorm_s(g_tmp1_d, L, weight_mean, weight_sd, "cm");
	g_tmp2_d = abs_v_d(g_tmp2_d, g_tmp1_d);
	// M<-matrix(0,ncol=N,nrow=N)
	// questa M non e` la stessa di prima, infatti sono di tipo diverso!
	CREAm_d(g_M, N, N);
	InitMatr_d(g_M, 0.0);
	// M[g_ind]<-aus
	assegna1_mv_indx_d(g_M, g_ind, g_tmp2_d);
	// return(list(M,Mdiscr))
	tipi[0] = MATRd;
	tipi[1] = MATRi;
	CreaLISTA(ris, tipi, 2);
	ris->dati[0].md = g_M;
	ris->dati[1].mi = g_Mdiscr;
	//~ CANCELLAv_i(g_scalare_i);
	//~ CANCELLAv_i(g_tmp1_i);
	//~ CANCELLAv_i(g_tmp2_i);
	//~ CANCELLAv_d(g_tmp1_d);
	//~ CANCELLAv_d(g_tmp2_d);
	//~ CANCELLAv_d(g_Prob);
	//~ CANCELLAv_d(g_Freq_out);
	//~ CANCELLAv_d(g_Freq_in);
	//~ CANCELLAv_d(g_STout);
	//~ CANCELLAv_d(g_STin);
	//~ CANCELLAv_d(g_toll);
	//~ CANCELLAv_d(g_p_out);
	//~ CANCELLAv_d(g_scalare_d);
	//~ CANCELLAv_i(g_Sr);
	//~ CANCELLAv_d(g_p);
	//~ CANCELLAv_i(g_h);
	//~ CANCELLAv_d(g_prob_mod);
	//~ CANCELLAm_d(g_aus);
	//~ CANCELLAv_d(g_Sc_v);
	//~ CANCELLAv_d(g_Cg);
	//~ CANCELLAv_i(g_Sin);
	//~ CANCELLAv_i(g_Sout);
	//~ CANCELLAv_i(g_a1);
	//~ CANCELLAv_i(g_a2);
	//~ CANCELLAv_i(g_ind);
	//~ CANCELLAv_i(g_ind_Sc);
	//~ CANCELLAv_i(g_h_new);
	//~ CANCELLAv_i(g_mod_type);
	//~ CANCELLAv_i(g_ind_s);
	//~ // CANCELLAv_i(hubs); hubs e` uno dei ris*.indices, che cancellero` tra un po'
	//~ CANCELLAv_i(g_aus0);
	//~ CANCELLAm_i(g_ris1.conn_matr);
	//~ CANCELLAm_i(g_ris2.conn_matr);
	//~ CANCELLAm_i(g_ris3.conn_matr);
	//~ CANCELLAm_d(g_ris1.score_matr);
	//~ CANCELLAm_d(g_ris2.score_matr);
	//~ CANCELLAm_d(g_ris3.score_matr);
	//~ CANCELLAv_i(g_ris1.indices);
	//~ CANCELLAv_i(g_ris2.indices);
	//~ CANCELLAv_i(g_ris3.indices);

	StrBilanciam();

#ifdef DET
	fprintf(fp_det, "cm output:\n");
	fprintf(fp_det, "\tlist(M,Mdiscr) = ");
	_StampaRawMatr_d(g_M);
	_StampaRawMatr_i(g_Mdiscr);
#endif

	return ris;
	// }
}

SEXP connectivity_modular(SEXP N, SEXP max_con, SEXP gamma, SEXP INdegree, SEXP Cf_cl, SEXP num_subnet, SEXP r_tol, SEXP a_tol, SEXP weight_mean, SEXP weight_sd)
{
	int nProtected = 0;
	int N1, max_con1;
	double gamma1, Cf_cl1, r_tol1, a_tol1, weight_mean1, weight_sd1;
	GString *INdegree1 = NULL;
	VETTOREi *num_subnet1;
	LISTA *l = NULL;
	SEXP ris;

	_InitDbg(false, false, false);

	_Intestazione("\n*** connectivity_modular ***\n");

	N1 = INTEGER_VALUE(N);
	max_con1 = INTEGER_VALUE(max_con);
	gamma1 = NUMERIC_VALUE(gamma);
	INdegree1 = inSTRINGA(INdegree, &nProtected, "INdegree");
	Cf_cl1 = NUMERIC_VALUE(Cf_cl);
	num_subnet1 = inVETTORE_i(num_subnet, &nProtected);
	r_tol1 = NUMERIC_VALUE(r_tol);
	a_tol1 = NUMERIC_VALUE(a_tol);
	weight_mean1 = NUMERIC_VALUE(weight_mean);
	weight_sd1 = NUMERIC_VALUE(weight_sd);

	l = connectivity_modular1(l, N1, max_con1, gamma1, INdegree1, Cf_cl1, num_subnet1, r_tol1, a_tol1, weight_mean1, weight_sd1);
	ris = daLISTA(l, &nProtected);

	CANCELLAstr(INdegree1);
	CANCELLAv_i(num_subnet1);

	StrBilanciam();
	ControllaCanc();

	UNPROTECT(nProtected);

	return ris;
}
