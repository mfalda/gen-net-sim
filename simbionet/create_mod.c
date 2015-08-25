#include "create_mod.h"

GList *createMOD1(int m, bool auto1, int *len)
{
	int i, j , b, L, k, aus;
	double CC;
	bool simm, fb;
	VETTOREi *codes = NULL, *tmp1_i = NULL, *tmp2_i = NULL;
	VETTOREd * ris_cc = NULL;
	MATRICEi *M = NULL, *tmpm1 = NULL, *mod = NULL;
	GString *etich = NULL, *path = NULL;
	LISTA *H = NULL;
	Mod *elem_ris;
	GList *lista_ris = NULL;

	_Intestazione("\n*** createMOD1 ***\n");
#ifdef DET
	fprintf(fp_det, "input:\n");
	fprintf(fp_det, "\tm = %d\n", m);
	if (auto1)
		fprintf(fp_det, "\tauto = TRUE\n");
	else
		fprintf(fp_det, "\tauto = FALSE\n");
#endif

//  MOD<-list()
//  b<-1
	b = 1;
//  path<-"modules/"
	CREAstr(path, "modules/");
	CREAstr(etich, "");
//  for (i in (2:m))
	for (i = 2; i <= m; i++) {
//   if (auto==FALSE) etich<-paste(path,"M",i,".txt",sep="")
		if (!auto1)
			g_string_printf(etich, "%s/library/SimBioNeT/%sM%d.txt", g_getenv("R_HOME"), path->str, i);
//    else  etich<-paste("autoreg_",path,"M",i,".txt",sep="")
		else
			g_string_printf(etich, "%s/library/SimBioNeT/autoreg_%sM%d.txt", g_getenv("R_HOME"), path->str, i);
//    M<-as.matrix(read.table(etich))
		M = read_m_i(M, etich->str);
//    codes<-as.vector(M[1,])
		codes = riga_i(codes, M, 1);
//    M<-M[-1,]
		elimina1_riga_i(M, 1);
//    L<-dim(M)[2]
		L = LENGTHm2_i(M);
//    for (j in (1:L))
		for (j = 1; j <= L; j++) {
//        mod<-matrix(M[,j],i,i)
			tmp1_i = colonna_i(tmp1_i, M, j);
			CREAm_i(mod, i, i);
			for (k = 1; k <= LENGTHv_i(tmp1_i); k++)
				ASSEGNAmv_i(mod, k, ACCEDIv_i(tmp1_i, k));
//         aus<-length(which((mod==t(mod))==0)) cioe` which(m != t(m))
			tmp2_i = which_m_indxnsimm_i(tmp2_i, mod);
			aus = LENGTHv_i(tmp2_i);
//         if (aus==0) simm<-TRUE
			if (aus == 0)
				simm = true;
//         else simm<-FALSE
			else
				simm = false;
//         H<-hubs(mod)
			H = hubs(H, mod, &fb);
			ris_cc = cluster_coeff2(ris_cc, mod, &CC);

//         MOD[[b]]<-list("code"=codes[j],"net"=mod,"hubs"=H[[1]],"CC"=cluster.coeff(mod)[[1]],"autoreg"=auto,"feedback"=H[[2]],"hubsio"=H[[3]],"SIMM"=simm,"dim.m"=i)
			elem_ris = mia_alloc(1, Mod);
			if (1 > 0 && elem_ris == NULL) {
				Rprintf("Not enough memory (create_mod # %d, elem_ris)", __LINE__ - 2);
				error("");
			}
			elem_ris->codice = ACCEDIv_i(codes, j);
			elem_ris->rete = NULL;
			elem_ris->rete = copia_m_i(elem_ris->rete, mod);
			CtrlLlst(H, 1);
			elem_ris->hubs = NULL;
			elem_ris->hubs = copia_v_i(elem_ris->hubs, ACCEDIlst(H, 1, vi), 1, LENGTHv_i(ACCEDIlst(H, 1, vi)));
			elem_ris->CC = CC;
			elem_ris->autoreg = auto1;
			elem_ris->feedback = fb;
			CtrlLlst(H, 2);
			elem_ris->hubsio = NULL;
			elem_ris->hubsio = copia_v_i(elem_ris->hubsio, ACCEDIlst(H, 2, vi), 1, LENGTHv_i(ACCEDIlst(H, 2, vi)));
			elem_ris->simm = simm;
			elem_ris->dim_m = i;
			lista_ris = g_list_append(lista_ris, elem_ris);
//         b<-b+1
			b++;
//       }
		}
//   }
	}

	*len = b;


	CANCELLAv_i(codes);
	CANCELLAv_i(tmp1_i);
	CANCELLAv_i(tmp2_i);
	CANCELLAv_d(ris_cc);
	CANCELLAm_i(M);
	CANCELLAm_i(mod);
	CANCELLAm_i(tmpm1);
	CancellaLISTA(H, false);
	CANCELLAstr(path);
	CANCELLAstr(etich);

	StrBilanciam();

	_Intestazione("\n*** Esco da createMOD1 ***\n");

// return(MOD)
	return lista_ris;
// }
}

// createMOD<-function(m=4,auto=FALSE)
SEXP createMOD(SEXP m, SEXP auto1)
{
	int i, nProtected = 0, m1, len;
	bool auto11;
	GList *ris, *li;
	Mod *elem_ris;
	SEXP rete1, hubs1, hubsio1;
	SEXP ret_lista, ret_lista1, nomi;

	_InitDbg(false, false, false);

	_Intestazione("\n*** createMOD ***\n");

	m1 = INTEGER_VALUE(m);
	auto11 = LOGICAL_VALUE(auto1);

	InitGlobali();

	ris = createMOD1(m1, auto11, &len);

	// alloca la lista per il risultato
	PROTECT(ret_lista = allocVector(VECSXP, len - 1));
	++nProtected;

	for (i = 0, li = ris; li != NULL; li = li->next, i++) {
		elem_ris = (Mod *) li->data;
		rete1 = daMATRICE_i(elem_ris->rete, &nProtected);
		hubs1 = daVETTORE_i(elem_ris->hubs, &nProtected);
		hubsio1 = daVETTORE_i(elem_ris->hubsio,
&nProtected);
		PROTECT(ret_lista1 = allocVector(VECSXP, 9));
		++nProtected;
		SET_VECTOR_ELT(ret_lista1, 0, Rf_ScalarInteger(elem_ris->codice));
		SET_VECTOR_ELT(ret_lista1, 1, rete1);
		SET_VECTOR_ELT(ret_lista1, 2, hubs1);
		SET_VECTOR_ELT(ret_lista1, 3, Rf_ScalarReal(elem_ris->CC));
		SET_VECTOR_ELT(ret_lista1, 4, Rf_ScalarLogical(elem_ris->autoreg));
		SET_VECTOR_ELT(ret_lista1, 5, Rf_ScalarLogical(elem_ris->feedback));
		SET_VECTOR_ELT(ret_lista1, 6, hubsio1);
		SET_VECTOR_ELT(ret_lista1, 7, Rf_ScalarLogical(elem_ris->simm));
		SET_VECTOR_ELT(ret_lista1, 8, Rf_ScalarInteger(elem_ris->dim_m));
		PROTECT(nomi = allocVector(STRSXP, 9));
		++nProtected;
		SET_STRING_ELT(nomi, 0, mkChar("code"));
		SET_STRING_ELT(nomi, 1, mkChar("net"));
		SET_STRING_ELT(nomi, 2, mkChar("hubs"));
		SET_STRING_ELT(nomi, 3, mkChar("CC"));
		SET_STRING_ELT(nomi, 4, mkChar("autoreg"));
		SET_STRING_ELT(nomi, 5, mkChar("feedback"));
		SET_STRING_ELT(nomi, 6, mkChar("hubsio"));
		SET_STRING_ELT(nomi, 7, mkChar("SIMM"));
		SET_STRING_ELT(nomi, 8, mkChar("dim.m"));
		setAttrib(ret_lista1, R_NamesSymbol, nomi);
		SET_VECTOR_ELT(ret_lista, i, ret_lista1);
	}

	CancGlobali();

	StrBilanciam();
	ControllaCanc();

	g_list_foreach(ris, (GFunc)g_free, NULL);
	g_list_free(ris);

	UNPROTECT(nProtected);

	return ret_lista;
}
