#include "next_op.h"

int next_op(const VETTOREi *r)
{
	int i, j, k, l, p;

	_Intestazione("\n***next_op***\n");
#ifdef DET
	fprintf(fp_det, "input:\n");
	fprintf(fp_det, "\tr = ");
	_StampaRawVett_i(r);
#endif

	//~ k<-(-1)
	k = -1;
	//~ l<-length(r);
	l = LENGTHv_i(r);
	//~ for (j in 0:(l-1)){
	for (j = 0; j <= l - 1; j++) {
		//~ i<- (l-j)
		i = (l - j);
		//~ p<-r[i]
		p = ACCEDIv_i(r, i);
		//~ if ((p<(-1))&(k==-1)) k<-i;
		if (p < -1 && k == -1)
			k = i;
	}
#ifdef DET
	fprintf(fp_det, "next_op output:\n\tk =  %d\n", k);
#endif

	return k;
}
