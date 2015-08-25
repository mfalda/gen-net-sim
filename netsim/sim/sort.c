#include "sort.h"


typedef struct {
	int val;
	int pos;
} Elem;

void mergesort_array(Elem a[], int size, Elem temp[])
{
	int i, j, i1, i2, tempi;
	Elem v;
	if (size < MIN_MERGESORT_LIST_SIZE) {
		/* Use insertion sort */
		for (i = 0; i < size; i++) {
			v = a[i];
			for (j = i - 1; j >= 0; j--) {
				if (a[j].val <= v.val) break;
				a[j + 1] = a[j];
			}
			a[j + 1] = v;
		}
		return;
	}

	mergesort_array(a, size / 2, temp);
	mergesort_array(a + size / 2, size - size / 2, temp);
	i1 = 0;
	i2 = size / 2;
	tempi = 0;
	while (i1 < size / 2 && i2 < size) {
		if (a[i1].val < a[i2].val) {
			temp[tempi] = a[i1];
			i1++;
		}
		else {
			temp[tempi] = a[i2];
			i2++;
		}
		tempi++;
	}

	while (i1 < size / 2) {
		temp[tempi] = a[i1];
		i1++;
		tempi++;
	}
	while (i2 < size) {
		temp[tempi] = a[i2];
		i2++;
		tempi++;
	}

	memcpy(a, temp, size * sizeof(Elem));
}


VETTOREi *ordine(VETTOREi *ris, VETTOREi *v)
{
	int size = LENGTHv_i(v);
	Elem *a    = mia_alloc(size, Elem);
	Elem *temp = mia_alloc(size, Elem);
	int i;

	CREAv_i(ris, size);
	for (i = 0; i < size; i++) {
		a[i].val = ACCEDIv_i(v, i + 1);
		a[i].pos = i + 1;
	}
	mergesort_array(a, size, temp);
#ifdef FDEBUG
	for (i = 1; i < size; i++) {
		if (!(a[i-1].val <= a[i].val))
			error("Error in MergeSort");
	}
#endif
	for (i = 0; i < size; i++) {
		ASSEGNAv_i(v, i + 1, a[i].val);
		ASSEGNAv_i(indici, i + 1, a[i].pos);
	}
	libera(a);
	libera(temp);

	return ris;
}