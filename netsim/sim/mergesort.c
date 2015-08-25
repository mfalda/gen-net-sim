#include <string.h>
#include <time.h>    /* time() */
#include <stdlib.h>  /* rand() */
#include <stdio.h>   /* puts() */

#define MIN_MERGESORT_LIST_SIZE    32

typedef struct {
	int val;
	int pos;
} Elem;

void mergesort_array(Elem a[], int size, Elem temp[], int decr)
{
	int i, j, i1, i2, tempi;
	Elem v;
	if (size < MIN_MERGESORT_LIST_SIZE) {
		/* Use insertion sort */
		for (i = 0; i < size; i++) {
			v = a[i];
			for (j = i - 1; j >= 0; j--) {
				if ((!decr && a[j].val <= v.val) || (decr && a[j].val >= v.val)) break;
				a[j + 1] = a[j];
			}
			a[j + 1] = v;
		}
		return;
	}

	mergesort_array(a, size / 2, temp, decr);
	mergesort_array(a + size / 2, size - size / 2, temp, decr);
	i1 = 0;
	i2 = size / 2;
	tempi = 0;
	while (i1 < size / 2 && i2 < size) {
		if ((!decr && a[i1].val < a[i2].val) || (decr && a[i1].val > a[i2].val)) {
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


int main(int argc, char* argv[])
{
	int size = atoi(argv[1]);
	int decr = atoi(argv[2]);
	Elem *a    = malloc(sizeof(Elem) * size);
	Elem *temp = malloc(sizeof(Elem) * size);
	int i;
	srand(time(NULL));
	for (i = 0; i < size; i++) {
		a[i].val = rand() % size;
		a[i].pos = i + 1;
		printf(" (v:%d, p:%d)", a[i].val, a[i].pos);
	}
	printf("\n");
	mergesort_array(a, size, temp, decr);
	for (i = 1; i < size; i++) {
		if (!((!decr && a[i-1].val <= a[i].val) || (decr && a[i-1].val >= a[i].val))) {
			puts("ERROR");
			return -1;
		}
	}
	puts("SUCCESS:\n\t");
	for (i = 0; i < size; i++) {
		printf(" (p:%d, v:%d)", a[i].pos, a[i].val);
	}
	printf("\n");
	return 0;
}