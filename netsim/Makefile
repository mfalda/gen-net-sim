CC=gcc
AR= ar

ifndef debug
	ifndef moduli
		CFLAGS= -DNDEBUG -c -std=c99 -IG:/C/glib2.18/include/glib-2.0 -IG:/R-2.8.1/include -fshort-enums -mno-cygwin -march=i686 -O2
	else
		CFLAGS= -g -c -std=c99 -DMDEBUG -DFDEBUG -DDET -IG:/C/glib2.18/include/glib-2.0 -IG:/R-2.8.1/include -mno-cygwin -Wconversion -Wfloat-equal -fshort-enums
	endif
else
	CFLAGS= -g -c -std=c99 -DMDEBUG -DFDEBUG -IG:/C/glib2.18/include/glib-2.0 -IG:/R-2.8.1/include -mno-cygwin -Wconversion -Wfloat-equal -fshort-enums
endif

LFLAGS= -lglib -lm -lR

LIBFLAGS= rcs

r2c: r_aux.o r_aux_i.o r_aux_d.o
	$(CC) $(CFLAGS) $(LFLAGS) r_aux.c r_aux_i.c r_aux_d.c
	$(AR) $(LIBFLAGS) libr2c.a r_aux.o r_aux_i.o r_aux_d.o

clean:
	rm -f *.o