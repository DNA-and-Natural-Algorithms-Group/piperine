
all: spuriousC shared

spuriousC: spuriousC.c
	gcc -Wall -g -O3 -o spuriousC spuriousC.c -I${VIENNAHOME}/H -L${VIENNAHOME}/lib -D'VIENNAHOME="${VIENNAHOME}"' -lRNA -lm

shared: spuriousC.c
	gcc -Wall -c -g -O3 -o spuriousC.o spuriousC.c -I${VIENNAHOME}/H -L${VIENNAHOME}/lib -D'VIENNAHOME="${VIENNAHOME}"'
	gcc -Wall -shared -g -O3 -o spuriousC.so spuriousC.o -I${VIENNAHOME}/H -L${VIENNAHOME}/lib -D'VIENNAHOME="${VIENNAHOME}"' -lRNA -lm

spuriousSSM:
	gcc -Wall -O3 spuriousSSM.c -o spuriousSSM -lm

clean:
	rm -f *.pyc
	rm -Rf *.dSYM
	rm -f spuriousC *.so *.o spuriousE spuriousD spuriousSSM
