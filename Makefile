#CC	= icc
CC	= gcc
MPICC 	= mpicc
CFLAG	= -O3

netconn:
	${CC} ${CFLAG} -o netconn randnum.c netconn.c

stngp:
	${MPICC} ${CFLAG} -o stngpMPI ion_ch.c neuronSTN.c neuronGP.c synAMPA.c synNMDA.c synGABAa.c synBg.c randnum.c stngpMPI.c
