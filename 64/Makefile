# Makefile for connecting many subroutines
#
#
#                TIANTIAN XU
#      DEPARTMENT OF COMPUTATIONAL SCIENCE AND ENGINEERING
#                YONSEI UNIVERSITY
#                    2018/01/02

EXE = a.out
RUNPROFILE = gprof

# The following information need to be modified for different directory
pwd =  $(shell pwd)/

main = $(pwd)IBM_Heat_A1.o
OBJS = $(main)

#compile
FCFLAG =  -O3 -W0 -r8 -extend-source -traceback

 
# library information
FFTWLIB = -L/shared/home/tian0917/fftw3/fftw_installed/lib
FFTWINC = -I/shared/home/tian0917/fftw3/fftw_installed/include
FFTWFLAGS  = -lfftw3 -lm
# ifort = /opt/intel/oneapi/compiler/2021.4.0/linux/bin/intel64/ifort

.SUFFIXES: .f90 .o
.f90.o:
	ifort -c -pg $(FFTWLIB) $(FFTWINC) ${FCFLAG} $*.f90 $(FFTWFLAGS)

all:${OBJS}
	ifort -pg ${OBJS} $(FFTWLIB) $(FFTWINC) $(FFTWFLAGS) -o ${EXE} 
#  ifort -pg ${OBJS} -o ${EXE}


#submit the job
exe:
	sbatch serial.sh

rm:
	clear
	rm *.o
	rm *.err
	rm *.out
