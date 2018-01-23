VPATH = ./src
INCFLAG := -I./include 

objtest=MyFunc.o Field.o Distribution.o LatticeMap.o test.o Mapmain.o  RFQMap.o QuadrupoleMap.o ElementMap.o BendingMap.o BeamRotationMap.o DriftMap.o ThinlensMap.o SolenoidMap.o EdgeangleMap.o GapMap.o BeamMap.o ParticleMap.o 

ifdef defines
DEF     = $(foreach def, $(defines) $(ADDONS), -D$(def))
override CPPFLAGS :=  $(DEF) $(CPPFLAGS)
endif

export
test:${objtest}
	mpicxx -o run ${objtest} #test
method: abprepose3d.o abspace3d.o
	${ARCHIVE} ${ARCHIVEFLAGS} ${LIBCORE} ${objtest}
	${RANLIB} ${LIBCORE}

ifeq ($(LTYPE),g)	# turn on debug flag
CXXFLAGS := -g 
INCFLAG := -I./include -I../../basic/include -I../../core/beam/include -I../solver/fft/include -I../../libs/include -I../../core/device/include
else			# Maximal optimization flags
ifeq ($(LTYPE),opt)
CXXFLAGS := -O3
INCFLAG := -I./include -I../../basic/include -I../../core/beam/include -I../solver/fft/include -I../../libs/include -I../../core/device/include
endif
endif

OBJS    = $(foreach module, $(bases) $(SPECIAL), $(module).o)
SRCE    = $(foreach module, $(bases) $(SPECIAL), $(module).C)
	
clean:
	rm -f *.o
	rm -f *.a
	rm -f run 
	rm -f *.dat
	rm -f ./particledata/*.dat
	rm -f ./temp/*


.C.o:
	$(CXX) $(CXXFLAGS) -w $(INCFLAG) -c $<
.c.o:
	$(CXX) $(CXXFLAGS) -w $(INCFLAG) -c $<
