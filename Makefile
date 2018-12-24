# makefile for version2_thread
# from fileDisposal.o , the command is automatically derived by make

Object =  coordinaryTransform.o    fileDisposal.o    initialization.o \
          main.o    matrixFunction.o    MOPSOAidFunction.o    \
          MOPSOFunction.o    parameterSettings.o \
	setTime.o     preDispose.o  multiplyBetterInput.o\
       	disposedCodes.o	 seiveRep.o  inputParticles.o  createRep.o\
	checkSimilarity.o	Debug.o	freeUpSpace.o	applyVariables.o \
	disposeTM_align.o DH_rotation.o inputFilter.o

AIR : $(Object)
	g++ -o AIR *o -lpthread

coordinaryTransform.o : MOPSO.h
	g++ -c coordinaryTransform.cpp

DH_rotation.o : MOPSO.h DH_rotation.h

fileDisposal.o : MOPSO.h

initialization.o : MOPSO.h

main.o : MOPSO.h DH_rotation.h TM_score_Combination_Order.h inputFilter.h

matrixFunction.o : MOPSO.h

MOPSOAidFunction.o : MOPSO.h

MOPSOFunction.o : MOPSO.h

parameterSettings.o : MOPSO.h

setTime.o : MOPSO.h

preDispose.o : MOPSO.h

multiplyBetterInput.o : MOPSO.h inputFilter.h

disposedCodes.o : MOPSO.h

seiveRep.o : MOPSO.h

inputParticles.o : MOPSO.h

updateRep.o : MOPSO.h

checkSimilarity.o : MOPSO.h

Debug.o : MOPSO.h

freeUpSpace.o : MOPSO.h

applyVariables.o : MOPSO.h

disposeTM_align.o : MOPSO.h TM_score_Combination_Order.h

inputFilter.o : inputFilter.h MOPSO.h

.PHONY : clear
clear :
	-rm *o AIR
.PHONY : delete
delete :
	-rm -r /home/genAIRing_new/mopso_v1/data/answer/newAnswer

