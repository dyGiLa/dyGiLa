# >>>>>>>>>>>  read me,plz <<<<<<<<<<<
# This is a short makefile for compiling phases_diagram.cpp.
# Because there are two *makefiles* in same directory, the -f or --file=
# flag must be used i.e.,
#      make --file=makefile [target]
# the *run* target will excute the binary ./pd.app and make plosts.
#
# 06.08.2022, elokuuta
# Quang. Zhang, timohyva@github

CPPc = g++
EXEC = pd.app
STD = -std=c++17
OBJS = phases_diagram.o matep.o

${EXEC}: ${OBJS}
	${CPPc} -g -o ${EXEC} ${OBJS}

matep.o: matep.cpp matep.hpp
	${CPPc} ${STD} -g -c matep.cpp matep.hpp

phases_diagram.o: phases_diagram.cpp matep.hpp
	${CPPc} ${STD} -g -c phases_diagram.cpp matep.hpp

.PHONY: clean
clean:
	rm -f ${OBJS}

.PHONY: run
run:
	./pd.app && python3 ./pd/sfdata.py
