.PHONY: all clean debug profile
ROOT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
export BIN_DIR = ${ROOT_DIR}/../bin

export MINIMAP2_DIR = ${ROOT_DIR}/../lib/minimap2

CXXFLAGS += -Wall -Wextra -pthread -std=c++11 -g
LDFLAGS += -pthread -std=c++11 -rdynamic

ASSEMBLE_BIN := ${BIN_DIR}/flye-assemble
REPEAT_BIN := ${BIN_DIR}/flye-repeat
CONTIGGER_BIN := ${BIN_DIR}/flye-contigger
POLISH_BIN := ${BIN_DIR}/flye-polish

profile: CXXFLAGS += -pg
profile: LDFLAGS +=  -lz -L${MINIMAP2_DIR} -lminimap2
profile: release

release: CXXFLAGS += -O3 -DNDEBUG
release: flye-polish flye-assemble

SANITIZE_FLAGS += -D_GLIBCXX_SANITIZE_VECTOR -U_FORTIFY_SOURCE -fsanitize=address -fno-omit-frame-pointer -fsanitize=undefined
#SANITIZE_FLAGS += -fsanitize=thread -fsanitize=leak -fsanitize=undefined
debug: CXXFLAGS += -Og ${SANITIZE_FLAGS}
#debug: CXXFLAGS += -D_GLIBCXX_DEBUG
debug: LDFLAGS += ${SANITIZE_FLAGS}
debug: flye-repeat flye-polish flye-assemble flye-contigger

#sequence module
sequence_obj := ${patsubst %.cpp,%.o,${wildcard sequence/*.cpp}}

sequence/%.o: sequence/%.cpp sequence/*.h common/*.h
	${CXX} -c ${CXXFLAGS} $< -o $@

#flye-assemble module
assemble_obj := ${patsubst %.cpp,%.o,${wildcard assemble/*.cpp}}
flye-assemble: ${assemble_obj} ${sequence_obj} bin/assemble.o
	${CXX} ${assemble_obj} ${sequence_obj} bin/assemble.o -o ${ASSEMBLE_BIN} ${LDFLAGS}

assemble/%.o: assemble/%.cpp assemble/*.h sequence/*.h common/*.h
	${CXX} -c ${CXXFLAGS} $< -o $@
bin/assemble.o: assemble/*.h sequence/*.h common/*.h bin/assemble.cpp
	${CXX} -c ${CXXFLAGS} bin/assemble.cpp -o $@

#flye-repeat module
repeat_obj := ${patsubst %.cpp,%.o,${wildcard repeat_graph/*.cpp}}
flye-repeat: ${repeat_obj} ${sequence_obj} bin/repeat.o
	${CXX} ${repeat_obj} ${sequence_obj} bin/repeat.o -o ${REPEAT_BIN} ${LDFLAGS}

repeat_graph/%.o: repeat_graph/%.cpp repeat_graph/*.h sequence/*.h common/*.h
	${CXX} -c ${CXXFLAGS} $< -o $@
bin/repeat.o: repeat_graph/*.h sequence/*.h common/*.h bin/repeat.cpp
	${CXX} -c ${CXXFLAGS} bin/repeat.cpp -o $@

#flye-contigger module
contigger_obj := ${patsubst %.cpp,%.o,${wildcard contigger/*.cpp}}
flye-contigger: ${contigger_obj} ${sequence_obj} ${repeat_obj} bin/contigger.o
	${CXX} ${contigger_obj} ${repeat_obj} ${sequence_obj} bin/contigger.o -o ${CONTIGGER_BIN} ${LDFLAGS}

contigger/%.o: contigger/%.cpp repeat_graph/*.h sequence/*.h common/*.h
	${CXX} -c ${CXXFLAGS} $< -o $@
bin/contigger.o: repeat_graph/*.h contigger/*.h sequence/*.h common/*.h bin/contigger.cpp
	${CXX} -c ${CXXFLAGS} bin/contigger.cpp -o $@


#flye-polish module
polish_obj := ${patsubst %.cpp,%.o,${wildcard polishing/*.cpp}}
flye-polish: ${polish_obj} bin/polisher.o
	${CXX} ${polish_obj} bin/polisher.o -o ${POLISH_BIN} ${LDFLAGS}

polishing/%.o: polishing/%.cpp bin/polisher.cpp polishing/*.h common/*h
	${CXX} -c ${CXXFLAGS} $< -o $@



clean:
	-rm ${repeat_obj}
	-rm ${sequence_obj}
	-rm ${assemble_obj}
	-rm ${polish_obj}
	-rm ${contigger_obj}
	-rm ${patsubst %.cpp,%.o,${wildcard bin/*.cpp}}
	-rm ${REPEAT_BIN}
	-rm ${ASSEMBLE_BIN}
	-rm ${POLISH_BIN}
	-rm ${CONTIGGER_BIN}
