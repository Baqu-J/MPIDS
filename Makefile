TARGET = greedy local_search metaheuristic
CXXFLAGS = -ansi -O3 -march=native -fpermissive -std=c++17 -g3 -pthread -lpthread
OBJS = Random.o Timer.o
CPLOBJS = Random.o Timer.o

SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic
GCC = gcc
CCC = g++
CCOPT = -O3 -march=native -g3 -fPIC -fexceptions -DNDEBUG -DIL_STD -std=c++17 -fpermissive -w

CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR)

all: ${TARGET}
greedy: greedy.cpp $(OBJS)
	${CCC} ${CXXFLAGS} -o $@ $^

local_search: local_search.cpp $(OBJS)
	${CCC} ${CXXFLAGS} -o $@ $^

metaheuristic: metaheuristic.cpp $(OBJS)
	${CCC} ${CXXFLAGS} -o $@ $^

clean:
	@rm -f *~ *.o ${TARGET} core


