GCC = g++
CXXFLAGS = -I./ -I/home/users/linacre/Analysis_2012_test/TopAFB/ma5_v1.0.22/madanalysis/SampleAnalyzer/ -lSampleAnalyzer
LIBFLAGS = -lGpad -lHist -lGraf -lGraf3d -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lRIO -lNet -lThread -lCore -lCint -pthread -lm -ldl -rdynamic -lSampleAnalyzer -lz

SRCS = $(wildcard */*.cpp)
OBJS = $(SRCS:.cpp=.o)
PROGRAM = SampleAnalyzer

all:	 compile link

compile:	$(OBJS)

link:	$(OBJS)
		$(GCC) $(CXXFLAGS) $(OBJS) $(LIBFLAGS) -o $(PROGRAM)
 
clean:;	@rm -f $(OBJS) $(PROGRAM) compilation.log linking.log *~ */*~ 
