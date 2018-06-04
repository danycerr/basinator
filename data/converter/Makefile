SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:.cpp=.o)

all: basin
	echo $(OBJS)

basin: $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<
