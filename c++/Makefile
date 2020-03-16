CXX=g++
RM=rm -f
CPPFLAGS=-O3 -I/usr/include/eigen3 -I. \
	-I./helpers \
	-I./floor_f1Hf2 \
	-I./floor_fHf \
	-I./floor_Hf
LDFLAGS=
LDLIBS=
DEPS=$(shell find . -name "*.h")

SRCS=$(shell find . -name "*.cpp")
OBJS=$(subst .cpp,.o,$(SRCS))

all: test

%.o: %.ccp $(DEPS)
	$(CXX) -c -o $@ $< $(CPPFLAGS)

test: $(OBJS)
	$(CXX) $(LDFLAGS) -o test $(OBJS) $(LDLIBS)

clean:
	$(RM) $(OBJS)
	$(RM) test
