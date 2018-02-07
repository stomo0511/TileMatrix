UNAME = $(shell uname)
  CXX = g++
ifeq ($(UNAME),Linux)
endif
ifeq ($(UNAME),Darwin)
  CXX = /usr/local/bin/g++
endif


#CXXFLAGS = -g -Wall -DDEBUG
CXXFLAGS = -O3 -Wall

LOBJS =		Matrix.o BMatrix.o TMatrix.o DMatrix.o

LIBS = libTileMatrix.a

TARGET =	test

$(LIBS):	$(LOBJS)
	$(AR) r $(LIBS) $(LOBJS)
	ranlib $(LIBS)

$(TARGET):	TileMatrixTest.o $(LIBS)
	$(CXX) $(CXXFLAGS) -o $@ TileMatrixTest.o $(LIBS)

all:	$(LIBS) $(TARGET)

clean:
	rm -f *.o $(LIBS) $(TARGET)

.cpp.o :
	$(CXX) $(CXXFLAGS) -c $<
