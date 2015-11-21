CXX = /usr/local/bin/g++
CXXFLAGS = -g -Wall -DDEBUG
#CXXFLAGS = -O3 -Wall

LOBJS =		Matrix.o BMatrix.o TMatrix.o

LIBS = libTileMatrix.a

TARGET =	test

$(LIBS):	$(LOBJS)
	$(AR) r $(LIBS) $(LOBJS)
	ranlib $(LIBS)

$(TARGET):	TileMatrixTest.o $(LIBS)
	$(CXX) $(CXXFLAGS) -o $@ TileMatrixTest.o $(LIBS)

all:	$(LIBS)

clean:
	rm -f $(LOBJS) $(LIBS)

.cpp.o :
	$(CXX) $(CXXFLAGS) -c $<
