# swap the following two lines to change between debug and release compilers
CXX=c++ -std=c++11 -g  # debug
CXX=c++ -std=c++11 -O3 # release

# compiler options
CXXFLAGS=`pkg-config opencv --cflags`
LDFLAGS=`pkg-config opencv --libs`

# no need to change anything below this line
OBJ=CradleFunctions.o DWT.o FDCT.o FFST.o HaarDWT.o MCA.o TextureRemoval.o
BIN=mainCradleRemoval mainTextureRemoval

all: $(BIN)

$(BIN) : % : %.o $(OBJ)
	$(CXX) $^ -o $@ $(LDFLAGS)

clean:
	rm -f *.o $(BIN)
