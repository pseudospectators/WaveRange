include config.mk

OUTPUTDIR = ./bin/
AUXDIR = ./libc/

MKDIR = mkdir -p  

SRCFLU = src/flusi/hdf5_interfaces.cpp
OBJFLU = hdf5_interfaces.o

SRCMSS = src/mssg/ctrl_aux.cpp 
OBJMSS = ctrl_aux.o

SRCGEN = src/generic/gen_aux.cpp
OBJGEN = gen_aux.o

all: generic mssg flusi

flusi: common
	$(MKDIR) $(OUTPUTDIR)flusi
	$(CXX) -c -I$(HDF_INC) $(CXXFLAGS) $(SRCFLU)
	$(CXX) -I$(HDF_INC) $(CXXFLAGS) ./src/flusi/main_enc.cpp $(OBJFLU) -L$(AUXDIR) -L$(HDF_LIB) -lwaverange -lhdf5 -o $(OUTPUTDIR)flusi/wrenc
	$(CXX) -I$(HDF_INC) $(CXXFLAGS) ./src/flusi/main_dec.cpp $(OBJFLU) -L$(AUXDIR) -L$(HDF_LIB) -lwaverange -lhdf5 -o $(OUTPUTDIR)flusi/wrdec
mssg: common
	$(MKDIR) $(OUTPUTDIR)mssg
	$(CXX) -c -I$(HDF_INC) $(CXXFLAGS) $(SRCMSS)
	$(CXX) $(CXXFLAGS) ./src/mssg/mssg_enc.cpp $(OBJMSS) -L$(AUXDIR) -lwaverange -o $(OUTPUTDIR)mssg/wrmssgenc
	$(CXX) $(CXXFLAGS) ./src/mssg/mssg_dec.cpp $(OBJMSS) -L$(AUXDIR) -lwaverange -o $(OUTPUTDIR)mssg/wrmssgdec
generic: common
	$(MKDIR) $(OUTPUTDIR)generic
	$(CXX) -c -I$(HDF_INC) $(CXXFLAGS) $(SRCGEN)
	$(CXX) $(CXXFLAGS) ./src/generic/gen_enc.cpp $(OBJGEN) -L$(AUXDIR) -lwaverange -o $(OUTPUTDIR)generic/wrenc
	$(CXX) $(CXXFLAGS) ./src/generic/gen_dec.cpp $(OBJGEN) -L$(AUXDIR) -lwaverange -o $(OUTPUTDIR)generic/wrdec
common:
	$(MKDIR) $(OUTPUTDIR)
	$(MKDIR) $(OUTPUTDIR)lib/
	$(MKDIR) $(AUXDIR)
	cd ./src/waveletcdf97_3d && $(MAKE) all
	cd ./src/rangecod && $(MAKE) all
	cd ./src/core && $(MAKE) all
	cp $(AUXDIR)libwaverange.a $(OUTPUTDIR)lib/ 
.PHONY: clean
clean:
	cd ./src/waveletcdf97_3d && $(MAKE) clean
	cd ./src/rangecod && $(MAKE) clean
	cd ./src/core && $(MAKE) clean
	$(RM) -rf $(OUTPUTDIR)
	$(RM) -rf $(AUXDIR)
	$(RM) ./*.gc??
	$(RM) ./*.o
