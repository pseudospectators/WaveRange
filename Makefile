include config.mk

OUTPUTDIR = ./bin/
AUXDIR = ./libc/
MKDIR = mkdir -p  

SRCCOM = src/core/wrappers.cpp
OBJCOM = wrappers.o  

SRCFLU = src/flusi/hdf5_interfaces.cpp
OBJFLU = hdf5_interfaces.o

SRCMSS = src/mssg/ctrl_aux.cpp 
OBJMSS = ctrl_aux.o


all: flusi mssg generic

flusi: common
	$(MKDIR) $(OUTPUTDIR)
	$(MKDIR) $(OUTPUTDIR)flusi
	$(CXX) -c -I$(HDF_INC) $(CXXFLAGS) $(SRCFLU)
	$(CXX) -I$(HDF_INC) $(CXXFLAGS) ./src/flusi/main_enc.cpp $(OBJCOM) $(OBJFLU) -L$(AUXDIR) -L$(HDF_LIB) -lcoders -lwcdf -lhdf5 -o $(OUTPUTDIR)flusi/wrenc
	$(CXX) -I$(HDF_INC) $(CXXFLAGS) ./src/flusi/main_dec.cpp $(OBJCOM) $(OBJFLU) -L$(AUXDIR) -L$(HDF_LIB) -lcoders -lwcdf -lhdf5 -o $(OUTPUTDIR)flusi/wrdec
mssg: common
	$(MKDIR) $(OUTPUTDIR)
	$(MKDIR) $(OUTPUTDIR)mssg
	$(CXX) -c -I$(HDF_INC) $(CXXFLAGS) $(SRCMSS)
	$(CXX) $(CXXFLAGS) ./src/mssg/mssg_enc.cpp $(OBJCOM) $(OBJMSS) -L$(AUXDIR) -lcoders -lwcdf -o $(OUTPUTDIR)mssg/wrmssgenc
	$(CXX) $(CXXFLAGS) ./src/mssg/mssg_dec.cpp $(OBJCOM) $(OBJMSS) -L$(AUXDIR) -lcoders -lwcdf -o $(OUTPUTDIR)mssg/wrmssgdec
generic: generic
	$(MKDIR) $(OUTPUTDIR)
	$(MKDIR) $(OUTPUTDIR)generic
	$(CXX) -c -I$(HDF_INC) $(CXXFLAGS) 
	$(CXX) $(CXXFLAGS) ./src/mssg/gen_enc.cpp $(OBJCOM) -L$(AUXDIR) -lcoders -lwcdf -o $(OUTPUTDIR)generic/wrenc
	$(CXX) $(CXXFLAGS) ./src/mssg/gen_dec.cpp $(OBJCOM) -L$(AUXDIR) -lcoders -lwcdf -o $(OUTPUTDIR)generic/wrdec
common:
	$(MKDIR) $(AUXDIR)
	cd ./src/waveletcdf97_3d && $(MAKE) all
	cd ./src/rangecod && $(MAKE) all
	$(CXX) -c -I$(HDF_INC) $(CXXFLAGS) $(SRCCOM)
.PHONY: clean
clean:
	cd ./src/waveletcdf97_3d && $(MAKE) clean
	cd ./src/rangecod && $(MAKE) clean
	$(RM) -rf $(OUTPUTDIR)
	$(RM) -rf $(AUXDIR)
	$(RM) ./*.gc??
	$(RM) ./*.o
