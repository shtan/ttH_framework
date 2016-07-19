# System type
UNAME := $(shell uname)

#############################################
### Project location and makefile commons ###
#############################################

PROJECT_DIR=./
BIN_DIR=./bin

include $(PROJECT_DIR)/Makefile_includes
include $(PROJECT_DIR)/Makefile_flags




#######################
### Target location ###
#######################

# Local BIN_DIR
ifeq ($(strip $(BIN_DIR)),)
	BIN_DIRL=.
else
	BIN_DIRL=$(BIN_DIR)
endif




################################
### Include objects to build ###
################################

include Makefile_objects




###################
### All targets ###
###################

test1=test1
test2=test2
test2p1=test2p1
test2p2=test2p2

all: $(test2)




######################
### Building block ###
######################

$(test1): $(objects) $(objects1) src/$(test1).o
	$(CXX) $(CXXFLAGS) -o $(BIN_DIRL)/$(test1) src/$(test1).o $(objects1) \
	$(objects) -L. $(LIBFLAGS) -l MEKD

$(test2): $(objects) $(objects2) src/$(test2).o
	$(CXX) $(CXXFLAGS) -o $(BIN_DIRL)/$(test2) src/$(test2).o $(objects2) \
	$(objects) -L. $(LIBFLAGS) -l ttbar_reco_C

$(test2p1): $(objects) $(objects2p1) src/$(test2p1).o
	$(CXX) $(CXXFLAGS) -o $(BIN_DIRL)/$(test2p1) src/$(test2p1).o \
	$(objects2p1) $(objects) -L. $(LIBFLAGS) -l MEKD

$(test2p2): $(objects) $(objects2p2) src/$(test2p2).o
	$(CXX) $(CXXFLAGS) -o $(BIN_DIRL)/$(test2p2) src/$(test2p2).o \
	$(objects2p2) $(objects) -L. $(LIBFLAGS) -l ttbar_reco_C -l MEKD


.PHONY: clean

clean_0:
	rm -f $(objects)

clean_1:
	rm -f $(objects1); \
	rm -f src/$(test1).o; \
	rm -f $(BIN_DIRL)/$(test1)

clean_2:
	rm -f $(objects2); \
	rm -f src/$(test2).o; \
	rm -f $(BIN_DIRL)/$(test2)

clean_2p1:
	rm -f $(objects2p1); \
	rm -f src/$(test2p1).o; \
	rm -f $(BIN_DIRL)/$(test2p1)

clean_2p2:
	rm -f $(objects2p2); \
	rm -f src/$(test2p2).o; \
	rm -f $(BIN_DIRL)/$(test2p2)

clean: clean_1 clean_2 clean_2p1 clean_2p2
