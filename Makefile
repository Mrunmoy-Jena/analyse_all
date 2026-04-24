CXX=g++
LINKER=g++
EXEC=analysis

ROOTCONFIG := root-config

CFLAGS := $(shell $(ROOTCONFIG) --cflags)
CFLAGS += -I${FAIRROOTPATH}/include
CFLAGS += -I$(SIMPATH)/include
CFLAGS += -I$(ROOT_INCLUDE_PATH)
CFLAGS += -I$(ROOT_INCLUDE_DIR)
CFLAGS += -I$(VMCWORKDIR)
CFLAGS += -I$(VMCWORKDIR)/r3bdata/footData
CFLAGS += -I$(VMCWORKDIR)/r3bdata/wrData
CFLAGS += -I$(VMCWORKDIR)/r3bdata/frsData
CFLAGS += -I$(VMCWORKDIR)/r3bdata/tofData
CFLAGS += -I$(VMCWORKDIR)/r3bdata/califaData
CFLAGS += -I$(VMCWORKDIR)/r3bdata/neulandData
CFLAGS += -I$(VMCWORKDIR)/r3bdata
CFLAGS += -I$(VMCWORKDIR)/r3bsource/base
CFLAGS += -I$(VMCWORKDIR)/r3bbase
CFLAGS += -I$(VMCWORKDIR)/tracking

#CFLAGS += -I$(UCESB_DIR)/hbook

CFLAGS += --std=c++17 -g -O0 -fexceptions

#CFLAGS += -g -Wall -W -Wconversion -Wshadow -Wcast-qual -Wwrite-strings 

LDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
LDFLAGS += -lEG $(shell $(ROOTCONFIG) --glibs)
LDFLAGS += -L$(ROOT_LIBRARY_DIR) -L$(FAIRROOTPATH)/lib
LDFLAGS += -g

#LDFLAGS += -L$(VMCWORKDIR)/../build/lib -lR3BSsd -lR3Bsource -lR3BBase

LDFLAGS += -L$(VMCWORKDIR)/../build/lib -lR3BSsd -lR3BBase -lR3BData -lR3BTracking
LDFLAGS += -L$(FAIRROOTPATH)/lib -lBase -lParBase -lFairTools
LDFLAGS += -lR3BData

INCLUDEDIR=include

DIR_INC=-I$(INCLUDEDIR)

SRC=analyse_all_WR.C
OBJ=$(SRC:.C=.o)
OBJ:=$(OBJ:.cxx=.o)

all: $(EXEC)

$(EXEC): $(OBJ)
	${LINKER} -o $@ $^ ${LDFLAGS}
	echo " COMP $@"

%.o : %.C
	$(MAKEDEPEND)
	${CXX} ${CFLAGS} $(DIR_INC) -c $< -o $@
	echo "	CXX $@"


clean:
	rm -f *.o $(EXEC)
