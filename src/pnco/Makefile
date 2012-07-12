PNCO_HOME:=/home/yifeng/Workspace/C/pnco
HDF5_HOME:=$(PNCO_HOME)/lib/hdf5
NETCDF4_HOME:=$(PNCO_HOME)/lib/netcdf4
PNETCDF_HOME:=$(PNCO_HOME)/lib/pnetcdf
ZLIB_HOME:=$(PNCO_HOME)/lib/zlib
#HEAD:=-I./ -I$(ZLIB_HOME)/include -I$(HDF5_HOME)/include -I$(NETCDF4_HOME)/include -I$(PNETCDF_HOME)/include
HEAD:=-I./ -I$(ZLIB_HOME)/include -I$(HDF5_HOME)/include -I$(NETCDF4_HOME)/include
LIB:=-L$(NETCDF4_HOME)/lib -lnetcdf -L$(PNETCDF_HOME)/lib -lpnetcdf -L$(HDF5_HOME)/lib -lhdf5_hl -lhdf5 -L$(ZLIB_HOME)/lib -lz -lcurl -lm 
# LIB:=-L$(NETCDF4_HOME)/lib -lnetcdf -L$(PNETCDF_HOME)/lib -L$(HDF5_HOME)/lib -lhdf5_hl -lhdf5 -L$(ZLIB_HOME)/lib -lz -lcurl -lm 
#SET LD_PATH:=-Wl rpath=$(NETCDF4_HOME)/lib -L$(NETCDF4_HOME)/lib
CC:=mpicc
#CC:=gcc
CFLAGS:= -g $(HEAD) $(LIB)
#Find files containing main function
GREPCMD:=grep -l 'int main\s*(.*)' *.c

EXESOURCES:=$(shell $(GREPCMD))
EXEOBJECTS:=$(patsubst %,c,%.o,$(EXEOBJECTS))
ALL_BIN:=$(patsubst %.c,%,$(EXESOURCES))
GREPCMD:=grep -L 'int main\s*(.*)' *.c
NONEXESOURCES:=$(shell $(GREPCMD))
NONEXEOBJECTS:=$(patsubst %.c,%.o,$(NONEXESOURCES))
ALLSOURCE:=$(wildcard *.c)
ALLOBJS:=$(patsubst %.c,%.o,$(SOURCE))
#ALL_BIN:=$(patsubst %.o,%,$(OBJS))


all:$(ALL_BIN)

$(ALL_BIN):$(NONEXEOBJECTS)
	$(CC) $(patsubst %,%.c,$@) $^ -o $@ $(CFLAGS) 
#	$(CC) $(patsubst %,%.o,$@) $^ -o $@ $(CFLAGS) 
%.o:%.c
	$(CC) -c $< -o $@ $(HEAD)
.PHONY:clean
LSCMD= ls *.o $(ALL_BIN)
EXISTOBJ=$(shell $(LSCMD))
clean:
	rm $(EXISTOBJ)




