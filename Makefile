IDIR=inc
VPATH=src
ODIR=src/obj
LDIR=lib
BINDIR=bin
MYBIN=/home/ymaan/pSofts/bin/

CC=gcc
CFLAGS=-I$(IDIR) -Wno-unused-result -O3 -march=native
LIBS=-lm

_DEPS = header.h  pcBeam_gmrt.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = median.o read_block.o strings_equal.o pack_unpack.o pcBeam_gmrt.o scaledata.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c  $(DEPS)
	@mkdir -p $(BINDIR)
	$(CC) -c -o $@ $< $(CFLAGS)

$(BINDIR)/pcBeam_gmrt: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: install
install:$(BINDIR)/pcBeam_gmrt
	install $(BINDIR)/pcBeam_gmrt $(MYBIN)/

.PHONY: clean
clean:
	rm -f $(ODIR)/*.o *~ core


