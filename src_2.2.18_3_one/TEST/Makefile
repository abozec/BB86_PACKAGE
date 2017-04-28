#
# --- HYCOM 1.1 TEST makefile, in subdirectory ./TEST
#
# --- primary HYCOM .o files are assumed to be up to date.
#
# --- Tunable parameters in ../../config/$(ARCH)_$(TYPE)
#

.SUFFIXES: 
.SUFFIXES: .c .F .f .o

include ../../config/$(ARCH)_$(TYPE)

default: test_xcl

# create handlinks to ../mod*.o and to machine-specific module info files
MODS =	mod_xc.o mod_za.o

# assumed up to date:
OBJS =	../machine.o ../machi_c.o ../wtime.o

xca:	test_xca.o
	$(LD)  $(LDFLAGS) -o test_xca  test_xca.o  $(OBJS) $(MODS)

xcl:	test_xcl.o
	$(LD)  $(LDFLAGS) -o test_xcl  test_xcl.o  $(OBJS) $(MODS)

xcs:	test_xcs.o
	$(LD)  $(LDFLAGS) -o test_xcs  test_xcs.o  $(OBJS) $(MODS)

xct:	test_xct.o
	$(LD)  $(LDFLAGS) -o test_xct  test_xct.o  $(OBJS) $(MODS)

xct_arctic: test_xct_arctic.o
	$(LD)  $(LDFLAGS) -o test_xct_arctic test_xct_arctic.o  $(OBJS) $(MODS)

zaio:	test_zaio.o
	$(LD)  $(LDFLAGS) -o test_zaio test_zaio.o $(OBJS) $(MODS)

test_xca.o:     	test_xca.f		dimensions.h
test_xcl.o:     	test_xcl.f		dimensions.h
test_xct.o:     	test_xct.f		dimensions.h
test_xct_arctic.o:	test_xct_arctic.f	dimensions.h
test_zaio.o:		test_zaio.f		dimensions.h
