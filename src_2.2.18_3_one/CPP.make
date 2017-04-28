#
# --- HYCOM 2.0 surce code makefile 
#
# --- Tunable parameters in ../config/$(ARCH)_$(TYPE)
#

.SUFFIXES: 
.SUFFIXES: .c .F .f .o

.F:
	@echo "Must have an explicit rule for" $*
.f:
	@echo "Must have an explicit rule for" $*
.c:
	@echo "Must have an explicit rule for" $*

include ../config/$(ARCH)_$(TYPE)

.F.f:
	$(RM) $<.f $<.C
	sed -e 's? */// *?/ / /?g' -e 's? *// *?/ /?g' $< >  $<.C
	$(CPP) $(CPPFLAGS) $<.C | sed -e '/^ *$$/d' > $<.f
	-\mv $<.f $*.f
	$(RM) $<.C

default: hycom

machine.f: machine.F
wtime.f:   wtime.F

mod_pipe.f: mod_pipe.F  
mod_xc.f:   mod_xc.F  mod_xc_sm.F mod_xc_mp.F
mod_za.f:   mod_za.F  mod_za_sm.F mod_za_mp.F mod_za_mp1.F
