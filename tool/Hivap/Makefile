OBJS =  alloch.o  barfit.o  beta.o  casout.o  cutoff.o  decodn.o  \
densty.o  dropb.o  dropl.o  engout.o  eva.o  fisrot.o  fuse.o  fusio.o  \
fusion.o  gammas.o  getput.o  grids.o  invpar.o  ip1.o  ip2.o  msben.o  \
myrd.o  om.o  out1.o  out2.o  over2.o  owkb.o  parabo.o  parap.o  popul.o  \
rot.o  rowout.o  sepenf.o  seqin.o  seqout.o  sumpop.o  tabip.o  tld.o  \
tot.o  transm.o  twkb.o  vlro.o  halloc.o 
# Fuer pgf77 VORHER: . pgilogin
# Fuer ifort VORHER: . /opt/intel-login 9.1
# F77 = f77 oder F77 = g77 oder F77 = pgf77 oder F77 = ifort
F77 = f77
# Debug
# FFLAGS = -g -O  -fpedantic  -Wunused  -Wuninitialized
# FFLAGS = -g -ffloat-store
# FFLAGS = -g -fno-trapping-math  -ffinite-math-only  (This option should never be turned on by any -O option )
# FFLAGS = -g -fbounds-check  -finit-local-zero
# FFLAGS = -g -O1
# FFLAGS = -g -O1 -ffloat-store
FFLAGS = -g
# Warnings unsichtbar
# FFLAGS = -g -w
# Optimize
# FFLAGS = -O
#
# Ohne Common-Fehlertexte:
# make 1> run.all 2>&1
# grep -i -v common run.all | grep -v "\^" | grep -v f77
#
.f.o:
	$(F77) -c $(FFLAGS) $<

hivapn: $(OBJS) hivap.f
	$(F77) $(FFLAGS) $(OBJS) hivap.f -o hivapn

clean: 
	rm -f $(OBJS) hivapn
