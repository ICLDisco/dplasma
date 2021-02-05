MPICC=mpicc
MPIF77=mpif77
CFLAGS+= #-O0 -g

CC=$(MPICC)
LD=mpif77

MKLDIR=-L${MKLROOT}/lib/intel64
DPLASMA_WRAPPER=../../tests/wrappers

SMKL=$(MKLDIR) -Wl,-Bstatic -Wl,--start-group -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 \
-lmkl_sequential -lmkl_gf_lp64 -lmkl_core -Wl,--end-group -Wl,-Bdynamic -ldl -lm

TMKL=$(MKLDIR) -Wl,-Bstatic -Wl,--start-group -lmkl_scalapack_lp64 -lmkl_gf_lp64 -lmkl_core \
-lmkl_gnu_thread -lmkl_blacs_openmpi_lp64 -Wl,--end-group -Wl,-Bdynamic -lm -fopenmp -ldl

TLIBS=$(TMKL)
SLIBS=$(SMKL)
WLIBS=-L$(DPLASMA_WRAPPER) -Wl,-rpath,$(DPLASMA_WRAPPER) -ldplasma_scalapack $(SLIBS)

TARGETS= \
pdpotrf tpdpotrf pdpotrf_wrapped \
pdgeqrf tpdgeqrf pdgeqrf_wrapped \
pdgetrf tpdgetrf pdgetrf_wrapped \
pdgemm tpdgemm pdgemm_wrapped \
pdtrmm tpdtrmm pdtrmm_wrapped \
pdtrsm tpdtrsm pdtrsm_wrapped \
pdsyev tpdsyev \
pdgebrd tpdgebrd \
pdpotri tpdpotri \


all: $(TARGETS)

common.o: common.c
	$(CC) ${CFLAGS}  -o $@ -c $<

p%.o: p%.c
	$(CC) ${CFLAGS}  -o $@ -c $<

p%: p%.o common.o
	$(LD) $(LDFLAGS) -o $@ $^ $(SLIBS)

t%: %.o common.o
	$(LD) $(LDFLAGS) -o $@ $^ $(TLIBS)


#WRAPPED VERSION
common_wrapped.o: common.c
	$(CC) ${CFLAGS} -DDPLASMA_WRAPPER_ON -o $@ -c $<

p%_wrapped.o: p%.c
	$(CC) ${CFLAGS} -DDPLASMA_WRAPPER_ON -o $@ -c $<

p%_wrapped: p%_wrapped.o common_wrapped.o
	$(LD) $(LDFLAGS) -o $@ $^ $(WLIBS)

clean:
	$(RM) *.o $(TARGETS)

