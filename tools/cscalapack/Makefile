MPICC=mpicc
MPIF77=mpif77
CFLAGS+=

CC=$(MPICC)
LD=mpif77

MKLDIR=-L${MKLROOT}/lib/intel64
DPLASMA_WRAPPER=../../tests/wrappers

SMKL=$(MKLDIR) -Wl,-Bstatic -Wl,--start-group -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 \
-lmkl_sequential -lmkl_gf_lp64 -lmkl_core -Wl,--end-group -Wl,-Bdynamic -ldl -lm

TMKL=$(MKLDIR) -Wl,-Bstatic -Wl,--start-group -lmkl_scalapack_lp64 -lmkl_gf_lp64 -lmkl_core \
-lmkl_gnu_thread -lmkl_blacs_openmpi_lp64 -Wl,--end-group -Wl,-Bdynamic -ldl -lm -fopenmp

TLIBS=$(TMKL)
SLIBS=$(SMKL)

ORIG_SMKL= $(MKLDIR) -Wl,-Bstatic -Wl,--start-group -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 \
-lmkl_sequential -lmkl_gf_lp64 -lmkl_core -Wl,--end-group -Wl,-Bdynamic -ldl
ORIG_TMKL= $(MKLDIR) -Wl,-Bstatic -Wl,--start-group -lmkl_scalapack_lp64 -lmkl_gf_lp64 -lmkl_core \
-lmkl_gnu_thread -lmkl_blacs_openmpi_lp64 -Wl,--end-group -Wl,-Bdynamic -lm -fopenmp -ldl

ORIG_SLIBS=$(ORIG_SMKL)
ORIG_TLIBS=$(ORIG_TMKL)

WRAP_SLIBS=-L$(DPLASMA_WRAPPER) -Wl,-rpath,$(DPLASMA_WRAPPER) -ldplasma_scalapack $(ORIG_SLIBS)
WRAP_TLIBS=-L$(DPLASMA_WRAPPER) -Wl,-rpath,$(DPLASMA_WRAPPER) -ldplasma_scalapack $(ORIG_TLIBS)

TARGETS= \
pdpotrf tpdpotrf \
pdpotrf_wrapped tpdpotrf_wrapped \
pdgeqrf tpdgeqrf \
pdgeqrf_wrapped tpdgeqrf_wrapped \
pdgetrf tpdgetrf \
pdgetrf_wrapped tpdgetrf_wrapped \
pdgemm tpdgemm \
pdgemm_wrapped tpdgemm_wrapped \
pdtrmm tpdtrmm \
pdtrmm_wrapped tpdtrmm_wrapped \
pdtrsm tpdtrsm \
pdtrsm_wrapped tpdtrsm_wrapped

all: $(TARGETS)

common.o: common.c
	$(CC) ${CFLAGS}  -o $@ -c $<

p%.o: p%.c
	$(CC) ${CFLAGS}  -o $@ -c $<

p%: p%.o common.o
	$(LD) $(LDFLAGS) -o $@ $^ $(ORIG_SLIBS)

t%: %.o common.o
	$(LD) $(LDFLAGS) -o $@ $^ $(ORIG_TLIBS)


#WRAPPED VERSION
common_wrapped.o: common.c
	$(CC) ${CFLAGS} -DDPLASMA_WRAPPER_ON -o $@ -c $<

p%_wrapped.o: p%.c
	$(CC) ${CFLAGS} -DDPLASMA_WRAPPER_ON -o $@ -c $<

p%_wrapped: p%_wrapped.o common_wrapped.o
	$(LD) $(LDFLAGS) -o $@ $^ $(WRAP_SLIBS)

t%_wrapped: %_wrapped.o common_wrapped.o
	$(LD) $(LDFLAGS) -o $@ $^ $(WRAP_TLIBS)

clean:
	$(RM) *.o $(TARGETS)

