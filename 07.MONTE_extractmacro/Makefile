# ------------------------------------------------------------
# This version of the Make file is for the fortran standard library implementation 
# Target CPU (Cray PE)
export CRAY_CPU_TARGET = x86-milan

# ---- COMPILER
FC      = ftn
FFLAGS  = -fopenmp -I$(FORTRANSTDLIB_INCLUDEDIR) -fcheck=all -g -O0

# ---- SUNDIALS (adjust if your install differs)
SUNDIALS_PREFIX     = /data/project/general/ans/COD/SUNDIALS/INSTDIR
SUNDIALS_INCLUDEDIR = $(SUNDIALS_PREFIX)/fortran
SUNDIALS_LIBDIR     = $(SUNDIALS_PREFIX)/lib64

INCLUDES = -I$(SUNDIALS_INCLUDEDIR)
LDFLAGS  = -L$(SUNDIALS_LIBDIR) -Wl,-rpath,$(SUNDIALS_LIBDIR) 
LIBS     = -lsundials_fida_mod -lsundials_ida -lsundials_fcore_mod -lm 

# ---- TARGET
TARGET = MONTE_extractmacro.x

# ---- SOURCES / OBJECTS
SRC = csv_nd_with_index.f90 globals.f90 MONTE_extractmacro.f90  main.f90
# If you have extra SUNDIALS helper modules/sources (e.g., RSDLS.f90), add them here:
# SRC += RSDLS.f90
OBJ = $(SRC:.f90=.o)

# ---- RULES
.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(LDFLAGS) $(LIBS)

%.o: %.f90
	$(FC) $(FFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)
