CC = g++

# Directories
IDIR_MKL = $(MKLROOT)/include
IDIR_HEADER = ../include
IDIR_EIGEN = ../
OBJDIR = obj
LIBDIR = ../lib

CPPFLAGS = -w -std=c++11 -DEIGEN_USE_MKL_ALL \
           -m64 -Ofast -DNDEBUG -fPIC -I$(IDIR_MKL) \
           -I$(IDIR_HEADER) -I$(IDIR_EIGEN)

LINK = -Wl,--no-as-needed -L${MKLROOT}/lib/intel64 \
       -lmkl_intel_lp64 \
       -lmkl_core \
       -lmkl_intel_thread \
       -liomp5 \
       -lpthread \
       -lm \
       -ldl

TARGET = $(LIBDIR)/libmodel.so
OBJ = $(addprefix $(OBJDIR)/, model.o)


$(TARGET): $(OBJ)
	$(CC) $(OBJ) $(LINK) -shared -o $(TARGET) 

$(OBJDIR)/%.o: %.cpp ../include/%.h
	$(CC) $(CPPFLAGS) -c -o $@ $<

.PHONY: clean
clean:
	rm $(TARGET) $(OBJ)