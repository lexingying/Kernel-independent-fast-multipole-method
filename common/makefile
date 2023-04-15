all: libkifmm3d.a

include ../makefile.opt

LIB_SRC = 	vecmatop.cpp	svdrep.cpp

LIB_OBJ = 	$(LIB_SRC:.cpp=.o)

DEP     = 	$(LIB_SRC:.cpp=.d)

libkifmm3d.a: 	$(LIB_OBJ)
	$(AR) $(ARFLAGS) ../libkifmm3d.a $(LIB_OBJ)
	$(RANLIB) ../libkifmm3d.a

-include $(DEP)

