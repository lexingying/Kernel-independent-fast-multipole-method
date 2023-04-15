all: libkifmm3d.a

include ../makefile.opt

LIB_SRC = 	kernel3d.cpp knlmat3d.cpp dense3d.cpp let3d.cpp fmm3d.cpp matmgnt3d.cpp fmm3d_setup.cpp fmm3d_eval.cpp fmm3d_check.cpp 

LIB_OBJ = 	$(LIB_SRC:.cpp=.o)

TST_SRC = 	tt.cpp

MEX_SRC = 	fmm3d_mex.cpp

DEP     = 	$(LIB_SRC:.cpp=.d)

libkifmm3d.a: 	$(LIB_OBJ)
	$(AR) $(ARFLAGS) ../libkifmm3d.a $(LIB_OBJ)
	$(RANLIB) ../libkifmm3d.a

tt: libkifmm3d.a tt.o
	${CXX} -o tt tt.o ../libkifmm3d.a ${LDFLAGS}

-include $(DEP)

