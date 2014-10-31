CC=/usr/bin/c++

INCLUDE =  -I/home/feixh/workspace/Theia\
	   -I/home/feixh/workspace/Theia/include\
	   -I/home/feixh/workspace/Theia/src\
	   -I/usr/local/include\
	   -I/usr/include/eigen3\
	   -I/home/feixh/workspace/Theia/libraries\
	   -I/home/feixh/workspace/Theia/libraries/agast\
	   -I/home/feixh/workspace/Theia/libraries/cimg\
	   -I/home/feixh/workspace/Theia/libraries/optimo\
           -I/home/feixh/workspace/Theia/libraries/statx\
	   -I/home/feixh/workspace/Theia/libraries/stlplus3\
	   -I/home/feixh/workspace/Theia/libraries/vlfeat


CFLAGS = -fopenmp -std=c++11 -Werror=all -Werror=extra -Wno-unknown-pragmas -Wno-sign-compare -Wno-unused-parameter -Wno-missing-field-initializers -Wno-comment -Wno-unused-variable -Wno-unused-result -Wno-unused-but-set-variable -O3 -DNDEBUG     


LIB_PREFIX = /home/feixh/workspace/theia-build/lib
LIBS =  -rdynamic $(LIB_PREFIX)/libgtest.a $(LIB_PREFIX)/libtheia.a /usr/local/lib/libceres.a -lgflags -lglog $(LIB_PREFIX)/libcimg.a $(LIB_PREFIX)/libstatx.a $(LIB_PREFIX)/libstlplus3.so $(LIB_PREFIX)/libvlfeat.so $(LIB_PREFIX)/libagast.a -lgomp -lpng -lz -ljpeg -ltiff -lz -ljpeg -ltiff /usr/local/lib/libceres.a -lspqr -ltbb -ltbbmalloc -lcholmod -lccolamd -lcamd -lcolamd -lamd -Wl,-Bstatic -lsuitesparseconfig -Wl,-Bdynamic -lrt -lpthread -lglog -llapack -lf77blas -latlas -lgomp -Wl,-rpath,/home/feixh/workspace/theia-build/lib
 
	



test: 
	$(CC) test.cpp $(INCLUDE) $(LIBS) $(CFLAGS) 

	

