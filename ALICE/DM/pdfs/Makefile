CFLAGS=`root-config --cflags`
LDFLAGS=`root-config --ldflags --glibs` -lRooFit -lRooFitCore -lMinuit -lopencv_flann

#temp : temp.cpp HWWLVJRooPdfs.cxx HWWLVJRooPdfsDict.cxx libHWWLVJRooPdfs.so
#	g++ -o temp temp.cpp HWWLVJRooPdfs.cxx HWWLVJRooPdfsDict.cxx $(CFLAGS) $(LDFLAGS)

libHWWLVJRooPdfs.so : libHWWLVJRooPdfs.so.1.0
	ln -sf libHWWLVJRooPdfs.so.1.0 libHWWLVJRooPdfs.so

libHWWLVJRooPdfs.so.1.0 : HWWLVJRooPdfs.o HWWLVJRooPdfsDict.o
	gcc -shared -Wl,-soname,libHWWLVJRooPdfs.so.1 -o libHWWLVJRooPdfs.so.1.0 HWWLVJRooPdfs.o HWWLVJRooPdfsDict.o

HWWLVJRooPdfsDict.o : HWWLVJRooPdfsDict.cxx
	g++ -c HWWLVJRooPdfsDict.cxx -fPIC $(CFLAGS) $(LDFLAGS)

HWWLVJRooPdfsDict.cxx : HWWLVJRooPdfs.h LinkDef.h
	rootcint -f HWWLVJRooPdfsDict.cxx -c HWWLVJRooPdfs.h LinkDef.h

HWWLVJRooPdfs.o : HWWLVJRooPdfs.cxx HWWLVJRooPdfs.h
	g++ -c HWWLVJRooPdfs.cxx -fPIC $(CFLAGS) $(LDFLAGS)
