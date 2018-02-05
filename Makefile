########################################################
#
#  Makefile for drsosc, drscl and drs_exam 
#  executables under linux
#
#  Requires wxWidgets 2.8.9 or newer
#
########################################################

# check if wxWidgets is installed
HAVE_WX       = $(shell which wx-config)
ifeq ($(HAVE_WX),)
$(error Error: wxWidgets required to compile "drsosc")
endif

# check for OS
OS            = $(shell uname)
ifeq ($(OS),Darwin)
DOS           = OS_DARWIN
else
DOS           = OS_LINUX
endif

#CC=g++
CFLAGS        = -g -O2 -Wall -Wuninitialized -fno-strict-aliasing -Iinclude -I/usr/local/include -D$(DOS) -DHAVE_USB -DHAVE_LIBUSB10 -DUSE_DRS_MUTEX  
LIBS          = -lpthread -lutil -lusb-1.0
LDFLAGS       = `root-config --glibs`
ROOTFLAGS     = `root-config --cflags`

ifeq ($(OS),Darwin)
CFLAGS        += -stdlib=libstdc++
endif         

# wxWidgets libs and flags
WXLIBS        = $(shell wx-config --libs)
WXFLAGS       = $(shell wx-config --cxxflags)

CPP_OBJ       = DRS.o averager.o ConfigDialog.o DOFrame.o DOScreen.o DRSOsc.o MeasureDialog.o Measurement.o Osci.o InfoDialog.o DisplayDialog.o AboutDialog.o EPThread.o TriggerDialog.o rb.o  functions.o
OBJECTS       = musbstd.o mxml.o strlcpy.o


ifeq ($(OS),Darwin)
all: drsosc drscl drs_exam drs_exam_multi DRSOsc.app drs_reading 
else
all: drsosc drscl drs_exam drs_exam_multi drs_reading 
endif

DRSOsc.app: drsosc
	-mkdir DRSOsc.app
	-mkdir DRSOsc.app/Contents
	-mkdir DRSOsc.app/Contents/MacOS
	-mkdir DRSOsc.app/Contents/Resources
	-mkdir DRSOsc.app/Contents/Resources/English.lproj
	echo 'APPL????' > DRSOsc.app/Contents/PkgInfo
	cp Info.plist DRSOsc.app/Contents/Info.plist
	cp DRSOsc.icns DRSOsc.app/Contents/Resources
	cp drsosc DRSOsc.app/Contents/MacOS/DRSOsc

drsosc: $(OBJECTS) $(CPP_OBJ) main.o
	$(CXX) $(CFLAGS)  $(ROOTFLAGS) $(OBJECTS) $(CPP_OBJ) main.o -o drsosc $(LIBS) $(WXLIBS) $(LDFLAGS)

drscl: $(OBJECTS) DRS.o averager.o drscl.o
	$(CXX) $(CFLAGS) $(OBJECTS) DRS.o averager.o drscl.o -o drscl $(LIBS) $(WXLIBS)

drs_exam: $(OBJECTS) DRS.o averager.o drs_exam.o
	$(CXX) $(CFLAGS) $(OBJECTS) DRS.o averager.o drs_exam.o -o drs_exam $(LIBS) $(WXLIBS)
	
drs_reading: $(OBJECTS) DRS.o averager.o drs_reading.o
	$(CXX) $(CFLAGS)  $(ROOTFLAGS) $(OBJECTS) DRS.o averager.o drs_reading.o  functions.o -o drs_reading $(LIBS) $(WXLIBS)  $(LDFLAGS)
		
drs_exam_multi: $(OBJECTS) DRS.o averager.o drs_exam_multi.o
	$(CXX) $(CFLAGS) $(OBJECTS) DRS.o averager.o drs_exam_multi.o -o drs_exam_multi $(LIBS) $(WXLIBS)

main.o: src/main.cpp include/mxml.h include/DRS.h
	$(CXX) $(CFLAGS) $(WXFLAGS) -c $<

drscl.o: src/drscl.cpp include/mxml.h include/DRS.h
	$(CXX) $(CFLAGS) -c $<

drs_exam.o: src/drs_exam.cpp include/mxml.h include/DRS.h
	$(CXX) $(CFLAGS) -c $<
	
drs_reading.o: src/drs_reading.cpp src/functions.cpp include/mxml.h include/DRS.h include/functions.h
	$(CXX) $(CFLAGS) $(ROOTFLAGS) -c $<

drs_exam_multi.o: src/drs_exam_multi.cpp include/mxml.h include/DRS.h
	$(CXX) $(CFLAGS) -c $<

$(CPP_OBJ): %.o: src/%.cpp include/%.h include/mxml.h include/DRS.h
	$(CXX) $(CFLAGS) $(WXFLAGS) $(ROOTFLAGS) -c $<

$(OBJECTS): %.o: src/%.c include/mxml.h include/DRS.h
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o drscl drsosc drs_reading
