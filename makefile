CC = g++ -std=c++14 -Wall
TARGET = Sibernetic
RM := rm -rf


#TEST_SOURCES = src/test/owPhysicTest.cpp

SRCEXT := cpp
SRCDIR := src
INCDIR := inc
BUILDDIR = ./Release
BINARYDIR = $(BUILDDIR)/obj
SOURCES = $(wildcard $(SRCDIR)/*.$(SRCEXT))
BINARYTESTDIR = $(BINARYDIR)/test
OBJECTS := $(patsubst $(SRCDIR)/%,$(BINARYDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
OBJECTS += $(BINARYTESTDIR)/owPhysicTest.o 

CPP_DEPS = $(OBJECTS:.o=.d)

LIBS := -lpython2.7 -lGL -lGLU -lOpenCL -lrt -lglut

CXXFLAGS = $(CC)
EXTRA_LIBS := -L/usr/lib64/OpenCL/vendors/amd/ -L/opt/AMDAPP/lib/x86_64/ -L/usr/lib/x86_64-linux-gnu/ 

ffmpeg: LIBS += -lavcodec -lswscale -lavutil
ffmpeg: CXXFLAGS += -DFFMPEG=1
ffmpeg: $(TARGET)

all: CXXFLAGS += -O3
all : $(TARGET)

debug: CXXFLAGS += -ggdb -O0
debug: $(TARGET)

$(TARGET):$(OBJECTS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker' 
	$(CXXFLAGS)  $(EXTRA_LIBS) -o $(BUILDDIR)/$(TARGET) $(OBJECTS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

$(BINARYDIR)/%.o: $(SRCDIR)/%.cpp 
	@mkdir -p $(BINARYDIR)
	@mkdir -p $(BINARYTESTDIR)
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXXFLAGS) -I/usr/include/python2.7 -I/opt/AMDAPPSDK-3.0/include/ -I/opt/AMDAPP/include/ -I$(INCDIR) -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean : 
	-$(RM) $(OBJECTS)$(CPP_DEPS) $(BUILDDIR)/$(TARGET)
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:
