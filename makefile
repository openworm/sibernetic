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
PYTHON_CONFIG ?= /usr/bin/python3.7-config

CPP_DEPS = $(OBJECTS:.o=.d)

LIBS := $(shell $(PYTHON_CONFIG) --embed --libs) -lGL -lGLU -lOpenCL -lrt -lglut

CXXFLAGS = $(CC) $(shell $(PYTHON_CONFIG) --embed --cflags) -fPIE
EXTRA_LIBS := -L/usr/lib64/OpenCL/vendors/amd/ -L/opt/AMDAPP/lib/x86_64/ -L/usr/lib/x86_64-linux-gnu/ 

ifeq ($(FFMPEG),true)
LIBS += -lavcodec -lswscale -lavutil
CXXFLAGS += -DFFMPEG=1 -I/usr/include/x86_64-linux-gnu
endif

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
	$(CXXFLAGS) -I/opt/AMDAPPSDK-3.0/include/ -I/opt/AMDAPP/include/ -I$(INCDIR) -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean : 
	-$(RM) $(OBJECTS)$(CPP_DEPS) $(BUILDDIR)/$(TARGET)
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:
