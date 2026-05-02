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

ARCH := $(shell uname -m)
ifeq ($(ARCH),arm64)
SOURCES := $(filter-out $(SRCDIR)/owOpenCLSolver.cpp,$(SOURCES))
endif

BINARYTESTDIR = $(BINARYDIR)/test
OBJECTS := $(patsubst $(SRCDIR)/%,$(BINARYDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
OBJECTS += $(BINARYTESTDIR)/owPhysicTest.o
PYTHON_VER_MAIN = $(shell python3 -c 'import sys; vv=sys.version_info[:];print(str(vv[0])+str(1.0)[1]+str(vv[1]))')
#PYTHON_VER_MAIN = 3.8 # Hardcode if necessary
PYTHON_CONFIG ?= /usr/bin/python$(PYTHON_VER_MAIN)-config

# NumPy include directory for fast array access (auto-detect if not set;
# can be overridden from the env to point at a venv install).
NUMPYHEADERDIR ?= $(shell python3 -c "import numpy; print(numpy.get_include())" 2>/dev/null)

CPP_DEPS = $(OBJECTS:.o=.d)

LIBS := -lGL -lGLU -lrt -lglut
ifeq ($(ARCH),arm64)
CXXFLAGS += -DOW_NO_OPENCL
else
LIBS += -lOpenCL
CXXFLAGS += -DCL_TARGET_OPENCL_VERSION=120 # Target OpenCL version is 1.2
CXXFLAGS += -DCL_HPP_TARGET_OPENCL_VERSION=120 # C++ bindings target OpenCL version 1.2
CXXFLAGS += -DCL_HPP_MINIMUM_OPENCL_VERSION=120 # Minimum OpenCL version supported is 1.2
CXXFLAGS += -DCL_HPP_ENABLE_PROGRAM_CONSTRUCTION_FROM_ARRAY_COMPATIBILITY # Enable compatibility for program construction from arrays
endif
PYTHON_CONFIG_BASENAME=$(basename $(PYTHON_CONFIG))
PYTHON_VERSION=$(patsubst python%-config,%,$(PYTHON_CONFIG_BASENAME))

ifneq (,$(findstring 3.6,$(PYTHON_CONFIG)))
LIBS += $(shell $(PYTHON_CONFIG) --libs)
CXXFLAGS += $(shell $(PYTHON_CONFIG) --cflags)
else ifneq (,$(findstring 3.7,$(PYTHON_CONFIG)))
LIBS += $(shell $(PYTHON_CONFIG) --libs)
CXXFLAGS += $(shell $(PYTHON_CONFIG) --cflags)
else
LIBS += $(shell $(PYTHON_CONFIG) --embed --libs)
CXXFLAGS += $(shell $(PYTHON_CONFIG) --embed --cflags)
endif

CXXFLAGS += -fPIE
EXTRA_LIBS := -L/usr/lib64/OpenCL/vendors/amd/ -L/opt/AMDAPP/lib/x86_64/ -L/usr/lib/x86_64-linux-gnu/

all: CXXFLAGS += -O3
all : $(TARGET)

debug: CXXFLAGS += -ggdb -O0
debug: $(TARGET)

$(TARGET):$(OBJECTS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	$(CC) $(CXXFLAGS)  $(EXTRA_LIBS) -o $(BUILDDIR)/$(TARGET) $(OBJECTS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

$(BINARYDIR)/%.o: $(SRCDIR)/%.cpp
	@echo 'Assuming Python: $(PYTHON_VER_MAIN)'
	@mkdir -p $(BINARYDIR)
	@mkdir -p $(BINARYTESTDIR)
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CC) $(CXXFLAGS) -I/opt/AMDAPPSDK-3.0/include/ -I/opt/AMDAPP/include/ -I$(NUMPYHEADERDIR) -I$(INCDIR) -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

clean :
	@echo 'Assuming Python: $(PYTHON_VER_MAIN)'
	-$(RM) $(OBJECTS)$(CPP_DEPS) $(BUILDDIR)/$(TARGET)
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:
