TARGET = Sibernetic
RM := rm -rf

SOURCES = src/owSignalSimulator.cpp \
src/owVtkExport.cpp \
src/main.cpp \
src/owHelper.cpp \
src/owOpenCLSolver.cpp \
src/owPhysicsFluidSimulator.cpp \
src/owWorldSimulation.cpp \
src/owNeuronSimulator.cpp

TEST_SOURCES = src/test/owPhysicTest.cpp

SRCEXT := cpp
SRCDIR := src
INCDIR := inc
BUILDDIR = ./Release
BINARYDIR = $(BUILDDIR)/obj
BINARYTESTDIR = $(BINARYDIR)/test
OBJECTS := $(patsubst $(SRCDIR)/%,$(BINARYDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
OBJECTS += $(BINARYTESTDIR)/owPhysicTest.o 

CPP_DEPS = $(OBJECTS:.o=.d)

LIBS := -lpython2.7 -lGL -lGLU -lOpenCL -lrt -lglut


CXXCOMPILER = g++
CXXFLAGS = $(CXXCOMPILER)

all : $(TARGET)
all:  CXXFLAGS += -O3 -Wall 
debug: CXXFLAGS += -g -O0
debug: $(TARGET)

$(TARGET):$(OBJECTS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	$(CXXCOMPILER) -L/usr/lib64/OpenCL/vendors/amd/ -o $(BUILDDIR)/$(TARGET) $(OBJECTS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

$(BINARYDIR)/%.o: $(SRCDIR)/%.cpp 
	@mkdir -p $(BINARYDIR)
	@mkdir -p $(BINARYTESTDIR)
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXXFLAGS) -I/usr/include/python2.7 -I/opt/AMDAPPSDK-3.0/include/ -I$(INCDIR) -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

clean : 
	-$(RM) $(OBJECTS) $(CPP_DEPS) $(BUILDDIR)/$(TARGET)
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:
