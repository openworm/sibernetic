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


all : $(TARGET)

$(TARGET):$(OBJECTS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ -L/opt/AMDAPP/lib/x86_64/  -o $(BUILDDIR)/$(TARGET) $(OBJECTS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

$(BINARYDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(BINARYDIR)
	@mkdir -p $(BINARYTESTDIR)
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/include/python2.7 -I$(INCDIR) -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean :
	-$(RM) $(OBJECTS)$(CPP_DEPS) $(BUILDDIR)/$(TARGET)
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:
