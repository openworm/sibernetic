################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/PyramidalSimulation.cpp \
../src/main.cpp \
../src/owHelper.cpp \
../src/owOpenCLSolver.cpp \
../src/owPhysicsFluidSimulator.cpp \
../src/owWorldSimulation.cpp 

OBJS += \
./src/PyramidalSimulation.o \
./src/main.o \
./src/owHelper.o \
./src/owOpenCLSolver.o \
./src/owPhysicsFluidSimulator.o \
./src/owWorldSimulation.o 

CPP_DEPS += \
./src/PyramidalSimulation.d \
./src/main.d \
./src/owHelper.d \
./src/owOpenCLSolver.d \
./src/owPhysicsFluidSimulator.d \
./src/owWorldSimulation.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/include/python2.7/ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


