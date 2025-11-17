################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../DDPModel/DDPDFNFunctions.cpp \
../DDPModel/DDPModel.cpp \
../DDPModel/DDPStudies.cpp \
../DDPModel/ExchangeFunctions.cpp 

OBJS += \
./DDPModel/DDPDFNFunctions.o \
./DDPModel/DDPModel.o \
./DDPModel/DDPStudies.o \
./DDPModel/ExchangeFunctions.o 

CPP_DEPS += \
./DDPModel/DDPDFNFunctions.d \
./DDPModel/DDPModel.d \
./DDPModel/DDPStudies.d \
./DDPModel/ExchangeFunctions.d 


# Each subdirectory must supply rules for building sources it contributes
DDPModel/%.o: ../DDPModel/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/Developer/Library/Platform/CGAL-3.2.1/include -I/Developer/Library/Platform/CGAL-3.2.1/include/CGAL/config/msvc7 -I/Developer/Library/Platform/boost_1_49_0 -I/Developer/Library/Platform/boost-numeric-bindings -I/Developer/Library/Platform/UMFPACK/Include -I/System/Library/Frameworks/GLUT.framework/Versions/A/Headers -I/System/Library/Frameworks/OpenGL.framework/Versions/A/Headers -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


