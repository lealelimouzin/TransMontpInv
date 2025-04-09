################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../DDPModel/DDPDFNFunctions.cpp \
../DDPModel/DDPModel.cpp \
../DDPModel/DDPStudies.cpp \
../DDPModel/ExchangeFunctions.cpp \
../DDPModel/SystLinearFunctions.cpp

OBJS += \
./DDPModel/DDPDFNFunctions.o \
./DDPModel/DDPModel.o \
./DDPModel/DDPStudies.o \
./DDPModel/ExchangeFunctions.o \
./DDPModel/SystLinearFunctions.o

CPP_DEPS += \
./DDPModel/DDPDFNFunctions.d \
./DDPModel/DDPModel.d \
./DDPModel/DDPStudies.d \
./DDPModel/ExchangeFunctions.d \
./DDPModel/SystLinearFunctions.d


# Each subdirectory must supply rules for building sources it contributes
DDPModel/%.o: ../DDPModel/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXX) $(CXXFLAGS) -I ../../Library/ -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


