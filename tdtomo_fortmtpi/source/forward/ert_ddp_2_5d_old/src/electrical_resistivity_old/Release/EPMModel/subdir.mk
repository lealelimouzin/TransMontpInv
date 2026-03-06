################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../EPMModel/EPMComputation.cpp \
../EPMModel/EPMDefinition.cpp 

OBJS += \
./EPMModel/EPMComputation.o \
./EPMModel/EPMDefinition.o 

CPP_DEPS += \
./EPMModel/EPMComputation.d \
./EPMModel/EPMDefinition.d 


# Each subdirectory must supply rules for building sources it contributes
EPMModel/%.o: ../EPMModel/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXX) $(CXXFLAGS) -I ../../Library/ -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


