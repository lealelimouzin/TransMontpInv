################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Utilitaries/Math/Laplace_Stehfest.cpp \
../Utilitaries/Math/LinearSystem.cpp \
../Utilitaries/Math/Operations.cpp \
../Utilitaries/Math/RandomNumber.cpp 

OBJS += \
./Utilitaries/Math/Laplace_Stehfest.o \
./Utilitaries/Math/LinearSystem.o \
./Utilitaries/Math/Operations.o \
./Utilitaries/Math/RandomNumber.o 

CPP_DEPS += \
./Utilitaries/Math/Laplace_Stehfest.d \
./Utilitaries/Math/LinearSystem.d \
./Utilitaries/Math/Operations.d \
./Utilitaries/Math/RandomNumber.d 


# Each subdirectory must supply rules for building sources it contributes
Utilitaries/Math/%.o: ../Utilitaries/Math/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXX) $(CXXFLAGS) -I ../../Library/ -I /usr/include/suitesparse/ -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


