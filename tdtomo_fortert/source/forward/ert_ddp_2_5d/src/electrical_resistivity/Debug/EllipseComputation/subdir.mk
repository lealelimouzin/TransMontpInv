################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../EllipseComputation/EllipseComputation.cpp 

OBJS += \
./EllipseComputation/EllipseComputation.o 

CPP_DEPS += \
./EllipseComputation/EllipseComputation.d 


# Each subdirectory must supply rules for building sources it contributes
EllipseComputation/%.o: ../EllipseComputation/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/Developer/Library/Platform/CGAL-3.2.1/include -I/Developer/Library/Platform/CGAL-3.2.1/include/CGAL/config/msvc7 -I/Developer/Library/Platform/boost_1_49_0 -I/Developer/Library/Platform/boost-numeric-bindings -I/Developer/Library/Platform/UMFPACK/Include -I/System/Library/Frameworks/GLUT.framework/Versions/A/Headers -I/System/Library/Frameworks/OpenGL.framework/Versions/A/Headers -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


