################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Utilitaries/Geometry/FluxPoint2D.cpp \
../Utilitaries/Geometry/Point_Cgal.cpp \
../Utilitaries/Geometry/Segment.cpp \
../Utilitaries/Geometry/scale.cpp 

OBJS += \
./Utilitaries/Geometry/FluxPoint2D.o \
./Utilitaries/Geometry/Point_Cgal.o \
./Utilitaries/Geometry/Segment.o \
./Utilitaries/Geometry/scale.o 

CPP_DEPS += \
./Utilitaries/Geometry/FluxPoint2D.d \
./Utilitaries/Geometry/Point_Cgal.d \
./Utilitaries/Geometry/Segment.d \
./Utilitaries/Geometry/scale.d 


# Each subdirectory must supply rules for building sources it contributes
Utilitaries/Geometry/%.o: ../Utilitaries/Geometry/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/Developer/Library/Platform/CGAL-3.2.1/include -I/Developer/Library/Platform/CGAL-3.2.1/include/CGAL/config/msvc7 -I/Developer/Library/Platform/boost_1_49_0 -I/Developer/Library/Platform/boost-numeric-bindings -I/Developer/Library/Platform/UMFPACK/Include -I/System/Library/Frameworks/GLUT.framework/Versions/A/Headers -I/System/Library/Frameworks/OpenGL.framework/Versions/A/Headers -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


