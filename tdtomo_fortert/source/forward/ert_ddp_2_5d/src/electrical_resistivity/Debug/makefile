################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include Utilitaries/Visualization/subdir.mk
-include Utilitaries/Storage/subdir.mk
-include Utilitaries/Math/subdir.mk
-include Utilitaries/Geometry/subdir.mk
-include RngStream/subdir.mk
-include InputOutput/subdir.mk
-include FourierModel/subdir.mk
-include EllipseComputation/subdir.mk
-include EPMModel/subdir.mk
-include DomainDefinition/subdir.mk
-include DFNModel/subdir.mk
-include DDPModel/subdir.mk
-include .metadata/.plugins/org.eclipse.cdt.make.core/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: ElectricalResistivity

# Tool invocations
ElectricalResistivity: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: MacOS X C++ Linker'
	g++ -L/Developer/Library/Platform/CGAL-3.2.1/lib/i386_Darwin-10.8_g++-4.2.1 -L/Developer/Library/Platform/UMFPACK/Lib -L/Developer/Library/Platform/CHOLMOD/Lib -L/Developer/Library/Platform/AMD/Lib -L/Developer/Library/Platform/COLAMD/Lib -L/Developer/Library/Platform/Atlas/lib/windows -framework GLUT -framework openGL -framework Cocoa -o "ElectricalResistivity" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(CC_DEPS)$(C++_DEPS)$(EXECUTABLES)$(OBJS)$(C_UPPER_DEPS)$(CXX_DEPS)$(CPP_DEPS)$(C_DEPS) ElectricalResistivity
	-@echo ' '

.PHONY: all clean dependents

-include ../makefile.targets
