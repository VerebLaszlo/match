
-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: test temp

# Tool invocations
test: main_New.o $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C Linker'
	gcc -std=gnu99 -L/usr/lib/ $(SHARED) -o"test" main_New.o $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

temp: main_Test.o $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C Linker'
	gcc -std=gnu99 -L/usr/lib/ $(SHARED) -o"temp" main_Test.o $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

doxy: ../*.h
	@echo 'Generating documentation'
	doxygen ../Doxyfile
	@echo 'Finished generating documentation'

# Other Targets
clean:
	-$(RM) $(OBJS)$(C_DEPS)$(EXECUTABLES) match main_Test.o main_New.o
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets
