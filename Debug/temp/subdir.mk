################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../temp/confuse-parser.c \
../temp/detector.c \
../temp/generator.c \
../temp/main_Generator.c \
../temp/main_Match.c \
../temp/main_Old_Match.c \
../temp/main_Spin_Match.c \
../temp/main_Test.c \
../temp/match.c \
../temp/match_Multi.c \
../temp/util.c \
../temp/util_math.c 

OBJS += \
./temp/confuse-parser.o \
./temp/detector.o \
./temp/generator.o \
./temp/main_Generator.o \
./temp/main_Match.o \
./temp/main_Old_Match.o \
./temp/main_Spin_Match.o \
./temp/main_Test.o \
./temp/match.o \
./temp/match_Multi.o \
./temp/util.o \
./temp/util_math.o 

C_DEPS += \
./temp/confuse-parser.d \
./temp/detector.d \
./temp/generator.d \
./temp/main_Generator.d \
./temp/main_Match.d \
./temp/main_Old_Match.d \
./temp/main_Spin_Match.d \
./temp/main_Test.d \
./temp/match.d \
./temp/match_Multi.d \
./temp/util.d \
./temp/util_math.d 


# Each subdirectory must supply rules for building sources it contributes
temp/%.o: ../temp/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


