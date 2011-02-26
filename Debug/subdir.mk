
# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../complex.o \
../complex_Test.o \
../windows.o 

C_SRCS += \
../detector.c \
../generator.c \
../main_New.c \
../match.c \
../match_qmss.c \
../util.c \
../util_math.c 

OBJS += \
./detector.o \
./generator.o \
./main_New.o \
./match.o \
./match_qmss.o \
./util.o \
./util_math.o 

C_DEPS += \
./detector.d \
./generator.d \
./main_New.d \
./match.d \
./match_qmss.d \
./util.d \
./util_math.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -std=gnu99 -I/usr/include/ $(INCLUDES) -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


