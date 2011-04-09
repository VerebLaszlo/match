
# Add inputs and outputs from these tool invocations to the build variables 
#O_SRCS += \
../complex.o \
../complex_Test.o \
../windows.o 

C_SRCS += \
../detector.c \
../parameters.c \
../io_handler.c \
../generator.c \
../main_Test.c \
../main_New.c \
../main_Time.c \
../main_SQT-ST.c \
../match.c \
../match_qmss.c \
../util.c \
../util_math.c \
../util_math_tensor.c \
../test.c 

OBJS += \
./parameters.o \
./io_handler.o \
./detector.o \
./generator.o \
./match.o \
./match_qmss.o \
./util.o \
./util_math.o \
./util_math_tensor.o \
./test.o

C_DEPS += \
./parameters.d \
./io_handler.d \
./detector.d \
./generator.d \
./main_Test.d \
./main_New.d \
./main_Time.d \
./main_SQT-ST.d \
./match.d \
./match_qmss.d \
./util.d \
./util_math.d \
./util_math_tensor.d \
./test.d

# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -std=gnu99 -I/usr/include/ $(INCLUDES) -O0 -g3 -Wall -Wextra -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


