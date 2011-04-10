/**
 * @file main.c
 *
 * @date Apr 8, 2011
 * @author vereb
 */
#include "match_qmss.h"

int main(int argc, char *argv[]) {
	if (argc < 4) {
		puts("Usage:\nprogram_file parameter_file -mode [options]\n");
		puts("modes:\n"
			"t for time measurement, with the samling interval as option, ex: -t 100\n"
			"d NOT IMPLEMENTED JET!!!!"
			"m NOT IMPLEMENTED JET!!!!"
			"g NOT IMPLEMENTED JET!!!!");
		printf("The %d used parameters were: ", argc - 1);
		fflush(stdout);
		while (argc-- > 1) {
			printf("%s, ", argv[argc]);
		}
		exit(EXIT_FAILURE);
	}
	Program_Parameters program_Parameters;
	System_Parameters parameters;
	short arg = 1;
	puts("Reading parameters ...");
	FILE*file = safely_Open_File_For_Reading(argv[arg++]);
	read_Program_Parameters(file, &program_Parameters);
	fclose(file);
	file = safely_Open_File_For_Reading(argv[arg++]);
	read_System_Parameters(file, &parameters);
	fclose(file);
	print_Program_Parameters(stdout, &program_Parameters);
	print_System_Parameters(stdout, &parameters);
	switch (argv[arg++][1]) {
	case 'd':
		puts("Starting difference ...");
		break;
	case 'm':
		puts("Starting match ...");
		break;
	case 't':
		if (argc != 5) {
			puts("Not enough parameter!!!");
			exit(EXIT_FAILURE);
		}
		puts("Starting time ...");
		calc_Time(&program_Parameters, &parameters, atoi(argv[arg++]));
		break;
	case 'g':
		puts("Starting generating ...");
		break;
	case 's':
		puts("Starting search for spin greter than 1 ...");
		find_Spin_Greater_Than1(&program_Parameters, &parameters);
	default:
		break;
	}
	LALCheckMemoryLeaks();
	puts("Done!");
	exit(EXIT_SUCCESS);
}
