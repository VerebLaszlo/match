/**
 * @file main.c
 *
 * @date Apr 8, 2011
 * @author vereb
 */
#include "match_qmss.h"

extern char apr[2][FILENAME_MAX];

static char*help = "Usage:\n"
	"-mode [options] parameter_file program_file\n"
	"mode option:\n"
	"g NOT IMPLEMENTED JET!!!!\n"
	"d NOT IMPLEMENTED JET!!!!\n"
	"t    sampling\n"
	"m NOT IMPLEMENTED JET!!!!\n"
	"s NOT IMPLEMENTED JET!!!!";

int main(int argc, char *argv[]) {
	if (argc < 2) {
		puts(help);
		exit(EXIT_FAILURE);
	}
	short arg = 0;
	argc--;
	if ((++argv)[arg][0] != '-') {
		printf("The %d used parameters were: ", argc);
		fflush(stdout);
		short i = 0;
		while (i < argc) {
			printf("%s, ", argv[i++]);
		}
		puts("\n");
		puts(help);
		exit(EXIT_FAILURE);
	}
	short sampling;
	char option = argv[arg++][1];
	switch (option) {
	case 'g':
		break;
	case 'd':
		break;
	case 't':
		sampling = atoi(argv[arg++]);
		break;
	case 'm':
		break;
	case 's':
		break;
	default:
		printf("The %d used parameters were: ", argc);
		fflush(stdout);
		short i = 0;
		while (i < argc) {
			printf("%s, ", argv[i++]);
		}
		puts("\n");
		puts(help);
		exit(EXIT_FAILURE);
		break;
	}
	if (argc - arg != 2) {
		printf("The %d used parameters were: ", argc);
		fflush(stdout);
		short i = 0;
		while (i < argc) {
			printf("%s, ", argv[i++]);
		}
		puts("\n");
		puts(help);
		exit(EXIT_FAILURE);
	}
	Program_Parameters program_Parameters;
	System_Parameters parameters;
	FILE *file = safely_Open_File_For_Reading(argv[arg++]);
	read_System_Parameters(file, &parameters);
	fclose(file);
	file = safely_Open_File_For_Reading(argv[arg]);
	read_Program_Parameters(file, &program_Parameters);
	fclose(file);
	strcpy(apr[0], parameters.approx[0]);
	strcpy(apr[1], parameters.approx[1]);
	switch (option) {
	case 'g':
	case 'd':
	case 't':
	case 'm':
	case 's':
		//print_System_Parameters(stdout, &parameters);
		print_Program_Parameters(stdout, &program_Parameters);
		break;
	default:
		break;
	}
	switch (option) {
	case 'g':
		puts("Starting generating ...");
		break;
	case 'd': {
		binary_System limits[2];
		memmove(limits, parameters.system, 2 * sizeof(binary_System));
		memset(parameters.system, 0, 2 * sizeof(binary_System));
		generate_Same_Parameters(&parameters, limits);
		print_System_Parameters(stdout, &parameters);
		signalStruct sig;
		puts("Starting difference ...");
		generate_Waveforms_For_Difference(&program_Parameters, &parameters, &sig);
		char file_Name[FILENAME_MAX];
		sprintf(file_Name, "%s/diff.txt", program_Parameters.folder);
		file = safely_Open_File_For_Writing(file_Name);
		print_Two_Signals_With_HPHC(file, &sig, parameters.time_Sampling, program_Parameters.width_Of_Number_To_Plot, program_Parameters.precision_To_Plot);
		fclose(file);
		break;
	}
	case 't':
		puts("Starting time ...");
		calc_Time(&program_Parameters, &parameters, sampling);
		break;
	case 'm':
		puts("Starting match ...");
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
