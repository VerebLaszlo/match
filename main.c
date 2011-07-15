/**
 * @file main.c
 *
 * @date Apr 8, 2011
 * @author vereb
 */
#include "match_qmss.h"

extern char apr[2][FILENAME_MAX];

int switchMode = LALSQTPN_PRECESSING;

static char*help = "Usage:\n"
	"-mode [options] program_file parameter_file\n"
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
	case 'e':
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
	read_Program_Parameters(file, &program_Parameters);
	fclose(file);
	strcpy(apr[0], parameters.approx[0]);
	strcpy(apr[1], parameters.approx[1]);
	switch (option) {
	case 'm':
	case 'e':
		break;
	case 'g':
	case 'd':
	case 't':
	case 's':
		file = safely_Open_File_For_Reading(argv[arg++]);
		read_System_Parameters(file, &parameters);
		fclose(file);
		//print_System_Parameters(stdout, &parameters);
		print_Program_Parameters(stdout, &program_Parameters);
		break;
	default:
		break;
	}
	switch (option) {
	case 'g': {
		binary_System limits[2];
		memmove(limits, parameters.system, 2 * sizeof(binary_System));
		memset(parameters.system, 0, 2 * sizeof(binary_System));
		srand(86);
		double percent;
		for (long i = 0; i < program_Parameters.number_Of_Runs; i++) {
			generate_Same_Parameters(&parameters, limits, ETAM);
			find_Waveform_Errors_At_Parameter(&program_Parameters, &parameters, i);
			//if (i == 0) {
			//	exit(EXIT_FAILURE);
			//}
			percent = 1000.0 * (double)i / (double)program_Parameters.number_Of_Runs;
			if (floor(percent) == ceil(percent)) {
				if (((int)percent + 1) % 10 == 0) {
					printf("%d%%o\n", (int)percent + 1);
					fflush(stdout);
				}
			}
		}
		break;
	}
	case 'd': {
		binary_System limits[2];
		memmove(limits, parameters.system, 2 * sizeof(binary_System));
		memset(parameters.system, 0, 2 * sizeof(binary_System));
		srand(86);
		for (int i = 0; i < program_Parameters.number_Of_Runs; i++) {
			generate_Same_Parameters(&parameters, limits, M1M2);
			print_System_Parameters(stdout, &parameters);
			signalStruct sig;
			puts("Starting difference ...");
			generate_Waveforms_For_Difference(&program_Parameters, &parameters, &sig);
			char file_Name[FILENAME_MAX];
			sprintf(file_Name, "%s/diff%d.txt", program_Parameters.folder,i);
			file = safely_Open_File_For_Writing(file_Name);
			//print_Two_Signals_With_HPHC(file, &sig, parameters.time_Sampling,
			//		program_Parameters.width_Of_Number_To_Plot, program_Parameters.precision_To_Plot);
			print_Two_Signals(file, &sig, parameters.time_Sampling,
					program_Parameters.width_Of_Number_To_Plot,
					program_Parameters.precision_To_Plot);
			destroy_Signal_Struct1(&sig);
		}
		fclose(file);
		break;
	}
	case 't':
		puts("Starting time ...");
		calc_Time(&program_Parameters, &parameters, sampling);
		break;
	case 'm': {
		FILE*fileIn = safely_Open_File_For_Reading(argv[arg++]);
		//int xx = 0;
		puts("Starting match ...");
		do {
			readExactParameters(fileIn, &parameters);
			//print_System_Parameters(stdout, &parameters);
			//if (xx == 0) {
			//parameters.system[0].incl = parameters.system[1].incl = 1.43;
			//}
			//print_System_Parameters(stdout, &parameters);
			convert_Masses(&parameters.system[0], FROM_M1M2);
			convert_Masses(&parameters.system[1], FROM_M1M2);
			convert_Spins(&parameters.system[0], FROM_KAPPA_PSI);
			convert_Spins(&parameters.system[1], FROM_KAPPA_PSI);
			printf("X: %d %d\n", parameters.amp_Code[0], parameters.amp_Code[1]);
			signalStruct sig;
			calc_Matches_For_ParameterPair(&program_Parameters, &parameters, &sig);
			char file_Name[FILENAME_MAX];
			sprintf(file_Name, "%s/%s.txt", program_Parameters.folder, parameters.name[0]);
			puts(file_Name);
			file = safely_Open_File_For_Writing(file_Name);
			print_System_Parameters_For_Plot(file, &parameters);
			print_Two_Signals(file, &sig, parameters.time_Sampling,
					program_Parameters.width_Of_Number_To_Plot,
					program_Parameters.precision_To_Plot);
			destroy_Signal_Struct(&sig);
			fclose(file);
			//if (xx++ == 0) {
			//	break;
			//}
		} while (!feof(fileIn));
		fclose(fileIn);
		break;
	}
	case 'e': {
		FILE*fileIn = safely_Open_File_For_Reading(argv[arg++]);
		puts("Starting exact ...");
		for (int i = 0; !feof(fileIn); i++) {
			readExactParameters(fileIn, &parameters);
			convert_Masses(&parameters.system[0], FROM_M1M2);
			convert_Masses(&parameters.system[1], FROM_M1M2);
			convert_Spins(&parameters.system[0], FROM_KAPPA_PSI);
			convert_Spins(&parameters.system[1], FROM_KAPPA_PSI);
			puts("X 2");
			signalStruct sig;
			generate_Waveforms_For_Difference(&program_Parameters, &parameters, &sig);
			puts("X 3");
			char file_Name[FILENAME_MAX];
			sprintf(file_Name, "%s/diff%d.txt", program_Parameters.folder,i);
			file = safely_Open_File_For_Writing(file_Name);
			print_Two_Signals(file, &sig, parameters.time_Sampling,
					program_Parameters.width_Of_Number_To_Plot,
					program_Parameters.precision_To_Plot);
			destroy_Signal_Struct1(&sig);
		}
		fclose(fileIn);
		break;
	}
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
