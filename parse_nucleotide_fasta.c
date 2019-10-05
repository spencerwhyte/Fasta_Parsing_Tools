#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>

typedef struct Flags {
    int flag_raised;  // 1 if raised, 0 if not.
    char flag_value[100];
} flag;

// can you initialize structs in a function and have them work??
/* int initialize_flags() {
    struct Flags filename;
    struct Flags others.
}*/

/******************************
*   ensure input it workable
*******************************/
int ensure_legal_arguments(int argcount, char** argvalues, struct Flags* flags) {
    if (argcount == 1) {
        printf("No arguments entered, please see maunal or enter ./parse_nuc_fasta -help\n");
        return -1;
    }

    static struct option long_options[] = {
        {"file", required_argument, NULL, 'f'},
        {"gc", optional_argument, NULL, 'c'},
        {"k_mers", required_argument, NULL, 'k'},
        {"match", required_argument, NULL, 'm'},
        {"help", no_argument, NULL, 'h'},
        {NULL, 0, NULL, 0}
        /*
        * Other future arguments:
        *   -out <output_file>
        *   -v <verbose_output>
        *   -seq <sequence_number> // only searches for the following sequences in the file (marked by >)
        *       example: -seq 2 would only work with the data chunk after the second > but before the third >
        *       -seq 2,4 should be allowed to be entered
        */
    };

    int opts;
    int argument_legality = -1;
    while ((opts = getopt_long_only(argcount, argvalues, "f:c::k:m:h", long_options, NULL)) != -1) {
        switch (opts) {
            case 'f':
                argument_legality = 0;
                flags[0].flag_raised = 1;
                //flags[0].flag_value = optarg;
                strcpy(flags[0].flag_value, optarg);
                printf("File %s passed in.\n", optarg);
                break;
            case 'c':
                argument_legality = 0;
                flags[1].flag_raised = 1;
                printf("optarg %s\n", optarg);
                if (optarg == NULL) {
                    printf("Calculating global GC content.\n");
                } else {
                strcpy(flags[1].flag_value, optarg);
                printf("Calculating GC content per %s nucleotides\n", optarg);
                }
                break;
            case 'k':
                argument_legality = 0;
                flags[2].flag_raised = 1;
                strcpy(flags[2].flag_value, optarg);
                printf("Determining all k-mers of size %s.\n", optarg);
                break;
            case 'm':
                argument_legality = 0;
                flags[3].flag_raised = 1;
                strcpy(flags[3].flag_value, optarg);
                printf("Searching for all instances of %s in file.\n", optarg);
                break;
            case 'h':
                printf("Parse_nuc_fasta Options: \n");
                printf("\t-help : displays this message\n");
                printf("\t-file <file> : Directs the program to the fasta file.\n");
                printf("\t-gc : global gc count\n");
                printf("\t-gc=<x> : gc count per x bases\n");
                printf("\t-k_mers <k> : lists all k-mers of size k\n");
                printf("\t-match <string> : searches for all instances of a string\n");
                argument_legality = -1;
                break;
            case '?':
                printf("Unknown or incorrect argument entered, exiting program. Please see maunal or enter ./parse_nuc_fasta -help\n");
                argument_legality = -1;
                break;
            default:
                printf("No legal options entered. Please see maunal or enter ./parse_nuc_fasta -help\n");
                argument_legality = -1;
                break;
        }
    }
    if (flags[0].flag_raised == 0) {
        printf("A fasta file must be passed in. Please see maunal or enter ./parse_nuc_fasta -help\n");
        argument_legality = -1;
    } else if (flags[1].flag_raised == 0 &&
               flags[2].flag_raised == 0 &&
               flags[3].flag_raised == 0) {
        printf("No process arguments entered, select a process for the program carry out. Please see maunal or enter ./parse_nuc_fasta -help\n");
        argument_legality = -1;
    }
    return argument_legality;
}

/******************************
*   read through fasta file
*******************************/
int parse_fasta(struct Flags* flags) {
    // this program opens the file and reads through each line, It will pass the lines to other functions which do stuff
    char prev_read[600];
    char curr_read[600];
    sequence_number = 0;
    FILE* file = fopen(flags[0].flag_value, "r");
    while (!feof(file)) {
        fscanf(file, "%s\n", read);
        // remember to skip the sequence identifier lines

        if (flags[])
    }

    fclose(file);
    return 0;
}
int skip_info_lines() {
    // skips the > lines in a fasta file
}

int GC_count() {
    // counts GC content per X bases
        site_index = 0;
    while (read[site_index] != '\0') {
        if(count == chunk_val) {
           printf("Chunk num: %d A: %d T: %d C: %d G:%d\n", chunk_num, count_a, count_t, count_c, count_g);
           printf("Chunk: %d GC content: %f\n", chunk_num, (((float)count_c + (float)count_g) / (float)count) * 100);
           count = 0;
           count_a = 0;
           count_t = 0;
           count_c = 0;
           count_g = 0;
           chunk_num++;
        }
        char site = read[site_index];
        if (site == 'A') {
            count_a++;
            count++;
        } else if (site == 'T') {
            count_t++;
            count++;
        } else if (site == 'C') {
            count_c++;
            count++;
        } else if (site == 'G') {
            count_g++;
            count++;
        } else {
          // Do Nothing
        }
        site_index++;
    }
}
    printf("Total size of final chunk: %d. Chunk num: %d A: %d T: %d C: %d G:%d\n", count, chunk_num, count_a, count_t, count_c, count_g);
    printf("Chunk: %d GC content: %f\n", chunk_num, (((float)count_c + (float)count_g) / (float)count) * 100);

    fclose(file);

}

int word_occurences_count() {
    // Counts the number of times a word is found in the fasta file
}

int find_all_kmer_permutations() {
    // Finds all the permutations for all K-mers of size K
}


int main(int argc, char** argv) {
    int number_of_flags = 4;
    struct Flags flags[number_of_flags];
    int error_flag = 0;

    int i;
    for (i = 0; i < number_of_flags; i++) {
        flags[i].flag_raised = 0;  // flag 0: file, flag 1: gc, flag 2: kmer, flag 3: match.
        strcpy(flags[i].flag_value, "");
    }

    error_flag = ensure_legal_arguments(argc, argv, flags);
    if (error_flag == -1) {
        return 0;
    }

    /* Check flags and flag values
    for (i = 0; i < number_of_flags; i++) {
        printf("is flag %d raised: %d\n", i, flags[i].flag_raised);  // flag 0: file, flag 1: gc , flag 2: kmer, flag 3: match
        printf("flag %i value: %s\n", i, flags[i].flag_value);
    }
    */

    read_through_file();


    long int chunk_val = strtol(chunk_size, NULL, 10);
    char read[600];
    int count = 0;
    FILE* file = fopen(fileName, "r");
    int site_index;
    int count_a = 0, count_t = 0, count_c = 0, count_g = 0, chunk_num = 1;
while (!feof(file)) {
    fscanf(file, "%s\n", read);
    site_index = 0;
    while (read[site_index] != '\0') {
        if(count == chunk_val) {
           printf("Chunk num: %d A: %d T: %d C: %d G:%d\n", chunk_num, count_a, count_t, count_c, count_g);
           printf("Chunk: %d GC content: %f\n", chunk_num, (((float)count_c + (float)count_g) / (float)count) * 100);
           count = 0;
           count_a = 0;
           count_t = 0;
           count_c = 0;
           count_g = 0;
           chunk_num++;
        }
        char site = read[site_index];
        if (site == 'A') {
            count_a++;
            count++;
        } else if (site == 'T') {
            count_t++;
            count++;
        } else if (site == 'C') {
            count_c++;
            count++;
        } else if (site == 'G') {
            count_g++;
            count++;
        } else {
          // Do Nothing
        }
        site_index++;
    }
}
    printf("Total size of final chunk: %d. Chunk num: %d A: %d T: %d C: %d G:%d\n", count, chunk_num, count_a, count_t, count_c, count_g);
    printf("Chunk: %d GC content: %f\n", chunk_num, (((float)count_c + (float)count_g) / (float)count) * 100);

    fclose(file);

    return 0;
}
