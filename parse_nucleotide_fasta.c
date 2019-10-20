#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <stdbool.h>

// typedef enum { false, true } bool;

typedef struct Flags {  // flag 0: file, flag 1: gc, flag 2: out file, flag 3: match, flag 4: merge, flag 5: seq.
    bool flag_raised;
    char flag_value[100];
} Flags;

typedef struct GC_data {
    int site_values[5];  // 0 = A; 1 = T; 2 = C, 3 = G, 4 = total
    int chunk_num;
    int chunk_size;
} GC_data;

typedef struct Reads {
    char prev_read[600];
    char curr_read[600];
} Reads;

typedef struct Files {
    FILE* fasta_file;
    FILE* merge_file;
    FILE* out_file;
} Files;


/******************************
*   Changes input into uppercase. 
*******************************/
void change_lowercase(char* text) {
    int i;
    for (i = 0; i < strlen(text); i++) {
        if (text[i] <= 'z' && text[i] >= 'a') {
            text[i] = text[i] - 32;
        }
    }
}

/******************************
*   Parses the input of the -seq option and splits the input by ,. 
*       If the sequence text passed in contains the same 
*******************************/
bool does_header_match(char* seq_search_string, char* curr_read, bool* prev_header) {
    if (curr_read[0] != '>') {
        return false;
    } else {
        char search_cpy[600];
        strcpy(search_cpy, seq_search_string);
        char* rest = search_cpy;
        char* token;
        while ((token = strtok_r(rest, ",", &rest))) {
            char cmp_string[600] = ">";
            if (!strcmp(curr_read, strcat(cmp_string, token))) {
                return true;
            }
        }
        return false;
    }
}
/******************************
*   Ensures input is workable and parses the input for the settings the program should use.
*******************************/
int ensure_legal_arguments(int argcount, char** argvalues, struct Flags* flags) {
    if (argcount == 1) {
        printf("No arguments entered, please see maunal or enter ./parse_nucleotide_fasta.out -help\n");
        return -1;
    }

    static struct option long_options[] = {
        {"file", required_argument, NULL, 'f'},
        {"gc", optional_argument, NULL, 'c'},
        {"match", required_argument, NULL, 'm'},
        {"merge", required_argument, NULL, 'g'},
        {"seq", required_argument, NULL, 's'},
        {"out", required_argument, NULL, 'o'},
        {"help", no_argument, NULL, 'h'},
        {NULL, 0, NULL, 0}
    };

    int opts;
    int argument_legality = -1;
    while ((opts = getopt_long_only(argcount, argvalues, "f:c::m:g:hs:o:", long_options, NULL)) != -1) {
        switch (opts) {
            case 'f':
                argument_legality = 0;
                flags[0].flag_raised = true;
                strcpy(flags[0].flag_value, optarg);
                printf("File %s passed in.\n", optarg);
                break;
            case 'c':
                argument_legality = 0;
                flags[1].flag_raised = true;
                if (optarg == NULL) {
                    printf("Calculating global GC content.\n");
                } else {
                    strcpy(flags[1].flag_value, optarg);
                    printf("Calculating GC content per %s nucleotides\n", optarg);
                }
                break;
            case 'o':
                argument_legality = 0;
                flags[2].flag_raised = true;
                strcpy(flags[2].flag_value, optarg);
                printf("Printing output into file %s.\n", flags[2].flag_value);
                break;
            case 'm':
                argument_legality = 0;
                flags[3].flag_raised = true;
                strcpy(flags[3].flag_value, optarg);
                change_lowercase(flags[3].flag_value);
                printf("Searching for all instances of %s in file.\n", flags[3].flag_value);
                break;
            case 'h':
                printf("Parse Nucleotide Fasta Options: \n");
                printf("\t-help : displays this message.\n");
                printf("\t-file <file> : directs the program to the fasta file.\n");
                printf("\t-gc : global gc count.\n");
                printf("\t\t-gc=<x> : gc count per x bases.\n");
                printf("\t-match <string> : searches for all instances of a string. (only words for same case searches)\n");
                printf("\t-merge <merged file name> : merges sequences into one and outputs result into file.\n");
                printf("\t-seq <options> : options must be separated by , or -\n");
                printf("\t\t-seq <seq_numbers> : performes operations on only the sequence numbers given.\n");
                printf("\t\t-seq <name> : performes operations on the named sequences.\n");
                printf("\t-out <filename> : outputs data to file.\n\n");
                argument_legality = -1;
                break;
            case 'g':
                argument_legality = 0;
                flags[4].flag_raised = true;
                strcpy(flags[4].flag_value, optarg);
                printf("Merging sequences into file: %s.\n", flags[4].flag_value);
                break;
            case 's':
                argument_legality = 0;
                flags[5].flag_raised = true;
                strcpy(flags[5].flag_value, optarg);
                printf("Carrying out operations only on selected sequences using sequence identifier string: %s.\n", flags[5].flag_value);
                break;
            case '?':
                printf("Unknown or incorrect argument entered, exiting program. Please see maunal or enter ./parse_nucleotide_fasta.out -help\n");
                argument_legality = -1;
                break;
            default:
                printf("No legal options entered. Please see maunal or enter ./parse_nucleotide_fasta.out -help\n");
                argument_legality = -1;
                break;
        }
    }
    if (flags[0].flag_raised == false) {
        printf("A fasta file must be passed in. Please see maunal or enter ./parse_nucleotide_fasta.out -help\n");
        argument_legality = -1;
    } else if (flags[1].flag_raised == false &&
               flags[2].flag_raised == false &&
               flags[3].flag_raised == false &&
               flags[4].flag_raised == false &&
               flags[5].flag_raised == false) {
        printf("No process arguments entered, select a process for the program carry out. Please see maunal or enter ./parse_nucleotide_fasta.out -help\n");
        argument_legality = -1;
    }
    return argument_legality;
}

/******************************
*   Prints output to either command line or to the file passed in from -out <file>.
*******************************/
void print_output(struct Flags* flags, struct GC_data* gc_data, int occurence_count, struct Files* files, bool called_from_gc_func) {
    if (flags[2].flag_raised != false) {
        if (called_from_gc_func != false) {
            fprintf(files->out_file, "Chunk number: %d A: %d T: %d C: %d G:%d\n", gc_data->chunk_num, gc_data->site_values[0], gc_data->site_values[1], gc_data->site_values[2], gc_data->site_values[3]);
            fprintf(files->out_file, "Chunk: %d GC content: %f\n\n", gc_data->chunk_num, (((float)gc_data->site_values[2] + (float)gc_data->site_values[3]) / (float)gc_data->site_values[4]) * 100);
        } else {
            if (flags[1].flag_raised != false) {
                fprintf(files->out_file, "Total size of final chunk: %d. Chunk num: %d A: %d T: %d C: %d G:%d\n", gc_data->site_values[4], gc_data->chunk_num, gc_data->site_values[0], gc_data->site_values[1], gc_data->site_values[2], gc_data->site_values[3]);
                fprintf(files->out_file, "Chunk: %d GC content: %f\n\n", gc_data->chunk_num, (((float)gc_data->site_values[2] + (float)gc_data->site_values[3]) / (float)gc_data->site_values[4]) * 100);
            }
            if (flags[3].flag_raised != false) { // need different file output here
                fprintf(files->out_file, "Total number of occurences of %s in the file was %d\n\n", flags[3].flag_value, occurence_count);
            }
        }
    } else {
        if (called_from_gc_func != false) {
            printf("Chunk number: %d A: %d T: %d C: %d G:%d\n", gc_data->chunk_num, gc_data->site_values[0], gc_data->site_values[1], gc_data->site_values[2], gc_data->site_values[3]);
            printf("Chunk: %d GC content: %f\n\n", gc_data->chunk_num, (((float)gc_data->site_values[2] + (float)gc_data->site_values[3]) / (float)gc_data->site_values[4]) * 100);
        } else {
            if (flags[1].flag_raised != false) {
                printf("Total size of final chunk: %d. Chunk num: %d A: %d T: %d C: %d G:%d\n", gc_data->site_values[4], gc_data->chunk_num, gc_data->site_values[0], gc_data->site_values[1], gc_data->site_values[2], gc_data->site_values[3]);
                printf("Chunk: %d GC content: %f\n\n", gc_data->chunk_num, (((float)gc_data->site_values[2] + (float)gc_data->site_values[3]) / (float)gc_data->site_values[4]) * 100);
            }
            if (flags[3].flag_raised != false) { // need different file output here
                printf("Total number of occurences of %s was %d\n\n", flags[3].flag_value, occurence_count);
            }
        }
    }
}

/******************************
*   Counts the GC content in the read passed into the function. 
*******************************/
void GC_count(struct Flags* flags, char* read, struct GC_data* gc_data, struct Files* files) {
    int site_index = 0;
    while (read[site_index] != '\0') {
        if(gc_data->site_values[4] == gc_data->chunk_size && gc_data->chunk_size != 0) {
            print_output(flags, gc_data, 0, files, true);
            int i;
            for (i = 0; i < 5; i++) {
                gc_data->site_values[i] = 0;
            }
            gc_data->chunk_num++;
        }
        char site = read[site_index];
        if (site == 'A') {
            gc_data->site_values[0]++;
            gc_data->site_values[4]++;
        } else if (site == 'T') {
            gc_data->site_values[1]++;
            gc_data->site_values[4]++;
        } else if (site == 'C') {
            gc_data->site_values[2]++;
            gc_data->site_values[4]++;
        } else if (site == 'G') {
            gc_data->site_values[3]++;
            gc_data->site_values[4]++;
        } else {
          // Do Nothing
        }
        site_index++;
    }
}

/******************************
*   Searches recursivly the bases in two reads (the current read that the program is on and the prevous read).
*   It will output 1 if the search gets to the first index in the word without matching errors (meaning the
*       read and the word match) and 0 otherwise.
*******************************/
int recursive_char_search(struct Flags* flags, int word_index, char* curr_read, char* prev_read, int index) {
    if (word_index < 0) {
        return 1;
    }
    if (index < 0) {
        int prev_read_size = strlen(prev_read);
        if (prev_read[prev_read_size + index] == flags[3].flag_value[word_index]) {
            return recursive_char_search(flags, word_index - 1, curr_read, prev_read, index - 1);
        } else {
            return 0;
        }
    } else {
        if (curr_read[index] == flags[3].flag_value[word_index]) {
            return recursive_char_search(flags, word_index - 1, curr_read, prev_read, index - 1);
        } else {
            return 0;
        }
    }
}

/******************************
*   Counts the number of times a word is found in the fasta file
*******************************/
void matching_occurences_count(struct Flags* flags, char* curr_read, char* prev_read, int* occurence_count) {
    int curr_read_size = strlen(curr_read);
    int word_size = strlen(flags[3].flag_value);
    int i;
    if (strcmp(prev_read, "")) {
        for (i = 0; i < curr_read_size; i++) {
            *occurence_count += recursive_char_search(flags, word_size - 1, curr_read, prev_read, i);
        }
    } else {
        for (i = word_size - 1; i < curr_read_size; i++) {
            *occurence_count += recursive_char_search(flags, word_size - 1, curr_read, prev_read, i);
        }
    }
}

/******************************
*   Iterates over the lines of the opened fasta file.
*******************************/
void iterate_over_lines(struct Files* files, struct Flags* flags, struct Reads* reads, struct GC_data* gc_data, int* occurence_count) {
    while (!feof(files->fasta_file)) {
        fscanf(files->fasta_file, "%s\n", reads->curr_read);
        bool prev_header = false;
        if (flags[5].flag_raised == false ||
            does_header_match(flags[5].flag_value, reads->curr_read, &prev_header) != false ||
            prev_header != false) {
            if (reads->curr_read[0] == '>') {
                strcpy(reads->prev_read, "");
                fscanf(files->fasta_file, "%s\n", reads->curr_read);
            }
            if (files->merge_file != NULL) {
                fprintf(files->merge_file, "%s\n", reads->curr_read);
            }
            change_lowercase(reads->curr_read);
            if (flags[1].flag_raised != false) {  // gc
                GC_count(flags, reads->curr_read, gc_data, files);
            }
            if (flags[3].flag_raised != false) {  // match
                matching_occurences_count(flags, reads->curr_read, reads->prev_read, occurence_count);
                strcpy(reads->prev_read, reads->curr_read);
            }
        }
    }
}

/******************************
*   Reads through fasta file and runs other functions depending on the input flags passed in.
*       This function opens the file and reads through each line, It will pass the lines to other functions which do stuff
*******************************/
void parse_fasta(struct Flags* flags) {
    struct Reads reads;
    strcpy(reads.curr_read, "");
    strcpy(reads.prev_read, "");

    struct GC_data gc_data;
    gc_data.chunk_num = 1;
    gc_data.chunk_size = strtol(flags[1].flag_value, NULL, 10);
    int i;
    for (i = 0; i < 5; i++) {
        gc_data.site_values[i] = 0;
    }

    int occurence_count = 0;
    struct Files files;
    files.fasta_file = fopen(flags[0].flag_value, "r");
    if (files.fasta_file == NULL) {
        printf("The fasta file %s does not appear to exist. Exiting program.\n", flags[0].flag_value);
        return;
    }

    if (flags[2].flag_raised != false) {
        files.out_file = fopen(flags[2].flag_value, "w");
        if (flags[4].flag_raised != false) {
            files.merge_file = fopen(flags[4].flag_value, "w");
            fprintf(files.merge_file, ">Sequences_Merged\n");
            iterate_over_lines(&files, flags, &reads, &gc_data, &occurence_count);
            print_output(flags, &gc_data, occurence_count, &files, false);
            fclose(files.merge_file);
        } else {
            iterate_over_lines(&files, flags, &reads, &gc_data, &occurence_count);
            print_output(flags, &gc_data, occurence_count, &files, false);
        }
        fclose(files.out_file);
    } else {
        if (flags[4].flag_raised != false) {
            files.merge_file = fopen(flags[4].flag_value, "w");
            fprintf(files.merge_file, ">Sequences_Merged\n");
            iterate_over_lines(&files, flags, &reads, &gc_data, &occurence_count);
            print_output(flags, &gc_data, occurence_count, &files, false);
            fclose(files.merge_file);
        } else {
            iterate_over_lines(&files, flags, &reads, &gc_data, &occurence_count);
            print_output(flags, &gc_data, occurence_count, &files, false);
        }
    }
    fclose(files.fasta_file);
}


int main(int argc, char** argv) {
    int number_of_flags = 6;
    struct Flags flags[number_of_flags];
    int error_flag = 0;

    int i;
    for (i = 0; i < number_of_flags; i++) {
        flags[i].flag_raised = false;
        strcpy(flags[i].flag_value, "");
    }

    error_flag = ensure_legal_arguments(argc, argv, flags);
    if (error_flag == -1) {
        return 0;
    }

    parse_fasta(flags);

    return 0;
}
