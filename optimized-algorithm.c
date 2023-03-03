#include <stdlib.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define CIRCULAR_BUFFER_SIZE 300
#define BUFFER_LOADING_MIN_THRESHOLD 128
#define SEQ_SEQUENCE_MEMORY_CHUNK 4096
#define SEQ_MICROSATELLITE_MEMORY_CHUNK 64

typedef struct {
    char *array;
    size_t used;
    size_t size;
} sequence;

void initSequence(sequence *a, size_t initialSize){
    // Structure memory allocation will be outside
    a->array = malloc(initialSize * sizeof(char) +1);
    a->used = 0;
    a->size = initialSize;
    a->array[a->used+1] == 0;
}

void insertSequence(sequence *a, char element){
    if (a->used+1 == a->size){
        a->size += SEQ_SEQUENCE_MEMORY_CHUNK;
        a->array = realloc(a->array, a->size * sizeof(char));
    }
    a->array[a->used++] = element;
}

void freeSequence(sequence *a){
    free(a->array);
    a->array = NULL;
    a->used = a->size = 0;
}

typedef struct {
    char *name;
    char *description;
    sequence *sequence;
} record;

typedef struct {
    char *sequence;
    char *motif;
    int period;
    int repeat;
    int start;
    int end;
    int length;
} microsatellite;

typedef struct {
    microsatellite *array;
    size_t used;
    size_t size;
} microsatelliteArray;

void initMicrosatelliteArray(microsatelliteArray *a, size_t initialSize){
    a->array = malloc(initialSize * sizeof(microsatellite));
    a->used = 0;
    a->size = initialSize;
}

void insertMicrosatelliteArray(microsatelliteArray *a, microsatellite *element){
    if (a->used == a->size){
        a->size += SEQ_MICROSATELLITE_MEMORY_CHUNK;
        a->array = realloc(a->array, a->size*sizeof(microsatellite));
    }
    a->array[a->used++] = *element;
}

void readConfig(int *output, FILE *f){
    char x[8];
    int cursor=0;
    while (fscanf(f, "%1023s", x) == 1){
        if(cursor%2==1)
            output[cursor/2]=atoi(x);
        cursor++;
    }
    int a = 0;
}


void search_perfect_microsatellites(microsatelliteArray *output, FILE *fptr, int *minRepeats, int *n_length,
                                    int *gc_length, int *total_length) {

    char sequenceName[64];
    int start;
    int length;
    int repeat;
    int i;
    int j;
    char motif[7] = "\0";


    char _buffer[4096]; // Buffer used to store temporary line
    int _buffer_offset=0;
    for (i=0; !feof(fptr); i++){
        if(_buffer[i-_buffer_offset] == 0){
            // Get length of the previous buffer
            int leap = strlen(_buffer);
            fgets(_buffer, sizeof(_buffer),  fptr);
            strtok(_buffer, "\n");

            if(_buffer[0] == '>'){
                sscanf(_buffer, "%63s", sequenceName);
                // empty buffer
                i--;
                _buffer[0] = (char)0;
                _buffer_offset=0;
                continue;
            }else{
                _buffer_offset+= leap;
            }
        }
        if(_buffer[i-_buffer_offset] == 'N')
            *n_length+=1;
        else if(_buffer[i-_buffer_offset] == 'C' || _buffer[i-_buffer_offset] == 'G')
            *gc_length+=1;

        *total_length+=1;
    }

    int buffer_offset = 0; // I think, no use
    char buffer[CIRCULAR_BUFFER_SIZE];
    char fgets_tmp[256];
    memset(buffer, 0, CIRCULAR_BUFFER_SIZE);
    rewind(fptr);
    int buffer_loading_cursor = 0;
    for (int i=0; !feof(fptr); i++){
        if(buffer[(i+BUFFER_LOADING_MIN_THRESHOLD) % CIRCULAR_BUFFER_SIZE] == 0){
            memset(fgets_tmp, 0, 256);
            fgets(fgets_tmp, sizeof(fgets_tmp), fptr),
            strtok(fgets_tmp, "\n");
            if(fgets_tmp[0] == '>'){
                i--;
                continue;
            }else{
                strncpy(buffer+(buffer_loading_cursor%CIRCULAR_BUFFER_SIZE), fgets_tmp, strlen(fgets_tmp));
                memset(buffer+(buffer_loading_cursor+strlen(fgets_tmp)%CIRCULAR_BUFFER_SIZE), 0, 128);
                if((buffer_loading_cursor%CIRCULAR_BUFFER_SIZE) > ((buffer_loading_cursor+strlen(fgets_tmp))%CIRCULAR_BUFFER_SIZE)){
                    int exceeding =buffer_loading_cursor+(int)strlen(fgets_tmp) - CIRCULAR_BUFFER_SIZE;
                    int loaded = (int)strlen(fgets_tmp) - exceeding;
                    strncpy(buffer, fgets_tmp+loaded, exceeding);
                    memset(buffer+exceeding, 0, 128);
                    printf("");
                }
                buffer_loading_cursor+=strlen(fgets_tmp);
                printf("Hello");
            }
        }

        // printf("%c", buffer[i%CIRCULAR_BUFFER_SIZE]);

        if(buffer[i%CIRCULAR_BUFFER_SIZE] == 78)
            continue;
        if(buffer[i%CIRCULAR_BUFFER_SIZE] != 'A' && buffer[i%CIRCULAR_BUFFER_SIZE] != 'C' && buffer[i%CIRCULAR_BUFFER_SIZE] != 'G' && buffer[i%CIRCULAR_BUFFER_SIZE] != 'T')
            continue;
        for (j=1; j<=6; j++) {
            start = i;
            length = j;
            while (buffer[i%CIRCULAR_BUFFER_SIZE] == buffer[(i+j)%CIRCULAR_BUFFER_SIZE] && buffer[i%CIRCULAR_BUFFER_SIZE] != 78) {
                i++;
                length++;
            }
            repeat=length/j;
            if(repeat >= minRepeats[j-1]){
                microsatellite *m = malloc(sizeof(microsatellite));
                m->motif = malloc(16*sizeof(char));
                strncpy(m->motif, buffer+(start%CIRCULAR_BUFFER_SIZE), j);
                m->motif[j] = 0;
                m->sequence=malloc(64*sizeof(char));
                strcpy(m->sequence, sequenceName+1);

                m->period=j;
                m->repeat=repeat;
                m->start=start+1;
                m->end=start+length;
                m->length=repeat*j;
                insertMicrosatelliteArray(output, m);
                free(m);

                i=start+length;
                j=0;
            }else{
                i=start;
            }
        }
    }
}

int main(int argc, char **argv){
    // Parse input
    const char *infile = NULL;
    const char *outfile = NULL;
    const char *configfile = NULL;
    int i;
    for(i = 1; i < argc; i++){
        if(strcmp(argv[i-1], "-i") == 0)
            infile = argv[i];
        if(strcmp(argv[i-1], "-o") == 0)
            outfile = argv[i];
        if(strcmp(argv[i-1], "-c") == 0)
            configfile = argv[i];
    }

    // Open the configuration file
    int *minRepeats;
    if(configfile != NULL){
        FILE *fptr;
        fptr = fopen(configfile, "r");
        minRepeats = malloc(6*sizeof(int));
        readConfig(minRepeats,  fptr);
    }


    int *N_length = malloc(sizeof(int));
    int *GC_length = malloc(sizeof(int));
    int *Total_length = malloc(sizeof(int));
    microsatelliteArray *microsatellites;
    if(infile != NULL) {
        FILE *fptr;
        fptr = fopen(infile, "r");
        microsatellites=malloc(sizeof(microsatelliteArray));
        initMicrosatelliteArray(microsatellites, SEQ_MICROSATELLITE_MEMORY_CHUNK);
        search_perfect_microsatellites(microsatellites, fptr, minRepeats, N_length
                , GC_length, Total_length);
    }

    if(outfile != NULL){
        FILE *fptr;
        fptr = fopen(outfile, "w");

        int i;
        for(i = 0; i < microsatellites->used; i++){
            fprintf(fptr, "%s\t%s\t%d\t%d\t%d\t%d\t%d\n",
                    microsatellites->array[i].sequence,
                    microsatellites->array[i].motif,
                    microsatellites->array[i].period,
                    microsatellites->array[i].repeat,
                    microsatellites->array[i].start,
                    microsatellites->array[i].end,
                    microsatellites->array[i].length);
        }
        fprintf(fptr, "# Total-length: %d\n", *Total_length);
        fprintf(fptr, "# GC-length: %d\n", *GC_length);
        fprintf(fptr, "# N-length: %d\n", *N_length);
        fclose(fptr);
        printf("File written\n");
    }

}

