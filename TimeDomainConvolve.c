#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*  Standard sample size in bits  */
#define BITS_PER_SAMPLE     16         

typedef struct FmtChunk {
    char chunkId[4];
    int chunkSize;
    int16_t audioFormat;
    int16_t numChannels;
    int sampleRate;
    int byteRate;
    int16_t blockAlign;
    int16_t bitsPerSample;
} FmtChunk;

typedef struct DataChunk {
    char chunkId[4];
    uint32_t chunkSize;
    int16_t data[];
    //int16_t *data;
} DataChunk;

typedef struct WavHeader {
    char chunkId[4];
    int chunkSize;
    char format[4];
    FmtChunk fmtChunk;
} WavHeader;

typedef struct WavFile {
    char *fileBuf;
    WavHeader *header;
    DataChunk *dataChunk;
} WavFile;

void convolve(float x[], int N, float h[], int M, float y[], int P);
WavFile* loadWav(char* fileName);
void scaleSamples(int16_t samples[], int numSamples, float scaled[]);
void unscaleSamples(float scaled[], int numSamples, int16_t samples[]);
void createWavFile(char* fileName, WavFile *wavFile);
void reverseEndianness(int size, void *value);

int maxValue;

int main(int argc, char *argv[]) {

    maxValue = (int)pow(2.0, (double)BITS_PER_SAMPLE - 1) - 1;
    
    char *inputFileName = NULL;
    char *IRfileName = NULL;
    char *outputFileName = NULL;
	
    /*  Process the command line arguments  */
    if (argc == 4) {
        /*  Set a pointer to the output filename  */
        inputFileName = argv[1];
        IRfileName = argv[2];
        outputFileName = argv[3];
    }
    else {
        /*  The user did not supply the correct number of command-line
            arguments.  Print out a usage message and abort the program.  */
        fprintf(stderr, "Usage: TimeDomainConvolve inputfile IRfile outputfile\n");
        exit(-1);
    }

    WavFile *inputWav = loadWav(inputFileName);
    WavFile *impulseWav = loadWav(IRfileName);

    int inputSamples = inputWav->dataChunk->chunkSize / (BITS_PER_SAMPLE / 8);
    float *x = malloc(sizeof(float) * inputSamples);
    scaleSamples(inputWav->dataChunk->data, inputSamples, x);

    int impulseSamples = impulseWav->dataChunk->chunkSize / (BITS_PER_SAMPLE / 8);
    float *h = malloc(sizeof(float) * impulseSamples);
    scaleSamples(impulseWav->dataChunk->data, impulseSamples, h);

    int convolvedSamples = inputSamples + impulseSamples - 1;
    float *y = malloc(sizeof(float) * convolvedSamples);

    convolve(x, inputSamples, h, impulseSamples, y, convolvedSamples);

    int addedBytes = (convolvedSamples - inputSamples) * (BITS_PER_SAMPLE / 8);
    inputWav->header->chunkSize += addedBytes;
    inputWav->dataChunk->chunkSize += addedBytes;
    inputWav->dataChunk = realloc(inputWav->dataChunk, inputWav->dataChunk->chunkSize + 8);
    unscaleSamples(y, convolvedSamples, inputWav->dataChunk->data);

    createWavFile(outputFileName, inputWav);

    return 0;
}

void convolve(float x[], int N, float h[], int M, float y[], int P) {
    int n, m;

    /*  Make sure the output buffer is the right size: P = N + M - 1  */
    if (P != (N + M - 1)) {
        printf("Output signal vector is the wrong size\n");
        printf("It is %-d, but should be %-d\n", P, (N + M - 1));
        printf("Aborting convolution\n");
        return;
    }

    /*  Clear the output buffer y[] to all zero values  */  
    for (n = 0; n < P; n++)
        y[n] = 0.0;

    /*  Do the convolution  */
    /*  Outer loop:  process each input value x[n] in turn  */
    uint64_t totalIterations = N * M;
    uint64_t performedIterations = 0;
    for (n = 0; n < N; n++) {
        /*  Inner loop:  process x[n] with each sample of h[]  */
        for (m = 0; m < M; m++) {
            y[n+m] += x[n] * h[m];
        }
        performedIterations += M;
        printf("Processing %llu/%llu...\r", performedIterations, totalIterations);
    }
    printf("\n");
}

WavFile* loadWav(char* fileName) {
    WavFile *wavFile = malloc(sizeof(WavFile));
    FILE *file = fopen(fileName, "rb");
    fseek(file, 0, SEEK_END);
    long size = ftell(file);
    wavFile->fileBuf = malloc(size);
    fseek(file, 0, SEEK_SET);
    fread(wavFile->fileBuf, 1, size, file);
    fclose(file);

    wavFile->header = (WavHeader *)wavFile->fileBuf;
    DataChunk* dataChunk = (DataChunk *)(wavFile->fileBuf + 20 + wavFile->header->fmtChunk.chunkSize);
    wavFile->dataChunk = malloc(dataChunk->chunkSize + 8);
    memcpy(wavFile->dataChunk, dataChunk, dataChunk->chunkSize + 8);
    wavFile->fileBuf = realloc(wavFile->fileBuf, wavFile->header->chunkSize + 8 - (dataChunk->chunkSize + 8));

    return wavFile;
}

void scaleSamples(int16_t samples[], int numSamples, float scaled[]) {
    for (int i = 0; i < numSamples; i++) {
        int16_t sample = *(samples + i);
        *(scaled + i) = sample / (float)(maxValue + (sample < 0 ? 1 : 0));
    }
}

void unscaleSamples(float scaled[], int numSamples, int16_t samples[]) {
    for (int i = 0; i < numSamples; i++) {
        float scaledSample = *(scaled + i);
        *(samples + i) = (int16_t)rintf(scaledSample * (maxValue + (scaledSample < 0 ? 1 : 0)));
    }
}

void createWavFile(char* fileName, WavFile *wavFile) {
    FILE *outputFile = fopen(fileName, "wb");

    //int dataChunkSize = wavFile->dataChunk.chunkSize / 2;

    // reverseEndianness(sizeof(int), &wavFile->header.chunkSize);
    // reverseEndianness(sizeof(int), &wavFile->header.fmtChunk.chunkSize);
    // reverseEndianness(sizeof(int16_t), &wavFile->header.fmtChunk.audioFormat);
    // reverseEndianness(sizeof(int16_t), &wavFile->header.fmtChunk.numChannels);
    // reverseEndianness(sizeof(int), &wavFile->header.fmtChunk.sampleRate);
    // reverseEndianness(sizeof(int), &wavFile->header.fmtChunk.byteRate);
    // reverseEndianness(sizeof(int16_t), &wavFile->header.fmtChunk.blockAlign);
    // reverseEndianness(sizeof(int16_t), &wavFile->header.fmtChunk.bitsPerSample); 
    // reverseEndianness(sizeof(int), &wavFile->dataChunk.chunkSize);

    fwrite(wavFile->header, sizeof(WavHeader), 1, outputFile);
    fwrite(wavFile->dataChunk->chunkId, sizeof(char), 4, outputFile);
    fwrite(&wavFile->dataChunk->chunkSize, sizeof(int), 1, outputFile);
    
    // for (int i = 0; i < numSamples; i++) {
    //     reverseEndianness(sizeof(int16_t), wavFile->dataChunk.data);
    // }
    fwrite(wavFile->dataChunk->data, sizeof(int8_t), wavFile->dataChunk->chunkSize, outputFile);
    fclose(outputFile);
}

void reverseEndianness(int size, void *value) {
    uint8_t reversed[4];
    for (int i = 0; i < size; i++) {
        reversed[i] = ((uint8_t *)value)[size - 1 - i];
    }
    for (int i = 0; i < size; i++) {
        ((uint8_t *)value)[i] = reversed[i];
    }
}