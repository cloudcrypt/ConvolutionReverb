#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
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
    int chunkSize;
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
    WavHeader header;
    DataChunk dataChunk;
} WavFile;

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

    FILE *inputFile = fopen(inputFileName, "rb");
    fseek(inputFile, 0, SEEK_END);
    long size = ftell(inputFile);
    WavFile *inputWav = malloc(size);
    fseek(inputFile, 0, SEEK_SET);
    fread(inputWav, 1, size, inputFile);

    int numSamples = inputWav->dataChunk.chunkSize / (BITS_PER_SAMPLE / 8);
    float *x = malloc(sizeof(float) * numSamples);
    scaleSamples(inputWav->dataChunk.data, numSamples, x);
    // for (int i = 0; i < numSamples; i++) {
    //     *(x + i) = *(x + i) / (float)2;
    // }
    unscaleSamples(x, numSamples, inputWav->dataChunk.data);

    createWavFile(outputFileName, inputWav);

    return 0;
}

void scaleSamples(int16_t samples[], int numSamples, float scaled[]) {
    for (int i = 0; i < numSamples; i++) {
        int16_t sample = *(samples + i);
        float test = sample / (float)(maxValue + (sample < 0 ? 1 : 0));
        *(scaled + i) = sample / (float)(maxValue + (sample < 0 ? 1 : 0));
    }
}

void unscaleSamples(float scaled[], int numSamples, int16_t samples[]) {
    for (int i = 0; i < numSamples; i++) {
        float scaledSample = *(scaled + i);
        int16_t test = (int16_t)rintf(scaledSample * maxValue);
        *(samples + i) = (int16_t)rintf(scaledSample * maxValue);
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

    fwrite(&wavFile->header, sizeof(WavHeader), 1, outputFile);
    fwrite(wavFile->dataChunk.chunkId, sizeof(char), 4, outputFile);
    fwrite(&wavFile->dataChunk.chunkSize, sizeof(int), 1, outputFile);
    
    // for (int i = 0; i < numSamples; i++) {
    //     reverseEndianness(sizeof(int16_t), wavFile->dataChunk.data);
    // }
    fwrite(wavFile->dataChunk.data, sizeof(int8_t), wavFile->dataChunk.chunkSize, outputFile);
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