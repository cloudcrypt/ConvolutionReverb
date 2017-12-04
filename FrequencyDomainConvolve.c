#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

/*  Standard sample size in bits  */
#define BITS_PER_SAMPLE     16       
#define SWAP(a,b)  tempr=(a);(a)=(b);(b)=tempr  

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

void four1(double data[], int nn, int isign);
WavFile* loadWav(char* fileName);
void scaleSamplesComplex(int16_t samples[], int numSamples, double scaled[]);
void unscaleSamplesComplex(double scaled[], int numSamples, int16_t samples[]);
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

    clock_t initial = clock();

    int inputSamples = inputWav->dataChunk->chunkSize / (BITS_PER_SAMPLE / 8);
    int impulseSamples = impulseWav->dataChunk->chunkSize / (BITS_PER_SAMPLE / 8);

    uint64_t neededLength = (inputSamples > impulseSamples ? inputSamples : impulseSamples) * 2;
    neededLength = pow(2, ceil(log(neededLength)/log(2)));

    double *x = calloc(neededLength, sizeof(double));
    scaleSamplesComplex(inputWav->dataChunk->data, inputSamples, x);
    four1(x - 1, neededLength / 2, 1);

    double *h = calloc(neededLength, sizeof(double));
    scaleSamplesComplex(impulseWav->dataChunk->data, impulseSamples, h);
    four1(h - 1, neededLength / 2, 1);

    double *f = calloc(neededLength, sizeof(double));
    for (int i = 0; i < neededLength; i+=2) {
        *(f + i) = (*(x + i) * *(h + i)) - (*(x + i + 1) * *(h + i + 1));
        *(f + i + 1) = (*(x + i + 1) * *(h + i)) + (*(x + i) * *(h + i + 1));
    }
    four1(f - 1, neededLength / 2, -1);

    for (int i = 0; i < neededLength; i+=2) {
        *(f + i) = *(f + i) / (double)neededLength;
        *(f + i + 1) = *(f + i + 1) / (double)neededLength;
    }

    int convolvedSamples = inputSamples + impulseSamples - 1;
    int addedBytes = (convolvedSamples - inputSamples) * (BITS_PER_SAMPLE / 8);
    inputWav->header->chunkSize += addedBytes;
    inputWav->dataChunk->chunkSize += addedBytes;
    inputWav->dataChunk = realloc(inputWav->dataChunk, inputWav->dataChunk->chunkSize + 8);
    unscaleSamplesComplex(f, convolvedSamples, inputWav->dataChunk->data);

    double elapsed = (clock() - initial) / (double)CLOCKS_PER_SEC;
    printf("Convolution completed in %.4f seconds\n", elapsed);

    createWavFile(outputFileName, inputWav);

    return 0;
}

void four1(double data[], int nn, int isign) {
    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;

    n = nn << 1;
    j = 1;

    for (i = 1; i < n; i += 2) {
        if (j > i) {
            SWAP(data[j], data[i]);
            SWAP(data[j+1], data[i+1]);
        }
        m = nn;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }

    mmax = 2;
    while (n > mmax) {
        istep = mmax << 1;
        theta = isign * (6.28318530717959 / mmax);
        wtemp = sin(0.5 * theta);
        wpr = -2.0 * wtemp * wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for (m = 1; m < mmax; m += 2) {
            for (i = m; i <= n; i += istep) {
                j = i + mmax;
                tempr = wr * data[j] - wi * data[j+1];
                tempi = wr * data[j+1] + wi * data[j];
                data[j] = data[i] - tempr;
                data[j+1] = data[i+1] - tempi;
                data[i] += tempr;
                data[i+1] += tempi;
            }
            wr = (wtemp = wr) * wpr - wi * wpi + wr;
            wi = wi * wpr + wtemp * wpi + wi;
        }
        mmax = istep;
    }
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

void scaleSamplesComplex(int16_t samples[], int numSamples, double scaled[]) {
    for (int i = 0; i < (numSamples * 2); i+=2) {
        int16_t sample = *(samples + (i / 2));
        *(scaled + i) = sample / (float)(maxValue + (sample < 0 ? 1 : 0));
        *(scaled + i + 1) = 0.0;
    }
}

void unscaleSamplesComplex(double scaled[], int numSamples, int16_t samples[]) {
    for (int i = 0; i < numSamples; i++) {
        float scaledSample = *(scaled + (i * 2));
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

    fwrite(wavFile->header, wavFile->header->fmtChunk.chunkSize + 20, 1, outputFile);
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