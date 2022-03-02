#include "header.h"
long int nsamp,naddc,bsize, bsize1, bsize2, naddt;
int headerless,obits;
long int off1,off2;
char inpfile1[128], inpfile2[128], outfile[128], gminfofile[128], gmhdrfile[128];
FILE *input, *input1, *input2, *output ;

float median(float arr[], long int n);
int strings_equal (char *string1, char *string2);

/* include these from sigproc-4.3 */
void slaCldj ( int iy, int im, int id, long double *djm, int *j );
double mjd(int year, int month, int day) ;
int read_block(FILE *input, int nbits, float *block, int nread) ;
int read_header(FILE *inputfile) ;
unsigned char charof2ints (int i, int j) ;
void char2ints (unsigned char c, int *i, int *j) ;
void char2fourints (unsigned char c, int *i, int *j, int *k, int *l);
void error_message(char *message) ;
void float2char(float *f, int n, float min, float max, unsigned char *c) ;
void float2four(float *f, int n, float min, float max, unsigned char *c) ;
void float2int(float *f, int n, int b, float min, float max, int *i) ;
void float2short(float *f, int n, float min, float max, unsigned short *s) ;
void int2float(int *i, int n, int b, float min, float max, float *f) ;
