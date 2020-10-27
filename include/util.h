#ifndef __UTIL
#define __UTIL
typedef enum {IntegerParm, FloatParm, DoubleParm, StringParm}
  CommandLineParmType;
void ErrorHandler(const char *msg, const char *item, int code);
void OpenFile(FILE **fp, char *name, const char *mode);
void AllocateByte(unsigned char **ptr, int len, const char *name);
void AllocateShort(short **ptr, int len, const char *name);
void AllocateInt(int **ptr, int len, const char *name);
void AllocateFloat(float **ptr, int len, const char *name);
void AllocateDouble(double **ptr, int len, const char *name);
void ReadByte(FILE *fp, unsigned char *data, int len, char *name);
void ReadShort(FILE *fp, short *data, int len, char *name);
void ReadInt(FILE *fp, int *data, int len, char *name);
void ReadFloat(FILE *fp, float *data, int len, char *name);
void ReadDouble(FILE *fp, double *data, int len, char *name);
void WriteByte(FILE *fp, unsigned char *data, int len, char *name);
void WriteShort(FILE *fp, short *data, int len, char *name);
void WriteInt(FILE *fp, int *data, int len, char *name);
void WriteFloat(FILE *fp, float *data, int len, char *name);
void WriteDouble(FILE *fp, double *data, int len, char *name);
void AverageByteToFloat(unsigned char *in, float *out, int tsize,
                        int xsize, int ysize);
void SaveFloatToImage(float *data, char *what, char *filename,
          int xsize, int ysize, int neg, int binary, int logflag);
void SaveByteToImage(unsigned char *im, const char *what, char *filename,
          int xsize, int ysize, int neg, int binary, int mask_code);
void SaveIntToImage(int *im, const char *what, char *filename,
                    int xsize, int ysize);
int Keyword(char *string, char *keyword);
int CommandLineParm(int argc, char *argv[], char *key,
  CommandLineParmType type, void *ptr, int required, char *usage);
void PrintMinAndMax(int w, int h, float *soln, char *keyword);
#endif
