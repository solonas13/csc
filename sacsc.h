#ifdef _USE_64
typedef int64_t INT;
#endif

#ifdef _USE_32
typedef int32_t INT;
#endif

unsigned int LCParray ( unsigned char *text, INT n, INT * SA, INT * ISA, INT * LCP );
void partitioning ( INT i, INT j, INT f, INT m, INT * mf, INT * ind );
unsigned int circular_sequence_comparison (  unsigned char * x, unsigned char * y, struct TSwitch  sw, unsigned int * rotation, unsigned int * distance );
