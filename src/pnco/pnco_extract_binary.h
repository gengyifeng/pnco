#ifndef _PNCO_EXTRACT_BINARY_H_
#define _PNCO_EXTRACT_BINRAY_H_
#include "pnco.h"
int extract_binary(BIN_OP op,char *var1,char* var2,char *varOut,char **outDims,int outDimNum,size_t *begins,int beginNum,size_t *ends,int endNum,ptrdiff_t *strides,int strideNum,char *fileIn1,char *fileIn2,char* fileOut);
#endif
