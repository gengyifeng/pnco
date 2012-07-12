#ifndef _PNCO_EXTRACT_CONCAT_H_
#define _PNCO_EXTRACT_CONCAT_H_
#include "pnco_init.h"
int extract_concat(char *var,char *varOut,char **outDims,int outDimNum,size_t *begins,int beginNum,size_t *ends,int endNum,ptrdiff_t *strides,int strideNum,char **fileIns,int fileInNum,char* fileOut);

#endif

