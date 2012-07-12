#ifndef _PNCO_EXTRACT_UNARY_H_
#define _PNCO_EXTRACT_UNARY_H_
#include "pnco.h"
int extract_unary(char *var,char* varOut,char** outputDims,int outputDimNum,size_t *begins,int beginNum,size_t *ends,int endNum,ptrdiff_t *strides,int strideNum,char *fileIn,char* fileOut);

int extract_unary_selector(int mpi_rank,int mpi_size, int ncid,int vlid,int ncidout,int vlidout,int ndims,nc_type vtype,size_t *shape,size_t *begins,size_t *ends, ptrdiff_t *strides,size_t preLen,size_t *outLen);

	int transfer_pos(size_t *poses,int ndims,ptrdiff_t *strides,size_t outlen,size_t *dividerOut,size_t *divider);
#endif

