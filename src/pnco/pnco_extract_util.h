#ifndef _PNCO_EXTRACT_UTIL_
#define _PNCO_EXTRACT_UTIL_
#include "netcdf.h"
#include "pnco.h"
#include "pnco_parser.h"
int getDivider(int ndims,const size_t *dims,size_t *divider);
int getStart(int ndims,size_t index,const size_t *divider,size_t *start);
int getIndex(int ndims,size_t *start,const size_t *divider,size_t *index);

int binary_op_uchar(BIN_OP op,unsigned char *data1,unsigned char*data2,unsigned char *dataOut);
int binary_op_schar(BIN_OP op,signed char *data1,signed char*data2,signed char *dataOut);
int binary_op_short(BIN_OP op,short *data1,short*data2,short *dataOut);
int binary_op_int(BIN_OP op,int *data1,int*data2,int *dataOut);
int binary_op_long(BIN_OP op,long *data1,long *data2,long *dataOut);
int binary_op_float(BIN_OP op,float *data1,float *data2,float *dataOut);
int binary_op_double(BIN_OP op,double *data1,double *data2,double *dataOut);
int copy_attrs(int ncid, int vlid, int ncidout, int vlidout, int nattrs);
int add_cmd_attr(int ncidout,int vlidout,int argc, char **argv,int mpi_size);

int dim_opts_handler(GLOBAL_ATTR * g,int *dimids,int ncid, int ndims,size_t *begins, size_t*ends,ptrdiff_t *strides);
int set_range_opts(char **dimList,int dimNum,int ncid,int vlid,int ndims,size_t*divider,size_t dlen,nc_type vtype,size_t *begins,size_t *ends,int index);


int init_begins_ends_strides(int ndims,size_t **begins,int *beginNum,int *is_begins_null, size_t **ends,int *endNum,int *is_ends_null,ptrdiff_t **strides,int *strideNum,int *is_strides_null,size_t *shape);
int check_begins_ends_strides(int ndims,size_t *begins,int beginNum,size_t *ends,int endNum,ptrdiff_t *strides,int strideNum,size_t *shape,char **dimsName);
int open_output_file_and_define_var(char *fileOut,int *ncidout,int *vlidout,int *dimidsOut,int append,int ndims,nc_type vtype,char **dimsName,size_t *begins,size_t *ends,ptrdiff_t *strides,size_t* shape,char *var,char *varOut,char **outDims,int outDimNum);
#endif
