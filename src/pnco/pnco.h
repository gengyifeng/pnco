#ifndef _PNCO_H_
#define _PNCO_H_
#include "netcdf.h"

typedef struct dims_t{
	int size;
	size_t *shape;
	char **names;
}DIMS;
typedef struct var_t{
	nc_type type;
	char *varName;
	DIMS dims;
}VAR;
typedef struct box_t{
	int size;
	size_t *start;
	size_t *count;
	size_t *stride;
}BOX;
typedef enum{
	ADD,SUB,DIV,MUL,MOD,AVG
}BIN_OP;

int init_dims(DIMS *dims,int size, size_t * arr,char **names);
int init_var(VAR *var,nc_type type,char *varName,DIMS *dims);

int init_box(BOX *box,int size,size_t *start,size_t *count);
int init_box_with_stride(BOX *box,int size,size_t *start,size_t *count,size_t *stride);

#endif
