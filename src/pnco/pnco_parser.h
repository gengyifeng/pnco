#ifndef _PNCO_PARSER_H_
#define _PNCO_PARSER_H_
#include "pnco.h"
#include "netcdf.h"
#define DELIM ","
#define DELIMA ','
typedef struct global_attr_t{
	int fileInNum;
	char **fileIns;
	int varNum;
	char **varList;
	size_t *begin;
	int beginNum;
	size_t *end;
	int endNum;
	char **dimList[NC_MAX_VAR_DIMS];
	int dimNum[NC_MAX_VAR_DIMS];
	int isRange[NC_MAX_VAR_DIMS];
	int listNum;
	ptrdiff_t *strides;
	int strideNum;
	char *fileOut;
	BIN_OP op;
	int append;
	int concat;
	int newSample;
	char *outputVar;
	char **outputDims;
	int outputDimNum;
	char **argv;
	int argc;
}GLOBAL_ATTR;


int pnco_parser(int argc,char **argv);

#endif
