#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
//#include "pnco_parser.h"
#include "pnco_create.h"
#include "netcdf.h"
#include "pnco_extract_util.h"
#include "pnco_extract_unary.h"
#include "pnco_extract_binary.h"
#include "pnco_extract_concat.h"
#include "pnco_parser.h"
extern GLOBAL_ATTR global_attr;
int test_create(){
	//create_sample("unlimited.nc");
	return 0;
}
int test_write(){
	return 0;
}
int test_read(){
	ptrdiff_t strides[2]={3,4};
	return 0;
}
int process(){
	GLOBAL_ATTR *g=&global_attr;
	//printf("fileIn[0]=%s,var1=%s,fileOut=%s\n",g->fileIns[0],g->varList[0],g->fileOut);
	if(g->fileOut==NULL)
		g->fileOut="out.nc";
	if(g->newSample){
		create_sample(g->fileOut);
		return 0;
	}
    if(g->varNum==0){
        printf("Error: There is no such specified variables!\n");
        exit(-1);
    }
	if(g->concat){
		extract_concat(g->varList[0],g->outputVar,g->outputDims,g->outputDimNum,g->begin,g->beginNum,g->end,g->endNum,g->strides,g->strideNum,g->fileIns,g->fileInNum,g->fileOut);
		return 0;
	}
	switch(g->fileInNum){
		case 1:
			extract_unary(g->varList[0],g->outputVar,g->outputDims,g->outputDimNum,g->begin,g->beginNum,g->end,g->endNum,g->strides,g->strideNum,g->fileIns[0],g->fileOut);
			break;
		case 2:
			extract_binary(g->op,g->varList[0],g->varList[0],g->outputVar,g->outputDims,g->outputDimNum,g->begin,g->beginNum,g->end,g->endNum,g->strides,g->strideNum,g->fileIns[0],g->fileIns[1],g->fileOut);
			
			break;
		default:
			printf("Please use -h to see the usage\n");
			//printf("fileInNum is %d\n",g->fileInNum);
	}
	return 0;
}
int main(int argc,char **argv){
	MPI_Init(&argc,&argv);
	pnco_parser(argc,argv);
	process();
	//test_unlimited();
	//test_create();
	//test_read();
	MPI_Finalize();
	return 0;
}
