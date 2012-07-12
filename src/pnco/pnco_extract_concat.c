#include "netcdf.h"
#include <mpi.h>
#include <assert.h>
#include "hdf5.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "pnco.h"
#include "pnco_extract_util.h"
#include "pnco_parser.h"
#include "pnco_extract_unary.h"
#define BAIL(e) do{printf("ERROR:: file: %s, line: %d, func: %s, code: %s.\n",__FILE__,__LINE__,__FUNCTION__, nc_strerror(1));return e;} while(0)

extern GLOBAL_ATTR global_attr;
/* concattenate a varaible in multiple netcdf files to a single netcdf file parallely. */
int extract_concat(char *var,char *varOut,char **outDims,int outDimNum,size_t *begins,int beginNum,size_t *ends,int endNum,ptrdiff_t *strides,int strideNum,char **fileIns,int fileInNum,char* fileOut){
	int mpi_size,mpi_rank;
	MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);

	MPI_Comm comm=MPI_COMM_WORLD;
	MPI_Info info=MPI_INFO_NULL;
	
	int ncid,vlid,ndims,dimids[NC_MAX_VAR_DIMS],nattrs;
	int ncidout,vlidout;
	int i,m,res;
	char varName[NC_MAX_NAME+1];
	nc_type vtype;	

	/* check input files start */
	if((res=nc_open_par(fileIns[0],NC_MPIIO,comm,info,&ncid))){
		printf("Cannot open input file %s\n",fileIns[0]);
		BAIL(res);
	}
	if((res=nc_inq_varid(ncid,var,&vlid))){
		printf("Cannot query variable %s in input file %s\n",var,fileIns[0]);
		BAIL(res);	
	}
	if((res=nc_inq_var(ncid,vlid,varName,&vtype,&ndims,dimids,&nattrs))){
		printf("Cannot query dimension info of varialbe %s in input file %s\n",var,fileIns[0]);
		BAIL(res);
	}
	char **dimsName=(char **)malloc(sizeof(char *)*ndims);
	size_t *shape=(size_t *)malloc(sizeof(size_t)*ndims);
	size_t *finalShape=(size_t *)malloc(sizeof(size_t)*ndims);

	for(i=0;i<ndims;++i){
		dimsName[i]=(char *)malloc(NC_MAX_NAME+1);
		nc_inq_dim(ncid,dimids[i],dimsName[i],&shape[i]);//get dimsName and dimsLen
		finalShape[i]=shape[i];
	}
	
	if((res=nc_close(ncid))){
		printf("Close input file %s failed\n",fileIns[0]);
		BAIL(res);
	}
	for(m=1;m<fileInNum;++m){
		int ndims_t;
		if((res=nc_open_par(fileIns[m],NC_MPIIO,comm,info,&ncid))){
			printf("Cannot open input file %s\n",fileIns[m]);
			BAIL(res);
		}
		if((res=nc_inq_varid(ncid,var,&vlid))){
			printf("Cannot query variable %s in input file %s\n",var,fileIns[m]);
			BAIL(res);	
		}
		if((res=nc_inq_var(ncid,vlid,varName,&vtype,&ndims_t,dimids,&nattrs))){
			printf("Cannot query dimension info of varialbe %s in input file %s\n",var,fileIns[m]);
			BAIL(res);
		}


		size_t *shape_t=(size_t *)malloc(sizeof(size_t)*ndims_t);


		for(i=0;i<ndims;++i){
			nc_inq_dim(ncid,dimids[i],NULL,&shape_t[i]);//get dimsName and dimsLen
		}
		if(ndims_t!=ndims){
			printf("Dimension size in file %s is %d, but in the first file %s is %d\n",fileIns[m],ndims_t,fileIns[0],ndims);
			return -1;
		}
		for(i=1;i<ndims;++i){
			if(shape_t[i]!=shape[i]){
				printf("Dimension shape in file %s is not equal to the first file %s\n",fileIns[m],fileIns[0]);
				return -1;
			}
		}
		finalShape[0]+=shape_t[0];
		free(shape_t);
		if((res=nc_close(ncid))){
			printf("Close input file %s failed\n",fileIns[m]);
			BAIL(res);
		}
	}
	/* check input files ends */

	int is_strides_null=0;
	int is_begins_null=0;
	int is_ends_null=0;
	
	/* initialize begins ends and strides*/
	init_begins_ends_strides(ndims,&begins,&beginNum,&is_begins_null,&ends,&endNum,&is_ends_null,&strides,&strideNum,&is_strides_null,shape);

	if(strides[0]!=1){
		printf("When cancating variables, first stride should be one\n");
		return -1;
	}
	if(!is_begins_null||!is_ends_null){
		printf("When cancat variables, begins and ends are not supported\n");
		return -1;
	}

	/*
		create output file and define dims and var
	*/
	int *dimidsOut=(int *)malloc(sizeof(int)*ndims);
	ends[0]=finalShape[0]-1;
if(-1==open_output_file_and_define_var(fileOut,&ncidout,&vlidout,dimidsOut,global_attr.append,ndims,vtype,dimsName,begins,ends,strides,shape,var,varOut,outDims,outDimNum)){return -1;}

	/* define attributes */
	copy_attrs(ncid, vlid, ncidout, vlidout, nattrs);
	add_cmd_attr(ncidout,vlidout,global_attr.argc,global_attr.argv,mpi_size);
	
	if((res=nc_enddef(ncidout))){
		printf("End define netcdf id %d in output file %s failed\n",ncidout,fileOut);
		BAIL(res);
	}

	if((res=nc_var_par_access(ncidout,vlidout,NC_INDEPENDENT))) {
		printf("Cannot set parallel access for variable %s in output file %s\n",var,fileOut);
		BAIL(res);
	}
	size_t preLen=0;
	size_t outLen;	
	for(m=0;m<fileInNum;++m){
		
		if((res=nc_open_par(fileIns[m],NC_MPIIO,comm,info,&ncid))) BAIL(res);
		if((res=nc_inq_varid(ncid,var,&vlid))) BAIL(res);
		if((res=nc_inq_var(ncid,vlid,varName,&vtype,&ndims,dimids,&nattrs))) BAIL(res);
		if((nc_var_par_access(ncid,vlid,NC_INDEPENDENT))) BAIL(res);
		for(i=0;i<ndims;++i){
		//	dimsName[i]=(char *)malloc(NC_MAX_NAME+1);
			nc_inq_dim(ncid,dimids[i],NULL,&shape[i]);//get dimsName and dimsLen
			begins[i]=0;
			ends[i]=shape[i]-1;
		}
		extract_unary_selector(mpi_rank,mpi_size,ncid,vlid,ncidout,vlidout,ndims,vtype,shape,begins,ends,strides,preLen,&outLen);
		preLen+=outLen;
/*        printf("mpi_rank %d, preLen %d, outLen %d\n",mpi_rank,preLen,outLen);*/
		if((res=nc_close(ncid))){
			printf("Close input file %s failed\n",fileIns[m]);
			BAIL(res);
		}
	}


	/* free resources */
	if(is_strides_null)
		free(strides);
	if(is_begins_null)
		free(begins);
	if(is_ends_null)
		free(ends);
	for(i=0;i<ndims;++i){
		free(dimsName[i]);
	}
	free(dimsName);
	free(dimidsOut);
	free(shape);
	free(finalShape);
	
	/*close file*/
	if((res=nc_close(ncidout))){
		printf("Close output file %s failed\n",fileOut);
		BAIL(res);
	}
	return 0;
}
