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
/*operator two varialbe in two netcdf file and write the result to output netcdf file.*/
int extract_binary(BIN_OP op,char *var1,char* var2,char *varOut,char **outDims,int outDimNum,size_t *begins,int beginNum,size_t *ends,int endNum,ptrdiff_t *strides,int strideNum,char *fileIn1,char *fileIn2,char* fileOut){
	int mpi_size,mpi_rank;
	MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);

	MPI_Comm comm=MPI_COMM_WORLD;
	MPI_Info info=MPI_INFO_NULL;
	
	int ncid1,vlid1,ndims1,dimids1[NC_MAX_VAR_DIMS],nattrs1;
	int ncid2,vlid2,ndims2,dimids2[NC_MAX_VAR_DIMS],nattrs2;
	int ncidout,vlidout;
	int i,res;
	char varName[NC_MAX_NAME+1];
	nc_type vtype1;
	nc_type vtype2;


	/* 
		open input file1
	*/
	if((res=nc_open_par(fileIn1,NC_MPIIO,comm,info,&ncid1))){
		printf("Cannot open input file %s\n",fileIn1);
		BAIL(res);
	}
	if((res=nc_inq_varid(ncid1,var1,&vlid1))){
		printf("Cannot query variable %s in input file %s\n",var1,fileIn1);
		BAIL(res);	
	}
	if((res=nc_inq_var(ncid1,vlid1,varName,&vtype1,&ndims1,dimids1,&nattrs1))){
		printf("Cannot query dimension info of varialbe %s in input file %s\n",var1,fileIn1);
		BAIL(res);
	}
	/* 
		open input file2
	*/
	if((res=nc_open_par(fileIn2,NC_MPIIO,comm,info,&ncid2))){
		printf("Cannot open input file %s\n",fileIn2);
		BAIL(res);
	}
	if((res=nc_inq_varid(ncid2,var2,&vlid2))){
		printf("Cannot query variable %s in input file %s\n",var2,fileIn2);
		BAIL(res);	
	}
	if((res=nc_inq_var(ncid2,vlid2,varName,&vtype2,&ndims2,dimids2,&nattrs2))){
		printf("Cannot query dimension info of varialbe %s in input file %s\n",var2,fileIn2);
		BAIL(res);
	}
	/* check dimension size*/
	if(ndims1!=ndims2){
		printf("Dimension size of file %s and file %s is not equal\n",fileIn1,fileIn2);
		return -1;
	}

	char **dimsName=(char **)malloc(sizeof(char *)*ndims1);
	size_t *shape1=(size_t *)malloc(sizeof(size_t)*ndims1);
	size_t *shape2=(size_t *)malloc(sizeof(size_t)*ndims1);

	for(i=0;i<ndims1;++i){
		dimsName[i]=(char *)malloc(NC_MAX_NAME+1);
		nc_inq_dim(ncid1,dimids1[i],dimsName[i],&shape1[i]);//get dimsName and dimsLen
		nc_inq_dim(ncid2,dimids2[i],NULL,&shape2[i]);//get dimsName and dimsLen
	}

	/* check dimension shape */ 
	for(i=0;i<ndims1;++i){
		if(shape1[i]!=shape2[i]){
			printf("Dimension shape of file %s and file %s is not equal\n",var1,var2);
			return -1;
		}
	}

	/* check variable data type */
	if(vtype1!=vtype2){
		printf("Data type of variable %s and variable %s is not equal\n",var1,var2);
		return -1;
	}
	
	int is_strides_null=0;
	int is_begins_null=0;
	int is_ends_null=0;


	/* initialize begins ends and strides*/
	init_begins_ends_strides(ndims1,&begins,&beginNum,&is_begins_null,&ends,&endNum,&is_ends_null,&strides,&strideNum,&is_strides_null,shape1);
	/* deal with -d and -r arguments */
	dim_opts_handler(&global_attr,dimids1,ncid1,ndims1,begins,ends,strides);

	/* check begins ends and strides */
	check_begins_ends_strides(ndims1,begins,beginNum,ends,endNum,strides,strideNum,shape1,dimsName);
	/* set independent parallel access for input files */
	if((res=nc_var_par_access(ncid1,vlid1,NC_INDEPENDENT))){
		printf("Cannot set parallel access for variable %s in input file %s\n",var1,fileIn1);
		BAIL(res);
	}
	if((res=nc_var_par_access(ncid2,vlid2,NC_INDEPENDENT))){
		printf("Cannot set parallel access for variable %s in input file %s\n",var2,fileIn2);
		BAIL(res);
	}

	/*
		create output file and define dims and var
	*/
	int *dimidsOut=(int *)malloc(sizeof(int)*ndims1);
	open_output_file_and_define_var(fileOut,&ncidout,&vlidout,dimidsOut,global_attr.append,ndims1,vtype1,dimsName,begins,ends,strides,shape1,var1,varOut,outDims,outDimNum);

	/* define attributes */
	copy_attrs(ncid1, vlid1, ncidout, vlidout, nattrs1);
	add_cmd_attr(ncidout,vlidout,global_attr.argc,global_attr.argv,mpi_size);
	
	if((res=nc_enddef(ncidout))){
		printf("End define netcdf id %d in output file %s failed\n",ncidout,fileOut);
		BAIL(res);
	}

	if((res=nc_var_par_access(ncidout,vlidout,NC_INDEPENDENT))) {
		printf("Cannot set parallel access for variable %s in output file %s\n",var1,fileOut);
		BAIL(res);
	}

	extract_binary_selector(op,mpi_rank,mpi_size,ncid1,vlid1,ncid2,vlid2,ncidout,vlidout,ndims1,vtype1,shape1,begins,ends,strides,0L,NULL);


	/* free resources */
	if(is_strides_null)
		free(strides);
	if(is_begins_null)
		free(begins);
	if(is_ends_null)
		free(ends);
	for(i=0;i<ndims1;++i){
		free(dimsName[i]);
	}
	free(dimsName);
	free(dimidsOut);
	free(shape1);
	free(shape2);

	/*close file*/
	if((res=nc_close(ncid1))){
		printf("Close input file %s failed\n",fileIn1);
		BAIL(res);
	}
	if((res=nc_close(ncid2))){
		printf("Close input file %s failed\n",fileIn2);
		BAIL(res);
	}
	if((res=nc_close(ncidout))){
		printf("Close output file %s fileOut failed\n",fileOut);
		BAIL(res);
	}
	return 0;
}


/* This function treats multi-dimension data as one dimension array, then distributes the array range
 * to each proccess according to the mpi_rank. Due to the underline pararllel netcdf library, only
 * nc_get_var1_* and nc_put_var2_* can be used to read or write single element. This function works
 * but has low performance, so it should be used properly when extract_binary_1D or
 * extract_binary_1D_all does not work parallelly well.*/
int extract_binary_single(BIN_OP op,int mpi_rank,int mpi_size,int ncid1,int vlid1, int ncid2,int vlid2,int ncidout,int vlidout,int ndims,nc_type vtype,size_t *shape,size_t *begins,size_t *ends, ptrdiff_t *strides,size_t preLen,size_t *outLen){
	
	int i,j,res;
	size_t *divider=(size_t *)malloc(sizeof(size_t)*ndims);    //input divider 
	size_t *dividerOut=(size_t *)malloc(sizeof(size_t)*ndims); // output divider
	size_t *start=(size_t*)malloc(sizeof(size_t)*ndims);     //start position for reading element from input file
	size_t *startOut=(size_t*)malloc(sizeof(size_t)*ndims);  //start position for writing element to output file
	size_t *shapeOut=(size_t*)malloc(sizeof(size_t)*ndims);  //output dimension shape

	int lenOut=1;
	for(i=0;i<ndims;++i){
		shapeOut[i]=(ends[i]-begins[i])/strides[i]+1;	
		lenOut*=shapeOut[i];
	}
	if(outLen!=NULL)
		*outLen=lenOut;
	getDivider(ndims,shape,divider);
	getDivider(ndims,shapeOut,dividerOut);

	/* decide element boundary for each mpi process */
	size_t beginOut;
	size_t endOut;
	if(lenOut>=mpi_size){                              
		beginOut=mpi_rank*(lenOut/mpi_size);
		if(mpi_rank!=mpi_size-1)
			endOut=(mpi_rank+1)*(lenOut/mpi_size);
		else
			endOut=lenOut;
	}else{ //mpi_size is bigger than lenOut
		if(mpi_rank<lenOut){
			beginOut=mpi_rank;
			endOut=mpi_rank+1;
		}else{
			beginOut=0;
			endOut=0;
		}
	}

	printf("mpi_rank %d, beginOut %d, endOut %d\n",mpi_rank,beginOut,endOut);
	void *data1=malloc(sizeof(double));
	void *data2=malloc(sizeof(double));
	void *dataOut=malloc(sizeof(double));

	size_t rem,remIn;
	for(i=beginOut;i<endOut;++i){
		rem=i+preLen;
		remIn=i;
		for(j=0;j<ndims;++j){
			startOut[j]=rem/dividerOut[j];
			start[j]=begins[j]+(remIn/dividerOut[j])*strides[j];
			rem=rem%dividerOut[j];
			remIn=remIn%dividerOut[j];
		}

	switch(vtype){
		case NC_BYTE:
			if((res=nc_get_var1_uchar(ncid1,vlid1,start,data1)))
				BAIL(res);
			if((res=nc_get_var1_uchar(ncid2,vlid2,start,data2)))
				BAIL(res);
			binary_op_uchar(op,data1,data2,dataOut);
			if((res=nc_put_var1_uchar(ncidout,vlidout,startOut,(unsigned char *)dataOut)))
				BAIL(res);
			break;
		case NC_CHAR:
			if((res=nc_get_var1_schar(ncid1,vlid1,start,(signed char *)data1)))
				BAIL(res);
			if((res=nc_get_var1_schar(ncid2,vlid2,start,(signed char *)data2)))
				BAIL(res);
			binary_op_schar(op,data1,data2,dataOut);
			if((res=nc_put_var1_schar(ncidout,vlidout,startOut,(signed char *)dataOut)))
				BAIL(res);
			break;
		case NC_SHORT:
			if((res=nc_get_var1_short(ncid1,vlid1,start,data1)))
				BAIL(res);
			if((res=nc_get_var1_short(ncid2,vlid2,start,data2)))
				BAIL(res);
			binary_op_short(op,data1,data2,dataOut);
			if((res=nc_put_var1_short(ncidout,vlidout,startOut,(short *)dataOut)))
				BAIL(res);
			break;
		case NC_INT:
			if((res=nc_get_var1_int(ncid1,vlid1,start,(int *)data1)))
				BAIL(res);
			if((res=nc_get_var1_int(ncid2,vlid2,start,(int *)data2)))
				BAIL(res);
			binary_op_int(op,data1,data2,dataOut);
			if((res=nc_put_var1_int(ncidout,vlidout,startOut,(int *)dataOut)))
				BAIL(res);
			break;
		case NC_FLOAT:
			if((res=nc_get_var1_float(ncid1,vlid1,start,data1)))
				BAIL(res);
			if((res=nc_get_var1_float(ncid2,vlid2,start,data2)))
				BAIL(res);
			binary_op_float(op,data1,data2,dataOut);
			if((res=nc_put_var1_float(ncidout,vlidout,startOut,(float *)dataOut)))
				BAIL(res);
			break;
		case NC_DOUBLE:
			if((res=nc_get_var1_double(ncid1,vlid1,start,data1)))
				BAIL(res);
			if((res=nc_get_var1_double(ncid2,vlid2,start,data2)))
				BAIL(res);
			binary_op_double(op,data1,data2,dataOut);
			if((res=nc_put_var1_double(ncidout,vlidout,startOut,(double *)dataOut)))
				BAIL(res);
			break;
		default:
			printf("Unknown data type\n");
			}
		}

	/* free resourses */
	free(divider);
	free(dividerOut);
	free(start);
	free(startOut);
	free(shapeOut);
	free(data1);
	free(data2);
	free(dataOut);
	return 0;
}

/* This function use the  first dimension to parallel. The difference between extract_binary_1D and
 * extract unary_1D_all is that extract_binary_1D in each iteration reads and writes a row like
 * v[1][?][?].
 * but extract_binary_1D_all reads and writes all needed data only once, there is no iteration.
 * */
int extract_binary_1D(BIN_OP op, int mpi_rank,int mpi_size, int ncid1,int vlid1, int ncid2, int vlid2, int ncidout,int vlidout,int ndims,nc_type vtype,size_t *shape,size_t *begins,size_t *ends, ptrdiff_t *strides,size_t preLen,size_t *outLen){
	
	int i,j,res;
	size_t *divider=(size_t *)malloc(sizeof(size_t)*ndims);    //input divider 
	size_t *dividerOut=(size_t *)malloc(sizeof(size_t)*ndims); // output divider
	size_t *start=(size_t*)malloc(sizeof(size_t)*ndims);     //start position for reading element from input file
	size_t *count=(size_t*)malloc(sizeof(size_t)*ndims);
	size_t *countOut=(size_t *)malloc(sizeof(size_t)*ndims);
	size_t *startOut=(size_t*)malloc(sizeof(size_t)*ndims);  //start position for writing element to output file
	size_t *shapeOut=(size_t*)malloc(sizeof(size_t)*ndims);  //output dimension shape
	size_t startOut0;
	size_t countOut0;
	int len0=1;
	int lenOut=1;
	for(i=0;i<ndims;++i){
		shapeOut[i]=(ends[i]-begins[i])/strides[i]+1;	
		lenOut*=shapeOut[i];
		if(i>0){
			startOut[i]=0;
			countOut[i]=shapeOut[i];
			start[i]=begins[i];
			count[i]=ends[i]-begins[i]+1;
			len0*=ends[i]-begins[i]+1;
		}
	}
	if(outLen!=NULL)
		*outLen=lenOut;
	if(shapeOut[0]>=mpi_size){
		startOut0=mpi_rank*(shapeOut[0]/mpi_size);
		//start[0]=begins[0]+startOut[0]*stride[0];
		if(mpi_rank!=mpi_size-1){
			countOut0=(shapeOut[0]/mpi_size);
		}else{
			countOut0=shapeOut[0]-startOut0;
		}
	}else{
		if(mpi_rank<shapeOut[0]){
			startOut0=mpi_rank;
			countOut0=1;
		}else{
			return 0;
		}
	}
	int dataEnd=lenOut/shapeOut[0];
	void* data1=(void*)malloc(sizeof(double)*len0);
	void* data2=(void*)malloc(sizeof(double)*len0);
	void* dataOut=(void*)malloc(sizeof(double)*dataEnd);
	size_t* poses=(size_t*)malloc(sizeof(size_t)*dataEnd);
/*    size_t *divider=(size_t*)malloc(sizeof(size_t)*(ndims-1));*/
/*    size_t *dividerOut=(size_t*)malloc(sizeof(size_t)*(ndims-1));*/
	getDivider(ndims,count,divider);
	getDivider(ndims,countOut,dividerOut);
	transfer_pos(poses,ndims-1,strides,dataEnd,dividerOut,divider);
	printf("mpi_rank %d,countOut0 %d\n",mpi_rank,countOut0);
	for(i=0;i<countOut0;++i){
		startOut[0]=startOut0+i+preLen/len0;
		countOut[0]=1;
		start[0]=begins[0]+(startOut0+i)*strides[0];
		count[0]=1;

	switch(vtype){
		case NC_BYTE:
			if((res=nc_get_vara_uchar(ncid1,vlid1,start,count,(unsigned char *)data1)))
				BAIL(res);
			if((res=nc_get_vara_uchar(ncid2,vlid2,start,count,(unsigned char *)data2)))
				BAIL(res);
			for(j=0;j<dataEnd;++j){
                binary_op_uchar(op,&((unsigned char*)data1)[j],&((unsigned char*)data2)[j],&((unsigned char *)dataOut)[j]);
/*                ((unsigned char *)dataOut)[j]=((unsigned char *)data1)[poses[j]]+((unsigned char *)data2)[poses[j]];*/
			}
			if((res=nc_put_vara_uchar(ncidout,vlidout,startOut,countOut,(unsigned char *)dataOut)))
				BAIL(res);

			break;
		case NC_CHAR:
			if((res=nc_get_vara_schar(ncid1,vlid1,start,count,(signed char *)data1)))
				BAIL(res);
			if((res=nc_get_vara_schar(ncid2,vlid2,start,count,(signed char *)data2)))
				BAIL(res);
			for(j=0;j<dataEnd;++j){
                binary_op_schar(op,&((signed char*)data1)[j],&((signed char*)data2)[j],&((signed char *)dataOut)[j]);
			}
			if((res=nc_put_vara_schar(ncidout,vlidout,startOut,countOut,(signed char *)dataOut)))
				BAIL(res);
		
			break;
		case NC_SHORT:
			if((res=nc_get_vara_short(ncid1,vlid1,start,count,(short *)data1)))
				BAIL(res);
			if((res=nc_get_vara_short(ncid2,vlid2,start,count,(short *)data2)))
				BAIL(res);
			for(j=0;j<dataEnd;++j){
                binary_op_short(op,&((short*)data1)[j],&((short*)data2)[j],&((short *)dataOut)[j]);
			}
			if((res=nc_put_vara_short(ncidout,vlidout,startOut,countOut,(short *)dataOut)))
				BAIL(res);
			
			break;
		case NC_INT:
/*            printf("mpi_rank %d,start[0] %d,start[1] %d,count[0] %d count[1] %d\n",mpi_rank,start[0],start[1],count[0],count[1]);*/
/*            printf("mpi_rank %d,startOut[0] %d,startOut[1] %d,countOut[0] %d countOut[1] %d\n",mpi_rank,startOut[0],startOut[1],countOut[0],countOut[1]);*/
			if((res=nc_get_vara_int(ncid1,vlid1,start,count,(int *)data1)))
				BAIL(res);
			if((res=nc_get_vara_int(ncid2,vlid2,start,count,(int *)data2)))
				BAIL(res);
			for(j=0;j<dataEnd;++j){
                binary_op_int(op,&((int*)data1)[j],&((int*)data2)[j],&((int *)dataOut)[j]);
			}
			if((res=nc_put_vara_int(ncidout,vlidout,startOut,countOut,(int *)dataOut)))
				BAIL(res);
			
			break;
		case NC_FLOAT:
			if((res=nc_get_vara_float(ncid1,vlid1,start,count,(float *)data1)))
				BAIL(res);
			if((res=nc_get_vara_float(ncid2,vlid2,start,count,(float *)data2)))
				BAIL(res);
			for(j=0;j<dataEnd;++j){
                binary_op_float(op,&((float*)data1)[j],&((float*)data2)[j],&((float *)dataOut)[j]);
			}
			if((res=nc_put_vara_float(ncidout,vlidout,startOut,countOut,(float *)dataOut)))
				BAIL(res);
			
			break;
		case NC_DOUBLE:
			if((res=nc_get_vara_double(ncid1,vlid1,start,count,(double *)data1)))
				BAIL(res);
			if((res=nc_get_vara_double(ncid2,vlid2,start,count,(double *)data2)))
				BAIL(res);
			for(j=0;j<dataEnd;++j){
                binary_op_double(op,&((double*)data1)[j],&((double*)data2)[j],&((double *)dataOut)[j]);
			}
			if((res=nc_put_vara_double(ncidout,vlidout,startOut,countOut,(double *)dataOut)))
				BAIL(res);
			
			break;
		default:
			printf("Unknown data type\n");
			}
		}

	/*free resourses*/
	free(divider);
	free(dividerOut);
	free(start);
	free(startOut);
	free(shapeOut);
	free(data1);
	free(data2);
	free(dataOut);
	free(poses);
	return 0;
}

/* This function use the  first dimension to parallel and reads and writes all needed data only once, there is no iteration. stride[0] is not supported.
 * */
int extract_binary_1D_all(BIN_OP op,int mpi_rank,int mpi_size, int ncid1,int vlid1,int ncid2, int vlid2,int ncidout,int vlidout,int ndims,nc_type vtype,size_t *shape,size_t *begins,size_t *ends,size_t preLen,size_t *outLen){
	
	int i,j,res;
	size_t *start=(size_t*)malloc(sizeof(size_t)*ndims);     //start position for reading element from input file
	size_t *countOut=(size_t *)malloc(sizeof(size_t)*ndims);
	size_t *startOut=(size_t*)malloc(sizeof(size_t)*ndims);  //start position for writing element to output file
	size_t *shapeOut=(size_t*)malloc(sizeof(size_t)*ndims);  //output dimension shape
	size_t startOut0;
	size_t countOut0;
	int len0=1;
	int lenOut=1;
	for(i=0;i<ndims;++i){
/*            shapeOut[0]=(ends[0]-begins[0])/strides[0]+1;*/
		shapeOut[i]=ends[i]-begins[i]+1;	
		lenOut*=shapeOut[i];
		if(i>0){
			startOut[i]=0;
			countOut[i]=shapeOut[i];
			start[i]=begins[i];
/*            count[i]=ends[i]-begins[i]+1;*/
			len0*=ends[i]-begins[i]+1;
		}
	}
	if(outLen!=NULL)
		*outLen=lenOut;
	if(shapeOut[0]>=mpi_size){
		startOut0=mpi_rank*(shapeOut[0]/mpi_size);
		//start[0]=begins[0]+startOut[0]*stride[0];
		if(mpi_rank!=mpi_size-1){
			countOut0=(shapeOut[0]/mpi_size);
		}else{
			countOut0=shapeOut[0]-startOut0;
		}
	}else{
		if(mpi_rank<shapeOut[0]){
			startOut0=mpi_rank;
			countOut0=1;
		}else{
			return 0;
		}
	}
	int dataEnd=countOut0*len0;
	printf("mpi_rank %d,countOut0 %d\n",mpi_rank,countOut0);
	startOut[0]=startOut0+preLen/len0;
	countOut[0]=countOut0;
/*        start[0]=begins[0]+startOut0*strides[0];*/
	start[0]=begins[0]+startOut0;
	switch(vtype){
		case NC_BYTE:
			{
			unsigned char* data1=(unsigned char *)malloc(sizeof(unsigned char)*dataEnd);
			unsigned char* data2=(unsigned char *)malloc(sizeof(unsigned char)*dataEnd);
			unsigned char* dataOut=(unsigned char *)malloc(sizeof(unsigned char)*dataEnd);
			if((res=nc_get_vara_uchar(ncid1,vlid1,start,countOut,data1)))
				BAIL(res);
			if((res=nc_get_vara_uchar(ncid2,vlid2,start,countOut,data2)))
				BAIL(res);
			for(j=0;j<dataEnd;++j){
				binary_op_uchar(op,&data1[j],&data2[j],&dataOut[j]);
			}
			if((res=nc_put_vara_uchar(ncidout,vlidout,startOut,countOut,(unsigned char *)dataOut)))
				BAIL(res);
			free(data1);
			free(data2);
			free(dataOut);
			}
			break;
		case NC_CHAR:
			{
			signed char* data1=(signed char *)malloc(sizeof(signed char)*dataEnd);
			signed char* data2=(signed char *)malloc(sizeof(signed char)*dataEnd);
			signed char* dataOut=(signed char *)malloc(sizeof(signed char)*dataEnd);
			if((res=nc_get_vara_schar(ncid1,vlid1,start,countOut,data1)))
				BAIL(res);
			if((res=nc_get_vara_schar(ncid2,vlid2,start,countOut,data2)))
				BAIL(res);
			for(j=0;j<dataEnd;++j){
				binary_op_schar(op,&data1[j],&data2[j],&dataOut[j]);
			}
			if((res=nc_put_vara_schar(ncidout,vlidout,startOut,countOut,(signed char *)dataOut)))
				BAIL(res);
			free(data1);
			free(data2);
			free(dataOut);
			}
			break;
		case NC_SHORT:
			{
			short* data1=(short *)malloc(sizeof(short)*dataEnd);
			short* data2=(short *)malloc(sizeof(short)*dataEnd);
			short* dataOut=(short *)malloc(sizeof(short)*dataEnd);
			if((res=nc_get_vara_short(ncid1,vlid1,start,countOut,data1)))
				BAIL(res);
			if((res=nc_get_vara_short(ncid2,vlid2,start,countOut,data2)))
				BAIL(res);
			for(j=0;j<dataEnd;++j){
				binary_op_short(op,&data1[j],&data2[j],&dataOut[j]);
			}
			if((res=nc_put_vara_short(ncidout,vlidout,startOut,countOut,(short *)dataOut)))
				BAIL(res);
			free(data1);
			free(data2);
			free(dataOut);
			}
			break;
		case NC_INT:
/*            printf("mpi_rank %d,start[0] %d,start[1] %d,count[0] %d count[1] %d\n",mpi_rank,start[0],start[1],count[0],count[1]);*/
/*            printf("mpi_rank %d,startOut[0] %d,startOut[1] %d,countOut[0] %d countOut[1] %d\n",mpi_rank,startOut[0],startOut[1],countOut[0],countOut[1]);*/
			{
			int* data1=(int *)malloc(sizeof(int)*dataEnd);
			int* data2=(int *)malloc(sizeof(int)*dataEnd);
			int* dataOut=(int *)malloc(sizeof(int)*dataEnd);
			if((res=nc_get_vara_int(ncid1,vlid1,start,countOut,data1)))
				BAIL(res);
			if((res=nc_get_vara_int(ncid2,vlid2,start,countOut,data2)))
				BAIL(res);
			for(j=0;j<dataEnd;++j){
				binary_op_int(op,&data1[j],&data2[j],&dataOut[j]);
			}
			if((res=nc_put_vara_int(ncidout,vlidout,startOut,countOut,(int *)dataOut)))
				BAIL(res);
			free(data1);
			free(data2);
			free(dataOut);
			}
			break;
		case NC_FLOAT:
			{
			float* data1=(float *)malloc(sizeof(float)*dataEnd);
			float* data2=(float *)malloc(sizeof(float)*dataEnd);
			float* dataOut=(float *)malloc(sizeof(float)*dataEnd);
			if((res=nc_get_vara_float(ncid1,vlid1,start,countOut,data1)))
				BAIL(res);
			if((res=nc_get_vara_float(ncid2,vlid2,start,countOut,data2)))
				BAIL(res);
			for(j=0;j<dataEnd;++j){
				binary_op_float(op,&data1[j],&data2[j],&dataOut[j]);
			}
			if((res=nc_put_vara_float(ncidout,vlidout,startOut,countOut,(float *)dataOut)))
				BAIL(res);
			free(data1);
			free(data2);
			free(dataOut);
			}
			break;
		case NC_DOUBLE:
			{
			double* data1=(double *)malloc(sizeof(double)*dataEnd);
			double* data2=(double *)malloc(sizeof(double)*dataEnd);
			double* dataOut=(double *)malloc(sizeof(double)*dataEnd);
			if((res=nc_get_vara_double(ncid1,vlid1,start,countOut,data1)))
				BAIL(res);
			if((res=nc_get_vara_double(ncid2,vlid2,start,countOut,data2)))
				BAIL(res);
			for(j=0;j<dataEnd;++j){
				binary_op_double(op,&data1[j],&data2[j],&dataOut[j]);
			}
			if((res=nc_put_vara_double(ncidout,vlidout,startOut,countOut,(double *)dataOut)))
				BAIL(res);
			free(data1);
			free(data2);
			free(dataOut);
			}
			break;
		default:
			printf("Unknown data type\n");
			}
	free(start);
	free(startOut);
	free(shapeOut);
	free(countOut);
/*    free(poses);*/
	return 0;
}
int extract_binary_selector(BIN_OP op,int mpi_rank,int mpi_size, int ncid1,int vlid1,int ncid2,int vlid2,int ncidout,int vlidout,int ndims,nc_type vtype,size_t *shape,size_t *begins,size_t *ends, ptrdiff_t *strides,size_t preLen,size_t *outLen){
	int i;
	size_t s=1;
	for(i=0;i<ndims;++i){
		s+=strides[i];
	}
	if(s==ndims){
		printf("use extract_binary_1D_all\n");
		extract_binary_1D_all(op,mpi_rank,mpi_size,ncid1,vlid1,ncid2,vlid2,ncidout,vlidout,ndims,vtype,shape,begins,ends,preLen,outLen);
	}
	else if(mpi_size<=(ends[0]-begins[0])/strides[0]+1||s/ndims<2){
		printf("use extract_binary_1D\n");
		extract_binary_1D(op,mpi_rank,mpi_size,ncid1,vlid1,ncid2,vlid2,ncidout,vlidout,ndims,vtype,shape,begins,ends,strides,preLen,outLen);
	}else{
		printf("use extract_binary_single\n");
		extract_binary_single(op,mpi_rank,mpi_size,ncid1,vlid1,ncid2,vlid2,ncidout,vlidout,ndims,vtype,shape,begins,ends,strides,preLen,outLen);
	}
	return 0;
}

