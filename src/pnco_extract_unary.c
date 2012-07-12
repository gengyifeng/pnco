#include "netcdf.h"
#include <mpi.h>
#include <assert.h>
#include "hdf5.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "pnco_init.h"
#include "pnco_extract_util.h"
#include "pnco_parser.h"

#define BAIL(e) do{printf("ERROR:: file: %s, line: %d, func: %s, code: %s.\n",__FILE__,__LINE__,__FUNCTION__, nc_strerror(1));return e;} while(0)

extern GLOBAL_ATTR global_attr;
/* extract a variable from a netcdf file ant write it to another netcdf parallely, the function
 * will create new netcdf file if the output file does not exist. If the output file exists and the
 * output variable doesn't exist in the output file, the function creates a new variable or aborts.*/
int extract_unary(char *var,char *varOut,char **outDims,int outDimNum,size_t *begins,int beginNum,size_t *ends,int endNum,ptrdiff_t *strides,int strideNum,char *fileIn,char* fileOut){
	int mpi_size,mpi_rank;
	MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);

	MPI_Comm comm=MPI_COMM_WORLD;
	MPI_Info info=MPI_INFO_NULL;
	
	int ncid,vlid,ndims,dimids[NC_MAX_VAR_DIMS],nattrs;
	int ncidout,vlidout;
	int i,res;
	char varName[NC_MAX_NAME+1];
	nc_type vtype;	

	/* 
		open input file
	*/
	if((res=nc_open_par(fileIn,NC_MPIIO,comm,info,&ncid))){
		printf("Cannot open input file %s\n",fileIn);
		BAIL(res);
	}
	if((res=nc_inq_varid(ncid,var,&vlid))){
		printf("Cannot query variable %s in input file %s\n",var,fileIn);
		BAIL(res);	
	}
	if((res=nc_inq_var(ncid,vlid,varName,&vtype,&ndims,dimids,&nattrs))){
		printf("Cannot query dimension info of varialbe %s in input file %s\n",var,fileIn);
		BAIL(res);
	}


	char **dimsName=(char **)malloc(sizeof(char *)*ndims);
	size_t *shape=(size_t *)malloc(sizeof(size_t)*ndims);


	for(i=0;i<ndims;++i){
		dimsName[i]=(char *)malloc(NC_MAX_NAME+1);
		nc_inq_dim(ncid,dimids[i],dimsName[i],&shape[i]);//get dimsName and dimsLen
	}

	int is_strides_null=0;
	int is_begins_null=0;
	int is_ends_null=0;
	
	/* initialize begins ends and strides*/
	init_begins_ends_strides(ndims,&begins,&beginNum,&is_begins_null,&ends,&endNum,&is_ends_null,&strides,&strideNum,&is_strides_null,shape);
	
	/* deal with -d and -r arguments */
	dim_opts_handler(&global_attr,dimids,ncid,ndims,begins,ends,strides);
	
	/* check begins ends and strides */
	if(-1==check_begins_ends_strides(ndims,begins,beginNum,ends,endNum,strides,strideNum,shape,dimsName)){return -1;}


	/* set independent parallel access for input file */
	if((res=nc_var_par_access(ncid,vlid,NC_INDEPENDENT))){
		printf("Cannot set parallel access for variable %s in input file %s\n",var,fileIn);
		BAIL(res);
	}
	/*
		create output file and define dims and var
	*/
	int *dimidsOut=(int *)malloc(sizeof(int)*ndims);
	if(-1==	open_output_file_and_define_var(fileOut,&ncidout,&vlidout,dimidsOut,global_attr.append,ndims,vtype,dimsName,begins,ends,strides,shape,var,varOut,outDims,outDimNum)){return -1;}
	
	
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
	extract_unary_selector(mpi_rank,mpi_size,ncid,vlid,ncidout,vlidout,ndims,vtype,shape,begins,ends,strides,0L,NULL);// select reasonable parallel funtion

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

	/*close file*/
	if((res=nc_close(ncid))){
		printf("Close input file %s failed\n",fileIn);
		BAIL(res);
	}
	if((res=nc_close(ncidout))){
		printf("Close output file %s failed\n",fileOut);
		BAIL(res);
	}
	return 0;
}


/* This function treats multi-dimension data as one dimension array, then distributes the array range
 * to each proccess according to the mpi_rank. Due to the underline pararllel netcdf library, only
 * nc_get_var1_* and nc_put_var2_* can be used to read or write single element. This function works
 * but has low performance, so it should be used properly when extract_unary_1D or
 * extract_unary_1D_all does not work parallelly well.*/

int extract_unary_single(int mpi_rank,int mpi_size, int ncid,int vlid,int ncidout,int vlidout,int ndims,nc_type vtype,size_t *shape,size_t *begins,size_t *ends, ptrdiff_t *strides,size_t preLen,size_t *outLen){
	
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
	void *data=malloc(sizeof(double));
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
			if((res=nc_get_var1_uchar(ncid,vlid,start,data)))
				BAIL(res);
			if((res=nc_put_var1_uchar(ncidout,vlidout,startOut,(unsigned char *)data)))
				BAIL(res);
			break;
		case NC_CHAR:
			if((res=nc_get_var1_schar(ncid,vlid,start,(signed char *)data)))
				BAIL(res);
			if((res=nc_put_var1_schar(ncidout,vlidout,startOut,(signed char *)data)))
				BAIL(res);
			break;
		case NC_SHORT:
			if((res=nc_get_var1_short(ncid,vlid,start,data)))
				BAIL(res);
			if((res=nc_put_var1_short(ncidout,vlidout,startOut,(short *)data)))
				BAIL(res);
			break;
		case NC_INT:
			if((res=nc_get_var1_int(ncid,vlid,start,(int *)data)))
				BAIL(res);
			if((res=nc_put_var1_int(ncidout,vlidout,startOut,(int *)data)))
				BAIL(res);
			break;
		case NC_FLOAT:
			if((res=nc_get_var1_float(ncid,vlid,start,data)))
				BAIL(res);
			if((res=nc_put_var1_float(ncidout,vlidout,startOut,(float *)data)))
				BAIL(res);
			break;
		case NC_DOUBLE:
			if((res=nc_get_var1_double(ncid,vlid,start,data)))
				BAIL(res);
			if((res=nc_put_var1_double(ncidout,vlidout,startOut,(double *)data)))
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
	free(data);
	return 0;
}

/* This function use the  first dimension to parallel. The difference between extract_unary_1D and
 * extract unary_1D_all is that extract_unary_1D in each iteration reads and writes a row like
 * v[1][?][?].
 * but extract_unary_1D_all reads and writes all needed data only once, there is no iteration.
 * */

int extract_unary_1D(int mpi_rank,int mpi_size, int ncid,int vlid,int ncidout,int vlidout,int ndims,nc_type vtype,size_t *shape,size_t *begins,size_t *ends, ptrdiff_t *strides,size_t preLen,size_t *outLen){
	
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
	void* data=(void*)malloc(sizeof(double)*len0);
	void* dataOut=(void*)malloc(sizeof(double)*dataEnd);
	size_t* poses=(size_t*)malloc(sizeof(size_t)*dataEnd);
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
			if((res=nc_get_vara_uchar(ncid,vlid,start,count,(unsigned char *)data)))
				BAIL(res);
			if(len0==dataEnd){
				if((res=nc_put_vara_uchar(ncidout,vlidout,startOut,countOut,(unsigned char *)data)))
					BAIL(res);
			}else{
				for(j=0;j<dataEnd;++j){
					((unsigned char *)dataOut)[j]=((unsigned char *)data)[poses[j]];
				}
				if((res=nc_put_vara_uchar(ncidout,vlidout,startOut,countOut,(unsigned char *)dataOut)))
					BAIL(res);
			}
			break;
		case NC_CHAR:
			if((res=nc_get_vara_schar(ncid,vlid,start,count,(signed char *)data)))
				BAIL(res);
			if(len0==dataEnd){
				if((res=nc_put_vara_schar(ncidout,vlidout,startOut,countOut,(signed char *)data)))
					BAIL(res);
			}else{
				for(j=0;j<dataEnd;++j){
					((signed char *)dataOut)[j]=((signed char *)data)[poses[j]];
				}
				if((res=nc_put_vara_schar(ncidout,vlidout,startOut,countOut,(signed char *)dataOut)))
					BAIL(res);
			}
			break;
		case NC_SHORT:
			if((res=nc_get_vara_short(ncid,vlid,start,count,data)))
				BAIL(res);
			if(len0==dataEnd){
				if((res=nc_put_vara_short(ncidout,vlidout,startOut,countOut,(short *)data)))
					BAIL(res);
			}else{
				for(j=0;j<dataEnd;++j){
					((short *)dataOut)[j]=((short *)data)[poses[j]];
				}
				if((res=nc_put_vara_short(ncidout,vlidout,startOut,countOut,(short *)dataOut)))
					BAIL(res);
			}
			break;
		case NC_INT:
			if((res=nc_get_vara_int(ncid,vlid,start,count,(int *)data)))
				BAIL(res);
			if(len0==dataEnd){
				if((res=nc_put_vara_int(ncidout,vlidout,startOut,countOut,(int *)data)))
					BAIL(res);
			}else{
				for(j=0;j<dataEnd;++j){
					((int *)dataOut)[j]=((int *)data)[poses[j]];
				}
				if((res=nc_put_vara_int(ncidout,vlidout,startOut,countOut,(int *)dataOut)))
					BAIL(res);
			}

			break;
		case NC_FLOAT:
			if((res=nc_get_vara_float(ncid,vlid,start,count,data)))
				BAIL(res);
			if(len0==dataEnd){
				if((res=nc_put_vara_float(ncidout,vlidout,startOut,countOut,(float *)data)))
					BAIL(res);
			}else{
				for(j=0;j<dataEnd;++j){
					((float *)dataOut)[j]=((float *)data)[poses[j]];
				}
				if((res=nc_put_vara_float(ncidout,vlidout,startOut,countOut,(float *)dataOut)))
					BAIL(res);
			}
			break;
		case NC_DOUBLE:
			if((res=nc_get_vara_double(ncid,vlid,start,count,(double *)data)))
				BAIL(res);
			if(len0==dataEnd){
				if((res=nc_put_vara_double(ncidout,vlidout,startOut,countOut,(double *)data)))
					BAIL(res);
			}else{
				for(j=0;j<dataEnd;++j){
					((double *)dataOut)[j]=((double *)data)[poses[j]];
				}
				if((res=nc_put_vara_double(ncidout,vlidout,startOut,countOut,(double *)dataOut)))
					BAIL(res);
			}
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
	free(data);
	free(dataOut);
	free(poses);
	return 0;
}


/* This function use the  first dimension to parallel and reads and writes all needed data only once, there is no iteration. Stride is not supported */
int extract_unary_1D_all(int mpi_rank,int mpi_size, int ncid,int vlid,int ncidout,int vlidout,int ndims,nc_type vtype,size_t *shape,size_t *begins,size_t *ends,size_t preLen,size_t *outLen){
	
	int i,res;
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
			len0*=ends[i]-begins[i]+1;
		}
	}
	if(outLen!=NULL)
		*outLen=lenOut;
	if(shapeOut[0]>=mpi_size){
		startOut0=mpi_rank*(shapeOut[0]/mpi_size);
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
			unsigned char* dataOut=(unsigned char *)malloc(sizeof(unsigned char)*dataEnd);
			if((res=nc_get_vara_uchar(ncid,vlid,start,countOut,dataOut)))
				BAIL(res);
			if((res=nc_put_vara_uchar(ncidout,vlidout,startOut,countOut,(unsigned char *)dataOut)))
				BAIL(res);
			free(dataOut);
			}
			break;
		case NC_CHAR:
			{
			char* dataOut=(char *)malloc(sizeof(char)*dataEnd);
			if((res=nc_get_vara_schar(ncid,vlid,start,countOut,(signed char *)dataOut)))
				BAIL(res);
			if((res=nc_put_vara_schar(ncidout,vlidout,startOut,countOut,(signed char *)dataOut)))
				BAIL(res);
			free(dataOut);
			}
			break;
		case NC_SHORT:
			{
			short *dataOut=(short *)malloc(sizeof(short)*dataEnd);
			if((res=nc_get_vara_short(ncid,vlid,start,countOut,(short *)dataOut)))
				BAIL(res);
			if((res=nc_put_vara_short(ncidout,vlidout,startOut,countOut,(short *)dataOut)))
				BAIL(res);
			free(dataOut);
			}
			break;
		case NC_INT:
			{
			int * dataOut=(int *)malloc(sizeof(int)*dataEnd);
			if((res=nc_get_vara_int(ncid,vlid,start,countOut,(int *)dataOut)))
				BAIL(res);
			if((res=nc_put_vara_int(ncidout,vlidout,startOut,countOut,(int *)dataOut)))
				BAIL(res);
			free(dataOut);
			}
			break;
		case NC_FLOAT:
			{
			float * dataOut=(float *)malloc(sizeof(float)*dataEnd);
			if((res=nc_get_vara_float(ncid,vlid,start,countOut,(float *)dataOut)))
				BAIL(res);
			if((res=nc_put_vara_float(ncidout,vlidout,startOut,countOut,(float *)dataOut)))
				BAIL(res);
			free(dataOut);
			}
			break;
		case NC_DOUBLE:
			{
			double* dataOut=(double *)malloc(sizeof(double)*dataEnd);
			if((res=nc_get_vara_double(ncid,vlid,start,countOut,dataOut)))
				BAIL(res);
			if((res=nc_put_vara_double(ncidout,vlidout,startOut,countOut,(double *)dataOut)))
				BAIL(res);
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
	return 0;
}
/* Pre-calculate the needed output data position in the input row. this is used in
 * extract_unary_1D*/
int transfer_pos(size_t *poses,int ndims,ptrdiff_t *strides,size_t outlen,size_t *dividerOut,size_t *divider){
	size_t i,j,rem,pos;
	for(i=0;i<outlen;++i){
		rem=i;
		pos=0;
		for(j=0;j<ndims;j++){
			pos+=(rem/dividerOut[j+1])*divider[j+1]*strides[j+1];
			rem=rem%dividerOut[j+1];
		}
		poses[i]=pos;
	}
	return 0;
}
/* This function selects the right parallel strategy according to the input arguments. Current
 * implementation is not well tested.*/
int extract_unary_selector(int mpi_rank,int mpi_size, int ncid,int vlid,int ncidout,int vlidout,int ndims,nc_type vtype,size_t *shape,size_t *begins,size_t *ends, ptrdiff_t *strides,size_t preLen,size_t *outLen){
	int i;
	size_t s=0;
	for(i=0;i<ndims;++i){
		s+=strides[i];
	}
	if(s==ndims){
		printf("use extract_unary_1D_all\n");
		extract_unary_1D_all(mpi_rank,mpi_size,ncid,vlid,ncidout,vlidout,ndims,vtype,shape,begins,ends,preLen,outLen);
	}
	else if(mpi_size<=(ends[0]-begins[0])/strides[0]+1||s/ndims<2){
		printf("use extract_unary_1D\n");
		extract_unary_1D(mpi_rank,mpi_size,ncid,vlid,ncidout,vlidout,ndims,vtype,shape,begins,ends,strides,preLen,outLen);
	}else{
		printf("use extract_unary_single\n");
		extract_unary_single(mpi_rank,mpi_size,ncid,vlid,ncidout,vlidout,ndims,vtype,shape,begins,ends,strides,preLen,outLen);
	}
	return 0;
}
