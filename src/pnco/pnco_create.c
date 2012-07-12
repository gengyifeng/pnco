#include "netcdf.h"
#include <mpi.h>
#include <assert.h>
#include <hdf5.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "pnco.h"
#define BAIL(e) do{ \
printf("ERROR:: file: %s, line: %d, func: %s, code: %s.\n",__FILE__,__LINE__,__FUNCTION__, nc_strerror(e)); \
return e;\
} while(0)

int create_sample_2D(char *filename){
	DIMS dims;
	int ndims=2;
	size_t arr[]={24,12};
	char *dim_names[]={"d1","d2"};
	init_dims(&dims,ndims,arr,dim_names);
	
	int mpi_size,mpi_rank;
	MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
	int block_size=dims.shape[0]/mpi_size*dims.shape[1];
	printf("mpi_size: %d mpi_rand: %d block_size: %d\n",mpi_size,mpi_rank,block_size);
	MPI_Comm comm=MPI_COMM_WORLD;
	MPI_Info info=MPI_INFO_NULL;
	int ncid,vlid,dimids[2];
	int d1_id,d2_id;
	
	int res;	
	if((res=nc_create_par(filename,NC_NETCDF4|NC_MPIIO,comm,info,&ncid)))
		BAIL(res);
/*    if((res=nc_def_dim(ncid,"d1",dims.shape[0],&dimids[0])))*/
/*        BAIL(res);*/
	if((res=nc_def_dim(ncid,"d1",NC_UNLIMITED,&dimids[0])))
		BAIL(res);
	if((res=nc_def_dim(ncid,"d2",dims.shape[1],&dimids[1])))
		BAIL(res);
/*    if((res=nc_def_dim(ncid,"d2",NC_UNLIMITED,&dimids[1])))*/
/*        BAIL(res);*/
	if((res=nc_def_var(ncid,"d1",NC_INT,1,&dimids[0],&d1_id)))
		BAIL(res);
	if((res=nc_def_var(ncid,"d2",NC_FLOAT,1,&dimids[1],&d2_id)))
		BAIL(res);
	if((res=nc_def_var(ncid,"v1",NC_DOUBLE,ndims,dimids,&vlid)))
		BAIL(res);
	if((res=nc_enddef(ncid)))
		BAIL(res);
	
	size_t start[2],count[2];
	start[0]=mpi_rank*(dims.shape[0]/mpi_size);
	start[1]=0;
	count[0]=dims.shape[0]/mpi_size;
	count[1]=dims.shape[1];
	size_t start1,start2,count1,count2;
	
	start1=mpi_rank*(dims.shape[0]/mpi_size);
	start2=mpi_rank*(dims.shape[1]/mpi_size);
	count1=dims.shape[0]/mpi_size;
	count2=dims.shape[1]/mpi_size;
	printf("mpi_rank=%d start[0]=%d start[1]=%d count[0]=%d count[1]=%d\n",mpi_rank, start[0],start[1],count[0],count[1]);

	printf("mpi_rank*block_size=%d (mpi_rank+1)*block_size-1=%d\n",mpi_rank*block_size,(mpi_rank+1)*block_size);
	int i;	
	int *data1=(int *)malloc(sizeof(int)*mpi_size*(dims.shape[0]/mpi_size));
	float*data2=(float *)malloc(sizeof(float)*mpi_size*(dims.shape[1]/mpi_size));

	double *data=(double *)malloc(sizeof(double)*mpi_size*block_size);
	for(i=mpi_rank*block_size;i<(mpi_rank+1)*block_size;++i)
		data[i]=287-i;
	for(i=mpi_rank*(dims.shape[0]/mpi_size);i<(mpi_rank+1)*(dims.shape[0]/mpi_size);++i)
		data1[i]=i*2;
	for(i=mpi_rank*(dims.shape[1]/mpi_size);i<(mpi_rank+1)*(dims.shape[1]/mpi_size);++i)
		data2[i]=(11-i)*2;
	if((res=nc_var_par_access(ncid,d1_id,NC_INDEPENDENT)))
		BAIL(res);
	if((res=nc_var_par_access(ncid,d2_id,NC_INDEPENDENT)))
		BAIL(res);
	if((res=nc_var_par_access(ncid,vlid,NC_INDEPENDENT)))
		BAIL(res);
	if((res=nc_put_vara_int(ncid,d1_id,&start1,&count1,&data1[mpi_rank*(dims.shape[0]/mpi_size)])))
		BAIL(res);
	if((res=nc_put_vara_float(ncid,d2_id,&start2,&count2,&data2[mpi_rank*(dims.shape[1]/mpi_size)])))
		BAIL(res);
	if((res=nc_put_vara_double(ncid,vlid,start,count,&data[mpi_rank*block_size])))
		BAIL(res);
	if((res=nc_close(ncid)))
		BAIL(res);
	return 0;

}
int create_sample(char *filename){
	DIMS dims;
	int ndims=3;
	size_t arr[]={24,12,6};
	char *dim_names[]={"d1","d2","d3"};
	init_dims(&dims,ndims,arr,dim_names);
	
	int mpi_size,mpi_rank;
	MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
	int block_size=(dims.shape[0]/mpi_size)*dims.shape[1]*dims.shape[2];
	printf("mpi_size: %d mpi_rand: %d block_size: %d\n",mpi_size,mpi_rank,block_size);
	MPI_Comm comm=MPI_COMM_WORLD;
	MPI_Info info=MPI_INFO_NULL;
	int ncid,vlid,dimids[3];
	int d1_id,d2_id,d3_id;
	
	int res;	
	if((res=nc_create_par(filename,NC_NETCDF4|NC_MPIIO,comm,info,&ncid)))
		BAIL(res);
/*    if((res=nc_def_dim(ncid,"d1",dims.shape[0],&dimids[0])))*/
/*        BAIL(res);*/
	if((res=nc_def_dim(ncid,"d1",NC_UNLIMITED,&dimids[0])))
		BAIL(res);
	if((res=nc_def_dim(ncid,"d2",dims.shape[1],&dimids[1])))
		BAIL(res);
	if((res=nc_def_dim(ncid,"d3",dims.shape[2],&dimids[2])))
		BAIL(res);
/*    if((res=nc_def_dim(ncid,"d2",NC_UNLIMITED,&dimids[1])))*/
/*        BAIL(res);*/
	if((res=nc_def_var(ncid,"d1",NC_INT,1,&dimids[0],&d1_id)))
		BAIL(res);
	if((res=nc_def_var(ncid,"d2",NC_FLOAT,1,&dimids[1],&d2_id)))
		BAIL(res);
	if((res=nc_def_var(ncid,"d3",NC_FLOAT,1,&dimids[2],&d3_id)))
		BAIL(res);
	if((res=nc_def_var(ncid,"v1",NC_DOUBLE,ndims,dimids,&vlid)))
		BAIL(res);
	if((res=nc_enddef(ncid)))
		BAIL(res);
	
	size_t start[3],count[3];
	start[0]=mpi_rank*(dims.shape[0]/mpi_size);
	start[1]=0;
	start[2]=0;
	count[0]=dims.shape[0]/mpi_size;
	count[1]=dims.shape[1];
	count[2]=dims.shape[2];
	size_t start1,start2,start3,count1,count2,count3;
	
	start1=mpi_rank*(dims.shape[0]/mpi_size);
	start2=mpi_rank*(dims.shape[1]/mpi_size);
	start3=mpi_rank*(dims.shape[2]/mpi_size);
	count1=dims.shape[0]/mpi_size;
	count2=dims.shape[1]/mpi_size;
	count3=dims.shape[2]/mpi_size;
	printf("mpi_rank=%d start[0]=%d start[1]=%d count[0]=%d count[1]=%d\n",mpi_rank, start[0],start[1],count[0],count[1]);

	printf("mpi_rank*block_size=%d (mpi_rank+1)*block_size-1=%d\n",mpi_rank*block_size,(mpi_rank+1)*block_size);
	int i;	
	int *data1=(int *)malloc(sizeof(int)*mpi_size*(dims.shape[0]/mpi_size));
	float*data2=(float *)malloc(sizeof(float)*mpi_size*(dims.shape[1]/mpi_size));
	float*data3=(float *)malloc(sizeof(float)*mpi_size*(dims.shape[2]/mpi_size));

	double *data=(double *)malloc(sizeof(double)*mpi_size*block_size);
	for(i=mpi_rank*block_size;i<(mpi_rank+1)*block_size;++i)
		data[i]=287-i;
	for(i=mpi_rank*(dims.shape[0]/mpi_size);i<(mpi_rank+1)*(dims.shape[0]/mpi_size);++i)
		data1[i]=i*2;
	for(i=mpi_rank*(dims.shape[1]/mpi_size);i<(mpi_rank+1)*(dims.shape[1]/mpi_size);++i)
		data2[i]=(11-i)*2;
	for(i=mpi_rank*(dims.shape[2]/mpi_size);i<(mpi_rank+1)*(dims.shape[2]/mpi_size);++i)
		data3[i]=i;
	if((res=nc_var_par_access(ncid,d1_id,NC_INDEPENDENT)))
		BAIL(res);
	if((res=nc_var_par_access(ncid,d2_id,NC_INDEPENDENT)))
		BAIL(res);
	if((res=nc_var_par_access(ncid,d3_id,NC_INDEPENDENT)))
		BAIL(res);
	if((res=nc_var_par_access(ncid,vlid,NC_INDEPENDENT)))
		BAIL(res);
	if((res=nc_put_vara_int(ncid,d1_id,&start1,&count1,&data1[mpi_rank*(dims.shape[0]/mpi_size)])))
		BAIL(res);
	if((res=nc_put_vara_float(ncid,d2_id,&start2,&count2,&data2[mpi_rank*(dims.shape[1]/mpi_size)])))
		BAIL(res);
	if((res=nc_put_vara_float(ncid,d3_id,&start3,&count3,&data3[mpi_rank*(dims.shape[2]/mpi_size)])))
		BAIL(res);
	if((res=nc_put_vara_double(ncid,vlid,start,count,&data[mpi_rank*block_size])))
		BAIL(res);
	if((res=nc_close(ncid)))
		BAIL(res);
	return 0;

}
