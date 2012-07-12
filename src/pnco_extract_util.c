#include "netcdf.h"
#include <mpi.h>
#include <assert.h>
#include "hdf5.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "pnco_init.h"
#include "pnco_parser.h"
#include "pnco_extract_util.h"


/* get divider for multi-dimension data*/
int getDivider(int ndims,const size_t *dims,size_t *divider){
	int i;
	divider[ndims-1]=1;
	if(ndims>1){
		divider[ndims-2]=dims[ndims-1];
	}
	for(i=ndims-3;i>=0;--i){
		divider[i]=dims[ndims-2-i]*divider[i+1];
	}
	return 0;
}
/* get multi-dimension start according to the one-dimension array index*/
int getStart(int ndims,size_t index,const size_t *divider,size_t *start){
	int i;
	size_t remain=index;
	for(i=0;i<ndims-1;++i){
		start[i]=remain/divider[i];
		remain=remain%divider[i];
	}
	start[ndims-1]=remain;
	return 0;
}
/* get one-dimension index according to the multi-dimension start*/
int getIndex(int ndims,size_t *start,const size_t *divider,size_t *index){
	int i;
	size_t sum=0;
	for(i=0;i<ndims-1;++i){
		sum+=start[i]*divider[i];
	}
	sum+=start[ndims-1];
	*index=sum;
	return 0;
}
/* binary_operator */
int binary_op_uchar(BIN_OP op,unsigned char *data1,unsigned char*data2,unsigned char *dataOut){
	switch(op){
		case ADD:
			*dataOut=(*data1)+(*data2);
			break;
		case SUB:
			*dataOut=(*data1)-(*data2);
			break;
		case DIV:
			if(*data2!=0)
				*dataOut=(*data1)/(*data2);
			else{
				printf("divider is zero!!!\n");
				*dataOut=0;
			}
			break;
		case MUL:
			*dataOut=(*data1)*(*data2);
			break;
		case MOD:
			if(*data2!=0)
				*dataOut=(*data1)%(*data2);
			else{
				printf("divider is zero!!!\n");
				*dataOut=0;
			}
			break;
		case AVG:
			*dataOut=((*data1)+(*data2))/2;
			break;
		default:
			printf("Unknown operation\n");
			*dataOut=0;
	}
	return 0;
}
int binary_op_schar(BIN_OP op,signed char *data1,signed char*data2,signed char *dataOut){
	switch(op){
		case ADD:
			*dataOut=(*data1)+(*data2);
			break;
		case SUB:
			*dataOut=(*data1)-(*data2);
			break;
		case DIV:
			if(*data2!=0)
				*dataOut=(*data1)/(*data2);
			else{
				printf("divider is zero!!!\n");
				*dataOut=0;
			}
			break;
		case MUL:
			*dataOut=(*data1)*(*data2);
			break;
		case MOD:
			if(*data2!=0)
				*dataOut=(*data1)%(*data2);
			else{
				printf("divider is zero!!!\n");
				*dataOut=0;
			}
			break;
		case AVG:
			*dataOut=((*data1)+(*data2))/2;
			break;
		default:
			printf("Unknown operation\n");
			*dataOut=0;
	}
	return 0;
}
int binary_op_short(BIN_OP op,short *data1,short*data2,short *dataOut){
	switch(op){
		case ADD:
			*dataOut=(*data1)+(*data2);
			break;
		case SUB:
			*dataOut=(*data1)-(*data2);
			break;
		case DIV:
			if(*data2!=0)
				*dataOut=(*data1)/(*data2);
			else{
				printf("divider is zero!!!\n");
				*dataOut=0;
			}
			break;
		case MUL:
			*dataOut=(*data1)*(*data2);
			break;
		case MOD:
			if(*data2!=0)
				*dataOut=(*data1)%(*data2);
			else{
				printf("divider is zero!!!\n");
				*dataOut=0;
			}
			break;
		case AVG:
			*dataOut=((*data1)+(*data2))/2;
			break;
		default:
			printf("Unknown operation\n");
			*dataOut=0;
	}
	return 0;
}
int binary_op_int(BIN_OP op,int *data1,int*data2,int *dataOut){
	switch(op){
		case ADD:
			*dataOut=(*data1)+(*data2);
			break;
		case SUB:
			*dataOut=(*data1)-(*data2);
			break;
		case DIV:
			if(*data2!=0)
				*dataOut=(*data1)/(*data2);
			else{
				printf("divider is zero!!!\n");
				*dataOut=0;
			}
			break;
		case MUL:
			*dataOut=(*data1)*(*data2);
			break;
		case MOD:
			if(*data2!=0)
				*dataOut=(*data1)%(*data2);
			else{
				printf("divider is zero!!!\n");
				*dataOut=0;
			}
			break;
		case AVG:
			*dataOut=((*data1)+(*data2))/2;
			break;
		default:
			printf("Unknown operation\n");
			*dataOut=0;
	}
	return 0;
}
int binary_op_long(BIN_OP op,long *data1,long *data2,long *dataOut){
	switch(op){
		case ADD:
			*dataOut=(*data1)+(*data2);
			break;
		case SUB:
			*dataOut=(*data1)-(*data2);
			break;
		case DIV:
			if(*data2!=0)
				*dataOut=(*data1)/(*data2);
			else{
				printf("divider is zero!!!\n");
				*dataOut=0;
			}
			break;
		case MUL:
			*dataOut=(*data1)*(*data2);
			break;
		case MOD:
			if(*data2!=0)
				*dataOut=(*data1)%(*data2);
			else{
				printf("divider is zero!!!\n");
				*dataOut=0;
			}
			break;
		case AVG:
			*dataOut=((*data1)+(*data2))/2;
			break;
		default:
			printf("Unknown operation\n");
			*dataOut=0;
	}
	return 0;
}
int binary_op_float(BIN_OP op,float *data1,float *data2,float *dataOut){
	switch(op){
		case ADD:
			*dataOut=(*data1)+(*data2);
			break;
		case SUB:
			*dataOut=(*data1)-(*data2);
			break;
		case DIV:
			if(*data2!=0)
				*dataOut=(*data1)/(*data2);
			else{
				printf("divider is zero!!!\n");
				*dataOut=0;
			}
			break;
		case MUL:
			*dataOut=(*data1)*(*data2);
			break;
		case MOD:
			printf("MOD is not available when data type is float!\n");
			*dataOut=0;
			break;
		case AVG:
			*dataOut=((*data1)+(*data2))/2;
			break;
		default:
			printf("Unknown operation!\n");
			*dataOut=0;
	}
	return 0;
}
int binary_op_double(BIN_OP op,double *data1,double *data2,double *dataOut){
	switch(op){
		case ADD:
			*dataOut=(*data1)+(*data2);
			break;
		case SUB:
			*dataOut=(*data1)-(*data2);
			break;
		case DIV:
			if(*data2!=0)
				*dataOut=(*data1)/(*data2);
			else{
				printf("divider is zero!!!\n");
				*dataOut=0;
			}
			break;
		case MUL:
			*dataOut=(*data1)*(*data2);
			break;
		case MOD:
			printf("MOD is not available when data type is double!\n");
			*dataOut=0;
			break;
		case AVG:
			*dataOut=((*data1)+(*data2))/2;
			break;
		default:
			printf("Unknown operation!\n");
			*dataOut=0;
	}
	return 0;
}
/* copy variable attributes from one variable to another*/
int copy_attrs(int ncid, int vlid, int ncidout, int vlidout, int nattrs){
	int i;
	char attrName[NC_MAX_NAME+1];
	for(i=0;i<nattrs;++i){
		nc_inq_attname(ncid,vlid,i,attrName);
		nc_copy_att(ncid,vlid,attrName,ncidout,vlidout);
	}
	return 0;
}
/* add pnco command to global attribute*/
int add_cmd_attr(int ncidout,int vlidout,int argc, char **argv,int mpi_size){
	char number[12];
	//itoa(mpi_size,number,10);
	sprintf(number,"%d",mpi_size);
	char cmd[512]="mpiexec -n ";
	strcat(cmd,number);
	int i,res;
	for(i=0;i<argc;++i){
		strcat(cmd," ");
		strcat(cmd,argv[i]);
	}
	cmd[511]='\0';
	i=0;
	do{
		++i;
		char name[50]="pnco_command_history_";
		sprintf(number,"%d",i);
		strcat(name,number);
		res=nc_put_att_text(ncidout,NC_GLOBAL,name,strlen(cmd),cmd);
	}while(res!=0);
	return 0;
}
/* find the position that is bigger or equal to the minimal value in '-d' argument in the ascending dimension variable*/
int find_var_min(int ncid, int vlid,int ndims,size_t *divider,size_t dlen,double des,nc_type vtype,size_t *pos){
	size_t *start=(size_t *)malloc(sizeof(size_t)*ndims);	
	size_t begin=0;
	size_t mid;
	size_t end=dlen-1;
	double min,val,max;
	getStart(ndims,begin,divider,start);
	get_var1_double(ncid,vlid,start,&min,vtype);
	getStart(ndims,end,divider,start);
	get_var1_double(ncid,vlid,start,&max,vtype);
	if(des>max){
		*pos=dlen;
		free(start);
		return -2;
	}
	if(des<min){
		*pos=0;
		free(start);
		return -2;
	}
	while(begin<=end){
		mid=begin+(end-begin)/2;
		getStart(ndims,mid,divider,start);
		get_var1_double(ncid,vlid,start,&val,vtype);
		if(val==des){
			while(val==des&&mid>0){
				mid--;
				getStart(ndims,mid,divider,start);
				get_var1_double(ncid,vlid,start,&val,vtype);
			}
			if(mid!=0)
				mid++;
			*pos=mid;
			free(start);
			return 0;
		}
		if(begin==end)
			break;
		if(val>des){
			end=mid-1;
		}else{
			begin=mid+1;
		}
	}
	if(mid==dlen-1){
		*pos=mid;
	}else{
		if(val>des){
			*pos=mid;
		}else{
			*pos=++mid;
		}
	}
	free(start);
	return -1;
}
/* find the position that is smaller or equal to the max value in '-d' argument in the ascending dimension variable*/
int find_var_max(int ncid, int vlid,int ndims,size_t *divider,size_t dlen,double des,nc_type vtype,size_t *pos){
	size_t *start=(size_t *)malloc(sizeof(size_t)*ndims);	
	size_t begin=0;
	size_t mid;
	size_t end=dlen-1;
	double val,min,max;
	getStart(ndims,begin,divider,start);
	get_var1_double(ncid,vlid,start,&min,vtype);
	getStart(ndims,end,divider,start);
	get_var1_double(ncid,vlid,start,&max,vtype);
	if(des<min){
		*pos=-1;
		free(start);
		return -2;
	}
	if(des>max){
		*pos=dlen-1;
		free(start);
		return -2;
	}
	while(begin<=end){
		mid=begin+(end-begin)/2;
		getStart(ndims,mid,divider,start);
		get_var1_double(ncid,vlid,start,&val,vtype);
		if(val==des){
			while(val==des&&mid<dlen-1){
				mid++;
				getStart(ndims,mid,divider,start);
				get_var1_double(ncid,vlid,start,&val,vtype);
			}
			if(mid!=dlen-1)
				mid--;
			*pos=mid;
			free(start);
			return 0;
		}
		if(begin==end){
			break;
		}
		if(val>des){
			end=mid-1;
		}else{
			begin=mid+1;
		}
	}
	
	if(mid==0){
		*pos=mid;
	}else{
		if(val<des){
			*pos=mid;
		}else{
			*pos=--mid;
		}
	}
	free(start);
	return -1;
}
/* find the position that is bigger or equal to the minimal value in '-d' argument in the descending dimension variable*/
int find_var_reverse_min(int ncid, int vlid,int ndims,size_t *divider,size_t dlen,double des,nc_type vtype,size_t *pos){
	size_t *start=(size_t *)malloc(sizeof(size_t)*ndims);	
	size_t begin=0;
	size_t mid;
	size_t end=dlen-1;
	double val,min,max;
	getStart(ndims,begin,divider,start);
	get_var1_double(ncid,vlid,start,&max,vtype);
	getStart(ndims,end,divider,start);
	get_var1_double(ncid,vlid,start,&min,vtype);
	if(des<min){
		*pos=dlen-1;
		free(start);
		return -2;
	}
	if(des>max){
		*pos=-1;
		free(start);
		return -2;
	}
	while(begin<=end){
		mid=begin+(end-begin)/2;
		getStart(ndims,mid,divider,start);
		get_var1_double(ncid,vlid,start,&val,vtype);
		if(val==des){
			while(val==des&&mid<dlen-1){
				mid++;
				getStart(ndims,mid,divider,start);
				get_var1_double(ncid,vlid,start,&val,vtype);
			}
			if(mid!=dlen-1)
				mid--;
			*pos=mid;
			free(start);
			return 0;
		}
		if(begin==end){
			break;
		}
		if(val<des){
			end=mid-1;
		}else{
			begin=mid+1;
		}
	}
	
	if(mid==0){
		*pos=mid;
	}else{
		if(val>des){
			*pos=mid;
		}else{
			*pos=--mid;
		}
	}
	free(start);
	return -1;
}
/* find the position that is smaller or equal to the maximal value in '-d' argument in the descending dimension variable*/
int find_var_reverse_max(int ncid, int vlid,int ndims,size_t *divider,size_t dlen,double des,nc_type vtype,size_t *pos){
	size_t *start=(size_t *)malloc(sizeof(size_t)*ndims);	
	size_t begin=0;
	size_t mid;
	size_t end=dlen-1;
	double min,val,max;
	getStart(ndims,begin,divider,start);
	get_var1_double(ncid,vlid,start,&max,vtype);
	getStart(ndims,end,divider,start);
	get_var1_double(ncid,vlid,start,&min,vtype);
	if(des>max){
		*pos=0;
		free(start);
		return -2;
	}
	if(des<min){
		*pos=dlen;
		free(start);
		return -2;
	}
	while(begin<=end){
		mid=begin+(end-begin)/2;
		getStart(ndims,mid,divider,start);
		get_var1_double(ncid,vlid,start,&val,vtype);
		if(val==des){
			while(val==des&&mid>0){
				mid--;
				getStart(ndims,mid,divider,start);
				get_var1_double(ncid,vlid,start,&val,vtype);
			}
			if(mid!=0)
				mid++;
			*pos=mid;
			free(start);
			return 0;
		}
		if(begin==end)
			break;
		if(val<des){
			end=mid-1;
		}else{
			begin=mid+1;
		}
	}
	if(mid==dlen-1){
		*pos=mid;
	}else{
		if(val<des){
			*pos=mid;
		}else{
			*pos=++mid;
		}
	}
	free(start);
	return -1;
}

int get_var1_double(int ncid, int vlid,size_t * start,double *data,nc_type dtype){
	switch(dtype){
		case NC_BYTE:{
				unsigned char val;
				nc_get_var1_uchar(ncid,vlid,start,&val);
				*data=val;
			 }
			break;
		case NC_CHAR:{
				signed char val;
				nc_get_var1_schar(ncid,vlid,start,&val);
				*data=val;
					 }
			break;
		case NC_SHORT:{
						  short val;
				nc_get_var1_short(ncid,vlid,start,&val);
				*data=val;	  
					  }
			break;
		case NC_INT:{
						int val;
				nc_get_var1_int(ncid,vlid,start,&val);
				*data=val;
					}
			break;
		case NC_FLOAT:{
						  float val;
			nc_get_var1_float(ncid,vlid,start,&val);
							*data=val;
					  }
			break;
		case NC_DOUBLE:{
			nc_get_var1_double(ncid,vlid,start,data);
					   }
			break;
		default:
			printf("Unknown data type in get_var1_double\n");
	}
	return 0;
}
/* detect the ascending or descending dimension varaible and set the begins and ends position*/
int set_range_opts(char **dimList,int dimNum,int ncid,int vlid,int ndims,size_t*divider,size_t dlen,nc_type vtype,size_t *begins,size_t *ends,int index){
	double min,max;
	size_t* start=(size_t*)malloc(sizeof(size_t)*ndims);
	size_t pos;
		double order1,order2;
		getStart(ndims,0,divider,start);
		get_var1_double(ncid,vlid,start,&order1,vtype);
		getStart(ndims,dlen-1,divider,start);
		get_var1_double(ncid,vlid,start,&order2,vtype);
		if(order2>=order1){
			if(dimNum>2){
				if(strlen(dimList[2])==0){
					//do nothing
				}else{
					sscanf(dimList[2],"%lf",&max);
					find_var_max(ncid,vlid,ndims,divider,dlen,max,vtype,&pos);
					if(pos>=dlen||pos<0){
						printf("max value for dimension %s is not correct!\n",dimList[0]);
						exit(1);
					}
					ends[index]=pos;
				}
			}
			if(dimNum>1){
				if(strlen(dimList[1])==0){
					//do nothing
				}else{
					sscanf(dimList[1],"%lf",&min);
					find_var_min(ncid,vlid,ndims,divider,dlen,min,vtype,&pos);
					if(pos>=dlen||pos<0){
						printf("min value for dimension %s is not correct!\n",dimList[0]);
						exit(1);
					}
					begins[index]=pos;
				}
			}
		}else{
			if(dimNum>2){
				if(strlen(dimList[2])==0){
					//don nothing
				}else{
					sscanf(dimList[2],"%lf",&max);
					find_var_reverse_max(ncid,vlid,ndims,divider,dlen,max,vtype,&pos);
					if(pos>=dlen||pos<0){
						printf("max value for dimension %s is not correct!\n",dimList[0]);
						exit(1);
					}
					begins[index]=pos;
				}
			}
			if(dimNum>1){
				if(strlen(dimList[1])==0){
					//do nothing
				}else{
					sscanf(dimList[1],"%lf",&min);
					find_var_reverse_min(ncid,vlid,ndims,divider,dlen,min,vtype,&pos);
					if(pos>=dlen||pos<0){
						printf("min value for dimension %s is not correct!\n",dimList[0]);
						exit(1);
					}
					ends[index]=pos;
				}
			}
		}
/*        printf("begins[%d] %d ends[%d] %d\n",index,begins[index],index,ends[index]);*/
	free(start);
	return 0;
}
int set_dim_opts(char **dimList,int dimNum,size_t *begins,size_t *ends,ptrdiff_t *strides,int index){
		//g->dimList[i],g->dimNum[i],begins,ends,j
		if(dimNum>3){
			if(strlen(dimList[3])==0){
					//do nothing
			}else{
				int stride;
				sscanf(dimList[3],"%d",&stride);
				strides[index]=stride;
			}
		}
		if(dimNum>2){
			if(strlen(dimList[2])==0){
					//do nothing
			}else{
				int end;
				sscanf(dimList[2],"%d",&end);
				ends[index]=end;
			}
		}
		if(dimNum>1){
			if(strlen(dimList[1])==0){
					//do nothing
			}else{
				int begin;
				sscanf(dimList[1],"%d",&begin);
				begins[index]=begin;
			}
		}
	return 0;	
}

/* del with '-d' and '-r' arguments*/
int dim_opts_handler(GLOBAL_ATTR * g,int *dimids,int ncid, int ndims,size_t *begins, size_t*ends,ptrdiff_t *strides){
	int i,j;
	int vlid,dimsNum,nattrs;
	int ids[NC_MAX_VAR_DIMS];
	int *isSet=(int *)malloc(sizeof(int)*ndims);
	memset(isSet,0,sizeof(int)*ndims);
	nc_type vtype;
	if(g->listNum){
		for(i=0;i<g->listNum;++i){
			int id,res;
			nc_inq_dimid(ncid,g->dimList[i][0],&id);
			int found=0;
            for(j=0;j<ndims;++j){
				if(id==dimids[j]){
                    found=1;
					if(isSet[j]==0){
						isSet[j]=1;
					}else{
						continue;
					}
					if(g->isRange[i]==1){//"-r "
					if((res=nc_inq_varid(ncid,g->dimList[i][0],&vlid))){
						printf("In argument '-r',variable %s is not found in the file\n",g->dimList[i][0]);
						continue;
					}
					nc_inq_var(ncid,vlid,NULL,&vtype,&dimsNum,ids,&nattrs);
					int m;
					size_t dlen=1;
					size_t *shape=(size_t *)malloc(sizeof(size_t)*dimsNum);
					size_t *divider=(size_t *)malloc(sizeof(size_t)*dimsNum);
					for(m=0;m<dimsNum;++m){
						nc_inq_dimlen(ncid,ids[m],&shape[m]);
						dlen*=shape[m];
					}
					getDivider(dimsNum,shape,divider);
					if(g->dimNum[i]>3){
						size_t s=atoi(g->dimList[i][3]);
						if(s>=dlen||s<0){
							printf("In argument '-r',stride for dimension %s is not correct in the file\n",g->dimList[i][0]);
							exit(1);
						}
						strides[j]=s;
					}
					set_range_opts(g->dimList[i],g->dimNum[i],ncid,vlid,dimsNum,divider,dlen,vtype,begins,ends,j);
					free(shape);
					free(divider);
				}else{//"-d"
					set_dim_opts(g->dimList[i],g->dimNum[i],begins,ends,strides,j);					
				}
				}
			}
            if(found==0){
                printf("ERROR: Dimension %s doesn't exist!\n", g->dimList[i][0]);
                exit(-1);
            }
		}
		
	}
	free(isSet);
	return 0;
}
int init_begins_ends_strides(int ndims,size_t **begins,int *beginNum,int *is_begins_null, size_t **ends,int *endNum,int *is_ends_null,ptrdiff_t **strides,int *strideNum,int *is_strides_null,size_t *shape){
	if(*begins==NULL){
		*begins=(size_t *)malloc(sizeof(size_t)*ndims);
		*beginNum=ndims;
		*is_begins_null=1;
	}
	if(*ends==NULL){
		*ends=(size_t *)malloc(sizeof(size_t)*ndims);
		*endNum=ndims;
		*is_ends_null=1;
	}
	if(*strides==NULL){
		*strides=(ptrdiff_t *)malloc(sizeof(ptrdiff_t)*ndims);
		*strideNum=ndims;
		*is_strides_null=1;
	}
	int i;
	for(i=0;i<ndims;++i){
		if(*is_begins_null||(*begins)[i]==-1){
			(*begins)[i]=0;
		}
		if(*is_ends_null||(*ends)[i]==-1){
			(*ends)[i]=shape[i]-1;
		}
		if(*is_strides_null||(*strides)[i]==-1){
			(*strides)[i]=1;
		}
	}
	return 0;
}
int check_begins_ends_strides(int ndims,size_t *begins,int beginNum,size_t *ends,int endNum,ptrdiff_t *strides,int strideNum,size_t *shape,char **dimsName){
	if(strideNum!=ndims){
		printf("Stride number is not correct, ndims:%d strideNum: %d\n",ndims,strideNum);
		return -1;
	}
	if(beginNum!=ndims||endNum!=ndims){
		printf("Boundary is not correct, ndims:%d beginNum:%d endNum:%d\n",ndims,beginNum,endNum);
		return -1;
	}
	int i;	
	for(i=0;i<ndims;++i){
		if(begins<0){
			printf("begin index should be  bigger than zero\n");	
			return -1;
		}
		if(ends[i]>shape[i]-1){
			printf("max index of dimension %s is %d, but user-set end is %d\n",dimsName[i],shape[i]-1,ends[i]);
			return -1;
		}
		if(begins[i]>ends[i]){
			printf("begin index shoulder not be bigger than end index\n");	
			return -1;
		}
		if(strides[i]<0||strides[0]>ends[i]-begins[i]+1){
			printf("stride of demension %s is not correct!\n",dimsName[i]);
			return -1;
		}
	}
	return 0;
}
int open_output_file_and_define_var(char *fileOut,int *ncidout,int *vlidout,int *dimidsOut,int append,int ndims,nc_type vtype,char **dimsName,size_t *begins,size_t *ends,ptrdiff_t *strides,size_t* shape,char *var,char *varOut,char **outDims,int outDimNum){
	int res;
	if (append){	
		if((res=nc_create_par(fileOut,NC_NOCLOBBER|NC_NETCDF4|NC_MPIIO,MPI_COMM_WORLD,MPI_INFO_NULL,ncidout))){
			if((res=nc_open_par(fileOut,NC_WRITE|NC_MPIIO,MPI_COMM_WORLD,MPI_INFO_NULL,ncidout))){
				printf("Try open output file %s\n failed!\n",fileOut);
				return -1;
			}
			if((res=nc_redef(*ncidout))){
				printf("Redefine netcdf id %d in output file %s failed!\n",*ncidout,fileOut);
				return -1;
			}
		}
	}else{
		if((res=nc_create_par(fileOut,NC_NETCDF4|NC_MPIIO,MPI_COMM_WORLD,MPI_INFO_NULL,ncidout))){
			printf("Cannot create output file %s\n",fileOut);
		}
	}
	int i;
	for(i=0;i<ndims;++i){
		char *name;
		if(outDims==NULL||i>outDimNum-1)
			name=dimsName[i];
		else
			name=outDims[i];
		if((res=nc_def_dim(*ncidout,name,(ends[i]-begins[i])/strides[i]+1,&dimidsOut[i]))){
			printf("Dimension %s exists in output file %s\n",name,fileOut);
			if((res=nc_inq_dimid(*ncidout,name,&dimidsOut[i]))){
				return -1;
			}
			size_t dlen=0;
			if((res=nc_inq_dimlen(*ncidout,dimidsOut[i],&dlen))){
				return -1;
			}
			if(!append&&dlen!=shape[i]){
				printf("Dimension %s conflicts, please retry --output_dims to set valid dimension names\n",name);
				return -1;
			}
		}
	}
	char *name;
	if(varOut==NULL)
		name=var;
	else
		name=varOut;
	if((res=nc_def_var(*ncidout,name,vtype,ndims,dimidsOut,vlidout))){
		printf("Define variable %s in output file %s failed, please retry --output_var to set valid variable name\n",var,fileOut);
		return -1;
	}
	return 0;
}
