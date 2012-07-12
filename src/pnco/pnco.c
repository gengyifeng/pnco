#include <string.h>
#include "pnco.h"
int init_dims(DIMS *dims,int size, size_t * arr,char **names){
	int i;
	dims->size=size;
	dims->shape=(size_t *)malloc(sizeof(size_t)*size);
	dims->names=(void *)malloc(sizeof(void *)*size);
	for(i=0;i<size;++i){
		dims->shape[i]=arr[i];
		dims->names[i]=(char*)malloc(strlen(names[i])+1);
		strcpy(dims->names[i],names[i]);
	}	
	return 0;
}
int init_var(VAR *var,nc_type type,char *varName,DIMS *dims){
	var->type=type;
	var->varName=(char *)malloc(strlen(varName)+1);
	strcpy(var->varName,varName);
	var->dims=*(DIMS*)malloc(sizeof(DIMS));
	init_dims(&var->dims,dims->size,dims->shape,dims->names);
	return 0;
}
int init_box(BOX *box,int size,size_t *start,size_t *count){
	int i;
	box->size=size;
	box->start=(size_t *)malloc(sizeof(size_t)*size);
	box->count=(size_t *)malloc(sizeof(size_t)*size);
	box->stride=(size_t *)malloc(sizeof(size_t)*size);
	for(i=0;i<size;++i){
		box->start[i]=start[i];
		box->count[i]=count[i];
		box->stride[i]=1;
	}
	return 0;
}
int init_box_with_stride(BOX *box,int size,size_t *start,size_t *count,size_t *stride){
	int i;
	box->size=size;
	box->start=(size_t *)malloc(sizeof(size_t)*size);
	box->count=(size_t *)malloc(sizeof(size_t)*size);
	box->stride=(size_t *)malloc(sizeof(size_t)*size);
	for(i=0;i<size;++i){
		box->start[i]=start[i];
		box->count[i]=count[i];
		box->stride[i]=stride[i];
	}
	return 0;
}
