#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include "pnco_parser.h"
GLOBAL_ATTR global_attr;
void init_global_attr(GLOBAL_ATTR *g){
	g->fileInNum=0;
	g->fileIns=NULL;
	g->varNum=0;
	g->varList=NULL;
	g->begin=NULL;
	g->beginNum=0;
	g->end=NULL;
	g->endNum=0;
	g->strides=NULL;
	g->strideNum=0;
	g->fileOut=NULL;
	g->op=ADD;
	g->append=0;
	g->concat=0;
	g->newSample=0;
	g->outputVar=NULL;
	g->outputDims=NULL;
	g->outputDimNum=0;
	g->argv=NULL;
	g->argc=0;
	g->listNum=0;
};
void display_usage(){
	printf("pnco [-h] [-c] [-b ,,,] [-d ,,,] [-e ,,,] [-i ,,,] [-v] [-s ,,,] [-m] [-o]\n");

	printf("  [-a] or [--append]\t set output file which is appended, if not set, out.nc is the default output file\n");
	printf("  [-h] or [--help]\t display usage\n");
	printf("  [-c] or [--concat]\t concat a variable in many files to a file, the variable should exits in each file and has the same dimension shape except the first dimension. When concating, the begin_dims and ends_dims are not valid, the first stride must be 1\n");
	printf("  [-b] or [--begin_dims]\t set variable begin positions which are seperated by ',' the length should be the same as the length of variable dimension. if not set, the default begin_dims are 0,0,0 when the varialbe is v[4][6][8]\n");
	printf("  [-d] or [--dim]\t set dimension's begin positon ,end position and stride with name, start and stop position which are seperated by ','.For varialbe v[d1][d2],use '-d d1,0,25,2': d1 is the name,0 is the start position, 25 is the stop position, 2 is the stride vaule. Attention:'-d' will override the result of -b -e and -s \n");
	printf("  [-r] or [--range]\t set dimension's begin positon ,end position and stride with name, minimal and maximal dimension variable's value which are seperated by ','.For varialbe v[d1][d2],use '-r d1,1,20,2': d1 is the name,1 is the minimal value, 20 is the maximal value, 2 is the stride vaule. This works when d1(d1) variable exists in the nc file and has ascending or descending values. Attention:'-d' will override the result of -b -e and -s \n");
	printf("  [-e] or [--end_dims]\t set variable ends positions which are seperated by ',' the length should be the same as the length of variable dimension. if not set, the default end_dims are 3,5,7 when the varialbe is v[4][6][8]\n");
	printf("  [-i] or [--input_files]\t set input files which are seperated by ','\n");
	printf("  [-v] or [--input_var]\t set input variable name, record variable(the first dimension is unlimited) is supported\n");
	printf("  [-s] or [--stride]\t set variable strides which are seperated by ',',the length should be the same as the length of variable dimension. if not set, the default strides are 1,1,1 when the varialbe is v[4][6][8], if strides are 2,2,2, the output variable is v[2][3][4], if strides are 3,3,3 the output variable is v[2][2][3]\n");
	printf("  [-m] or [--method]\t set the binary operator method. Binary operation is used only when two input files are set and -m arguments exists.Currently support ADD SUB MUL DIV MOD AVG five methods.\n");
	printf("  [-o] or [--override]\t set output file which is overwritten, if not set, out.nc is the default output file\n");
	printf("  [--output_var]\t set output variable name in output file, if not set, use varialbe name in the first input file\n");
	printf("  [--output_dims]\t set output dimension names in output file, if not set, use dimension names in the first input file. Besides unlimited dimension is not supported in the output file\n");
	printf(" exapmle:\n");
	printf("\t pnco -i file.nc -v v1 -dim d1,0,12,2 -o result.nc   or   pnco --input_files file.nc --input_var v1 --dim d1,0,12,2 --override result.nc\n");
	printf("\t   extract two dimension variable v1 with the strides 2,2 from file.nc to result.nc\n");
	printf("\t pnco -i file1.nc,file2.nc -v v1 -s 2,2 -m avg -o result.nc --output_var v2 --output_dims dim1,dim2\n");
	printf("\t   get the average result of variable v1 with strides 2,2 from file1.nc and file2.nc to result.nc, the output variable and dimension names are set\n");
	printf("\t pnco -i file1.nc,file2.nc -v v1 -s 2,2 -m avg -o result.nc\n");
	printf("\t   get the average result of variable v1 with strides 2,2 from file1.nc and file2.nc to result.nc\n");
	printf("\t pnco -c -i file1.nc,file2.nc,file3.nc,file4.nc -v v1 -o result.nc\n");
	printf("\t   concat variable v1 from file1.nc file2.nc file3.nc and file4.nc to result.nc\n");
			
}
/*int parse_names(char ***list,int *len,char *value){*/
/*    char *p;*/
/*    *len=1;*/
/*    p=value; */
/*    while(p){// calculate len*/
/*        p=strpbrk(p,DELIM);*/
/*        if(p){*/
/*            p=p+1;*/
/*            (*len)++;*/
/*            //printf("aaa:%s\n",p);*/
/*        }*/
/*    }*/
/*    //printf("len:%d\n",(*len));*/
/*    int i=1;*/
/*    *list=(char **)malloc(sizeof(char *)*(*len));*/
/*    p=strtok(value,DELIM);*/
/*    (*list)[0]=(char *)malloc(strlen(p)+1);*/
/*    strcpy((*list)[0],p);*/
/*    //printf("strcpy %s\n",(*list)[0]);*/
/*    while((p=strtok(NULL,DELIM))){*/
/*        //	printf("%s ",p);*/
/*        (*list)[i]=(char *)malloc(strlen(p)+1);*/
/*        strcpy((*list)[i],p);*/
/*        //printf("!!:%s\n",(*list)[i]);*/
/*        ++i;*/
/*    }*/
/*    //print_info("",*list,*len);*/
/*    return 0;*/
/*}*/
int parse_names(char ***list,int *len,char *value){
	char *p;
	char *start;
	*len=1;
	p=value;
	while(*p!='\0'){
		if(*p==DELIMA){
			++(*len);
		}
		++p;
	}
	*list=(char **)malloc(sizeof(char *)*(*len));
	
	int i=0;
	p=value;
	start=value;
	while(*p!='\0'){
		if(*p==DELIMA){
			if(*start!=DELIMA){
				*p='\0';
				(*list)[i]=(char *)malloc(strlen(start)+1);
				strcpy((*list)[i],start);
			}else{
				*p='\0';
				*start='\0';
				(*list)[i]=(char *)malloc(1);
				(*list)[i][0]='\0';
			}
			++i;
			++p;
			if(*p!='\0')
				start=p;
			continue;
		}
		++p;
	}
	(*list)[i]=(char *)malloc(strlen(start)+1);
	strcpy((*list)[i],start);
/*    for(i=0;i<(*len);++i){*/
/*        printf("%d %s\n",i,(*list)[i]);*/
/*    }*/
	return 0;
}
int parse_arr(size_t **arr,int *len,char *value){
	char *p;
	char *start;
	*len=1;
	p=value;
	while(*p!='\0'){
		if(*p==DELIMA){
			++(*len);
		}
		++p;
	}
	*arr=(size_t *)malloc(sizeof(size_t *)*(*len));
	
	int i=0;
	p=value;
	start=value;
	while(*p!='\0'){
		if(*p==DELIMA){
			//printf("%d\n",*start);
			if(*start!=DELIMA){
				*p='\0';
				(*arr)[i]=atoi(start);
			}else{
				*p='\0';
				*start='\0';
				(*arr)[i]=-1;
			}
			++i;
			++p;
			if(*p!='\0')
				start=p;
			continue;
		}
		++p;
	}
	(*arr)[i]=atoi(start);
/*    for(i=0;i<(*len);++i){*/
/*        printf("%d %d\n",i,(*arr)[i]);*/
/*    }*/
	return 0;
}

/*int parse_arr(size_t **arr,int *len,char *value){*/
/*    char *p;*/
/*    *len=1;*/
/*    p=value; */
/*    while(p){// calculate len*/
/*        p=strpbrk(p,DELIM);*/
/*        if(p){*/
/*            p=p+1;*/
/*            (*len)++;*/
/*            //printf("aaa:%s\n",p);*/
/*        }*/
/*    }*/
/*    *arr=(int *)malloc(sizeof(int)*(*len));*/
/*    //printf("len:%d\n",(*len));*/
/*    int i=1;*/
/*    p=strtok(value,DELIM);*/
/*    (*arr)[0]=atoi(p);*/
/*    if((*arr)[0]==0){*/
/*        printf("arr[0] is zero!\n");*/
/*        exit(1);*/
/*    }*/
/*    while((p=strtok(NULL,DELIM))){*/
/*        (*arr)[i]=atoi(p);*/
/*        if((*arr)[i]==0){*/
/*            printf("arr[%d] is zero\n",i);*/
/*            exit(1);*/
/*        }*/
/*        ++i;*/
/*    }*/
/*|+    for(i=0;i<*len;++i)+|*/
/*|+        printf("%d ",(*strides)[i]);+|*/
/*|+    printf("\n");+|*/
/*    return 0;*/
/*}*/
int parse_op(BIN_OP *op,char *value){
	if(strcasecmp(value,"add")==0){
		*op=ADD;
	}else if(strcasecmp(value,"sub")==0){
		*op=SUB;
	}else if(strcasecmp(value,"mul")==0){
		*op=MUL;
	}else if(strcasecmp(value,"div")==0){
		*op=DIV;
	}else if(strcasecmp(value,"mod")==0){
		*op=MOD;
	}else if(strcasecmp(value,"avg")==0){
		*op=AVG;
	}
	return 0;
}
const struct option longopt[]={
	{"append",1,NULL,'a'},
	{"input_files",1,NULL,'i'},
	{"input_var",1,NULL,'v'},
	{"override",1,NULL,'o'},
	{"dim",1,NULL,'d'},
	{"range",1,NULL,'r'},
	{"begin_dims",1,NULL,'b'},
	{"end_dims",1,NULL,'e'},
	{"stride",1,NULL,'s'},
	{"concat",0,NULL,'c'},
	{"help",0,NULL,'h'},
	{"method",1,NULL,'m'},
	{"output_var",1,NULL,1},
	{"output_dims",1,NULL,2},
	{NULL,0,NULL,0}
};
int pnco_parser(int argc,char **argv){
	init_global_attr(&global_attr);
	global_attr.argc=argc;
	global_attr.argv=(char **)malloc(sizeof(char *)*argc);
	int i;
	for(i=0;i<argc;++i){
		global_attr.argv[i]=(char *)malloc(strlen(argv[i])+1);
		strcpy(global_attr.argv[i],argv[i]);
	}
	int ch;
	//opterr=0;
	while((ch=getopt_long(argc,argv,"a:chni:v:d:s:m:o:b:e:r:1:2",longopt,NULL))!=-1){
/*    while((ch=getopt(argc,argv,"achni:v:d:s:m:o:b:e:"))!=-1){*/
		switch(ch){
			case 'a':
				global_attr.append=1;
				global_attr.fileOut=(char *)malloc(strlen(optarg)+1);
				strcpy(global_attr.fileOut,optarg);
				break;
			case 'b':
				parse_arr(&global_attr.begin,&global_attr.beginNum,optarg);
				break;
			case 'e':
				parse_arr(&global_attr.end,&global_attr.endNum,optarg);
				break;
			case 'c':
				global_attr.concat=1;
				break;
			case 'd':
				parse_names(&global_attr.dimList[global_attr.listNum],&global_attr.dimNum[global_attr.listNum],optarg);
				global_attr.isRange[global_attr.listNum]=0;
				global_attr.listNum++;
				break;
			case 'h':
				display_usage();
				break;
			case 'n':
				global_attr.newSample=1;
				break;
			case 'i':
				parse_names(&global_attr.fileIns,&global_attr.fileInNum,optarg);
				break;
			case 'v':
				parse_names(&global_attr.varList,&global_attr.varNum,optarg);
				break;
			case 'm':
				parse_op(&global_attr.op,optarg);
				break;
			case 'r':
				parse_names(&global_attr.dimList[global_attr.listNum],&global_attr.dimNum[global_attr.listNum],optarg);
				global_attr.isRange[global_attr.listNum]=1;
				global_attr.listNum++;
				break;
			case 's':
				parse_arr(&global_attr.strides,&global_attr.strideNum,optarg);
				break;
			case 'o':
				global_attr.append=0;
				global_attr.fileOut=(char *)malloc(strlen(optarg)+1);
				strcpy(global_attr.fileOut,optarg);
				break;
			case 1:
				global_attr.outputVar=(char *)malloc(strlen(optarg)+1);
				strcpy(global_attr.outputVar,optarg);
				break;
			case 2:
				parse_names(&global_attr.outputDims,&global_attr.outputDimNum,optarg);
				break;
			default:	
				printf("Unknown option :%c\n",ch);
		}
	}
	return 0;
}
