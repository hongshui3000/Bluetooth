#include <stdio.h>
#include <stdlib.h>
//#include <iostream.h>

// write dec data to target file
void WriteDecFile(char *filename, int n, int *ptr)
{
	FILE *fp = fopen(filename,"wt");
	int idx;
	for (idx=0; idx<n; idx++)
	{
		fprintf(fp,"%d\n",ptr[idx]);
	}
	fclose(fp);
}

void WriteFltFile(char *filename, int n, float *ptr)
{
	FILE *fp = fopen(filename,"wt");
	int idx;
	for (idx=0; idx<n; idx++)
	{
		fprintf(fp,"%f\n",ptr[idx]);
	}
	fclose(fp);
}

// read dec data from target file
void ReadDecFile(char *filename, int n, int *ptr)
{
	FILE *fp = fopen(filename,"r");
	int idx;
	char temp[16];
	for(idx=0; idx<n; idx++)
	{
         fscanf(fp,"%s",&temp);
		 ptr[idx] = atoi(temp);
	}
}

void ReadFltFile(char *filename, int n, float *ptr)
{
	FILE *fp = fopen(filename,"r");
	int idx;
	char temp[16];
	for(idx=0; idx<n; idx++)
	{
         fscanf(fp,"%s",&temp);
		 ptr[idx] = atof(temp);
	}
}
