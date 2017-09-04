#include <stdio.h>
#include <stdlib.h>
#include "CommFun.h"


// write data to target file
void WriteIntFile(char *filename, int n, int *ptr)
{
	char filename1[128];
	sprintf(filename1,"%s.txt",filename);
	FILE *fp = fopen(filename1,"wt");
	int idx;
	for (idx=0; idx<n; idx++)
	{
		fprintf(fp,"%d\n",ptr[idx]);
	}
	fclose(fp);
}

void WriteFltFile(char *filename, int n, float *ptr)
{
	char filename1[128];
	sprintf(filename1,"%s.txt",filename);
	FILE *fp = fopen(filename1,"wt");
	int idx;
	for (idx=0; idx<n; idx++)
	{
		fprintf(fp,"%f\n",ptr[idx]);
	}
	fclose(fp);
}

void WriteCompxFile(char *filename, int n, compx *ptr)
{
	char filename1[128];
	char filename2[128];
	sprintf(filename1,"%s_real.txt",filename);
    sprintf(filename2,"%s_imag.txt",filename);

	FILE *fp1 = fopen(filename1,"wt");
	FILE *fp2 = fopen(filename2,"wt");
	int idx;
	for (idx=0; idx<n; idx++)
	{
		fprintf(fp1,"%f\n",ptr[idx].real);
		fprintf(fp2,"%f\n",ptr[idx].imag);
	}
	fclose(fp1);
	fclose(fp2);
}

// read dec data from target file
void ReadIntFile(char *filename, int n, int *ptr)
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

