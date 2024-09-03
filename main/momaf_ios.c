#include <stdio.h>
#include <stdlib.h>

void write_momaf_magfile_(char *file, float *mag, float *dmag, int *ng)
{
  int       i;
  FILE      *fd;

  if(!(fd = fopen(file,"w")))
    {
      printf("can't open file `%s'\n", file);
      exit(1);
    }
  fwrite(ng,sizeof(int),1,fd);
  for (i = 0; i < *ng; i++) 
    {
      fwrite(&mag[i],sizeof(float),1,fd);
      fwrite(&dmag[i],sizeof(float),1,fd);
    }
  fclose(fd);
}
 

void write_momaf_empty_magfile_(char *file)
{
  int       i;
  FILE      *fd;

  if(!(fd = fopen(file,"w")))
    {
      printf("can't open file `%s'\n", file);
      exit(1);
    }
  i = 0;
  fwrite(&i,sizeof(int),1,fd);
  fclose(fd);
}

void write_momaf_propfile_(char *file, int *ng, long long *GalID, long long *HaloID, float *Inclination, float *Pos, float *Vel)
{
  int       i;
  FILE      *fd;

  if(!(fd = fopen(file,"w")))
    {
      printf("can't open file %s \n", file);
      exit(1);
    }
  fwrite(ng,sizeof(int),1,fd);
  for (i = 0; i < *ng; i++) 
  {
      fwrite(&GalID[i],      sizeof(long long),1,fd);
      fwrite(&HaloID[i],     sizeof(long long),1,fd);
      fwrite(&Inclination[i],sizeof(float),1,fd);
      fwrite(&Pos[3*i],      sizeof(float),1,fd);
      fwrite(&Pos[3*i+1],    sizeof(float),1,fd);
      fwrite(&Pos[3*i+2],    sizeof(float),1,fd);
      fwrite(&Vel[3*i],      sizeof(float),1,fd);
      fwrite(&Vel[3*i+1],    sizeof(float),1,fd);
      fwrite(&Vel[3*i+2],    sizeof(float),1,fd);
  }
  fclose(fd);
}

void write_momaf_empty_propfile_(char *file)
{
  int       i;
  FILE      *fd;

  if(!(fd = fopen(file,"w")))
    {
      printf("can't open file %s \n", file);
      exit(1);
    }
  i = 0;
  fwrite(&i,sizeof(int),1,fd);
  fclose(fd);
}

void update_momaf_with_ids_(char *file, long long *GalID, long long *HaloID)
{
  int  i,ng,skip;
  FILE *fd;
  if(!(fd = fopen(file,"r+")))
    {
      printf("can't open file %s \n", file);
      exit(1);
    }
  fread(&ng,sizeof(int),1,fd);
  skip = 7 * sizeof(float);
  for (i=0; i<ng; i++) 
    {
      fwrite(&GalID[i], sizeof(long long),1,fd);
      fwrite(&HaloID[i],sizeof(long long),1,fd);
      fseek(fd,skip,SEEK_CUR); /* skip incl, pos, and vel*/
    }
  fclose(fd);
}

