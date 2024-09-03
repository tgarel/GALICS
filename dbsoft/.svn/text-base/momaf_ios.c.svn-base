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

int get_ng_(char *file)
{
  int  ng;
  FILE *fd;
  if(!(fd = fopen(file,"r")))
    {
      printf("can't open file %s \n", file);
      exit(1);
    }
  fread(&ng,sizeof(int),1,fd);
  fclose(fd);
  return(ng);
}

void read_momaf_propfile_(char *file, long long *GalID, long long *HaloID, float *Inclination, float *Pos, float *Vel)
{
  int       i,ng;
  FILE      *fd;

  if(!(fd = fopen(file,"r")))
    {
      printf("can't open file %s \n", file);
      exit(1);
    }
  fread(&ng,sizeof(int),1,fd);
  for (i = 0; i < ng; i++) 
  {
      fread(&GalID[i],      sizeof(long long),1,fd);
      fread(&HaloID[i],     sizeof(long long),1,fd);
      fread(&Inclination[i],sizeof(float),1,fd);
      fread(&Pos[3*i],      sizeof(float),1,fd);
      fread(&Pos[3*i+1],    sizeof(float),1,fd);
      fread(&Pos[3*i+2],    sizeof(float),1,fd);
      fread(&Vel[3*i],      sizeof(float),1,fd);
      fread(&Vel[3*i+1],    sizeof(float),1,fd);
      fread(&Vel[3*i+2],    sizeof(float),1,fd);
  }
  fclose(fd);
}

void read_momaf_magfile_(char *file, float *mag, float *dmag)
{
  int       i,ng;
  FILE      *fd;

  if(!(fd = fopen(file,"r")))
    {
      printf("can't open file `%s'\n", file);
      exit(1);
    }
  fread(&ng,sizeof(int),1,fd);
  for (i = 0; i < ng; i++) 
    {
      fread(&mag[i],sizeof(float),1,fd);
      fread(&dmag[i],sizeof(float),1,fd);
    }
  fclose(fd);
}


void write_halo_momaf_file_(char *file, int *nh, long long *HaloID, float *Pos, float *Vel)
{
  int       i;
  FILE      *fd;

  if(!(fd = fopen(file,"w")))
    {
      printf("can't open file %s \n", file);
      exit(1);
    }
  fwrite(nh,sizeof(int),1,fd);
  for (i = 0; i < *nh; i++) 
  {
      fwrite(&HaloID[i],     sizeof(long long),1,fd);
      fwrite(&Pos[3*i],      sizeof(float),1,fd);
      fwrite(&Pos[3*i+1],    sizeof(float),1,fd);
      fwrite(&Pos[3*i+2],    sizeof(float),1,fd);
      fwrite(&Vel[3*i],      sizeof(float),1,fd);
      fwrite(&Vel[3*i+1],    sizeof(float),1,fd);
      fwrite(&Vel[3*i+2],    sizeof(float),1,fd);
  }
  fclose(fd);
}

void write_empty_halo_momaf_file_(char *file)
{
  int       nh;
  FILE      *fd;

  if(!(fd = fopen(file,"w")))
    {
      printf("can't open file %s \n", file);
      exit(1);
    }
  nh = 0;
  fwrite(&nh,sizeof(int),1,fd);
  fclose(fd);
}
