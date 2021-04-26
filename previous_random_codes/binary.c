#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

struct io_header_1
{
  unsigned int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  int npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8];	/* fills to 256 Bytes */
} header1;



int NumPart, Ngas;

struct particle_data
{
  float Pos[3];
  float Vel[3];
  float Mass;
  int Type;

  float Rho, U, Temp, Ne;
} *P;

int *Id;

double Time, Redshift;

int main()
{
  int i,j,n;
  char posfile[80],velfile[80], outfile[80];
  FILE *out,*pos,*vel;

  for(i=0;i<6;i++) header1.npart[i]=0;
  for(i=0;i<6;i++) header1.mass[i]=0.0;

  printf("Enter number of halo and disk particles (comma separated) >");
  scanf("%d,%d",&header1.npart[1],&header1.npart[2]);

  printf("Enter mass of halo and disk particles (comma separated) >");
  scanf("%lf,%lf",&header1.mass[1],&header1.mass[2]);

  header1.time = 0.0;
  header1.redshift = 0.0;
  header1.flag_sfr = 0;
  header1.flag_feedback = 0;

  for(i=0;i<6;i++){
    header1.npartTotal[i] = header1.npart[i];
  }

  header1.flag_cooling = 0;
  header1.num_files = 0;
  header1.BoxSize = 0.0;
  header1.Omega0 = 0.0;
  header1.OmegaLambda = 0.0;
  header1.HubbleParam = 0.0;

  n = 0;
  for(i=0;i<6;i++){
    n = n + header1.npartTotal[i];
  }

  float Pos[n][3], Vel[n][3];
  int Id[n];
  float Mass[n];
  int Type[n];
  float Rho[n],U[n],Temp[n],Ne[n];

  printf("Enter name of position file (format x y z) >");
  scanf("%s",posfile);

  printf("Enter name of velocity file (format vx vy vz) >");
  scanf("%s",velfile);

  if((pos = fopen(posfile,"r")) == NULL) {
    printf("\nCannot find position file\n");
  }

  i=0;
  while(fscanf(pos,"%f %f %f",&Pos[i][0],&Pos[i][1],&Pos[i][2]) != EOF) i++;

  fclose(pos);

  if((vel = fopen(velfile,"r")) == NULL) {
    printf("\nCannot find position file\n");
  }

  i=0;
  while(fscanf(vel,"%f %f %f",&Vel[i][0],&Vel[i][1],&Vel[i][2]) != EOF) i++;

  fclose(vel);

  for(i=0;i<n;i++){
    Id[i] = i+1;
    Rho[i] = 0.0;
    U[i] = 0.0;
    Temp[i] = 0.0;
    Ne[i] = 0.0;
  }

  for(i=0;i<6;i++){
    for(j=0;j<header1.npartTotal[i];j++){
      Mass[j] = header1.mass[i];
      Type[j] = j;
    }
  }

  int dummy1, dummy2;
  dummy1=256;
  dummy2=dummy1;

  printf("Enter name of initial condition file to make >");
  scanf("%s",outfile);

  if((out = fopen(outfile,"wb")) == NULL) {
    printf("\nCannot find output file\n");
  }

  fwrite(&dummy1, sizeof(dummy1), 1, out);
  fwrite(&header1, sizeof(struct io_header_1), 1, out);
  fwrite(&dummy2, sizeof(dummy2), 1, out);
  dummy1=n*3*4;
  dummy2=dummy1;
  fwrite(&dummy1, sizeof(dummy1), 1, out);
  fwrite(&Pos, sizeof(Pos), 1, out);
  fwrite(&dummy2, sizeof(dummy2), 1, out);
  fwrite(&dummy1, sizeof(dummy1), 1, out);
  fwrite(&Vel, sizeof(Vel), 1, out);
  fwrite(&dummy2, sizeof(dummy2), 1, out);
  dummy1=n*4;
  dummy2=dummy1;
  fwrite(&dummy1, sizeof(dummy1), 1, out);
  fwrite(&Id, sizeof(Id), 1, out);
  fwrite(&dummy2, sizeof(dummy2), 1, out);

  fclose(out);

  return 0;
}
