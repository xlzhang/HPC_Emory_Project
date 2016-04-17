
//Domer_v2.0
//SRMC regular potential
//Ziwei Guo  Feb-22-2016


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define  PI 3.1415926
#define  NA 6.022e23
#define  K (8.9875517873681764e9)*(1.6021892e-19)*(1.6021892e-19)*pow(10,9)*NA/1000 /* unit: N·m2/C2 */
#define  EPSILON 5.0 /*unit:kJ/mol*/
#define  SIGMA 4.7 /*unit:A*/
#define  L_A 4*EPSILON*pow(SIGMA,12)
#define  L_B 4*EPSILON*pow(SIGMA,6)
#define  RC_LJ 12.0 /*unit:A*/
#define  RC_SW 9.0 /*unit:A*/
#define  frand() ((double) rand() / (RAND_MAX+1.0))

double L6[3],L12[3];

double min(double a, double b)
	{return ((a)>(b)) ? (b):(a);}

void copy_array(double origin[][3], double new[][3], int n)
{
	int x,y;
	for(x=0;x<n;x++)
    {
    	for(y=0;y<3;y++)
    	{
    		new[x][y] = origin[x][y];
    	}
    }
}

void copy_array_2(double origin[3], double new[3])
{
	int x;
	for(x=0;x<3;x++)
    	new[x] = origin[x];
}

void read_input(char address_input[], char InputConfigName[], double sys_para[], double sim_para[])
{
	FILE *fp;
	char line[100];

	fp = fopen(address_input, "r");
	fgets(line, 100, fp);
	fgets(line, 100, fp);
	sscanf(line,"%s", InputConfigName);

	fgets(line, 100, fp);
	fgets(line, 100, fp);
	sscanf(line,"%lf%lf%lf",&sys_para[0],&sys_para[1],&sys_para[2]);

	fgets(line, 100, fp);
	fgets(line, 100, fp);
	fgets(line, 100, fp);
	sscanf(line,"%lf%lf%lf%lf",&sim_para[0],&sim_para[1],&sim_para[2],&sim_para[3]);
	fgets(line, 100, fp);
	fgets(line, 100, fp);
	sscanf(line,"%lf",&sim_para[4]);
}

void read_coordinate(double coordinate[][3], double box[], char address[], char Title[], int natom)
{
	int i;
	char line[100];
	FILE *fp;
	fp = fopen(address, "r");
	fgets(line, 100, fp);
	sscanf(line,"%s", Title);
	fgets(line, 100, fp);
	sscanf(line,"%5d", &natom);
	for (i=0; i<natom; i++)
	{
			fgets(line, 100, fp);
			sscanf(&line[22],"%lf%lf%lf",&coordinate[i][0],&coordinate[i][1],&coordinate[i][2]);
	}
	fscanf(fp,"%lf%lf%lf\n",&box[0],&box[1],&box[2]);
	fclose(fp);
}

double distance(double i[], double j[], double L)
{
	double rx = i[0]-j[0], ry = i[1]-j[1], rz = i[2]-j[2];
	rx = rx - L*round(rx/L);
	ry = ry - L*round(ry/L);
	rz = rz - L*round(rz/L);
	return sqrt(rx*rx + ry*ry + rz*rz);
}

double distance_sqr(double i[], double j[], double L)
{
	double rx = i[0]-j[0], ry = i[1]-j[1], rz = i[2]-j[2];
	rx = rx - L*round(rx/L);
	ry = ry - L*round(ry/L);
	rz = rz - L*round(rz/L);
	return rx*rx + ry*ry + rz*rz;
}

void distance_pbc(double i[], double j[], double L, double r[])  /*rji  note the order of i,j*/
{
	double rx = i[0]-j[0], ry = i[1]-j[1], rz = i[2]-j[2];
	if(rx > L/2)
		{rx = rx - L;}
	else if(rx < -L/2)
		{rx = rx + L;}
	if(ry > L/2)
		{ry = ry - L;}
	else if(ry < -L/2)
		{ry = ry + L;}
	if(rz > L/2)
		{rz = rz - L;}
	else if(rz < -L/2)
		{rz = rz + L;}
	r[0] = rx;
	r[1] = ry;
	r[2] = rz;
}

void nojump(double L, double r[])  /*ensure this site does not jump out of box !!!Set the origin at the center of the box!!!*/
{
	if(r[0]<-L/2)
		{r[0] += L;}
	else if(r[0]>L/2)
		{r[0] -= L;}
	if(r[1]<-L/2)
		{r[1] += L;}
	else if(r[1]>L/2)
		{r[1] -= L;}
	if(r[2]<-L/2)
		{r[2] += L;}
	else if(r[2]>L/2)
		{r[2] -= L;}
}

void shift_potential(int a, double LJ_COEFF[])
{
	LJ_COEFF[0] = -a*((a+4)*RC_LJ-(a+1)*RC_SW)/(pow(RC_LJ,a+2)*pow((RC_LJ-RC_SW),2));
	LJ_COEFF[1] = a*((a+3)*RC_LJ-(a+1)*RC_SW)/(pow(RC_LJ,a+2)*pow((RC_LJ-RC_SW),3));
	LJ_COEFF[2] = 1/pow(RC_LJ,a)-LJ_COEFF[0]*pow((RC_LJ-RC_SW),3)/3-LJ_COEFF[1]*pow((RC_LJ-RC_SW),4)/4;
}

double LJ(double r) /* unit: kJ/mol */
{
	if(r <= RC_SW/10)
		return (L_A*pow(1/(10*r),12)-L_B*pow(1/(10*r),6))+(L_B*L6[2]-L_A*L12[2]);
	else
		return -L_B*(1/pow((10*r),6)- L6[0]*pow((r*10-RC_SW),3)/3-L6[1]*pow((r*10-RC_SW),4)/4-L6[2]) + L_A*(1/pow((10*r),12)- L12[0]*pow((r*10-RC_SW),3)/3-L12[1]*pow((r*10-RC_SW),4)/4-L12[2]);
}

void translate(double L, double OS[], double coordinate[3], double New[3])
{
	int i;
	double x,y,z,X,Y,Z;
	x = OS[0];
	y = OS[1];
	z = OS[2];
	for(i=0;i<3;i++)
	{
		X = coordinate[0]-x;
		Y = coordinate[1]-y;
		Z = coordinate[2]-z;
		New[0] = X;
		New[1] = Y;
		New[2] = Z;
		nojump(L, New);
	}
}

void rand_sphere(double center[], double coordinate[], double r)  /*generate a vector in a sphere*/
{
	double rsq = 1000.0;
	while(rsq >= r*r)
	{
		coordinate[0] = 2*r*(frand()-0.5);
		coordinate[1] = 2*r*(frand()-0.5);
		coordinate[2] = 2*r*(frand()-0.5);
		rsq = coordinate[0]*coordinate[0]+coordinate[1]*coordinate[1]+coordinate[2]*coordinate[2];
	}
	coordinate[0] += center[0];
	coordinate[1] += center[1];
	coordinate[2] += center[2];
}

double Usite(double coordinate[][3], int imol, int N, double L)
{
	int tmol;
	double ULJ=0.0,rc;
	for(tmol=0;tmol<N;tmol++)
	{
		if(tmol==imol)
			continue;
		else
		{
			rc = distance(coordinate[imol], coordinate[tmol], L);
//			printf("imol:%f %f %f\n",coordinate[imol][0],coordinate[imol][1],coordinate[imol][2]);
//			printf("jmol:%f %f %f\n",coordinate[tmol][0],coordinate[tmol][1],coordinate[tmol][2]);
//			printf("???tmol%d  imol%d  rrrrr:%f L%f\n",tmol,imol,rc,L);
			if(rc < RC_LJ/10)
				ULJ += LJ(rc);
//			printf("LJ(rc):%f\n",LJ(rc));
		}
	}
//	printf("ULJ:%.20f   UC:%.20f   Usite:%.20f\n", ULJ,Ucoulomb,ULJ+Ucoulomb);
	return ULJ;
}

double Utotal(double coordinate[][3], int N, double L)
{
	int imol;
	double tot = 0.0;
	for(imol=0;imol<N;imol++)
	{
		tot += Usite(coordinate, imol, N, L);
//		printf("imol:%d  check:%f\n",imol,Usite(coordinate, imol, N, L));
	}

	return tot/2;
}

int serach_number(int n,int array[],int len)       /*used for check a int number is in a array or not,and return the index of the number in the array*/
{
	int i,judge = -1;
	for(i=0;i<len;i++)
	{
		if(n == array[i])
		{
			judge = i;
			break;
		}
	}
	return judge;
}

int Select(double w[], double sumw)
{
	int n;
	double ws,cumw;
	ws = frand()*sumw;
	cumw = w[0];
	n = 0;
	while(cumw < ws)
	{
		n++;
		cumw += w[n];
	}
	return n;
}

double Usite_cav(double coordinate[][3], int imol, int N, double L, int index_sel[], int nsel, int isel)           /*used for calculating a particle in cavity*/
{           															/*index_sel takes the index of molecule to be deleted, isel is the order for insertion*/
	int tmol, s;
	double ULJ=0.0,rc;
	for(tmol=0;tmol<N;tmol++)
	{
		s = serach_number(tmol,index_sel,nsel);
		if(tmol==imol)
		{
			continue;
		}
		else if((s!=-1)&&(s>isel))
		{
			continue;
		}
		else
		{
			rc = distance(coordinate[imol], coordinate[tmol], L);
			if(rc < (RC_LJ/10))
				ULJ += LJ(rc);
		}
	}
	return ULJ;
}

double gen_config(char F, double coordinate[][3],double cps[][3],double temp_cps[][3], double r_cav, int num_cps,int index_sel[],int index_sel_cps[],int nsel,int k,double beta,double center[],double L,int imol, FILE* fplog)
{
	int j,jbook=0;

	double alpha=1, sumw=0, select[5000][3], w[5000] = {0}, w_t[5000] = {0};


	for(j=0;j<k;j++)
	{
		if(j==0 && F == 'O')
			copy_array_2(temp_cps[index_sel_cps[imol]], cps[index_sel_cps[imol]]);
		else
		{
			rand_sphere(center, cps[index_sel_cps[imol]], r_cav);
			nojump(L,cps[index_sel_cps[imol]]);
		}
		copy_array_2(cps[index_sel_cps[imol]],select[j]);
		w[j] =  Usite_cav(cps, index_sel_cps[imol], num_cps, L, index_sel_cps, nsel,imol);
		w_t[j] = w[j];
		if(j==0 && F == 'O')
		{
			printf("w(old)%d:%.2f\n",imol,w[j]);
			fprintf(fplog,"w(old)%d:%.2f\n",imol,w[j]);
		}
		w[j] = exp(-beta*w[j]*1000/NA);
		sumw += w[j];
	}

	if(F == 'O')
		jbook = 0;
	else if(F == 'N')
	{
		jbook = Select(w,sumw);
		printf("w(new)%d:%.2f\n",imol,w_t[jbook]);
		fprintf(fplog,"w(new)%d:%.2f\n",imol,w_t[jbook]);
	}

	alpha = w[jbook]/sumw;
	printf("%c: alpha:%.5e   wi:%.2f  W:%.2f\n",F,alpha,w[jbook],sumw);
	fprintf(fplog,"%c: alpha:%.5e   wi:%.2f  W:%.2f\n",F,alpha,w[jbook],sumw);

	copy_array_2(select[jbook], cps[index_sel_cps[imol]]);
	if(F == 'N')
		copy_array_2(select[jbook], coordinate[index_sel[imol]]);

	if(F == 'N' && w[jbook] < 10)
		alpha = 0;

	return alpha;
}

void F_vector(int imol, int isite, int jmol, int jsite, double L, double coordinate[][3], double f_vector[])
{
	/*f_vector[fx, fy, fz, rx, ry, rz]   fij  */
	double rc,mag_lj=0.0,mag,d[3];

	rc = distance(coordinate[imol], coordinate[jmol], L);
	if(rc < RC_LJ/10)
		mag_lj = (L_A*pow(1/(10*rc),13)-L_B*pow(1/(10*rc),7))*(10e10);
	distance_pbc(coordinate[jmol], coordinate[imol], L, d);
	mag = -(mag_lj)*1000/NA;
	f_vector[0] = (d[0])*mag/rc;      /*unit: N*/
	f_vector[1] = (d[1])*mag/rc;
	f_vector[2] = (d[2])*mag/rc;
	f_vector[3] = d[0];      /*unit: nm*/
	f_vector[4] = d[1];
	f_vector[5] = d[2];
}

double Pressure(int N, double L, double beta, double coordinate[][3])
{
	int imol=0,isite=0,jmol=0,jsite=0;
	double vir = 0.0, fij[6];
	for(imol=0;imol<N-1;imol++)
		for(jmol=imol+1;jmol<N;jmol++)
			for(isite=0;isite<3;isite++)
				for(jsite=0;jsite<3;jsite++)
				{
					F_vector(imol, isite, jmol, jsite, L, coordinate, fij);
					vir += fij[0]*fij[3]+fij[1]*fij[4]+fij[2]*fij[5];     /*fij*rij*/
				}
	vir = vir/3;
	return (N*(10e9)/beta+vir)*(10e9)*(10e9)/(L*L*L)/101325;  /*unit: atm*/
}

void output_coord(double coordinate[][3], double box[], char address1[], char Title[], int natom)
{
	int i;
	char residuename[] = "W";
	char atomname = 'W';
	FILE *fp1;  /*fp1 generate a gro file, fp2 generate a psf file*/
	fp1 = fopen(address1, "w");
	fprintf(fp1,"%s\n", Title);
	fprintf(fp1,"%5d\n", natom);

	for (i=0; i<natom; i++)
	{
		fprintf(fp1,"     %-5s%5c%5d%8.3lf%8.3lf%8.3lf\n",residuename,atomname,i+1,coordinate[i][0],coordinate[i][1],coordinate[i][2]);
	}
	fprintf(fp1,"%10.5lf%10.5lf%10.5lf\n",box[0],box[1],box[2]);
	fclose(fp1);
}

void output_utot(char address1[], char address2[], double u[],int n, int N)
{
	int i;
	FILE *fp1,*fp2;
	fp1 = fopen(address1, "w");
	fp2 = fopen(address2, "w");
	fprintf(fp1,"@   title\"Total Energy\"\n");
	fprintf(fp1,"@   xaxis label \"Index of output\"\n");
	fprintf(fp1,"@   yaxis label \"Utot(*e4 kJ/mol)\"\n");
	fprintf(fp1,"@TYPE xy\n\n");
	fprintf(fp2,"@   title\"Molecular Energy\"\n");
	fprintf(fp2,"@   xaxis label \"Index of output\"\n");
	fprintf(fp2,"@   yaxis label \"Umol(kJ/mol)\"\n");
	fprintf(fp2,"@TYPE xy\n\n");
	for(i=0;i<n;i++)
	{
		fprintf(fp1,"%10d%10.3f\n",i,u[i]/10000);
		fprintf(fp2,"%10d%10.3f\n",i,u[i]/N);
	}
	fclose(fp1);
	fclose(fp2);
}

void output_pressure(char address[], double p[],int n)
{
	int i;
	FILE *fp;
	fp = fopen(address, "w");
	fprintf(fp,"@   title\"Pressure of System\"\n");
	fprintf(fp,"@   xaxis label \"Index of output\"\n");
	fprintf(fp,"@   yaxis label \"Pressure(atm)\"\n");
	fprintf(fp,"@TYPE xy\n\n");
	for(i=0;i<n;i++)
		fprintf(fp,"%10d%10.3f\n",i,p[i]/10000);
	fclose(fp);
}

void output_rdf(double coordinate[][3], int N, double L, char directory_address[])
{
	int i,ig, nhis,imol,jmol;
	double delg = 0.01 ,r,vb,nid_WW;
	double gWW[10000];
	char address1[100];
	FILE *fp1;

	sprintf(address1,"%srdfWW.xvg",directory_address);
	fp1 = fopen(address1,"w");
	fprintf(fp1,"@   title\"Radial distribution\"\n");
	fprintf(fp1,"@   xaxis label \"r\"\n");
	fprintf(fp1,"@   yaxis label \"\"\n");
	fprintf(fp1,"@TYPE xy\n");
	fprintf(fp1,"@ subtitle \"W - W\"\n");

	nhis = (int)((L/2)/delg)+1;
	for(i=0;i<nhis;i++)
	{
		gWW[i] = 0;
	}
	for(imol=0;imol<N-1;imol++)
		for(jmol=imol+1;jmol<N;jmol++)
		{
			r = distance(coordinate[imol],coordinate[jmol],L);
			if(r<L/2)
			{
				ig = (int)(r/delg);
				gWW[ig] += 2;
			}
		}
	for(i=0;i<nhis;i++)
	{
		r = delg*(i+0.5);
		vb = (pow(i+1,3)-pow(i,3))*pow(delg,3);
		nid_WW = (4./3.)*PI*vb*(N/pow(L,3));
		gWW[i]=gWW[i]/(N*nid_WW);
		fprintf(fp1,"%10.3f%10.4f\n",r,gWW[i]);
	}
	fclose(fp1);
}

void output_report(char address[], int N, int steps_energy, double u[], int steps_acc, double acc[])
{
	int i,count=0;
	double utot_avg = 0,umol_avg = 0, acc_avg = 0, acc_rmsd, utot_rmsd,umol_rmsd,max=-10e6,min=0;
	FILE *fp;


	fp = fopen(address, "w");
	fprintf(fp,"*******************FINAL REPORT*******************\n\n");

	for(i=0;i<steps_energy;i++)
	{
		utot_avg += u[i];
		if(u[i] > max)
			max = u[i];
		if(u[i] < min)
			min = u[i];
		umol_avg += u[i]/N;
	}
	utot_avg /= steps_energy;
	umol_avg /= steps_energy;

	utot_rmsd = 0;
	umol_rmsd = 0;
	for(i=0;i<steps_energy;i++)
	{
		utot_rmsd += (u[i]-utot_avg)*(u[i]-utot_avg);
		umol_rmsd += (u[i]/N-umol_avg)*(u[i]/N-umol_avg);
	}
	utot_rmsd /= steps_energy;
	umol_rmsd /= steps_energy;
	utot_rmsd = sqrt(utot_rmsd);
	umol_rmsd = sqrt(umol_rmsd);

	fprintf(fp,"-----------  POTENTIAL  -----------\n");
	fprintf(fp,"molecular potential:  max:%f  min:%f  avg:%f  rmsd:%f\n",max/N,min/N,umol_avg,umol_rmsd);
	fprintf(fp,"    total potential:  max:%f  min:%f  avg:%f  rmsd:%f\n",max,min,utot_avg,utot_rmsd);

	min = 1;
	max = 0;
	for(i=steps_acc/2;i<steps_acc;i++)
	{
		acc_avg += acc[i];
		if(acc[i] > max)
			max = acc[i];
		if(acc[i] < min)
			min = acc[i];
		count ++;
	}
	acc_avg /= count;

	acc_rmsd = 0;
	for(i=steps_acc/2;i<steps_acc;i++)
		acc_rmsd += (acc[i]-acc_avg)*(acc[i]-acc_avg);

	acc_rmsd /= count;
	acc_rmsd = sqrt(acc_rmsd);


	fprintf(fp,"--------------  ACC  --------------\n");
	fprintf(fp,"acceptance probability:  max:%f  min:%f  avg:%f  rmsd:%f\n",max,min,acc_avg,acc_rmsd);
	fprintf(fp,"##Note: just use the data from the second half of outputs##\n");

	fclose(fp);
}

/********************************************
**********MAIN PROGRAM STARTS HERE***********
********************************************/

int main(void)
{
	/*N:moleclue number; nmove:number of moves in trajectory; fout:frequency of output calculation
	 * L:box dimension; T:temperature; max_angle:max angle of rotation move(unit:angle! not radial);
	 *  coordinate:[molecule][site][x,y,z,charge] is used to store the  current new structure
	 * temp[number of molecule][site][x,y,z,charge] is used to store the  structure
	 *  for the last move; save_ener[nmove][1] is used to store the total energy*/
	int N=50000 , nmove=100, fout=1, fout_acc, acc_count = 0;
	int n, imol=0;         /*loop control*/
	int nsel, k, num_cps, rand_seed;
	int index_cps[N], index_sel[100], index_sel_cps[100];
	double L, T, beta, coordinate[N][3], temp[N][3],cps[N][3], temp_cps[N][3], center[3], box[3], sys_para[3], sim_para[5];
	double duration,utot[100000],upress[100000],save_acc[1500],acc;
	double alpha_n, alpha_o, uact_n, uact_o, r, r_cav;
	char Title[100], address_out1[100],address_out3[100], address_out4[100],address_out5[100],address_in[100], address_log[100], address_report[100], address_acc[100];
	char directory_in[100] = "/Users/Ziwei/Desktop/Mar_04_001/", open_name[100] = "input.txt", InputConfigName[100], directory_out[100] = "";
	FILE *fp_log;
	clock_t start, finish;

	start = clock();
	printf("ProName:Domer_v1.0  Job starts\n");
	strcpy(directory_out,directory_in);
	n = frand();    /*discard the first random*/

	sprintf(address_log, "%slog.txt", directory_out);
	fp_log = fopen(address_log, "w");
	fprintf(fp_log,"ProName:Domer_v1.0     JOB STARTS\n");

	sprintf(address_in, "%s%s", directory_in, open_name);
	read_input(address_in, InputConfigName, sys_para, sim_para);

    shift_potential(6, L6);
	shift_potential(12, L12);

	N = (int)sys_para[0];
	L = sys_para[1];
	T = sys_para[2];
	nmove = (int)sim_para[0];
	fout = (int)sim_para[1];
	r_cav = sim_para[2]/10;
	k = (int)sim_para[3];
	rand_seed = (int)sim_para[4];
	beta = 1.0/(1.3806488*pow(10,-23)*T);  /*unit for Boltzmann constant is J/K */

	for(n=0;n<rand_seed;n++)
		r = frand();

	fprintf(fp_log,"N:%d L:%.5f T:%.2f \nnmove:%d fout:%d r_cav:%.3f k:%d\n\n",N,L,T,nmove,fout,r_cav,k);
	printf("N:%d L:%.5f T:%.2f \nnmove:%d fout:%d r_cav:%.3f k:%d\n\n",N,L,T,nmove,fout,r_cav,k);

	sprintf(address_in, "%s%s", directory_in, InputConfigName);
	printf("%s\n",address_in);

	read_coordinate(coordinate, box, address_in, Title, N);
	utot[0] = Utotal(coordinate, N, L);
	upress[0] = Pressure(N, L, beta, coordinate);
	sprintf(address_out3, "%stotal_energy.xvg", directory_out);
	sprintf(address_out4, "%smol_energy.xvg", directory_out);
	sprintf(address_out5, "%spressure.xvg", directory_out);
	sprintf(address_report, "%sREPORT.txt", directory_out);
	sprintf(address_acc, "%sacc.xvg", directory_out);

	for(imol=0;imol<N;imol++)
			nojump(L, coordinate[imol]);

	fprintf(fp_log,"*******************The original Total Potential is %.3f******************\n", Utotal(coordinate, N, L));

	copy_array(coordinate,temp,N);

	for(n=1;n<(nmove+1);n++)
	{
		fprintf(fp_log,"***************************************************step:%d*************************************************\n",n);
		nsel = 0;
		num_cps = 0;
		alpha_n = 1;
		alpha_o = 1;
		uact_n = 0;
		uact_o = 0;

		center[0] = L*(frand()-0.5);
		center[1] = L*(frand()-0.5);
		center[2] = L*(frand()-0.5);

		for(imol=0;imol<N;imol++)
		{
			r = distance(coordinate[imol],center,L);
			if(r < (r_cav+RC_LJ/10))
			{
				index_cps[num_cps] = imol;
				copy_array_2(coordinate[imol],cps[num_cps]);
				num_cps++;
				if(r < r_cav)
				{
					index_sel[nsel] = imol;
					index_sel_cps[nsel] = num_cps-1;
					nsel++;
				}
			}
		}
		copy_array(cps,temp_cps,num_cps);

		printf("nsel:%d\n",nsel);
		fprintf(fp_log,"nsel:%d\n",nsel);

		for(imol=0;imol<nsel;imol++)
			alpha_o *= gen_config('O', coordinate, cps, temp_cps, r_cav, num_cps, index_sel, index_sel_cps, nsel, k, beta, center, L, imol, fp_log);

		printf("\n");
		fprintf(fp_log,"\n");

		for(imol=0;imol<nsel;imol++)
			alpha_n *= gen_config('N', coordinate, cps, temp_cps, r_cav, num_cps, index_sel, index_sel_cps, nsel, k, beta, center, L, imol, fp_log);

		for(imol=0;imol<nsel;imol++)
		{
			uact_n += Usite_cav(cps, index_sel_cps[imol], num_cps, L,index_sel_cps, nsel, imol);
			uact_o += Usite_cav(temp_cps, index_sel_cps[imol], num_cps, L,index_sel_cps, nsel, imol);
		}

		/*start to calculate the acc*/
		acc = (alpha_o/alpha_n)*(exp(-beta*(uact_n-uact_o)*1000/NA));

		if(nsel == 0 || alpha_n < 10e-40)
			acc = 0;

		if(frand() < min(1,acc))
		{
			for(imol=0;imol<nsel;imol++)
				copy_array_2(coordinate[index_sel[imol]], temp[index_sel[imol]]);
			acc_count += 1;
			fprintf(fp_log,"ACCEPT\n");
			printf("ACCEPT  ");
		}
		else
		{
			for(imol=0;imol<nsel;imol++)
				copy_array_2(temp[index_sel[imol]], coordinate[index_sel[imol]]); /*set coordinate back to the last one*/
			fprintf(fp_log,"REJECT\n");
			printf("REJECT  ");
		}
		if(n%fout==0)
		{
			sprintf(address_out1, "%swater_out_%d.gro", directory_out,n/fout);
			utot[n/fout] = Utotal(coordinate, N, L);
			upress[n/fout] = Pressure(N, L, beta, coordinate);
			fprintf(fp_log,"The %dth total U is %.3fkJ/mol\n", n/fout, utot[n/fout]);
			fprintf(fp_log,"The %dth system pressure is %.3fatm\n", n/fout, upress[n/fout]);
			output_coord(temp, box, address_out1, Title, N);
		}

		printf("n_sel:%d	",nsel);
		printf("a_n:%.5e	",alpha_n);
		printf("a_o:%.5e	",alpha_o);
		printf("uact_n：%.3f	",uact_n);
		printf("uact_o：%.3f",uact_o);
		printf("acc: %.5f\n",acc);

		fprintf(fp_log,"n_sel:%d	",nsel);
		fprintf(fp_log,"a_n:%.5e	",alpha_n);
		fprintf(fp_log,"a_o:%.5e	",alpha_o);
		fprintf(fp_log,"uact_n：%.3f	",uact_n);
		fprintf(fp_log,"uact_o：%.3f",uact_o);
		fprintf(fp_log,"acc: %.5f\n",acc);

		printf("The %dth trial is performed\n",n);
		printf("prob %.4f\n\n**************************\n",((double)acc_count/(double)n));
	}
	output_utot(address_out3, address_out4, utot, nmove/fout+1, N);
	output_pressure(address_out5, upress, nmove/fout+1);
	output_rdf(coordinate, N, L, directory_out);
	output_report(address_report, N, nmove/fout+1, utot, nmove/fout_acc+1,  save_acc);

	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	fprintf(fp_log,"\nJob ended, it lasts %.2f seconds, acceptance probability is %.4f\n",duration,((double)acc_count/(double)nmove));
	printf("Job quits normally\n");
	fclose(fp_log);
}


