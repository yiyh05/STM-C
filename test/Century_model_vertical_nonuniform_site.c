/* 
    Solve the SOC vertical profile based on the Century SOM module
    with prescribed litterfall and soil temperature;
    for a single-pixel or site.  
*/

#include  "/home/yonghong.yi/Clib/include/function.h"
#include  "/home/yonghong.yi/Clib/include/constants.h"  

/* there are some missing data from 2001-2002 in MODIS GPP... */
/* looks like I have not produced 0.5 deg MODIS GPP */
#define year_start   2008
#define year_end     2010

#define NCLASS  9     /* total 9 land cover classes */

#define NCOMP  46

/* the structure of SOC pools */
#define nlayer 15     /* the soil layers */
#define npool  7      /* the number of pools */

typedef struct bplut
{
	/* parameters for the soil decomposition model */
    double  fLabile; 
	double  fCellulose;
	double  fLignin;
	double  fwoody; 
}BPLUT;

/* the soil carbon model function */
void Cal_TWConstraint(float tsoil, float *twmult_out);
void soc_model(float32 *NPP_arry, float32 *TWmult_arry, float32 *soc_arry, float32 *Rh_out, BPLUT *bplut, float32 *AD);
void CalSOCContent(float32 *SOC_arry, float *totsoc);
void TDMA(int ndim, double *a, double *b, double *c, double *d, double *x);

/* ancillary function */
void FILL_BPLUT(char* BPLUT_file_name, BPLUT* bplut);
int ReadTsoilData(char *filepath, char *var_name, char *site_name, int year, float32 *data_arry);
int ReadTsoilData_v2(char *filepath, char *var_name, char *site_name, int year, float32 *data_arry);
int ReadCFluxData(char *filepath, char *var_name, char *site_name, int year, int var_index, float32 *data_arry);
int ReadTsurfData(char *filepath, char *var_name, char *site_name, int year, int var_index, float32 *data_arry);
void WriteSOCData(char *filepath, char *, char *, int year, int max_iter, float32 *data_arry);
void WriteCFluxData(char *filepath, char *, char *var_name, int year, float32 *NEE_arry, float32 *NPP_arry, float32 *Rh_arry);

int main(int argc, char** argv)
{
	char filename[255];
	char temp_char[255];
	char year_char[5];
	char *pch;

	int i, ic, iyr, icell;
	int imth, ipool;
	int flag;
	int nmonth;
	int nyear;
	int iter, max_iter;
	int id, ilayer;

	float32 *Rh_arry;
	float32 *Tsoil_arry;
	float32 *TWmult_arry;

	float32 NPP_arry[NCOMP];
	float32 NEE_arry[NCOMP];

	float32 *SOC_out;
	float32 *SOC_arry;

	float tsoil, twmult;
	float Rh, totsoc0, totsoc1;
	float32 AD[npool];   /* accelerating rates */

	/* the soil depth and node */
	float znode[nlayer] = {0.01,0.03,0.08,0.13,0.23,0.33,0.45,0.55,0.7,1.05,1.4,1.75,2.25,2.75,3.25};
	float dzmm[nlayer];

	BPLUT bplut[NCLASS];

	FILE *sfp_input;
	FILE *sfp_output;

	if(argc != 5)
	{
		printf("exe: <data direct> <output filepath> <site name> <landcover type>!\n");
		exit(-1);
	}

	/* land cover type: open water(1), Ice/snow(2), Developed(3), Barren land(4), 
	   Forest(5), Shrub(6), Grassland/tundra(7), Cultivated(8), Wetlands(9) */
	ic = atoi(argv[4]);
	if((ic<=4) || (ic>9))
	{
		printf("unknown vegetation types!\n");
		exit(-1);
	}

	/* used to test the algorithm */
	/*for(ilayer=0; ilayer<nlayer; ilayer++)
	{
		znode[ilayer] = (ilayer+0.5)*0.1;
	}*/

	/* the depth of each layer; is it consistent with the Tsoil layer allocation?? */
	for(ilayer=nlayer-1; ilayer>=0; ilayer--)
	{
		if(ilayer>0)
			dzmm[ilayer] = znode[ilayer] - znode[ilayer-1];
		else
			dzmm[ilayer] = znode[ilayer+1] - znode[ilayer]; /* a virtual node in the surface*/

		printf("layer: %d, soil depth: %f\n", ilayer+1, dzmm[ilayer]);
	}

	/* Fill the BPLUT */
	filename[0] = 0;
        strcat(filename, argv[2]);
        strcat(filename, "som_bplut.txt");
        FILL_BPLUT(filename, bplut);	

	/* initialization */
	nyear = year_end - year_start + 1;

	/* for data reading */
	Tsoil_arry = (float32*)calloc(nlayer*NCOMP, sizeof(float32)); /* soils below */

	/* temporal arrays*/
	TWmult_arry = (float32*)calloc(nlayer*NCOMP, sizeof(float32));

	/* output data */
	Rh_arry = (float32*)calloc(NCOMP*nlayer, sizeof(float32));
	
	/* SOC pools include: 3  litter pools, CWD pool, and 3 SOM pools; allocated for each depth  */
	max_iter = 50000;
	SOC_out = (float32*)calloc(npool*nlayer*max_iter, sizeof(float32));
	/* SOC pools at a specific time step */
	SOC_arry = (float32*)calloc(npool*nlayer, sizeof(float32));

	for(i=0; i<npool*nlayer; i++)
	{
		SOC_arry[i] = 0.;
	}

	/* read the NPP dataset */
	for(iyr=year_start; iyr<=year_end; iyr++)
	{
		/* initialization */
		for(id=0; id<NCOMP; id++)
		{
			NPP_arry[id] = Nodata;
			NEE_arry[id] = Nodata;
		}

		for(i=0; i<nlayer*NCOMP; i++)
		{
			Tsoil_arry[i] = Nodata;
		}

		for(i=0; i<nlayer*NCOMP; i++)
		{
			TWmult_arry[i] = Nodata;
		}

		for(i=0; i<npool*nlayer*max_iter; i++)
		{
			SOC_out[i] = Nodata;
		}

		/* Read NPP data */		
		flag = ReadCFluxData(argv[1], "C_fluxes.", argv[3], iyr, 6, NPP_arry);
		if(flag == 0)
		{
			printf("NPP < 0.!\n");
			continue;
		}

		/* Read the Tsoil data */
		flag = ReadTsoilData(argv[1], "soilt.", argv[3], iyr, Tsoil_arry); 
		/* flag = ReadTsoilData_v2(argv[1], ".Tsoil_interp.", argv[3], iyr, Tsoil_arry); */
		if(flag == 0)
		{
			printf("soil temperature <= Nodata\n");
			continue;
		}

		/* calculate the TWmult */
		for(id=0; id<NCOMP; id++)
		{
			for(ilayer=0; ilayer<nlayer; ilayer++)
			{
				tsoil = Tsoil_arry[ilayer*NCOMP+id];

				twmult = Nodata;
				Cal_TWConstraint(tsoil, &twmult);
				TWmult_arry[ilayer*NCOMP+id] = twmult;
			}
		}

		if(iyr==year_start)
			max_iter = 50000;
		else
			max_iter = 0;

		/* the total SOC data */
		totsoc0 = Nodata;
		CalSOCContent(SOC_arry, &totsoc0);
		
		for(ipool=0; ipool<npool; ipool++)
		{
			AD[ipool] = 1.;
		}

		iter = 0;
		if(max_iter > 1) /* spin up the SOC pool */
		{
			/* first spin-up the model using accelerating rates */
			AD[4] = 5.;    /* the slow pool accelerating rates */
			AD[5] = 100.;  /* the passive pool */
			AD[6] = 1.;    /* the cwd pool */

			while(iter < max_iter)
			{			
				soc_model(NPP_arry, TWmult_arry, SOC_arry, Rh_arry, &bplut[ic-1], AD);

				for(ipool=0; ipool<npool; ipool++)
				for(ilayer=0; ilayer<nlayer; ilayer++)			
				{
					SOC_out[iter*npool*nlayer+ipool*nlayer+ilayer] = SOC_arry[ipool*nlayer+ilayer];
					if(SOC_arry[ipool*nlayer+ilayer] < 0.)
					{
						printf("soc < 0.!\n");
						exit(-1);
					}
				}

				/* check whether it should stop spin-up process */
				totsoc1 = Nodata;
				CalSOCContent(SOC_arry, &totsoc1);
			 
				if(fabs(totsoc1 - totsoc0) < 0.5)  /* unit: gC/m2 */
				{
					break;
				}

				totsoc0 = totsoc1;
			
				iter ++;
			} /* while(iter < max_iter) */

			/* then scale the SOC_arry */
			for(ipool=0; ipool<npool; ipool++)
			{
				for(ilayer=0; ilayer<nlayer; ilayer++)
				{
					SOC_arry[ipool*nlayer+ilayer]  = SOC_arry[ipool*nlayer+ilayer] * AD[ipool];
				}
			}

			totsoc0 = Nodata;
			CalSOCContent(SOC_arry, &totsoc0);	

			/* then do normal spinup process */
			for(ipool=0; ipool<npool; ipool++)
			{
				AD[ipool] = 1.;
			}

			while(iter < max_iter)
			{
				soc_model(NPP_arry, TWmult_arry, SOC_arry, Rh_arry, &bplut[ic-1], AD);

				for(ipool=0; ipool<npool; ipool++)
				for(ilayer=0; ilayer<nlayer; ilayer++)			
				{
					SOC_out[iter*npool*nlayer+ipool*nlayer+ilayer] = SOC_arry[ipool*nlayer+ilayer];
					if(SOC_arry[ipool*nlayer+ilayer] < 0.)
					{
						printf("soc < 0.!\n");
						exit(-1);
					}
				}

				/* check whether it should stop spin-up process */
				totsoc1 = Nodata;
				CalSOCContent(SOC_arry, &totsoc1);
			 
				if(fabs(totsoc1 - totsoc0) < 0.1)  /* unit: gC/m2 */
				{
					break;
				}

				totsoc0 = totsoc1;
			
				iter ++;
			} /* while(iter < max_iter) */
		} /* if(max_iter > 1)  */
		
		/* followed by transit run */
		soc_model(NPP_arry, TWmult_arry, SOC_arry, Rh_arry, &bplut[ic-1], AD);
		
		for(ipool=0; ipool<npool; ipool++)
		for(ilayer=0; ilayer<nlayer; ilayer++)			
		{
			SOC_out[iter*npool*nlayer+ipool*nlayer+ilayer] = SOC_arry[ipool*nlayer+ilayer];
			if(SOC_arry[ipool*nlayer+ilayer] < 0.)
			{
				printf("soc < 0.!\n");
				exit(-1);
			}
		}

		/* the C flux */
		for(id=0; id<NCOMP; id++)
		{
			Rh = 0.;
			for(ilayer=0; ilayer<nlayer-1; ilayer++)
			{
				if(ilayer ==0)
					Rh = Rh + dzmm[ilayer] * Rh_arry[id*nlayer+ilayer];
				else
					Rh = Rh + dzmm[ilayer+1] * (Rh_arry[id*nlayer+ilayer] + Rh_arry[id*nlayer+ilayer+1])/2.;
			}

			if(Rh < 0.)
			{
				printf("Rh < 0.!\n");
				exit(-1);
			}

			NEE_arry[id] = Rh - NPP_arry[id];
		}


		/* output the data */
		WriteSOCData(argv[2], argv[3], "soc_spinup.", iyr, 50000, SOC_out);
		WriteCFluxData(argv[2], argv[3], "CFlux_mod.", iyr, NEE_arry, NPP_arry, Rh_arry);

		printf("finish processing year %d\n", iyr);		
	} /* loop for each year */	
	free(SOC_out);	
	free(SOC_arry);	
	free(Tsoil_arry);
	free(TWmult_arry);
	free(Rh_arry);
	
	return 1;
}

void soc_model(float32 *NPP_arry, float32 *TWmult_arry, float32 *soc_arry, float32 *Rh_out, BPLUT *bplut, float32 *AD)
{
	int i, imth;
	int ipool;
	int id, num;
	int ilayer;
	int numdays;	
	int iday;
	int ndim = nlayer;

	float ann_npp;
	float litterfall, litter_above;
	float litter_below, litter_wood;
	float litter, influx, fs;
	float Fs, Fd, Fs_wood, Fd_wood;
	float cwd_pool, lit1_pool;
	float lit2_pool, lit3_pool;
	float Rh_cwd;

	/* temporay array */
	float  TWmult[nlayer], Rh_arry[nlayer];
	float  Rh_pool[npool], krate[npool];
	float  frac_co2[npool];

	/* to solve the tridiagonal matrix */
	double ATR[nlayer],BTR[nlayer],CTR[nlayer];
	double DTR[nlayer],XTR[nlayer];
	double *ATR_arry;
	double *BTR_arry;
	double *CTR_arry;
	double *DTR_arry;

	float fwood_cellu = 0.7;
	float fwood_lignin = 0.3;

	float fsom3_co2 = 1.0 ;/* which is the fraction of respired carbon */

	/* the diffusion rates in permafrost areas, ranging from 1-5 cm2 yr-1 */
	float  D0 = 20.;  /* cm2 yr-1*/
	float  D0_daily;
	double D_daily[nlayer];

	/* the percentage of aboveground litterfall */
	float fabove = 0.4; 
	float zr = 0.25;   /* m, the efolding depth of root litter */
	float zs = 0.1; /*m, the efolding depth of surface litter */

	/* constant for the matrix computation */
	float delt_t  = 1.;   /* in days */
	float delt_z2_up, delt_z2_down;

	/* the soil depth and node */
	float znode[nlayer] = {0.01,0.03,0.08,0.13,0.23,0.33,0.45,0.55,0.7,1.05,1.4,1.75,2.25,2.75,3.25}; 
	float dzmm[nlayer];

	/* allocate memory for this */
	ATR_arry = (double *)calloc(nlayer*npool, sizeof(double));
	BTR_arry = (double *)calloc(nlayer*npool, sizeof(double));
	CTR_arry = (double *)calloc(nlayer*npool, sizeof(double));
	DTR_arry = (double *)calloc(nlayer*npool, sizeof(double));

	/* the maximum decomposition rates*/
	krate[0] = 0.040303;
	krate[1] = 0.01064;
	krate[2] = 0.01064;
	krate[3] = 0.015647;
	krate[4] = 0.000443;
	krate[5] = 0.00001;
	krate[6] = 0.000659;

	/* accelerate the SOC decomposition and diffusion rates as well */
	for(ipool=4; ipool<npool; ipool++)
	{
		krate[ipool] = krate[ipool] * AD[ipool];
	}

	/* the fraction of respiration lost as CO2  */
	frac_co2[0] = 0.55; /*litterpool 1*/
	frac_co2[1] = 0.5;
	frac_co2[2] = 0.5;
	frac_co2[3] = 0.4;
	frac_co2[4] = 0.55;
	frac_co2[5] = fsom3_co2;
	frac_co2[6] = 0.; /* no CO2 respiration loss for the CWD pool */

	/* temporary */
	/*for(ilayer=0; ilayer<nlayer; ilayer++)
	{
		znode[ilayer] = (ilayer+0.5)*0.1;
	} */

	/* the depth of each layer */
	for(ilayer=nlayer-1; ilayer>=0; ilayer--)
	{
		if(ilayer>0)
			dzmm[ilayer] = znode[ilayer] - znode[ilayer-1];
		else
			dzmm[ilayer] = znode[ilayer+1] - znode[ilayer]; /* a virtual node in the surface*/
	}

	/* the diffusion rates converted to m2 day-1 */
	/* D0_daily = D0/(365.*100.*100.); */
	D0_daily = D0;

	/* modify the diffusion rates for each depth */
	for(ilayer=0; ilayer<nlayer; ilayer++)
	{
		/* set the threshold as 3.0m */
		if(znode[ilayer] <= 1.0)
			D_daily[ilayer] = D0_daily;
		else if(znode[ilayer] >= 3.0)
			D_daily[ilayer] = 0.;
		else
			D_daily[ilayer] = D0_daily*(1.-(znode[ilayer] - 1.0)/(2*1.0));

		/* printf("ilayer: %d, depth: %f, D_daily: %lf\n", ilayer+1, znode[ilayer], D_daily[ilayer]); */

		if(D_daily[ilayer] < 0.)
		{
			printf("D_daily<0.!\n");
			exit(-1);
		}
	}

	ann_npp = 0.;
	for(id=0; id<NCOMP; id++)
	{
		ann_npp += NPP_arry[id];
	}

	if(ann_npp <= 0.)
	{
		printf("ann_npp: %f\n", ann_npp);
		exit(-1);
	}

	/* daily litterfall */
	litterfall   = ann_npp/365.;
	litter_above = litterfall * fabove;
	litter_below = litterfall - litter_above;
	if((litter_above<0.) || (litter_below<0.))
	{
		printf("litter<0.!\n");
		exit(-1);
	}

	/* make sure there is no missing data in the input data  */
	for(i=0; i<nlayer*npool; i++)
	{
		if(soc_arry[i] <0.)
		{
			printf("soc: %f\n", soc_arry[i]);
			exit(-1);
		}
	}

	for(i=0; i<nlayer*NCOMP; i++)
	{
		if((TWmult_arry[i] < 0.) || (TWmult_arry[i] > 1.))
		{
			printf("TWmult: %f\n", TWmult_arry[i]);
			exit(-1);
		}
	}

	for(i=0; i<NCOMP*nlayer; i++)
	{
		Rh_out[i] = 0.;
	}

	for(id=0; id<NCOMP; id++)
	{
		for(ilayer=0; ilayer<nlayer; ilayer++)
		{
			TWmult[ilayer] = TWmult_arry[ilayer*NCOMP+id];
		}

		if(id<NCOMP-1)
		{
			numdays = 8;
		}
		else
		{
			numdays = 5;
		}

		for(iday=0; iday<numdays; iday++)
		{
			/* initialize the arrays */
			for(i=0; i<nlayer*npool; i++)
			{
				ATR_arry[i] = Nodata;
				BTR_arry[i] = Nodata;
				CTR_arry[i] = Nodata;
				DTR_arry[i] = Nodata;
			}

			/************** the soil at depth (0-3m) ****************************/
			for(ilayer=0; ilayer<nlayer; ilayer++)
			{
				/* initialization */
				for(ipool=0; ipool<npool; ipool++)
				{
					Rh_pool[ipool] = 0.;
				}

				Rh_arry[ilayer] = 0.;

				/* aboveground and belowground litterfall, gC/m3 for this depth*/
				Fd = litter_below/zr * exp(-1.*znode[ilayer]/zr);
				Fs = litter_above/zs * exp(-1.*znode[ilayer]/zs);
				Fd_wood = Fd * bplut->fwoody;
				Fs_wood = Fs * bplut->fwoody;

				/* the CWD pool */
				cwd_pool = soc_arry[nlayer*6+ilayer];
				cwd_pool += (Fd_wood + Fs_wood);
				Rh_cwd    = cwd_pool * krate[6] * TWmult[ilayer];	
				cwd_pool -= Rh_cwd;
				soc_arry[nlayer*6+ilayer] = cwd_pool;

				/* respiration from the litterfall and soc pools */
				for(ipool=0; ipool<6; ipool++)
				{
					Rh_pool[ipool] = soc_arry[ipool*nlayer+ilayer] * krate[ipool] * TWmult[ilayer];
				}

				/* the total soil respiration at this soil depth */
				for(ipool=0; ipool<6; ipool++)
				{
					Rh_arry[ilayer] += Rh_pool[ipool] * frac_co2[ipool];
				}

				for(ipool=0; ipool<6; ipool++)
				{
					if(ilayer == 0) /* the surface layer */
					{
						delt_z2_up = 2./(dzmm[ilayer+1]*(dzmm[ilayer]+dzmm[ilayer+1]));

						ATR_arry[ipool*nlayer+ilayer] = 0.;
						BTR_arry[ipool*nlayer+ilayer] = 1. + 2.*D_daily[ilayer]*AD[ipool]/(365.*100.*100.)*delt_t*delt_z2_up;
						CTR_arry[ipool*nlayer+ilayer] = -2.*D_daily[ilayer]*AD[ipool]/(365.*100.*100.)*delt_t*delt_z2_up;
					}
					else if(ilayer < nlayer-1)
					{
						delt_z2_up = 2./(dzmm[ilayer+1]*(dzmm[ilayer]+dzmm[ilayer+1]));
						delt_z2_down = 2./(dzmm[ilayer]*(dzmm[ilayer]+dzmm[ilayer+1]));

						ATR_arry[ipool*nlayer+ilayer] = -1.* D_daily[ilayer]*AD[ipool]/(365.*100.*100.) *delt_t*delt_z2_down;
						BTR_arry[ipool*nlayer+ilayer] =  1. + D_daily[ilayer]*AD[ipool]/(365.*100.*100.) *delt_t*(delt_z2_up+delt_z2_down);
						CTR_arry[ipool*nlayer+ilayer] = -1.*D_daily[ilayer]*AD[ipool]/(365.*100.*100.) * delt_t*delt_z2_up;

					}
					else   /* the bottom layer */
					{
						delt_z2_down = 2./(dzmm[ilayer]*(dzmm[ilayer]+dzmm[ilayer]));

						ATR_arry[ipool*nlayer+ilayer] = -2.*D_daily[ilayer]*AD[ipool]/(365.*100.*100.) * delt_t * delt_z2_down;
						BTR_arry[ipool*nlayer+ilayer] = 1. + 2.*D_daily[ilayer]*AD[ipool]/(365.*100.*100.) * delt_t * delt_z2_down; 
						CTR_arry[ipool*nlayer+ilayer] = 0.;
					}

					if(ipool == 0) /* litterfall pool 1*/
					{
						influx = (Fs+Fd-Fs_wood-Fd_wood) * bplut->fLabile;
					}
					else if(ipool == 1) /* litterfall pool 2*/
					{
						influx = (Fs+Fd-Fs_wood-Fd_wood) * bplut->fCellulose + Rh_cwd * fwood_cellu;
					}
					else if(ipool == 2) /* litterfall pool 3 */
					{
						influx = (Fs+Fd-Fs_wood-Fd_wood) * bplut->fLignin + Rh_cwd * fwood_lignin;
					}
					else if(ipool==3) /* the active SOM pool */
					{
						influx = Rh_pool[0]*(1.-frac_co2[0]) + Rh_pool[1]*(1.-frac_co2[1]) + Rh_pool[4]*(1.-0.03-frac_co2[4]) + (1.-frac_co2[5])*Rh_pool[5]; 
					}
					else if(ipool == 4) /* the slow SOM pool */
					{
						influx = Rh_pool[2]*(1.-frac_co2[2]) + (1.-0.004-frac_co2[3])*Rh_pool[3]; 
					}
					else  /* the passive SOM pool */
					{
						influx = Rh_pool[3]*0.004 + Rh_pool[4]*0.03;
					}

					DTR_arry[ipool*nlayer+ilayer] = soc_arry[ipool*nlayer+ilayer] + (influx - Rh_pool[ipool]) * delt_t;
				} /* for(ipool=3; ipool<6; ipool++) */
			} /* for(ilayer=0; ilayer<nlayer; ilayer++) */

			for(ipool=0; ipool<6; ipool++)
			{
				for(ilayer=0; ilayer<nlayer; ilayer++)
				{
					ATR[ilayer] = ATR_arry[ipool*nlayer+ilayer];
					BTR[ilayer] = BTR_arry[ipool*nlayer+ilayer];
					CTR[ilayer] = CTR_arry[ipool*nlayer+ilayer];
					DTR[ilayer] = DTR_arry[ipool*nlayer+ilayer];

					XTR[ilayer] = Nodata;
				}

				TDMA(ndim, ATR, BTR, CTR, DTR, XTR);

				for(ilayer=0; ilayer<nlayer; ilayer++)
				{
					soc_arry[ipool*nlayer+ilayer] = XTR[ilayer];
					if(XTR[ilayer] < 0.)
					{
						printf("XTR < 0.: %f!\n", XTR[ilayer]);
						exit(-1);
					}
				}
			}

			/* update the C fluxes */
			for(ilayer=0; ilayer<nlayer; ilayer++)
			{
				if(Rh_arry[ilayer] >= 0.)
				{
					Rh_out[id*nlayer+ilayer] += Rh_arry[ilayer];
				}
				else
				{
					printf("Rh_arry (%f) < 0.!\n", Rh_arry[ilayer]);
					exit(-1);
				}
			}
		} /*  for(iday=0; iday<numdays; iday++) */
	} /* for(id=0; id<NCOMP; id++) */
}

void CalSOCContent(float32 *SOC_arry, float *totsoc)
{
	float znode[nlayer] = {0.01,0.03,0.08,0.13,0.23,0.33,0.45,0.55,0.7,1.05,1.4,1.75,2.25,2.75,3.25}; 
	float sum, soc;

	int ipool, ilayer;
	
	/* temporary */
	/*for(ilayer=0; ilayer<nlayer; ilayer++)
	{
		znode[ilayer] = (ilayer+0.5)*0.1;
	} */

	sum = 0.;
	for(ipool=0; ipool<npool; ipool++)
	for(ilayer=0; ilayer<nlayer-1; ilayer++)
	{
		soc = (SOC_arry[ipool*nlayer+ilayer] + SOC_arry[ipool*nlayer+ilayer+1])/2.;
		if(soc < 0.)
		{
			printf("SOC<0.!\n");
			exit(-1);
		}
		else
		{
			sum = sum + soc * (znode[ilayer+1] - znode[ilayer]);
		}
	}

	(*totsoc) = sum;
}

void TDMA(int ndim, double *a, double *b, double *c, double *d, double *x)
{
	/* use the Thomas algorithm to solve the tridiagonal algorithm (TDMA) */
	int i;
	double m;

	/* make sure the matrix is diagonally dominate; otherwise the solution is not stable */
	/*for(i=0; i<ndim; i++)
	{
		printf("a: %lf, b: %lf, c: %lf, d: %lf\n", a[i], b[i], c[i], d[i]);
	} */
	
	for(i=0; i<ndim; i++)
	{
		if( (fabs(b[i]) <= (fabs(a[i]) + fabs(c[i]))) || (d[i] < 0.))
		{
			printf("cannot find stable solution for that!\n");
			printf("a: %lf, b: %lf, c: %lf, d:%lf\n", a[i], b[i], c[i], d[i]);
			exit(-1);
		}

		x[i] = Nodata;
	}
	
	for(i=1; i<ndim; i++)
	{
		m = a[i]/b[i-1];
		b[i] = b[i] - m*c[i-1];
		d[i] = d[i] - m*d[i-1];
	}

	x[ndim-1] = d[ndim-1]/b[ndim-1];

	for(i=ndim-2; i>=0; i--)
	{
		x[i] = (d[i] - c[i]*x[i+1])/b[i];
	}

	for(i=0; i<ndim; i++)
	{
		if(x[i] < 0.)
		{
			printf("x[%d]: %lf\n", i+1, x[i]);
			exit(-1);
		}
	}
}


void Cal_TWConstraint(float tsoil, float *twmult_out)
{
	float Tmult, Wmult;
	float TWmult;

	/* temperature constraint parameters */
	float topt = 25.; 
	float a1, a2;
	
	/* initialize the TCF temperature parameters */
	/* the optimal air temperature (topt) only affects the coefficient "a2" */
	a1 = 308.56;
	a2 = 66.02+topt-20.;

	if(tsoil <= Nodata)
	{
		printf("temperature data invalid!\n");
		exit(-1);
	}

	if(tsoil >= topt)
		Tmult = 1.;
	else if(tsoil <= -10.)
		Tmult = 0.0;
	else
	{
		Tmult = exp(a1*(1./a2-1./(a2+tsoil-topt))); /* this is the right expression */
	}

	/* the unfrozen water content */
	if(tsoil >= 0.)
		Wmult = 1.0;
	else
	{
		Wmult = pow((0.1-tsoil)/0.01, -0.5)/pow(10.,-0.5);
	}

	TWmult = Tmult * Wmult;

	/* printf("tsoil: %f, tmult: %f, wmult: %f, TWmult: %f\n", tsoil, Tmult, Wmult, TWmult); */

	if((TWmult < 0.) ||(TWmult > 1.))
	{
		printf("tmult: %f\n", TWmult);
		exit(-1);
	}

	(*twmult_out) = TWmult;
}

void FILL_BPLUT(char* BPLUT_file_name, BPLUT* bplut)
{
/* land cover type: open water(1), Ice/snow(2), Developed(3), Barren land(4), 
   Forest(5), Shrub(6), Grassland/tundra(7), Cultivated(8), Wetlands(9) */
	char  temp_str[255];
	int   ic;

	FILE  *sfp_input;

	if((sfp_input = fopen(BPLUT_file_name, "r")) == NULL)
	{
	    printf("cannot open MOD17_BPLUT text file %s, exit...\n",BPLUT_file_name);
	    exit(0);
	}

	/* skip the first two lines */
	fgets(temp_str, 255, sfp_input);

	/* the labile fraction of leaf/fine root litterfall */
	fscanf(sfp_input, "%s", temp_str);
	for(ic=0; ic<NCLASS; ic++)
	{
		fscanf(sfp_input, "%lf", &bplut[ic].fLabile);
		printf("%lf, ", bplut[ic].fLabile);
	}
	printf("\n");

	/* the cellulose fraction of leaf/root litterfall */
	fscanf(sfp_input, "%s", temp_str);
	for(ic=0; ic<NCLASS; ic++)
	{
		fscanf(sfp_input, "%lf", &bplut[ic].fCellulose);
		printf("%lf, ", bplut[ic].fCellulose);
	}
	printf("\n");

	fscanf(sfp_input, "%s", temp_str);
	for(ic=0; ic<NCLASS; ic++)
	{
		fscanf(sfp_input, "%lf", &bplut[ic].fLignin);
		printf("%lf, ", bplut[ic].fLignin);
	}
	printf("\n");

	/* the fraction of woody litterfall */
	fscanf(sfp_input, "%s", temp_str);
	for(ic=0; ic<NCLASS; ic++)
	{
		fscanf(sfp_input, "%lf", &bplut[ic].fwoody);
		printf("%lf, ", bplut[ic].fwoody); 
	}
	printf("\n");
	
	fclose(sfp_input);
}

int ReadCFluxData(char *filepath, char *var_name, char *site_name, int year, int var_index, float32 *data_arry)
{
	char filename[255];
	char temp_char[1000];
	char year_char[5];
	char *pch;

	int i, id, flag;
	int ilayer;

	FILE *sfp_input;

	filename[0] = 0;
	strcat(filename, filepath);
	strcat(filename, var_name);
	strcat(filename, site_name);
	/*strcat(filename, ".");*/
	sprintf(year_char, "%d",  year);
	strcat(filename, year_char);
	strcat(filename, ".txt");
	if(!(sfp_input = fopen(filename, "r")))
	{
		printf("cannot open file %s for read!\n", filename);
		exit(-1);
	}

	fgets(temp_char, 1000, sfp_input);
	id = 0;
	flag = 1;
	while(fgets(temp_char, 1000, sfp_input) != NULL)
	{
		pch = strtok(temp_char, " ");

		for(i=0; i<var_index-1; i++)
		{
			pch = strtok(NULL, " ");
		}

		data_arry[id] = atof(pch);
		/* printf("id %d, data: %f\n", id+1, data_arry[id]); */ 
		if(data_arry[id] < 0.)
			flag = 0;

		id++;
	}
	fclose(sfp_input);
	
	if(id != NCOMP)
	{
		printf("id: %d\n", id);
		exit(-1);
	}

	return flag;
}

int ReadTsoilData(char *filepath, char *var_name, char *site_name, int year, float32 *data_arry)
{
	char filename[255];
	char temp_char[5000];
	char year_char[5];
	char *pch;

	int  i, id, ilayer, flag;

	FILE *sfp_input;

	filename[0] = 0;
	strcat(filename, filepath);
	strcat(filename, var_name);
	strcat(filename, site_name);
	strcat(filename, ".");
	sprintf(year_char, "%d", year);
	strcat(filename, year_char);
	strcat(filename, ".txt");
	if(!(sfp_input = fopen(filename, "r")))
	{
		printf("cannot open file %s for read!\n", filename);
		exit(-1);
	}

	fgets(temp_char, 5000, sfp_input);
	id=0; 
	while(fgets(temp_char, 5000, sfp_input) != NULL)
	{
		pch = strtok(temp_char, " ");

		for(i=0; i<3; i++)
		{
			pch = strtok(NULL, " ");
		}

		for(ilayer=0; ilayer<nlayer; ilayer++)
		{
			pch = strtok(NULL, " ");
			data_arry[ilayer*NCOMP+id] = atof(pch);
		}

		id++;
	}
	if(id != NCOMP)
	{
		printf("id: %d\n", id);
		exit(-1);
	}
	fclose(sfp_input);

	/* the deepest layer has no data */
	flag = 1;
	for(id=0; id<NCOMP; id++)
	{
		/*data_arry[(nlayer-1)*NCOMP+id] = data_arry[(nlayer-2)*NCOMP+id]; */

		for(ilayer=0; ilayer<nlayer; ilayer++)
		{
			if(data_arry[ilayer*NCOMP+id] <= Nodata)
				flag = 0;
		}
	}

	return flag;
}

void WriteSOCData(char *filepath, char *site_name, char *var_name, int year, int max_iter, float32 *data_arry)
{
	char filename[255];
	char year_char[5];

	int iter, flag, ipool;
	int ilayer;

	float soc;

	FILE *sfp_output;

	filename[0] = 0;
	strcat(filename, filepath);
	strcat(filename, site_name);
	strcat(filename, "//");
	strcat(filename, site_name);
	strcat(filename, ".");
	strcat(filename, var_name);
	sprintf(year_char, "%d", year);
	strcat(filename, year_char);
	strcat(filename, ".century_SOM_vertical.non-uniform.test.csv");
	if(!(sfp_output = fopen(filename, "w")))
	{
		printf("cannot open file %s for write!\n", filename);
		exit(-1);
	}

	fprintf(sfp_output, "spinup_year, ");
	for(ipool=0; ipool<npool; ipool++)
	for(ilayer=0; ilayer<nlayer; ilayer++)
	{
		fprintf(sfp_output, "pool%d_layer%d,", ipool+1, ilayer+1);
	}
	fprintf(sfp_output, "\n");

	for(iter=0; iter<max_iter; iter++)
	{
		flag = 1;
		fprintf(sfp_output, "%d, ", iter+1);
		for(ipool=0; ipool<npool; ipool++)
		for(ilayer=0; ilayer<nlayer; ilayer++)
		{
			soc = data_arry[iter*npool*nlayer+ipool*nlayer+ilayer];
			fprintf(sfp_output, "%f, ", soc);
			if(soc < 0.)
				flag = 0;
		}
		fprintf(sfp_output, "\n");
		if(flag == 0)
			break;
	}
	fclose(sfp_output);
}


void WriteCFluxData(char *filepath, char *site_name, char *var_name, int year, float32 *NEE_arry, float32 *NPP_arry, float32 *Rh_arry)
{
	char filename[255];
	char year_char[5];

	int id, ipool, ilayer;

	FILE *sfp_output;

	filename[0] = 0;
	strcat(filename, filepath);
	strcat(filename, site_name);
	strcat(filename, "//");
	strcat(filename, site_name);
	strcat(filename, ".");
	strcat(filename, var_name);
	sprintf(year_char, "%d", year);
	strcat(filename, year_char);
	strcat(filename, ".century_SOM_vertical.non-uniform.test.csv");
	if(!(sfp_output = fopen(filename, "w")))
	{
		printf("cannot open file %s for write!\n", filename);
		exit(-1);
	}

	fprintf(sfp_output, "id, NPP, NEE, Rh, ");
	for(ilayer=0; ilayer<nlayer; ilayer++)
	{
		fprintf(sfp_output, "Rh_layer%d,", ilayer+1);
	}
	fprintf(sfp_output, "\n");

	for(id=0; id<NCOMP; id++)
	{
		fprintf(sfp_output, "%d, %f, %f, %f,", id+1, NPP_arry[id], NEE_arry[id], NPP_arry[id]+NEE_arry[id]);
		for(ilayer=0; ilayer<nlayer; ilayer++)
		{
			fprintf(sfp_output, "%f, ", Rh_arry[id*nlayer+ilayer]);
		}
		fprintf(sfp_output, "\n");
	}
	fclose(sfp_output);
}


int ReadTsoilData_v2(char *filepath, char *var_name, char *site_name, int year, float32 *data_arry)
{
	char filename[255];
	char temp_char[5000];
	char year_char[5];
	char *pch;

	int  i, id, ilayer, flag;

	FILE *sfp_input;

	filename[0] = 0;
	strcat(filename, filepath);
	strcat(filename, site_name);
	strcat(filename, var_name);
	sprintf(year_char, "%d", year);
	strcat(filename, year_char);
	strcat(filename, ".csv");
	if(!(sfp_input = fopen(filename, "r")))
	{
		printf("cannot open file %s for read!\n", filename);
		exit(-1);
	}

	fgets(temp_char, 5000, sfp_input);
	id=0; 
	flag = 1;
	while(fgets(temp_char, 5000, sfp_input) != NULL)
	{
		pch = strtok(temp_char, ",");

		for(i=0; i<3; i++)
		{
			pch = strtok(NULL, ",");
		}

		for(ilayer=0; ilayer<nlayer; ilayer++)
		{
			pch = strtok(NULL, ",");
			data_arry[ilayer*NCOMP+id] = atof(pch);
			if(data_arry[ilayer*NCOMP+id] <= Nodata)
			{
				flag = 0;
			}
		}

		id++;
	}
	if(id != NCOMP)
	{
		printf("id: %d\n", id);
		exit(-1);
	}
	fclose(sfp_input);

	return flag;
}











