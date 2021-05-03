/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%								       %%
%% Program : history_v15_060Mo.c      				       %%
%% REFUGEE ( Rotational historiEs oF yoUnG stEllar objEcts )           %%
%% -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-..-.-.-..-.-.-..-.-.-.-%%
%% Rotation evolution for a young M*=0.60Mo star                       %%       
%% -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-..-.-.-..-.-.-..-.-.-.-%%
%% APR 2021  						               %%
%% gapinzone@unal.edu.co                                               %%
%%								       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <ctype.h>
#include <unistd.h>
#define PRECISION 0.00001  /* in Sean's Rt/Rco Solver */
//#define TEFF 4000.0        /* not used, tg used instead */
//#define TACC 3e6   /* escala de tiempo de acrecion */
#define macclimit 2e-14  /* lower limit to accretion rate in Msun/yr */
//#define tg 5.17e15 /*  (44 Myr), Teff = 4277K, to "fit" Seiss model radii */
//#define tg 3.15576e13
//#define tg 1.4e15 /*  (44 Myr), Teff = 4277K, to "fit" Seiss model radii */
#define WINDSTART 44.0 /*Myr*/ 
#define omega_sat 8.0  /* (8.0 Barnes & Sofia 2001) */
#define solarperiod 24.47/*Mean Rotation Period in days Mamajek */ 
#define GAMMAC 1.0
#define BETA 0.01 
#define RCSTART 44.0
#define ksm 2.11  /* MP05 constant for stellar wind rA equation */
#define mm 0.223  /* MP05 exponent for stellar wind rA equation */
#define TINIT 0.03  /* Myr Beginning age of star   */

#define cl 29979245800.0  /*ligh speed in cm / s  */
#define PI 3.14159265
#define SIGMA 5.67e-5 /* cgs    */
#define G 6.67259e-8 /*  cgs   */
#define NSA 31557600e0
#define NSMA 3.15576e13
#define MSUN 1.98855e33 /* g    */
#define RSUN 6.95508e10 /* cm    */
#define omegasun 2.86e-06

#define NPMODEL 300 /* Number of lines Siess Model files for 1Mo*/

#define STRSIZ 100
#define DV 2 /*step size for derivative in Siess Model*/


/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zona de declaracion de subprogramas  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/
double dmui_dtau(double taoi,double mucero, double TACC,double tg);
double dalpha_dtau(double muiprima, double dzdt, double alphai,double zeti,double mui,double torque, double kdos,double tg, double dinertiacoredt, double inertiaenv);
double dalpha_dtaucore(double torqueoncore, double taoc,double inertiacore, double radc, double muc, double alphai, double dinertiacoredt, double alphaicore, double kdosrad, double mc, double dzdtcore, double tg, double zeti);
double corotation(double alphai, double mui);
double BE(double alphai,double zeti,double BCERO,double mamon);
double laqueesconciso(double first_guess,double psi,double ef);
double verladiferencia(double x, double RHS);
double lapendiente(double x);
float polintrad(float exis);
float polintradDEV(float exis);
float polintmascoreDEV(float exis);
float rcore(float exis);
float kadosconv(float exis);
float kadosrad(float exis);
float mrad(float exis);
float didtcore(float exis);
float polintradcoreDEV(float exis);
float litio7(float exis);
float tece(float exis);
float turnovernata(float exis);
float labolo(float jk);
float calculin(float MACC, float mui,float TFINAL,float FINCTTS,float taoc, float BCERO, float TACC, float tg, float zeti, float TEFF, float FRAC);

main ()
{

float endtt,tfinalsim,tdeco,masainidisk,MACC,MACC1,BCERO,FRAC,mui,FINCTTS,TFINAL,modelspeed;
float taoc,TACC,TEFF,tg,zeti0,zeti;
int jj=1;
FILE  *readt;  
readt=fopen("input.dat","r");

//  3 < endtt < 30 Myr 
//  1 < tdeco  < 10           [i.e 10-100Myr]
//  0.01 < masainidisk < 0.1Mo
for(jj=1;jj<=1;jj++)

{

fscanf(readt,"%f\t""%f\t""%f\t""%f\n",&endtt,&tfinalsim,&tdeco,&masainidisk);

  FINCTTS=endtt;

  
  TFINAL=tfinalsim;
  
 // TFINAL=5.0;
  
  
 // TACC=3e6;
  taoc=1e7; /*taoc in yr*/
  BCERO=masainidisk;

         TACC=tdeco;  /* lee el archivo input y asigna TACC en y. (tercera columna input.dat)*/

}

  FRAC=0.6;
 
   mui=FRAC; 
 

   MACC= ( 2.51342 * log10(FRAC) ) - 8.179010;
   
   
   

   
//Radio inicial de la estrella
  zeti =5.66 ;  
 //tg =  58*3.15576e13;
  TEFF=4011.0;
  
  
//tg=3.0*G*mui*mui*MSUN*MSUN / (28.0*PI*SIGMA*zeti*zeti*zeti*RSUN*RSUN*RSUN*pow(TEFF,4.0)) ;
 
 //tg=(3e7*mui*MSUN/(4.0*PI*SIGMA*pow(TEFF,4.0)*pow(zeti*RSUN,2.0)*zeti/3.844e33))*NSA;

  //if (FINCTTS < 0.8)
  //{
  
    //tg=(80e6)*NSA;
  //}
  
  //tg=(60e6)*NSA;   /* for star with disk */ 
  
  
  tg=(147e6)*NSA;

  
//tg=3.0*G*0.6*0.6*MSUN*MSUN / (28.0*PI*SIGMA*pow(0.512, 3.0)*RSUN*RSUN*RSUN*pow(TEFF,4.0)) ;



//tg=NSMA;
 modelspeed=calculin(MACC,mui,TFINAL,FINCTTS,taoc,BCERO,TACC,tg,zeti,TEFF, FRAC);
 
 printf("\n");
 printf("Tiempo(Myr)=%6.4f\n""Velocidad(km/s)=%6.4f\n""Radio(Ro)=%6.4f\n""Masa(Mo)=%6.4f\n""Teff(K)=%6.4f\n",TFINAL,modelspeed,zeti0,mui,TEFF);
 printf("\n");
 
 return 0;
 
}
  
  
  
 
/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zona definicion de subprogramas   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/  

  
  
float calculin (float MACC, float mui,float TFINAL,float FINCTTS,float taoc, float BCERO, float TACC, float tg, float zeti, float TEFF, float FRAC)
 
{
  
/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zona declaracion de variables  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/
 double taoi,efecin,taofinal,mucero,turnover,mamon;
  double m1,m2,m3,m4,k1,k2,k3,k4,l1,l2,l3,l4;
  double alphai,alphaicore,muiprima,rat1,torque,torque1,wm;
  double omegm,acc_core,kmn;
  double queca,kbd,B,Bconv,x,funcioncita,zero;
  double Tc,rhomean,rhoc,pc,tg1,k_grande,dinertiacoredt;
  double timeconvectivo,alphainit,deltaj;
  double gamin,factor,muc;
  double  age,iconv,dostercios,age1;
  double xc,omrt,cte,cte1,rart;
  double rt,rin,rout,torqueacc_mp,torquemag_mp;
  double ef,efcore,dip,psi,omrin;
  double ere,vamola,torquedec_gp,kdos,kdosconv,kdosrad;
  double coisa,veqstar,radmodel,imodel;
  double radialdisci,omegadisc,period,vk,vd,gamma;
  double torque_wind,alphacore, ksk,ksk_new,temp_eff,torqueneto[10000],torque_winds;
  double espesor,beto,mwind,aefe,nicoret,nicoret1,kei;
  double omega_crit,torque_wind_sku,torque_wind_mayor,torque_wind_sat,torqueoncore;
  double kmm,drdt,ceqprepare,feq,omeq,BTT;
  double INICIAL;
  double ceq,velo,velocore,dzdtcore,xeq,xequil;xeq,xequil;
  double tover,tdeco;
  float inertiacore,exis,lithiumseven;
  double rossby,inertiaenv;
  double mc,radc=0.0,imc,ce;  //core mass and other core properties
  double muiprimas, muiprimass; // placeholders for RK
  double dzdt,equisss,turno,rhk;            // dR/dt (normalized units)
  double tcentral,peq; //Tc in K Siess
  double zc,oi;
  double swfrac; // fraction of accretion rate that goes to stellar wind
  double sgfrac; // 1-swfrac, net fraction gained by star
  double mdot_sw; // stellar wind mass loss rate 
  double torque_sw,torqueoncore1,torqueoncore2,cujn,cumn,enne; // stellar wind torque
  double vesc;     // escape speed from stellar surface
  double rarstar,gilmour;  // stellar wind alfven radius/stellar radius
  double rarstarms; 
  double feqsw,veloeq;    // equilibrium f for APSW
  double omeqsw,eletot,eleacc,ele;   // equilibrium spin rate for APSW
  double step, stepm, stepz, stepa; // timestep variables
  double C_M = 0.01;  // timestep factor, sets max change in Mstar
  double C_Z = 0.001;  // timestep factor, sets max change in Rstar.
  double C_A = 0.01;  // timestep factor, sets max change in Omegastar.
  double omegaeq,a0,a1,qu,erre,de,ese,tee,eleque,eleque1,eleque2,eleque3;
  float periodo,jk,jk1,mk,mv,dist,propermotion,pm1,pm2,pm3,prlx,prob,agess,bvdelorme;
  float bmenosv,bece,flux;
  double a,b,c,d,term1,runotres;
  char   datfile[STRSIZ],datfiledisk[STRSIZ],ids[50];
  char   id [50];
  int    s=1,eni=10;
  int    i=1,ii=1,jota=1;          // loop variable
  int    state;        // magnetic conneB500GINFM1e-2A03BT1.datction state of MP05
  FILE   *write,*writeaqui,*escriba,*lea,*writeunin,*readlabolo;
   
  
/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!! zona lectura de archivo de entrada   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/  
 
//lee archivo que contiene B-V (12 datos representativos entre 0.36 y  1.48
//este archivo se usa pues permite conocer la respectiva BC usando el
//subprograma labolo()
  readlabolo = fopen("labolo.txt","r"); 

  
     
/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zona condiciones iniciales   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/  
 INICIAL=0.3;
 BTT=2000.0;
 
//fraccion de la tasa de acrecion que se va en el viento estelar
  swfrac = 0.1;  
 
//fraccion neta de acrecion que se va a la estrella  
  sgfrac = 1.0 - swfrac; 

//tao final de la simualcion; (tao es el tiempo normalizado por tg)
//taofinal = TFINAL*NSMA/tg;
taofinal=TFINAL*NSMA/tg;
//tao inicial de la simulacion;
efecin   = TINIT*NSMA/tg;
taoi     = efecin;


//dmu/dtao incial i.e. en el tao inicial
//mucero=(MACC/(TACC*(exp(-efecin/(TACC*NSA/tg)) - exp(-(44*NSMA/tg)/(TACC*NSA/tg))) ))*(tg/NSA); /*Msun/yr*/

 mucero = ( tg /  MSUN ) * pow(10, MACC) * ( MSUN / NSA ) * exp(TFINAL*1e6/TACC);
//mucero = MACC / TACC ;
//mucero =0.0000001*pow(mui,2)*(tg/NSA);
//mucero =1e1*(exp((-8) + ((3e6)/TACC)));
//-taoi*(tg/NSA)/TACC
// Frecuencia angular inicial de rotacion (alphai):
alphai=INICIAL*sqrt(G*FRAC*MSUN/pow(RSUN*zeti,3.0))/omegasun;

//BCERO=pow(6.28/(alphai*omegasun*86400),-1.32);
  
 /*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zona de configuracion de usuario   (desactivada)   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  printf("\n");printf("\n");
  printf(":: ROTATIONAL EVOLUTION OF A SOLAR-TYPE STAR :: \n");
  printf("-.-.-.-.-.-.-.-.-.-.-.-.-.-..-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n");
  printf("\n");printf("\n");
  printf("Please, set the model to your interest: \n");
  printf("\n");printf("\n");
  printf("Initial Stellar Mass (in Solar Units):""\n");
  scanf("%lf",&FRAC);
  printf("Initial Disk Mass (in Solar Units):""\n");
  scanf("%lf",&MACC);
  mui=FRAC;  BCero// Initial mass
  //MACC = 1.0 - FRAC;
  printf("Duration of CTTS stage (in Myr):""\n");
  scanf("%lf",&FINCTTS);
  printf("Simulation ends at (in Myr) [4900 for run up to solar age]:""\n");
  scanf("%lf",&TFINAL);
  printf("Initial Magnetic Field Strenght (in G) :""\n");
  scanf("%lf",&BCERO);
  printf("Initial Angular Velocity (in Solar Units):""\n");
  scanf("%lf",&INICIAL);
  printf("Decoupling time-scale (in yr):""\n");
  scanf("%lf",&taoc);
 */

  
/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!! zona nombre del archivo de salida   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/  
 sprintf(datfile, "output");
 
  write=fopen(datfile,"w");

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zona comienzo del loop principal    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/  


// taoi es adimensional  taoi = Millones de años x NSMA / tg  con tg en segundos
  while (taoi <= taofinal) {
tcentral=tece((taoi*(tg/NSMA))*1e6);
//  Tasa de acrecion de disco
    muiprima=dmui_dtau(taoi,mucero,TACC,tg);  
    if (muiprima*NSA/tg <= macclimit )  /* Si multiplico muiprima*NSA/tg entonces tengo muiprima en Mo/y */
    {
      muiprima=macclimit*tg/NSA;
    }
 
    /* Si multiplico muiprima*NSA/tg entonces tengo muiprima en Mo/y */
   
    
//  Radio estelar 
 zeti=polintrad((taoi*(tg/NSMA))*1e6);

 ele=SIGMA*pow(TEFF,4.0);
  eleacc=G*mui*MSUN*(muiprima/(tg/NSA))*(MSUN/NSA)/(4.0*PI*zeti*RSUN*zeti*RSUN*zeti*RSUN);
  eletot=ele+eleacc;
  
  
// derivada del radio respecto al tiempo
  dzdt=polintradDEV((taoi*(tg/NSMA))*1e6); /* in Ro/yr*/

//  Stellar k^2 Siess et al.2000 
    kdosconv=kadosconv((taoi*(tg/NSMA))*1e6);
    kdosrad=kadosrad((taoi*(tg/NSMA))*1e6);
    kdos=kdosconv;

//  Mass of the radiative core in Mo and its derivative in Mo/yr   
    muc=mrad((taoi*(tg/NSMA))*1e6);
    mc=-polintmascoreDEV((taoi*(tg/NSMA))*1e6)/(tg/NSA);
    lithiumseven=litio7((taoi*(tg/NSMA))*1e6);
    
// taoi es adimensional  taoi = Millones de años x NSMA / tg  con tg en segundos

    
//  Radius of the radiative core in Restelares
    radc=rcore((taoi*(tg/NSMA))*1e6)*zeti;

//  turnover convectivo de Landin et al.2010 A&A 510,A46 table 1
//  edad, luminosidad, teff, logg, logtc, tg, Ro para 1Mo PMS models con rotacion
 //turnover=turnovernata((taoi*(tg/NSMA))*1e6);

  turnover=41.3960*pow(taoi*(tg/NSMA),2.0)-419.985*taoi*(tg/NSMA)+1338.7278;
     
//  Moment of inertia of the star Siess et al.2000 normalized to MSUN x RSUN**2    
    inertiacore=kdosrad*muc*pow(radc,2.0);
    inertiaenv=kdos*mui*MSUN*pow(zeti*RSUN,2.0);
//     
//  Derivative of the moment of inertia of the core normalized and per year
    dinertiacoredt=(didtcore((taoi*(tg/NSMA))*1e6));
   //BCERO=8.94-4.849*log(taoi*(tg*1e6/NSMA))+0.624*log(taoi*(tg*1e6/NSMA))*log(taoi*(tg*1e6/NSMA));
//BCERO=1000.0;
//      Stellar Dipolar Magnetic field. Saturation limit = 3000G
   
   B=BE(alphai,zeti,BCERO,turnover); 
    
   // taoi es adimensional  taoi = Millones de años x NSMA / tg  con tg en segundos

    
//printf("OIIII=%6.3e\n",oi);

   
enne=2.0;
    cumn=pow((20.0/B),(enne-1.0)/20.0)*pow(2.71,(1.22-1.42*enne+0.19*enne*enne+0.01*enne*enne*enne));
    
    cujn=4.05*pow(2.71,-1.4*enne)+(enne-1)/(60.0*B*enne);
    
// f
    ef=(alphai*omegasun)/sqrt(G*mui*MSUN/pow(zeti*RSUN,3.0));

    
// Enforce maximum spin rate.
if (ef > 1.0) {
    ef = 1.0;
 	alphai = ef *sqrt(G*mui*MSUN/pow(zeti*RSUN,3.0))/omegasun;
    }

    
// Corotation Radius
  queca=corotation(alphai,mui);

    
// Dipolar moment of the star
//dip=B*0.5*pow(zeti*RSUN,3.0);
  dip=B*pow(zeti*RSUN,3.0);

  /* Si multiplico muiprima*NSA/tg entonces tengo muiprima en Mo/y */
  
/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zona parametros de interaccion estrella disco   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/  
// If star is accreting, calculate star-disk interaction variables.

    if(muiprima > 0.0) {

      // Psi = (internal disk torque)/(magnetic torque per gamma)
      psi=2.0*pow(dip,2.0)*pow((muiprima/(tg/NSA))*(MSUN),-1.0) * 
	pow(G*mui*MSUN,-0.5)*pow(zeti*RSUN,-3.5);

      // Computation of Rout 
      rout=pow(1.+(GAMMAC*BETA),(2./3.))*queca;

      // Computation of Rin 
      rin=pow(1.-(GAMMAC*BETA),(2./3.))*queca;
 // printf("Minit(Msun)=%6.2f\n""Omega_init(SolarUnits)=%6.4f\n""tacc(yr)=%6.2e\n",mui,alphai,TACC);
 // printf("tKH(s)=%6.4e\n""Macc(Mo)=%6.2f\t""R*=%f\n",tg,MACC,zeti);
 // printf("Teff(K)=%6.2f\n""tINIT(Myr)=%6.4f\n""tFINAL(Myr)=%6.2f\n",TEFF,efecin*(tg/NSMA),taofinal*(tg/NSMA));
  // printf("Rrad(Rsun)=%6.2f\n",radc);

      // Determine the MP05 state of the system.
      if (ef < (1-BETA*GAMMAC)*pow(GAMMAC*psi,-3./7.))
	state = 1; else state = 2;

      /***** Calculate Rt and torques. *****/
      if (state == 1) {

	rt = pow(GAMMAC*psi, 2./7.)*zeti;  // Rt/Rsun
	
	// Check to see if disk reaches stellar surface.
	if (rt < zeti) rt = zeti;

	zero = rt/queca;  // Rt/Rco

	// Magnetic torque
	torquemag_mp = 0.0;  // no magnetic torque in state 1

      } else {

	zero=laqueesconciso(zero, psi, ef);

	// Check to see if disk reaches stellar surface.
	if (zero < zeti/queca) zero = zeti/queca;

	rt=zero*queca;

        
/*
%%%%%%%%%%%%%%%%%%%%%%%%%%
zona torque magnetico   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
*/  
	torquemag_mp=(pow(dip,2.0)/(3.0*BETA*pow(queca*RSUN,3.0))) * 
	  (2.0*pow(queca/rout,1.5)-pow(queca/rout,3.0) - 
	   2.0*pow(queca/rt,1.5)+pow(queca/rt,3.0));
      }
      
      /* end if state 1, 2 */
if ( taoi*(tg/NSMA) > FINCTTS   )
{
    
    rt=zeti;
 
    
    
}

/* Si multiplico muiprima*NSA/tg entonces tengo muiprima en Mo/y */
    
/*
%%%%%%%%%%%%%%%%%%%%%%%%
zona torque acrecion  %%
%%%%%%%%%%%%%%%%%%%%%%%%
*/  
      torqueacc_mp=(muiprima*(NSA/tg)*(MSUN/NSA))*pow(G*mui*MSUN*rt*RSUN,0.5); 
    } else {  // If not accreting, set star-disk variables to zeros.

      rout = 0.0;
      rin  = 0.0;
      rt   = 0.0;
      zero = 0.0;
      torquemag_mp = 0.0;
      torqueacc_mp = 0.0;
    }

    /* Si multiplico muiprima*NSA/tg entonces tengo muiprima en Mo/y */
/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zona torque de viento estimulado por acrecion  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/  
    // Stellar wind mass loss rate
    mdot_sw = swfrac * muiprima*(NSA/tg)*MSUN/NSA;

    // Surface escape speed
    vesc = sqrt(2.0*G*mui*MSUN/(zeti*RSUN));

    // Calculate Alfven radius in wind
    if(mdot_sw != 0.0) {
      rarstar = ksm * pow(B*B*zeti*zeti*RSUN*RSUN/(mdot_sw*vesc),mm);
    } else{
      rarstar = 0.0;
    }

    // Stellar wind torque
    gilmour = 0.09 * rt/zeti;
    torque_sw = -mdot_sw*alphai*omegasun*zeti*zeti*RSUN*RSUN*rarstar*rarstar;

    //rastar es el radio de alfven de la estrella em Ro
    
    
    
   
   
    
    
/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zona torque de viento estimulado por actividad %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/


   ksk=9e48; 
 if(alphai*pow(omega_sat,2.0)<=pow(alphai,3.0))
 {
torque_winds = -ksk*pow(zeti/mui,0.5)*omegasun*alphai*pow(omegasun*omega_sat,2.0);
}
else
{
torque_winds = -ksk*pow(zeti/mui,0.5)*pow(alphai*omegasun,3.0);

} 


    
    //if(alphai*pow(omega_sat,2.0)<=pow(alphai,3.0))
  // {
  //torque_winds = -9.5e30*pow(turnover/12.9,2.0)*pow(alphai,3.0)*cujn;
    
//}
 //else
 //{
  //torque_winds = -9.5e30*pow(zeti,3.1)*pow(mui,0.5)*pow(10.0,2.0)*alphai*cujn;
  //} 
   
    
/*
%%%%%%%%%%%%%%%%%%%%%%
zona torque  total  %%
%%%%%%%%%%%%%%%%%%%%%%

// taoi es adimensional  taoi = Millones de años x NSMA / tg  con tg en segundos

*/  
if((taoi*(tg/NSMA)) <= FINCTTS  )
   {
    alphaicore=alphai;
    efcore=ef;
    deltaj=0.0;
    torqueoncore=0.0; 
    torque_winds=0.0;

//torque durante la etapa T Tauri 
 torque = torqueacc_mp + torquemag_mp + torque_sw;  
   
//  torque = 0.0;
    torque1=torque/1e37;
  
       
   
}

  
  
if((taoi*(tg/NSMA)) > FINCTTS && (taoi*(tg/NSMA)) <= WINDSTART )
  
//  if((taoi*(tg/NSMA)) > FINCTTS && (zeti-radc)/zeti >= 0.3 )
  {

   alphaicore=alphai;


 deltaj=0.0;
 efcore=ef;


 
 //AQUICITO DESCOMENTE NOMAS
 torqueoncore=0.0;
torque_winds=0.0;

//torque = torqueacc_mp;
torque = 0.0;

} 

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zona torque transferencia interna de momento angular %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/  




if( (taoi*(tg/NSMA)) > WINDSTART)
 // if ( (zeti-radc)/zeti < 0.3 )  
{

dzdtcore=polintradcoreDEV((taoi*(tg/NSMA))*1e6);

deltaj=(inertiaenv*inertiacore*MSUN*RSUN*RSUN*omegasun)*(alphaicore-alphai)/(inertiacore+inertiaenv);

efcore=(alphaicore*omegasun)/sqrt(G*muc*MSUN/pow(radc*RSUN,3.0));

torqueoncore1=deltaj/(taoc*NSA);
torqueoncore2=-(2.0/3.0)*pow(radc*RSUN,2.0)*alphai*omegasun*mc;

torqueoncore=(torqueoncore1+torqueoncore2);

 // torque =  torqueacc_mp;
torque = 0.0 ;
  
}
  
/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zona Integracion de las ecuaciones    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/  

    /***************************************/
    /***** Calculate optimal timestep. *****/
    /***************************************/


    // Calculate timestep for changes in Mass.
    stepm = (muiprima == 0.0) ? 1.0e20 :
      fabs(C_M * mui / (muiprima*sgfrac));
//     Calculate timestep for changes in radius.
     stepz = (dzdt == 0.0) ? 1.0e20 :
      fabs(C_Z * zeti / dzdt);

    // Calculate timestep for changes in rotation rate.
    stepa = (torque == 0.0) ? 1.0e20 :
      fabs(C_A * alphai * omegasun * (kdos)*mui*MSUN*pow(RSUN*zeti,2.0) / 
	   torque)/tg;

    // Set best timestep as minimum of above steps.
    step = (stepm <= stepz) ? stepm : stepz;
    step = (step  <= stepa) ? step  : stepa;

    // Be sure timestep does not go beyond final time.
    if (taoi+step >= taofinal) {
      step = taofinal-taoi;
      // Do one last calculation for output of final variables.
      if (step == 0.0) step = 1.0e-8;
    }
    //    printf("timesteps %e %e %e %e\n",stepm, stepz, stepa, step);

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zona parametros de salida    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/  
    /***** Calculate quantities used only for output. *****/
    
    // Spin period of the envelope
    period=2.0*3.14159/(alphai*omegasun*86400e0);

    velo=(alphai*omegasun*(zeti)*RSUN)/(100000.0);

   // Speed of the core
    
    if((taoi*(tg/NSMA)) <= RCSTART)
    {
      velocore=velo;
    }
    else
    {
    velocore=(alphaicore*omegasun*radc*RSUN)/100000.0;
    }
    
    
    
   // Equilibrium disk locked spin rate (MP05 eq. 27)
    kei=2.*pow(1.+BETA*GAMMAC,-1.)-pow(1.+BETA*GAMMAC,-2.);
    xeq = (8.0 + sqrt(64. - 28.*kei))/14.;
    ceqprepare=2.*pow(BETA,-1.)*pow(xeq,4./3.)*(xeq-1.);
    ceq=pow(ceqprepare,-3.0/7.0);

    //    omeq = pow(ceq*pow(muiprima*MSUN/tg,3./7.)*pow(G*mui*MSUN,5./7.) * 
    //	       pow(BCERO/omegasun,-6./7.)*pow(RSUN,-12./7.) * 
    //	       pow(zeti*RSUN,-6./7.),7./13.)/omegasun;


    //CAUTION!!!! Testing constant mangetic field.

    omeq = ceq*pow(muiprima*MSUN/tg,3./7.)*pow(G*mui*MSUN,5./7.) * 
	       pow(B,-6./7.)*pow(zeti*RSUN,-18./7.) / omegasun;
	       
feq = (omeq*omegasun)/pow(G*mui*MSUN/pow(zeti*RSUN,3.0),0.5);

    // Equilibrium APSW spin rate, assuming Rt=Rco (MP08III eq. 18)

    if(mdot_sw != 0.0) {
      feqsw = pow(ksm,-1.5)/pow(swfrac,.75-1.5*mm)*
	pow(psi*pow(2.0,-1.5),-1.5*mm);
    } else{
      feqsw = 0.0;
    }
    
      omeqsw = feqsw*sqrt(G*mui*MSUN)*pow(zeti*RSUN,-1.5)/omegasun;

 
  omegaeq = alphai;
 
veloeq=200.0;
peq=10000.0;
    

for(jota=1;jota<=12;jota++)
		
	{
		
		fscanf(readlabolo,"%f\n",&bmenosv);
		//  fprintf(write,"%6.3f\t""%6.3f\n",-vh+vk,-vj+vh);
	//	
                
                equisss=1.0-bmenosv;
                if(equisss>0)
                {
                turno=1.362-0.166*pow(equisss,2.0)-5.323*pow(equisss,3.0);
                }
                else
                 {
                turno=1.362-0.14*equisss;
                }   
                
                rhk=log10((6e-5)*exp(-0.9*period/turnover));
                
                bece=labolo(bmenosv);  
                                                rart=(rarstar-rt)/zeti;

	//	printf("%6.3f\n",bece);   
		
	}

    /***** Output to screen. *****/
    
    
    printf("\n""\n""\n""t(Myr)=%6.2f\t""alpha=%6.3f\t""alphacore=%6.3f\n""mdot=%6.3e\t""R*(Rsun)=%6.3f\t""Rco(R*)=%6.3f\n" "Rt(R*)=%6.3f\t""Rout(R*)=%6.3f\t""B=%6.3f\t""T=%6.3e\n""T_sw=%6.3e\t""Ta=%6.3e\t""Twinds=%6.3e\t""dzdtcore=%6.3e\n""root=%6.3f\t""dip=%6.3e\t""stepz=%8.6f\t""stepm=%8.6f\n""stepa=%6.6f\t""f=%6.3f\t"
"kei=%6.3f\t""mui=%6.3f\n""dzdt=%6.4f\t""Rrad(Rsun)=%6.3f\n""step(yr)=%6.3f\t""BCERO=%6.3f\t""tg(Ma)=%6.4f\t""turnover=%6.4f\n""rhk=%6.3f\t""period=%6.3f\n", taoi*(tg/NSMA),alphai,alphaicore,muiprima*NSA/tg,zeti,queca/zeti,rt/zeti,rout/zeti,B,torque,torque_sw,torqueacc_mp,torque_winds,dzdtcore,zero,dip,stepz*(tg/NSA),stepm*(tg/NSA),stepa*(tg/NSA),ef,kei,mui,dzdt,radc,step*(tg/NSA),BCERO,tg/NSMA,turnover,rhk,period);



    /***** Print to file. *****/

   fprintf(write,"%e\t%e\t%e\t",taoi*(tg/NSA), alphai*2*3.14159,muiprima*NSA/tg);
    fprintf(write,"%e\t%e\t%e\t",zeti, velo, period);
    fprintf(write,"%e\t%e\t%e\t",turnover, rt*0.00465047, queca);
    fprintf(write,"%e\t%e\t%e\t",torque_sw*1e-38, torqueacc_mp*1e-38, torquemag_mp*1e-38);
    //ef estaba antes de de
    fprintf(write,"%e\t%e\t%e\t",vesc, inertiacore,torque);
    fprintf(write,"%e\t%e\t%e\t",FRAC, mui, mdot_sw);
    fprintf(write,"%e\t%e\t%e\t",dzdt, dinertiacoredt, mdot_sw*NSA/MSUN);
    fprintf(write,"%e\t""%e\t""%e\t",dip,kdos,torque);
     fprintf(write,"%e\t""%e\n",B,rarstar*0.00465047);
 //    fprintf(write,"%e\t%e\t%e\t",);

 /*
 printf("\n");printf("\n");
 printf("-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.\n");
  printf(" :: Used parameters  :: \n");
  printf("\n");
  printf("tINIT(yr)=%6.2f\n""tFINAL(Myr)=%6.2f\n",efecin*(tg/NSA),taofinal*(tg/NSMA));
  printf("tDISK(Myr)=%6.2f\n",FINCTTS);

  printf("tKH(Myr)=%6.2f\n""tacc(Myr)=%6.1f\n",tg/NSMA,TACC/1e6);
  printf("tao_DEC(yr)=%6.2e\n",taoc);
  
  printf("R*(Ro) [at the end of simulation!]=%6.3f\n",zeti);
     printf("Rrad(Ro)[at the end of simulation!]=%6.2f\n",radc);

  printf("B(G)[at t=0!]=%6.1f\n""B(G)[at the end of simulation!]=%6.1f\n""Teff(K)=%6.2f\n""BETA=%6.4f\n""GAMMAC=%6.1f\n",BCERO,B,TEFF,BETA,GAMMAC);
  printf("M* at t=0 (in Mo)=%6.2f\n""Mass of the Disk at t=0 (in Mo)=%6.2f\n""Omega_init (in Solar Units)=%6.2f\n",FRAC,MACC,INICIAL);
  printf("KSK=%6.2e\n""KMN=%6.2e\n""Omega_Sat (in Solar Units)=%6.2f\n",ksk,kmn,omega_sat);
  printf("\n");printf("\n");
 printf("\n");printf("\n");
 */
 
    /****************************************/
    /***** Do fourth-order Runge-Kutta. *****/
    /****************************************/

    muiprimas=dmui_dtau(taoi+step/2.0,mucero,TACC,tg);
    muiprimass=dmui_dtau(taoi+step,mucero,TACC,tg);

    k1=step*dmui_dtau(taoi,mucero,TACC,tg)*sgfrac;
    l1=step*dalpha_dtau(muiprima,dzdt,alphai,zeti,mui,torque,kdos,tg, dinertiacoredt, inertiaenv);
    m1=step*dalpha_dtaucore(torqueoncore,taoc,inertiacore, radc,muc,alphai,dinertiacoredt,alphaicore,kdosrad,mc,dzdtcore,tg, zeti);
			
		
    k2=step*dmui_dtau(taoi+step/2.0,mucero,TACC,tg)*sgfrac;
    l2=step*dalpha_dtau(muiprimas,dzdt,alphai+l1/2.0,zeti,mui+k1/2.0,torque,kdos,tg,dinertiacoredt, inertiaenv);
    m2=step*dalpha_dtaucore(torqueoncore, taoc,inertiacore, radc, muc, alphai+l1/2.0,dinertiacoredt,alphaicore+m1/2.0,kdosrad,mc,dzdtcore,tg,zeti);			
    
    k3=step*dmui_dtau(taoi+step/2.0,mucero,TACC,tg)*sgfrac;
    l3=step*dalpha_dtau(muiprimas,dzdt,alphai+l2/2.0,zeti,mui+k2/2.0,torque,kdos,tg, dinertiacoredt, inertiaenv);
    m3=step*dalpha_dtaucore(torqueoncore, taoc,inertiacore, radc, muc, alphai+l2/2.0,dinertiacoredt,alphaicore+m2/2.0,kdosrad,mc,dzdtcore,tg,zeti);
    
    k4=step*dmui_dtau(taoi+step,mucero,TACC,tg)*sgfrac;
    l4=step*dalpha_dtau(muiprimass,dzdt,alphai+l3,zeti,mui+k3,torque,kdos,tg, dinertiacoredt, inertiaenv);
    m4=step*dalpha_dtaucore(torqueoncore, taoc,inertiacore, radc, muc, alphai+l3, dinertiacoredt,alphaicore+m3,kdosrad,mc,dzdtcore,tg,zeti);
   
    mui+=(k1/6.0+k2/3.0+k3/3.0+k4/6.0);
    alphai+=(l1/6.0+l2/3.0+l3/3.0+l4/6.0);
    alphaicore+=(m1/6.0+m2/3.0+m3/3.0+m4/6.0);
 
    // Increment time
   taoi += step;

  }  /* end main time loop */
 printf("\n");
 printf(" ## computation was done succesfully ##  \n");
 printf(" ## please check the file output.dat ##  \n");
 printf("\n");

  
 
 
 
 
 
  // Close input file.
 fclose(write);

  return velo;

} /* end main */

//************************************************
//***** Calculates Accretion Rate, dmui/dtau *****
//************************************************
double dmui_dtau(double taoi,double mucero, double TACC,double tg)
{
  double dmu_dtao;
  // Exponential accretion decay, from CC93
  dmu_dtao=mucero*exp(-taoi*(tg/NSA)/TACC);
  // Hyperbolic tangent function from ?? Yi??
    

//  dmu_dtao=mucero*(1.+tanh(1.-5.*log(taoi/0.18)));
  
  //dmu_dtao=mucero*(1.+tanh(1.-5.*log(taoi/0.04)));
  
  
  
  
  // Power-law accretion rate
   //   dmu_dtao=1e-8*tg/NSA*pow(TACC/taoi/(tg/NSA),2.0);
  // Constant accretion rate
  //    dmu_dtao = 1e-8/NSA*tg;
  // If less than lower limit, set to zero.
  //  dmu_dtao = (dmu_dtao > macclimit*tg/NSA) ? dmu_dtao : 0.0;
  return dmu_dtao;
}


//**********************************************
//***** Calculates dalpha/dtau for envelope*****
//**********************************************

double dalpha_dtau(double muiprima, double dzdt, double alphai,double zeti,double mui,double torque,double kdos, double tg, double dinertiacoredt, double inertiaenv)
{
    

 //return  tg*torque/(kdos*MSUN*pow(RSUN,2.0)*omegasun) - (alphai/NSA)*(((muiprima*NSA/tg)/mui)+(2.0*dzdt/zeti));

    return (tg*torque/(mui*MSUN*pow(zeti*RSUN,2.0)*omegasun)) +
   ( (alphai/mui) * (  ( 2.0*zeti*zeti*zeti/mui )   - 5.0*muiprima) );
    
  
}

//**********************************************
//***** Calculates dalpha/dtau for envelope*****
//**********************************************
double dalpha_dtaucore(double torqueoncore, double taoc,double inertiacore, double radc, double muc, double alphai, double dinertiacoredt, double alphaicore,  double kdosrad, double mc, double dzdtcore, double tg, double zeti)
{
return -taoc*NSA*torqueoncore/(kdosrad*muc*MSUN*pow(RSUN*radc,2.0)*omegasun)+
    (alphaicore/NSA)*(mc/muc+(2.0*mc/radc));
    
}

/****************************************/
/***** Calculates corotation radius *****/
/****************************************/
 double corotation(double alphai, double mui)
{
  double radio;
  radio=pow((G*mui*MSUN/(pow(alphai*omegasun,2.0))),1./3.)/RSUN;
  return radio;
}

/*************************************/
/***** Calculates Magnetic Field *****/
/*************************************/
double BE(double alphai,double zeti,double BCERO, double turnover)
{
  double magnetic;
  //magnetic=BCERO*alphai/pow(zeti,2.0);
 // magnetic=BCERO*turnover*(alphai/6.28)/120.0;
 magnetic = BCERO;
   

  return magnetic;
}

/***************************************/
/***** Sean's Rt/Rco (zero) solver *****/
/***************************************/
double laqueesconciso(double first_guess, double psi, double ef)
{
  double valu, new_guess, fval, slope, falta;
  valu = BETA / psi / pow(ef,7./3.);
  new_guess = first_guess;
  fval = verladiferencia(new_guess, valu);
  falta = fval / valu;
  while((falta > PRECISION) || (falta < -PRECISION)) {
    slope = lapendiente(new_guess);
    new_guess = new_guess - fval/slope;

    // Check new_guess for negatives.
    if(new_guess <= 0.0) new_guess = PRECISION;
    fval = verladiferencia(new_guess, valu);
    falta = fval / valu;
  }
 return new_guess;
}

/****************************************************************/
/***** For Rt solver: calculates LHS - RHS of MP05 eq. (15) *****/
/****************************************************************/

double verladiferencia(double x, double RHS)
{
  return (pow(x, -3.5) - 1.0/(x*x) - RHS);
}

/*******************************************************************/
/***** For Rt solver: calculates slope of LHS of MP05 eq. (15) *****/
/*******************************************************************/

double lapendiente(double x)
{
  return (2.0/(x*x*x) - 3.5 * pow(x,-4.5));
}

/****************************************************************************/
/***** For Radius from Siess(2000) : calculates the radius at any time *****/
/**************************************************************************/


float tece(float exis)

{

 float xa[1000],ya[1000],c[1000],d[1000],den,dif,dift,ho,hp,w;
 float y,dy,age1,radius,age2,k2cmk2r,xx,xxin,yy,yyin,cados;
 int j=1,i,m,ns=1;
 float model,phase,lum,mbol,Reff,teff,rhoeff,logg, masa,idade; 
 float tc, roc, etac,mcore, rcore, lnuc, lgrav, k2conv, k2rad,irad,iconv;
 FILE *readmdtres,*readthisdtres,*write3tres,radius1;
 
 // readm = fopen("k2ov1.csv","r"); /* 209 puntos */
  readmdtres = fopen("m0.6z02.hrd","r"); 
  //ojo:
 readthisdtres = fopen("edad.dat","r"); 
 write3tres=fopen("mcore.txt","w"); 
 while(j<=NPMODEL)

{
fscanf(readmdtres,"%d\t""%d\t""%f\t""%f\t""%f\t""%f\t""%f\t""%e\t""%e\t""%f\t""%e\n",&model, &phase, &lum, &mbol,&Reff, &radius1, &teff, &rhoeff, &logg,&masa,&age1);
  
 fscanf(readthisdtres,"%f\n",&idade);
   xx=idade;
   yy=teff;
   xa[j]=xx;
   ya[j]=yy;

  j++;
}

 fclose(readmdtres);
 fclose(readthisdtres);
 fclose(write3tres);
dif=fabs(exis-xa[1]);
 for(i=1;i<NPMODEL;i++){
   if((dift=fabs(exis-xa[i]))<dif) {
     ns=i;
     dif=dift;
   }
   c[i]=ya[i];
   d[i]=ya[i];
}
 y=ya[ns--];
 return y;

}





float polintrad(float exis)

{

 float xa[5000],ya[5000],c[5000],d[5000],den,dif,dift,ho,hp,w;
 float viejin[5000]; 
float y,dy,age1,radius1,age2,k2cmk2r,xxx,xxin,yyy,yyin;
 int s=1,i,m,ns=1,model,phase;
 float lum,mbol,Reff,teff,rhoeff,logg, masa; 

 FILE *read11;
 read11 = fopen("m0.6z02.hrd","r");
 //read11 = fopen("vallejo.dat","r");


while(s<=NPMODEL)
//while (s<=10)
{


fscanf(read11,"%d\t""%d\t""%f\t""%f\t""%f\t""%f\t""%f\t""%e\t""%e\t""%f\t""%e\n",&model, &phase, &lum, &mbol,&Reff, &radius1, &teff, &rhoeff, &logg,&masa,&age1);

//fscanf(read11,"%f\t""%f\n",&age1,&radius1);


   xxx=age1;
   yyy=radius1;

//printf("%6.2f\n", xxx);
// printf("\n");printf("\n");
//system("pause"); 

//   viejin[j]=xx;


// printf("\n");printf("\n");
 

   xa[s]=xxx;
   ya[s]=yyy;

 

  s++;


}

fclose(read11);

dif=fabs(exis-xa[1]);

 for(i=1;i<=NPMODEL;i++){
   if((dift=fabs(exis-xa[i]))<dif) {
     ns=i;
     dif=dift;
   }

   c[i]=ya[i];
   d[i]=ya[i];

 }
 
 y=ya[ns--];


 return y;

}

/*******************************************************************/
/***** For Radius from Siess(2000) : calculates derivative of the radius at any time *****/
/*******************************************************************/

float polintradDEV(float exis)

{

 float xa[5000],ya[5000],c[5000],d[5000],deriv[5000],den,dif,dift,ho,hp,w;
 float viejin[5000];
 float y,dy,age1,radius,age2,k2cmk2r,xx,xxin,yy,yyin;
 int j=1,i,m,ns=1;
 float lum,mbol,Reff,teff,rhoeff,logg, masa,model,phase; 

 FILE *readm33,*writethis33;
 readm33 = fopen("m0.6z02.hrd","r");
 writethis33 = fopen("edad.dat","w");
  while(j<=NPMODEL)


{
  fscanf(readm33,"%f\t""%f\t""%f\t""%f\t""%f\t""%f\t""%f\t""%e\t""%e\t""%f\t""%f\n",&model, &phase, &lum, &mbol, &Reff, &radius, &teff, &rhoeff, &logg, &masa, &age1);
   xx=age1;
   yy=radius;


   xa[j]=xx;
   ya[j]=yy;
   viejin[j]=age1;

  // deriv[j]=((ya[j+ DV ]-ya[j-DV])/(xa[j+DV]-xa[j-DV]))*(tg/NSA);
  deriv[j]=((ya[j+ DV ]-ya[j-DV])/(xa[j+DV]-xa[j-DV]));  
    if(j==1){
    
    
 //   deriv[j]=((ya[j+DV]-ya[j])/(xa[j+DV]-xa[j]))*(tg/NSA);
 deriv[j]=((ya[j+DV]-ya[j])/(xa[j+DV]-xa[j]));
  }
  
   if(j==NPMODEL){
    
 //   deriv[j]=((ya[j]-ya[j-DV])/(xa[j]-xa[j-DV]))*(tg/NSA);
 deriv[j]=((ya[j]-ya[j-DV])/(xa[j]-xa[j-DV]));
  }
   //    printf("%d\t""%e\t""%f\n",j,xx,deriv[j]);
     fprintf(writethis33,"%6.3e\n",age1);

  j++;
}

 fclose(readm33);
 fclose(writethis33);

dif=fabs(exis-xa[1]);

 for(i=1;i<NPMODEL;i++){
   if((dift=fabs(exis-xa[i]))<dif) {
     ns=i;
     dif=dift;
   }

   c[i]=deriv[i];
   d[i]=deriv[i];
 }

 y=deriv[ns--];
 return y;
}

float kadosconv(float exis)

{

 float xa[5000],ya[5000],c[5000],d[5000],den,dif,dift,ho,hp,w;
 float y,dy,age1,radius,age2,k2cmk2r,xx,xxin,yy,yyin,cados;
 int j=1,i,m,ns=1;
 float model,phase,lum,mbol,Reff,teff,rhoeff,logg, masa,idade; 
 float tc, roc, etac,mcore, rcore, lnuc, lgrav, k2conv, k2rad,irad,iconv;
 FILE *readmd,*readthisd,*write3;
 
 // readm = fopen("k2ov1.csv","r"); /* 209 puntos */
  readmd = fopen("m0.6z02.var1","r"); 
  //ojo:
 readthisd = fopen("edad.dat","r"); 
 write3=fopen("mcore.txt","w"); 
 while(j<=NPMODEL)

{
  fscanf(readmd,"%f\t""%f\t""%f\t""%f\t""%f\t""%f\t""%f\t""%f\t""%f\t""%f\n",&model, &tc,&roc,&etac,&mcore,&rcore,&lnuc,&lgrav, &k2conv, &k2rad);
  
 fscanf(readthisd,"%f\n",&idade);
   xx=idade;
   yy=k2conv;
   xa[j]=xx;
   ya[j]=yy;

  j++;
}

 fclose(readmd);
 fclose(readthisd);
 fclose(write3);
dif=fabs(exis-xa[1]);
 for(i=1;i<NPMODEL;i++){
   if((dift=fabs(exis-xa[i]))<dif) {
     ns=i;
     dif=dift;
   }
   c[i]=ya[i];
   d[i]=ya[i];
}
 y=ya[ns--];
 return y;

}


float kadosrad(float exis)

{

 float xa[5000],ya[5000],c[5000],d[5000],den,dif,dift,ho,hp,w;
 float y,dy,age1,radius,age2,k2cmk2r,xx,xxin,yy,yyin,cados;
 int j=1,i,m,ns=1;
 float model,phase,lum,mbol,Reff,teff,rhoeff,logg, masa,idade; 
 float tc, roc, etac,mcore, rcore, lnuc, lgrav, k2conv, k2rad,irad,iconv;
 FILE *readm4,*readthis4,*write4;
 
 // readm = fopen("k2ov1.csv","r"); /* 209 puntos */
  readm4 = fopen("m0.6z02.var1","r"); 
  
  //ojo:
 readthis4 = fopen("edad.dat","r"); 
 write4=fopen("mcore.txt","w"); 
 while(j<=NPMODEL)


{
  fscanf(readm4,"%f\t""%f\t""%f\t""%f\t""%f\t""%f\t""%f\t""%f\t""%f\t""%f\n",&model, &tc,&roc,&etac,&mcore,&rcore,&lnuc,&lgrav, &k2conv, &k2rad);
  
cados = k2rad;
 fscanf(readthis4,"%f\n",&idade);
   xx=idade;
   yy=cados;
   xa[j]=xx;
   ya[j]=yy;
   fprintf(write4,"%6.3e\t""%6.3e\n",xx,mcore);

  j++;
}

 fclose(readm4);
 fclose(readthis4);
 fclose(write4);
dif=fabs(exis-xa[1]);
 for(i=1;i<NPMODEL;i++){
   if((dift=fabs(exis-xa[i]))<dif) {
     ns=i;
     dif=dift;
   }

   c[i]=ya[i];
   d[i]=ya[i];
 }
 y=ya[ns--];
 return y;

}


/*******************************************************************/
/***** For Radius from Siess(2000) : calculates derivative of the radius at any time *****/
/*******************************************************************/


float polintmascoreDEV(float exis)

{

 float xa[5000],ya[5000],c[5000],d[5000],deriv[5000],den,dif,dift,ho,hp,w;
 float viejin[5000];
 float y,dy,age1,radius,age2,k2cmk2r,xx,xxin,yy,yyin;
 int j=1,i,m,ns=1;
 float lum,mbol,Reff,teff,rhoeff,logg, masa,model,phase; 
 float te,masacore;
 FILE *readmy;
 readmy = fopen("mradi.txt","r");
 // writethis = fopen("edad.dat","w");
  while(j<=NPMODEL)

{
 fscanf(readmy,"%f\t""%f\n",&te, &masacore);
   xx=te;
   yy=masacore;

   xa[j]=xx;
   ya[j]=yy;

   deriv[j]=((ya[j+ DV ]-ya[j-DV])/(xa[j+DV]-xa[j-DV]));

   if(j==1){
    deriv[j]=((ya[j+DV]-ya[j])/(xa[j+DV]-xa[j]));
  }

   if(j==NPMODEL){
    deriv[j]=((ya[j]-ya[j-DV])/(xa[j]-xa[j-DV]));
  }
  j++;
}
 fclose(readmy);
 
 
dif=fabs(exis-xa[1]);
 for(i=1;i<NPMODEL;i++){
   if((dift=fabs(exis-xa[i]))<dif) {
     ns=i;
    dif=dift;
   }
   c[i]=deriv[i];
   d[i]=deriv[i];
 }
 y=deriv[ns--];
 return y;

}


float rcore(float exis)

{

 float xa[5000],ya[5000],c[5000],d[5000],den,dif,dift,ho,hp,w;
 float y,dy,age1,radius,age2,k2cmk2r,xx,xxin,yy,yyin,cados;
 int j=1,i,m,ns=1;
 float model,phase,lum,mbol,Reff,teff,rhoeff,logg, masa,idade; 
 float mconv1,rconv1,muc,renv,tenv,rhoenv;
 FILE *readm7,*readthis7;
 
 // readm = fopen("k2ov1.csv","r"); /* 209 puntos */
  readm7 = fopen("m0.6z02.var2","r"); 
 readthis7 = fopen("edad.dat","r"); 
 // write3=fopen("mcore.txt","w"); 
 while(j<=NPMODEL)

{fscanf(readm7,"%f\t""%f\t""%f\t""%f\t""%f\t""%f\t""%f\n",&model, &mconv1,&rconv1,&muc,&renv,&tenv,&rhoenv);
cados = renv;
 fscanf(readthis7,"%f\n",&idade);
   xx=idade;
   yy=cados;
   xa[j]=xx;
   ya[j]=yy;
  j++;
}
 fclose(readm7);
 fclose(readthis7);
dif=fabs(exis-xa[1]);

 for(i=1;i<NPMODEL;i++){
   if((dift=fabs(exis-xa[i]))<dif) {
     ns=i;
     dif=dift;
   }
   c[i]=ya[i];
   d[i]=ya[i];
 }
 y=ya[ns--];
 return y;
}



float turnovernata(float exis)

{

 float xa[5000],ya[5000],c[5000],d[5000],den,dif,dift,ho,hp,w;
 float y,dy,age1,radius,age2,k2cmk2r,xx,xxin,yy,yyin,cados;
 int j=1,i,m,ns=1;
 float model,natatime,to,loglumla,logteff,loggrav,taog,nrossby,masas;
 
 FILE *readm264;
 
 // readm = fopen("k2ov1.csv","r"); /* 209 puntos */
  readm264 = fopen("natalia06.dat","r"); 
  
  //ojo:

  while(j<=20)

{
  fscanf(readm264,"%f\t""%f\t""%f\t""%f\t""%f\t""%f\t""%f\t""%f\n",&masas,&natatime,&loglumla,&logteff,&loggrav,&to,&taog,&nrossby);
   xx=pow(10.0,natatime);
   yy=taog;
   xa[j]=xx;
   ya[j]=yy;
   
   printf("%6.3e\t""%6.3e\n",xx,yy);

  j++;
}

 fclose(readm264);

 dif=fabs(exis-xa[1]);
 for(i=1;i<20;i++){
   if((dift=fabs(exis-xa[i]))<dif) {
     ns=i;
     dif=dift;
   }

   c[i]=ya[i];
   d[i]=ya[i];
 }
 y=ya[ns--];
 return y;

}

   

float mrad(float exis)

{

 float xa[5000],ya[5000],c[5000],d[5000],den,dif,dift,ho,hp,w;
 float y,dy,age1,radius,age2,k2cmk2r,xx,xxin,yy,yyin,cados;
 int j=1,i,m,ns=1;
 float model,phase,lum,mbol,Reff,teff,rhoeff,logg, masa,idade; 
 float mconv1,rconv1,muc,renv,tenv,rhoenv;
 FILE *readm66,*readthis66,*write66;
 
  readm66 = fopen("m0.6z02.var2","r"); 
 readthis66 = fopen("edad.dat","r"); 
  write66=fopen("mradi.txt","w"); 
 while(j<=NPMODEL)

 {fscanf(readm66,"%f\t""%f\t""%f\t""%f\t""%f\t""%f\t""%f\n",&model, &mconv1,&rconv1,&muc,&renv,&tenv,&rhoenv);
cados = muc;
 fscanf(readthis66,"%f\n",&idade);
   xx=idade;
   yy=cados;
   xa[j]=xx;
   ya[j]=yy;
    fprintf(write66,"%6.3e\t""%6.3e\n",xx,yy);
  j++;
}
 fclose(readm66);
 fclose(readthis66);
 fclose(write66);
dif=fabs(exis-xa[1]);
 for(i=1;i<NPMODEL;i++){
   if((dift=fabs(exis-xa[i]))<dif) {
     ns=i;
     dif=dift;
   }
   c[i]=ya[i];
   d[i]=ya[i];
 }
 y=ya[ns--];
 return y;

}

float didtcore(float exis)

{

 float xa[5000],ya[5000],c[5000],d[5000],deriv[5000],den,dif,dift,ho,hp,w;
 float viejin[5000];
 float y,dy,age1,radius,age2,k2cmk2r,xx,xxin,yy,yyin;
 int j=1,i,m,ns=1,model;
 float lum,mbol,Reff,teff,rhoeff,logg, masa,phase; 
 float te,masacore,idade;
 float tc,roc,etac,mcore,rcore,lnuc,lgrav, k2conv, k2rad;
 FILE *readm4p, *readthisplis;

 
 
 readthisplis = fopen("edad.dat","r");
 readm4p = fopen("m0.6z02.var1","r"); 

 while(j<=NPMODEL)

  {
  
  fscanf(readm4p,"%d\t""%e\t""%e\t""%f\t""%f\t""%f\t""%f\t""%f\t""%f\t""%f\n",&model, &tc,&roc,&etac,&mcore,&rcore,&lnuc,&lgrav, &k2conv, &k2rad);
   fscanf(readthisplis,"%f\n",&idade);

   yy = k2conv;
xx=idade;
  j++;
}
fclose(readm4p);
fclose(readthisplis);
dif=fabs(exis-xa[1]);
 for(i=1;i<NPMODEL;i++){
   if((dift=fabs(exis-xa[i]))<dif) {
     ns=i;
     dif=dift;
   }
   c[i]=deriv[i];
   d[i]=deriv[i];
 }
 y=deriv[ns--];

 return y;

}

float polintradcoreDEV(float exis)
/*******************************************************************/
/***** For Radius from Siess(2000) : calculates derivative of the radius at any time *****/
/*******************************************************************/
{
 
 float viejin[5000];
 float xa[5000],ya[5000],c[5000],d[5000],deriv[5000],den,dif,dift,ho,hp,w;
 float y,dy,age1,radius,age2,k2cmk2r,xx,xxin,yy,yyin,cados;
 int j=1,i,m,ns=1;
 float model,phase,lum,mbol,Reff,teff,rhoeff,logg, masa,idade; 
 float mconv1,rconv1,muc,renv,tenv,rhoenv;
 FILE *readm1,*readthis1,*writethis1;
 
 // readm = fopen("k2ov1.csv","r"); /* 209 puntos */
  readm1 = fopen("m0.6z02.var2","r"); 
 readthis1 = fopen("edad.dat","r"); 
 // write3=fopen("mcore.txt","w"); 
 while(j<=NPMODEL)

{fscanf(readm1,"%f\t""%f\t""%f\t""%f\t""%f\t""%f\t""%f\n",&model, &mconv1,&rconv1,&muc,&renv,&tenv,&rhoenv);
cados = renv;
 fscanf(readthis1,"%f\n",&idade);
   xx=idade;
   yy=cados;
   xa[j]=xx;
   ya[j]=yy;
 
   
   
   // deriv[j]=((ya[j+ DV ]-ya[j-DV])/(xa[j+DV]-xa[j-DV]))*(tg/NSA);
  deriv[j]=((ya[j+ DV ]-ya[j-DV])/(xa[j+DV]-xa[j-DV]));  
    if(j==1){
    
    
 //   deriv[j]=((ya[j+DV]-ya[j])/(xa[j+DV]-xa[j]))*(tg/NSA);
 deriv[j]=((ya[j+DV]-ya[j])/(xa[j+DV]-xa[j]));
  }
  
   if(j==NPMODEL){
    
 //   deriv[j]=((ya[j]-ya[j-DV])/(xa[j]-xa[j-DV]))*(tg/NSA);
 deriv[j]=((ya[j]-ya[j-DV])/(xa[j]-xa[j-DV]));
  }
   //    printf("%d\t""%e\t""%f\n",j,xx,deriv[j]);
  //  fprintf(writethis1,"%6.3e\n",age1);
//
  j++;
}
   
    fclose(readthis1);
   fclose(readm1);
   
   
   
   
dif=fabs(exis-xa[1]);
 for(i=1;i<NPMODEL;i++){
   if((dift=fabs(exis-xa[i]))<dif) {
     ns=i;
    dif=dift;
   }
   c[i]=deriv[i];
   d[i]=deriv[i];
 }
 y=deriv[ns--];
 return y;

}
 
float litio7(float exis)

{

 float xa[5000],ya[5000],c[5000],d[5000],den,dif,dift,ho,hp,w;
 float y,dy,age1,radius,age2,k2cmk2r,xx,xxin,yy,yyin,cados;
 int j=1,i,m,ns=1;
 float model,surfacemassfraction2H,He3,Li6,Li7,Be9,B10,B11,idade; 
 FILE *readm24,*readthis24,*write24;
 
 // readm = fopen("k2ov1.csv","r"); /* 209 puntos */
  readm24 = fopen("m0.6z02.xsurf","r"); 
  
  //ojo:
 readthis24 = fopen("edad.dat","r"); 
 write24=fopen("mcore.txt","w"); 
 while(j<=NPMODEL)


{
  fscanf(readm24,"%f\t""%e\t""%e\t""%e\t""%e\t""%e\t""%e\t""%e\n",&model,&surfacemassfraction2H,&He3,&Li6,&Li7,&Be9,&B10,&B11);
  
cados = Li7;
 fscanf(readthis24,"%f\n",&idade);
   xx=idade;
   yy=cados;
   xa[j]=xx;
   ya[j]=yy;
//   fprintf(write24,"%6.3e\t""%6.3e\n",xx,Li7);

  j++;
}

 fclose(readm24);
 fclose(readthis24);
 fclose(write24);
dif=fabs(exis-xa[1]);
 for(i=1;i<NPMODEL;i++){
   if((dift=fabs(exis-xa[i]))<dif) {
     ns=i;
     dif=dift;
   }

   c[i]=ya[i];
   d[i]=ya[i];
 }
 y=ya[ns--];
 return y;

}

float labolo(float jk)

{
	float visual[170],exit[170],uband[170],mod_d[170];
	
	float uvkh,vikh,vickh,bvkh,den,dif,dift,ho,hp,w,bv;
	float xa[1000],ya[1000],c[1000],d[1000];
	
	float y,dy;
	float h,citauno=0.999,citainicial=0,xx,yy;
	float bc, tempe, vic,vij,vj,vh,vk;
	int i,m,ns=1;
	char id [40];
	int j=1,k=1;
	
	FILE *readss5;
	
	readss5 = fopen("KH95_BC_TEFF_VI.dat","r");
	
	for(j=1;j<=55;j++)
		
	{
		
		fscanf(readss5,"%f\t""%f\t""%f\t""%f\t""%f\t""%f\t""%f\t""%f\n",&bc, &tempe, &vic, &vij, &vj, &vh, &vk, &bv);
		//  fprintf(write,"%6.3f\t""%6.3f\n",-vh+vk,-vj+vh);
	//	printf("%6.3f\t""%6.3f\n",-vj+vk,bv);
		xa[j]=bc;  
		ya[j]=bv;   
		
	}
	
	
	dif=fabs(jk-xa[1]);
	
	for(i=1;i<55;i++){
		if((dift=fabs(jk-xa[i]))<dif) {
			ns=i;
			dif=dift;
		}
		
		c[i]=ya[i];
		d[i]=ya[i];
		
	}
	
	y=ya[ns--];
	
	
	fclose(readss5);
	
	return y;
}







   


// End of the program
// Giovanni Pinzon Estrada gapinzone@unal.edu.co
// 2021
