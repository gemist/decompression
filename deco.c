/*
Program uses Buhlman ZH-L16 algorithm and Baker's gradient factors.
Program is written for education purpose only. Originally written in Fortran by Baker

Things that are (not) included:
- no repetitive dives
- only dives in the sea water
- oxygen toxicity is not included
- ... 
Missing things will be included in the future. Program is free 
and further modifications and improvements are welcome. 
For further questions contact me.

Last modified: 1. november 2010
email: alan.bizjak@gmail.com

#example of input.txt file:
### START file
# SAMPLE DIVE TO 90 METERS OF SEAWATER GAUGE (MSWG) FOR 20 MINUTES
# 
# Everything behind character # is comment and it is ignored. Comments 
# can be only between sections GASES, DIVE, DECO
# 
# GASES: firstly you specify the number of gases that are going 
# to be used during the dive. After that are fractions of gases (f) 
# separeted by comma in following order fO2, fHe, fN2 
# 
GASES
4
0.13, 0.50, 0.37
0.36,0.00,0.64
0.50,0.00,0.50
0.80,0.00,0.20
#
# DIVE: 
#  first row starting_depth=0 meters final_depth=90 rate=23 meter/min number of gas mixture 1
#  
DIVE
0 90 23 1
90 20 1
# 
DECO
1 -10 3 0.75 0.3
33
2 -10 3
21
3 -10 3
9
4 -10 3
### END of file

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void ascdec(float sdepth, float fdepth,float rate, int mixnum);
void cdepth(float depth, float srtime,int mixnum);
float safasc(void);
float Mvcalc(float depth);
void decostop(float stopd,float nxstop,int mixnum);
float deepest_deco_stop(float depth,float rate);
void initialize(void);

const float eps=1E-4;
float pH2O=0.567; 
float pCO2=0.534;
int nummix,mixnum;
float factor;
//float fO2[10];
//float fHe[10];
//float fN2[10];
float *fO2;
float *fHe;
float *fN2;
float pHe[16],pN2[16];
float kHe[16];
float kN2[16];
float aHe[16];
float bHe[16];
float aN2[16];
float bN2[16];
float halftHe[16];
float halftN2[16];
float runtime,segtime;
int segnum;
int main(int argc, char *argv[]){
 
  initialize();
 
  char line[128];
  
  runtime =0.0;
  segnum = 0;
  
  FILE *outFile=fopen("output.txt","wt");
  if (outFile == NULL) {
    printf("Error opening file output.txt\n");
    exit(1);
  }
  
  int i;
  FILE *inFile=fopen("input.txt","rt");
  if (inFile == NULL) {
    printf("Error opening file input.txt\n");
    exit(1);
  }
  float depth,sdepth,fdepth,rate;
  while (fgets(line,sizeof(line),inFile)!=NULL){
      if (line[0]!='#'){ //skip comments 
	char tmp[8];
	sscanf(line,"%d",&nummix);
	//allocate dynamic memory
	  fO2 = (float *)malloc(sizeof(float)*nummix); 
	  fHe = (float *)malloc(sizeof(float)*nummix); 
	  fN2 = (float *)malloc(sizeof(float)*nummix); 
	  for (i=0;i<nummix;i++){
	    fgets(line, sizeof(line), inFile);
	    sscanf(line,"%f,%f,%f",&fO2[i],&fHe[i],&fN2[i]);
	    float chksum = fO2[i]+fHe[i]+fN2[i];
	    
	    if (chksum != 1.0f){ 
	      printf("Error: In input file sum of " 
		     "coefficients must be 1.0\n");
	      exit(1);
	    }
	  }
	  
	  fprintf(outFile,"DISCLAIMER:   anybody who uses this program to calculate actual dives takes\n"
		  "full responsability for their own actions and understand/accept all the risks\n"
		  "inherent in technical diving.\n"
		  "This dive schedule is experimental, under no circumstances there is any\n"
		  "implication that is safe or that you wont get bent, use it at your own risk.\n\n\n");
	  
	  fprintf(outFile,"DIVE PROFILE\n\n");
	  fprintf(outFile,"%4s %5s %5s | %6s | %7s %6s  %6s %9s %8s\n",
		  "Seg-","Segm.","Run","Gasmix","Ascent","From","To",
		  "Rate","Constant");
	  fprintf(outFile,"%4s %5s %5s | %6s | %7s %6s  %6s %9s %8s\n",
		  "ment","Time","Time","Used", "or","Depth","Depth",
		  "+Dn/-Up","Depth");
	  fprintf(outFile,"%4s %5s %5s | %6s | %7s %6s  %6s %9s %8s\n",
		  "#","(min)","(min)","#","Descent","(mswg)","(mswg)",
		  "(msw/min)","(mswg)");
	  fprintf(outFile,"%4s %5s %5s | %6s | %7s %6s  %6s %9s %8s\n",
		  "----","----","----","----","----","----","----",
		  "----","----");

	  char profil[8];
	  do {
	    fgets(line, sizeof(line), inFile);
	    sscanf(line,"%s",profil);
	    if (strcmp(profil,"ascent")==0 || strcmp(profil,"descent")==0){ 
	      fgets(line, sizeof(line), inFile);
	      sscanf(line,"%f %f %f %d",&sdepth,&fdepth,&rate,&mixnum);
	      mixnum--;
	      ascdec(sdepth,fdepth,rate,mixnum);
	      fprintf(outFile,"%-4d %-5.1f %5.1f | %6d | %7s %6.1f %6.1f %9.1f  %8s\n",segnum,segtime,runtime,mixnum+1,profil,sdepth,fdepth,rate," ");
	    }else if (strcmp(profil,"const_depth")==0){
	      float srtime;
	    fgets(line, sizeof(line), inFile);
	    sscanf(line,"%f %f %d",&depth,&srtime,&mixnum);
	    mixnum--;
	    cdepth(depth,srtime,mixnum);
	    fprintf(outFile,"%-4d %5.1f %5.1f | %6d | %7s %6s %6s %9s %8.1f   \n",segnum,segtime,runtime,mixnum+1," "," "," "," ",depth);
	    }else if (strcmp(profil,"deco")==0){

 fprintf(outFile,"\nDECOMPRESSION PROFILE\n\n");
	  fprintf(outFile,"%4s %5s %5s | %6s | %7s %9s  %5s | %6s %5s %5s %8s\n",
		  "Seg-","Segm.","Run","Gasmix","Ascent","Ascent","Max",
		  "DECO","STOP","RUN","Gradient");
	  fprintf(outFile,"%4s %5s %5s | %6s | %7s %9s  %5s | %6s %5s %5s %8s\n",
		  "ment","Time.","Time","Used","To","Rate","%M-",
		  "STOP","TIME","TIME","Factor");
	  fprintf(outFile,"%4s %5s %5s | %6s | %7s %9s  %5s | %6s %5s %5s %8s\n",
		  "#","(min)","(min)","#","(mswg)","(msw/min)","Value",
		  "(mswg)","(min)","(min)","(GF)");
	  fprintf(outFile,"%4s %5s %5s | %6s | %7s %9s  %5s | %6s %5s %5s %8s\n",
		  "---","---","---","---","---","---","---",
		  "---","---","---","---");
	  
	    }else{
	      printf("Error: word for profil must be ascent, descent, const_depth or deco\n");
	      exit(1);
	    }
	      

	  }while(strcmp(profil,"deco")!=0);
	  
	  

	  float triald;
	  float stepsz,fctrhi,fctrlo,change,change_time;
	  fgets(line, sizeof(line), inFile);
	  sscanf(line,"%d %f %f %f %f",&mixnum,&rate,&stepsz,&fctrhi,&fctrlo);
	  mixnum--;
	  fgets(line, sizeof(line), inFile);
	  sscanf(line,"%f %f",&change,&change_time);

	  printf("deepest deco stop %f\n",deepest_deco_stop(depth,rate));	  
	  factor=fctrlo;
	  float fctrsl=0.0;
	  float triald_tmp,stopd;
	  do{
	    triald = fdepth - stepsz;
	    if (triald <= change) {triald_tmp = triald; triald =change;} 
	    float ceiling=safasc();
	    if (ceiling > triald){
	      float nxstop = triald;
	      if (triald==change)
	      	stopd = triald_tmp + stepsz;
	      else
		stopd = triald+stepsz;
	      if (stopd > 0.0f && fctrsl ==0.0f){
		fctrsl = (fctrhi-fctrlo)/(0.0-stopd);
	      }
	      float stopgf = factor;
	      factor = nxstop*fctrsl+fctrhi;
	      
	      decostop(stopd,nxstop,mixnum);
	      fprintf(outFile,"%-4d %5.1f %5.1f | %6d | %7s %9s  %5s | %6.1f %5.1f %5.1f %8.3f\n",
		      segnum,segtime,runtime,mixnum+1," "," "," ",stopd,segtime,runtime,stopgf);
	      
	    }
	    
	    sdepth = fdepth;
	    fdepth = triald;
	    
	    ascdec(sdepth,fdepth,rate,mixnum);
	    float pmvmax=Mvcalc(fdepth);
	    fprintf(outFile,
		    "%-4d %5.1f %5.1f | %6d | %7.1f %9.1f  %5.3f | %6s %5s  %5s %8s\n",
		    segnum,segtime,runtime,mixnum+1,triald,rate,pmvmax," "," "," "," ");
	    	    
	    if (change==triald) {
	      //cdepth(change,change_time,mixnum);
	      fgets(line, sizeof(line), inFile);
	      sscanf(line,"%d %f %f",&mixnum,&rate,&stepsz);
	      mixnum--;
	      fgets(line, sizeof(line), inFile);
	      sscanf(line,"%f",&change);
	    }
	    
	  }while(triald > 0.0);
	  
      }     

    }


  fclose(inFile);
  
  
  fclose(outFile);
  free(fO2);
  free(fHe);
  free(fN2);

  return 0;
}

void initialize(void){
  //haltime constants for tisues
  halftN2[0]=5.0;      halftHe[0]=1.88;   
  halftN2[1]=8.0;      halftHe[1]=3.02;  
  halftN2[2]=12.5;     halftHe[2]=4.72;  
  halftN2[3]=18.5;     halftHe[3]=6.99;  
  halftN2[4]=27.0;     halftHe[4]=10.21; 
  halftN2[5]=38.3;     halftHe[5]=14.48; 
  halftN2[6]=54.3;     halftHe[6]=20.53; 
  halftN2[7]=77.0;     halftHe[7]=29.11; 
  halftN2[8]=109.0;    halftHe[8]=41.2;  
  halftN2[9]=146.0;    halftHe[9]=55.19; 
  halftN2[10]=187.0;   halftHe[10]=70.69; 
  halftN2[11]=239.0;   halftHe[11]=90.34; 
  halftN2[12]=305.0;   halftHe[12]=115.29;
  halftN2[13]=390.0;   halftHe[13]=147.42;
  halftN2[14]=498.0;   halftHe[14]=188.24;
  halftN2[15]=635.0;   halftHe[15]=240.02;
  
  // constants 
  aHe[0]=16.189; aHe[1]=13.830; aHe[2]=11.919; aHe[3]=10.458;
  aHe[4]=9.220;  aHe[5]=8.205;  aHe[6]=7.305;  aHe[7]=6.502;
  aHe[8]=5.950;  aHe[9]=5.545;  aHe[10]=5.333; aHe[11]=5.189; 
  aHe[12]=5.181; aHe[13]=5.176; aHe[14]=5.172; aHe[15]=5.119;
  
  bHe[0]=0.4770;  bHe[1]=0.5747;  bHe[2]=0.6527;  bHe[3]=0.7223;
  bHe[4]=0.7582;  bHe[5]=0.7957;  bHe[6]=0.8279;  bHe[7]=0.8553;
  bHe[8]=0.8757;  bHe[9]=0.8903;  bHe[10]=0.8997; bHe[11]=0.9073; 
  bHe[12]=0.9122; bHe[13]=0.9171; bHe[14]=0.9217; bHe[15]=0.9267;

  aN2[0]=11.696;  aN2[1]=10.000;  aN2[2]=8.618;   aN2[3]=7.562;
  aN2[4]=6.667;   aN2[5]=5.600;   aN2[6]=4.947;   aN2[7]=4.500;
  aN2[8]=4.187;   aN2[9]=3.798;   aN2[10]=3.497;  aN2[11]=3.223;
  aN2[12]=2.850;  aN2[13]=2.737;  aN2[14]=2.523;  aN2[15]=2.327;
  
  bN2[0]=0.5578;  bN2[1]=0.6514;  bN2[2]=0.7222;   bN2[3]=0.7825;
  bN2[4]=0.8126;  bN2[5]=0.8434;  bN2[6]=0.8693;   bN2[7]=0.8910;
  bN2[8]=0.9092;  bN2[9]=0.9222;  bN2[10]=0.9319;  bN2[11]=0.9403;
  bN2[12]=0.9477; bN2[13]=0.9544; bN2[14]=0.9602;  bN2[15]=0.9653;


  int i;
  for (i=0;i<16;i++){
    kHe[i]=log(2.0f)/halftHe[i];
    kN2[i]=log(2.0f)/halftN2[i];
    pHe[i]=0.0;
    pN2[i]=7.452;
  }


}


void ascdec(float sdepth, float fdepth,float rate, int mixnum){ //ascent or descent calculations
  segtime=(fdepth-sdepth)/rate; //segment time
  runtime += segtime;
  segnum++;
  float sp_amb = sdepth+10.0; //start pressure
  // float fp_amb = fdepth+10.0; //final pressure
  float pHe_alv=(sp_amb-pH2O)*fHe[mixnum];
  float pN2_alv=(sp_amb-pH2O)*fN2[mixnum];
  float He_rate=rate*fHe[mixnum];
  float N2_rate=rate*fN2[mixnum];
  
  int i;
  float pHe0[16],pN20[16];
  for (i=0;i<16;i++){
    pHe0[i]=pHe[i];
    pN20[i]=pN2[i];

    pHe[i]=pHe_alv + He_rate*(segtime-1.0/kHe[i])
      - (pHe_alv-pHe0[i]-He_rate/kHe[i])*exp(-kHe[i]*segtime);

    pN2[i]=pN2_alv + N2_rate*(segtime-1.0/kN2[i])
      - (pN2_alv-pN20[i]-N2_rate/kN2[i])*exp(-kN2[i]*segtime);

  }
}  


void cdepth(float depth, float srtime,int mixnum){ //constant depth

  segtime = srtime-runtime;
  runtime+=segtime;
  segnum++;
  float p_amb= depth +10.0;
  float pHe_alv=(p_amb-pH2O)*fHe[mixnum];
  float pN2_alv=(p_amb-pH2O)*fN2[mixnum];
  
  int i;
  float pHe0[16],pN20[16];
  for (i=0;i<16;i++){
    pHe0[i]=pHe[i];
    pN20[i]=pN2[i];
    pHe[i]=pHe0[i] + (pHe_alv -pHe0[i])*(1.0-exp(-kHe[i]*segtime));
    pN2[i]=pN20[i] + (pN2_alv -pN20[i])*(1.0-exp(-kN2[i]*segtime));
  }
}

float safasc(void){ //safe ascent

  float ceiling = 0.0;
  int i;
  float pHeN2[16],aHeN2[16],bHeN2[16],p_amb[16];
  float safeadd[16];
  for (i=0;i<16;i++){
    pHeN2[i]=pHe[i]+pN2[i];
    aHeN2[i]=(pHe[i]*aHe[i]+pN2[i]*aN2[i])/pHeN2[i];
    bHeN2[i]=(pHe[i]*bHe[i]+pN2[i]*bN2[i])/pHeN2[i];
    p_amb[i]=(pHeN2[i]-aHeN2[i]*factor)/(factor/bHeN2[i]-factor+1.0);
    safeadd[i]=p_amb[i]-10.0;
    if (safeadd[i] > ceiling) ceiling =safeadd[i];
  }
  return ceiling;
}

float Mvcalc(float depth){ //M value
  float pmvmax =0.0;
  int i;
  float Mvalue[16],percMv[16];
  float pHeN2[16],aHeN2[16],bHeN2[16];
  float p_amb = depth+10.0;
  for (i=0;i<16;i++){
    pHeN2[i]=pHe[i]+pN2[i];
    aHeN2[i]=(pHe[i]*aHe[i]+pN2[i]*aN2[i])/pHeN2[i];
    bHeN2[i]=(pHe[i]*bHe[i]+pN2[i]*bN2[i])/pHeN2[i];
    
    Mvalue[i]=p_amb/bHeN2[i]+aHeN2[i];
    percMv[i]=pHeN2[i]/Mvalue[i];
    if (percMv[i] > pmvmax) pmvmax=percMv[i];
  }
  return pmvmax;
}


void decostop(float stopd,float nxstop,int mixnum){ //decostop

  float ceiling;
  //  float round = roundf(runtime + 0.5);
  //segtime = round - runtime;
  segtime=1.0;
  float p_amb = stopd + 10.0;
  float pHe_alv=(p_amb-pH2O)*fHe[mixnum];
  float pN2_alv;
  
  pN2_alv=(p_amb-pH2O)*fN2[mixnum];
  int i;
  float pHe0[16],pN20[16];

  for (i=0;i<16;i++){
  pHe0[i]=pHe[i];
  pN20[i]=pN2[i];
      
  }

  do {
    for (i=0;i<16;i++){
      pHe[i]=pHe0[i]+(pHe_alv-pHe0[i])*(1.0-exp(-kHe[i]*segtime));
      pN2[i]=pN20[i]+(pN2_alv-pN20[i])*(1.0-exp(-kN2[i]*segtime));
    }
    ceiling=safasc();
    if (ceiling > nxstop){
      segtime+=1.0;
    }
  }while(ceiling > nxstop);
  runtime+=segtime;
  segnum++;
}

float deepest_deco_stop(float depth,float rate){
  /*
    Function calculate deepest deco stop. Deepest deco stop is defined
    as a depth where pressure in one compartment equals ambient
    pressure.  Here the decompression begins. This value tells the
    diver that at this dpth he needs to break the ascent. The starting
    deco depth is calculated using bisection method. With bisection
    method we are searching for the time during ascent when this happen.
    One bound is set to zero and the other to at time when it would take 
    to zero ambient pressure (absolute).
   */
  float first_time,second_time,mid_time;
  float deepest_deco =0.0; 
  float p_amb = depth +10.0;
  float pHe_alv=(p_amb-pH2O)*fHe[mixnum];
  float pN2_alv=(p_amb-pH2O)*fN2[mixnum];
  float He_rate=rate*fHe[mixnum];
  float N2_rate=rate*fN2[mixnum];
  int i;
  for (i=0;i<16;i++){
    float pHe0=pHe[i];
    float pN20=pN2[i];
    float pHeN20=pHe0+pN20-p_amb;
    first_time = 0.0;
    second_time = -p_amb/rate; //multiply by -1 (negative rate)
  
    float pHe=pHe_alv + He_rate*(second_time-1.0/kHe[i])
      - (pHe_alv-pHe0-He_rate/kHe[i])*exp(-kHe[i]*second_time);

    float pN2=pN2_alv + N2_rate*(second_time-1.0/kN2[i])
      - (pN2_alv-pN20-N2_rate/kN2[i])*exp(-kN2[i]*second_time);
    float pHeN2 = pHe+pN2;
    if (pHeN20*pHeN2 > 0.0){ //bisection method fails in that case
      printf("Error in deepest decompression stop calculation. Root is not withing brackets\n");
      exit(1);
    }

  //start bisection method
  do{
    mid_time = (first_time+second_time)/2.0;
    float mid_pHe =pHe_alv + He_rate*(mid_time-1.0/kHe[i])
      - (pHe_alv-pHe0-He_rate/kHe[i])*exp(-kHe[i]*mid_time);
    float mid_pN2=pN2_alv + N2_rate*(mid_time-1.0/kN2[i])
      - (pN2_alv-pN20-N2_rate/kN2[i])*exp(-kN2[i]*mid_time);
    float mid_pHeN2=mid_pHe+mid_pN2-(p_amb+rate*mid_time);
    if ( mid_pHeN2< 0.0){
       pHeN20 = mid_pHeN2;
       first_time = mid_time;
    }else{
       pHeN2 = mid_pHeN2;
       second_time=mid_time;
    }
  }while(fabs(pHeN20-pHeN2)>eps);
  float cpt_deepest_deco = p_amb+rate*mid_time-10.0;
  if (cpt_deepest_deco > deepest_deco) deepest_deco=cpt_deepest_deco; 
  }
  return deepest_deco;
}

