/*
 * Copyright (c) 2022  Yogesh Maan <ymaan@ncra.tifr.res.in>
 * 
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of version 2 of the GNU General Public License as
 * published by the Free Software Foundation.
 * 
 * This program is distributed in the hope that it would be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * 
 * Further, this software is distributed without any warranty that it is
 * free of the rightful claim of any third person regarding infringement
 * or the like.  Any license provided herein, whether implied or
 * otherwise, applies only to this software file.  Patent licenses, if
 * any, provided herein do not apply to combinations of this program with
 * other software, or any other product whatsoever.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write the Free Software Foundation, Inc., 59
 * Temple Place - Suite 330, Boston MA 02111-1307, USA.
 * 
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "pcBeam_gmrt.h"
#include <math.h>
#include <unistd.h>
//#include "utils.h"



int help_required(char *string)
{
  if (strings_equal(string,"--help")) return(1);
  if (strings_equal(string,"-h")) return(1);
  return(0);
}

void pcbeam_gmrt_help()
{
  puts("");
  puts("pcBeam_gmrt - Construct post-correlation beam using GMRT's IA and PA/CDPA beam data\n");
  puts("usage: pcBeam_gmrt -{options} \n");
  puts("options:\n");
  puts(" -f1 <file-name> : GMRT data file1 (PA/CDPA)");
  puts(" -f2 <file-name> : GMRT data file2 (IA)");
  //puts("            ***(result will be f1-f2)***");
  puts(" ");
  //puts(" -block <bytes>  : block-size to read at once in bytes");
  puts(" -nch1 <nchan>   : No. of channels in file1 ");
  puts(" -ts1  <tsamp>   : Sampling time (in ms) in file1 ");
  puts(" -nch2 <nchan>   : No. of channels in file2 ");
  puts(" -ts2  <tsamp>   : Sampling time (in ms) in file2 ");

//  puts("-t numsamps - number of time samples to fetch (def=1024)");
//  puts("-T dur.(s)  - duration to fetch ");
//  puts("-gm hdr-file- GMRT header info file ");
//  puts("-n numbits  - specify output number of bits (def=input-size)");
//  puts("-white      - whitten each spectrum");
  puts(" -scale <fact>   : scale input-2 by 'fact' before subtraction ");
  puts("-o filename - specify output filename (def=stdout)");
//  puts("-headerless - do not broadcast resulting header (def=broadcast)");
  puts("");
  puts(" -startch <nch>   : Starting chan-num to be output ");
  puts(" -outnch  <nch>   : Total channels to be output ");
  puts("");
}




int file_exists(char *filename)
{
  if ((fopen(filename,"rb"))==NULL) { return(0);}
  else { return(1); }
}

FILE *open_file(char *filename, char *descriptor)
{
  FILE *fopen(), *fptr;
  if ((fptr=fopen(filename,descriptor)) == NULL) {
    fprintf(stderr,"Error in opening file: %s\n",filename);
    exit(1);
  }
  return fptr;
}

void error_message(char *message)
{
  fprintf(stderr,"ERROR: %s\n",message);
  exit(1);
}
//==============
//============== robust median rms ===============================
/* TO COMPUTE MEDIAN and RMS OF 'arr' BY EXCLUDING WHAT MAY BE
SOME CONTRIBUTION FROM INTERFERENCE */
void medrms(float *arr, int np)
{
    int i,iter,maxiter,k;
    long int nn;
    float *barr,med,rms,rms0,diff,atemp,an;

    nn = np;
    barr=(float *) malloc(np*sizeof(float));
    for (i=0;i<np;i++) barr[i] = arr[i];
    med = median(barr,nn);
    free(barr);

    atemp = 0.0;
    for (i=0;i<np;i++) atemp = atemp + (arr[i]-med)*(arr[i]-med);
    rms = sqrt(atemp/np);
    iter=0;
    maxiter=10;
    rms0=rms*2.0;
    while (iter<maxiter && fabs((rms0/rms)-1.0) > 0.05){
      rms0 = rms;
      atemp = 0.0;
      an = 0.0;
      for (i=0;i<np;i++){
        diff = arr[i]-med;
        if(fabs(diff) <= 3.5*rms0){
          an = an + 1.0;
          atemp = atemp + diff*diff;
        }
      }
      rms = sqrt(atemp/an);
      iter = iter+1;
    }
    arr[np] = med;
    arr[np+1] = rms;

}
//============================================================================

/* read in the general header info and timestamp for GMRT-obs from text files*/
long double gm_mjd(char gmhdrfile[]) /* includefile */
{
  char string[80], str1[32],str2[32],str3[32];
  int hh,mm,day,month,year,nyear,istat;
  long double seconds,mjd_fracday;
  long double mjd_day,gm_mjd;
  FILE *timefile;

  timefile = fopen(gmhdrfile,"rb");
  fscanf(timefile, "%s %s %s %s\n", string,string,string,string);
  fscanf(timefile, "%s %s %d:%d:%Lf  \n", str1,str2,&hh,&mm,&seconds);
  fscanf(timefile, "%s  %d:%d:%d  \n", str1,&day,&month,&year);
  fclose(timefile);

  slaCldj(year,month,day, &mjd_day, &istat);
  mjd_fracday = (hh + (mm + (seconds / 60.0)) / 60.0) / 24.0;
  gm_mjd = mjd_day + mjd_fracday - 5.5/24.0 ;  // correct for 5.5 hours diff between IST and UT
  //printf("mjd_day, mjd_fracday, tstart  %lf  %.15Lf  %.15lf\n",mjd_day, mjd_fracday, tstart);

  return gm_mjd;
}

void slaCldj ( int iy, int im, int id, long double *djm, int *j ) /*includefile*/
/*
**  - - - - - - - -
**   s l a C l d j
**  - - - - - - - -
**  Gregorian calendar to Modified Julian Date.
**  Given:
**     iy,im,id     int    year, month, day in Gregorian calendar
**  Returned:
**     *djm         double Modified Julian Date (JD-2400000.5) for 0 hrs
**     *j           int    status:
**                           0 = OK
**                           1 = bad year   (MJD not computed)
**                           2 = bad month  (MJD not computed)
**                           3 = bad day    (MJD computed)
**  The year must be -4699 (i.e. 4700BC) or later.
**  The algorithm is derived from that of Hatcher 1984 (QJRAS 25, 53-55).
**  Last revision:   29 August 1994
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   long iyL, imL;
/* Month lengths in days */
   static int mtab[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
/* Validate year */
   if ( iy < -4699 ) { *j = 1; return; }
/* Validate month */
   if ( ( im < 1 ) || ( im > 12 ) ) { *j = 2; return; }
/* Allow for leap year */
   mtab[1] = ( ( ( iy % 4 ) == 0 ) &&
             ( ( ( iy % 100 ) != 0 ) || ( ( iy % 400 ) == 0 ) ) ) ?
             29 : 28;
/* Validate day */
   *j = ( id < 1 || id > mtab[im-1] ) ? 3 : 0;
/* Lengthen year and month numbers to avoid overflow */
   iyL = (long) iy;
   imL = (long) im;
/* Perform the conversion */
   *djm = (double)
        ( ( 1461L * ( iyL - ( 12L - imL ) / 10L + 4712L ) ) / 4L
        + ( 306L * ( ( imL + 9L ) % 12L ) + 5L ) / 10L
        - ( 3L * ( ( iyL - ( 12L - imL ) / 10L + 4900L ) / 100L ) ) / 4L
        + (long) id - 2399904L );
}



//============================================================================
void main (int argc, char *argv[])
{
  int i, j, k, jt1,jt2, nc, headersize, headerless=0,gm=0;
  int ii,jj,kk, nchans1, nchans2, nchans, itemp, startch, outnch;
  float afact,anorm=-9999.0, med_fact,std;
  long double tsamp1, tsamp2, tsamp, mjd1, mjd2, mjd0, dtemp;
  float tstart=0.0,tend,dur1,dur2,dt1,dt2;
  char string[128],infoFile[128],ctt;
  float *fblock, *fblock1, *fblock2, *tempblock, min,max, atemp;
  unsigned short *sblock;
  unsigned char  *cblock;
  int ns=0, ns1=0, ns2=0,nsblk,nsblk1,nsblk2,t1,t2,f1,f2,opened=0,nout,iter;
  FILE *hdr1, *hdr2, *infohdr;
  

  /* set up default global variables */
  obits=headerless=naddt=nsamp=off1=off2=bsize=nchans1=nchans2=tsamp1=tsamp2=0;
  startch=outnch=0;
  //strcpy(inpfile,"stdin");
  //output=stdout;
  strcpy(outfile,"dummy");
  strcpy(inpfile1,"dummy");
  strcpy(inpfile2,"dummy");

  if (argc >= 1) {
    i=1;
    while (i<argc) {
      if (strings_equal(argv[i],"-o")) {
	i++;
        strcpy(outfile,argv[i]);
	output=fopen(outfile,"wb");
      } else if (strings_equal(argv[i],"-f1")) {
	i++;
        if((file_exists(argv[i]))) {
          strcpy(inpfile1,argv[i]);
          input1=open_file(inpfile1,"rb");
        } else {
          printf("specified input file does not exist!\n");
          exit(0);
        }
      } else if (strings_equal(argv[i],"-f2")) {
	i++;
        if((file_exists(argv[i]))) {
          strcpy(inpfile2,argv[i]);
          input2=open_file(inpfile2,"rb");
        } else {
          printf("specified input file does not exist!\n");
          exit(0);
        }
      } else if (strings_equal(argv[i],"-gm")) {
	i++;
	gm=1;
	strcpy(gminfofile,argv[i]);
      } else if (strings_equal(argv[i],"-block")) {
	i++;
	bsize=atoi(argv[i]);
      } else if (strings_equal(argv[i],"-nch1")) {
	i++;
	nchans1=atoi(argv[i]);
      } else if (strings_equal(argv[i],"-nch2")) {
	i++;
	nchans2=atoi(argv[i]);
      } else if (strings_equal(argv[i],"-startch")) {
	i++;
	startch=atoi(argv[i]);
      } else if (strings_equal(argv[i],"-outnch")) {
	i++;
	outnch=atoi(argv[i]);
      } else if (strings_equal(argv[i],"-ts1")) {
	i++;
	tsamp1=atof(argv[i]);
      } else if (strings_equal(argv[i],"-ts2")) {
	i++;
	tsamp2=atof(argv[i]);
      } else if (strings_equal(argv[i],"-headerless")) {
	headerless=1;
      } else if (strings_equal(argv[i],"-scale")) {
	i++;
	anorm=atof(argv[i]);
      } else if (help_required(argv[1])) {
	pcbeam_gmrt_help();
	exit(0);
      } else {
	pcbeam_gmrt_help();
	sprintf(string,"unknown argument (%s) passed to pcbeam_gmrt",argv[i]);
	error_message(string);
      }
      i++;
    }
  } else {
      pcbeam_gmrt_help();
      exit(0);
  }

  if ((strings_equal(outfile,"dummy")) || (strings_equal(inpfile1,"dummy")) || (strings_equal(inpfile2,"dummy")) || (nchans1==0) || (nchans2==0) || (tsamp1==0) || (tsamp2==0) ) {
    pcbeam_gmrt_help();
    exit(0);
  }
  printf("\n");
//=====================================

  nbits = 16;
  if (obits == 0) obits=nbits;

  if (naddt <= 1) naddt=4096;
  if (obits == 0) obits=nbits;

  // deduce factors by which down-sampling in time/freq is needed.
  // and some sanity checks...
  if (nchans1 == nchans2){
     nchans = nchans1;
     atemp = 1.0;
     f1 = 1;
     f2 = 1;
  } else if (nchans1 > nchans2){
     nchans = nchans2;
     atemp = (float)nchans1/(float)nchans2;
     f1 = (int)atemp;
     f2 = 1;
  } else {
     nchans = nchans1;
     atemp = (float)nchans2/(float)nchans1;
     f2 = (int)atemp;
     f1 = 1;
  }
  if (atemp-(int)atemp != 0.0){
	  printf("No. of chans in two files are not integer multiple of each other!");
	  exit(0);
  }
  if (outnch<=0) outnch = nchans;
  if ( (startch+outnch-1)>=nchans ) {
	  printf("Requested no. of output chans not possible!");
	  exit(0);
  }

  if (tsamp1 == tsamp2){
     tsamp = tsamp1;
     atemp = 1.0;
     t1 = 1;
     t2 = 1;
  } else if (tsamp1 > tsamp2){
     tsamp = tsamp1;
     atemp = (float)tsamp1/(float)tsamp2;
     t2 = (int)atemp;
     t1 = 1;
  } else {
     tsamp = tsamp2;
     atemp = (float)tsamp2/(float)tsamp1;
     t1 = (int)atemp;
     t2 = 1;
  }
  if (atemp-(int)atemp != 0.0){
	  printf("Sampling times in two files are not integer multiple of each other!");
	  exit(0);
  }
  bsize = nchans*4096*nbits;

  // read the timestamp files and get the starting offset
  strcpy(gmhdrfile, inpfile1);
  strcat(gmhdrfile,".hdr");
  mjd1 = gm_mjd(gmhdrfile);
  strcpy(gmhdrfile, inpfile2);
  strcat(gmhdrfile,".hdr");
  mjd2 = gm_mjd(gmhdrfile);
  if (mjd2==mjd1){
	  dtemp = 0.0;
	  itemp = 0;
	  off1 = 0;
	  off2 = 0;
          strcpy(string, inpfile1);
          printf("Starting timestamps are same in both the files.\n");
  }
  else if (mjd2>mjd1){
	  dtemp = (mjd2-mjd1)*86400.0*1000.0/tsamp1;
	  dtemp = round(dtemp*100.0)/100.0;
	  itemp = (int)dtemp;
	  off1 = itemp*nchans1*nbits/8; // offset in bytes
	  off2 = 0;
          printf("File1 has %Lf extra samples (%f seconds) in the beginning.\n",dtemp,(float)(dtemp*tsamp1/1000.0));
          strcpy(string, inpfile2);
  } else {
	  dtemp = (mjd1-mjd2)*86400.0*1000.0/tsamp2;
	  dtemp = round(dtemp*100.0)/100.0;
	  itemp = (int)dtemp;
	  off2 = itemp*nchans2*nbits/8; // offset in bytes
	  off1 = 0;
          printf("File2 has %Lf extra samples (%f seconds) in the beginning.\n",dtemp,(float)(dtemp*tsamp2/1000.0));
          strcpy(string, inpfile1);
  }
  if (dtemp-itemp != 0.0){
	  printf("No. of offset samples are not integer!");
	  exit(0);
  }
  strcat(string,".hdr");
  hdr1 = fopen(string, "r");
  strcpy(gmhdrfile, outfile);
  strcat(gmhdrfile,".hdr");
  hdr2 = fopen(gmhdrfile, "w");
  // copy contents from hdr1 to hdr2
  ctt = fgetc(hdr1);
  while (ctt != EOF)
  {
      fputc(ctt, hdr2);
      ctt = fgetc(hdr1);
  }
  fclose(hdr1);
  fclose(hdr2);
  // also write the info-header
  strcpy(infoFile,outfile);
  strcat(infoFile,".info");
  infohdr = fopen(infoFile, "w");
  fprintf(infohdr,"%Lf\n",tsamp);
  fprintf(infohdr,"freq-ch1\n");
  fprintf(infohdr,"band-width\n");
  fprintf(infohdr,"%d\n",outnch);
  fprintf(infohdr,"src-name\n");
  fclose(infohdr);


  // Now the down-sampling in time/freq., if needed, and actual subtraction of data
  nsblk1=bsize*f1*t1/nbits;
  fblock1=(float *) malloc(nsblk1*sizeof(float));
  nsblk2=bsize*f2*t2/nbits;
  fblock2=(float *) malloc(nsblk2*sizeof(float));
  nsblk=bsize/nbits;
  tempblock=(float *) malloc((nsblk+5)*sizeof(float));
  nsblk=(bsize*outnch/nchans)/nbits;
  fblock=(float *) malloc(nsblk*sizeof(float));
  sblock=(unsigned short *) malloc(nsblk*sizeof(unsigned short));
  cblock=(unsigned char *) malloc(nsblk*sizeof(unsigned short));
  min=0.0;
  max=(float) pow(2.0,(double)obits) -1.0; 

  dur1=dur2=0.0;

  fseek(input1, off1, SEEK_SET); 
  fseek(input2, off2, SEEK_SET); 
  while ((ns1=read_block(input1,nbits,fblock1,nsblk1))>0 && (ns2=read_block(input2,nbits,fblock2,nsblk2))>0){
      ns = nsblk;
      dt1 = ns1*tsamp1/(nchans1*1000.0);
      dt2 = ns2*tsamp2/(nchans2*1000.0);
      // do the required averaging in data from first file
      //printf ("ns1, ns2: %d %d\n",ns1,ns2);
      if (f1>1){  // simple averaging of consecutive points wil do
	      nout = ns1/f1;
	      for (i=0; i<nout; i++){
		      atemp = 0.0;
		      itemp = i*f1;
		      for (k=0; k<f1; k++) atemp = atemp + fblock1[itemp+k];
		      fblock1[i] = atemp/(float)f1;
	      }
	      ns1 = nout;
      }
      if (t1>1){  // downsampling in time
	      nout = ns1/t1;
              kk = (ns1/nchans)/t1;
      	      for (j=0;j<kk;j++) {
                      jt1 = j*t1;
                      jt2 = j*nchans;
                      for (i=0; i<nchans; i++){
                         atemp = 0.0;
                         for (jj=0;jj<t1;jj++){
                           k = (jt1+jj)*nchans + i;
                           atemp = atemp + fblock1[k];
                         }
                         k = jt2 + i;
                         fblock1[k] = atemp/(float)t1;
                      }
              }
	      ns1 = nout;
      }
      // do the required averaging in data from second file
      if (f2>1){  // simple averaging of consecutive points wil do
	      nout = ns2/f2;
	      for (i=0; i<nout; i++){
		      atemp = 0.0;
		      itemp = i*f2;
		      for (k=0; k<f2; k++) atemp = atemp + fblock2[itemp+k];
		      fblock2[i] = atemp/(float)f2;
	      }
	      ns2 = nout;
      }
      if (t2>1){  // downsampling in time
	      nout = ns2/t2;
              kk = (ns2/nchans)/t2;
      	      for (j=0;j<kk;j++) {
                      jt1 = j*t2;
                      jt2 = j*nchans;
                      for (i=0; i<nchans; i++){
                         atemp = 0.0;
                         for (jj=0;jj<t2;jj++){
                           k = (jt1+jj)*nchans + i;
                           atemp = atemp + fblock2[k];
                         }
                         k = jt2 + i;
                         fblock2[k] = atemp/(float)t2;
                      }
              }
	      ns2 = nout;
      }
      //printf ("ns1, ns2: %d %d\n",ns1,ns2);
      if (ns2 != ns1){
	      printf ("\nUnequal durations of the two files? ");
	      printf ("(Potentially reached EOF for one of the files) \n");
	      printf ("Total durations (s) read, processed and unloaded: %f %f\n\n",dur1,dur2);
	      ns = 0;
	      dt1=dt2=0.0;
	      //exit(0);
	      break;
      } else {
	      ns = ns1;
      }
      dur1 = dur1 + dt1;
      dur2 = dur2 + dt2;
      //printf ("Durations processed: %f   %f   %d  %d\n",dur1,dur2,ns,nsblk);


      if(anorm==-9999.0){ // first block and uder didn't specify scaling
                        // check if scaling is needed
        for (i=0; i<ns; i++) {
		if (fblock2[i] != 0.0) tempblock[i] = fblock1[i]/fblock2[i];
		else tempblock[i] = 0.0;
	}
        medrms(tempblock,ns);
        med_fact = tempblock[ns];
        std = tempblock[ns+1];
        printf ("Median-scaling-factor, rms:  %f   %f  \n",med_fact,std);
	anorm = round(med_fact*10.0)/10.0;
        printf ("Rounded-off Scaling-factor, rms:  %f   %f  \n",anorm,std);
        printf ("Assumed scaling factor (input-1/input-2): %f\n",anorm);
      }
      for (i=0; i<ns; i++) tempblock[i] = fblock1[i] - anorm*fblock2[i];
      if (startch==0 && outnch==nchans){
	      for (i=0; i<ns; i++) fblock[i] = tempblock[i];
      } else {
              kk = ns/nchans; // number of samples
	      k = 0;
      	      for (j=0;j<kk;j++) {
                      jt1 = j*nchans;
                      for (i=startch; i<startch+outnch; i++){
			      fblock[k] = tempblock[jt1+i];
			      k = k+1;
                      }
              }
	      ns = k;

      }
      nout=ns;
      switch (obits) {
      case 32:
        fwrite(fblock,sizeof(float),nout,output);
        break;
      case 16:
        float2short(fblock,nout,min,max,sblock);
        fwrite(sblock,sizeof(unsigned short),nout,output);
        break;
      case 8:
        float2char(fblock,nout,min,max,cblock);
        fwrite(cblock,sizeof(unsigned char),nout,output);
        break;
      } 
  } 
  printf("\nSampling time (ms) in the PC-beam data: %Lf \n",tsamp);
  printf("No. of Freq-channels in the PC-beam data (out of %d): %d \n\n",nchans,outnch);

  free (tempblock);
  free (fblock);
  free (fblock1);
  free (sblock);
  free (cblock);
  free (fblock2);
  fclose(input1);
  fclose(input2);
  fclose(output);
  exit(0);
}
