/*Authors: Klaus Seywald (klaus.seywald@mytum.de) 
         and Simon Binder (simon.binder@tum.de)

This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
*/
#include "mex.h"
#include <math.h>
/*
 * xtimesy.c - example found in API guide
 *
 * multiplies an input scalar times an input matrix and outputs a
 * matrix
 *
 * This is a MEX-file for MATLAB.
 * Copyright 1984-2011 The MathWorks, Inc.
 */

/* $Revision: 1.10.6.4 $ */

void compute_influence_vring(double *grid,double *panels, double *colloc,double *w1,double *w2,double *w3,mwSize n_colloc,mwSize n_panels)
{
    mwSize i,j,count=0;
    
    double wAB_1,wAB_2,wAB_3;
    double r1xr0_1,r1xr0_2,r1xr0_3;
    double norm_r0;
    double r01,r02,r03,r11,r12,r13,r21,r22,r23,r1xr2_1,r1xr2_2,r1xr2_3,r1dr2_1,r1dr2_2,r1dr2_3,dotp_r0_r1dr2,norm_r1xr2,normr1,normr2;
    double coredist;
    double lambda,xi,sp1,sp2,sp3;
    
    double p1,p2,p3,a1,a2,a3,n1,n2;
    double vortex[12];
    int k;
    
    int doomed;
    doomed=0;
    for (k=0;k<n_panels;k++){
	//for (k=0;k<10;k++){
        
        vortex[0]=*(grid+((int)*(panels+4*k)-1)*3);
        vortex[1]=*(grid+((int)*(panels+4*k)-1)*3+1);
        vortex[2]=*(grid+((int)*(panels+4*k)-1)*3+2);
        
        vortex[3]=*(grid+((int)*(panels+4*k+1)-1)*3);
        vortex[4]=*(grid+((int)*(panels+4*k+1)-1)*3+1);
        vortex[5]=*(grid+((int)*(panels+4*k+1)-1)*3+2);
		
		vortex[6]=*(grid+((int)*(panels+4*k+2)-1)*3);
        vortex[7]=*(grid+((int)*(panels+4*k+2)-1)*3+1);
        vortex[8]=*(grid+((int)*(panels+4*k+2)-1)*3+2);
        
        vortex[9]=*(grid+((int)*(panels+4*k+3)-1)*3);
        vortex[10]=*(grid+((int)*(panels+4*k+3)-1)*3+1);
        vortex[11]=*(grid+((int)*(panels+4*k+3)-1)*3+2);
        
      // mexPrintf("x   %f  %f  %f %f  \n",vortex[0],vortex[4],vortex[7],vortex[10]);
      // mexPrintf("y   %f  %f  %f %f  \n",vortex[1],vortex[5],vortex[8],vortex[11]);
      // mexPrintf("z   %f  %f  %f %f  \n",vortex[3],vortex[6],vortex[9],vortex[12]);
      // mexPrintf("hallooo   %d, %d \n",(int)*(te_idx+(((int)*(panels+4*k+3))-1)-1)-1,(int)*(te_idx+(((int)*(panels+4*k+2))-1)-1)-1);
      // mexPrintf("hallooo2   %d, %d \n",((int)*(panels+4*k+3)-1),(((int)*(panels+4*k+2)-1)-1)-1);
        for (j=0;j<n_colloc;j++){
		//for (j=0;j<10;j++){
            wAB_1=0;
            wAB_2=0;
            wAB_3=0;
			doomed=0;
            
            p1=*(colloc+j*3);
            p2=*(colloc+1+j*3);
            p3=*(colloc+2+j*3);
                    
            for (i=0;i<4;i++){
				r11=*(colloc+j*3)-*(vortex+i*3);
				r12=*(colloc+1+j*3)-*(vortex+i*3+1);
				r13=*(colloc+2+j*3)-*(vortex+i*3+2);
                
                a1=*(vortex+i*3);
                a2=*(vortex+i*3+1);
                a3=*(vortex+i*3+2);

				if (i==3){
					r21=*(colloc+j*3)-*(vortex+(0)*3);
					r22=*(colloc+1+j*3)-*(vortex+(0)*3+1);
					r23=*(colloc+2+j*3)-*(vortex+(0)*3+2);
                
					r01=*(vortex+0)-*(vortex+i*3);
					r02=*(vortex+1)-*(vortex+i*3+1);
					r03=*(vortex+2)-*(vortex+i*3+2);
				}
				else{
					r21=*(colloc+j*3)-*(vortex+(i+1)*3);
					r22=*(colloc+1+j*3)-*(vortex+(i+1)*3+1);
					r23=*(colloc+2+j*3)-*(vortex+(i+1)*3+2);
                
					r01=*(vortex+(i+1)*3)-*(vortex+i*3);
					r02=*(vortex+(i+1)*3+1)-*(vortex+i*3+1);
					r03=*(vortex+(i+1)*3+2)-*(vortex+i*3+2);
				}
                
                n1=r01*r03*r13-r03*r03*r11-r02*r02*r11+r01*r02*r12;
                n2=-r03*r03*r12+r02*r03*r13+r01*r02*r11-r01*r01*r12;
                
                r1xr2_1=r12*r23-r13*r22;
                r1xr2_2=r13*r21-r11*r23;
                r1xr2_3=r11*r22-r12*r21;
                
                norm_r1xr2=sqrt(r1xr2_1*r1xr2_1+r1xr2_2*r1xr2_2+r1xr2_3*r1xr2_3);
                norm_r0=sqrt(r01*r01+r02*r02+r03*r03);
                
                normr1=sqrt(r11*r11+r12*r12+r13*r13);
                normr2=sqrt(r21*r21+r22*r22+r23*r23);
                
                r1xr0_1=r12*r03-r13*r02;
                r1xr0_2=r13*r01-r11*r03;
                r1xr0_3=r11*r02-r12*r01;
                
                coredist=sqrt(r1xr0_1*r1xr0_1+r1xr0_2*r1xr0_2+r1xr0_3*r1xr0_3)/norm_r0;
                
                if (coredist<0.19){
                    // improve for vertical elements(r03!!)
                    if (abs(r01)>0.1){
                        lambda=(r11)/r01;
                       // lambda=(p1/r01-a1/r01+a2/n2*n1-p3/n2*n1)/(1-r02/n2*n1);
                    }
                    else{
                        //lambda=(p3/r02-a2/r02+a1/n1*n2-p1/n1*n2)/(1-r01/n1*n2);
                        lambda=r12/r02;
                    }
                    
                    sp1=*(vortex+i*3)+r01*lambda;
                    sp2=*(vortex+i*3+1)+r02*lambda;
                    sp3=*(vortex+i*3+2)+r03*lambda;
                    
                    if (abs(r01)>0.1){
                        lambda=(sp1-*(vortex+i*3))/r01;
                    }
                    else{
                        lambda=(sp2-*(vortex+i*3+1))/r02;
                    }
                    
                    if (i==3){
                        if (abs(r01)>0.1){
                            xi=(sp1-*(vortex+(0)*3))/r01;
                        }
                        else{
                            xi=(sp2-*(vortex+(0)*3+1))/r02;
                        }
                    }
                    else
                    {
                        if (abs(r01)>0.1){
                            xi=(sp1-*(vortex+(i+1)*3))/r01;
                        }
                        else{
                            xi=(sp2-*(vortex+(i+1)*3+1))/r02;
                        }
                    }
                    
                    if ((xi<0.19) && (lambda>-0.19)){
                        doomed=1;
                    }
                    if ((xi>-0.19) && (lambda<0.19)){
                        doomed=1;
                    }
                }

                if  ( (!((normr1<=1E-7)||(normr2<=1E-7))) && coredist>0.01)
                {
                    r1dr2_1=r11/normr1-r21/normr2;
                    r1dr2_2=r12/normr1-r22/normr2;
                    r1dr2_3=r13/normr1-r23/normr2;
                    
                    dotp_r0_r1dr2=r1dr2_1*r01+r1dr2_2*r02+r1dr2_3*r03;
                    
                    wAB_1=wAB_1-r1xr2_1/(norm_r1xr2*norm_r1xr2)*dotp_r0_r1dr2;
                    wAB_2=wAB_2-r1xr2_2/(norm_r1xr2*norm_r1xr2)*dotp_r0_r1dr2;
                    wAB_3=wAB_3-r1xr2_3/(norm_r1xr2*norm_r1xr2)*dotp_r0_r1dr2;
                }
                else
                {
                mexPrintf("no influence i j k %d %d %d\n",i,j,k);
                }
            }
            if (doomed==1){
             //wAB_1=0;
             //wAB_2=0;
             //wAB_3=0;         
            }
           //  *(w+n_colloc*k+j)=*(colloc_nvec+j*3)*wAB_1+*(colloc_nvec+1+j*3)*wAB_2+*(colloc_nvec+2+j*3)*wAB_3;
            //if ((wAB_1>-20)&&(wAB_1<20)){
                *(w1+n_colloc*k+j)=wAB_1;
           // }
           // else{
           //     *(w1+n_colloc*k+j)=0;
           // } 
         //  if ((wAB_2>-20)&&(wAB_2<20)){
           *(w2+n_colloc*k+j)=wAB_2;
         //  }
         //  else{
         //       *(w2+n_colloc*k+j)=0;
         //  }
         //  if ((wAB_3>-20)&&(wAB_3<20)){  
           *(w3+n_colloc*k+j)=wAB_3;
         //  }
         //  else
         //  {
          //    *(w3+n_colloc*k+j)=0;  
          // }
           // mexPrintf("   %f \n",*(w+n_colloc*k+j)/(4*3.141592));
           // *(wind+n_colloc*k+j)=*(colloc_nvec+j*3)*wAB_1+*(colloc_nvec+1+j*3)*wAB_2+*(colloc_nvec+2+j*3)*wAB_3;
        }
    }
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    double *grid,*colloc,*colloc_nvec,*w1,*w2,*w3,*wind,*Uinf;
    double *panels,*te_idx;
    size_t n_colloc;
    size_t n_panels;
    
    /*  check for proper number of arguments */
    /* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks you out of
     the MEX-file) */
    if(nrhs!=6)
        mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumInputs",
                "Two inputs required.");
    if(nlhs!=3)
        mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumOutputs",
                "One output required.");
    
    /* check to make sure the first input argument is a scalar */
    //if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
    //    mxGetN(prhs[0])*mxGetM(prhs[0])!=1 ) {
    //  mexErrMsgIdAndTxt( "MATLAB:xtimesy:xNotScalar",
    //          "Input x must be a scalar.");
    //}
    
    /*  get the scalar input vortex */
    grid =(double *)mxGetPr(prhs[0]);
   
    panels =(double *)mxGetPr(prhs[1]);
    n_panels= mxGetN(prhs[1]);
    te_idx =(double *)mxGetPr(prhs[2]);
    /*  create a pointer to the input matrix colloc */
    
    colloc =(double *)mxGetPr(prhs[3]);
    n_colloc= mxGetN(prhs[3]);
    
    
    colloc_nvec =(double *)mxGetPr(prhs[4]);
    Uinf =(double *)mxGetPr(prhs[5]);
    /*  get the dimensions of the matrix input y */
    /* mrows = mxGetM(prhs[0]);
 /* ncols = mxGetN(prhs[0]);
  
  /*  set the output pointer to the output matrix */
    plhs[0]=mxCreateDoubleMatrix(n_colloc,n_panels,mxREAL);
    plhs[1]=mxCreateDoubleMatrix(n_colloc,n_panels,mxREAL);
    plhs[2]=mxCreateDoubleMatrix(n_colloc,n_panels,mxREAL);
     //plhs[1]=mxCreateDoubleMatrix(n_colloc,n_colloc,mxREAL);
    /*  create a C pointer to a copy of the output matrix */
    w1 = mxGetPr(plhs[0]);
    w2 = mxGetPr(plhs[1]);
    w3 = mxGetPr(plhs[2]);
    //wind = mxGetPr(plhs[1]);
    //mexPrintf("bla %i \n ",n_colloc);
    // mexPrintf("Vortex %f,%f,%f \n %f,%f,%f \n %f,%f,%f \n ", *(vortex),*(vortex+1),*(vortex+2),*(vortex+3),*(vortex+4),*(vortex+5),*(vortex+6),*(vortex+7),*(vortex+8));
    /*  call the C subroutine */
    compute_influence_vring(grid,panels,colloc,w1,w2,w3,n_colloc,n_panels);
}
