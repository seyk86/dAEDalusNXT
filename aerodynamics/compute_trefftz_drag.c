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

void compute_trefftz_drag(double* grid,double* panels,double* te_idx,double* Gamma,double* Uinf,double* n_pan,double *CD,int n,double * coords_y,double * coords_z,double *vel_y, double *vel_z)
{
    
    mwSize i,j,k,count=0;
    const double pi=3.1415926535897932384626433832795;
    double dist,u,v,w;
    double te_x,te_y,te_z;
    double lambda;
    double  *trefftz_coords_in,*trefftz_coords_out;
    double max_y,min_y,max_z,min_z;
    double start_y,start_z;
         double trefftz_width,trefftz_height;   
    double dl_y,dl_z,lr;
    double v1_y,v2_y,v1_z,v2_z;
    double r_y,r_z;
    
    double core_radius;
        
    mxArray *zw1,  *zw2;
    zw1=mxCreateDoubleMatrix(n,3,mxREAL);
    zw2=mxCreateDoubleMatrix(n,3,mxREAL);
    
    trefftz_coords_in=(double *)mxGetPr(zw1);
    trefftz_coords_out=(double *)mxGetPr(zw2);
    
    dist=-200;
    
    u=Uinf[0];
    v=Uinf[1];
    w=Uinf[2];
    
    core_radius=0.01;
    
  for (i=0;i<n;i++){
        
        te_x=*(grid+((int)*(te_idx+((int)*(panels+4*i+3)-1)-1)-1)*3);
        te_y=*(grid+((int)*(te_idx+((int)*(panels+4*i+3)-1)-1)-1)*3+1);
        te_z=*(grid+((int)*(te_idx+((int)*(panels+4*i+3)-1)-1)-1)*3+2);
        
        lambda=(dist*u-(te_x*u+te_y*v+te_z*w))/(u*u+v*v+w*w);
        
        *(trefftz_coords_in+i)=te_x+lambda*u;
        *(trefftz_coords_in+1*(n)+i)=te_y+lambda*v;
        *(trefftz_coords_in+2*(n)+i)=te_z+lambda*w;
        
        te_x=*(grid+((int)*(te_idx+((int)*(panels+4*i+2)-1)-1)-1)*3);
        te_y=*(grid+((int)*(te_idx+((int)*(panels+4*i+2)-1)-1)-1)*3+1);
        te_z=*(grid+((int)*(te_idx+((int)*(panels+4*i+2)-1)-1)-1)*3+2);
        
        lambda=(dist*u-(te_x*u+te_y*v+te_z*w))/(u*u+v*v+w*w);
        
        *(trefftz_coords_out+i)=te_x+lambda*u;
        *(trefftz_coords_out+1*(n)+i)=te_y+lambda*v;
        *(trefftz_coords_out+2*(n)+i)=te_z+lambda*w;
    }
    
    min_y=*(trefftz_coords_in+(n)*1);
    min_z=*(trefftz_coords_in+(n)*2);
    max_y=*(trefftz_coords_in+(n)*1);
    max_z=*(trefftz_coords_in+(n)*2);
    
    

    
    for(i=1;i<n;i++)
    {
        if(*(trefftz_coords_in+(n)*1+i)<min_y){
            min_y=*(trefftz_coords_in+(n)*1+i);
        }
        if(*(trefftz_coords_in+(n)*2+i)<min_z){
            min_z=*(trefftz_coords_in+(n)*2+i);
        }
        if(*(trefftz_coords_in+(n)*1+i)>max_y){
            max_y=*(trefftz_coords_in+(n)*1+i);
        }
        if(*(trefftz_coords_in+(n)*2+i)>max_z){
            max_z=*(trefftz_coords_in+(n)*2+i);
        }
    }
    
    trefftz_width=3.0*abs(max_y-min_y);
    trefftz_height=3.0*abs(max_z-min_z);
    
    if( trefftz_height<trefftz_width){
           trefftz_height= trefftz_width;
    }

    start_y=min_y-(trefftz_width/2-(max_y-min_y)/2);
    start_z=min_z-(trefftz_height/2-(max_z-min_z)/2);

    dl_y=trefftz_width/ *n_pan;
    dl_z=trefftz_height/ *n_pan;
    
    *CD=0.0000000;

    for (i=0;i<(int)*n_pan;i++)
    {
        for(j=0;j<(int)*n_pan;j++)
        {
            v1_y=0.000000000000;
            v1_z=0.000000000000;
            v2_y=0.000000000000;
            v2_z=0.000000000000;
            
            for(k=0;k<n;k++)
            {
                r_y=(start_y+dl_y/2)+i*dl_y-*(trefftz_coords_in+(n)*1+k);
                r_z=(start_z+dl_z/2)+j*dl_z-*(trefftz_coords_in+(n)*2+k);
                
                lr=sqrt(r_y*r_y+r_z*r_z);
                
                if(lr>core_radius){
                    v1_y=v1_y-*(Gamma+k)/(2*pi*lr)*(r_z/lr);
                    v1_z=v1_z-*(Gamma+k)/(2*pi*lr)*(-r_y/lr);
                }
                else{
                    v1_y=v1_y-*(Gamma+k)/(2*pi*core_radius)*(r_z/core_radius);
                    v1_z=v1_z-*(Gamma+k)/(2*pi*core_radius)*(-r_y/core_radius);  
                }
                
                r_y=(start_y+dl_y/2)+i*dl_y-*(trefftz_coords_out+(n)*1+k);
                r_z=(start_z+dl_z/2)+j*dl_z-*(trefftz_coords_out+(n)*2+k);
                
                lr=sqrt(r_y*r_y+r_z*r_z);
                
                if(lr>core_radius){
                    v2_y=v2_y+*(Gamma+k)/(2*pi*lr)*(r_z/lr);
                    v2_z=v2_z+*(Gamma+k)/(2*pi*lr)*(-r_y/lr);
                }
                else{
                    v2_y=v2_y+*(Gamma+k)/(2*pi*core_radius)*(r_z/core_radius);
                    v2_z=v2_z+*(Gamma+k)/(2*pi*core_radius)*(-r_y/core_radius); 
                    
                }
            }
            *(vel_y+i*((int)*n_pan)+j)=v1_y+v2_y;
            *(vel_z+i*((int)*n_pan)+j)=v1_z+v2_z;
            *(coords_y+i*((int)*n_pan)+j)=(start_y+dl_y/2)+i*dl_y;
            *(coords_z+i*((int)*n_pan)+j)=(start_z+dl_z/2)+j*dl_z;
            
            *CD=*CD+(*(vel_y+i*((int)*n_pan)+j)+*(vel_z+i*((int)*n_pan)+j))*(*(vel_y+i*((int)*n_pan)+j)+*(vel_z+i*((int)*n_pan)+j))*dl_y*dl_z;
        }
    } 
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    double *grid,*panels,*te_idx,*Gamma,*Uinf;
    double *CD,*n_pan,*treffz_width,*treffz_height;
    
    double *coords_y,*coords_z,*vel_y,*vel_z;
        size_t n;
    size_t n_colloc;
    
    /*  check for proper number of arguments */
    /* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
     * within an if statement, because it will never get to the else
     * statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks you out of
     * the MEX-file) */
    if(nrhs!=6)
        mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumInputs",
                "Eight inputs required.");
    if(nlhs!=5)
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
    te_idx =(double *)mxGetPr(prhs[2]);
    Gamma  =(double *)mxGetPr(prhs[3]);
    Uinf =(double *)mxGetPr(prhs[4]);
    n_pan =(double *)mxGetPr(prhs[5]);
    
    /*  get the dimensions of the matrix input y */

    n= mxGetM(prhs[3]);
    
    /*  set the output pointer to the output matrix */
    plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
    plhs[1]=mxCreateDoubleMatrix((int)(*n_pan),(int)(*n_pan),mxREAL);
    plhs[2]=mxCreateDoubleMatrix((int)(*n_pan),(int)(*n_pan),mxREAL);
    plhs[3]=mxCreateDoubleMatrix((int)(*n_pan),(int)(*n_pan),mxREAL);
    plhs[4]=mxCreateDoubleMatrix((int)(*n_pan),(int)(*n_pan),mxREAL);
    
    CD=mxGetPr(plhs[0]);
    coords_y = mxGetPr(plhs[1]);
    coords_z = mxGetPr(plhs[2]);
    vel_y = mxGetPr(plhs[3]);
    vel_z = mxGetPr(plhs[4]);
    //mexPrintf("bla %i \n ",n_colloc);
    // mexPrintf("Vortex %f,%f,%f \n %f,%f,%f \n %f,%f,%f \n ", *(vortex),*(vortex+1),*(vortex+2),*(vortex+3),*(vortex+4),*(vortex+5),*(vortex+6),*(vortex+7),*(vortex+8));
    /*  call the C subroutine */
    compute_trefftz_drag(grid,panels,te_idx,Gamma,Uinf,n_pan,CD,n,coords_y,coords_z,vel_y,vel_z);
}
