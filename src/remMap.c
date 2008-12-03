#include <stdlib.h>
#include <math.h>
#include <stdio.h>

/////////   c code for joint modeling CGH and expression data using group lasso  ////////////
// Regression model;
// using original shooting methods;

/////////  accessary functions
double SoftShrink(double y, double lam);
void Assign(int P, int Q, double *data_old, double *data_copy);
void CalBnorm(int P, int Q, double *Beta, int *C_m, double *Bnorm);
void Update(int cur_p, int N, int P, int Q, double lambda1, double lambda2, int *C_m, double *X_m, double *Xnorm, double *E, double *Beta, double *Bnorm, double *phi_old, double *phi);
double Dist(int P, int Q, double *phi_last, double * phi);

////////   main function
void MultiRegGroupLassoDegree (int *NN, int *PP, int *QQ, double *X_m, double *Y_m, int *C_m, int *NN1, double *lam1, int *NN2,
 double *lam2, double *Degree, double *Debug);

void MultiRegGroupLassoIni (int *NN, int *PP, int *QQ, double *X_m, double *Y_m, int *C_m, double *lam1,  double *lam2,
                 double *Phi_initial, double *Phi_output, double *Phi_debug, int *N_iter, double *RSS, double *E_debug);


//////////////////////////////////////////////////////
////////////////////////// no initial value

void MultiRegGroupLasso (int *NN, int *PP, int *QQ, double *X_m, double *Y_m, int *C_m, double *lam1,  double *lam2,
                         double *Phi_output, double *Phi_debug, int *N_iter, double *RSS, double *E_debug)
{
	/// Variables:
	/// N: Number of observations;
	/// P: Number of DNA;
	/// Q: Number of mRNA;

	/// X_m: matrix of N by P, input DNA data; each column has mean 0.
	/// Y_m: matrix of N by Q, input RNA data; each column has mean 0.
	/// C_m: matrix of P by Q, input data; indicator for "self relationship";
	///      C_m[i,j]=1 means the corresponding beta[i,j] will be penalized;	
      ///      C_m[i,j]=2 means the corresponding beta[i,j] will not be penalized.
     ///       C_m[i,j]=0 means the corresponding beta[i,j] will not be considered in the model.

	/// lambda1: group panelty;
	/// lambda2: lasso panelty;

	/// Phi_output: matrix of P by Q, estimated coefficients.
	/// Phi_debug: matrix of P by Q, for debug
	/// N_iter:   indicate the total iteration
	/// RSS: the sum of squre error

	int N, P, Q;
	int n,p,q;
	int i,j;
	int n_iter;
	int cur_p;
	int *pick;
	int n_pick;
	int *unpick;
	int n_unpick;

        double rss;
	double temp,temp1, temp2;
	double lambda1, lambda2;
	double *Xnorm;
	double *Beta;
	double *phi;
	double *phi_old, *phi_last;
	double *E;
	double *Bnorm;
	double eps;
	double flag;
	double flag_a;

	lambda1=*lam1;
	lambda2=*lam2;
	N=*NN;
	P=*PP;
	Q=*QQ;

        n_iter=0;
	eps=1e-6;

	pick=(int *) malloc (P*sizeof(int));
	unpick=(int *) malloc (P*sizeof(int));
	Xnorm=(double *) malloc (P*sizeof(double));
	Beta=(double *) malloc (P*Q*sizeof(double));
	phi=(double *) malloc (P*Q*sizeof(double));
	phi_old=(double *) malloc (P*Q*sizeof(double));
	phi_last=(double *) malloc (P*Q*sizeof(double));
	E=(double *) malloc (N*Q*sizeof(double));
	Bnorm=(double *) malloc (P*sizeof(double));

///// initial value
	for(p=0;p<P;p++)
	 {
		 Xnorm[p]=0;
		 for(n=0;n<N;n++)
		   Xnorm[p]=Xnorm[p]+X_m[n*P+p]*X_m[n*P+p];
	 }

////////(1) calculate lasso solution Beta
	for(p=0;p<P;p++)
	 for(q=0;q<Q;q++)
	  {
         
            if(C_m[p*Q+q]==0) /// not considered in the model
               Beta[p*Q+q]=0;
            else
            {
		  temp=0;
		  for(n=0;n<N;n++)
		    temp=temp+X_m[n*P+p]*Y_m[n*Q+q];
		  if(C_m[p*Q+q]==1) /// penalized
				Beta[p*Q+q]=SoftShrink(temp, lambda2)/Xnorm[p];
                   else /// not penalized
                                Beta[p*Q+q]=temp/Xnorm[p];
            }
	  }
 ///////// (2) calculate Bnorm;
    CalBnorm(P, Q, Beta, C_m, Bnorm);

 //////// (3) calculate phi

    for(p=0;p<P;p++)
       for(q=0;q<Q;q++)
         {
             
            if(C_m[p*Q+q]==0) /// not considered in the model
               phi[p*Q+q]=0;
            else if(C_m[p*Q+q]==1 && Bnorm[p]>1e-6) /// penalized			  
                  {
			   temp=Xnorm[p]*Bnorm[p];
			   temp1=SoftShrink(temp, lambda1);
    		           phi[p*Q+q]=temp1*Beta[p*Q+q]/temp;
		       }            else
		phi[p*Q+q]=Beta[p*Q+q]; /// not penalized
         }


////////// (4) Derive active set
     CalBnorm(P, Q, phi, C_m, Bnorm);
///////// (5) Residue
    for(n=0;n<N;n++)
     for(q=0;q<Q;q++)
      {
		  temp=0;
		  for(p=0;p<P;p++)
		    temp=temp+phi[p*Q+q]*X_m[n*P+p];
 		  E[n*Q+q]=Y_m[n*Q+q]-temp;
      }
        /////////////////////////////////////////////////
        ///////////// (6) begin update
        ////////////////////////////////////////////////

             flag=100;
             while(flag>1e-6 && n_iter<1e+10)
             {
			   //printf("iter= %d", n_iter);

			   //////////// derive active set
			   n_pick = 0;
			   n_unpick=0;
			   for(p = 0; p<P;p++)
			    {
			      if( Bnorm[p]>1e-6)
				  {
				 	 pick[n_pick] =p;
				 	 n_pick = n_pick + 1;
		           }
		           else
		           {
					 unpick[n_unpick] =p;
				 	 n_unpick = n_unpick + 1;
				   }
			     }

			   //////////// prepare for active set
			   flag_a=100;
			   //if(n_pick==0)
			   //    flag.a=0;

			   while(flag_a>1e-6 && n_iter<1e+10)/////////////////1) Active set
			    {

				 ///phi_last=phi
				 Assign(P, Q, phi, phi_last);				 
                                 Assign(P, Q, phi, phi_old);                
                               /// update coefficients that is penalized                 
                               for(i=0;i<n_pick;i++)
                                  {	
                        	    cur_p=pick[i];
        		                Update(cur_p, N, P, Q, lambda1, lambda2, C_m, X_m, Xnorm, E, Beta, Bnorm, phi_old, phi);
			                n_iter=n_iter+1;				  
                                  }

				  /// update coefficients that is not penalized
				  for(i=0; i<n_unpick; i++)
				  {
					  p=unpick[i];
					  for(q=0; q<Q; q++)
						   if(C_m[p*Q+q]==2)
						    {
								temp=0;
								for(n=0;n<N;n++)
								  temp=temp+E[n*Q+q]*X_m[n*P+p];
								phi[p*Q+q]=temp/Xnorm[p]+phi_old[p*Q+q];	
                                                ///update residue
								for(n=0; n<N;n++)
								     E[n*Q+q]=E[n*Q+q]+(phi_old[p*Q+q]-phi[p*Q+q])*X_m[n*P+p];
								///update phi_old
								phi_old[p*Q+q]=phi[p*Q+q];                 
                                                                n_iter=n_iter+1;
					         } /// end if
				  }              
                         flag_a=Dist(P, Q, phi_last, phi);
			    }///end of active set while(flag.a>1e-6)
               //////////////// 2) Full loop
               Assign(P, Q, phi, phi_last);
               Assign(P, Q, phi, phi_old);
               for(p=0;p<P;p++)               
                {
		        Update(p, N, P, Q, lambda1, lambda2, C_m, X_m, Xnorm, E, Beta, Bnorm, phi_old, phi);	
                 	n_iter=n_iter+1;	
            	}	
       	flag=Dist(P, Q, phi_last, phi);	
            }//end of all loop while(flag>1e-6)

  //////////// calculate Residue
    rss=0;  
    for(n=0;n<N;n++) 
     for(q=0;q<Q;q++)
	   {
		  rss=rss+E[n*Q+q]*E[n*Q+q];
	   }
    Assign(N, Q, E, E_debug);

//////////// return phi
Assign(P, Q, phi, Phi_output);
Assign(P, Q, phi_last, Phi_debug);
*N_iter=n_iter;
*RSS=rss;

//////// free allocated variables
free(pick);
free(unpick);
free(Xnorm);
free(Beta);
free(phi);
free(phi_old);
free(phi_last);
free(E);
free(Bnorm);
}///end MultiRegGroupLasso function 


///////////////////////////////////////////////////////////////
/////////////////////////// with initial value
void MultiRegGroupLassoIni (int *NN, int *PP, int *QQ, double *X_m, double *Y_m, int *C_m, double *lam1,  double *lam2, double *Phi_initial, double *Phi_output, double *Phi_debug, int *N_iter, double *RSS, double *E_debug)
{
	/// Variables:
	/// N: Number of observations;
	/// P: Number of DNA;
	/// Q: Number of mRNA;
	/// X_m: matrix of N by P, input DNA data; each column has mean 0.
	/// Y_m: matrix of N by Q, input RNA data; each column has mean 0.
	/// C_m: matrix of P by Q, input data; indicator for "self relationship";

	
       ///      C_m[i,j]=1 means the corresponding beta[i,j] will be penalized;
	///      C_m[i,j]=2 means the corresponding beta[i,j] will not be penalized.
        ///       C_m[i,j]=0 means the corresponding beta[i,j] will not be considered in the model.	
       /// lambda1: group panelty;
	/// lambda2: lasso panelty;
        /// Phi_inital: inital value;
	/// Phi_output: matrix of P by Q, estimated coefficients.
	/// Phi_debug: matrix of P by Q, for debug
	/// N_iter:   indicate the total iteration
	/// RSS: the sum of squre error
	int N, P, Q;
	int n,p,q;
	int i,j;
	int n_iter;
	int cur_p;
	int *pick;
	int n_pick;
	int *unpick;
	int n_unpick;
        double rss;
	double temp,temp1, temp2;
	double lambda1, lambda2;
	double *Xnorm;
	double *Beta;
	double *phi;
	double *phi_old, *phi_last;
	double *E;
	double *Bnorm;
	double eps;
	double flag;
	double flag_a;
	lambda1=*lam1;
	lambda2=*lam2;
	N=*NN;
	P=*PP;
	Q=*QQ;        
      n_iter=0;
	eps=1e-6;
	pick=(int *) malloc (P*sizeof(int));
	unpick=(int *) malloc (P*sizeof(int));
	Xnorm=(double *) malloc (P*sizeof(double));
	Beta=(double *) malloc (P*Q*sizeof(double));
	phi=(double *) malloc (P*Q*sizeof(double));
	phi_old=(double *) malloc (P*Q*sizeof(double));
	phi_last=(double *) malloc (P*Q*sizeof(double));
	E=(double *) malloc (N*Q*sizeof(double));
	Bnorm=(double *) malloc (P*sizeof(double));

///// initial value
	for(p=0;p<P;p++)
	 {
		 Xnorm[p]=0;
		 for(n=0;n<N;n++)
		   Xnorm[p]=Xnorm[p]+X_m[n*P+p]*X_m[n*P+p];
	 }

 //////// (1) initial phi

    for(p=0;p<P;p++)
       for(q=0;q<Q;q++)
         {
   if(C_m[p*Q+q]==0)
         phi[p*Q+q]=0;
    else
	   phi[p*Q+q]=Phi_initial[p*Q+q];
         }
///////// (3) initial residue
    for(n=0;n<N;n++)
     for(q=0;q<Q;q++)
      {
		  temp=0;
		  for(p=0;p<P;p++)
		    temp=temp+phi[p*Q+q]*X_m[n*P+p];
 		  E[n*Q+q]=Y_m[n*Q+q]-temp;
      }

        /////////////////////////////////////////////////
        ///////////// (4) begin update
        ////////////////////////////////////////////////

       flag=100;
       while(flag>1e-6 && n_iter<1e+10)
       {
			   //////////// derive active set
			   n_pick = 0;
			   n_unpick=0;
			   for(p = 0; p<P;p++)
			    {
			      if( Bnorm[p]>1e-6)
				  {
				 	 pick[n_pick] =p;
				 	 n_pick = n_pick + 1;
		                   }
		              else
		                  {
					 unpick[n_unpick] =p;
				 	 n_unpick = n_unpick + 1;
				   }
			     }

			   //////////// prepare for active set
			   flag_a=100;
			   //if(n_pick==0)
			   //    flag.a=0;

			   while(flag_a>1e-6 && n_iter<1e+10)/////////////////1) Active set
			    {
				 
                                ///phi_last=phi
				
                                Assign(P, Q, phi, phi_last);
			        Assign(P, Q, phi, phi_old);
                 
                               if(n_pick>0)
                                 for(i=0;i<n_pick;i++)
                                   {
					    cur_p=pick[i];
				          Update(cur_p, N, P, Q, lambda1, lambda2, C_m, X_m, Xnorm, E, Beta, Bnorm, phi_old, phi);
					    n_iter=n_iter+1;
				           }

				  /// update coefficients that is not penalized
				 if(n_unpick>0)
				  for(i=0; i<n_unpick; i++)
				  {
                                    p=unpick[i];
				    for(q=0; q<Q; q++)
					   if(C_m[p*Q+q]==2)
						    {
							temp=0;
							for(n=0;n<N;n++)
								  temp=temp+E[n*Q+q]*X_m[n*P+p];
							phi[p*Q+q]=temp/Xnorm[p]+phi_old[p*Q+q];

								///update residue
							for(n=0; n<N;n++)
								 E[n*Q+q]=E[n*Q+q]+(phi_old[p*Q+q]-phi[p*Q+q])*X_m[n*P+p];

								///update phi_old
							phi_old[p*Q+q]=phi[p*Q+q];
                                                         n_iter=n_iter+1;
					        }
				  }
                             flag_a=Dist(P, Q, phi_last, phi);
			    }///end of active set while(flag.a>1e-6)               
 
               //////////////// 2) Full loop
               Assign(P, Q, phi, phi_last);
               Assign(P, Q, phi, phi_old);
               for(p=0;p<P;p++)
               {
		   		Update(p, N, P, Q, lambda1, lambda2, C_m, X_m, Xnorm, E, Beta, Bnorm, phi_old, phi);
		   		n_iter=n_iter+1;
		}
		   flag=Dist(P, Q, phi_last, phi);
		     }//end of all loop while(flag>1e-6)

//////////// calculate Residue

  rss=0;
  for(n=0;n<N;n++)
     for(q=0;q<Q;q++)
	   {
		  rss=rss+E[n*Q+q]*E[n*Q+q];
	   }

Assign(N, Q, E, E_debug);

//////////// return phi
Assign(P, Q, phi, Phi_output);
Assign(P, Q, phi_last, Phi_debug);
*N_iter=n_iter;
*RSS=rss;

//////// free allocated variables
free(pick);
free(unpick);
free(Xnorm);
free(Beta);
free(phi);
free(phi_old);
free(phi_last);
free(E);
free(Bnorm);

}///end MultiRegGroupLasso function




//////////////////////////////////////
////// function to calculate degree
//////////////////////////////////////

void MultiRegGroupLassoDegree (int *NN, int *PP, int *QQ, double *X_m, double *Y_m, int *C_m, int *NN1, double *lam1, int *NN2, double *lam2,
                     double *Degree, double *Debug)
{
		/// Variables:
		/// N: Number of observations;
		/// P: Number of DNA;
		/// Q: Number of mRNA;

		/// X_m: matrix of N by P, input DNA data; each column has mean 0.
		/// Y_m: matrix of N by Q, input RNA data; each column has mean 0.
		/// C_m: matrix of P by Q, input data; indicator for "self relationship";
		///      C_m[i,j]=1 means the corresponding beta[i,j] will be penalized;		
            ///      C_m[i,j]=2 means the corresponding beta[i,j] will not be penalized.

		
                ///      C_m[i,j]=0 means the corresponding beta[i,j] will not be considered in the model.


		/// lambda1: group panelty; vector of lenght N1;
		/// lambda2: lasso panelty; vector of length N2;

		/// Degree: matrix of N1 by N2, estimated coefficients.

    int N, P, Q, N1, N2;
	int n,p,q, n1, n2;
	int i,j, k;
	int cur_p;
	int FIX;
	int *countN;

	double lambda1, lambda2;
	double temp;
	double *Xnorm;
	double *XYnorm;
	double *BetaLasso;
	double *Bnorm;


	N=*NN;
	P=*PP;
	Q=*QQ;
	N1=*NN1;
	N2=*NN2;

	countN=(int *) malloc (P*sizeof(int));
	Xnorm=(double *) malloc (P*sizeof(double));
	XYnorm=(double *) malloc (P*Q*sizeof(double));
	BetaLasso=(double *) malloc (P*Q*sizeof(double));
	Bnorm=(double *) malloc (P*sizeof(double));	

     //// initial value
		for(p=0;p<P;p++)
		 {
			 Xnorm[p]=0;
			 for(n=0;n<N;n++)
			   Xnorm[p]=Xnorm[p]+X_m[n*P+p]*X_m[n*P+p];
	      }

	    for(p=0;p<P;p++)
		   for(q=0;q<Q;q++)
		    {
				XYnorm[p*Q+q]=0;
				for(n=0;n<N;n++)
				   XYnorm[p*Q+q]=XYnorm[p*Q+q]+X_m[n*P+p]*Y_m[n*Q+q];
			}
        FIX=0;
        for(p=0;p<P;p++)
			 for(q=0;q<Q;q++)
				    {
                                     if(C_m[p*Q+q]==2)
						FIX=FIX+1;
					}	
        //// Calculate df for each pair of lam1 and lam2
	    for(n2=0;n2<N2;n2++)
	     {
			 lambda2=lam2[n2];		     
              //// 1. calculate lasso solution, and count non zero beta for each lambda2
	         for(p=0;p<P;p++)
	           {
				
                     countN[p]=0;
	             for(q=0;q<Q;q++)
		          {
				    if(C_m[p*Q+q]==1)
				     {
                                   BetaLasso[p*Q+q]=SoftShrink(XYnorm[p*Q+q], lambda2); 
                                   BetaLasso[p*Q+q]=BetaLasso[p*Q+q]/Xnorm[p];
	               		      if(BetaLasso[p*Q+q]!=0)
                                            countN[p]=countN[p]+1;
					 }
			      }
			    }			 
                   //// 2. calculate beta norm.
		     
                      CalBnorm(P, Q, BetaLasso, C_m, Bnorm);
			 //// 3. for each lam1, calculate df
			 
                     for(n1=0;n1<N1;n1++)
			   {
                      lambda1=lam1[n1];
			    Degree[n1*N2+n2]=0;
			    for(p=0;p<P;p++)
		   		{		 
                           if(Bnorm[p]*Xnorm[p]>lambda1)
				    {
 				     temp=countN[p]-lambda1*(countN[p]-1)/(Xnorm[p]*Bnorm[p]);
   			           Degree[n1*N2+n2]=temp+Degree[n1*N2+n2];		                    
                             }
		               }			
                          Degree[n1*N2+n2]= Degree[n1*N2+n2]+FIX;				
                    } ///end of n1 loop

	     }///end of n2 loop

*Debug=temp;

free(countN);
free(Xnorm);
free(BetaLasso);
free(XYnorm);
free(Bnorm);
}



////////////////////////////////////////
//////////////////  other functions
/////////////////////////////////////////

double SoftShrink(double y, double lam)
{
	double temp;
	double result;
	if(y>0)
	 temp=y-lam;
	else
	 temp=-y-lam;
	if(temp<=0)
	 result=0;
	else
	 {
		 result=temp;
		 if(y<0)
		  result=temp*(-1);
	  }
	return(result);
}


//////////////////
void CalBnorm(int P, int Q, double *Beta, int *C_m, double *Bnorm)
{
	int p,q;
    for(p=0;p<P;p++)
     {
		 Bnorm[p]=0;
         for(q=0;q<Q;q++)
         {
		   if(C_m[p*Q+q]==1)
            Bnorm[p]=Bnorm[p]+Beta[p*Q+q]*Beta[p*Q+q];
		 }
         Bnorm[p]=sqrt(Bnorm[p]);
     }
}

//////////////////
void Assign(int P, int Q, double *data_old, double *data_copy)
{
	int p, q;

	for(p=0;p<P;p++)
	 for(q=0;q<Q;q++)
	  data_copy[p*Q+q]=data_old[p*Q+q];
}

///////////////////
void Update(int cur_p, int N, int P, int Q, double lambda1, double lambda2, int *C_m, double *X_m, double *Xnorm, double *E, double *Beta, double *Bnorm, double *phi_old, double *phi)
{
	int n, p, q;
	double temp, temp1;
	p=cur_p;	
      ////// calculate lasso solution Beta
	for(q=0;q<Q;q++)
	 {
       if(C_m[p*Q+q]==0)
            Beta[p*Q+q]=0;
       else
       {
	    temp=0;
	    for(n=0;n<N;n++)
	      temp=temp+E[n*Q+q]*X_m[n*P+p];
    	temp1=temp+phi_old[p*Q+q]*Xnorm[p];
	    if(C_m[p*Q+q]==1)
	      Beta[p*Q+q]=SoftShrink(temp1, lambda2)/Xnorm[p];
	    else
		  Beta[p*Q+q]=temp1/Xnorm[p];
       }
	 }    
    ////// calculate norm of lasso solution Beta    
    Bnorm[p]=0;
    for(q=0;q<Q;q++)
    {
	  if(C_m[p*Q+q]==1)
        Bnorm[p]=Bnorm[p]+Beta[p*Q+q]*Beta[p*Q+q];
    }
    Bnorm[p]=sqrt(Bnorm[p]);
    ////// update phi    
    for(q=0;q<Q;q++)
    {
      if(C_m[p*Q+q]==0)
            phi[p*Q+q]=0;
	  else if(C_m[p*Q+q]==1 && Bnorm[p]>1e-6)
	   {
		temp=Xnorm[p]*Bnorm[p];
	    phi[p*Q+q]=Beta[p*Q+q]*SoftShrink(temp, lambda1)/temp;
	   }
	  else
	    phi[p*Q+q]=Beta[p*Q+q];
	}    
    ////// update residue    
    for(q=0;q<Q;q++)
      for(n=0;n<N;n++)
       {
		E[n*Q+q]=E[n*Q+q]+(phi_old[p*Q+q]-phi[p*Q+q])*X_m[n*P+p];
	   }
   /////// update phi_old    
   for(q=0;q<Q;q++)
      phi_old[p*Q+q]=phi[p*Q+q];
   /////// update Bnorm
    Bnorm[p]=0;
	    for(q=0;q<Q;q++)
	    {
		  if(C_m[p*Q+q]==1)
	        Bnorm[p]=Bnorm[p]+phi[p*Q+q]*phi[p*Q+q];
	    }
    Bnorm[p]=sqrt(Bnorm[p]);
}
//////////////////////////////
double Dist(int P, int Q, double *phi_last, double * phi)
{
	int p,q;
	double temp,result;	
      result=0;
	for(p=0;p<P;p++)
	 for(q=0;q<Q;q++)
	  {
		  temp=phi_last[p*Q+q]-phi[p*Q+q];
		  temp=fabs(temp);
		  if(temp>result)
		    result=temp;
      }
    return(result);
}
