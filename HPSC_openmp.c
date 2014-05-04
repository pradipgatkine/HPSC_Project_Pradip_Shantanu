//To set number of threads, in terminal: export OMP_NUM_THREADS=8
//Reference: https://computing.llnl.gov/tutorials/openMP/samples/C/omp_hello.c
//To compile: gcc -fopenmp HPSC_openmp.c
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <malloc.h>
#include <time.h>
double max(double a, double b, double c){
	double ans=a;
	if (ans<b){
		ans=b;
	}
	if (ans<c){
		ans=c;
	}
	return ans;
}

double absolute(double a){
	double ans;
	if (a>=0){
		ans=a;
	}
	else{
		ans=-a;
	}
	return ans;
}

int main(){
	double Lx = 0.40; // length along X (m)
	double Ly = 0.40; // length along Y (m)
	double delta_x = 0.002; // length of CV along X(m)
	double delta_y = 0.002; // length of CV along Y(m)
	int nx = Lx/delta_x; // no of CVs along X
	int ny = Ly/delta_y; // no of CVs along Y
	double nu = 0.004; // kinematic viscosity (m^2/s)
	double mu = nu; // dynamic viscosity ()
	double rho = 1; // density (kg/m^3)
	double rel_u=0.5;
	double rel_v=0.5;
	double rel_p=0.5;
	int Nx=nx+2;	
	int Ny=ny+2;
	int i,j,k,u_no,v_no,b_count,p_no;
	double max_b, delx_e, delx_w, dely_n, dely_s, De_u, Dw_u, Dn_u, Ds_u, De_v, Dw_v, Dn_v, Ds_v;
	double Fe_u, Fw_u, Fn_u, Fs_u, Fe_v, Fw_v, Fn_v, Fs_v, ae_u, aw_u, an_u, as_u, ae_v, aw_v, an_v, as_v, ap_u, ap_v;
	double ae_p, aw_p, an_p, as_p, ap_p;
	// Initialise with 0
	double p[ny][nx], u[Ny][Nx-1], v[Ny-1][Nx];
	double Ap_u[Ny][Nx-1], Ap_v[Ny-1][Nx], b[ny][nx];
	double uc[Ny][Nx-1], vc[Ny-1][Nx], pc[ny][nx];
	double converge=1e-6;//Convergance criterion
	
	for (i=0;i<ny;i++){
		for (j=0;j<nx;j++){
			p[i][j]=0;
			b[i][j]=0;
		}
	}
	for (i=0;i<Ny;i++){
		for (j=0;j<Nx-1;j++){
			u[i][j]=0;
			uc[i][j]=0;
			Ap_u[i][j]=0;
		}
	}
	for (i=0;i<Ny-1;i++){
		for (j=0;j<Nx;j++){
			v[i][j]=0;
			vc[i][j]=0;
			Ap_v[i][j]=0;
		}
	}

	// Boundary Condition
	for (j=0;j<Nx-1;j++){
		u[0][j]=4;
	}
	int count=0;
	while (1){
		count=count+1;
	    // X mom GS
	    for (u_no=1;u_no<6;u_no++){
	    	for (i=1;i<Ny-1;i++){ 
	    		#pragma omp parallel for private(j)
				for (j=1;j<Nx-2;j++){
					delx_e=delta_x;
	                delx_w=delta_x;
	                dely_n=delta_y;
	                dely_s=delta_y;
	                if (i==1){ 
	                    dely_n=delta_y/2;
	                }
	                else if (i==Ny-2){ 
	                    dely_s=delta_y/2;
	                }
	                De_u=mu/delx_e;
	                Dw_u=mu/delx_w;
	                Dn_u=mu/dely_n;
	                Ds_u=mu/dely_s;
	                Fe_u=rho*(u[i][j]+u[i][j+1])/2;
	                Fw_u=rho*(u[i][j]+u[i][j-1])/2;
	                Fn_u=rho*(v[i-1][j]+v[i-1][j+1])/2;
	                Fs_u=rho*(v[i][j]+v[i][j+1])/2;
	                ae_u=max(-Fe_u,De_u-Fe_u/2,0.0);
	                aw_u=max(Fw_u,Dw_u+Fw_u/2,0);	
	                an_u=max(-Fn_u,Dn_u-Fn_u/2,0);
	                as_u=max(Fs_u,Ds_u+Fs_u/2,0);
	                ap_u=ae_u+aw_u+an_u+as_u;
	                Ap_u[i][j]=ap_u;
	                if (ap_u==0){//This is just error condition, put to avoid division by zero
	                	printf("Error: ap_u is zero %i %i \n",i,j);
	                	//break;
	                	//goto end;
	                }
	                u[i][j]=rel_u*(ae_u*u[i][j+1]+aw_u*u[i][j-1]+an_u*u[i-1][j]+as_u*u[i+1][j]+(p[i-1][j-1]-p[i-1][j])*delta_y)/ap_u+(1-rel_u)*u[i][j];
				}
			}
		}
                
	    // Y mom GS
	    for (v_no=1;v_no<6;v_no++){
	    	for (i=1;i<Ny-2;i++){
	    		for (j=1;j<Nx-1;j++){
					delx_e=delta_x;
	                delx_w=delta_x;
	                dely_n=delta_y;
	                dely_s=delta_y;
	                if (j==1){
	                    delx_w=delta_x/2;
	                }
	                else if (j==Nx-1){ 
	                    delx_e=delta_x/2;
	                }
	                De_v=mu/delx_e;
	                Dw_v=mu/delx_w;
	                Dn_v=mu/dely_n;
	                Ds_v=mu/dely_s;
	                Fe_v=rho*(u[i][j]+u[i+1][j])/2;
	                Fw_v=rho*(u[i][j-1]+u[i+1][j-1])/2;
	                Fn_v=rho*(v[i][j]+v[i-1][j])/2;
	                Fs_v=rho*(v[i][j]+v[i+1][j])/2;
	                ae_v=max(-Fe_v,De_v-Fe_v/2,0);
	                aw_v=max(Fw_v,Dw_v+Fw_v/2,0);
	                an_v=max(-Fn_v,Dn_v-Fn_v/2,0);
	                as_v=max(Fs_v,Ds_v+Fs_v/2,0);
	                ap_v=ae_v+aw_v+an_v+as_v;
	                Ap_v[i][j]=ap_v;
	                if (ap_v==0){//This is just error condition, put to avoid division by zero
	                	printf("Error: ap_v is zero %i %i \n",i,j);	                
	                	//goto end;
	                	break;
	                }
	                v[i][j]=rel_v*(ae_v*v[i][j+1]+aw_v*v[i][j-1]+an_v*v[i-1][j]+as_v*v[i+1][j]+(p[i][j-1]-p[i-1][j-1])*delta_x)/ap_v+(1-rel_v)*v[i][j];
				}
			}
		}	               
     
		b_count=0;
		max_b=0;
		for (i=0;i<ny;i++){
			for (j=0;j<nx;j++){
				b[i][j]=(u[i+1][j]-u[i+1][j+1])*delta_y+(v[i+1][j+1]-v[i][j+1])*delta_x;
				if (max_b<absolute(b[i][j])){
					max_b=absolute(b[i][j]);
				}
			}
		}
		printf("Current value: %f ", max_b);
		printf("Will stop when the value is <=%f ",converge);
		printf("Iteration no.: %i\n", count);
		if (max_b<converge){//Convergence Criterion
			break;
		}
		else{
			for (i=0;i<ny;i++){
				for (j=0;j<nx;j++){ //initialise as zeros in each iteration
					pc[i][j]=0;
				}
			} 
	        // Pressure correction GS
			for (p_no=1;p_no<6;p_no++){
				for (i=0;i<ny;i++){
					for (j=0;j<nx;j++){ 					
    					ae_p=delta_y*delta_y/Ap_u[i+1][j+1];
    					aw_p=delta_y*delta_y/Ap_u[i+1][j];
    					an_p=delta_x*delta_x/Ap_v[i][j+1];
    					as_p=delta_x*delta_x/Ap_v[i+1][j+1];
    					if (i==0){
    						an_p=0;
    					}
    					else if (i==ny-1){
    						as_p=0;
    					}
    					if (j==0){
    						aw_p=0;
    					}
    					else if (j==nx-1){
    						ae_p=0;
    					}
    					ap_p=ae_p+aw_p+an_p+as_p;
    					if (i==0){
    						if (j==0){
    							pc[i][j]=(ae_p*pc[i][j+1]+as_p*pc[i+1][j]+b[i][j])/ap_p;
    						}
    						else if (j==nx-1){
    							pc[i][j]=(aw_p*pc[i][j-1]+as_p*pc[i+1][j]+b[i][j])/ap_p;
    						}
    						else{
    							pc[i][j]=(ae_p*pc[i][j+1]+aw_p*pc[i][j-1]+as_p*pc[i+1][j]+b[i][j])/ap_p;
    						}
    					}
    					else if (i==ny-1){
    						if (j==0){
    							pc[i][j]=(ae_p*pc[i][j+1]+an_p*pc[i-1][j]+b[i][j])/ap_p;
    						}
    						else if (j==nx-1){
    							pc[i][j]=(aw_p*pc[i][j-1]+an_p*pc[i-1][j]+b[i][j])/ap_p;
    						}
    						else{
    							pc[i][j]=(ae_p*pc[i][j+1]+aw_p*pc[i][j-1]+an_p*pc[i-1][j]+b[i][j])/ap_p;
    						}
    					}
    					else{
    						if (j==0){
    							pc[i][j]=(ae_p*pc[i][j+1]+an_p*pc[i-1][j]+as_p*pc[i+1][j]+b[i][j])/ap_p;
    						}
    						else if (j==nx-1){
    							pc[i][j]=(aw_p*pc[i][j-1]+an_p*pc[i-1][j]+as_p*pc[i+1][j]+b[i][j])/ap_p;
    						}
    						else{
    							pc[i][j]=(ae_p*pc[i][j+1]+aw_p*pc[i][j-1]+an_p*pc[i-1][j]+as_p*pc[i+1][j]+b[i][j])/ap_p;
    						}
    					}
    				}
    			}
    		}
    	}
                                
        // corrections
        for (i=0;i<ny;i++){
        	for (j=0;j<nx;j++){
        		p[i][j]=p[i][j]+rel_p*pc[i][j];
        	}
        }
        for (i=1;i<Ny-1;i++){
        	for (j=1;j<Nx-2;j++){
        		uc[i][j]=delta_y*(pc[i-1][j-1]-pc[i-1][j])/Ap_u[i][j];
                u[i][j]=u[i][j]+uc[i][j];
            }
        }
        for (i=1;i<Ny-2;i++){
        	for (j=1;j<Nx-1;j++){
        		vc[i][j]=delta_x*(pc[i][j-1]-pc[i-1][j-1])/Ap_v[i][j];
                v[i][j]=v[i][j]+vc[i][j];    
            }
        }
	}
	printf("Total Number of Iterations: %i\n", count);
	/*FILE *fp;
    fp = fopen("u_velocity.txt", "w+");
	for (i=0;i<Ny;i++){
		for (j=0;j<Nx-1;j++){
			fprintf(fp,"%f ", u[i][j]);
		}
		fprintf(fp,"\n");
	}
    fclose(fp);
    fp = fopen("v_velocity.txt", "w+");
	for (i=0;i<Ny-1;i++){
		for (j=0;j<Nx;j++){
			fprintf(fp,"%f ", v[i][j]);
		}
		fprintf(fp,"\n");
	}
	fp = fopen("max_b.txt", "w+");
	for (i=0;i<count;i++){
			fprintf(fp,"%f \n", max_b_record[i]);
	}*/
	end://Refering to "goto" in X-mom and Y_mom loops; just in case there is error
		return 0;
}     
