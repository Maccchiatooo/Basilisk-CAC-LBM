#include "run.h"
#include "view.h"

double rho_1=0.1;
double rho_2=1;  

event init(t=0){
    foreach(){
        double dist =sq(x-0.5*L0)+sq(y-0.5*L0)-sq(0.25*L0);
        phi[]=0.5-0.5*tanh(2.0*dist/D_delta);
        rho[]=rho_1*phi[]+rho_2*(1.0-phi[]);
        t_p[]=t_1*phi[]+t_2*(1.0-phi[]);
        nu[]=t_p[]*cs2;
        p[]=0.0;
        pp[]=0.0;
        foreach_dimension(){
            u.x[]=u.y[]=0;
        }
        
    }
    foreach(){

        d_phi.x[]=iso_x(point,phi,Delta);
        d_phi.y[]=iso_y(point,phi,Delta);
    }
    foreach(){
        double sqd=sqrt(sq(d_phi.x[])+sq(d_phi.y[]))+1e-15;
        div_phi_x[]=d_phi.x[]/sqd;
        div_phi_y[]=d_phi.y[]/sqd;
    }
    foreach(){
        div_phi[]=iso_x(point,div_phi_x,Delta)+iso_y(point,div_phi_y,Delta);
    }

    foreach(){
        int i=0;
        for(scalar f_temp in f){
            double edu = e[i].x*u.x[]+e[i].y*u.y[];
            double udu = u.x[]*u.x[]+u.y[]*u.y[];
            double gamma = wt[i]*(1.+edu/cs2+edu*edu/2/cs2/cs2-udu/2/cs2);
            
            f_temp[]=wt[i]*pp[]+(gamma-wt[i])*cs2;
            
            i++;
        }

        i=0;
        for(scalar g_temp in g){
            double edu = e[i].x*u.x[]+e[i].y*u.y[];
            double udu = u.x[]*u.x[]+u.y[]*u.y[];
            double gamma = wt[i]*(1.+edu/cs2+edu*edu/2/cs2/cs2-udu/2/cs2);

            g_temp[]=gamma*phi[];
            i++;
        }
    }


}


event collision_f(i++){

    if(i>0){
    foreach(){
        double sqd=sqrt(sq(d_phi.x[])+sq(d_phi.y[]))+1e-15;

        double fors_sx=-3.0/2.0*sigma*D_delta*div_phi[]*sqd*d_phi.x[];
        double fors_sy=-3.0/2.0*sigma*D_delta*div_phi[]*sqd*d_phi.y[];

        double fors_px=-iso_x(point,p,Delta);
        double fors_py=-iso_y(point,p,Delta);

        double fors_ppx=iso_x(point,pp,Delta);
        double fors_ppy=iso_y(point,pp,Delta);

        double fors_visx=nu[]*(2.0*iso_x(point,u.x,Delta)*iso_x(point,rho,Delta)+
                                  (iso_y(point,u.x,Delta)+iso_x(point,u.y,Delta))*iso_y(point,rho,Delta));

        double fors_visy=nu[]*(2.0*iso_y(point,u.y,Delta)*iso_y(point,rho,Delta)
                                 +(iso_y(point,u.x,Delta)+iso_x(point,u.y,Delta))*iso_x(point,rho,Delta));


        
        int ii=0;

        for(scalar f_temp in f){
            double edu = e[ii].x*u.x[]+e[ii].y*u.y[];
            double udu = u.x[]*u.x[]+u.y[]*u.y[];
            double gamma = wt[ii]*(1.+edu/cs2+edu*edu/2/cs2/cs2-udu/2/cs2);

            double fors=(e[ii].x-u.x[])*gamma/rho[]*(fors_px+fors_visx+fors_sx)+
                        (e[ii].y-u.y[])*gamma/rho[]*(fors_py+fors_visy+fors_sy)+
                        (e[ii].x-u.x[])*wt[ii]*fors_ppx+
                        (e[ii].y-u.y[])*wt[ii]*fors_ppy;

            double feq=gamma*pp[]+(gamma-wt[ii])*cs2-0.5*fors;
            f_temp[]-=(f_temp[]-feq)/(t_p[]+0.5)+fors;
            ii++;
        }

    }
    }
}


event collision_g(i++){
    if(i>0){
    foreach(){

        int ii=0;

        for(scalar g_temp in g){

            double gamma = wt[ii]*(1.+3.0*(e[ii].x*u.x[]+e[ii].y*u.y[])+4.5*sq(e[ii].x*u.x[]+e[ii].y*u.y[])-1.5*sq(sq(u.x[])+sq(u.y[])));


            double fors_phi= gamma*((e[ii].x-u.x[])*4.0*phi[]*(1.0-phi[])/D_delta*div_phi_x[]+
                                    (e[ii].y-u.y[])*4.0*phi[]*(1.0-phi[])/D_delta*div_phi_y[]);

            double geq=gamma*phi[]-0.5*fors_phi;

            g_temp[]=g_temp[]-(g_temp[]-geq)/(t_phi+0.5)+fors_phi;
            ii++;
        }
    }
    }
}

// // event boundary(i++){

// // }
scalar f_new[];
event streaming_f(i++){
    if(i>0){
    int ii=0;
    for(scalar f_old in f){
        foreach(){
            double cfl=1./Delta;
            if(cfl<=0.99)
            f_new[]=cfl / 2. * (1. + cfl) * f_old[-(int)e[ii].x, -(int)e[ii].y] + (1. - cfl * cfl) * f_old[] - cfl / 2. * (1. - cfl) * f_old[(int)e[ii].x, (int)e[ii].y];
                else
                f_new[] = f_old[-(int)e[ii].x, -(int)e[ii].y];
        }
            foreach(){
            f_old[]=f_new[];
        }


        ii++;
    }

    }
}

scalar g_new[];
event streaming_g(i++){
    if(i>0){
    int ii=0;
    for(scalar g_old in g){
        foreach(){
            double cfl=1./Delta;
            if(cfl<=0.99)
            g_new[]=cfl / 2. * (1. + cfl) * g_old[-(int)e[ii].x, -(int)e[ii].y] + (1. - cfl * cfl) * g_old[] - cfl / 2. * (1. - cfl) * g_old[(int)e[ii].x, (int)e[ii].y];
                else
                g_new[] = g_old[-(int)e[ii].x, -(int)e[ii].y];
        }

            foreach(){
            g_old[]=g_new[];}
        ii++;
    }
    }

}



event update(i++){

    if(i>0){
    foreach(){
        phi[]=0.0;
        for(scalar g_temp in g){
            phi[]+=g_temp[];
        }
        
    }

    foreach(){
        rho[]=rho_1*phi[]+rho_2*(1.0-phi[]);
        t_p[]=t_1*phi[]+t_2*(1.0-phi[]);
        nu[]=t_p[]*cs2;
    }

    foreach(){

        d_phi.x[]=iso_x(point,phi,Delta);
        d_phi.y[]=iso_y(point,phi,Delta);
    }
    foreach(){
        double sqd=sqrt(sq(d_phi.x[])+sq(d_phi.y[]))+1e-15;
        div_phi_x[]=d_phi.x[]/sqd;
        div_phi_y[]=d_phi.y[]/sqd;

    }
    foreach(){
        div_phi[]=iso_x(point,div_phi_x,Delta)+iso_y(point,div_phi_y,Delta);
    }

    foreach(){
        pp[]=0.0;
        for(scalar f_temp in f)
            pp[]+=f_temp[];
        
        

        pp[]-=u.x[]*iso_x(point,pp,Delta)/2.0-u.y[]*iso_y(point,pp,Delta)/2.0;    
        
        p[]=pp[]*rho[];
    }

    foreach(){
            u.y[]=0.0;
            u.x[]=0.0;
            int ii=0;
            for(scalar f_temp in f){
               
                u.x[]+=f_temp[]*e[ii].x/cs2;
                u.y[]+=f_temp[]*e[ii].y/cs2;
                ii++;
            }
    }

    foreach(){

        double sqd=sqrt(sq(d_phi.x[])+sq(d_phi.y[]))+1e-15;

        double fors_sx=-3.0/2.0*sigma*D_delta*div_phi[]*sqd*d_phi.x[];
        double fors_sy=-3.0/2.0*sigma*D_delta*div_phi[]*sqd*d_phi.y[];

        double fors_px=-iso_x(point,p,Delta);
        double fors_py=-iso_y(point,p,Delta);

        double fors_ppx=iso_x(point,pp,Delta)*rho[];
        double fors_ppy=iso_y(point,pp,Delta)*rho[];

        double fors_visx=nu[]*(2.0*iso_x(point,u.x,Delta)*iso_x(point,rho,Delta)+
                                  (iso_y(point,u.x,Delta)+iso_x(point,u.y,Delta))*iso_y(point,rho,Delta));

        double fors_visy=nu[]*(2.0*iso_y(point,u.y,Delta)*iso_y(point,rho,Delta)
                                 +(iso_y(point,u.x,Delta)+iso_x(point,u.y,Delta))*iso_x(point,rho,Delta));


        u.x[]+=(fors_px+fors_ppx+fors_visx+fors_sx)/2.0/rho[];
        u.y[]+=(fors_py+fors_ppy+fors_visy+fors_sy)/2.0/rho[];
        }
    }
}
