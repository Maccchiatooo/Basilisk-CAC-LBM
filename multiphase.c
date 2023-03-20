//#include "grid/multigrid.h"
#include "view.h"
#include "utils.h"
#include "vtk.h"
#include "isotropic.h"
#include "macro.h"
#include "lbm.h"

#define lmin (6)
#define lmax (8)
int main(){

    L0=pow(2,lmax);
    init_grid(1<<lmin);

    run();

}


 event dump_out(i+=1000){
    char ti[100];
    sprintf(ti,"%d",i/1000);  
    dump(ti);
 }

event adapt(i++)
{
    scalar cri[];
    foreach(){
        cri[]=sqrt(sq(d_phi.x[])+sq(d_phi.y[]));
    }
    adapt_wavelet({cri}, (double[]){1.e-6}, maxlevel = (lmax), minlevel = (lmin));
}
event end (i = 10000){}
