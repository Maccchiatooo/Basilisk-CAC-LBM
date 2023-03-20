#include "grid/multigrid.h"
#include "view.h"
#include "utils.h"
#include "vtk.h"
#include "isotropic.h"
#include "macro.h"
#include "lbm.h"

#define lmax ((int)log2(L0))

int main(){

    L0=128;
    init_grid(1<<7);

    run();

}


 event dump_out(i+=1000){
    char ti[100];
    sprintf(ti,"%d",i/1000);  
    dump(ti);
 }

event end (i = 100000){}
