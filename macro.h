#define cs2 1./3.
#define cs 1./sqrt(3.)
#define q 9
#define D_delta 8


double sigma=0.0001;


double t_1=0.1;
double t_2=0.1;
double t_phi=0.5;

scalar rho[];
scalar phi[];

scalar t_p[];

scalar nu[];
scalar div_phi_x[];
scalar div_phi_y[];
scalar div_phi[];

scalar div_[];
scalar p[];
scalar pp[];

vector u[];
vector d_phi[];
vector d_p[];
vector d_pp[];
vector d_u[];


scalar f0[], f1[], f2[], f3[], f4[], f5[], f6[], f7[], f8[];
scalar * f = {f0, f1, f2, f3, f4, f5, f6, f7, f8};

scalar g0[], g1[], g2[], g3[], g4[], g5[], g6[], g7[], g8[];
scalar * g = {g0, g1, g2, g3, g4, g5, g6, g7, g8};
const double wt[q]={4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.};
const coord e[q]={{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};
//////////////////////////////////////////////////////////

///rising bubble simulation///
