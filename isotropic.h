inline double iso_x(Point point, scalar phi,int del){
    return(((phi[1,0]-phi[-1,0])/3.0+(phi[1,1]-phi[-1,-1])/12.0+(phi[1,-1]-phi[-1,1])/12.0));

} 
inline double iso_y(Point point, scalar phi,int del){
    return(((phi[0,1]-phi[0,-1])/3.0+(phi[1,1]-phi[-1,-1])/12.0+(phi[-1,1]-phi[1,-1])/12.0));

} 

inline double iso_laplace(Point point, scalar phi,int i, int j ){
    return((phi[i+1,j+1]+phi[i-1,j+1]+phi[i-1,j-1]+phi[i+1,j-1]+4.0*phi[i+1,j]+4.0*phi[i-1,j]+4.0*phi[i,j+1]+4.0*phi[i,j-1]-20*phi[i,j])/6.0);

} 