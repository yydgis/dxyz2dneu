
#include <iostream>
#include <vector>
#include <cstring>
#include <cmath>

#define PI          3.1415926535897932  /* pi */
#define D2R         (PI/180.0)          /* deg to rad */
#define R2D         (180.0/PI)          /* rad to deg */
#define RE_WGS84    6378137.0           /* earth semimajor axis (WGS84) (m) */
#define FE_WGS84    (1.0/298.257223563) /* earth flattening (WGS84) */

/* inner product ---------------------------------------------------------------
* inner product of vectors
* args   : double *a,*b     I   vector a,b (n x 1)
*          int    n         I   size of vector a,b
* return : a'*b
*-----------------------------------------------------------------------------*/
extern double dot(const double *a, const double *b, int n)
{
    double c=0.0;
    
    while (--n>=0) c+=a[n]*b[n];
    return c;
}
/* transform ecef to geodetic postion ------------------------------------------
* transform ecef position to geodetic position
* args   : double *r        I   ecef position {x,y,z} (m)
*          double *pos      O   geodetic position {lat,lon,h} (rad,m)
* return : none
* notes  : WGS84, ellipsoidal height
*-----------------------------------------------------------------------------*/
extern void ecef2pos(const double *r, double *pos)
{
    double e2=FE_WGS84*(2.0-FE_WGS84),r2=dot(r,r,2),z,zk,v=RE_WGS84,sinp;
    
    for (z=r[2],zk=0.0;fabs(z-zk)>=1E-4;) {
        zk=z;
        sinp=z/sqrt(r2+z*z);
        v=RE_WGS84/sqrt(1.0-e2*sinp*sinp);
        z=r[2]+v*e2*sinp;
    }
    pos[0]=r2>1E-12?atan(z/sqrt(r2)):(r[2]>0.0?PI/2.0:-PI/2.0);
    pos[1]=r2>1E-12?atan2(r[1],r[0]):0.0;
    pos[2]=sqrt(r2+z*z)-v;
}

/* scale for lat to north, lon to east */
static double lat2local(double lat, double* lat2north)
{
	double f_WGS84 = FE_WGS84;
	double e2WGS84 = (2.0 * f_WGS84 - f_WGS84 * f_WGS84);
	double slat = sin(lat);
	double clat = cos(lat);
	double one_e2_slat2 = 1.0 - e2WGS84 * slat * slat;
	double Rn = RE_WGS84 / sqrt(one_e2_slat2);
	double Rm = Rn * (1.0 - e2WGS84) / (one_e2_slat2);
	*lat2north = Rm;
	return Rn * clat;
}

int main(int argc, const char *argv[])
{
  if (argc > 6)
  {
    double xyz1[3] = { 0 };
    double xyz2[3] = { 0 };
    xyz1[0] = atof(argv[1]);
    xyz1[1] = atof(argv[2]);
    xyz1[2] = atof(argv[3]);
    xyz2[0] = atof(argv[4]);
    xyz2[1] = atof(argv[5]);
    xyz2[2] = atof(argv[6]);

    double blh1[3] = { 0 };
    double blh2[3] = { 0 };
    ecef2pos(xyz1,blh1);
    ecef2pos(xyz2,blh2);
    double l2n = 0;
    double l2e = lat2local(blh1[0], &l2n);
    printf("dXYZ = %10.4f,%10.4f,%10.4f\n", xyz2[0]-xyz1[0], xyz2[1]-xyz1[1], xyz2[2]-xyz1[2]);
    printf("dNEU = %10.4f,%10.4f,%10.4f\n", (blh2[0]-blh1[0])*l2n, (blh2[1]-blh1[1])*l2e, (blh2[2]-blh1[2]));

  }
  return 0;
}
