#include "math.h"
#include "float.h"
#include "common.h"

/*
   [x1,y1,z1] dot [x2,y2,z2]
 */
double dot(double x1, double y1, double z1, double x2, double y2, double z2)
{
	return (x1 * x2 + y1 * y2 + z1 * z2);
}

double norm(double x1, double y1, double z1)
{
	return sqrt(dot(x1, y1, z1, x1, y1, z1));
}
/*
	perform cross product of vector a with vector b
*/
void cross(double a1, double a2, double a3, double b1, double b2, double b3, double *c1, double *c2, double *c3)
{
	*c1 = a2 * b3 - a3 * b2;
	*c2 = a3 * b1 - a1 * b3;
	*c3 = a1 * b2 - a2 * b1;
}
