#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "float.h"
#include "common.h"

/*
  Set up buffers for different mosiacing routines - e.g. speckleTrack
*/
void setupBuffers(outputImageStructure *outputImage, float ***vXimage, float ***vYimage, float ***vZimage,float ***scaleX,float ***scaleY, float ***scaleZ,
		  float ***vxTmp, float ***vyTmp, float ***vzTmp,float ***sxTmp,float ***syTmp, float ***fScale,float ***errorX,float ***errorY) {
	*vXimage = (float **)outputImage->image;
	*vYimage = (float **)outputImage->image2;
	*vZimage = (float **)outputImage->image3;
	*scaleX = (float **)outputImage->scale;
	*scaleY = (float **)outputImage->scale2;
	*scaleZ = (float **)outputImage->scale3;
	*vxTmp = (float **) outputImage->vxTmp;
	*vyTmp = (float **) outputImage->vyTmp;
	*vzTmp = (float **) outputImage->vzTmp;
	*sxTmp = (float **) outputImage->sxTmp;
	*syTmp = (float **) outputImage->syTmp;
	*fScale= (float **) outputImage->fScale; 
	*errorX=(float **)outputImage->errorX;
	*errorY=(float **)outputImage->errorY;
}

/*
  Compute tide error and add to phase error, twok=2pi/lambda for phase, 1.0 for offsets
*/
void interpTideError(double *phaseError,  inputImageStructure *phaseImage,	vhParams *params,double x,double y, double psi, double twok) {
	double tdTemp,tideError;

	if(phaseImage->tideDiffFlag == TRUE ) {
		tdTemp=interpTideDiff(x,y,phaseImage->tideDiff);
		if(tdTemp > MINELEVATION) phaseImage->tideCorrection=tdTemp*365.25/(double)params->nDays;
		else phaseImage->tideCorrection=0.0;
		interpTideDiff(x,y,phaseImage->tideDiff)*365.25/(double)params->nDays;
		tideError= 0.1 * cos(psi) *twok;
		*phaseError=sqrt((*phaseError) * (*phaseError) + tideError*tideError);
	}
}

/* 
	Compute theta and other info: not for now compute thetaCfixed and thetaC the same way 
*/
conversionDataStructure *setupGeoConversions (inputImageStructure *currentImage, float *azSLPixSize, float *rSLPixSize, double *Re, double *ReH,
					      double *thetaC,double *ReHfixed, double *thetaCfixedReH) {
	conversionDataStructure *cP;

	initllToImageNew( currentImage);
	*azSLPixSize=currentImage->azimuthPixelSize/currentImage->nAzimuthLooks;
	*rSLPixSize=currentImage->rangePixelSize/currentImage->nRangeLooks;
	cP = &(currentImage->cpAll);
	*Re = cP->Re;
	if(cP->ReH != NULL) *ReH =cP->ReH[(int) (currentImage->azimuthSize)/2]; 
	else { error(" setup setupGeoConversion *ReH = cP->Re + cP->H"); }
	/* This uses the central value state-vector determined value */
	*thetaC = thetaRReZReH( cP->RCenter, (*Re+0), *ReH);
	/* This uses the nominal value from the geodat and corresponds to what was used by changeflat
	It is needed to compensate for the changeflat corrections in phase calculations */
	*ReHfixed=*Re + currentImage->par.H;		
	*thetaCfixedReH=thetaRReZReH(cP->RCenter, (*Re+0), *ReHfixed);
	return(cP);
}

/*
  Compute corrections based on shelf mask
*/
float shelfMaskCorrection(inputImageStructure *currentImage, vhParams *currentParams, unsigned char sMask, double x, double y, double psi, double*sigmaR) {
	double tdTemp,drCorrect,tideError;
	if(sMask == SHELF) {
		/* DeltaTideError=-vz * cospsi * deltaT */   
		if(currentImage->tideDiffFlag == TRUE ) {
			tdTemp=interpTideDiff(x,y,currentImage->tideDiff);
			if(tdTemp > MINELEVATION) currentImage->tideCorrection=tdTemp*365.25/(double)currentParams->nDays;
			else currentImage->tideCorrection=0.0;
		}
		drCorrect = - currentImage->tideCorrection * cos(psi) * (double)currentParams->nDays/365.25;
		/* Added tide error 5/21/03 */
		tideError= 0.1 * cos(psi);
		*sigmaR=sqrt( (*sigmaR)*(*sigmaR) + tideError*tideError);
	} else { currentImage->tideCorrection=0; drCorrect=0.0;}
	return( drCorrect);
}

/*
  Compute Range,theta,thetatd, and ReH for azimuth, range coords
*/
void geometryInfo(conversionDataStructure *cP,inputImageStructure *currentImage,double azimuth,double range, double z, 
					double thetaC,double *ReH,double *Range,double *theta,double *thetaD, double *psi,double zSp) {
	
	if(cP->ReH == NULL) error("geometryInfo : ReH[] not defined");
	if(azimuth < cP->azSize && azimuth >= 0) *ReH = cP->ReH[(int)azimuth];
	else *ReH = cP->ReH[(int) cP->azSize/2]; /* this should not actually get used, since no valid azimuth */
	*Range = cP->RNear + range * currentImage->rangePixelSize;
	/* 
	   Compute look and inc angles
	*/
	*theta = thetaRReZReH(*Range, (cP->Re+zSp), *ReH);
	*thetaD = *theta - thetaC;
	*psi = psiRReZReH(*Range,(cP->Re+zSp),*ReH);
}

/*
  Compute error, return based on covariance matrix
*/
double computeSig2Base(double sinThetaD,double cosThetaD, double azimuth, inputImageStructure *inputImage,Offsets *offsets)
{
	double v[7], tmpV[7];
	double sig2Base;
	double imageLength, normAzimuth, xsq;
	int i,j;
	/*
	  Compute azimuth in image coord stuff
	*/
	imageLength = (double)inputImage->azimuthSize;   
	normAzimuth = (azimuth-0.5*imageLength)/imageLength;
	xsq=normAzimuth * normAzimuth;
	/*
	  Vector used to comptue baseline error 
	*/
	v[1]= - sinThetaD;
	v[2]= - sinThetaD* normAzimuth;	v[3]= - cosThetaD * normAzimuth;
	v[4]= - sinThetaD* xsq;			v[5]= - cosThetaD * xsq;
	if(offsets->deltaB == DELTABNONE) v[6]= 1; /* const */
	else v[6]= -cosThetaD;  /* with state vectors correction to Bp */
	/*
	  C*v
	*/
	for(i=1; i <= 6; i++) {
		tmpV[i]=0;
		for(j=1; j<=6; j++)  {tmpV[i] += offsets->Cr[i][j] * v[j]; }
	}
	/*
	  sigma^2= vT * C*v
	*/
	sig2Base=0.0;
	for(j=1; j<=6; j++)  sig2Base += tmpV[j] * v[j];
	return(sig2Base);
}

/* 
   Compute errors based on azimuth fit covariance matrix.
 */
double computeSig2AzParam(double sinTheta,double cosTheta, double azimuth, double Range, inputImageStructure *inputImage,Offsets *offsets)
{
	double v[5], tmpV[5];
	double sig2Off;
	double imageLength, normAzimuth;
	int i,j;
	/*
	  Compute azimuth in image coord stuff
	*/
	imageLength = (double)inputImage->azimuthSize;   
	normAzimuth = (azimuth-0.5*imageLength)/imageLength;
	/*
	  Vector used to comptue baseline error 
	*/
	v[1]= 1;
	v[2]= Range * sinTheta;	v[3] = Range * cosTheta; 
	v[4]= normAzimuth;
	/*
	  C*v
	*/
	for(i=1; i <= 4; i++) {
		tmpV[i]=0;
		for(j=1; j<=4; j++)  tmpV[i] += offsets->Ca[i][j] * v[j];
	}
	/*
	  sigma^2= vT * C*v
	*/
	sig2Off=0.0;
	for(j=1; j<=4; j++)  sig2Off += tmpV[j] * v[j];
	return(sig2Off);
}

/*
  Compute velocity determination matrix
*/
void computeA(double lat,double lon,double x,double y,    inputImageStructure *aPhaseImage,  inputImageStructure *dPhaseImage,  double A[2][2])
{
	extern int HemiSphere;	
	double alpha, beta;
	double aHAngle, dHAngle,xyAngle;
	conversionDataStructure *aCp,*dCp;
	double invSin2;
	double twok;
	double cosAlphaBeta, cosAlpha, cosBeta, sinBeta, sinAlphaBeta;
	/* Get geometry params */
	aCp = &(aPhaseImage->cpAll);
	dCp = &(dPhaseImage->cpAll);
	/* Heading angles */
	aHAngle = computeHeading(lat,lon,0.0,aPhaseImage,aCp);
	dHAngle = computeHeading(lat,lon,0.0,dPhaseImage,dCp);
	/* other angles */
	alpha = aHAngle - dHAngle;
	xyAngle = atan2(-y,-x);
	if(HemiSphere == SOUTH) xyAngle += PI;
	beta = xyAngle - aHAngle;
	/*  Compute A matrix	*/
	invSin2 = 1.0/pow(sin(alpha),2.0);
	cosAlphaBeta = cos(alpha+beta);
	cosAlpha = cos(alpha);
	cosBeta = cos(beta);
	sinBeta = sin(beta);
	sinAlphaBeta = sin(alpha+beta);
	A[0][0] = invSin2*(cosBeta       -cosAlpha * cosAlphaBeta);
	A[0][1] = invSin2*(cosAlphaBeta - cosAlpha * cosBeta);
	A[1][0] = invSin2*(sinBeta       -sinAlphaBeta*cosAlpha);
	A[1][1] = invSin2*(sinAlphaBeta - sinBeta*cosAlpha);
	/* make sure sufficient difference   in angles for 3d solution	*/
	if(fabs(alpha) < 0.8 || fabs(alpha) > (2*PI-0.8)) A[0][0]=-LARGEINT;
}

/*
  Return FALSE if z value not in reasonable range.
 */
static int badZ(double z1) {
	if(z1 <= MINELEVATION || z1 > 10000.0) return TRUE;
	return(FALSE);
}

/*
Avoid overly large slopes
*/
double limitSlope(double slope, double maxSlope) {
	if(fabs(slope) <= maxSlope) return slope;
	return copysign(maxSlope, slope);
}

/*
  Computer B matrix for 3D vel solution
 */
void computeB(double x, double y, double z, double B[2][2], double *dzdx, double *dzdy, double aPsi, double dPsi,xyDEM *xydem)
{
	double  zx1,zx2,zy1,zy2;
	double x1,y1,x2,y2;
	double  dx,dy;
	double tanAPsi, tanDPsi;
	xyDEM dem;
	dem = *xydem;
	/* Use 90 meter spacing */
	dx = max(dem.deltaX,0.09);
	dy = max(dem.deltaY,0.09);
	/*  Compute gradient	*/
	x1 = x + 0.5*dx;	x2 = x - 0.5*dx;
	y1 = y + 0.5*dy;	y2 = y - 0.5*dy;
	zx1 = interpXYDEM(x1,y,dem);
	zx2 = interpXYDEM(x2,y,dem);
	zy1 = interpXYDEM(x,y1,dem);
	zy2 = interpXYDEM(x,y2,dem);
	/* Zero slopes if bad data */		
	if( badZ(zx1)==TRUE || badZ(zx2) == TRUE || badZ(zy1)==TRUE || badZ(zy2) == TRUE) {
		*dzdx = 0.0; *dzdy = 0.0;
	} else {
		*dzdx = limitSlope((zx1-zx2) / (KMTOM * dx), 0.1);
		*dzdy = limitSlope((zy1-zy2) / (KMTOM * dy), 0.1);
	}
	/* Form B matrix */
	tanAPsi = tan(aPsi);
	tanDPsi = tan(dPsi);
	B[0][0] = *dzdx/tanAPsi;
	B[0][1] = *dzdy/tanAPsi;
	B[1][0] = *dzdx/tanDPsi;
	B[1][1] = *dzdy/tanDPsi;
	return;
}

/*
  Compute 3D velocities
*/
void computeVxy(double aP, double dP,double aPe, double dPe,double A[2][2], double B[2][2], double *vx, double *vy, double *scaleX,double *scaleY)
{
	double C[2][2], Cinv[2][2], D[2][2];
	double detC;
	double ex,ey;
	/*
	  C = I - AB
	*/
	C[0][0] = 1.0 - (A[0][0] * B[0][0] + A[0][1] * B[1][0]);
	C[0][1] =     - (A[0][0] * B[0][1] + A[0][1] * B[1][1]);
	C[1][0] =     - (A[1][0] * B[0][0] + A[1][1] * B[1][0]);
	C[1][1] = 1.0 - (A[1][0] * B[0][1] + A[1][1] * B[1][1]);
	/*
	  Invert C
	*/
	detC = C[0][0] * C[1][1] - C[0][1]*C[1][0];
	Cinv[0][0] = C[1][1]/detC;
	Cinv[1][1] = C[0][0]/detC;
	Cinv[0][1] = -C[0][1]/detC;
	Cinv[1][0] = -C[1][0]/detC;
	/*
	  (I-AB)^-1 * A
	*/        
	D[0][0] = Cinv[0][0] * A[0][0] + Cinv[0][1] * A[1][0];
	D[0][1] = Cinv[0][0] * A[0][1] + Cinv[0][1] * A[1][1];
	D[1][0] = Cinv[1][0] * A[0][0] + Cinv[1][1] * A[1][0];
	D[1][1] = Cinv[1][0] * A[0][1] + Cinv[1][1] * A[1][1];
	/* 
	   Vxy = D * [px,py] (column)
	*/
	*vx = D[0][0] * aP + D[0][1] * dP;
	*vy = D[1][0] * aP + D[1][1] * dP;
	/*
	  Compute errors
	  This computes the diagnal elements of the covariance matrix
	  (D Cxx)D'
	  Note these are variances
	*/
	ex = (D[0][0]*D[0][0]*aPe*aPe + D[0][1]*D[0][1]*dPe*dPe);
	ey = (D[1][0]*D[1][0]*aPe*aPe+  D[1][1]*D[1][1]*dPe*dPe);
	/*
	  Compute scale.
	*/
	*scaleX = 1.0/ex;
	*scaleY = 1.0/ey;
	/* 
	   Added this to avoid large errors with wonky slopes 
	   Forces no solution for low detC
	*/
	if( detC < 0.25) {
		*vx=-LARGEINT;	*vy=-LARGEINT;
		*scaleX=0.0;   	*scaleY=0.0;
	} 
	return;
}

/*
   Convert WGS84 z at lat, to sperical with Re as reference
 */
double sphericalElev(double z,double lat,double Re)
{
	z +=  (earthRadius(lat*DTOR,EMINOR,EMAJOR)*KMTOM - Re);
	return z;
}

/*
  Convert spherical elevation at lat with radius Re to WGS
 */
double sphericalToWGSElev(double z,double lat,double Re)
{
	z +=  (Re-earthRadius(lat*DTOR,EMINOR,EMAJOR)*KMTOM);
	return z;
}

/*
   Evaluate quadratic baseline polynomial
 */
double bPoly(double b0,double b1,double b2,double x) {
	return  (b0 + b1 * x +  b2 * x*x);
}

/*
  Rotate flow direction vector to xy (PS)  coordinates
*/
void errorsToXY(double er,double ea, double *ex, double *ey, double xyAngle, double hAngle)
{
	double rotAngle,cosRot,sinRot;
	/*
	  Compute angle for CW rotation of ps to radar
	*/
	rotAngle = hAngle - xyAngle;
	/*
	  Compute range displacement vyra by ccw rot from ps to ra
	*/      
	sinRot =  sin(rotAngle);
	cosRot =  cos(rotAngle);
	*ex =  (er * er) * (cosRot * cosRot) + (ea * ea) * (sinRot * sinRot);
	*ey =  (er * er) * (sinRot  * sinRot) +  (ea * ea) * (cosRot * cosRot);                  
}

/*
  Returns positive value if a point (x0,y0) is to the left of the line formed by two points (x1,y1) and(x2,y2)
*/
static float isLeft(double x0, double y0,double x1,double y1,double x2,double y2){
	return (x1-x0) * (y2-y0) - (x2 -x0) * (y1-y0);
}

/*
	Check if a point is in a closed rect with ccw points.
*/
static int inRect(double xP, double yP,double *xR,double *yR){
	int i;
	/* Any point not left indicates point outside rectangle for CCW check */
	for(i=0; i < 4; i++) if(isLeft(xP, yP, xR[i], yR[i], xR[i+1], yR[i+1]) < 0) return(FALSE);
	return TRUE;
}

static void xyMinMax(double *xR, double *yR, double *minX, double *minY, double *maxX, double *maxY){
	int i;
	*minX = 1e20; *minY = 1e20; *maxX = -1e20; *maxY = -1e20;
	for(i=0; i < 4; i++) {
		*minX = min(xR[i], *minX);  *minY = min(yR[i], *minY);
		*maxX = max(xR[i], *maxX);  *maxY = max(yR[i], *maxY);
	}
}

/*
  Compute bounding box for area of intersection of and ascending and descending pass
 */
void getIntersect(inputImageStructure *dPhaseImage,inputImageStructure *aPhaseImage,
		  int *iMin,int *iMax,int *jMin,int *jMax, outputImageStructure *outputImage) 
{
	extern double Rotation;
	double minX,maxX,minY,maxY;
	double minXa,maxXa,minYa,maxYa;
	double minXd,maxXd,minYd,maxYd;
	double pad;
	double xaP[6],yaP[6];
	double xdP[6],ydP[6];
	int inA,inD, intersect;
	int i,j;
	*iMin=0; *iMax=0; *jMin=0; *jMax=0;
	/*
	 Compute xy coords of rectangles, with CCW order. Tack center point on the end.
	*/
    lltoxy1(aPhaseImage->latControlPoints[1],aPhaseImage->lonControlPoints[1],xaP,yaP,Rotation,outputImage->slat);
	lltoxy1(aPhaseImage->latControlPoints[2],aPhaseImage->lonControlPoints[2],xaP+1,yaP+1,Rotation,outputImage->slat);
	lltoxy1(aPhaseImage->latControlPoints[4],aPhaseImage->lonControlPoints[4],xaP+2,yaP+2,Rotation,outputImage->slat);
	lltoxy1(aPhaseImage->latControlPoints[3],aPhaseImage->lonControlPoints[3],xaP+3,yaP+3,Rotation,outputImage->slat);
	lltoxy1(aPhaseImage->latControlPoints[0],aPhaseImage->lonControlPoints[0],xaP+5,yaP+5,Rotation,outputImage->slat);
	xaP[4] = xaP[0]; yaP[4] = yaP[0];  /* Complete polygon */
	lltoxy1(dPhaseImage->latControlPoints[1],dPhaseImage->lonControlPoints[1],xdP,ydP,Rotation,outputImage->slat);
	lltoxy1(dPhaseImage->latControlPoints[2],dPhaseImage->lonControlPoints[2],xdP+1,ydP+1,Rotation,outputImage->slat);
	lltoxy1(dPhaseImage->latControlPoints[4],dPhaseImage->lonControlPoints[4],xdP+2,ydP+2,Rotation,outputImage->slat);
	lltoxy1(dPhaseImage->latControlPoints[3],dPhaseImage->lonControlPoints[3],xdP+3,ydP+3,Rotation,outputImage->slat);
	lltoxy1(dPhaseImage->latControlPoints[0],dPhaseImage->lonControlPoints[0],xdP+5,ydP+5,Rotation,outputImage->slat);
	xdP[4] = xdP[0]; ydP[4] = ydP[0];
	intersect = FALSE;
	/* If any point from one rect falls in the other they intersect */
	minX = 1e20; minY = 1e20; maxX = -1e20; maxY = -1e20;
	for(i=0; i < 6; i++) { /* This checks corners and centers */
		inA = inRect(xdP[i], ydP[i], xaP, yaP);
		inD = inRect(xaP[i], yaP[i], xdP, ydP);
		/*fprintf(stderr, "%i %i %f %f %f %f %f\n", inA, inD, xdP[i],ydP[i], xaP[i],yaP[i], dPhaseImage->latControlPoints[i]);*/
		if (inA == TRUE || inD  == TRUE) {
			intersect = TRUE;
		}
	}
	/* Compute xy bounds */
	xyMinMax(xdP, ydP, &minXd, &minYd, &maxXd, &maxYd);
	xyMinMax(xaP, yaP, &minXa, &minYa, &maxXa, &maxYa);
	minX = max(minXa, minXd); /* Take inner bounds to restrict to just overlap */
	maxX = min(maxXa, maxXd);
	minY = max(minYa, minYd); 
	maxY = min(maxYa, maxYd);
	/* Compute image bounds */
	if(intersect == FALSE) { *iMin=0; *iMax=0; *jMin=0; *jMax=0; return;}
	pad=15000.; 
	*iMin=(int)((minY*KMTOM-outputImage->originY - pad)/outputImage->deltaY);
	*jMin=(int)((minX*KMTOM-outputImage->originX - pad)/outputImage->deltaX);
	*iMax=(int)((maxY*KMTOM-outputImage->originY + pad)/outputImage->deltaY);
	*jMax=(int)((maxX*KMTOM-outputImage->originX + pad)/outputImage->deltaX);
	*iMin=max(*iMin,0); *jMin=max(*jMin,0);
	*iMax=min(outputImage->ySize,*iMax); *jMax=min(outputImage->xSize,*jMax);
	return;
}



/*
  OBSOLETE REMOVE AT SOME POINT
 */
void getIntersectOld(inputImageStructure *dPhaseImage,inputImageStructure *aPhaseImage,
		  int *iMin,int *iMax,int *jMin,int *jMax, outputImageStructure *outputImage) 
{
	extern double Rotation;
	double a,b,c,d;
	double xa1,xa2,ya1,ya2;
	double xd1,xd2,yd1,yd2;
	double xi,yi;
	double minX,maxX,minY,maxY;
	double minXa,maxXa,minYa,maxYa;
	double minXd,maxXd,minYd,maxYd;
	double pad;
	double xiP[4],yiP[4];
	int inA,inD;
	int i,j;
	*iMin=0; *iMax=0; *jMin=0; *jMax=0;
	/*
	  Loop through asc/des and find intersections of image lines along track
	*/    
	for(i=0; i < 2; i++) {
		lltoxy1(aPhaseImage->latControlPoints[1+i],aPhaseImage->lonControlPoints[1+i],&xa1,&ya1,Rotation,outputImage->slat);
		lltoxy1(aPhaseImage->latControlPoints[3+i],aPhaseImage->lonControlPoints[3+i],&xa2,&ya2,Rotation,outputImage->slat);
		/*
		  Compute line params
		*/
		a=(ya2-ya1)/(xa2-xa1);
		b=ya1-a*xa1;
		/*
		  Save min/max to define image bounding box.
		*/
		if(i==0) {
			minXa=min(xa1,xa2); minYa=min(ya1,ya2);
			maxXa=max(xa1,xa2); maxYa=max(ya1,ya2);
		} else {
			minXa=min(minXa,min(xa1,xa2)); minYa=min(minYa,min(ya1,ya2));
			maxXa=max(maxXa,max(xa1,xa2)); maxYa=max(maxYa,max(ya1,ya2));
		}   

		for(j=0; j < 2; j++) {
			lltoxy1(dPhaseImage->latControlPoints[1+j],dPhaseImage->lonControlPoints[1+j],&xd1,&yd1,Rotation,outputImage->slat);
			lltoxy1(dPhaseImage->latControlPoints[3+j],dPhaseImage->lonControlPoints[3+j],&xd2,&yd2,Rotation,outputImage->slat);
			/*
			  Compute line params
			*/
			c=(yd2-yd1)/(xd2-xd1);
			/*           
				     Check for parallel tracks. Added 05/31/07
			*/
			if(fabs(a-c) < .2) { *iMin=0; *iMax=0; *jMin=0; *jMax=0; return;}
			d=yd1-c*xd1;
			/*
			  Compute intersection
			*/
			xi = (b-d)/(c-a);
			yi = (b*c - a*d)/(c-a);
			/*
			  Compute bounding box of intersecting points
			*/
			if(i==0 && j==0) {minX=xi; minY=yi; maxX=xi; maxY=yi;}
			else {
				minX=min(xi,minX);  minY=min(yi,minY);
				maxX=max(xi,maxX);  maxY=max(yi,maxY);
			}
			/*
			  Save min/max to define image bounding box.
			*/
			if(j==0) {
				minXd=min(xd1,xd2); minYd=min(yd1,yd2);
				maxXd=max(xd1,xd2); maxYd=max(yd1,yd2);
			} else {
				minXd=min(minXd,min(xd1,xd2)); minYd=min(minYd,min(yd1,yd2));
				maxXd=max(maxXd,max(xd1,xd2)); maxYd=max(maxYd,max(yd1,yd2));
			}   
			/*
			  Save intersection 
			*/
			yiP[i*2+j]=yi;
			xiP[i*2+j]=xi;
		}  /* End j */
	}
	/*
	  Compute i,j min,max with pad
	*/ 
	pad=15000.; 
	*iMin=(int)((minY*KMTOM-outputImage->originY - pad)/outputImage->deltaY);
	*jMin=(int)((minX*KMTOM-outputImage->originX - pad)/outputImage->deltaX);
	*iMax=(int)((maxY*KMTOM-outputImage->originY + pad)/outputImage->deltaY);
	*jMax=(int)((maxX*KMTOM-outputImage->originX + pad)/outputImage->deltaX);
	*iMin=max(*iMin,0); *jMin=max(*jMin,0);
	*iMax=min(outputImage->ySize,*iMax); *jMax=min(outputImage->xSize,*jMax);
	/*
	  Final check. At least one intersection point must fall in bounding box for 
	  each image.
	*/
	inA=FALSE; inD=FALSE;
	for(i=0; i < 4; i++) {
		if(xiP[i] >= minXa && xiP[i]< maxXa && yiP[i] >=minYa && yiP[i] < maxYa)
			inA=TRUE;
		if(xiP[i] >= minXd && xiP[i]< maxXd && yiP[i] >=minYd && yiP[i] < maxYd)
			inD=TRUE;
		fprintf(stderr,"**** %f %f %i %i\n", xiP[i], yiP[i], inA, inD);
	}
	return;
	if(inA==FALSE && inD == FALSE) { *iMin=0; *iMax=0; *jMin=0; *jMax=0;}
}

/*
  Computes the overall heading angles for two scenes and checks if the angular diffrence sufficient to continue. 
  if it returns *iMax=jMax=0 then it will skip this pair. 
*/
void computeSceneAlpha(outputImageStructure *outputImage,	inputImageStructure *aOffImage, inputImageStructure *dOffImage,
		       conversionDataStructure *aCp,  conversionDataStructure *dCp, xyDEM *dem,  int *iMin,int *iMax,int *jMin,int *jMax) {
	extern int HemiSphere;
	extern double Rotation;	
	double aHAngle,dHAngle,x,y;
	double aTmp,dTmp,alpha;
	double lat,lon;
	double aAzimuth, dAzimuth,drange,arange;
	int i,j;
	aHAngle=9999.0; dHAngle=9999.0;
	if(*iMax > 0 && *jMax > 0) {
		for(i=*iMin; i < *iMax; i+=max((*iMax-*iMin)/4,1)) {
			for(j=*jMin; j < *jMax; j+=max((*jMax-*jMin)/4,1)) {
				y = (outputImage->originY + i*outputImage->deltaY) * MTOKM;
				x = (outputImage->originX + j*outputImage->deltaX) * MTOKM;
				xytoll1(x,y,HemiSphere,&lat,&lon,Rotation,dem->stdLat);
				llToImageNew(lat,lon,0,&arange,&aAzimuth,aOffImage);
				aTmp = computeHeading(lat,lon,0.0,aOffImage,aCp);
				llToImageNew(lat,lon,0,&drange,&dAzimuth,dOffImage);
				dTmp = computeHeading(lat,lon,0.0,dOffImage,dCp);
				if(aTmp < 999 && dTmp < 999) { aHAngle=aTmp; dHAngle=dTmp; }
			}
		}
		alpha=0.0;
		if(aHAngle < 999 && dHAngle < 999) alpha = aHAngle - dHAngle;
		/* Changed from 0.45 to 0.6 06/07/07 */
		if(fabs(alpha) < 0.7 || fabs(alpha) > (2*PI-0.7)) {
			*iMax=0; *jMax=0;
			fprintf(stderr,"---No solution for headings %f %f %f\n",aHAngle*RTOD,dHAngle*RTOD,alpha*RTOD);
		} /* else fprintf(stderr,"*** Solution for headings %f %f %f \n",aHAngle*RTOD,dHAngle*RTOD,alpha*RTOD); */
	}
}

/*
  Malloc a floating point image
 */
float **mallocImage(int nr,int nc)
{
	float **tmp,*tmp1;
	int i;

	tmp = (float **) malloc((size_t)(nr*sizeof(float *)));

	for(i=0; i < nr; i++ ) {
		tmp1=(float *) malloc((size_t)(nc*sizeof(float)));
		if(tmp1 == NULL) 
			error("mallocImage(getGeoCodeInputImage.c) - unable to malloc data for %i x %i %i\n",nr,nc,nc*sizeof(float));
		tmp[i] = tmp1;
	}
	return tmp;
}


void setTiePointsMapProjectionForHemisphere(tiePointsStructure *tiePoints)
{
	extern int HemiSphere;
	extern double Rotation;	

	if(tiePoints->lat[0] < 0) {
		fprintf(stderr,"**** SOUTHERN HEMISPHERE ****");
		HemiSphere=SOUTH;
		tiePoints->stdLat=71;
		Rotation=0.0;
	} else{
		fprintf(stderr,"**** NORTHERN HEMISPHERE ****");
		HemiSphere=NORTH;
		tiePoints->stdLat=70;
		Rotation=45.0;
	}
	tiePoints->imageCoords = FALSE;	
}

/* 
append a suffix to a file name, use buff for space. 
*/
char *appendSuffix(char *file, char *suffix,char *buf ) {
	char *datFile;
	datFile=&(buf[0]);
	buf[0]='\0';
	datFile= strcat(datFile,file);   
	datFile = strcat(datFile,suffix);
	return datFile;
}

double secondForSAR(SARData *par) {
	return par->hr*3600.0 + par->min*60. + par->sec;
}
double getReH( conversionDataStructure *cP, inputImageStructure *inputImage,double azimuth) {
	double ReH;
	int index;
	/* for  index inbounds */
	index=(int)min(max(0,azimuth),cP->azSize -1);
	if(cP->ReH != NULL ) ReH = cP->ReH[index];
	else { ReH = cP->Re+inputImage->par.H; fprintf(stderr,"*** Warning using fixed ReH ***");} /* Should be obsolete*/
		
	return(ReH);
}
/*
 Compute rho given R, Re+z,Re+H
 */
double rhoRReZReH(double R,double ReZ, double ReH) {
	return acos( (R*R - ReH * ReH - ReZ*ReZ)/ (-2.0 * ReZ * ReH) );
}

/*
 Compute look angle theta  given R, Re+z,Re+H
 */
double thetaRReZReH(double R,double ReZ, double ReH) {
	return acos( (R*R + ReH*ReH - ReZ*ReZ) / (2.0*R*ReH)  );
}

/*
 Compute incidence angle psi  given R, Re+z,Re+H
 */
double psiRReZReH(double R,double ReZ, double ReH) {
	double theta;
	theta = thetaRReZReH( R,ReZ,ReH );
	return asin((ReH/ReZ)*sin(theta));
}

/*
  Compute slant range give rho, Re+z, and Re+H
 */
double slantRange(double rho,double ReZ, double ReH) {
	return sqrt( ReH*ReH + ReZ*ReZ - 2.0*ReZ*ReH*cos(rho) );
}


/*
  Init output image
*/
void initOutputImage(outputImageStructure *outputImage, inputImageStructure inputImage)
{    
	extern double Rotation;
	double pixelSize;
	double x,y;
	double minX=FLT_MAX,maxX=-FLT_MAX;
	double minY=FLT_MAX,maxY=-FLT_MAX;
	int i;
	fprintf(stderr,"InitOutput Image %f\n",Rotation);
	/*
	  Determine pixels Size
	*/
	pixelSize = AzimuthPixelSize * inputImage.nAzimuthLooks;
	if(outputImage->deltaX > 0.0) pixelSize = outputImage->deltaX;
	outputImage->deltaX = pixelSize;
	outputImage->deltaY = pixelSize;
	fprintf(stderr,"pixelSize = %f\n",pixelSize);
	fprintf(stderr,"--x0,y0 %lf %lf\n\n",outputImage->originX,outputImage->originY);
	/*
	  Determin lower left and upper right corner of image
	*/
	for(i = 1; i < 5; i++ ) {
		lltoxy(inputImage.latControlPoints[i],inputImage.lonControlPoints[i],  &x, &y, Rotation);
		fprintf(stderr," %f %f %f %f\n",  x,y,inputImage.latControlPoints[i],inputImage.lonControlPoints[i]);
		x = x * KMTOM; /* Conversion from km to m */
		y = y * KMTOM; /* Conversion from km to m */
		minX = min(minX,x);
		minY = min(minY,y);
		maxX = max(maxX,x);
		maxY = max(maxY,y);
	}
	if(outputImage->originX > (LARGEINT-1))  outputImage->originX = minX; /* else outputImage->originX=0.0;*/
	if(outputImage->originY > (LARGEINT-1))   outputImage->originY = minY; /* else outputImage->originY=0.0;*/
	fprintf(stderr,"--x0,y0 %lf %lf\n\n",outputImage->originX,outputImage->originY);
	/*
	  Determine image type float/complex
	*/
	if(inputImage.imageType == COMPLEX) outputImage->imageType = COMPLEX;
	else outputImage->imageType = FLOAT;
	/*
	  Determine pixel size of image
	*/    
	if(outputImage->xSize <= 0.0)
		outputImage->xSize = (int)((maxX-minX)/outputImage->deltaX) + 1;
	if(outputImage->ySize <= 0.0)
		outputImage->ySize = (int)((maxY-minY)/outputImage->deltaY) + 1;
	fprintf(stderr,"xSize,ySize = %i %i\n",outputImage->xSize,outputImage->ySize);
	/*
	  Init space for image 
	*/
	if(outputImage->noMem == FALSE) {
		outputImage->image=(void **)   mallocImage(outputImage->ySize,outputImage->xSize);
	}

	return;
}
