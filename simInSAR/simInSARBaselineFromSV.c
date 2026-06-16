#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "mosaicSource/common/common.h"
#include "simInSARInclude.h"

/*
  Allocate and fill scene->bnArray and scene->bpArray with per-azimuth-line
  baseline components derived from the state vectors of the primary scene
  (scene->I) and a secondary geodat file.  scene->aSize must already be set
  (i.e. call after parseSceneFile).
*/
void simInSARBaselineFromSV(char *geodat2File, sceneStructure *scene)
{
	inputImageStructure inputImage2;
	double Re, ReH, thetaC, dt1t2;
	double azTime, bn, bp;
	int32_t i;

	parseInputFile(geodat2File, &inputImage2);
	initllToImageNew(&inputImage2);

	Re     = scene->I.cpAll.Re;
	ReH    = getReH(&(scene->I.cpAll), &(scene->I), scene->I.azimuthSize / 2);
	thetaC = thetaRReZReH(scene->I.cpAll.RCenter, Re, ReH);
	dt1t2  = scene->I.cpAll.sTime - inputImage2.cpAll.sTime;

	scene->bnArray = (float *)malloc(scene->aSize * sizeof(float));
	scene->bpArray = (float *)malloc(scene->aSize * sizeof(float));
	if (scene->bnArray == NULL || scene->bpArray == NULL)
		error("simInSARBaselineFromSV: malloc failed\n");

	for (i = 0; i < scene->aSize; i++)
	{
		/* Convert output azimuth index to SAR single-look time */
		azTime = scene->I.cpAll.sTime +
		         (scene->aO + i * scene->dA) * scene->I.nAzimuthLooks / scene->I.par.prf;
		svBnBp(azTime, thetaC, dt1t2,
		       &(scene->I.sv), &inputImage2.sv,
		       &bn, &bp, scene->I.lookDir);
		scene->bnArray[i] = (float)bn;
		scene->bpArray[i] = (float)bp;
	}
	fprintf(stderr, "SV baseline: bn[0]=%f bp[0]=%f  bn[end]=%f bp[end]=%f\n",
	        scene->bnArray[0], scene->bpArray[0],
	        scene->bnArray[scene->aSize - 1], scene->bpArray[scene->aSize - 1]);
}
