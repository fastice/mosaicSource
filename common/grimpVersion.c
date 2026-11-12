#include <stdio.h>
#include "mosaicSource/common/common.h"
#include "mosaicSource/common/grimpVersionInfo.h"

/*
  Software version of a mosaicSource binary, so every product can be traced to the commit
  that made it. The values come from common/grimpVersionInfo.h, which the Makefile regenerates
  on every make (scripts/makeVersionHeader.sh): the VERSION number, the mosaicSource commit and
  the commits of the linked library repos, each marked -dirty.<diff hash> when built from
  uncommitted changes.
*/

/* One-line version, e.g. "mosaicSource 1.0.3 (a60e090; clib 5590ec2, gdalIO c9a8945, ...)" */
const char *mosaicSourceVersionString(void)
{
	static char versionString[1024];
	snprintf(versionString, sizeof(versionString), "mosaicSource %s (%s; %s)", MOSAICSOURCE_VERSION,
			 MOSAICSOURCE_COMMIT, MOSAICSOURCE_LIBCOMMITS);
	return versionString;
}

/* Add the version keys to a product's metadata (written to every GTiff/VRT that gets it). */
void addMosaicSourceVersion(dictNode **metaData, char *program)
{
	insert_node(metaData, "mosaicSourceProgram", program);
	insert_node(metaData, "mosaicSourceVersion", MOSAICSOURCE_VERSION);
	insert_node(metaData, "mosaicSourceCommit", MOSAICSOURCE_COMMIT);
	insert_node(metaData, "mosaicSourceLibCommits", MOSAICSOURCE_LIBCOMMITS);
}
