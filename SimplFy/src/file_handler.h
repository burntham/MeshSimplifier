/*
 * FileHandler.h
 *
 *  Created on: Sep 18, 2014
 *      Author: brunt
 */

#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export_ply.h>
//#include "simplifier/mesh_properties.h"

#ifndef FILEHANDLER_H_
#define FILEHANDLER_H_

namespace brndan022
{
	template <typename MeshType> int readFile(MeshType &m, char * filePath)
	{
		printf("Attempted to read Mesh%s: \n ",filePath);

		int err=vcg::tri::io::Importer<MeshType>::Open(m,filePath);

		if (err)
		{
			printf("Unable to open mesh at %s\n", filePath);
			exit(-1);
		}else
			printf("%s Was read successfully\n", filePath);

		return 0;
	}

	template <typename MeshType> int writeMesh(MeshType &m, char * outfilePath)
	{
		printf("Writing simplified mesh to %s\n", outfilePath);
		vcg::tri::io::ExporterPLY<MeshType>::Save(m,outfilePath);
		return 0;
	}
}





#endif /* FILEHANDLER_H_ */
