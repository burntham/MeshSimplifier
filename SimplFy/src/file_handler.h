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
int mask=0;

	template <typename MeshType> int readFile(MeshType &m, char * filePath)
	{
        printf("Reading mesh at: %s \n",filePath);

        int err=vcg::tri::io::Importer<MeshType>::Open(m,filePath,mask);

		if (err)
		{
			printf("Unable to open mesh at %s\n", filePath);
			exit(-1);
		}else
        {
            printf("mesh loaded with %d vertices and %d faces \n\n",m.vn,m.fn);
        }

		return 0;
	}

	template <typename MeshType> int writeMesh(MeshType &m, char * outfilePath)
    {
        printf("Saving Mesh\n");
        vcg::tri::io::ExporterPLY<MeshType>::Save(m,outfilePath,mask);
        printf("Simplfied mesh saved to: %s\n", outfilePath);
		return 0;
	}
}





#endif /* FILEHANDLER_H_ */
