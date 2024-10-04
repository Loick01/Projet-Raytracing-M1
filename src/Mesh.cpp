#include "Mesh.h"
#include <iostream>
#include <fstream>

void Mesh::loadOFF (const std::string & filename) {
    std::ifstream in (filename.c_str ());
    if (!in)
        exit (EXIT_FAILURE);
    std::string offString;
    unsigned int sizeV, sizeT, tmp;
    in >> offString >> sizeV >> sizeT >> tmp;
    vertices.resize (sizeV);
    triangles.resize (sizeT);
    
    float pX, pY, pZ, cR, cG, cB;
    for (unsigned int i = 0; i < sizeV; i++){
        in >> pX >> pY >> pZ >> cR >> cG >> cB;
        
        MeshVertex mv;
        mv.position = Vec3(pX,pY,pZ);
        mv.color = Vec3(cR,cG,cB);
        vertices[i] = mv;
    }
    int nbSommet, ind1, ind2, ind3;
    for (unsigned int i = 0; i < sizeT; i++){
        in >> nbSommet >> ind1 >> ind2 >> ind3;
        
        MeshTriangle mt;
        mt.v[0] = ind1;
        mt.v[1] = ind2;
        mt.v[2] = ind3;
        
        for (int k = 0 ; k < 3 ; k++){
        	vertices[mt.v[k]].listeTriangles.push_back(mt); // Ajoute le triangle à ses 3 sommets, ça servira pour le KD Tree
        }
        
       	triangles[i] = mt;
    }
    
    in.close ();
}

void Mesh::recomputeNormals () {
    for (unsigned int i = 0; i < vertices.size (); i++)
        vertices[i].normal = Vec3 (0.0, 0.0, 0.0);
    for (unsigned int i = 0; i < triangles.size (); i++) {
        Vec3 e01 = vertices[triangles[i].v[1]].position -  vertices[triangles[i].v[0]].position;
        Vec3 e02 = vertices[triangles[i].v[2]].position -  vertices[triangles[i].v[0]].position;
        Vec3 n = Vec3::cross (e01, e02);
        n.normalize ();
        for (unsigned int j = 0; j < 3; j++)
            vertices[triangles[i].v[j]].normal += n;
    }
    for (unsigned int i = 0; i < vertices.size (); i++)
        vertices[i].normal.normalize ();
}

void Mesh::centerAndScaleToUnit () {
    Vec3 c(0,0,0);
    for  (unsigned int i = 0; i < vertices.size (); i++)
        c += vertices[i].position;
    c /= vertices.size ();
    float maxD = (vertices[0].position - c).length();
    for (unsigned int i = 0; i < vertices.size (); i++){
        float m = (vertices[i].position - c).length();
        if (m > maxD)
            maxD = m;
    }
    for  (unsigned int i = 0; i < vertices.size (); i++)
        vertices[i].position = (vertices[i].position - c) / maxD;
}
