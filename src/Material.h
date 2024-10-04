#ifndef MATERIAL_H
#define MATERIAL_H

#include "imageLoader.h"
#include "Vec3.h"
#include <cmath>

#include <GL/glut.h>

enum MaterialType {
    Material_Diffuse_Blinn_Phong ,
    Material_Glass,
    Material_Mirror
};


struct Material {
    Vec3 ambient_material;
    Vec3 diffuse_material;
    Vec3 specular_material;
    double shininess;

    float index_medium;
    float transparency;
    
    bool isTextured; // true si l'objet aura une texture
    ppmLoader::ImageRGB texture; // Si l'objet a besoin d'une texture, elle sera spécifiée ici

    MaterialType type;

    Material() {
        type = Material_Diffuse_Blinn_Phong;
        transparency = 0.0; // Si différent de 0.0 alors le matériaux est transparent
        index_medium = 1.0;
        ambient_material = Vec3(0.2, 0.2, 0.2);
        isTextured = false; // Par défaut pas de texture
    }
    
};



#endif // MATERIAL_H
