#ifndef SQUARE_H
#define SQUARE_H
#include "Vec3.h"
#include <vector>
#include "Mesh.h"
#include <cmath>

struct RaySquareIntersection{
    bool intersectionExists;
    float t;
    float u,v;
    Vec3 intersection;
    Vec3 normal;
};


class Square : public Mesh {
public:
    Vec3 m_normal;
    Vec3 m_bottom_left;
    Vec3 m_right_vector;
    Vec3 m_up_vector;
    ppmLoader::ImageRGB squareTexture;

    Square() : Mesh() {}
    Square(Vec3 const & bottomLeft , Vec3 const & rightVector , Vec3 const & upVector , float width=1. , float height=1. ,
           float uMin = 0.f , float uMax = 1.f , float vMin = 0.f , float vMax = 1.f) : Mesh() {
        setQuad(bottomLeft, rightVector, upVector, width, height, uMin, uMax, vMin, vMax);
    }

    void setQuad( Vec3 const & bottomLeft , Vec3 const & rightVector , Vec3 const & upVector , float width=1. , float height=1. ,
                  float uMin = 0.f , float uMax = 1.f , float vMin = 0.f , float vMax = 1.f) {
        m_right_vector = rightVector;
        m_up_vector = upVector;
        m_normal = Vec3::cross(rightVector , upVector);
        m_bottom_left = bottomLeft;

        m_normal.normalize();
        m_right_vector.normalize();
        m_up_vector.normalize();

        m_right_vector = m_right_vector*width;
        m_up_vector = m_up_vector*height;

        vertices.clear();
        vertices.resize(4);
        vertices[0].position = bottomLeft;                                      vertices[0].u = uMin; vertices[0].v = vMin;
        vertices[1].position = bottomLeft + m_right_vector;                     vertices[1].u = uMax; vertices[1].v = vMin;
        vertices[2].position = bottomLeft + m_right_vector + m_up_vector;       vertices[2].u = uMax; vertices[2].v = vMax;
        vertices[3].position = bottomLeft + m_up_vector;                        vertices[3].u = uMin; vertices[3].v = vMax;
        vertices[0].normal = vertices[1].normal = vertices[2].normal = vertices[3].normal = m_normal;
        triangles.clear();
        triangles.resize(2);
        triangles[0][0] = 0;
        triangles[0][1] = 1;
        triangles[0][2] = 2;
        triangles[1][0] = 0;
        triangles[1][1] = 2;
        triangles[1][2] = 3;


    }
	
	
    RaySquareIntersection intersect(const Ray &ray) const {
        RaySquareIntersection intersection;
        intersection.intersectionExists = false;
        
        // En suivant la formule du cours
    	Vec3 d = ray.direction(); // Vecteur directeur du rayon
    	Vec3 o = ray.origin(); // Origine du rayon
    	Vec3 n = this->m_normal; // Normale à la surface
    	n.normalize(); // Ne pas oublier de normaliser n
    	
    	Vec3 a = this->m_bottom_left;
    	Vec3 r_vector = this->m_right_vector;
    	Vec3 u_vector = this->m_up_vector;
    	
    	
    	float D = Vec3::dot(a,n);
    	
    	float denominateur = Vec3::dot(d,n);
    	
    	if (std::abs(denominateur) < 2*FLT_MIN){ // Dans ce cas, rayon parallèle au plan donc ne pas chercher d'intersection (2 fois FLT_MIN à cause des imprécisions numériques éventuelles)
    		return intersection;
    	}
    	
    	float t = (D - Vec3::dot(o,n)) / denominateur;
    	if (t >= 0) {
            Vec3 pointIntersection = o + t * d;

            // u et v sur le plan
            Vec3 tr = pointIntersection - a;
            float u = Vec3::dot(tr, r_vector) / r_vector.squareLength();
            float v = Vec3::dot(tr, u_vector) / u_vector.squareLength();

            // Si l'intersection est dans les limites du plan
            if (u >= 0 && u <= 1 && v >= 0 && v <= 1) {
		        intersection.intersectionExists = true;
		        intersection.intersection = pointIntersection;
		        intersection.t = t;
		        intersection.normal = n; // Normale au point d'un plan est la normale du plan
		        // n doit bien être normalisé (c'est le cas un peu plus haut)
		        
		        intersection.u = u;
		        intersection.v = v;
            }
        }
    	
    	/* Ancienne version qui ne fonctionnait pas pour les plans alignés sur l'axe z
    	intersection.intersection = o + t * d;
        Vec3 top_right = m_bottom_left + m_right_vector + m_up_vector;
       	
		if(intersection.intersection[0] < m_bottom_left[0] ||//en dehors du square
			intersection.intersection[1] < m_bottom_left[1] ||
			//intersection.intersection[2] < m_bottom_left[2] || //Inutile de regarder pour z apparemment, c'est ça qui causait ce problème de texture bizarre que j'obtennais pour mes squares
			intersection.intersection[0] > top_right[0] ||
			intersection.intersection[1] > top_right[1]){ // ||
			//intersection.intersection[2] > top_right[2]){
				return intersection;
        	}
        
    	
    	if (t > 0){ // Il y a une intersection
    		intersection.intersectionExists = true;
    		intersection.t = t;
    		intersection.normal = n; // La normale à l'intersection est la normale du plan
    		
    		//float u,v;	
    	}
    	*/

        return intersection;
    }
    
};
#endif // SQUARE_H
