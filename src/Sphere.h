#ifndef Sphere_H
#define Sphere_H
#include "Vec3.h"
#include <vector>
#include "Mesh.h"
#include <cmath>

struct RaySphereIntersection{
    bool intersectionExists;
    float t;
    float theta,phi;
    Vec3 intersection;
    Vec3 secondintersection;
    Vec3 normal;
};

static
Vec3 SphericalCoordinatesToEuclidean( Vec3 ThetaPhiR ) {
    return ThetaPhiR[2] * Vec3( cos(ThetaPhiR[0]) * cos(ThetaPhiR[1]) , sin(ThetaPhiR[0]) * cos(ThetaPhiR[1]) , sin(ThetaPhiR[1]) );
}
static
Vec3 SphericalCoordinatesToEuclidean( float theta , float phi ) {
    return Vec3( cos(theta) * cos(phi) , sin(theta) * cos(phi) , sin(phi) );
}

static
Vec3 EuclideanCoordinatesToSpherical( Vec3 xyz ) {
    float R = xyz.length();
    float phi = asin( xyz[2] / R );
    float theta = atan2( xyz[1] , xyz[0] );
    return Vec3( theta , phi , R );
}



class Sphere : public Mesh {
public:
    Vec3 m_center;
    float m_radius;
    ppmLoader::ImageRGB sphereTexture;

    Sphere() : Mesh() {}
    Sphere(Vec3 c , float r) : Mesh() , m_center(c) , m_radius(r) {}

    void build_arrays(){
        unsigned int nTheta = 20 , nPhi = 20;
        positions_array.resize(3 * nTheta * nPhi );
        normalsArray.resize(3 * nTheta * nPhi );
        uvs_array.resize(2 * nTheta * nPhi );
        for( unsigned int thetaIt = 0 ; thetaIt < nTheta ; ++thetaIt ) {
            float u = (float)(thetaIt) / (float)(nTheta-1);
            float theta = u * 2 * M_PI;
            for( unsigned int phiIt = 0 ; phiIt < nPhi ; ++phiIt ) {
                unsigned int vertexIndex = thetaIt + phiIt * nTheta;
                float v = (float)(phiIt) / (float)(nPhi-1);
                float phi = - M_PI/2.0 + v * M_PI;
                Vec3 xyz = SphericalCoordinatesToEuclidean( theta , phi );
                positions_array[ 3 * vertexIndex + 0 ] = m_center[0] + m_radius * xyz[0];
                positions_array[ 3 * vertexIndex + 1 ] = m_center[1] + m_radius * xyz[1];
                positions_array[ 3 * vertexIndex + 2 ] = m_center[2] + m_radius * xyz[2];
                normalsArray[ 3 * vertexIndex + 0 ] = xyz[0];
                normalsArray[ 3 * vertexIndex + 1 ] = xyz[1];
                normalsArray[ 3 * vertexIndex + 2 ] = xyz[2];
                uvs_array[ 2 * vertexIndex + 0 ] = u;
                uvs_array[ 2 * vertexIndex + 1 ] = v;
            }
        }
        triangles_array.clear();
        for( unsigned int thetaIt = 0 ; thetaIt < nTheta - 1 ; ++thetaIt ) {
            for( unsigned int phiIt = 0 ; phiIt < nPhi - 1 ; ++phiIt ) {
                unsigned int vertexuv = thetaIt + phiIt * nTheta;
                unsigned int vertexUv = thetaIt + 1 + phiIt * nTheta;
                unsigned int vertexuV = thetaIt + (phiIt+1) * nTheta;
                unsigned int vertexUV = thetaIt + 1 + (phiIt+1) * nTheta;
                triangles_array.push_back( vertexuv );
                triangles_array.push_back( vertexUv );
                triangles_array.push_back( vertexUV );
                triangles_array.push_back( vertexuv );
                triangles_array.push_back( vertexUV );
                triangles_array.push_back( vertexuV );
            }
        }
    }


    RaySphereIntersection intersect(const Ray &ray) const {
        RaySphereIntersection intersection;
        intersection.intersectionExists = false;// Par défaut, dans le cas où il n'y a pas d'intersection
    	
    	// En suivant la formule du cours
    	Vec3 d = ray.direction(); // Vecteur directeur du rayon
    	Vec3 o = ray.origin(); // Origine du rayon
    	Vec3 center = this->m_center; // Centre de la sphère dont on cherche l'intersection avec le rayon
    	float r = this->m_radius; // Rayon de la sphère
    	
    	// On cherche le plus petit t >= 0 tel que : t^2 * dd + 2t * d(o-c) + ||o - c||^2 - r^2 = 0
    	// Calcul du discriminant (b^2 - 4ac) :
    	float a = Vec3::dot(d,d);
    	float b = (2 * Vec3::dot(d,(o-center)));
    	float c = ((o-center).squareLength()-r*r);
    	float delta = pow(b,2) - 4 * a * c;
    	if (delta >= 0){ // Il existe au moins une solution
    		float sol;
    		if (delta > 0){ // 2 solutions
    			sol = std::min( ((-1)*b - sqrt(delta)) / (2*a) , ((-1)*b + sqrt(delta)) / (2*a) );
    		}
    	 	if (delta == 0){ // 1 solution
    			sol = (-1)*b / (2*a); // -b / 2a
    		}
    		intersection.intersectionExists = true;
    		intersection.t = sol;
    		intersection.intersection = o + d * sol;
    		intersection.normal = intersection.intersection - center;
    		intersection.normal.normalize(); // Ne pas oublier de normaliser les normales
    		/* Pour l'instant, je ne met rien dans ces variables
    		float theta,phi;
    		Vec3 secondintersection;
    		*/
    	} // Si pas d'intersection, alors intersection.intersectionExists = false
    	
        return intersection;
    }
};
#endif
