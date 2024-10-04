#ifndef TRIANGLE_H
#define TRIANGLE_H
#include "Vec3.h"
#include "Ray.h"
#include "Plane.h"
#include <cfloat>

struct RayTriangleIntersection{
    bool intersectionExists;
    float t;
    float w0,w1,w2;
    unsigned int tIndex;
    Vec3 intersection;
    Vec3 normal;
};

class Triangle {
private:
    Vec3 m_normal;
    float area;
public:
    Triangle() {}
    Triangle( Vec3 const & c0 , Vec3 const & c1 , Vec3 const & c2 ) {
        m_c[0] = c0;
        m_c[1] = c1;
        m_c[2] = c2;
        updateAreaAndNormal();
    }
    
    Vec3 m_c[3];
    
    void updateAreaAndNormal() {
        Vec3 nNotNormalized = Vec3::cross( m_c[1] - m_c[0] , m_c[2] - m_c[0] );
        float norm = nNotNormalized.length();
        m_normal = nNotNormalized / norm;
        area = norm / 2.f;
    }
    void setC0( Vec3 const & c0 ) { m_c[0] = c0; } // remember to update the area and normal afterwards!
    void setC1( Vec3 const & c1 ) { m_c[1] = c1; } // remember to update the area and normal afterwards!
    void setC2( Vec3 const & c2 ) { m_c[2] = c2; } // remember to update the area and normal afterwards!
    Vec3 const & normal() const { return m_normal; }
    Vec3 projectOnSupportPlane( Vec3 const & p ) const {
        Vec3 result;
        //TODO completer
        return result;
    }
    float squareDistanceToSupportPlane( Vec3 const & p ) const {
        float result;
        //TODO completer
        return result;
    }
    float distanceToSupportPlane( Vec3 const & p ) const { return sqrt( squareDistanceToSupportPlane(p) ); }
    bool isParallelTo( Line const & L ) const {
        bool result;
        //TODO completer
        return result;
    }
    Vec3 getIntersectionPointWithSupportPlane( Line const & L ) const {
        // you should check first that the line is not parallel to the plane!
        Vec3 result;
        //TODO completer
        return result;
    }
    void computeBarycentricCoordinates( Vec3 const & p , float & u0 , float & u1 , float & u2 ) const {
        //TODO Complete
    }

    RayTriangleIntersection intersect( Ray const & ray ) const {
        RayTriangleIntersection intersection;
        intersection.intersectionExists = false;
        
        Vec3 d = ray.direction(); // Vecteur directeur du rayon
    	Vec3 o = ray.origin(); // Origine du rayon
    	Vec3 n = this->m_normal; // Normale à la surface
    	n.normalize(); // Ne pas oublier de normaliser n
    	
    	Vec3 a = m_c[0];
    	Vec3 r_vector = m_c[1] - m_c[0];
    	Vec3 u_vector = m_c[2] - m_c[0];
    	
    	
    	float D = Vec3::dot(a,n);
    	
    	float denominateur = Vec3::dot(d,n);
    	
    	if (std::abs(denominateur) < 2*FLT_MIN){ // Dans ce cas, rayon parallèle au plan auquel appartient le triangle donc ne pas chercher d'intersection (2 fois FLT_MIN à cause des imprécisions numériques éventuelles)
    		return intersection;
    	}
    	
    	float t = (D - Vec3::dot(o,n)) / denominateur;
    	if (t >= 0) {
            Vec3 prime = o + t * d;
            
            // On vérifie si le rayon est à l'intérieur du triangle
            Vec3 p0 = this->m_c[0];
            Vec3 p1 = this->m_c[1];
            Vec3 p2 = this->m_c[2];
            
            Triangle t1 = Triangle(p0,p1,prime);
            Triangle t2 = Triangle(p1,p2,prime);
            Triangle t3 = Triangle(p2,p0,prime);
            float w0 = t1.area / this->area;
            float w1 = t2.area / this->area;
            float w2 = t3.area / this->area;
            float total = w0 + w1 + w2;
            if (total >= 1. - 1e-4 && total <= 1 + 1e-4){
		        intersection.intersectionExists = true;
				intersection.intersection = prime;
				intersection.t = t;
				intersection.normal = n;
				intersection.w0 = w0;
				intersection.w1 = w1;
				intersection.w2 = w2;
    			//unsigned int tIndex;
		    }
        }
        // 1) check that the ray is not parallel to the triangle:

        // 2) check that the triangle is "in front of" the ray:

        // 3) check that the intersection point is inside the triangle:
        // CONVENTION: compute u,v such that p = w0*c0 + w1*c1 + w2*c2, check that 0 <= w0,w1,w2 <= 1

        // 4) Finally, if all conditions were met, then there is an intersection! :

        return intersection;
    }
};
#endif
