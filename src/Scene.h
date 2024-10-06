#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include <string>
#include "Mesh.h"
#include "Sphere.h"
#include "Square.h"
#include "kdNode.h"


#include <GL/glut.h>

// Caractéristiques des sources de lumières (toutes les lumières ont pour l'instant les mêmes caractéristiques)
#define AMBIANT_INTENSITY 1. // Isa dans diapo -> 1. pour Vec3(1.,1.,1.)
#define DIFFUSE_INTENSITY 1. // Isd dans diapo -> 1. pour Vec3(1.,1.,1.)
#define SPECULAR_INTENSITY 1. // Iss dans diapo -> 1. pour Vec3(1.,1.,1.)

#define INDICE_REFRACTION_VIDE 1.

enum LightType {
    LightType_Spherical,
    LightType_Quad
};


struct Light {
    Vec3 material;
    bool isInCamSpace;
    LightType type;

    Vec3 pos;
    float radius;

    Square quad; // J'ai changé le type de Mesh vers Square
    int nb_points; // Nombre de point qui seront considérés (aléatoirement) sur le carré

    float powerCorrection;

    Light() : powerCorrection(1.0) {}

};

struct RaySceneIntersection{
    bool intersectionExists;
    unsigned int typeOfIntersectedObject;
    unsigned int objectIndex;
    float t;
    RayTriangleIntersection rayMeshIntersection;
    RaySphereIntersection raySphereIntersection;
    RaySquareIntersection raySquareIntersection;
    RayTriangleIntersection rayTriangleIntersection;
    MeshTriangle mt;
    RaySceneIntersection() : intersectionExists(false) , t(FLT_MAX) {}
};



class Scene {
    std::vector< Mesh > meshes;
    std::vector< Sphere > spheres;
    std::vector< Square > squares;
    std::vector< Light > lights;
    
    KdNode *kdRacine; // Racine du KD Tree

public:


    Scene() {
    }
    
    void buildKDTree(){
    	if (meshes.size() == 0){
    		return;
    	}
    	
		std::vector<MeshVertex> sommets; // On récupère tous les sommets de la scène pour construire le KD Tree
		
		/* 
		for (int i = 0 ; i < meshes.size() ; i++){
			for (int j = 0 ; j < meshes[i].vertices.size() ; j++){
				sommets.push_back(meshes[i].vertices[j]); // Récupère tous les MeshVertex de tous les sommets de tous les mesh
			}
		}
		*/
		
		for (int j = 0 ; j < meshes[0].vertices.size() ; j++){
				sommets.push_back(meshes[0].vertices[j]); // Récupère tous les MeshVertex de tous les sommets de meshes[0] uniquement
		}
		
		std::vector<std::vector<MeshTriangle>> lesListesTriangles;
		for (int i = 0 ; i < sommets.size() ; i++){ 
			//sommets[i].listeTriangles = meshes[0].vertices[i].listeTriangles;
			lesListesTriangles.push_back(meshes[0].vertices[i].listeTriangles);
		}
		//std::cout << "Nb de sommets = " << sommets.size() << "\n";
		
		Vec3 minPoint = sommets[0].position;
		Vec3 maxPoint = sommets[0].position;
		for (int i = 1 ; i < sommets.size() ; i++){
			float x_pos = sommets[i].position[0];
			float y_pos = sommets[i].position[1];
			float z_pos = sommets[i].position[2];
			if (x_pos < minPoint[0]){
				minPoint[0] = x_pos;
			}else if (x_pos > maxPoint[0]){
				maxPoint[0] = x_pos;
			}
			if (y_pos < minPoint[1]){
				minPoint[1] = y_pos;
			}else if (y_pos > maxPoint[1]){
				maxPoint[1] = y_pos;
			}
			if (z_pos < minPoint[2]){
				minPoint[2] = z_pos;
			}else if (z_pos > maxPoint[2]){
				maxPoint[2] = z_pos;
			}
		}
		for (int k = 0 ; k < 3 ; k++){ // On met un peu d'écart pour prendre une grille englobante un peu plus grande que le mesh
			minPoint[k] -= 0.01;
			maxPoint[k] += 0.01;
		}
		
		Vec3 vectors; // Contient les dimensions de la cellule en x, y et z (attention, les dimensions sont divisées par 2, ça nous arrange)
		Vec3 refPoint;// Point de référence de la racine. Les points de références sont les points aux centres de chaques cellules
		for (int k = 0 ; k < 3 ; k++){
			vectors[k] = (maxPoint[k] - minPoint[k]) / 2; // Ne surtout pas normaliser
			refPoint[k] = (minPoint[k] + maxPoint[k]) / 2; // Centre de la cellule sur les 3 dimensions
		}
		
		// sommets.size() est la taille du tableau lesListesTriangles
		KdNode *k = new KdNode(sommets,0,refPoint,vectors, lesListesTriangles); // Pour l'instant, la première division se fait sur l'axe x
		//k->showSommet(); // Avant tri (montre juste la coordonée x des sommets)

		k->divide(); // Construit le KD Tree en faisant la répartition des sommets

		//k->showSommet(); // Après tri (montre juste la coordonée x des sommets)
		
		this->kdRacine = k; // KdNode qui est la racine de l'arbre
    }

    void draw() {
        // iterer sur l'ensemble des objets, et faire leur rendu :
        for( unsigned int It = 0 ; It < meshes.size() ; ++It ) {
            Mesh const & mesh = meshes[It];
            mesh.draw();
        }
        for( unsigned int It = 0 ; It < spheres.size() ; ++It ) {
            Sphere const & sphere = spheres[It];
            sphere.draw();
        }
        for( unsigned int It = 0 ; It < squares.size() ; ++It ) {
            Square const & square = squares[It];
            square.draw();
        }
    }
    
    Material materialFromIntersection(RaySceneIntersection rsi){
    	Material m;
    	int n = rsi.typeOfIntersectedObject;
    	int ind = rsi.objectIndex;
    	
    	if (n == 0){
    		m  = spheres[ind].material;
    		
    		if (m.isTextured){ // Si la sphère doit avoir une texture, on va récupérer le texel correspondant aux bonnes coordonnées u,v de la sphère 
    			// m.texture contient le chemin vers la texture
    			Vec3 d = rsi.raySphereIntersection.intersection - spheres[ind].m_center;
    			
    			// Calcul des coordonnées u,v (voir sur diapo https://physique.cmaisonneuve.qc.ca/svezina/nyc/note_nyc/NYC_XXI_Chap%206.7.pdf)
    			// Dans ce cas (pour les objets de type Sphere, elles sont forcément entre 0 et 1)
    			float u = 0.5 + (atan2(d[2],d[0])/(2.*M_PI));
    			float v = 0.5 - (asin(d[1]/spheres[ind].m_radius)/M_PI);
    			//std::cout << "u = " << u << ", v = " << v << std::endl;
    			
    			// Utiliser u et v pour aller lire le bon texel dans la texture
    			
    			// ppmLoader::ImageRGB texture = spheres[ind].sphereTexture; Très bizarre --> Ca ralentit énormément le rendu de faire cette ligne
    			// Du coup je récupère les champs un à un
    			int width = spheres[ind].sphereTexture.w;
    			int height = spheres[ind].sphereTexture.h;
    			// std::vector<ppmLoader::RGB> data = spheres[ind].sphereTexture.data; Pareil, ça ralentit beaucoup de faire ça
    			
    			int nLigne = floor(v * (height - 1));
    			int nColonne = floor(u * (width - 1));
    			ppmLoader::RGB texel = spheres[ind].sphereTexture.data[nLigne * width + (width - 1 - nColonne)]; // On fait (width - 1 - nColonne) car sinon la texture se lisait à l'envers
    			m.diffuse_material = Vec3(texel.r/255.0,texel.g/255.0,texel.b/255.0); // Attention à bien diviser par 255 pour avoir des valeurs entre 0 et 1			
    			
    		}
    		
    	}
    	else if (n == 1){
    		m = squares[ind].material;
    		
    		if (m.isTextured){ // Si le square doit avoir une texture, on va récupérer le texel correspondant aux bonnes coordonnées u,v du square 
    			// m.texture contient le chemin vers la texture
    			Vec3 pLocal = rsi.raySquareIntersection.intersection - squares[ind].m_bottom_left;
    			
    			// Calcul des coordonnées u,v (formule trouvé sur internet)
    			float u = Vec3::dot(pLocal,squares[ind].m_right_vector) / squares[ind].m_right_vector.squareNorm();
    			float v = Vec3::dot(pLocal,squares[ind].m_up_vector) / squares[ind].m_up_vector.squareNorm();
    			//std::cout << "u = " << u << ", v = " << v << std::endl;
    			
    			// Utiliser u et v pour aller lire le bon texel dans la texture
    			
    			int width = squares[ind].squareTexture.w;
    			int height = squares[ind].squareTexture.h;
    			
    			int nLigne = floor(v * (height - 1));
    			int nColonne = floor(u * (width - 1));
    			
    			ppmLoader::RGB texel = squares[ind].squareTexture.data[(height - 1 - nLigne) * width + nColonne]; // On fait (width - 1 - nColonne) car sinon la texture se lisait à l'envers
    			m.diffuse_material = Vec3(texel.r/255.0,texel.g/255.0,texel.b/255.0); // Attention à bien diviser par 255 pour avoir des valeurs entre 0 et 1			
    			
    		}
    		
    	}
    	else if (n == 2){
    		m = meshes[ind].material; // A noter que tous les triangles d'un même mesh ont le même material
    		
    		if (m.isTextured){
    			// Je défini ici les coordonnées uv des MeshVertex --> (0,0) pour le sommet 1 ; (1,0) pour le sommet 2 ; (0,1) pour le sommet 3
				meshes[ind].vertices[rsi.mt.v[0]].u = 0.0f;
				meshes[ind].vertices[rsi.mt.v[0]].v = 0.0f;
				meshes[ind].vertices[rsi.mt.v[1]].u = 1.0f;
				meshes[ind].vertices[rsi.mt.v[1]].v = 0.0f;
				meshes[ind].vertices[rsi.mt.v[2]].u = 0.0f;
				meshes[ind].vertices[rsi.mt.v[2]].v = 1.0f;
		        
    			RayTriangleIntersection rtriangle = rsi.rayTriangleIntersection;
    			float w0 = rtriangle.w0;
        		float w1 = rtriangle.w1;
        		
        		// Même principe que l'interpolation des normales
        		float u = meshes[ind].vertices[rsi.mt.v[1]].u * (1 - w0 - w1) + meshes[ind].vertices[rsi.mt.v[2]].u * w0 + meshes[ind].vertices[rsi.mt.v[0]].u * w1 ;
        		float v = meshes[ind].vertices[rsi.mt.v[1]].v * (1 - w0 - w1) + meshes[ind].vertices[rsi.mt.v[2]].v * w0 + meshes[ind].vertices[rsi.mt.v[0]].v * w1 ;
        		
        		// Utiliser u et v pour aller lire le bon texel dans la texture
    			
    			int width = meshes[ind].meshTexture.w;
    			int height = meshes[ind].meshTexture.h;
    			
    			int nLigne = floor(v * (height - 1));
    			int nColonne = floor(u * (width - 1));
    			
    			ppmLoader::RGB texel = meshes[ind].meshTexture.data[nLigne * width + nColonne];
    			m.diffuse_material = Vec3(texel.r/255.0,texel.g/255.0,texel.b/255.0); // Attention à bien diviser par 255 pour avoir des valeurs entre 0 et 1	
				
    		}
    		
    	}
    	return m;
    }

    std::vector<Vec3> getInfoFromIntersection(RaySceneIntersection rsi){ // Selon le type de l'objet intersecté ayant crée l'intersection en paramètre, retourne un vector contenant dans l'ordre l'intersection, la normale, la couleur de l'intersection
    	std::vector<Vec3> res;
    	Vec3 intersection;
    	Vec3 normal;
    	Vec3 intersectionColor = materialFromIntersection(rsi).diffuse_material;
    	
    	 if (rsi.typeOfIntersectedObject == 0){ // Sphère
        	RaySphereIntersection rsphere = rsi.raySphereIntersection;
        	normal = rsphere.normal; // Déjà normalisé
        	intersection = rsphere.intersection;        	
        }else if (rsi.typeOfIntersectedObject == 1){ // Plan
        	RaySquareIntersection rsquare = rsi.raySquareIntersection;
        	normal = rsquare.normal; // Déjà normalisé
        	intersection = rsquare.intersection;
        }else if (rsi.typeOfIntersectedObject == 2){ // Triangle
        	RayTriangleIntersection rtriangle = rsi.rayTriangleIntersection;
        	intersection = rtriangle.intersection;
        	
        	// On a besoin de ça pour les interpolations (couleurs et normales)
        	MeshTriangle mt = rsi.mt;
        	int indMesh = rsi.objectIndex;
        	float w0 = rtriangle.w0;
        	float w1 = rtriangle.w1;
        	
        	// Sans normale interpolée (À COMMENTER/DÉCOMMENTER)
        	normal = rtriangle.normal; // Déjà normalisé	
        	
        	
        	/*
        	// Avec normales interpolées (À COMMENTER/DÉCOMMENTER)
        	Vec3 n0 = meshes[indMesh].vertices[mt.v[0]].normal;
        	Vec3 n1 = meshes[indMesh].vertices[mt.v[1]].normal;
        	Vec3 n2 = meshes[indMesh].vertices[mt.v[2]].normal;
        	normal = n1 * (1 - w0 - w1) + n2 * w0 + n0 * w1; // Attention à l'ordre des n0, n1 et n2
        	*/
        	
        	/*
        	// Avec couleurs interpolée (À COMMENTER/DÉCOMMENTER)
        	Vec3 c0 = meshes[indMesh].vertices[mt.v[0]].color;// Couleur du sommet 0
        	Vec3 c1 = meshes[indMesh].vertices[mt.v[1]].color;// Couleur du sommet 1
        	Vec3 c2 = meshes[indMesh].vertices[mt.v[2]].color;// Couleur du sommet 2
        	intersectionColor = c1 * (1 - w0 - w1) + c2 * w0 + c0 * w1;
        	*/
      		
      		  	
        } // Rentre forcément dans l'un de ces 3 cas sinon il y a un problème
    	
    	res.push_back(intersection);
    	res.push_back(normal);
    	res.push_back(intersectionColor);
    	return res;
    }


    RaySceneIntersection computeIntersection(Ray const & ray) {
        RaySceneIntersection result;
        
        for( unsigned int It = 0 ; It < spheres.size() ; ++It ) { // Calcule les intersections pour les sphères
        	Sphere const & sphere = spheres[It];
            RaySphereIntersection rsi = sphere.intersect(ray);
            if (rsi.intersectionExists && rsi.t >= 0 && rsi.t < result.t){ // On s'assure de prendre l'intersection la plus proche, et positive (ou nulle)
				result.intersectionExists = true;
				result.typeOfIntersectedObject = 0; // 0 pour les sphères, 1 pour les square (utilisé dans materialFromIntersection)
				result.objectIndex = It;
				result.t = rsi.t;
				result.raySphereIntersection = rsi;
            }
        }
        
        for( unsigned int It = 0 ; It < squares.size() ; ++It ) { // Calcule les intersections pour les squares
        	Square const & square = squares[It];
            RaySquareIntersection rsi = square.intersect(ray);
            /* On s'assure de prendre l'intersection la plus proche, et supérieur à une très petite valeur
               Auparavant, on pouvait se permettre de prendre l'intersection dès que t était positif ou nulle, mais depuis l'ajout des surfaces réfléchissantes, il a fallu
               restreindre à une valeur suffisament petite, ici 1e-4
            */
            if (rsi.intersectionExists && rsi.t > 1e-4 && rsi.t < result.t){ 
				result.intersectionExists = true;
				result.typeOfIntersectedObject = 1; // 0 pour les sphères, 1 pour les square (utilisé dans materialFromIntersection)
				result.objectIndex = It;
				result.t = rsi.t;
				result.raySquareIntersection = rsi;
            }
        }
        
        
        
        // À COMMENTER/DÉCOMMENTER
        // Sans KD Tree (dans ce cas on peut enlever si on veut l'appel à buildKDTree qui est dans main.cpp il ne sert à rien) ----------------------------------------------------------
        for( unsigned int It = 0 ; It < meshes.size() ; ++It ) { // Calcule les intersections pour les triangles (qui forment les meshes)
        	std::vector<MeshVertex> vertices = meshes[It].vertices; // Ensemble des MeshVertex du mesh
        	
        	for( unsigned int It_tr = 0 ; It_tr < meshes[It].triangles.size() ; ++It_tr ) {
        		
        		MeshTriangle const & mt = meshes[It].triangles[It_tr];
        		Triangle triangle = Triangle(vertices[mt.v[0]].position, vertices[mt.v[1]].position, vertices[mt.v[2]].position);
        		
		        RayTriangleIntersection rsi = triangle.intersect(ray);
		        if (rsi.intersectionExists && rsi.t > 1e-4 && rsi.t < result.t){ // On s'assure de prendre l'intersection la plus proche,
					result.intersectionExists = true;
					result.typeOfIntersectedObject = 2; // 0 pour les sphères, 1 pour les square (utilisé dans materialFromIntersection)
					result.objectIndex = It;
					result.t = rsi.t;
					result.rayTriangleIntersection = rsi;
					result.mt = mt; // On conserve le triangle dans lequel se trouve l'intersection, on en aura besoin pour l'interpolation des normales
		        }
		    }
        }
        // Fin Sans KD Tree ------------------------------------------------------------------------------------------------------------------------------------------------------------
        
        
        
        /*
        À COMMENTER/DÉCOMMENTER si vous voulez essayer la version kd tree
        // Avec KD Tree ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        
        
        if (meshes.size() == 0){ // On s'arrête là s'il n'y a pas de mesh
    		return result;
    	}
    	
        std::vector<MeshTriangle> *allTriangle = new std::vector<MeshTriangle>(); // Pour conserver les indices des triangles pour lesquels on testera l'intersection avec le rayon une fois l'exploration de l'arbre fini UNIQUEMENT
        
        exploreNodes(kdRacine, ray, allTriangle);
		
		
		      
        // Temporaire, pour l'instant avec un seul mesh
        std::vector<MeshVertex> vertices = meshes[0].vertices; // Ensemble des MeshVertex du mesh
        
        for (int i = 0 ; i < allTriangle->size() ; i++){
        	MeshTriangle const & mt = allTriangle->at(i);
        	Triangle triangle = Triangle(vertices[mt.v[0]].position, vertices[mt.v[1]].position, vertices[mt.v[2]].position);
        	
        	RayTriangleIntersection rsi = triangle.intersect(ray);
		    if (rsi.intersectionExists && rsi.t > 1e-4 && rsi.t < result.t){ // On s'assure de prendre l'intersection la plus proche,
				result.intersectionExists = true;
				result.typeOfIntersectedObject = 2; // 0 pour les sphères, 1 pour les square (utilisé dans materialFromIntersection)
				result.objectIndex = 0; // Ici 0 pour le numéro du mesh, comme je fais avec un seul mesh pour l'instant c'est forcément 0
				result.t = rsi.t;
				result.rayTriangleIntersection = rsi;
				result.mt = mt; // On conserve le triangle dans lequel se trouve l'intersection, on en aura besoin pour l'interpolation des normales
			}
        }
        
        
        */
        // Fin Avec KD Tree ----------------------------------------------------------------------------------------------------------------------------------------------------------------
        
        
        
        return result;
    }
    
    void exploreNodes(KdNode *cn, Ray const & ray, std::vector<MeshTriangle> *acc){
        if (cn->intersectNode(ray)){
        	//std::cout << "Intersecte\n";
        	if (cn->isLeaf()){
        		//std::cout << "C'est une feuille de l'arbre\n";
        		std::vector<MeshVertex> sommetsCellule = cn->getSommets();
        		std::vector<std::vector<MeshTriangle>> llt = cn->getLesListesTriangles();
        		for (int i = 0 ; i < sommetsCellule.size() ; i++){
        			//std::cout << "Taille = " << sommetsCellule[i].listeTriangles.size() << "\n";
        			for (int j =0 ; j < llt[i].size() ; j++){
        				
        				acc->push_back(llt[i][j]);
        			}
        		}
        	} else { // Sinon explorer fils gauche et droit (s'il existe)
        		if (cn->hasLeftChild()){ 
        			exploreNodes(cn->getLeftChild(), ray, acc);
        		}
        		if (cn->hasRightChild()){
        			exploreNodes(cn->getRightChild(), ray, acc);
        		}
        	}
        }
        
    }
    
    float getDist(Vec3 a, Vec3 b){
    	return sqrt(pow(b[0]-a[0],2) + pow(b[1]-a[1],2) + pow(b[2]-a[2],2)) ;
    }
    
    bool computeShadeIntersection(Ray const & i_to_light, float dist_light) { // Retourne true s'il y a une intersection créeant une ombre pour le point d'origine du rayon en paramètre
    	/*
    	Pour vérifier les intersections créeant des ombres, ce serait mieux de créer des nouvelles 
    	fonctions intersect pour Sphere, Square et Triangle, on a pas besoin de toutes les infos que 
    	retournent les fonctions actuellement utilisées 
    	*/
    	
    	Vec3 o = i_to_light.origin(); // Origine du rayon
    	Vec3 d = i_to_light.direction(); // Vecteur directeur du rayon
    	for( unsigned int It = 0 ; It < spheres.size() ; ++It ) { // Calcule les intersections pour les sphères
        	Sphere const & sphere = spheres[It];
            RaySphereIntersection rsi = sphere.intersect(i_to_light);
            if (rsi.intersectionExists && rsi.t > 1e-4 && getDist(o, o + rsi.t * d) < dist_light){ // 10^-4 pour qu'il n'y ait pas une ombre créé par une intersection avec lui-même
				return true;	
            }
        }
              
        for( unsigned int It = 0 ; It < squares.size() ; ++It ) { // Calcule les intersections pour les squares
        	Square const & square = squares[It];
            RaySquareIntersection rsi = square.intersect(i_to_light);
            if (rsi.intersectionExists && rsi.t > 1e-4 && getDist(o, o + rsi.t * d) < dist_light){ // 10^-4 pour qu'il n'y ait pas une ombre créé par une intersection avec lui-même
				return true;
            }
        }
        
        for( unsigned int It = 0 ; It < meshes.size() ; ++It ) { // Calcule les intersections pour les triangles (qui forment les meshes)
        	std::vector<MeshVertex> vertices = meshes[It].vertices; // Ensemble des MeshVertex du mesh
        	for( unsigned int It_tr = 0 ; It_tr < meshes[It].triangles.size() ; ++It_tr ) {
        		MeshTriangle const & mt = meshes[It].triangles[It_tr];
        		Triangle triangle = Triangle(vertices[mt.v[0]].position, vertices[mt.v[1]].position, vertices[mt.v[2]].position); 
        		RayTriangleIntersection rsi = triangle.intersect(i_to_light);
        		if (rsi.intersectionExists && rsi.t > 1e-4 && getDist(o, o + rsi.t * d) < dist_light){ // 10^-4 pour qu'il n'y ait pas une ombre créé par une intersection avec lui-même
					return true;	
           		}
		}
        }
        
        return false;
    }


	Vec3 computeLightPhong(RaySceneIntersection rsi, Ray ray, int NRemainingBounces){
		Material m = materialFromIntersection(rsi);
		
		std::vector<Vec3> info_intersection = getInfoFromIntersection(rsi); // Contient l'intersection, la normale, et la couleur de l'intersection
        Vec3 normal_intersection = info_intersection[1];
        Vec3 intersection = info_intersection[0];
        m.diffuse_material = info_intersection[2];
		
		if (m.type == Material_Mirror && NRemainingBounces > 0){ // Pour les matériaux réfléchissant
	    		
	    	// Calcul du rayon réfléchi à partir du rayon initial
			Vec3 new_dir=ray.direction() - 2 * (Vec3::dot(ray.direction(),normal_intersection)) * normal_intersection;
            new_dir.normalize();
			
			Ray new_ray = Ray(intersection, new_dir);
			m.diffuse_material = rayTraceRecursive(new_ray, NRemainingBounces - 1); // Appel récursif puisque computeLightPhong est appelé dans rayTraceRecursive
			
			return m.diffuse_material; // ATTENTION : Dans le cas d'un matériau réfléchissant, je ne calcule pas les ombres. A voir s'il ne faut quand même pas prendre en compte certain paramètre, notamment la couleur de la lumière
    	}
    	
    	
    	// Pour les sphères transparentes, avec rayons réfractés
        if (m.transparency != 0.0){ // Je teste pas si l'objet intersecté est une sphère, je considère que c'est toujours le cas si m.transparency est != de 0
	    	
        	Vec3 vecteur_observ = ray.origin() - intersection; // Vecteur allant de l'intersection à l'origine du rayon
        	Vec3 normal_utilise; // La normale dépend si on rentre ou sort de la sphère
	    	vecteur_observ.normalize();
	    	
	    	// n1 * sin(i1) = n2 * sin(i2)
	    	double i1; // Angle (en radians) entre le vecteur de l'intersection à la caméra et le vecteur normal de l'intersection. ATTENTION si on rentre/sort de la sphère ça change le calcul
	    	double i2; // Angle (en radians) du rayon réfracté par rapport à la normale de l'intersection
	    	float n1; // Indice de réfraction du milieu du rayon incident
	    	float n2;  // Indice de réfraction du milieu du rayon réfracté
	    	
	    	if (Vec3::dot(vecteur_observ,normal_intersection) >= 0){ // On rentre dans la sphère, donc l'indice de réfraction n2 dans la formule est celui du matériau de la sphère
	    		n1 = INDICE_REFRACTION_VIDE; // = 1.0
	    		n2 = m.index_medium;
	    		normal_utilise = normal_intersection;
	    	} else { // On sort de la sphère donc l'indice de réfraction n2 dans la formule est celui de l'extérieur (on prendra l'indice de réfraction du vide donc 1.0)
	    		n1 = m.index_medium;
	    		n2 = INDICE_REFRACTION_VIDE; // = 1.0
	    		normal_utilise = normal_intersection * (-1);
	    	}
	    	
	    	i1 = acos(Vec3::dot(vecteur_observ,normal_utilise) / ( vecteur_observ.length() * normal_utilise.length()) );
	    	i2 = asin((n1 * sin(i1)) / n2);
	    	
	    	Vec3 M = (normal_utilise * cos(i1) - vecteur_observ ) / sin(i1);// Voir diapo récupéré sur internet
	    	
	    	M.normalize(); // A priori c'est déjà normalisé mais dans le doute je le fais quand même
	    	Vec3 T = sin(i2) * M - cos(i2) * normal_utilise; // Vecteur directeur du rayon réfracté
	    	T.normalize();
	    	
	    	Ray new_ray = Ray(intersection + T * 1e-3, T); // Le rayon réfracté ne doit surtout pas démarrer là où a été trouvé l'intersection précédente, sinon il y aura un problème d'auto-intersection et donc il y aura des valeurs nan. Pour empecher ça, j'ajoute T * 1e-3 à l'origine du rayon
	    	
	    	return rayTraceRecursive(new_ray, NRemainingBounces); 
            // Pour l'instant, uniquement des sphères parfaitement transparentes, donc on ne prends pas en compte la couleur
            // de la sphère. Si on voulait que ce soit le cas, il faudrait modifier m.diffuse_material au lieu de
            // directement retourner le résultat (comme ci-dessus) et continuer la fonction computeLightPhong avec cette valeur.
            // ATTENTION --> Je pense qu'il faudra modifier m.diffuse_material uniquement au moment de la première réfraction 
            // du rayon, pour les autres réfractions (s'il y en a) il faudrait là par contre directement retourner le résultat 
            // (comme ci-dessus)
	    	
	    	/*
	    	Normalement pour la réfraction c'est bon, voir juste s'il ne faut pas s'occuper
	    	du problème de réfraction interne
	    	*/
	   		
        }
        
        
        
        
        
		Vec3 vecteur_observ = ray.origin() - intersection; // V dans diapo (vecteur allant de l'intersection à la caméra càd l'origine du rayon)		
		vecteur_observ.normalize(); // C'est le même pour toutes les lumières donc on le calcule qu'une seule fois
		
		Vec3 ambientLight = AMBIANT_INTENSITY * m.ambient_material; // Attention : Par défaut, ambient_material est à Vec3(0.,0.,0.) dans le code
		Vec3 diffuse_and_specular_light = Vec3(0.,0.,0.); // Combiné à cause du calcul des ombres diffuses
		
		
        for (unsigned int i = 0 ; i < lights.size() ; i++){
        	
        	if (lights[i].type == LightType_Spherical){ // Lumière en un point, ne crée que des ombres simples
        		Vec3 vecteur_lumiere = lights[i].pos - intersection; // Ll dans le diapo
		    	vecteur_lumiere.normalize();
		    	
		    	Ray i_to_light = Ray(intersection, vecteur_lumiere); // Rayon du point vers la lumière
		    	float dist_light = getDist(intersection, lights[i].pos); // Distance entre le point et la lumière 
		    	
		    	if (computeShadeIntersection(i_to_light, dist_light)){ // Si on trouve une intersection entre le point et la lumière lights[i]
					continue; // Passe à la lumière suivante
		    	}
		    	
				Vec3 vecteur_reflexion = 2 * Vec3::dot(normal_intersection,vecteur_lumiere) * normal_intersection - vecteur_lumiere; // R dans diapo --> R = 2 * (N . L) * N - L
				vecteur_reflexion.normalize();
				
				diffuse_and_specular_light += DIFFUSE_INTENSITY * m.diffuse_material * std::max(Vec3::dot(vecteur_lumiere,normal_intersection),0.f) 
											+ SPECULAR_INTENSITY * m.specular_material * pow(std::max(Vec3::dot(vecteur_reflexion, vecteur_observ), 0.f),m.shininess); 
											// On utilise max (avec 0) pour avoir une valeur positive
											
        	}else if (lights[i].type == LightType_Quad){ // Lumière créeant des ombres douces
        		// Echantilloner le carré de la lumière
        		// Pour chaque échantillon, lancer un rayon depuis l'intersection (dont on calcule la lumière) vers le point de l'échantillon
        		// Si aucune intersection pour ce rayon --> +1 à nb_light
        		
        		float nb_light = 0.; // Nombre de point sur le carré de lumière qui éclaire (sans intersection) le point courant
        		float nb_point_light = lights[i].nb_points;; // Nombre de point total sur le carré de lumière
        		
        		Square carre = lights[i].quad;
        		
        		for(unsigned int p = 0 ; p < nb_point_light ; p++){
        			Vec3 courant = carre.m_bottom_left + ((rand()%100) / 99.) *carre.m_right_vector + ((rand()%100) / 99.) * carre.m_up_vector; // Donne un nouveau point aléatoirement sur le carré à chaque itération
        		 		
        		 	Vec3 vecteur_lumiere = courant - intersection; // Ll dans le diapo
					vecteur_lumiere.normalize();
						
					Ray i_to_light = Ray(intersection, vecteur_lumiere);
					float dist_light = getDist(intersection, courant); 
        		 		
        		 	if (!computeShadeIntersection(i_to_light, dist_light)){ // S'il n'y a pas d'intersection entre le point courant sur le carré et l'intersection dont on calcule l'ombre
						nb_light += 1.; // Ajoute 1 au nombre de lumière éclairant l'intersection
					}
        		}
        		
        		
        		if (nb_light != 0.){
        			Vec3 vecteur_lumiere = lights[i].pos - intersection; // Ll dans le diapo
		    		vecteur_lumiere.normalize();
		    		Vec3 vecteur_reflexion = 2 * Vec3::dot(normal_intersection,vecteur_lumiere) * normal_intersection - vecteur_lumiere; // R dans diapo --> R = 2 * (N . L) * N - L
					vecteur_reflexion.normalize();
				
					diffuse_and_specular_light += DIFFUSE_INTENSITY * m.diffuse_material * std::max(Vec3::dot(vecteur_lumiere,normal_intersection),0.f) 
											+ SPECULAR_INTENSITY * m.specular_material * pow(std::max(Vec3::dot(vecteur_reflexion, vecteur_observ), 0.f),m.shininess); 
											// On utilise max (avec 0) pour avoir une valeur positive
											
					diffuse_and_specular_light *= (nb_light/nb_point_light);
        		} 
        	}
        	
        }
        
        return Vec3::compProduct(ambientLight + diffuse_and_specular_light, lights[0].material); // Pour prendre en compte la couleur de la lumière (dans lights[i].material) 
	}


    Vec3 rayTraceRecursive( Ray ray , int NRemainingBounces ) {

        //TODO RaySceneIntersection raySceneIntersection = computeIntersection(ray);
        Vec3 color = Vec3(0.,0.,0.); // Retournera du noir si aucune intersection n'est trouvé (ATTENTION : valeur entre 0 et 1, et non 0 et 255)
        RaySceneIntersection rsi = computeIntersection(ray);
        if (rsi.intersectionExists){
        	
        	
        	/* Pour tester la phase 1
        	Material m = materialFromIntersection(rsi); // La couleur retourné sera celle de la première intersection rencontrée par le rayon 
        	color = m.diffuse_material;
        	*/
        	
        	// Pour tester la phase 2
        	// Pour vérifier de quel côté on regarde les plans (1) et les triangles (2) (affichés en noir si on les regarde dans le mauvais sens)
        	std::vector<Vec3> info_intersection = getInfoFromIntersection(rsi); // Normale de l'objet intersecté
        	if (rsi.typeOfIntersectedObject == 1 || rsi.typeOfIntersectedObject == 2){
        		Vec3 versOrigine = ray.origin() - info_intersection[0];
		        versOrigine.normalize();
		        if (Vec3::dot(info_intersection[1], versOrigine) < 0){
		        	return color; // On retourne du noir quand on regarde pas du bon côté du square/triangle
		        }
        	}
        	
        	color = computeLightPhong(rsi, ray, NRemainingBounces);
        	
        }      
        return color;
    }


    Vec3 rayTrace( Ray const & rayStart ) {
        Vec3 color = rayTraceRecursive(rayStart,2);
        return color; 
    }

    void setup_single_sphere() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(5,5,5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }
        {
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(0. , 0. , 0.);
            s.m_radius = 1.f;
            
            s.material.isTextured = true;
            // load_ppm dans imageLoader.h
            load_ppm(s.sphereTexture,"img/sphereTextures/s7.ppm"); 
            //std::cout << "w = " << s.sphereTexture.w << ", h = " << s.sphereTexture.h << std::endl;
            
            s.build_arrays();
            //s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,1.,1 );
            s.material.specular_material = Vec3( 0.2,0.2,0.2 );
            s.material.shininess = 20;
        }
    }

    void setup_single_square() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(-5,5,5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }

        {
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            
            s.material.isTextured = true;
            // load_ppm dans imageLoader.h
            load_ppm(s.squareTexture,"img/squareTextures/paysage.ppm");
            //std::cout << "w = " << s.squareTexture.w << ", h = " << s.squareTexture.h << std::endl;
            
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 0.8,0.8,0.8 );
            s.material.shininess = 20;
        }
    }
    
    void setup_my_scene(){ 
    	meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(-5,5,5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }
        
        { // Tétraèdre mesh
            meshes.resize( meshes.size() + 1 );
            Mesh & m = meshes[meshes.size() - 1];
            m.loadOFF("src/data/model/my_mesh.off");
            m.build_arrays();
            
	    	m.material.diffuse_material = Vec3( 0.5,0.8,0.1 );
	    	m.material.specular_material = Vec3( 0.8,0.8,0.8 );
	    	m.material.shininess = 20;
	    	
	    	m.material.isTextured = true;
            // load_ppm dans imageLoader.h
            load_ppm(m.meshTexture,"img/squareTextures/paysage.ppm");
        } 
    }
    
    void setup_cornell_box(){
    	meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();
		
	// Lumière utilisé pour les ombres simples
        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3( 0.0, 1.5, 0.0 );
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }
        
        
        
        /*
        // Une deuxième lumière (ombre simple) si besoin
        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3( 0.0, -1., 2.0 );
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }
        */
        
        
        /*
        // Source lumineuse = un carré pour faire des ombres douces
        
        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3( 0.0, 1.5, 0.0 );
            light.powerCorrection = 2.f;
            
            light.type = LightType_Quad;
            Square s;
            s.setQuad(Vec3(-1., 1.5, 1.), Vec3(1., 0, 0.), Vec3(0., 0., -1.), 2., 2.); // Carré autour du point de la lumière
            light.quad = s;
            light.nb_points = 40; // Attention : Plus le nombre de point sur le carré est élevé, plus le rendu prendra du temps
            
            light.material = Vec3(1.,1.,1.); // Lumière blanche
            //light.material = Vec3(0.6,0.7,1.); // Lumière légèrement bleu
            light.isInCamSpace = false;
        }
        */
        
       	
        { // Mur de fond
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-2., -2., -2.), Vec3(1., 0., 0.), Vec3(0., 1., 0.), 4., 4.);
            s.build_arrays();
            //s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 0.,0.,1. ); // Plan bleu
            s.material.specular_material = Vec3( 0.8,0.8,0.8 );
            s.material.shininess = 20;
        }
        { // Plafond
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-2., 2., -2.), Vec3(1., 0., 0.), Vec3(0., 0, 1.), 4., 4.);
            s.build_arrays();
            //s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 0.8,0.2,0.8 ); // Plan rose
            s.material.specular_material = Vec3( 0.8,0.8,0.8 );
            s.material.shininess = 20;
        }
        { // Sol 
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-2., -2., 2.), Vec3(1., 0, 0.), Vec3(0., 0., -1.), 4., 4.);
            s.build_arrays();
            //s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 0.,1.,0. ); // Plan vert
            s.material.specular_material = Vec3( 0.8,0.8,0.8 );
            s.material.shininess = 20;
        }
        { // Mur de gauche 
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-2., -2., 2.), Vec3(0. , 0, -1.), Vec3(0., 1., 0.), 4., 4.);
            s.build_arrays();
            //s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,1.,0. ); // Plan blanc
            s.material.specular_material = Vec3( 0.8,0.8,0.8 );
            s.material.shininess = 20;
        }
        { // Mur de droite
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(2., -2., -2.), Vec3(0., 0, 1.), Vec3(0., 1, 0.), 4., 4.);
            s.build_arrays();
            //s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,0.,0. ); // Plan rouge
            s.material.specular_material = Vec3( 0.8,0.8,0.8 );
            s.material.shininess = 20;
        }
        
        /*
        { // Tétraèdre mesh
            meshes.resize( meshes.size() + 1 );
            Mesh & m = meshes[meshes.size() - 1];
            std::cout << "Chargement du mesh" << std::endl;
            m.loadOFF("src/data/model/my_mesh.off");
            std::cout << "Fin du chargement du mesh" << std::endl;
            m.build_arrays();   
			//m.material.type = Material_Mirror;
			m.material.diffuse_material = Vec3( 0.5,0.8,0.1 );
			m.material.specular_material = Vec3( 0.8,0.8,0.8 );
			m.material.shininess = 20;
        } 
        */
        
        /*
        { // Mur de devant
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-2., -2., 2.), Vec3(0., 1., 0.), Vec3(1., 0., 0.), 4., 4.);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.1,0.3, 0.6 ); // Plan bleu
            s.material.specular_material = Vec3( 0.8,0.8,0.8 );
            s.material.shininess = 20;
        }
        */
	
	
        { //GLASS Sphere (J'ai enlevé le Mirroir sur cette sphère)
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1.0, -1.25, 0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.ambient_material = s.material.diffuse_material / 6.f;
            s.material.shininess = 16;
            s.material.transparency = 0.0;
            s.material.index_medium = 1; // Si on mettait 1, les rayons réfractés seraient les mêmes que les rayons incidents, donc on verrait parfaitement à travers de la sphère
        } 
        { //MIRRORED Sphere
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-1.0, -1.25, -0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            //s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,1.,1. ); // Sphère blanche
            s.material.specular_material = Vec3(  1.,1.,1. );
            s.material.ambient_material = s.material.diffuse_material / 6.f;
            s.material.shininess = 16;
            s.material.transparency = 0.;
            s.material.index_medium = 1;
        }
                
        /*
        { //Une sphère derrière la caméra pour voir si on la voit dans le reflet des surfaces réfléchissantes
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1.0, -1.25, 5.);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 0.,1.,0. ); // Sphère verte
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.ambient_material = s.material.diffuse_material / 6.f;
            s.material.shininess = 16;
            s.material.transparency = 0.0;
            s.material.index_medium = 0.0;
        }
        */ 
    }
    
    void setup_my_mesh(){
    	meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();
        
    	 {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3( -5.0, 5.0, 5.0 );
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }
        
        { //Sphere mesh
            meshes.resize( meshes.size() + 1 );
            Mesh & m = meshes[meshes.size() - 1];
            m.loadOFF("src/data/model/sphere.off");
            m.build_arrays();
			m.material.diffuse_material = Vec3( 0.5,0.8,0.1 ); 
			m.material.specular_material = Vec3( 0.8,0.8,0.8 );
			m.material.shininess = 20;
			
			m.material.isTextured = true;
		    // load_ppm dans imageLoader.h
		    load_ppm(m.meshTexture,"img/squareTextures/nourriture.ppm");
        } 
            
    }
    
    void setup_complete_scene(){
    	meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        { // Ombres simple
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3( 0.0, 1.5, 0.0 );
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }
        
       	
        { // Mur de fond
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-2., -2., -2.), Vec3(1., 0., 0.), Vec3(0., 1., 0.), 4., 4.);
            s.material.isTextured = true;
            load_ppm(s.squareTexture,"img/squareTextures/paysage.ppm");
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.,0.,1. ); // Plan bleu
            s.material.specular_material = Vec3( 0.8,0.8,0.8 );
            s.material.shininess = 20;
        }
        { // Plafond
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-2., 2., -2.), Vec3(1., 0., 0.), Vec3(0., 0, 1.), 4., 4.);
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 0.8,0.2,0.8 ); // Plan rose
            s.material.specular_material = Vec3( 0.8,0.8,0.8 );
            s.material.shininess = 20;
        }
        { // Sol 
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-2., -2., 2.), Vec3(1., 0, 0.), Vec3(0., 0., -1.), 4., 4.);
            s.build_arrays();
            //s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 0.,1.,0. ); // Plan vert
            s.material.specular_material = Vec3( 0.8,0.8,0.8 );
            s.material.shininess = 20;
        }
        { // Mur de gauche 
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-2., -2., 2.), Vec3(0. , 0, -1.), Vec3(0., 1., 0.), 4., 4.);
            s.material.isTextured = true;
            load_ppm(s.squareTexture,"img/squareTextures/nourriture.ppm");
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,1.,0. ); // Plan blanc
            s.material.specular_material = Vec3( 0.8,0.8,0.8 );
            s.material.shininess = 20;
        }
        { // Mur de droite
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(2., -2., -2.), Vec3(0., 0, 1.), Vec3(0., 1, 0.), 4., 4.);
            s.material.isTextured = true;
            load_ppm(s.squareTexture,"img/squareTextures/nourriture.ppm");
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,0.,0. ); // Plan rouge
            s.material.specular_material = Vec3( 0.8,0.8,0.8 );
            s.material.shininess = 20;
        }
        
        /*
        { // Tétraèdre mesh
            meshes.resize( meshes.size() + 1 );
            Mesh & m = meshes[meshes.size() - 1];
            std::cout << "Chargement du mesh" << std::endl;
            m.loadOFF("src/data/model/my_mesh.off");
            std::cout << "Fin du chargement du mesh" << std::endl;
            m.build_arrays();   
			//m.material.type = Material_Mirror;
			m.material.diffuse_material = Vec3( 0.5,0.8,0.1 );
			m.material.specular_material = Vec3( 0.8,0.8,0.8 );
			m.material.shininess = 20;
        } 
        */
	
        { //GLASS Sphere (J'ai enlevé le Mirroir sur cette sphère)
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1.0, -1.25, 0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            //s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.ambient_material = s.material.diffuse_material / 6.f;
            s.material.shininess = 16;
            s.material.transparency = 1.0;
            s.material.index_medium = 1.08; // Si on mettait 1, les rayons réfractés seraient les mêmes que les rayons incidents, donc on verrait parfaitement à travers de la sphère
        } 
        { //MIRRORED Sphere
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-1.0, -1.25, -0.5);
            s.m_radius = 0.75f;
            s.material.isTextured = true;
            load_ppm(s.sphereTexture,"img/sphereTextures/s7.ppm");
            s.build_arrays();
            //s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,1.,1. ); // Sphère blanche
            s.material.specular_material = Vec3(  1.,1.,1. );
            s.material.ambient_material = s.material.diffuse_material / 6.f;
            s.material.shininess = 16;
            s.material.transparency = 0.;
            s.material.index_medium = 0.;
        }
    }

};

#endif
