#ifndef KdNode_H
#define KdNode_H

#include "Mesh.h"

#define SEUIL 22 // La division d'une cellule continue tant que ce seuil n'est pas atteint

class KdNode { // Objet pour une cellule du Kd tree / Noeud d'un arbre ayant 0 ou 2 fils 
	private:
		std::vector<MeshVertex> sommets;
		KdNode *left;
		KdNode *right;
		int split_dim; // 0 pour x, 1 pour y, 2 pour z
		std::vector<std::vector<MeshTriangle>> lesListesTriangles;
		
		Vec3 refPoint; // Point au centre de la cellule
		
		// J'ai mit directement les 3 dimensions dans un seul paramètre de type de Vec3 (vectors) pour avoir un seul paramètre, mais sinon il faudrait 3 floats
		// Ca permit aussi de ne pas avoir à testé la valeur de split_dim dès qu'on veut faire quelque chose
		Vec3 vectors; // Permet de connaitre les dimensions de la cellule, vectors[0] pour dimension en x, vectors[1] pour dimension en y, vectors[2] pour dimension en 2,  
		
	public:
		KdNode(std::vector<MeshVertex> sommets, int split_dim, Vec3 refPoint, Vec3 vectors, std::vector<std::vector<MeshTriangle>> lesListesTriangles){
			// sommets.size() est la taille du tableau lesListesTriangles
			
			left = nullptr; // Pas de fils gauche par défaut
			right = nullptr; // Pas de fils droit par défaut
			this->sommets = sommets;
			this->split_dim = split_dim; // La racine commence avec 0 puis on alterne 0,1,2
			
			this->refPoint = refPoint;
			this->vectors = vectors;
			this->lesListesTriangles = lesListesTriangles;
		}
		
		std::vector<MeshVertex> getSommets(){
			//std::cout << "Nb sommet dans getSommets = " << (this->sommets).size() << "\n";
			return this->sommets;
		}
		
		std::vector<std::vector<MeshTriangle>> getLesListesTriangles(){
			return this->lesListesTriangles;
		}
		
		bool hasLeftChild(){ // true si le noeud a un fils gauche
			return (left != nullptr);
		}
		
		bool hasRightChild(){ // true si le noeud a un fils droit
			return (right != nullptr);
		}
		
		bool isLeaf(){ // true si le noeud est une feuille
			return !(this->hasLeftChild()) && !(this->hasRightChild());
		}
		
		KdNode* getLeftChild(){
			return (*this).left;
		}
		
		KdNode* getRightChild(){
			return (*this).right;
		}
		
		void divide(){ // Répartition des sommets dans le Kd-tree
			if (this->sommets.size() < SEUIL){
				//std::cout << "Nb Limite\n";
				return; // Fin de la division de cette cellule
			}
			
			this->trier();

			// On découpe à la médiane
			int half_size = this->sommets.size() / 2;
			std::vector<MeshVertex> premiereMoitie(this->sommets.begin(), this->sommets.begin() + half_size);
			std::vector<MeshVertex> deuxiemeMoitie(this->sommets.begin() + half_size, this->sommets.end());
			
			std::vector<std::vector<MeshTriangle>> leftListesTriangles;
			std::vector<std::vector<MeshTriangle>> rightListesTriangles;

			for (int i = 0 ; i < half_size ; i++){
				leftListesTriangles.push_back(this->lesListesTriangles[i]);
			}

			for (int i = half_size ; i < this->sommets.size() ; i++){
				rightListesTriangles.push_back(this->lesListesTriangles[i]);
			}

			Vec3 leftRefPoint = refPoint; // Point de référence pour le fils gauche (ne restera pas identique à celui du noeud parent)
			Vec3 rightRefPoint = refPoint; // Point de référence pour le fils droit (ne restera pas identique à celui du noeud parent)
			
			// On fait la moyenne entre le dernier de la première moitié et le premier de la deuxième moitié
			// Et on utilise split_dim pour savoir dans quel dimension on fait la division. On peut se permettre ça uniquement parceque les cellules sont alignées sur les axes.
			float separation = (this->sommets.at(half_size - 1).position[this->split_dim] + this->sommets.at(half_size).position[this->split_dim]) / 2;
			
			float minp = refPoint[this->split_dim] - this->vectors[this->split_dim]; // Point au minimum de la cellule sur la dimension split_dim 
			float maxp = refPoint[this->split_dim] + this->vectors[this->split_dim]; // Point au maximum de la cellule sur la dimension split_dim 
			
			// Par rapport au point de référence de la cellule parente, il n'y a que sur la dimension split_dim que ça change
			leftRefPoint[this->split_dim] = (minp + separation) / 2;
			rightRefPoint[this->split_dim] = (separation + maxp) / 2;
			
			Vec3 left_vectors = this->vectors; // Paramètre vectors pour le fils du noeud left (intialement le même que le noeud parent)
			Vec3 right_vectors = this->vectors; // Paramètre vectors pour le fils du noeud right (intialement le même que le noeud parent)
			
			// Par rapport au vectors de la cellule parente, il n'y a que sur la dimension split_dim que ça change
			left_vectors[this->split_dim] = leftRefPoint[this->split_dim] - minp; // On pourrait aussi faire : separation - leftRefPoint[this->split_dim]
			right_vectors[this->split_dim] = maxp - rightRefPoint[this->split_dim]; // On pourrait aussi faire : rightRefPoint[this->split_dim] - separation
			
			
			// On fait alterner split_dim avec 0 --> 1 --> 2 --> 0 --> ...
			left = new KdNode(premiereMoitie, this->split_dim == 2 ? 0 : this->split_dim + 1, leftRefPoint, left_vectors, leftListesTriangles);
			right = new KdNode(deuxiemeMoitie, this->split_dim == 2 ? 0 : this->split_dim + 1, rightRefPoint, right_vectors, rightListesTriangles);
			//std::cout << "On divise\n";
			left->divide();
			right->divide();
			
		}
		
		// Comparaison de 2 points, selon leur dimension x,y ou z (j'ai pas trouvé comment faire automatiquement 
		// une division selon split_dim, c'est plus dur que ça en a l'air, le problème c'est que ces fonctions doivent obligatoirement être static donc on peut pas utiliser this->split_dim)
		
		static bool comparaisonSommetX(MeshVertex a, MeshVertex b){
			return a.position[0] < b.position[0];
		}
		static bool comparaisonSommetY(MeshVertex a, MeshVertex b){
			return a.position[1] < b.position[1];
		}
		static bool comparaisonSommetZ(MeshVertex a, MeshVertex b){
			return a.position[2] < b.position[2];
		}
		
		
		void trier(){ // Trie les sommets qui sont dans la cellule, selon leur dimension split_dim (j'ai pas trouvé autrement que de tester split_dim)
			if (this->split_dim == 0) {
				std::sort(this->sommets.begin(), this->sommets.end(), comparaisonSommetX); 
			} else if (this->split_dim == 1) {
				std::sort(this->sommets.begin(), this->sommets.end(), comparaisonSommetY); 
			} else { // split_dim = 2
				std::sort(this->sommets.begin(), this->sommets.end(), comparaisonSommetZ); 
			}
		}
		
		std::vector<Square> get6Square(){ // Retourne le vector contenant les 6 Squares formant la cellule, grâce au point de référence (càd au milieu de la cellule) et ses dimensions (dans vectors)
			std::vector<Square> res;
			
			Square s1;
			
			// 1 (Face gauche)
            s1.setQuad(Vec3(this->refPoint[0] - this->vectors[0], this->refPoint[1] - this->vectors[1], this->refPoint[2] - this->vectors[2]),
            				 Vec3(0., 0., 1.), Vec3(0., 1., 0.), 2*this->vectors[2], 2*this->vectors[1]);
            res.push_back(s1);
            
            // 2 (Face droite)
            s1.setQuad(Vec3(this->refPoint[0] + this->vectors[0], this->refPoint[1] - this->vectors[1], this->refPoint[2] + this->vectors[2]),
            				 Vec3(0., 0., -1.), Vec3(0., 1., 0.), 2*this->vectors[2], 2*this->vectors[1]);
            res.push_back(s1);
            
            // 3 (Face haut)
            s1.setQuad(Vec3(this->refPoint[0] - this->vectors[0], this->refPoint[1] + this->vectors[1], this->refPoint[2] + this->vectors[2]),
            				 Vec3(1., 0., 0.), Vec3(0., 0., -1.), 2*this->vectors[0], 2*this->vectors[2]);
            res.push_back(s1);
            
            // 4 (Face bas)
            s1.setQuad(Vec3(this->refPoint[0] - this->vectors[0], this->refPoint[1] - this->vectors[1], this->refPoint[2] - this->vectors[2]),
            				 Vec3(1., 0., 0.), Vec3(0., 0., 1.), 2*this->vectors[0], 2*this->vectors[2]);
            res.push_back(s1);
            
            // 5 (Face devant)
            s1.setQuad(Vec3(this->refPoint[0] - this->vectors[0], this->refPoint[1] - this->vectors[1], this->refPoint[2] + this->vectors[2]),
            				 Vec3(1., 0., 0.), Vec3(0., 1., 0.), 2*this->vectors[0], 2*this->vectors[1]);
            res.push_back(s1);
            
            // 6 (Face derrière)
            s1.setQuad(Vec3(this->refPoint[0] + this->vectors[0], this->refPoint[1] - this->vectors[1], this->refPoint[2] - this->vectors[2]),
            				 Vec3(-1., 0., 0.), Vec3(0., 1., 0.), 2*this->vectors[0], 2*this->vectors[1]);
            res.push_back(s1);
           
            
			return res;
		}
		
		bool intersectNode(Ray const & ray){// true si le rayon en paramètre passe par la cellule
			 std::vector<Square> listeFaces = this->get6Square(); // Formation des 6 faces formant la cellule courante
			 
			 for (int i = 0 ; i < 6 ; i++){
			 	if (Vec3::dot(listeFaces[i].m_normal, ray.direction()) < 0){ // On regarde uniquement si le Square fait face au rayon (ça parait contre-intuitif mais c'est bien dot < 0)
			 		RaySquareIntersection rsi = listeFaces[i].intersect(ray);
				    if (rsi.intersectionExists && rsi.t > 1e-4){
						return true;
				    }
			 	}
			 }
			 return false; // Aucune intersection trouvé avec l'une des faces de la cellule
		}
		
		// Pour les tests
		void showSommet(){// Affiche la coordonnée x de tous les sommets pour vérifier s'ils sont triés
			for (int i = 0; i < this->sommets.size() ; i++){
				std::cout << (this->sommets)[i].position[0] << "\n";
			}
		}
};

#endif
