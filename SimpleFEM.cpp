#include <iostream>
#include <cstdlib>

#ifndef ELEMENT_COUNT
#define ELEMENT_COUNT 50
#endif

#ifndef ENVIRONMENT_TEMPERATURE
#define ENVIRONMENT_TEMPERATURE 40.0f
#endif

#ifndef HEAT_SOURCE_DENSITY
#define HEAT_SOURCE_DENSITY -150
#endif

#ifndef CONVECTION_COEFFICIENT
#define CONVECTION_COEFFICIENT 10.0f
#endif

#ifndef HEAT_TRANSFER_COEFFICIENT
#define HEAT_TRANSFER_COEFFICIENT 75.0f
#endif

#ifndef CROSSECTION_AREA
#define CROSSECTION_AREA 1.0f
#endif

#ifndef ROD_LENGTH
#define ROD_LENGTH 5.0f
#endif

#if HEAT_SOURCE_DENSITY > 0
#error HEAT_SOURCE_DENSITY must be greater than 0
#endif

//---- basic types

struct MeshNode {
	constexpr MeshNode(float x = 0, unsigned boundary = -1): m_X(x), m_BoundaryIndex(boundary) { }
	
	float m_X;
	unsigned m_BoundaryIndex;
	
	constexpr float Distance(MeshNode &n) {
		float r = m_X - n.m_X;
		if(r < 0)
			r *= -1.0f;
		return r;
	}
};

struct MeshElement {
	constexpr MeshElement(int a = 0, int b = 0) : m_Node1(a), m_Node2(b) { }

	int m_Node1;
	int m_Node2;
};

enum class BoundaryType {
	None,
	Flow,
	Convection,
};

struct Boundary {
	BoundaryType m_Type;
	
	static constexpr Boundary CreateConvection(float alpha, float envTemp) {
		return Boundary(BoundaryType::Convection, alpha, envTemp);
	}
	
	static constexpr Boundary CreateFlow(float q) {
		return Boundary(BoundaryType::Flow, q);
	}	
	struct Convection_t {
		float alpha;
		float envTemp;
	};
	struct Flow_t {
		float q;
	};
	
	Flow_t Flow;
	Convection_t Convection;

	constexpr Boundary(): m_Type(), Flow(), Convection() { }
private:
	constexpr Boundary(BoundaryType type, float alpha, float env) :
		m_Type(type), Flow(), Convection({alpha, env}) { }
	constexpr Boundary(BoundaryType type, float q) :
		m_Type(type), Flow({q}), Convection()  { }
};

//---- configuration

constexpr float Length = (float)(ROD_LENGTH);
constexpr float vk = (float)(HEAT_TRANSFER_COEFFICIENT);
constexpr float vS = (float)(CROSSECTION_AREA);

//---- code begins here

template<size_t ElementCount>
struct StaticProblem {
	enum {
		NodeCount = ElementCount + 1,
		BoundaryCount = 2,
	};

	constexpr StaticProblem():
			m_Elements(),
			m_Nodes(),
			m_Boundaries() {
			
		for(size_t i = 0; i < ElementCount; ++i) {
			m_Elements[i].m_Node1 = i;
			m_Elements[i].m_Node2 = i + 1;
		}
			
		for(size_t i = 0; i < NodeCount; ++i)
			m_Nodes[i] = MeshNode((float)i * (Length / (float)(ElementCount)), -1);		
			
		m_Nodes[0].m_BoundaryIndex = 0;
		m_Nodes[NodeCount - 1].m_BoundaryIndex = 1;
		
		m_Boundaries[0] = Boundary::CreateFlow((float)(HEAT_SOURCE_DENSITY));
		m_Boundaries[1] = Boundary::CreateConvection((float)(CONVECTION_COEFFICIENT), (float)(ENVIRONMENT_TEMPERATURE));
	}

	MeshElement m_Elements[ElementCount];
	MeshNode m_Nodes[NodeCount];
	Boundary m_Boundaries[BoundaryCount];
private:
};

template<size_t W, size_t H>
struct StaticMatrix {
	constexpr StaticMatrix(): m_T() { }
	constexpr StaticMatrix(const StaticMatrix<W,H> &other): m_T() { 
		for(size_t i = 0; i < W * H; ++i)
			m_T[i] = other.m_T[i];
	}
	
	constexpr void fill(float val) {
		for(size_t i = 0; i < W * H; ++i)
			m_T[i] = val;
	}
	
	constexpr float& at(int w, int h) { 
		return m_T[w + W * h]; 
	}
	float at(int w, int h) const { 
		return m_T[w + W * h]; 
	}
	
	void print() const {
		for(size_t i = 0; i < W; ++i){
			for(size_t j = 0; j < H; ++j)
				printf("%6.2f  ", at(j, i));
			std::cout << std::endl;
		}
	}	

#define xat(w, h) m_T[w + W * h]
	constexpr void Invert(StaticMatrix<W,H> &out) const {
		static_assert(W == H, "Invalid matrix size!");
		StaticMatrix<W,H> self(*this);
		out.fill(0.0f);
		for(size_t i = 0; i < W; i++)
		   out.at(i, i) = 1.0f;    
		   
		for(size_t k = 0; k < W; k++) {            
			float temp = self.xat(k, k);            
			for(size_t j = 0; j < W; j++) {
				self.xat(k, j) /= temp;    
				out.xat(k, j) /= temp;
			}                                

			for(size_t i=0; i<W; i++) {
				float temp = self.xat(i, k); 
				for(size_t j = 0; j < W; j++) {
					if(i == k) break;
					self.xat(i, j) -= self.xat(k, j) * temp;
					out.xat(i, j) -= out.xat(k, j) * temp;
				}
			}
		}
	}	
#undef xat
private:
	float m_T[W * H];
};

	
template<size_t aW, size_t aH, size_t W>
constexpr void MultiplyMatrices(StaticMatrix<aW, aH>& a, StaticMatrix<W, aH>& b, StaticMatrix<W, aH> &out) {
	for(size_t i = 0; i < W; ++i){
		for(size_t j = 0; j < aH; ++j){
			float res = 0;
			for(size_t k = 0; k < aW; ++k){
				res += b.at(i, k) * a.at(k, j);
			}
			out.at(i, j) = res;
		}
	}
}

template<size_t ElementCount>
struct StaticSolver {
	constexpr StaticSolver(): 
			m_Problem(),
			m_Matrix(),
			m_Inverted(),
			m_Solution(),	
			m_P() {
	
		solve();
	}
	
	constexpr void solve() {
		m_Matrix.fill(0.0f);
		m_Solution.fill(0.0f);
		m_P.fill(0.0f);

		for(size_t i = 0; i < ElementCount; ++i) {
			MeshElement &e = m_Problem.m_Elements[i];
			MeshNode &n1 = m_Problem.m_Nodes[e.m_Node1];
			MeshNode &n2 = m_Problem.m_Nodes[e.m_Node2];
			int p1 = e.m_Node1;
			int p2 = e.m_Node2;
		
			float dist = n1.Distance(n2);
			
			float C = vk * vS / dist;	
			if(n1.m_BoundaryIndex < Problem_t::BoundaryCount)
				ProcessBoundary(p1, &n1);
			if(n2.m_BoundaryIndex < Problem_t::BoundaryCount)
				ProcessBoundary(p2, &n2);				
			
			AddToGlobal(p1, p2, C);
			
		}
		
		m_Matrix.Invert(m_Inverted);
		MultiplyMatrices(m_Inverted, m_P, m_Solution);
	}	
	
	using Problem_t = StaticProblem<ElementCount>;
	Problem_t m_Problem;

	StaticMatrix<Problem_t::NodeCount, Problem_t::NodeCount> m_Matrix;
	StaticMatrix<Problem_t::NodeCount, Problem_t::NodeCount> m_Inverted;
	StaticMatrix<1, Problem_t::NodeCount> m_Solution;
	StaticMatrix<1, Problem_t::NodeCount> m_P;
	
private:
	constexpr void AddToGlobal(int p1, int p2, float C) {
		m_Matrix.at(p1, p1) += C;
		m_Matrix.at(p2, p2) += C;
		m_Matrix.at(p2, p1) += -C;
		m_Matrix.at(p1, p2) += -C;	
	}

	constexpr void ProcessBoundary(int nodeid, MeshNode *n) {
		auto *b = &m_Problem.m_Boundaries[n->m_BoundaryIndex];
		if(b->m_Type == BoundaryType::Flow) {
			double d = b->Flow.q * vS;
			m_P.at(0, nodeid) += -d;
			return;
		}
		if(b->m_Type == BoundaryType::Convection) {
			float d = b->Convection.alpha * vS;
			m_P.at(0, nodeid) += ( b->Convection.envTemp * d);
			m_Matrix.at(nodeid, nodeid) += d;
			return;
		}		
	}
};

constexpr StaticSolver<ELEMENT_COUNT> Solution = StaticSolver<ELEMENT_COUNT>();


#ifndef ELEMENT_COUNT
#define ELEMENT_COUNT 50
#endif

#ifndef ENVIRONMENT_TEMPERATURE
#define ENVIRONMENT_TEMPERATURE 40.0f
#endif

#ifndef HEAT_SOURCE_DENSITY
#define HEAT_SOURCE_DENSITY -150
#endif

#ifndef CONVECTION_COEFFICIENT
#define CONVECTION_COEFFICIENT 1.0f
#endif

#ifndef HEAT_TRANSFER_COEFFICIENT
#define HEAT_TRANSFER_COEFFICIENT 75.0f
#endif

#ifndef CROSSECTION_AREA
#define CROSSECTION_AREA 1.0f
#endif

#ifndef ROD_LENGTH
#define ROD_LENGTH 5.0f
#endif

int main() {
	std::cout << "Compile-time conditions:\n"
		<< "ENVIRONMENT_TEMPERATURE: " << ENVIRONMENT_TEMPERATURE << "\n"
		<< "HEAT_SOURCE_DENSITY: " << HEAT_SOURCE_DENSITY << "\n"
		<< "CONVECTION_COEFFICIENT: " << CONVECTION_COEFFICIENT << "\n"
		<< "HEAT_TRANSFER_COEFFICIENT: " << HEAT_TRANSFER_COEFFICIENT << "\n"
		<< "CROSSECTION_AREA: " << CROSSECTION_AREA << "\n"
		<< "ROD_LENGTH: " << ROD_LENGTH << "\n"
		<< "\n";
	
    std::cout << "Solution for " << ELEMENT_COUNT << " elements: \n";
    Solution.m_Solution.print();		
    return 0;
}