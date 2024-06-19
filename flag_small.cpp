/* 
Small flag algebra's program written by Bernard Lidicky

 Files and file formats:
 
 Problem:
    min/max c_1xG_1 + c_2xG_2 + c_3xG_3
 subject to
    c_4xG_4 + c_5xG_5 + c_6 \geq 0
    c_7xG_7 + c_8xG_8 + c_9 \geq 0
 
 Objective file: (F_edges2__objective.txt)
   c_1 G_1
   c_2 G_2
   c_3 G_3
 0
   c_6
   c_4 G_5
   c_5 G_5
 0
   c_9
   c_7 G_7
   c_8 G_8
 0
 Optional program description (that does not start with a number)
 
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <istream>
#include <iterator>
#include <vector>
#include <set>
#include <utility>
#include <assert.h>
#include <cstring>
#include <algorithm>
#include <cstdarg>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <limits>
#include <limits.h>
#include <unistd.h>
#include <iomanip>  
#include <errno.h>
#include <time.h>
#include <unordered_map>
#include <sys/types.h> 
#include <sys/wait.h> 
//#include <filesystem>

#ifdef _USING_OMP_
#include <omp.h>
#endif

using namespace std;


#ifndef G_COLORED_EDGES
#define G_COLORED_EDGES 2
#endif





#ifndef V
#define V      11  // Maximum number of vertices of flags. If you know you need smaller, use smaller - will be faster and less memory consuming
#endif

#define G_NOT_ALL_FLAGS_USED // Some extra warnings if there are unexpected flags. This should be enabled

#define G_USING_ZERO_AS_ANY_COLOR // If something has color 0, it is considered to have anycolor or being uncolored. Usefull when trying to find extensions of some graph

#define G_BE_BRAVE_AND_IGNORE_SAFETY_ASSERTS // Will be slower, but if not defined. But if you don't trust the program, comment it :-)

#define USE_REDUCED_TYPES // 4 4 1:123  and 4 4 1:234 are considered redundant and only one is used. This you want to have on. Not having it is stupid unless you have ordered vertices


#define MAX_FLAG_TYPES  300000 // Just a constant - raise if neded - means different theta



bool g_use_simple_linear_constraints = true;
vector<int>  g_use_simple_linear_constraints_list;

bool g_use_product_linear_constraints = false;
bool g_use_product_linear_constraints_by_expanding = false;
bool g_use_square_linear_constraints = false;
bool g_use_linear_constraints_self_products = false;



// Copy to +1 counterparts.... to make life easier
#ifdef G_COLORED_EDGES
#define COLORS_EDGES (G_COLORED_EDGES+1)
#endif


#ifndef G_BLOW_UP_COLOR_EDGES
#define G_BLOW_UP_COLOR_EDGES  1
#endif

int g_blow_up_color_edges[V]; // adds more precision in specifying the blow-up


#ifndef G_PRECISION
#define G_PRECISION 16
#endif





#ifdef DONT_USE_C11
double stod(const string&  str)
{
    stringstream ss(str);
    double d;
    ss >> d;
    return d;
}

long stol(const string&  str)
{
    stringstream ss(str);
    long l;
    ss >> l;
    return l;
}

string to_string (long long val)
{
    stringstream ss;
    ss << val;
    return ss.str();
}
#endif



// This is the main class that holds a flag. 
class flag
{
public:
    flag()
    {
		m_Theta = 0;
        m_Theta_class = 0;
		m_vertices = 0;


    }
    
    virtual ~flag()
    {
        
    }
    
    
    
    
    void make_all_vertices_labeled()
    {
        m_Theta = m_vertices;
    }
    
    bool is_labeled(const int v) const
    {
        assert(v < m_vertices);
        return  v < m_Theta;
    }
    
    // returns the number of labeled vertices
    int labeled_vertices_cnt() const
    {
        return m_Theta;
    }    

    
    // First m_Theta vertices are considered labeled in that order
    void set_Theta(int Theta, int Theta_class = 0)
    {
        // if !G_ORDERED_VERTICES, then the labeled vertices are first. Not true for ordered.
        //cerr << m_vertices << " " << Theta << endl;
        if (Theta > m_vertices)
        {
            cerr << "FAIL: Attempting to load grap on "<< m_vertices << " vertices with Theta=" << Theta << endl;
        }
        assert(Theta <= m_vertices);
        m_Theta = Theta;
        m_Theta_class = Theta_class;
    }
    
	
    int get_Theta()
    {
        return m_Theta;
    }
    
    void set_vertices_and_Theta(int vertices, int Theta, int Theta_class = 0)
    {
        set_vertices(vertices);
        set_Theta(Theta, Theta_class);
    }
	
	
    // set the nymber of vertices and clears the graph
    void set_vertices(int vertices)
    {
        assert(vertices <= V);
        m_vertices = vertices;
        

		
#ifdef G_COLORED_EDGES
        for (int u = 0; u < V; u++)
        {
	        for (int v = 0; v < V; v++)
	        {
	            m_color_edge[u][v] = 0;
	        }
	    }
        
        m_colored_edges[0] = e();
        for (int c = 1; c < COLORS_EDGES; c++)  m_colored_edges[c] = 0;
#endif

        
        

        
        
        
        
        

    }



    

#ifdef G_COLORED_EDGES
    void color_edge(int u, int v, int color)
    {

    	assert(u != v);
    	assert(u < m_vertices);
    	assert(v < m_vertices);
    	assert(0 <= u);
    	assert(0 <= v);
    	assert(color >= 0);
        assert(color < COLORS_EDGES);
        
    	int old_color = m_color_edge[u][v];
        
    	m_colored_edges[old_color]--;
    	
        m_color_edge[u][v] = m_color_edge[v][u] = color;
    	
    	m_colored_edges[color]++;

    }
#endif
    
 

#ifdef G_COLORED_EDGES
    // Makes a permutation of edges of this flag
    // Maybe be rewritten for higher effectivity
    void permute_edge_colors(const vector<int> &color_permutation)
    {
        assert (color_permutation.size() >= COLORS_EDGES);
        
        for (int u = 0; u < m_vertices-1; u++)
            for (int v = u+1; v < m_vertices; v++)
        {
            color_edge(u,v,color_permutation[m_color_edge[u][v]]);      
        }
    }
#endif 
    
    
    
// Sequence of functions for checking if two flags have the same type //
private:
    // TODO: This sequence could be rewritten by using early check and no need
    // for creating copies of F - could be more effective.
    bool have_same_type_colorblind_colored_3edges(const flag &h) const
    {
        return have_same_type_colorblind_colored_edges(h);
    }
    
    
    bool have_same_type_colorblind_colored_edges(const flag &h) const
    {
#ifdef G_COLORED_EDGES_BLIND
        flag F;
        for (int p = 0; p < (int)g_allowed_edgecolor_permutations.size(); p++)
        {
            F = h; 
            F.permute_edge_colors(g_allowed_edgecolor_permutations[p]);
            if (have_same_type_colorblind_oriented_edges(F)) return true;            
        }
        return false;
#else
        return have_same_type_colorblind_oriented_edges(h);
#endif
    }
    
    bool have_same_type_colorblind_oriented_edges(const flag &h) const
    {
        return have_same_type_colorblind_vertices(h);
    }
    
    bool have_same_type_colorblind_vertices(const flag &h) const
    {
        return have_same_type_leftright_blind(h);
    }


    bool have_same_type_leftright_blind(const flag&h) const
    {
        return have_same_type_rotation_reverse_blind(h); 
    }

    bool have_same_type_rotation_reverse_blind(const flag &h) const
    {
        return have_same_type_not_colorblind(h);
    }

    // See if the lists are identical up to a rotation
    bool same_lists_up_to_rotation(const int *L1, int L1sz, const int *L2, int L2sz) const
    {
        
        if (L1sz != L2sz)
        {
            return false;
        }

        if (L1sz == 0)
        {
            return true;
        }
        
        // Finds the smallest element in both lists
        int min1 = 0;
        int min2 = 0;
        for (int i = 1; i < L1sz; i++)
        {
            if (L1[i] < L1[min1]) min1 = i;
            if (L2[i] < L2[min2]) min2 = i;
        }
        
        // And now go from the min element and check all entries
        for (int i = 0; i < L1sz; i++)
        {
            if (L1[(i+min1)%L1sz] != L2[(i+min2)%L1sz]) return false;
        }
        
        return true;
    }

    bool have_same_type_not_colorblind(const flag &f) const
    {
        //if (m_Theta != f.m_Theta) return false;
        //if (m_Theta == 0) return true;
        
        
        
#ifdef G_COLORED_EDGES
        assert(m_Theta < V);
        for (int u = 0; u < m_Theta; u ++)
            for (int v = u+1; v < m_Theta; v++)
            {
#ifdef G_USING_ZERO_AS_ANY_COLOR
                if (m_color_edge[u][v] == 0 || f.m_color_edge[u][v] == 0) continue;
#endif            
                if (m_color_edge[u][v] != f.m_color_edge[u][v]) return false;
            }
#endif
        
        
        
        
        
        
        
        
        return true;
    }
    

public:    
    // This is the function to call in general program
    bool have_same_type(const flag &f) const
    {
        if (m_Theta_class != 0 && f.m_Theta_class != 0 && m_Theta_class != f.m_Theta_class) return false;
        if (m_Theta == 0 && f.m_Theta == 0) return true;
        if (labeled_vertices_cnt() != f.labeled_vertices_cnt()) return false;
        
        return have_same_type_colorblind_colored_3edges(f);
    }
    

    template <bool verbose_output=false> 
    bool is_isomorphic_to_colorblind_colored_3edges(const flag &h) const
    {
        return is_isomorphic_to_colorblind_colored_edges<verbose_output>(h);
    }
    
    template <bool verbose_output=false> 
    bool is_isomorphic_to_colorblind_colored_edges(const flag &h) const
    {
#ifdef G_COLORED_EDGES_BLIND
        flag F;
        for (int p = 0; p < (int)g_allowed_edgecolor_permutations.size(); p++)
        {
            // A noticable speed-up for big ones...
            bool good_permuation = true;
            for (int c = 0; c < COLORS_EDGES; c++)
            {
                if (m_colored_edges[g_allowed_edgecolor_permutations[p][c]] != h.m_colored_edges[c])
                {
                    good_permuation = false;
                    break;
                }
            }
            if (!good_permuation) continue;
            F = h; 
            F.permute_edge_colors(g_allowed_edgecolor_permutations[p]);
            if (is_isomorphic_to_colorblind_oriented_edges<verbose_output>(F)) 
            {
                return true;
            }
        }
        assert(g_allowed_edgecolor_permutations.size() != 0);
        return false;  
#else
        return is_isomorphic_to_colorblind_oriented_edges<verbose_output>(h);
#endif
    }
    
    template <bool verbose_output=false> 
    bool is_isomorphic_to_colorblind_oriented_edges(const flag &h) const
    {
        return is_isomorphic_to_colorblind_vertices<verbose_output>(h);
    }

    template <bool verbose_output=false> 
    bool is_isomorphic_to_colorblind_vertices(const flag &h) const
    {
        return is_isomorphic_to_reverseblind_leftright_system<verbose_output>(h);
    }

    template <bool verbose_output=false> 
    bool is_isomorphic_to_reverseblind_leftright_system(const flag &h) const
    {
        return is_isomorphic_to_reverseblind_rotation_system(h);
    }

    
    bool is_isomorphic_to_reverseblind_rotation_system(const flag &h) const
    {
        return is_isomorphic_to_not_colorblind(h);
    }

    template <bool verbose_output=false>    
    bool is_isomorphic_to(const flag &h) const
    {
        // Just a quick kill
        if (m_vertices != h.m_vertices) 
        {
            if (verbose_output) cerr << "The number of vertices is different" << endl;
            return false;
        }
        if (labeled_vertices_cnt() != h.labeled_vertices_cnt())
        {
            if (verbose_output) cerr << "The number of labeled vertices is different" << endl;
            return false;
        }
        
        if (m_Theta_class != 0 && h.m_Theta_class != 0  && m_Theta_class != h.m_Theta_class)
        {
            if (verbose_output) cerr << "The m_Theta_class is different" << endl;
            return false;
        } 

       
        
        
        //cerr << "Testing " << h.print() << " and " << print() << endl;
        
        return is_isomorphic_to_colorblind_colored_3edges<verbose_output>(h);
    }
	
	
	// perm is mapping vertices of h to vertices of this. pv is a candidate for
	// mapping of v, all vertices before v are already mapped.
	bool is_map_up_to_v_correct(int v, int pv, const int *perm, const flag &h) const
	{
			
						
#ifdef G_COLORED_EDGES
      // check edges, 0 as joker
      for (int u = 0; u < v; u++)
         if (m_color_edge[perm[u]][pv] != h.m_color_edge[u][v] && m_color_edge[perm[u]][pv] != 0 && h.m_color_edge[u][v] != 0) return false;        
#endif


        
        
        

        
        
        
		return true;
	}

    // perm is saying how are vertices of h mapped to *this
    // it means color of x in h is the same as color of perm[x] in *this
    //
    //  v in h is mapped to perm[v] in *this
    //
    bool is_mapping_an_equality(int *perm, const flag &h) const
    {
#ifdef G_COLORED_EDGES
        //	check coloring
        for (int u = 0; u < m_vertices; u++)
        {
            for (int x = u+1; x < m_vertices; x++)
            {
                if (m_color_edge[perm[u]][perm[x]] == 0 || h.m_color_edge[u][x]==0) continue;
                if (m_color_edge[perm[u]][perm[x]] != h.m_color_edge[u][x]) return false;
            }
        }
#endif
        
        
        
        

        
        
        
        
        

        return true;
    }
    
    // Real isomorphism testing
	bool make_perm_isomorphic(int v, int *perm, bool *used, const flag &h) const
	{
		// Test if the permutation is isomorphism
		if (v >= m_vertices)
		{
#ifndef G_BE_BRAVE_AND_IGNORE_SAFETY_ASSERTS
            assert(is_mapping_an_equality(perm,h) == true);
#endif
     
            
            
			return true;
		}

		for (int pv = m_Theta; pv < m_vertices; pv++)
		{
			// check if used
			if (used[pv]) continue;

			if (!is_map_up_to_v_correct(v, pv, perm, h)) continue;

            
			perm[v] = pv;
			used[pv]  = true;
			
			if (make_perm_isomorphic(v+1, perm, used, h)) return true;
			
			used[pv]  = false;
		}
		
		return false;
	}
	
	
    template <bool verbose_output=false> 
	bool is_isomorphic_to_not_colorblind(const flag &h) const
	{
        //if (labeled_vertices_cnt() > 0)
        //    cerr << "Starting testing not colorblind " << print() <<" and "<< h.print() << endl;

        // cerr << endl << "Testing " << endl << h.print() << endl << print() << endl;

        
        // quick check
        if (m_vertices != h.m_vertices) 
        {
            //if (labeled_vertices_cnt() > 0)
            //cerr << "Wrong number of vertices " << endl;
            return false;
        }
        if (labeled_vertices_cnt() != h.labeled_vertices_cnt())
        {
            //if (labeled_vertices_cnt() > 0)
            //cerr << "Wrong number of labeled vertices " << labeled_vertices_cnt() << " " << h.labeled_vertices_cnt() << endl;
            return false;
        }
        
		
		
		

		
#ifdef G_COLORED_EDGES
        for (int c = 0; c < COLORS_EDGES; c++)
        {
            //if (m_vertices==0)
            //    cerr << c << "  " << m_colored_edges[c] << " " << h.m_colored_edges[c] << endl;
	        if (m_colored_edges[c] >  h.m_colored_edges[c]+h.m_colored_edges[0]) return false;
            if (m_colored_edges[c]+h.m_colored_edges[0] <  h.m_colored_edges[c]) return false;
	    }
#endif

		

        // Check for same type
        if (!have_same_type_not_colorblind(h)) return false;

		// try check all permutations --- WASTING HERE!!!!!!!!
		int perm[V];
		bool used[V];

		// Theta is preserved by the mapping
        for (int i = 0; i < m_Theta; i++)
        {
            perm[i] = i;
            used[i] = true;
        }
        
        for (int i = m_Theta; i < V; i++)
        {
            perm[i] = 0;
            used[i] = false;
        }
		

        bool rv = make_perm_isomorphic(m_Theta, perm, used, h);
        if (verbose_output==true && rv==true)
        {
            cerr << "Isomorphic with permutation ";
            for (int i = 0; i < m_vertices; i++)
            { 
                cerr << i << "->" << perm[i] << " ";
            }
            cerr << endl;
            cerr << "This maps the second graph to the first one" << endl;
        }
        return rv;        
	}


    bool has_as_notinduced_subflag(const flag &h) const
    {
        // Not implemented for other
        assert(0);
        return false;
    }

    bool is_identical_to(const flag &h, bool write_why_not=true) const
    {
        if (m_vertices != h.m_vertices)
        {
            if (write_why_not)
                cerr << "Number of vertices is different: " << m_vertices << " vs " << h.m_vertices << endl;
            return false;
        }

        if (m_Theta != h.m_Theta)
        {
            if (write_why_not)
                cerr << "Thetas are different: " << m_Theta << " vs " << h.m_Theta  << endl;
            return false;
        }
        
        int perm[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
        assert(m_vertices < 20);
        return is_mapping_an_equality(perm, h);
    }
    
    
    
    void create_minlex_signature_blind_3edges(const flag &h, string &lex_min) const
    {
        string Hprint = h.print("");
        if (Hprint.compare(lex_min) < 0)
        {
            lex_min = Hprint;
        }
    }
    
    void create_minlex_signature_blind_edges(const flag &h, string &lex_min) const
    {
#ifdef G_COLORED_EDGES_BLIND
        flag F;
        for (int p = 0; p < (int)g_allowed_edgecolor_permutations.size(); p++)
        {
            F = h; 
            F.permute_edge_colors(g_allowed_edgecolor_permutations[p]);
            create_minlex_signature_blind_3edges(F, lex_min);            
        }
#else
        create_minlex_signature_blind_3edges(h, lex_min);
#endif
    }
    
    void create_minlex_signature_oriented_edges(const flag &h, string &lex_min) const
    {
        create_minlex_signature_blind_edges(h,lex_min);
    }
    
    
    void create_minlex_signature_blind_vertices(const flag &h, string &lex_min) const
    {
        create_minlex_signature_oriented_edges(h, lex_min);
    }

        
    void create_minlex_signature()
    {
    }


    void get_type_subflag(flag &f) const
    {
        const int mapping[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
        f.as_subflag(*this,mapping,m_Theta,m_Theta);        
        f.m_Theta_class = m_Theta_class;
    }
    
    
    // removes two vertices
    void remove_vertices(int u, int v)
    {
        int mapping[V];
        
        int id = 0;
        for (int i = 0; i < m_vertices;i++)
        {
            if (i == u || i == v) continue;
            mapping[id++] = i;
        }
        
        flag h;
        h.as_subflag(*this, mapping, m_vertices-2, 0);
        
        *this = h;
    }

    void remove_labeled_vertices()
    {
        if (m_Theta == 0)
            return;
        
        int mapping[V];
        
        for (int i = 0; i < m_vertices-m_Theta;i++)
        {
            mapping[i] = i+m_Theta;
        }
        
        flag h;
        h.as_subflag(*this, mapping, m_vertices-m_Theta, 0);
        
        *this = h;
    }



    void as_subflag(const flag &h,const int *mapping, int vertices, int theta, bool create_signature=true)
    {
        set_vertices_and_Theta(vertices,theta);
		
#ifdef G_COLORED_EDGES
        for (int x = 0; x < m_vertices-1; x++)
            for (int y = x+1; y < m_vertices; y++)
            {
                color_edge(x,y,h.m_color_edge[mapping[x]][mapping[y]]);
            }
#endif
        

        

        
        

		
        
        
        
    }
    

    void as_subflag_in_blowup(const flag &h,const int *mapping, int vertices, int theta, bool create_signature=true)
    {
        set_vertices_and_Theta(vertices,theta);
        
#ifdef G_COLORED_EDGES
        for (int x = 0; x < m_vertices-1; x++)
            for (int y = x+1; y < m_vertices; y++)
            {
                if (mapping[x] != mapping[y])
                {
                    color_edge(x,y,h.m_color_edge[mapping[x]][mapping[y]]);
                }
                else
                {
                    // color_edge(x,y,G_BLOW_UP_COLOR_EDGES);
                    color_edge(x,y,g_blow_up_color_edges[mapping[x]]);
                }
            }
#endif
        
        
        
    }
    
    
        
    int e()
    {
        return (m_vertices*(m_vertices-1))/2;
    }
    
    int e3()
    {
        return (m_vertices*(m_vertices-1)*(m_vertices-2))/6;
    }	
	
	
    bool contains_as_subflag_map_rest(const flag &g, int* mapping, int next_to_map, int first_available) const
    {
        if (next_to_map >= g.m_vertices)
        {
            flag h;
            h.as_subflag(*this, mapping,g.m_vertices, g.m_Theta);
            
            //This may happen!!
            // assert(g.m_Theta == 0);
         
            //cerr << "Testing " << h.print() << " and " << g.print() << endl;

            
            if (h.is_isomorphic_to(g))
            {
                return true;
            }
                         
            return false;
        }
        
        if (first_available < g.m_Theta)
        {
            mapping[next_to_map] = first_available;
            if (contains_as_subflag_map_rest(g,mapping,next_to_map+1,first_available+1)) return true;            
        }
        else
        {
        for (int i = first_available; i < m_vertices; i++)
        {
            mapping[next_to_map] = i;
            if (contains_as_subflag_map_rest(g,mapping,next_to_map+1,i+1)) return true;
        }
        }
        
        return false;
    }
    
    // Count density of g as a subflags
    bool contains_as_subflag(const flag &g) const
    {
        if (g.m_vertices > m_vertices) return false;
        
        int mapping[g.m_vertices];
                
        return contains_as_subflag_map_rest(g,mapping,0,0);        
    }
    
    
	void density_subflag_map_rest(const flag &g, int &total, int &good, int* mapping, int next_to_map, int first_available) const
	{
		if (next_to_map >= g.m_vertices)
		{
			flag h;
			h.as_subflag(*this, mapping,g.m_vertices, g.m_Theta);
			
			assert(g.m_Theta == 0);
			
			if (h.is_isomorphic_to(g)) good++;
			
			total++;
			return;
		}
		
		for (int i = first_available; i < m_vertices; i++)
		{
			mapping[next_to_map] = i;
			density_subflag_map_rest(g,total,good,mapping,next_to_map+1,i+1);
		}
	}
	
	// Count density of g as a subflags
    double density_subflag(const flag &g) const
    {
		if (g.m_vertices > m_vertices) return 0;
		
		int mapping[g.m_vertices];
		
		int total = 0;
		int good = 0;
		
		density_subflag_map_rest(g,total,good,mapping,0,0);
		
        return (double)good/(double)total;
	}

	int count_subflag(const flag &g) const
	{
		if (g.m_vertices > m_vertices) return 0;
		
		int mapping[g.m_vertices];
		
		int total = 0;
		int good = 0;
		
		density_subflag_map_rest(g,total,good,mapping,0,0);
		
      return good;
	}
    

    
    
    
    
    void count_labeled_copies_of_map_rest(const flag &g_labeled, int &good_maps, int* mapping, int *used, int next_to_map) const
    {
        if (next_to_map >= g_labeled.m_vertices)
        {
            flag h;
            h.as_subflag(*this, mapping,g_labeled.m_vertices, g_labeled.m_Theta);
            
            assert(g_labeled.m_Theta == g_labeled.m_vertices);
            
            if (h.have_same_type(g_labeled)) good_maps++;

            return;
        }
        
        for (int i = 0; i < m_vertices; i++)
        {
            if (used[i]) continue;
            mapping[next_to_map] = i;
            used[i] = 1;
            count_labeled_copies_of_map_rest(g_labeled,good_maps,mapping,used,next_to_map+1);
            used[i] = 0;
        }
    }
    

    int count_labeled_copies_of(const flag &g) const
    {
        if (g.m_vertices > m_vertices) return 0;
       
        // This does not work correctly for unlabeled graphs
        assert(g.m_Theta == 0 && m_Theta == 0);
        
        
        // We use have_same_type to check correct mappings...
        flag g_labeled = g;
        g_labeled.m_Theta = g.m_vertices;
        
        int mapping[g.m_vertices];
        int used[m_vertices];
        for (int i = 0; i < m_vertices; i++)
        {
            used[i] = 0;
        }
        
        int good_maps = 0;
        
        count_labeled_copies_of_map_rest(g_labeled,good_maps,mapping,used,0);
        
        return good_maps;
    }
    
    
    
    void generate_subflags_of_size_n_rest(int n, vector<flag> &subflags, int* mapping, int next_to_map, int first_available, const vector<int> &available_to_map) const
    {
        //cerr << "Runnig n=" << n << " mapping=" << mapping[0] << "," << mapping[1] << "," << mapping[2] << " next_to_map=" << next_to_map << endl; 

        if (next_to_map >= n)
        {
            flag h;
            h.as_subflag(*this, mapping, n, 0);
            bool h_is_new = true;
            for (int i = 0; i < (int)subflags.size(); ++i)
            {
                if (h.is_isomorphic_to(subflags[i])) 
                {
                    h_is_new = false;
                    break;
                } 
            }
            if (h_is_new) subflags.push_back(h);            
            return;
        }

        //cerr << "first_available=" << first_available  << " (n-next_to_map)=" << (n-next_to_map) << endl;
        //cerr << "(int)available_to_map.size()-(n-next_to_map)=" << (int)available_to_map.size()-(n-next_to_map) << endl;
        for (int i = first_available; i < (int)available_to_map.size()-(n-next_to_map-1); i++)
        {
            //cerr << "i=" << i << endl;
            mapping[next_to_map] = available_to_map[i];
            generate_subflags_of_size_n_rest(n,subflags,mapping,next_to_map+1,i+1,available_to_map);
        }
    }
   
    
    // Could be written as a template with size of in_subflag - should be faster
    // an/or with arrays
    // if in_subrgaph is specified, then only subflags containing in_subrgaph are generated
    int generate_subflags_of_size_n(int n, vector<flag> &subflags, const vector<int> &in_subflag = vector<int>()) const
    {
        if (n > m_vertices) return 0;
        
        // This does not work correctly for unlabeled graphs
        assert(m_Theta == 0);
                
        int mapping[n];
        
        vector<int> available_to_map;
        for (int i = 0; i < m_vertices; i++)
            if (find(in_subflag.begin(), in_subflag.end(), i) == in_subflag.end()) 
                available_to_map.push_back(i);
        
        for (int i = 0; i < (int)in_subflag.size(); i++)
            mapping[i] = in_subflag[i];
        
        //cerr << "Working on it... " << available_to_map.size() << endl;
        //cerr << "Starting with " << (int)in_subflag.size() << endl;

        generate_subflags_of_size_n_rest(n, subflags, mapping, (int)in_subflag.size(),0,available_to_map);
        
        return (int)subflags.size();
    }
    
    
    
    
    // Copy data from g if the number of vertices of g is <= current number
    void copy_from(const flag &g)
    {
        assert(m_vertices >= g.m_vertices);
     
        

        



        

        
        
#ifdef G_COLORED_EDGES
        for (int x = 0; x < g.m_vertices-1; x++)
            for (int y = x+1; y < g.m_vertices; y++)
            {
                color_edge(x,y,g.m_color_edge[x][y]);
            }        
#endif        
    }
    
/***************************************************************************************************************** Printing ****************/
	void print_for_human() const
	{
		// Prints the number of vertices, zero as not type and then part of the diagonal matrix
		cout << m_vertices << " " << m_Theta << " ";
		
		
		cout << endl;
	}

    template<typename T>
    string print_latex(bool use_label, const T &graph_label, bool color_1_nonedge = false)
    {
        if (m_vertices == 0)
        {
            return "";
        }
        
        stringstream ss;
        
        // Prints the number of vertices, zero as not type and then part of the diagonal matrix
        ss << "\n\\begin{tikzpicture}[flag_pic]";
        ss << "\\outercycle{" << m_vertices << "}{" << m_Theta << "}\n";

        const string edgestr = "--";        
        
#ifdef G_COLORED_EDGES
        for (int u = 0; u < m_vertices; u++)
        {
            for (int v = u+1; v < m_vertices; v++)
            {
                ss << "\\draw[edge_color" << m_color_edge[u][v] << "] (x"<<u<<")" << edgestr << "(x"<<v<<");"; 
            }
            ss << "  ";
        }
        ss << endl;
#else  // Edges have no colors
#endif
        
      
        for (int i = 0; i < m_vertices; i++)
        {
            if (i < m_Theta)
                ss << "\\draw (x"<<i<<") node[labeled_vertex]{};";
            else
                ss << "\\draw (x"<<i<<") node[unlabeled_vertex]{};";
        }
        ss << endl;
        for (int i = 0; i < m_vertices; i++)
        {
            ss << "\\labelvertex{" << i << "}";
        }
        
        
        
        


        
        
        if (use_label)
        {
            ss << "\\draw (labelpoint) node{"<<graph_label<<"};" << endl;
        }
        ss << "\\end{tikzpicture} " << endl;
        
        return ss.str();
    }
    
    
	
	string print(const string &delimeter=" ") const
	{
		stringstream ss;
		
		// Prints the number of vertices, zero as not type and then part of the diagonal matrix
        ss << m_vertices << delimeter;
        if (m_Theta_class != 0)
            ss << m_Theta << "." <<  m_Theta_class << delimeter << delimeter;
        else 
            ss << m_Theta << delimeter << delimeter;
        
        

        

#ifdef G_COLORED_EDGES
  		for (int u = 0; u < m_vertices; u++)
  		{
    		for (int v = u+1; v < m_vertices; v++)
				ss << delimeter << m_color_edge[u][v];
            ss << delimeter;
		}
#endif
                
		return ss.str();
	}
    
	
	
    
/***************************************************************************************************************** Loading ****************/	
	void load_from_string(const char *str)
	{
        if (strcmp(str,"cin") == 0)
        {
            load_from_stream(cin, -1, -1);            
        }
        else
        {
            stringstream s(str);
            load_from_stream(s, -1, -1);
        }
	}
		
	
   bool load_from_stream(istream &stream, int assumed_vertices, int assumed_theta)
   {
        int vertices, theta, theta_class=0; //, theta;
        //string theta_dot_class; //  (theta loaded as string since it may be  THETA.CLASS)
		stream >> vertices;

		if (stream.fail()) return false;
        
        stream >> theta;

		if (stream.fail()) return false;
       
       /*
//       ss << m_Theta << "." <<  m_Theta_class << delimeter << delimeter;
       string only_theta="", only_theta_class="";
       
       std::size_t found_dot = theta_dot_class.find(".");
       if (found_dot != std::string::npos)
       {
           only_theta = theta_dot_class.substr(0,found_dot);
           only_theta_class = theta_dot_class.substr(found_dot+1);
#ifndef G_FLAG_PRODUCTS_SLOW_LOW_MEMORY
           cerr << "The program was not compiled with support for repeated flag types.\n";
           cerr << "Use -DG_FLAG_PRODUCTS_SLOW_LOW_MEMORY when compiling the program." << endl;
           assert(0);
#endif           
       }
       else
       {
           only_theta = theta_dot_class;
           only_theta_class = "";
       }
       theta = stol(only_theta);
       if (only_theta_class.length() > 0)
       {
           theta_class = stol(only_theta_class);
       }
       */

        if (assumed_theta != -1 && theta != assumed_theta)
        {
            cerr << "Err: assumed theta " << assumed_theta << " but loaded " <<  theta << " for " << vertices << " vertices " << endl;
        }      
		assert(assumed_theta == -1 || theta == assumed_theta);
        assert(assumed_vertices == -1 || assumed_vertices == vertices);
		set_vertices_and_Theta(vertices, theta, theta_class);
       
       //cerr << "Set " <<vertices << " " << theta << " " << theta_class << endl;
		
		
#ifdef G_COLORED_EDGES
		int color;
		for (int u = 0; u < vertices; u++)
		{
			for (int v = u+1; v < vertices; v++)
			{
				stream >> color;
				color_edge(u,v,color);
			}
		}
#endif




       


       

       



       
       if (!stream) 
       {
           cerr << "Stream failed while loading flag " << print() << endl;
           assert(0);
       }  


    //std::cout << " good()=" << stream.good() << endl;
    //std::cout << " eof()=" << stream.eof() << endl;
    //std::cout << " fail()=" << stream.fail() << endl;
    //std::cout << " bad()=" << stream.bad() << endl;       
       //cerr << print() << endl;
       
       create_minlex_signature();

		return true;
	}
	
    


	

	  

public:
    
    flag( const flag& f)
    {
        *this = f;
    }

    flag& operator=(const flag& f)
    {
        m_vertices = f.m_vertices;
        m_Theta    = f.m_Theta;
        m_Theta_class = f.m_Theta_class;


        
        
        
#ifdef G_COLORED_EDGES
        for (int i = 0; i < m_vertices; i++)
            for (int j = 0; j < m_vertices; j++)
                {
                    m_color_edge[i][j] = f.m_color_edge[i][j];
                }
        
        for (int i = 0; i < COLORS_EDGES; i++)
        {
            m_colored_edges[i] = f.m_colored_edges[i];
        }         
        
#endif
        
        
        
        
        

        
        return *this;
    }
    
		
	int m_vertices;  // number of vertices
    int m_Theta;     // number of labeled vertices in case of not ordered. Binary string which are labeled in ordered
    int m_Theta_class; // Optional thing - same theta can have several classes - this should help with doing fancier rounding.
                       // m_Theta_class == 0 means any class. If not specified, any class is the default


//	int m_id; // id of the original flag

    
	
	
#ifdef G_COLORED_EDGES
	int m_color_edge[V][V];  // adjacency matrix
//	int m_colored_deg[V][COLORS_EDGES]; // number of edges from each vertex of a particular color
	int m_colored_edges[COLORS_EDGES]; // number of edges of each color
//	vector<int> m_deg_sorted;
//	bool m_deg_sorted_valid;
#endif
    


};


class flag_and_coefficient
{
public:

    bool operator==(const flag_and_coefficient& fc) const
    {
        if (fc.coefficient != coefficient) return false;
        return g.is_isomorphic_to(fc.g);
    }

	flag   g;
	double coefficient;
};

// Try if the number should be actuallly rounded
double smart_round(double d, double precision)
{
    double dr = round(d);
    if (abs(dr-d) < precision) return dr;
    return d;
}

double g_smart_round_precision = 0.00000001;
double smart_round(double d)
{
    return smart_round(d, g_smart_round_precision);
}



std::ostream& operator<< (std::ostream& stream, const flag_and_coefficient& fc)
{
    if (fc.coefficient == 0) return stream;
    stream.precision(G_PRECISION);
    stream << smart_round(fc.coefficient) << "  " << fc.g.print() << endl;

    return stream;
}




// constraint is of the form...  m_flag + m_constant >= 0
class linear_constraint
{
public:
    linear_constraint()
    {
        m_checked = false; 
        m_labeled_vertices_in_type_cnt = 0;

        m_required_coefficients_sum_for_elcp_constraints = false;
        m_required_coefficients_sum_for_elcp_constraints_value = 0;
    }
    
    
    bool check_constraint()
    {
        if (m_entries.size() == 0) 
        {
            m_checked = false;
            //m_same_types = false;
            return false;
        }
        
        m_entries_max_size = m_entries[0].g.m_vertices;
        
        for (int i = 1; i < (int)m_entries.size(); i++)
        {
            if (! m_entries[0].g.have_same_type(m_entries[i].g))
            {
                m_checked = false;
                return false;
            }
            if (m_entries_max_size < m_entries[i].g.m_vertices)
            {
                m_entries_max_size = m_entries[i].g.m_vertices;
            }
        }
        
        m_labeled_vertices_in_type_cnt = m_entries[0].g.labeled_vertices_cnt();
        
        m_entries[0].g.get_type_subflag(m_type);
        
        m_checked = true;

//        m_same_types = true;
        
        if (m_constant != 0)
        {
            flag_and_coefficient fc;
            fc.g = m_type;
            fc.coefficient = m_constant;
            m_entries.push_back(fc);
            m_constant = 0;
        }
        
        return true;
    }
    
    string print()
    {
        stringstream ss;
        ss << "0 " << m_constant << endl;
        for (int i = 0; i < (int)m_entries.size(); i++)
        {
            ss.precision(G_PRECISION);
            ss << m_entries[i].coefficient << " " << m_entries[i].g.print() << endl;
        }
        return ss.str();
    }
    
    void add_entry(const flag_and_coefficient& fc)
    {
        for (int i = 0; i < (int)m_entries.size(); i++)
        {
            if (m_entries[i].g.is_isomorphic_to(fc.g))
            {
                m_entries[i].coefficient += fc.coefficient;
                if (m_entries[i].coefficient == 0)
                {
                    m_entries.erase(m_entries.begin()+i);
                }
                return;
            }
        }
        m_entries.push_back(fc);
    }
    
    bool operator==(const linear_constraint &lc) const
    {
        if (lc.m_labeled_vertices_in_type_cnt != m_labeled_vertices_in_type_cnt) return false;
        if (lc.m_entries.size() != m_entries.size()) return false;
        if (lc.m_constant != m_constant) return false;
        if (lc.m_type.is_isomorphic_to(m_type) == false) return false;
        if (lc.m_entries_max_size != m_entries_max_size) return false;

        for (flag_and_coefficient entry : m_entries)
        {
            // Test if it is the same as some other entry in lc
            bool found_match = false;
            for (flag_and_coefficient entry2 : lc.m_entries)
            {
                if (entry == entry2)
                {
                    found_match = true;
                    break;
                }
            }
            if (found_match == false) return false;
        }

        return true;
    }

    bool is_identical_after_type_permutation(const linear_constraint &lc) const
    {
        if (lc == *this) return true;

        /*
        if (lc.m_labeled_vertices_in_type_cnt != m_labeled_vertices_in_type_cnt) return false;
        if (lc.m_entries.size() != m_entries.size()) return false;
        if (lc.m_constant != m_constant) return false;
        if (lc.m_type.is_isomorphic_to(m_type) == false) return false;
        if (lc.m_entries_max_size != m_entries_max_size) return false;
        if (m_entries.size() == 0) return true;
         */
 
        
        int permutation[m_type.m_vertices];
        for (int i = 0; i < m_type.m_vertices; i++) permutation[i] = i;

        // This will be used to permute vertices of flags
        int permutation_fc[V];
        for (int i = 0; i < V; i++) permutation_fc[i] = i;


        flag permuted_type;

        // Skip the identity permutation, already tested at the beginning
        while ( std::next_permutation(permutation,permutation+m_type.m_vertices) )
        {
            permuted_type.as_subflag(m_type, permutation, m_type.m_vertices, m_type.m_Theta);
            if (!m_type.is_isomorphic_to(permuted_type)) continue;

            // Now we got an automorphism of the type
            // so we permute the whole constraint - it is wasting a little but will
            // be easier to read
            linear_constraint lc_tmp;
            
            for (const flag_and_coefficient &fc : m_entries)
            {
                // make permutation of type, the rest is identity
                // from the beginning
                for (int i = 0; i < m_type.m_vertices; i++)
                    permutation_fc[i] = permutation[i];

                flag_and_coefficient new_fc;
                new_fc.coefficient = fc.coefficient;
                new_fc.g.as_subflag(fc.g, permutation_fc, fc.g.m_vertices, fc.g.m_Theta);

                lc_tmp.add_entry(new_fc);
            }

            assert(lc_tmp.check_constraint());

            if (lc == lc_tmp)
                return true;

        }

        return false;
    }

public:
	vector<flag_and_coefficient> m_entries;
	double m_constant;
    
    bool   m_checked;
    //bool   m_same_types;  // true if the entries are of same type and size 
    int    m_entries_max_size;
    int    m_labeled_vertices_in_type_cnt;

    // This is a special bonus thing
    bool    m_required_coefficients_sum_for_elcp_constraints;
    double  m_required_coefficients_sum_for_elcp_constraints_value;

    flag   m_type;
};

std::ostream& operator<< (std::ostream& stream, const linear_constraint& lc)
{
    stream << "0 ";

    // Constant is made into a labeled flag, so this tries to recover it
    // if possible. Maybe there should be a switch for this?
    bool constant_found = false;
    if (lc.m_constant == 0)
    {
        for (const auto & fc :  lc.m_entries)
        {
            if (fc.g.m_vertices == fc.g.m_Theta)
            {
                assert(constant_found == false);
                stream.precision(G_PRECISION);
                stream << smart_round(fc.coefficient) << endl;
                constant_found = true;
            }
        }
    }
    if (!constant_found)
    {
        stream.precision(G_PRECISION);
        stream << smart_round(lc.m_constant) << endl;
    }

    
    for (const auto & fc :  lc.m_entries)
    {
        if (fc.g.m_vertices != fc.g.m_Theta)
        {
            stream << fc;
        }
    }
    return stream;
}




vector<flag> g_unlabeled_flags[V+1];
vector<flag> g_types[V+1]; // types - it is just same as g_unlabeled_flags, but all vertices are labeled :-) Used for reduced types
vector<vector<flag> >  g_flags; // flags to process using multiplication - every type&size is a separate list
vector<vector<flag> >  g_flag_square_linear_constraints; // flag used only for linear constraints


vector<linear_constraint> g_linear_constraints;
vector<flag_and_coefficient> g_objective_combination;
vector<flag_and_coefficient> g_objective_divisor; // may be nothing...
vector<flag> g_forbidden_subflags;
vector<flag> g_forbidden_subflags_by_size[V+1];
vector<flag> g_no_slack_flags; // means no slack in CSDP
string       g_program_description;  // loaded from the objective function
vector<string> g_additional_csdp_blocks;




// When using the program for crossing number, it may happen that there is a fixed
// number of vertices in each color (for the -n ? we try to compute). Then it makes
// no sense to use some flags that will 'clearly' never be used because their sqaure
// is never used. This helps avoid such situations.


string time_to_str(time_t time_taken)
{

    int days_t  = time_taken/60/60/24;
    int hours_t = (time_taken%(60*60*24))/60/60;
    int min_t = (time_taken%(60*60))/60;
    int sec_t = time_taken%(60);

    stringstream ss;
    ss << days_t << "d+"
         << hours_t << "h+" 
         << min_t << "m+" 
         << sec_t << "s";

    return ss.str();
}

class mini_timer
{
public:
    mini_timer()
    {
        m_time_start = 0;
    }

    void start()
    {
        m_time_start = time(NULL);
    }

    string report(int processing_id, int total_ids)
    {
        if (m_time_start == 0)
        {
            m_time_start = time(NULL);
            return "";
        }

        time_t time_taken = time(NULL) - m_time_start;

        #ifdef _USING_OMP_
            int adjust = omp_get_num_threads();
        #else
            int adjust = 1;
        #endif

        stringstream ss;
        if (time_taken > 0 && processing_id > adjust)
        {
            time_t total_time  = (time_taken*total_ids)/(processing_id-adjust);
            ss << " took " << time_to_str(time_taken);
            ss << " will take " << time_to_str(total_time-time_taken);
            return ss.str();
        }

        return "";
    } 
    
private:
    time_t m_time_start;
};


int binomial(int n, int k) {
    int b = 1;
    for (int i = 0; i < k; ++i) {
        b *= (n - i);
        b /= (i + 1);
    }
    return b;
}

int factorial(int n) {
    int b = 1;
    for (int i = 2; i <= n; i++) {
        b *= i;
    }
    return b;
}



string filename_prefix()
{
    stringstream filename;
    
    filename << "F"
#ifdef G_COLORED_EDGES
    << "_edges" << COLORS_EDGES-1
#endif
#ifdef G_COLORED_EDGES_BLIND
    << "blind"
#endif
    ;
    return filename.str();
}



inline int find_flag_in_list_nonparalel(const flag &f, const vector<flag> &fv)
{
    for (int i = 0; i < (int)fv.size(); ++i)
    {
        if (f.is_isomorphic_to(fv[i]))
        {
            return i;
        }
    }
    return -1;
}


inline int find_flag_in_list(const flag &f, const vector<flag> &fv)
{
    // Not worth paralelizing the search.
    if (fv.size() < 500)
    {
        return find_flag_in_list_nonparalel(f,fv);
    }

    // attempt to paralelize the search if it was not parallel yet
#ifdef _USING_OMP_
    if (omp_in_parallel() == false)
    {        
        //cerr << "Trying in parallel" << endl;
        volatile int location = -1;
        
        #pragma omp parallel for shared(location)    
        for (int i = 0; i < (int)fv.size(); ++i)
        {
            if (location != -1) continue;
            if (f.is_isomorphic_to(fv[i]))
            {
                location = i;
            }
        }
        return location;
    }
#endif

    return find_flag_in_list_nonparalel(f,fv);
}

inline int get_flag_type_in_list(const flag& f, vector< vector<flag> > &flags)
{
#ifdef USE_REDUCED_TYPES
    const int identity[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30};
    flag f_type;
    f_type.as_subflag(f,identity,f.m_Theta,0);
#endif 
    
    // Find correct type
    for (int i = 0; i < (int)flags.size(); i++)
    {
        if (flags[i].size() == 0)
        {
            cerr << "Some bug" << endl;
        }
        if (f.have_same_type(flags[i][0])) return i;
#ifdef USE_REDUCED_TYPES
        if (f.m_Theta == flags[i][0].m_Theta && f.m_Theta_class == flags[i][0].m_Theta_class)
        {
            flag g_flags_type;
            g_flags_type.as_subflag(flags[i][0],identity,flags[i][0].m_Theta,0);
            if (g_flags_type.is_isomorphic_to(f_type))
                return -2;
        }
#endif        
        
    }
    
    return -1;
}

inline void include_flag_in_list(const flag& f, vector< vector<flag> > &flags, int id = -3, bool ignore_reduced_types = true)
{
    if (id == -3)
    {
        id = get_flag_type_in_list(f, flags);
    }
#ifdef  USE_REDUCED_TYPES   
    if (id == -2 && ignore_reduced_types) return;
#endif   
    if (id == -1 || id == -2)
    {
        vector<flag> new_type_list;
        new_type_list.push_back(f);
        flags.push_back(new_type_list);
        return;
    }
    
    flags[id].push_back(f);
}


inline void include_flag_in_list_if_new(const flag& f, vector< vector<flag> > &flags, bool ignore_reduced_types = true)
{
    int id = get_flag_type_in_list(f, flags);
#ifdef USE_REDUCED_TYPES
    if (id == -2 && ignore_reduced_types) 
    {
        return;
    }
#endif
    if (id == -1 || id == -2)
    {
//        cerr << "new type" << endl;
        include_flag_in_list(f, flags, id,ignore_reduced_types);
        return;
    }
    
    if (find_flag_in_list(f,flags[id]) != -1)
    {
//        cerr << "duplicate" << endl;
        return;   
    }
//    cerr << "new for known type" << f.print() << endl;
    flags[id].push_back(f);
}

inline bool is_flag_forbidden_noninduced(const flag &g)
{
    for (int i = 0; i < (int)g_forbidden_subflags.size(); i++)
    {
        if (g.has_as_notinduced_subflag(g_forbidden_subflags[i])) return true;
    }
    return false;
}


inline bool is_flag_forbidden(const flag &g, int verbose_output=0)
{
    for (int i = 0; i < (int)g_forbidden_subflags.size(); i++)
    {
        //cerr << "Testing " << g.print() << " and " << g_forbidden_subflags[i].print() << endl;
        if (g.contains_as_subflag(g_forbidden_subflags[i])) 
        {
            if (verbose_output > 0)
            {
                cerr << "Forbidding " << g.print() << " because of " << g_forbidden_subflags[i].print() << endl;
            }
            return true;
        }
    }
    return false;
}

//#ifdef G_USE_PERMITTED_SUBFLAGS  not used in early generation - maybe should be... ?
inline bool is_flag_forbidden_avoid_zeros(const flag &g, const vector<int> &in_subrgaph = vector<int>())
{
    vector<flag> subflags_list;
    for (int i = 1; i <= g.m_vertices; i++)
    {
        //if (g_permitted_subflags[i].size() == 0) continue;        
        if (g_forbidden_subflags_by_size[i].size() == 0) continue;
        
        subflags_list.clear();
        g.generate_subflags_of_size_n(i, subflags_list, in_subrgaph);
                
        for (int j = 0; j < (int)subflags_list.size(); j++)
        {
// Do not check subgraphs that are not fully colored
#ifdef G_COLORED_EDGES            
            if (subflags_list[j].m_colored_edges[0] > 0) continue;
#endif            
            if (find_flag_in_list(subflags_list[j], g_forbidden_subflags_by_size[i]) != -1) return true;
        }
    }
    return false;
}


inline bool g_already_in_known_flags(flag &g, vector<flag> &flag_list)
{
    for (unsigned int i = 0; i < flag_list.size(); i++)
    {
        if (flag_list[i].is_isomorphic_to(g))
        {
            return true;
        }
    }
    return false;
}

inline void add_g_to_known_flags(flag &g, vector<flag> &flag_list)
{
    flag_list.push_back(g);
    //if (g.m_3edges_cnt >= 20)
    {
    //    cerr << g.print() << endl;
    }
}

inline void add_g_to_flags_list_if_new(flag &g, vector<flag> &flag_list)
{
    if (g_already_in_known_flags(g, flag_list)) return;
    add_g_to_known_flags(g, flag_list);
}


#ifdef G_COLORED_EDGES
void try_color_edge(flag &g, int u, int v, vector<flag> &flag_list)
{
    if (v >= g.m_vertices || v == u)
    {
        u++;
        v = u+1;
    }
    if (v >= g.m_vertices)
    {
        if (is_flag_forbidden(g)) return;
        
        g.create_minlex_signature();
        
        if (!g_already_in_known_flags(g,flag_list)) add_g_to_known_flags(g,flag_list);
        //		add_g_to_known_flags(flagh_list);
        return;
    }
    
    
    // Color only uncolored edges.....
    if (g.m_color_edge[u][v] == 0)
    {
// Either color 1 us used for blow-up or for tournaments it is not used at all
// since color 1 does not change orientation
        for (int color = 1; color < COLORS_EDGES; color++)
        {
            g.color_edge(u,v,color);
            try_color_edge(g,u,v+1,flag_list);
        }
        g.color_edge(u,v,0);
    }
    else
    {
        try_color_edge(g,u,v+1,flag_list);
    }
}
#endif








void extensions_of_g_edges(flag &g, vector<flag> &flag_list)
{
    //cerr << "Generating extensions of " << g.print() << endl;
    
#ifdef G_COLORED_EDGES
    try_color_edge(g,0,1,flag_list);
    return;
#endif


}





void extensions_of_g(flag &g, vector<flag> &flag_list)
{
    extensions_of_g_edges(g,flag_list);
}



// u is original, v will be a copy
// edges containing both will have color 1
void duplicate_vertex_u_to_v(flag &g, int u, int v)
{

    
#ifdef G_COLORED_EDGES
    
    for (int x = 0; x < g.m_vertices-1; x++)
    {
        if (x == u) 
            g.color_edge(v,x,1);
        else
            g.color_edge(v,x,g.m_color_edge[u][x]);
    }
#endif


    
//    assert(0);
}




///////////////////////////////////////////////////////////////////////////////////////////////////////////////// P(F1,F2,H)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////// P(F1,F2,H)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////// P(F1,F2,H)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////// P(F1,F2,H)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////// P(F1,F2,H)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////// P(F1,F2,H)


// Rest of F2 is also piskced as a subset
void pick_F2_mapping(const flag &F1, const flag &F2, const flag &H, int *mapping_F2, int next_to_map,  bool *used_from_H, int &good_maps)
{
    if (next_to_map >= F2.m_vertices)
    {
        flag H_F2;
        H_F2.as_subflag(H, mapping_F2, F2.m_vertices, F2.m_Theta);
        if (!H_F2.is_isomorphic_to(F2)) return;
		
        good_maps++;
        
        return;
    }
    
    int min_vertex = 0;
    if (next_to_map > F2.m_Theta) min_vertex = mapping_F2[next_to_map-1]+1;
    
    for (int v = min_vertex; v < H.m_vertices; v++)
    {
        if (used_from_H[v]) continue;
        
        used_from_H[v] = true;
        mapping_F2[next_to_map] = v;
        pick_F2_mapping(F1,F2,H,mapping_F2,next_to_map+1,used_from_H,good_maps);
        used_from_H[v] = false;
    }
    
}

// Rest of F1 is picked only as a subset
void pick_F1_mapping(const flag &F1, const flag &F2, const flag &H, int *mapping_F1, int *mapping_F2, int next_to_map,  bool *used_from_H, int &good_maps)
{
    if (next_to_map >= F1.m_vertices)
    {
		//cout << "." << endl;

        // Check if the mapping is corect!
        flag H_F1;
        H_F1.as_subflag(H, mapping_F1, F1.m_vertices,F1.m_Theta);
		
        if (!H_F1.is_isomorphic_to(F1)) return;
        
        pick_F2_mapping(F1,F2,H,mapping_F2,F2.m_Theta,used_from_H,good_maps);
        return;
    }
    
    int min_vertex = 0;
    if (next_to_map > F1.m_Theta) min_vertex = mapping_F1[next_to_map-1]+1;
    
    for (int v = min_vertex; v < H.m_vertices; v++)
    {
        if (used_from_H[v]) continue;
        
        used_from_H[v] = true;
        mapping_F1[next_to_map] = v;
        pick_F1_mapping(F1,F2,H,mapping_F1,mapping_F2,next_to_map+1,used_from_H,good_maps);
        used_from_H[v] = false;
    }
}

// Theta mapping is picked as all permutations...
void pick_theta_mapping(const flag &F1, const flag &F2, const flag &H, int *mapping_theta, int next_to_map, bool *used_from_H, int &good_maps)
{
    if (next_to_map >= F1.m_Theta)
    {
        int mapping_F1[F1.m_vertices];
        int mapping_F2[F2.m_vertices];
        for (int i = 0; i < F1.m_Theta; i++)
        {
            mapping_F2[i] = mapping_F1[i] = mapping_theta[i];
        }
				
        pick_F1_mapping(F1,F2,H,mapping_F1,mapping_F2, F1.m_Theta, used_from_H,good_maps);
        return;
    }
	
    for (int v = 0; v < H.m_vertices; v++)
    {
        if (used_from_H[v]) continue;
        
        if (!H.is_map_up_to_v_correct(next_to_map, v, mapping_theta, F1)) continue;
       
         used_from_H[v] = true;
         mapping_theta[next_to_map] = v;
         pick_theta_mapping(F1,F2,H,mapping_theta,next_to_map+1,used_from_H,good_maps);
		 used_from_H[v] = false;
    }
}

inline double P_F1_F2_IN_H(const flag &F1, const flag &F2, const flag &H, bool return_only_good_maps=false)
{
    assert(F1.m_Theta == F2.m_Theta);

#if defined(G_COLORED_EDGES) && !defined(G_COLORED_EDGES_BLIND)
	for (int c = 1; c < COLORS_EDGES; c++)
	{
        //		cout << "e" << endl;
        if (F1.m_colored_edges[c] > H.m_colored_edges[c]) return 0;
        if (F2.m_colored_edges[c] > H.m_colored_edges[c]) return 0;
	}
#endif

    

    
    
    
	
		
    int mapping_theta[F1.m_Theta];
    bool used_from_H[H.m_vertices];
    
    for (int i = 0; i < H.m_vertices; i++)
    {
        used_from_H[i] = false;
    }
	
    int good_maps = 0;
    pick_theta_mapping(F1,F2,H,mapping_theta,0,used_from_H,good_maps);
	
    if (return_only_good_maps) return good_maps;
    
    int all_maps = 1;
    for (int i = 0; i < F1.m_Theta; i++) all_maps *= (H.m_vertices-i);
    all_maps *= binomial(H.m_vertices-F1.m_Theta,F1.m_vertices-F1.m_Theta);
    all_maps *= binomial(H.m_vertices-F1.m_vertices,F2.m_vertices-F1.m_Theta);
    
	//    if (good_maps > 0)
	//        cout << "good_maps/all_maps: " << good_maps << "/" << all_maps << endl;
    
    
    return (double)good_maps/(double)all_maps;
}



inline double P_F1_F2_F3_IN_H(const flag &F1, const flag &F2, const flag &F3, const flag &H)
{
    assert(F1.m_Theta == F2.m_Theta);
    assert(F1.m_Theta == F3.m_Theta);
    
    assert(F1.m_vertices+F2.m_vertices+F3.m_vertices - 2*F1.labeled_vertices_cnt() <= H.m_vertices);
    
    flag type;
    F1.get_type_subflag(type);    

    
    // Some easy checks...
#if defined(G_COLORED_EDGES) && !defined(G_COLORED_EDGES_BLIND)
    for (int c = 1; c < COLORS_EDGES; c++)
    {
        if (F1.m_colored_edges[c] + F2.m_colored_edges[c] + F3.m_colored_edges[c] - 2*type.m_colored_edges[c] > H.m_colored_edges[c]) return 0;
    }
#endif
    
    
    
    


    
    


    

    int good_maps = 0;
    int all_maps = 0;
    
    int vertex_split[V];

    int t = 0;
    for (int i = 0; i < type.m_vertices; i++) vertex_split[t++] = i;
    for (int i = type.m_vertices; i < F1.m_vertices; i++) vertex_split[t++] = -1;
    for (int i = type.m_vertices; i < F2.m_vertices; i++) vertex_split[t++] = -2;
    for (int i = type.m_vertices; i < F3.m_vertices; i++) vertex_split[t++] = -3;
    for (;t < H.m_vertices; t++) vertex_split[t] = -4;
    
    std::sort (vertex_split,vertex_split+H.m_vertices);
    
    flag tmp;
    int mapping[V];
    int next_to_map; // temporary
    do {
        all_maps++;
        
        
        //cerr << "split: ";
        //for (int i = 0; i < H.m_vertices; i++)
        //    cerr << vertex_split[i] << " " ;
        //cerr << endl;        
        
        
        // check if Theta correct
        for (int j = 0; j < H.m_vertices; j++)
        {
            if (vertex_split[j] < 0) continue;
            if (vertex_split[j] < type.m_vertices) mapping[vertex_split[j]] = j;
        }
        tmp.as_subflag(H,mapping,type.m_vertices,type.m_Theta);
        //cerr << tmp.print() << " vs " << type.print();
        if (!tmp.is_isomorphic_to(type)) continue;

        
        // check if F1 correct
        next_to_map = type.m_vertices;
        for (int j = 0; j < H.m_vertices; j++)  if (vertex_split[j] == -1) mapping[next_to_map++] = j;
        tmp.as_subflag(H,mapping,F1.m_vertices,type.m_Theta);
        if (!tmp.is_isomorphic_to(F1)) continue;


        // check if F2 correct
        next_to_map = type.m_vertices;
        for (int j = 0; j < H.m_vertices; j++)  if (vertex_split[j] == -2) mapping[next_to_map++] = j;
        tmp.as_subflag(H,mapping,F2.m_vertices,type.m_Theta);
        if (!tmp.is_isomorphic_to(F2)) continue;


        
        // check if F3 correct
        next_to_map = type.m_vertices;
        for (int j = 0; j < H.m_vertices; j++)  if (vertex_split[j] == -3) mapping[next_to_map++] = j;
        tmp.as_subflag(H,mapping,F3.m_vertices,type.m_Theta);
        if (!tmp.is_isomorphic_to(F3)) continue;
        
        good_maps++;
        
    } while ( std::next_permutation(vertex_split,vertex_split+H.m_vertices) ); 

    //cerr << "Good/All: " << good_maps << "/" << all_maps << endl;
    
    return (double)good_maps/(double)all_maps;
}


inline double P_F1_F2_IN_labeled_H(const flag &F1, const flag &F2, const flag &H)
{

    //cerr << F1.m_Theta << " " << F2.m_Theta << endl;
    assert(F1.m_Theta == F2.m_Theta);
    if (F1.m_vertices+F2.m_vertices - F1.labeled_vertices_cnt() > H.m_vertices)
    {
        cerr << F1.m_vertices<< "+" << F2.m_vertices << "-" << F1.labeled_vertices_cnt() << " > " << H.m_vertices << endl;
    }
    assert(F1.m_vertices+F2.m_vertices - F1.labeled_vertices_cnt() <= H.m_vertices);
    assert(F1.labeled_vertices_cnt() >= H.labeled_vertices_cnt());
    
    flag type;
    F1.get_type_subflag(type);    
    
    if (H.labeled_vertices_cnt() > 0)
    {
        flag typeH;
        H.get_type_subflag(typeH);    
        if (!type.contains_as_subflag(typeH)) return 0;
    }
    
    // Some easy checks...
#if defined(G_COLORED_EDGES) && !defined(G_COLORED_EDGES_BLIND)
    for (int c = 1; c < COLORS_EDGES; c++)
    {
        if (F1.m_colored_edges[c] + F2.m_colored_edges[c] - type.m_colored_edges[c] > H.m_colored_edges[c]) return 0;
    }
#endif
    
    
    
    int good_maps = 0;
    int all_maps = 0;
    
    int vertex_split[V];
    
    int t = 0;
    for (int i = 0; i < type.m_vertices; i++) vertex_split[t++] = i;
    for (int i = type.m_vertices; i < F1.m_vertices; i++) vertex_split[t++] = -1;
    for (int i = type.m_vertices; i < F2.m_vertices; i++) vertex_split[t++] = -2;
    //for (int i = type.m_vertices; i < F3.m_vertices; i++) vertex_split[t++] = -3;
    for (;t < H.m_vertices; t++) vertex_split[t] = -4;
    
    std::sort (vertex_split+H.labeled_vertices_cnt(),vertex_split+H.m_vertices);
    
    flag tmp;
    int mapping[V];
    int next_to_map; // temporary
    do {
        all_maps++;
        
        
        
        // check if Theta correct
        for (int j = 0; j < H.m_vertices; j++)
        {
            if (vertex_split[j] < 0) continue;
            if (vertex_split[j] < type.m_vertices) mapping[vertex_split[j]] = j;
        }
        tmp.as_subflag(H,mapping,type.m_vertices,type.m_Theta);
        //cerr << tmp.print() << " vs " << type.print();
        if (!tmp.is_isomorphic_to(type)) continue;
        
        
        // check if F1 correct
        next_to_map = type.m_vertices;
        for (int j = 0; j < H.m_vertices; j++)  if (vertex_split[j] == -1) mapping[next_to_map++] = j;
        tmp.as_subflag(H,mapping,F1.m_vertices,type.m_Theta);
        if (!tmp.is_isomorphic_to(F1)) continue;
        
        
        // check if F2 correct
        next_to_map = type.m_vertices;
        for (int j = 0; j < H.m_vertices; j++)  if (vertex_split[j] == -2) mapping[next_to_map++] = j;
        tmp.as_subflag(H,mapping,F2.m_vertices,type.m_Theta);
        if (!tmp.is_isomorphic_to(F2)) continue;
        
        
        good_maps++;
        
    } while ( std::next_permutation(vertex_split+H.labeled_vertices_cnt(),vertex_split+H.m_vertices) ); 
    
    //cerr << "Good/All: " << good_maps << "/" << all_maps << endl;
    
    return (double)good_maps/(double)all_maps;
}



inline double P_F1_IN_H(const flag &F1, const flag &H, bool density=true)
{
    
    assert(F1.m_vertices <= H.m_vertices);


    flag type;
    F1.get_type_subflag(type);    

    if (H.labeled_vertices_cnt() > 0)
    {
        flag typeH;
        H.get_type_subflag(typeH);    
        if (!type.contains_as_subflag(typeH)) return 0;
    }    
    
    // Some easy checks...
#if defined(G_COLORED_EDGES) && !defined(G_COLORED_EDGES_BLIND)
        
    for (int c = 1; c < COLORS_EDGES; c++)
    {
        if (F1.m_colored_edges[c] > H.m_colored_edges[c]) return 0;
    }
#endif
    
    
    
    
    
    
    
    

    int good_maps = 0;
    int all_maps = 0;
    
    int vertex_split[V];
    
    int t = 0;
    for (int i = 0; i < type.m_vertices; i++) vertex_split[t++] = i;
    for (int i = type.m_vertices; i < F1.m_vertices; i++) vertex_split[t++] = -1;
    for (;t < H.m_vertices; t++) vertex_split[t] = -4;
    
    std::sort (vertex_split+H.labeled_vertices_cnt(),vertex_split+H.m_vertices);
    
    flag tmp;
    int mapping[V];
    int next_to_map; // temporary
    do {
        all_maps++;
        
        
        // check if Theta correct
        for (int j = 0; j < H.m_vertices; j++)
        {
            if (vertex_split[j] < 0) continue;
            if (vertex_split[j] < type.m_vertices) mapping[vertex_split[j]] = j;
        }
        tmp.as_subflag(H,mapping,type.m_vertices,type.m_Theta);
        //cerr << tmp.print() << " vs " << type.print();
        if (!tmp.is_isomorphic_to(type)) continue;
        
        
        // check if F1 correct
        next_to_map = type.m_vertices;
        for (int j = 0; j < H.m_vertices; j++)  if (vertex_split[j] == -1) mapping[next_to_map++] = j;
        tmp.as_subflag(H,mapping,F1.m_vertices,type.m_Theta);
        if (!tmp.is_isomorphic_to(F1)) continue;
        
        good_maps++;
        
    } while ( std::next_permutation(vertex_split+H.labeled_vertices_cnt(),vertex_split+H.m_vertices) ); 
    
    //cerr << "Good/All: " << good_maps << "/" << all_maps << endl;
    if (density)
    {
        if (all_maps == 0)
        {
            cerr << "all_maps is zero" << endl;
            assert(0);
        }
        return (double)good_maps/(double)all_maps;
    }
    else
        return good_maps;    
}


//  The following is by James Kanze
// http://gabisoft.free.fr/articles/fltrsbf1.html
// https://openclassrooms.com/forum/sujet/fluxtampon-gestion-31091?page=1
// https://lists.boost.org/Archives/boost/att-49459/fltrsbf1.htm
struct FilteringInputStreambuf : public std::streambuf
{
    FilteringInputStreambuf(streambuf * source, bool deleteWhenFinished = false) : mySource(source), myDeleteWhenFinished(deleteWhenFinished)
    {
    }
    
    virtual ~FilteringInputStreambuf()
    {
        resetSource(NULL);
    }
    
    void resetSource(streambuf * newSource, bool deleteWhenFinished = false)
    {
        sync();
        
        if (myDeleteWhenFinished)
        delete mySource;
        
        mySource = newSource;
        myDeleteWhenFinished = deleteWhenFinished;
        setg(NULL, NULL, NULL);
    }
    
    virtual int underflow()
    {
        int result(EOF);
        
        if ( gptr() < egptr() )
        result = static_cast<unsigned char>(*gptr());
        else if ( mySource != NULL )
        {
            result = extract(*mySource);
            
            if (result != EOF)
            {
                assert(result >= 0); // && result <= UCHAR_MAX); // does not work on older machines
                myPushbackBuffer = result;
                setg(&myPushbackBuffer, &myPushbackBuffer, &myPushbackBuffer + 1);
            }
        }
        
        return result;
    }
    virtual int sync()
    {
        int result(EOF);
        
        if (mySource != NULL)
        {
            if (gptr() == egptr() || mySource->sputbackc(*gptr()) != EOF)
            result = mySource->pubsync();
            
            setg( NULL, NULL, NULL ) ;
        }
        
        return result;
    }
    virtual std::streambuf * setbuf(char * buffer, std::streamsize length)
    {
        return mySource == NULL
        ?   static_cast<std::streambuf*>(NULL)
        :   mySource->pubsetbuf(buffer, length);
    }
    
    virtual int extract(std::streambuf & source)
    {
        int ch( source.sbumpc() ) ;
        if ( ch == '#' )
        {
            while ( ch != EOF && ch != '\n' && ch != '\r' )
            {
                ch = source.sbumpc() ;
            }
        }
        if (ch >= 126)
        {
            cerr << "WARNING: Reading from input: " << ch << " which is '" << char(ch) << "'" << endl;
        }
        return ch ;
    }
    
    std::streambuf * mySource;
    char myPushbackBuffer;
    bool myDeleteWhenFinished;    
} ;


class FilteringIstream 
:   private FilteringInputStreambuf
,   public istream
{
    public:

    FilteringIstream( istream& source ):
    FilteringInputStreambuf(source.rdbuf()), istream(this)
    {        
    }
    virtual             ~FilteringIstream(){} ;    
    virtual FilteringInputStreambuf* rdbuf() 
    {
        return this;
    }
     void changeSource(streambuf * newSource)
     {
         resetSource(newSource);
     }
} ;


bool starts_with_flag_expression(string filename)
{
    stringstream ss(filename);

    char ch;
    ss >> ch;
    
    return isdigit(ch) || ch=='-' || ch == '(' || ch == '.' || ch == '{';

    // Below is an old version, we need ( for flag calculator
    double number;
    ss >> number;
    return !ss.fail();
}




// This is a macro that creates  istream *istr variable that opened the filename in 
// a smart and nice way. 
//
//     istream *istr = &std::cin; 
#define OPEN_FILE_SMARTLY_BASE(istr, filename, fail_operation) \
    FilteringIstream fscin_XXX(std::cin); \
    istream *istr = &fscin_XXX;\
    stringstream ss_filename_XXX(filename);\
    ifstream infile_XXX; \
    bool file_exists_XXX = false; \
    if (filename != "cin") \
    { \
        infile_XXX.open (filename.c_str(), ifstream::in); \
        if (!infile_XXX.good()) \
        { \
            std::cerr << "Failed opening file " << filename << endl; \
            if (starts_with_flag_expression(filename)) \
            { \
                std::cerr << "Using it as the string itself" << endl; \
                istr = &ss_filename_XXX; \
            } \
            else \
            { \
                std::cerr << "Failed interpreting as flags without comments" << endl; \
                fail_operation; \
            } \
        } \
        else \
        { \
            istr = &infile_XXX; \
            file_exists_XXX = true; \
        } \
    } \
    FilteringIstream infileFilter_XXX(infile_XXX); \
    if (filename != "cin" && file_exists_XXX) \
    { \
        istr = &infileFilter_XXX; \
    } 


#define OPEN_FILE_SMARTLY_RETURN_FALSE_ON_FAIL(istr, filename)  OPEN_FILE_SMARTLY_BASE(istr, filename, return false)
#define OPEN_FILE_SMARTLY(istr, filename)  OPEN_FILE_SMARTLY_BASE(istr, filename, assert(0))

bool load_flags_from_file(string filename, vector<flag> &flag_list,int verbose_output=0)
{
    OPEN_FILE_SMARTLY_RETURN_FALSE_ON_FAIL(istr, filename);
    
 
    if (verbose_output)
        cerr << "Loading labeled flags from file " << filename << endl;
    
    
    flag f;
    while (f.load_from_stream((*istr),-1,-1))
    {
        flag_list.push_back(f);
    }
    
    ///infile.close();
    
    return true;    
}

bool dump_flags_to_file(string filename, vector<flag> &flag_list)
{
    ofstream outfile;
    outfile.open(filename.c_str(), ofstream::out);
    if (!outfile.good())
    {
        cerr << "Failed opening file " << filename << endl;
        return false;
    }
    
    cerr << "Writing flags to file " << filename << endl;
    
    for (unsigned int x = 0; x < flag_list.size(); x++)
    {
        outfile << flag_list[x].print() << endl;
    }
    
    outfile.close();
    
    return true;
}



bool load_flags_and_coefficients_from_file(string filename, vector<flag_and_coefficient> &flag_list, int verbose_output=0)
{

    OPEN_FILE_SMARTLY_RETURN_FALSE_ON_FAIL(istr, filename);

    if (verbose_output)
        cerr << "Loading labeled flags from file " << filename << endl;
    
    
    flag_and_coefficient fc;
    (*istr) >> fc.coefficient;
    cerr << "Loaded coefficient " << fc.coefficient << endl;
    while ((*istr) && fc.g.load_from_stream((*istr),-1,-1))
    {
        flag_list.push_back(fc);
        (*istr) >> fc.coefficient;
    }
    
    //infile.close();
    
    return true;
}

void dump_flag_and_coefficient(const flag_and_coefficient &fc, bool use_smart_round = false)
{
        if (fc.coefficient == 0) return;
        
        cout.precision(G_PRECISION);
        if (use_smart_round)
        {
            if (smart_round(fc.coefficient) == 0) return;
            cout << smart_round(fc.coefficient);
        }
        else
            cout << fc.coefficient;
        cout << "  " << fc.g.print() << endl;
}

void dump_flags_and_coefficients(vector<flag_and_coefficient> &flag_list, bool use_smart_round = false)
{
    for (unsigned int x = 0; x < flag_list.size(); x++)
    {
        dump_flag_and_coefficient(flag_list[x], use_smart_round);

    }
}



// final_flags means that the flags will not be used any further
// for genertin other flags. Important when combined with G_USE_FIRST_EDGE_COLOR_FOR_BLOWUP_ONLY
bool load_unlabeled_flags_from_file(int sizeKn, bool final_flags = true, bool remove_duplicates_while_loading=false, bool remove_forbidden_wile_loading=false)
{
    stringstream filename;
    
    filename << filename_prefix() << "__n" << sizeKn << "_unlabeled.txt";
    
    ifstream infile;
    infile.open (filename.str().c_str(), ifstream::in);
    if (!infile.good())
    {
        cerr << "Failed opening file with unlabeled flags " << filename.str() << endl;
        return false;
    }
    
    cerr << "Loading unlabeled flags from file " << filename.str() << endl;
    
    // For large instances, it really helps to reserve the space in containers
    //  before loading... so we read the file twice
    int sizes_cnt[V], duplicates_cnt[V], forbidden_cnt[V];
    for (int i = 0; i < V; i++)
    {
        sizes_cnt[i] = 0;
        duplicates_cnt[i] = 0;
        forbidden_cnt[i] = 0;
    }
    flag g;
    while (g.load_from_stream(infile,-1,0))
    {
        sizes_cnt[g.m_vertices]++;
    }
    for (int i = 0; i < V; i++) g_unlabeled_flags[i].reserve(sizes_cnt[i]);
    
    //    infile.seek(0);
    infile.close();
    infile.open (filename.str().c_str(), ifstream::in);
    
    
    int tested=0;
    
    while (g.load_from_stream(infile,-1,0))
    {
        
        tested++;
        assert(sizeKn >= g.m_vertices);
        
        if (remove_forbidden_wile_loading && is_flag_forbidden(g))
        {
            forbidden_cnt[g.m_vertices]++;
            continue;
        }
        
        if (remove_duplicates_while_loading)
        {
            volatile bool new_flag = true;
            
            #pragma omp parallel for shared(new_flag)
            for (int i = 0; i < (int)g_unlabeled_flags[g.m_vertices].size();i++)
            {
                
                if (new_flag == false)
                    continue;
                
                if (g.is_isomorphic_to(g_unlabeled_flags[g.m_vertices][i]))
                {
                        #pragma omp atomic write
                    new_flag = false;
                    // This allows to break when a duplicate found - may be usefull to check if the
                    // duplicates are indeed duplicates or if it is a mistake in the program....
                    //if (duplicates_cnt[g.m_vertices] == 2)
                    //{
                    //    cerr << "Found duplicates" << endl;
                    //    cout << g.print() << endl;
                    //    cout << g_unlabeled_flags[g.m_vertices][i].print() << endl;
                    //    exit(0);
                    //}
#ifndef _USING_OMP_                   
                    break;
#endif                    
                }
            }
            if (!new_flag)
            {
                duplicates_cnt[g.m_vertices]++;
                continue;
            }
        }
        g_unlabeled_flags[g.m_vertices].push_back(g);
        //g.reverse_rotation_system();
        //g_unlabeled_flags[g.m_vertices].push_back(g);
        //cerr << g.print() << " " << g_unlabeled_flags[g.m_vertices].size() << " tested " << tested << endl;
    }
    
    infile.close();
    
    for (int i = 0; i <= sizeKn; i++)
    {
        cerr << "Loaded # of unlabeled flags of size " << i << " is " << g_unlabeled_flags[i].size();
        if (remove_forbidden_wile_loading)
        {
            cerr << " and " << forbidden_cnt[i] << " forbiddens";            
        }
        if (remove_duplicates_while_loading)
        {
            cerr << " and " << duplicates_cnt[i] << " duplicates";
        }
        cerr << endl;
    }
    
    
    return true;
}



bool dump_unlabeled_flags(int sizeKn)
{
    stringstream filename;
    filename << filename_prefix() << "__n" << sizeKn << "_unlabeled.txt";
    
    ofstream outfile;
    outfile.open (filename.str().c_str(), ofstream::out);
    if (!outfile.good())
    {
        cerr << "Failed opening file " << filename.str() << endl;
        return false;
    }
    
    cerr << "Writing unlabeled flags to file " << filename.str() << endl;
    
    for (int f = 0; f < V; f++)
    {
        for (unsigned int x = 0; x < g_unlabeled_flags[f].size(); x++)
        {
            outfile << g_unlabeled_flags[f][x].print() << endl;
        }
    }
    
    outfile.close();
    
    return true;
}


// Creates all subflags of flags of size sizeKn
void generate_all_unlabeled_subflags_from_size(int sizeKn, int verbose_output)
{
    cerr << "Generating smaller unlabeled flags from size " << sizeKn << endl;

    int mapping[V];

    for (int n = sizeKn; n > 0; n--)
    {
        // Delete the ones we are just about to start generating
        g_unlabeled_flags[n-1].clear();
        
        // graph to be processed
        for (int id = 0; id < (int)g_unlabeled_flags[n].size(); id++)
        {
            // vertex to be skipped
            for (int skip=0; skip < n; skip++)
            {
                int mapID = 0;
                for (int u = 0; u < n; u++)
                {
                    if (u == skip) continue;
                    mapping[mapID++] = u;                    
                }
                
                flag F;
                F.as_subflag(g_unlabeled_flags[n][id], mapping, n-1, 0);
                
                volatile bool new_flag = true;
                
                   #pragma omp parallel for shared(new_flag)
                for (int i = 0; i < (int)g_unlabeled_flags[F.m_vertices].size();i++)
                {
                    if (new_flag == false) continue;
                    
                    if (F.is_isomorphic_to(g_unlabeled_flags[F.m_vertices][i]))
                    {
                        #pragma omp atomic write                        
                        new_flag = false;
                        #ifndef _USING_OMP_                   
                        break;
                        #endif 
                    }
                }
                if (!new_flag) continue;
                
                g_unlabeled_flags[F.m_vertices].push_back(F);
            } 
        }
        if (verbose_output)
        {
            cerr << "# of unlabeled flags of size " << n-1 << " is " << g_unlabeled_flags[n-1].size() << endl;
        }
    }
}




void try_extensions_of_g_to_last_vertex_edges(flag &g, vector<flag> &flag_list)
{
    
    
    int v = g.m_vertices;
#ifdef G_COLORED_EDGES
    try_color_edge(g,v-1,0,flag_list);
    return;
#endif




    v++; // this is to avoid warning of unused v
    assert(0);
}


void try_extensions_of_g_to_last_vertex(flag &g, vector<flag> &flag_list)
{
    try_extensions_of_g_to_last_vertex_edges(g, flag_list);
}




void dump_vector_of_flags(vector <flag> *flag_lists, int previous_size, int sizeKn)
{
    stringstream filename;
    filename << filename_prefix() << "__n" << sizeKn << "_unlabeled_dump.txt";
    
    ofstream outfile;
    outfile.open (filename.str().c_str(), ofstream::out);
    if (!outfile.good())
    {
        cerr << "Failed opening file " << filename.str() << endl;
        return;
    }
    
    // dumping of flag lists...
    for (int i = 0; i < previous_size; i++)
    {
        for (int f = 0; f < (int)flag_lists[i].size(); f++)
        {
            outfile << flag_lists[i][f].print() << endl;
        }
        outfile << endl;
    }
    
    outfile.close();
}



// The following function megres dest_flags with src_flags and clears src_flags
void merge_flags_lists_parallel(vector<flag> &dest_flags, vector<flag> &src_flags)
{
    
    //    cout << flag_list.size() << endl;
    
    vector<flag> new_flags;
    new_flags.reserve(src_flags.size());
    
#pragma omp parallel for
    for (int k = 0; k < (int)src_flags.size(); k++)
    {
        if (!g_already_in_known_flags(src_flags[k],dest_flags))
        {
#pragma omp critical (merging_lists)
            {
                new_flags.push_back(src_flags[k]);
            }
        }
    }
    //#ifdef DONT_USE_C11
    dest_flags.insert(dest_flags.end(), new_flags.begin(), new_flags.end() );
    src_flags.clear();
    vector<flag>().swap( src_flags ); // free the memory
    //#else
    //    dest_flags.insert(dest_flags.end(), std::make_move_iterator(new_flags.begin()), std::make_move_iterator(new_flags.end()) );
    //    src_flags.clear();
    //    src_flags.shrink_to_fit(); // free the memory - does not work well on icc
    //#endif  
}


// The following function megres dest_flags with src_flags and clears src_flags
template< class data_type >
void merge_vectors_parallel(vector<data_type> &dest_vector, vector<data_type> &src_vector)
{    
    vector<data_type> new_vector;
    new_vector.reserve(src_vector.size());
    
#pragma omp parallel for
    for (int k = 0; k < (int)src_vector.size(); k++)
    {
        bool already_exists = false;
        for (const auto& test : dest_vector)
        {
            if (test == src_vector[k])
            {
               already_exists = true; 
            }
        }
        if (already_exists == false)
        {
#pragma omp critical (merging_lists)
            {
                new_vector.push_back(src_vector[k]);
            }
        }
    }
    //#ifdef DONT_USE_C11
    dest_vector.insert(dest_vector.end(), new_vector.begin(), new_vector.end() );
    src_vector.clear();
    vector<data_type>().swap( src_vector ); // free the memory
    //#else
    //    dest_flags.insert(dest_flags.end(), std::make_move_iterator(new_flags.begin()), std::make_move_iterator(new_flags.end()) );
    //    src_flags.clear();
    //    src_flags.shrink_to_fit(); // free the memory - does not work well on icc
    //#endif  
}

//#define CROSSING_EXTENSIONS_COLOR 1

// We expect v has highest index of all
void try_crossings_extension_of_g_add_rest(flag &g, vector<flag> &flag_list, int v, int u, int x, int y)
{
    if (u >= v)
    {
        // test, we generated it all
        return;
    } 

    if (x >= v)
    {
        try_crossings_extension_of_g_add_rest(g, flag_list, v, u+1, 1, 0);
        return;
    }

    if (x == u)
    {
        try_crossings_extension_of_g_add_rest(g, flag_list, v, u, x+1, 0);
        return;
    }

    if (x == u)
    {
        try_crossings_extension_of_g_add_rest(g, flag_list, v, u, x+1, 0);
        return;
    }

    //TODO!!!
    assert(0);
}

 
 
void generate_unlabeled_flags_of_size(int i, int verbose_output, bool dump_unlabeled_while_generating)
{

    
    int previous_size = (int)g_unlabeled_flags[i-1].size();
    
    //vector <flag> flag_lists[previous_size]; // does not work on mac (llvm :-( )
    vector <flag> *flag_lists;
    flag_lists = new vector<flag>[previous_size];
    if (flag_lists == NULL)
    {
        cerr << "Memory allocation failed" << endl;
        exit(1);
    }
    
    if ( (int)g_unlabeled_flags[i-1].size() < 1 )
    {
        cerr << "No flags of size " << i-1 << endl;
        cerr << "Generating flags of size " << i << " failed." << endl;
        exit(0);
    }
    
    //            #pragma omp parallel for
#pragma omp parallel for ordered schedule(dynamic)
    for (int j = 0; j < (int)g_unlabeled_flags[i-1].size(); j++)
    {
        flag_lists[j].reserve(previous_size);

        
        flag g;
        g.set_vertices(i);
        g.copy_from(g_unlabeled_flags[i-1][j]);
        //vector <flag> valid_extensions_of_g;
        //cout << "going for " << j << endl;
        try_extensions_of_g_to_last_vertex(g, flag_lists[j]);
        if (verbose_output)
        {
            cerr << "generated size " << i << ": " << j  << "/" << (int)g_unlabeled_flags[i-1].size() << " found " << flag_lists[j].size() << " flags " << endl;
        }
    }
    if (verbose_output)
    {
        cerr << "Flags generated. Merging part starts." << endl;
    }
    
    if (dump_unlabeled_while_generating)
        dump_vector_of_flags(flag_lists, previous_size, i);
    
    for (int power = 1; power < previous_size; power *= 2)
    {
        if (verbose_output)
        {
            cerr << "Trying " << power << " in " << previous_size << endl;
        }
        for (int j = 0; j < previous_size; j += 2*power)
        {
            if (j+power < previous_size)
            {
                if (verbose_output)
                {
                    cerr << "Merging  " << j << "<-" << j+power << "  of sizes " <<  flag_lists[j].size() << " " << flag_lists[j+power].size() << endl;
                }
                merge_flags_lists_parallel(flag_lists[j],flag_lists[j+power]);
                //flag_lists[j+power].clear();
            }
        }
        if (dump_unlabeled_while_generating)
            dump_vector_of_flags(flag_lists, previous_size, i);
    }
    
    g_unlabeled_flags[i].swap(flag_lists[0]);
    
    delete[] flag_lists;
}


void get_unlabeled_flags_of_size(int Kn, bool force_generate_flags = false, bool remove_duplicates_while_loading = false, bool remove_forbidden_wile_loading = false, int verbose_output = 0, bool dump_unlabeled_while_generating= false)
{
    // Already loaded
    if (g_unlabeled_flags[Kn].size() > 0)
    {
        return;
    }
    
    // Getting unlabeled flags...
    if (force_generate_flags || !load_unlabeled_flags_from_file(Kn,true,remove_duplicates_while_loading, remove_forbidden_wile_loading))
    {
        // Try to reuse previous graphs...
        int loaded_size = 0;
        
        // blow-upping and generating would not work correctly...
        if (!force_generate_flags)
        {
            loaded_size = Kn-1;
            cerr << "Trying to load file with smaller flags..." << endl;
            while (loaded_size > 0 && !load_unlabeled_flags_from_file(loaded_size,false)) loaded_size--;
        }
        
        if (loaded_size == 0)
        {
            cerr << "Unable to load any smaller graphs. Starting from empty graph." << endl;
            flag g_zero;
            g_zero.set_vertices_and_Theta(0,0);
            g_unlabeled_flags[0].push_back(g_zero);
        }
        
        cerr << "Generating unlabeled flags of sizes " << loaded_size+1 << " to " << Kn << endl;
        
        for (int i = loaded_size+1; i <= Kn; i++)
        {
            generate_unlabeled_flags_of_size(i, verbose_output, dump_unlabeled_while_generating);
            cerr << "Unlabbeled flags of size "<< i << ": " << g_unlabeled_flags[i].size() << endl;
        }
        
        
        dump_unlabeled_flags(Kn);
    }

    if ( (int)g_unlabeled_flags[Kn].size() < 1 )
    {
        cerr << "No flags of size " << Kn << endl;
        cerr << "Generating flags of size " << Kn << " failed." << endl;
        exit(0);
    }

}


// This should make a cash of loaded flags and have the things a little faster
unordered_map<string, vector<flag> > g_labeled_flags_of_one_type_map;

inline string get_filename_for_labeled_flags(int flag_size, const flag &type)
{
    // Saving of the flags to a separate file and/or loading them
    return filename_prefix() + "__size"+to_string((long long)flag_size)+"_type"+ type.print("")+".txt";
}

bool labeled_flags_of_one_already_exist(int flag_size, const flag &type)
{
    string filename = get_filename_for_labeled_flags(flag_size,type);

    // if it is already in the map, it exists
    if (g_labeled_flags_of_one_type_map.find(filename) != g_labeled_flags_of_one_type_map.end())
        return true;

    vector<flag> flag_list;
    if (load_flags_from_file(filename,flag_list))
    {
        g_labeled_flags_of_one_type_map.insert(make_pair(filename, flag_list));
        return true;
    }
    
    return false;
}


void generate_next_map(int type_size, int flag_size, int id, vector<int> &current_map,  vector<bool> &used_in_map, vector<vector<int> > &mappings_to_try)
{
    // mapping found
    if (id >= flag_size)
    {
        mappings_to_try.push_back(current_map);
        return;
    }

    // finish mapping of unlabeled vertices in any order
    if (id >= type_size)
    {
        for (int i = 0; i < flag_size; i++)
        {
            if (used_in_map[i] == false)
            {
                used_in_map[i] = true;
                current_map[id] = i; 
                generate_next_map(type_size, flag_size, id+1, current_map,  used_in_map, mappings_to_try);
                used_in_map[i] = false;
                return;
            }
        }
        assert(0);
    }

    // Try all possibilities for labeled vertices
    for (int i = 0; i < flag_size; i++)
    {
        if (used_in_map[i] == false)
        {
            used_in_map[i] = true;
            current_map[id] = i; 
            generate_next_map(type_size, flag_size, id+1, current_map,  used_in_map, mappings_to_try);
            used_in_map[i] = false;
        }
    }
}

void generate_maps(int type_size, int flag_size, vector<vector<int> > &mappings_to_try)
{
    vector<int> current_map;
    vector<bool> used_in_map;
    for (int i = 0; i < flag_size; i++)
    {
        current_map.push_back(i);
        used_in_map.push_back(false);        
    }
    generate_next_map(type_size, flag_size, 0, current_map,  used_in_map, mappings_to_try);
}


// Not paralel version - very slow and not completely correct - in particular for ordered version.
void generate_labeled_flags_of_one_type(int flag_size, const flag &type, vector<flag> &flag_list)
{

    cerr << "Generating labeled flags of size " << flag_size << " of type " << type.print() << endl;

    string filename = get_filename_for_labeled_flags(flag_size,type);
                    
    // Make sure we have the unlabeled flags of right size loaded
    get_unlabeled_flags_of_size(flag_size);
    
    int to_label_size = (int)g_unlabeled_flags[flag_size].size();
    //cerr << "Type size " << type_size << " for " << g_unlabeled_flags[flag_size][i].print() << endl;
    int mapping[flag_size];
    flag F;
    
    int type_size = type.labeled_vertices_cnt();
    
    if (type_size != 0)
    {
    
        vector<vector<int> > mappings_to_try;
        generate_maps(type_size, flag_size, mappings_to_try);

        vector<flag> found_already;
        for (int i = 0; i < to_label_size; i++)
        {
            found_already.clear();
            //if (verbose_output)
            //     cerr << "Flag size "<<flag_size << " type size " << type_size << " labeling " << i << "/" << to_label_size << endl;
        
            //for (int j = 0; j < flag_size; j++) mapping[j]=j;
        
            //do {
            for (int mapID = 0; mapID < (int)mappings_to_try.size(); mapID++)
            {
                for (int j = 0; j < flag_size; j++) mapping[j]=mappings_to_try[mapID][j];

                F.as_subflag(g_unlabeled_flags[flag_size][i], mapping, flag_size, type_size); // taking a subflag
                //cerr << "Found " << F.print() << " X " << F.have_same_type(type) << " " << !find_flag_in_list(F, flag_list) << endl;
                //if (F.have_same_type(type) && find_flag_in_list(F, flag_list) == -1)
                //if (F.have_same_type(type) && find_flag_in_list_nonparalel(F, flag_list) == -1)
                if (F.have_same_type(type) && find_flag_in_list_nonparalel(F, found_already) == -1)
                {
                    found_already.push_back(F);
                    flag_list.push_back(F);
                }
                //cerr << "x";
            }
            //} while ( std::next_permutation(mapping,mapping+flag_size) );
            //cerr << endl;
            //cerr << "found_already.size()=" << found_already.size() << " flag_list.size()=" << flag_list.size() << endl;
        }
    }
    else
    {
        //cerr << "Bubu " << flag_size << " "  << g_unlabeled_flags[flag_size].size() << endl;
        //flag_list.erase();
        flag_list = g_unlabeled_flags[flag_size];
    }    
    
    if (dump_flags_to_file(filename,flag_list))
    {
        cerr << "Generated and written flags to " << filename << endl;
    } 
}

// Not paralel version - very slow and not completely correct - in particular for ordered version.
void get_labeled_flags_of_one_type(int flag_size, const flag &type, vector<flag> &flag_list)
{
    
    //cerr << "XXXXX" << endl;

    // Saving of the flags to a separate file and/or loading them
    string filename = get_filename_for_labeled_flags(flag_size,type);

    //umap.insert(make_pair("e", 2.718)); 
 
    if (g_labeled_flags_of_one_type_map.find(filename) == g_labeled_flags_of_one_type_map.end()) 
    {
        #pragma omp critical (generating_special_labeled_flags)
        {   
            if (g_labeled_flags_of_one_type_map.find(filename) == g_labeled_flags_of_one_type_map.end())
            {
                //cerr << "About to start loading " << filename << endl;
                if (load_flags_from_file(filename,flag_list))
                {
                    cerr << "Loaded " <<  flag_list.size() << " flags from " << filename << endl;
                }
                else
                {
                    cerr << "Generating flags for file " << filename << endl;
                    generate_labeled_flags_of_one_type(flag_size, type, flag_list);   
                }
                    
                g_labeled_flags_of_one_type_map.insert(make_pair(filename, flag_list));
            }
        }   
    }
    else
    {
        //cerr << "Type already exists" << endl;
    }


    assert(g_labeled_flags_of_one_type_map.find(filename) != g_labeled_flags_of_one_type_map.end());

    flag_list = g_labeled_flags_of_one_type_map[filename];
}



double unlabel_F1(const flag &F1, int remaining_labeled)
{
    
    flag H = F1;
    H.m_Theta = remaining_labeled;
    
    return P_F1_IN_H(F1, H);
}

vector<flag_and_coefficient> F1_times_F2(const flag &F1, const flag &F2, const flag &type)
{
    // do something with ID's
    assert(F1.have_same_type(F2));
    
    int labeled_vertices_cnt = F1.labeled_vertices_cnt();
    int result_size = F1.m_vertices + F2.m_vertices - labeled_vertices_cnt;

    //cout << labeled_vertices_cnt << " " << result_size << endl;
    
    vector<flag> big_list;

    get_labeled_flags_of_one_type(result_size, type, big_list);
    
    //cout << big_list.size() << endl;
    
    vector<flag_and_coefficient> combination;
    
    for (int i = 0; i < (int)big_list.size(); i++)
    {
        double d = P_F1_F2_IN_labeled_H(F1,F2,big_list[i]);
        if (d != 0)
        {
            flag_and_coefficient fc;
            fc.g = big_list[i];
            fc.coefficient = d;
            combination.push_back(fc);
        }
    }

    return combination;
}





void fc_add_FC_to_first(vector<flag_and_coefficient> &FC, vector<flag_and_coefficient> &FC_add, double coeff = 1)
{
    for (int j = 0; j < (int)FC_add.size(); j++)
    {
        int found_ID = -1;
        for (int k = 0; k < (int)FC.size(); k++)
        {
            if (FC_add[j].g.is_isomorphic_to(FC[k].g))
            {
                found_ID = k;
                break;
            }
        }
        if (found_ID == -1)
        {
            flag_and_coefficient fc = FC_add[j];
            fc.coefficient *= coeff;
            FC.push_back(fc);
            
        }
        else
        {
            FC[found_ID].coefficient += coeff*FC_add[j].coefficient;
        }
    }
}

void multiply_FC_by_C(vector<flag_and_coefficient> &FC, double c)
{
    for (int i = 0; i < (int)FC.size(); i++)
    {
        FC[i].coefficient *= c;
    }
}

void fc_F1_times_F2(vector<flag_and_coefficient> &F1, vector<flag_and_coefficient> &F2, vector<flag_and_coefficient> &F1F2sum)
{
    F1F2sum.clear();

    //cerr << "Sizes for product " << F1.size() << " "  << F2.size() << endl;

    for (int i = 0; i < (int)F1.size(); i++)
    {
        for (int j = 0; j < (int)F2.size(); j++)
        {
            //cerr << "Multiplying " << i << " " << j << endl;
            flag type;
            F1[i].g.get_type_subflag(type);
            vector<flag_and_coefficient> F1F2 = F1_times_F2(F1[i].g, F2[j].g, type);
            
            //cerr << F1F2.size() << endl;
            
            multiply_FC_by_C(F1F2, F1[i].coefficient*F2[j].coefficient);
            
            fc_add_FC_to_first(F1F2sum,F1F2);
        }
    }
}


/*
Expression := Term { ("+" | "-") Term }
Term       := Factor { ( "*" | "/" ) Factor }
Factor     := Power [^2]
Power      := CFlag | "{" Expression "}"
CFlag      := [-]CExpression  Flag
Flag       := Bunch of Digits
Digit      := "0" | "1" | "2" | "3" | "4" | "5" | "6" | "7" | "8" | "9"

CExpression := CTerm { ("+" | "-") CTerm }
CTerm       := CFactor { ( "*" | "/" ) CFactor }
CFactor     := CPower [^ CExpression]
CPower      := CNumber 
CNumber      := [-]Bunch of Digits inclduing .

https://codereview.stackexchange.com/questions/54273/simple-c-calculator-which-follows-bomdas-rules
*/


bool fc_cexpression(istream *istr, double &result, ostream *ostr);
void fc_expression(istream *istr, vector<flag_and_coefficient> &result, ostream *ostr);

char fc_token(istream *istr) 
{ 
    char ch;
    (*istr) >> ch;
    return ch;
}

bool fc_cfactor(istream *istr, double &result, ostream *ostr) 
{ 
    result = 0;
    char ch = fc_token(istr);
    if (ch == '{') 
    {
        if (ostr != NULL) (*ostr) << " \\left( ";
        fc_cexpression(istr, result, ostr);
        ch = fc_token(istr);
        if (ch != '}') 
        {
            std::string error = std::string("Expected '}', got: ") + ch;
            throw std::runtime_error(error.c_str());
        }
        if (ostr != NULL) (*ostr) << " \\right) ";
    }
    else if (isdigit(ch) || ch=='-' || ch == '.') {
        istr->unget();

        (*istr) >> result;

        if (ostr != NULL) 
        {
            (*ostr) << result << " ";
        }
    }
    else 
    {
        std::string error = std::string("Unexpected character fc_cfactor got: ") + ch;
        throw std::runtime_error(error);
    }
    return true;
}

bool fc_cpower(istream *istr, double &result, ostream *ostr) { 
    char ch;
    fc_cfactor(istr, result, ostr);
    ch = fc_token(istr);
    if (ch == '^') 
    {
        double power;
        if (ostr != NULL) (*ostr) << "^{";
        fc_cexpression(istr, power, ostr);
        if (ostr != NULL) (*ostr) << "}";
        result = pow(result,power);
    }
    else istr->unget();
    return true;
}

bool fc_cterm(istream *istr, double &result, ostream *ostr) { 
    char ch;
    fc_cpower(istr, result, ostr);

    while(true)
    {
        ch = fc_token(istr);
        if (ch == '*' || ch == '/')
        { 
            if (ch == '*') 
            {
                if (ostr != NULL) (*ostr) << " \\times ";
                double b;
                fc_cterm(istr, b, ostr);

            result *= b;
            }
            if (ch == '/') 
            {
                if (ostr != NULL) (*ostr) << " / ";
                double b;
                fc_cterm(istr, b, ostr);

                if (b == 0)
                {
                    throw std::runtime_error("Division by 0 in flag coefficient");
                }

            result /= b;
            }
        }
        else 
        {
            istr->unget();
            break;
        }
    }
    return true;
}

// One has to be careful with expresion because of -
bool fc_cexpression(istream *istr, double &result, ostream *ostr) 
{
    fc_cterm(istr, result, ostr);

    while(true)
    {
        char ch = fc_token(istr);
        if (ch == '-' || ch=='+') 
        {
            if (ostr != NULL) (*ostr) << ch;
            double b;
            fc_cterm(istr, b, ostr);
            if (ch == '+')
            {
                result += b;
            }
            else
            {
                result -= b;
            }
        }
        else 
        {
            istr->unget(); 
            break;
        }
    }
    return true;
}


void fc_factor(istream *istr, vector<flag_and_coefficient> &result, ostream *ostr) 
{ 
    result.clear();
    char ch = fc_token(istr);
    if (ch == '(') 
    {
        if (ostr != NULL) (*ostr) << " \\left( ";
        fc_expression(istr, result, ostr);
        ch = fc_token(istr);
        if (ch != ')') 
        {
            std::string error = std::string("Expected ')', got: ") + ch;
            throw std::runtime_error(error.c_str());
        }
        if (ostr != NULL) (*ostr) << " \\right) ";
    }
    else if (isdigit(ch) || ch=='-' || ch=='.' || ch=='{')  
    {
        istr->unget();

        flag_and_coefficient fc;
//        fc_cexpression(istr, fc.coefficient, ostr);
        fc_cexpression(istr, fc.coefficient, NULL);
        //(*istr) >> fc.coefficient;
        fc.g.load_from_stream((*istr),-1,-1);

        if (ostr != NULL) 
        {
            //(*ostr) << " [" << fc.coefficient << "]" << " \\vc{" << fc.g.print_latex(false, 0) << "} ";
            (*ostr) << " " << fc.coefficient  << " \\vc{" << fc.g.print_latex(false, 0) << "} ";
        }

        result.push_back(fc);
    }
    else 
    {
        std::string error = std::string("Unexpected character in fc_factor got: ") + ch;
        throw std::runtime_error(error);
    }
}

void fc_power(istream *istr, vector<flag_and_coefficient> &result, ostream *ostr) { 
    char ch;
    fc_factor(istr, result, ostr);
    ch = fc_token(istr);
    if (ch == '^') 
    {
        ch = fc_token(istr);
        if (isdigit(ch))
        {
            int power=(int)(ch-'0');
            if (power < 1 || power > 10)
            {
                throw std::runtime_error("Unexpected power ");
            }
            if (ostr != NULL) (*ostr) << " ^{"<<power<<"} ";
            vector<flag_and_coefficient> base = result;
            for(;power > 1;power--)
            {
                vector<flag_and_coefficient> tmp = result;
                fc_F1_times_F2(tmp,base,result);
            }
        }
        else
        {
            throw std::runtime_error("^ must be followed by a digit");
        }
    }
    else istr->unget();
}

void fc_term(istream *istr, vector<flag_and_coefficient> &result, ostream *ostr) { 
    char ch;
    fc_power(istr, result, ostr);
    ch = fc_token(istr);
    //if (ch == '*' || ch == '/') 
    if (ch == '*') 
    {
        if (ostr != NULL) (*ostr) << " \\times ";
        vector<flag_and_coefficient> b;
        fc_term(istr, b, ostr);

        vector<flag_and_coefficient> tmp=result;

        fc_F1_times_F2(tmp,b,result);
    }
    else istr->unget();
}

// One has to be careful with expresion because of -
void fc_expression(istream *istr, vector<flag_and_coefficient> &result, ostream *ostr) 
{
    fc_term(istr, result, ostr);

    while(true)
    {
        char ch = fc_token(istr);
        if (ch == '-' || ch=='+') 
        {
            if (ostr != NULL) (*ostr) << ch;
            vector<flag_and_coefficient> b;
            fc_term(istr, b, ostr);
            if (ch == '+')
            {
                fc_add_FC_to_first(result, b);
            }
            else
            {
                multiply_FC_by_C(b, -1);
                fc_add_FC_to_first(result, b);
            }
        }
        else 
        {
            istr->unget(); 
            break;
        }
    }
}

bool flag_calculator_simple(vector<flag_and_coefficient> &result, string& inputname, ostream *ostr = NULL, int verbose_output=0)
{
   OPEN_FILE_SMARTLY_RETURN_FALSE_ON_FAIL(istr, inputname);

    //ostr = &cout;

    if (ostr != NULL)
    {
        (*ostr) << " $ ";
    }

    fc_expression(istr, result, ostr);

    if (ostr != NULL)
    {
        (*ostr) << " $ ";
    }

    return true;
}


    
// f is an unlabeled flag
void count_flag_products(stringstream &ss, int matrixID, const flag &f)
{    
    int subset_split[V];
    int mapping_F1[V];
    int mapping_F2[V];
    flag F1, F2; 
 
    // Try to find decomposition of vertices of f into 3 sets.
    // first is of size typesize and the other two are of sizes privatesize
    // Then typesize+privatesize creates one flag.
    
    // We go through all possible type sizes
    for (int typesize = f.m_vertices-2; typesize >= 0; typesize -= 2)
    {
        vector<pair<string,int> > found_pairs;

        int privatesize = (f.m_vertices-typesize)/2; 
        int totalsize = typesize+privatesize;  // size of F1 and F2
        int label_f1 = typesize;
        int label_f2 = typesize+1;

        for (int i = 0; i < typesize; i++) subset_split[i] = i;
        assert(typesize+2*privatesize <= V);
        for (int i = 0; i < privatesize; i++)
        {   
            subset_split[typesize+i] = label_f1;
            subset_split[typesize+i+privatesize] = label_f2; // I don't know why the warning is here...
        }
        
        do {
            //std::cout << subset_split[0] << ' ' << subset_split[1] << ' ' << subset_split[2] << ' ' << subset_split[3] << '\n';
        
            // Construction of inverz mapping
            int inF1 = typesize;
            int inF2 = typesize;            
            for (int i = 0; i < f.m_vertices; i++)
            {
                if (subset_split[i] < typesize) 
                {
                    mapping_F1[subset_split[i]] = i;
                    mapping_F2[subset_split[i]] = i;
                }
                if (subset_split[i] == label_f1) mapping_F1[inF1++] = i; 
                if (subset_split[i] == label_f2) mapping_F2[inF2++] = i; 
            }

            F1.as_subflag(f,mapping_F1, totalsize, typesize);
        
            int typeID = get_flag_type_in_list(F1, g_flags);
            if (typeID == -2) continue;
#ifdef G_NOT_ALL_FLAGS_USED
            if (typeID == -1) continue;
#endif            
            assert (typeID != -1);
            
            F2.as_subflag(f,mapping_F2, totalsize, typesize);

            int idF1 = find_flag_in_list(F1,g_flags[typeID]);
            int idF2 = find_flag_in_list(F2,g_flags[typeID]);
            if (idF1 > idF2 ) continue; // count only once when idF1 <= idF2         
            if (idF1 == -1)
            {
#ifdef G_NOT_ALL_FLAGS_USED
                continue;
#endif                
                cerr << "Missing " << F1.print() << endl;
            }
            assert(idF1 >= 0 && idF2 >= 0); 
            
            stringstream what_found_ss;
            what_found_ss << typeID+1 << " " << idF1+1 << " " << idF2+1;
            string what_found = what_found_ss.str();
            
            bool known = false;
            for(int i = 0; i < (int)found_pairs.size(); i++)
            {
                if (found_pairs[i].first == what_found) 
                {
                    found_pairs[i].second++;
                    known = true;
                    break;
                }
            }
            if (known == false)
            {
                found_pairs.push_back(make_pair(what_found,1));
            }
            
        } while ( std::next_permutation(subset_split,subset_split+f.m_vertices) );
        
        for(int i = 0; i < (int)found_pairs.size(); i++)
        {
            ss << matrixID << " " << found_pairs[i].first << " " <<  found_pairs[i].second << endl; 
        }
    }
}



double compute_linear_combination_in_g(const flag &g, const vector<flag_and_coefficient> &lc)
{
    double total_density = 0;    
    
    for (int i = 0; i < (int)lc.size(); i++)
    {
        total_density += lc[i].coefficient * P_F1_IN_H(lc[i].g, g);
    }
    
    return total_density;
}



double compute_objective_function_for_SDP(const flag &g)
{

    return compute_linear_combination_in_g(g, g_objective_combination);
    
	double total_density = 0;	

	for (int i = 0; i < (int)g_objective_combination.size(); i++)
	{
        total_density += g_objective_combination[i].coefficient * P_F1_IN_H(g_objective_combination[i].g, g);
        
        // HACK
        //cerr << P_F1_IN_H(g_objective_combination[i].g, g) << " " << g.m_rotation_system_noncrosssings << endl;
    }
    
	return total_density;
}


int print_CSDP_simple_linear_constraints(ostream &ostr, int Kn, int i, int matrixID, int blockID, bool print_blocks, bool print_products, int verbose_output)
{
    
    if (g_use_simple_linear_constraints_list.size() == 0)
    {
        for (int c = 0; c < (int)g_linear_constraints.size(); c++)
        {
            g_use_simple_linear_constraints_list.push_back(c);
        }
    }
    
    int constraint_id = 0;
//    for (int j = 0; j < (int)g_linear_constraints.size(); j++)
    for (int c = 0; c < (int)g_use_simple_linear_constraints_list.size(); c++)
    {
        int j = g_use_simple_linear_constraints_list[c];
        assert(g_linear_constraints[j].m_constant == 0);
        
        constraint_id++;
        
        if (print_products)
        {
            double d = 0;
            for (int k = 0; k < (int)g_linear_constraints[j].m_entries.size(); k++)
            {
                assert (g_linear_constraints[j].m_entries[k].coefficient != 0);
                d += g_linear_constraints[j].m_entries[k].coefficient * P_F1_IN_H(g_linear_constraints[j].m_entries[k].g, g_unlabeled_flags[Kn][i]);
            }
            if (d != 0)
            {
                ostr.precision(G_PRECISION);
                ostr <<  matrixID << " " << blockID << " " << constraint_id << " " << constraint_id << " " << smart_round(d) << endl;
            }
        }
    } 
    
    if (print_blocks)
    {
        ostr << -constraint_id << " ";
        if (verbose_output)
        {
            cerr << "Simple linear are used for " << constraint_id << "/" << g_linear_constraints.size() << " constraints" << endl;
        }
    }
    
    return constraint_id;
}


// We assume print_products may run in parallel
int print_CSDP_product_linear_constraints(ostream &ostr, int Kn, int i, int matrixID, int blockID, bool print_blocks, bool print_products, int verbose_output)
{
    // g_flag_product_linear_constraints[ constraint ] is a vector of flags that can be used to multiply  [constraint]
    static vector<vector<flag> >  g_flag_product_linear_constraints;
    
    if (g_flag_product_linear_constraints.size() != g_linear_constraints.size())
        g_flag_product_linear_constraints.resize(g_linear_constraints.size());
    
    
    int constraints_possible = 0;
    
    int constraint_id = 0;
    for (int j = 0; j < (int)g_linear_constraints.size(); j++)
    {
         // The linear constraints need the same types
        //if (!g_linear_constraints[j].m_same_types) continue;
        
                
        // Size of the other flags in multiplication
        int type_size = g_linear_constraints[j].m_labeled_vertices_in_type_cnt;         
        int flag_size = (Kn - g_linear_constraints[j].m_entries_max_size) + type_size;

        if (flag_size < 2)   // multiplying by just 1 vertex does not do much (unless the vertex has a color)
        {
            if (print_blocks)
                cerr << "Size of big graphs is too small to use products linear constraint " << j+1 << endl;
            continue;
        }
        
        if (print_products && g_flag_product_linear_constraints[j].size() == 0)
        {
            continue;
        }
                
        // getting labeled flags of the right size
        if (g_flag_product_linear_constraints[j].size() == 0)
        {
            flag type;
            type = g_linear_constraints[j].m_type;
            get_labeled_flags_of_one_type(flag_size, type, g_flag_product_linear_constraints[j]);
            
            if (g_flag_product_linear_constraints[j].size() == 0)
            {
                cerr << "Linear constraint " << j+1 << " cannot be used since generating flags of type " <<  type.print() << " and size " <<  flag_size  << " gave zero flags" << endl;
                continue;
            }
            
            if (verbose_output)
                cerr << "Constraint " << j << "/" << g_linear_constraints.size() << " can be multiplied by " <<  g_flag_product_linear_constraints[j].size()  << " other flags " << endl;
            //for (int i = 0; i < g_flag_product_linear_constraints[j].size(); i++)
            //{
            //    cerr << g_flag_product_linear_constraints[j][i].print() << endl;
            //}
            
        }
        
        constraints_possible++;

        for (int z = 0; z < (int)g_flag_product_linear_constraints[j].size(); z++)
        {
            constraint_id++;
            
            if (print_products)
            {
                // constant
                double d = 0; //g_linear_constraints[j].m_constant * P_F1_IN_H(g_flag_product_linear_constraints[j][z], g_unlabeled_flags[Kn][i]);
                
                // products with other flags in the linear constraint
                for (int k = 0; k < (int)g_linear_constraints[j].m_entries.size(); k++)
                {
                    assert (g_linear_constraints[j].m_entries[k].coefficient != 0);
                    double PF = P_F1_F2_IN_labeled_H(g_linear_constraints[j].m_entries[k].g, g_flag_product_linear_constraints[j][z], g_unlabeled_flags[Kn][i]);
                    
                    d += PF * g_linear_constraints[j].m_entries[k].coefficient;
                }
                if (d != 0)
                {
                    ostr.precision(G_PRECISION);
                    ostr << matrixID << " " << blockID << " " << constraint_id << " " << constraint_id << " " << smart_round(d) << endl;
                    
                }
            }
        }
    } 
    
    
    if (print_blocks && constraint_id != 0)
    {
        ostr << -constraint_id << " ";
        //if (verbose_output)
        {
            cerr << "Product linear constraints are working for " << constraints_possible << "/" << g_linear_constraints.size() << " constraints" << endl;
        }
    }   
    
    return constraint_id;
}

// This is not yet written
int print_CSDP_product_linear_constraints_extra_lines_in_csdp_nontrivial_sum(ostream &ostr, int Kn, int matrixID, int blockID, bool print_constraints, bool  print_constants, int verbose_output)
{
    return 0;
}


void generate_all_types_containing_one_subtype(int big_type_size, const flag &type, vector<flag> &flag_list)
{
    flag g;
    g.set_vertices_and_Theta(big_type_size,big_type_size);
    g.copy_from(type);
    
    extensions_of_g(g,flag_list);
}




/*
 This is pocessing constraints like  F1+F2 >= 0 
 It tries
 (F1+F2)  [(F^T M F)]
 Let the type of the constraint be t. We take all types t_i that contain t as a subtype.
 The we compute  F_i^T M_i F_i  and do partial unlabeling down to t and finally multiply with (F1+F2).
So we have a set of constraints of form
 (F1+F2) [(F_i^T M_i F_i)]_t
 */
int print_CSDP_square_linear_constraints_v2(ostream &ostr, int Kn, int i, int matrixID, int blockID, bool print_blocks, bool print_products, int verbose_output)
{  
    static vector<vector< vector<vector<flag> > > > g_flag_square_linear_constraints_v2; // flag used only for linear constraints
    
    // For a constraint i, it remembers all 
    static vector< vector<vector<flag> > >  g_flag_square_linear_constraints_v2_typelist; // flag used only for linear constraints

    
    if (g_flag_square_linear_constraints_v2.size() != g_linear_constraints.size())
    {
        g_flag_square_linear_constraints_v2.resize(g_linear_constraints.size());
        g_flag_square_linear_constraints_v2_typelist.resize(g_linear_constraints.size());
    }
    
    
    int constraint_id = 0;
    int constraints_applicable = 0;
    
    for (int j = 0; j < (int)g_linear_constraints.size(); j++)
    {
        // TODO: Shuld be done sooner
        //g_linear_constraints[j].check_if_same_types();
        
        // The linear constraints need the same types
        // if not, just skip this shit
        assert(g_linear_constraints[0].m_checked);
        //if (!g_linear_constraints[j].m_same_types) continue;
         
        // Size of the other flags in multiplication
        int type_size = g_linear_constraints[j].m_labeled_vertices_in_type_cnt;
        
        flag type;
        g_linear_constraints[j].m_entries[0].g.get_type_subflag(type);
        
        int free_vertices = Kn - g_linear_constraints[j].m_entries_max_size; // vertices not used by the constraint
        if (free_vertices < 2)
        {
                if (print_blocks)
                    cerr << "Size of unlabeled flags is too small to use squares in constraint " << j+1 << endl;
                continue;
        }
        
        constraints_applicable++;

        // generating additional flags
        // int current_type_id=0;
        //int max_inner_type_size = type_size+free_vertices-2;
        int max_inner_type_size = type_size;

        
        if (g_flag_square_linear_constraints_v2_typelist[j].size() == 0)
        {
            if (print_products) continue;
            g_flag_square_linear_constraints_v2_typelist[j].resize(max_inner_type_size+1);
            g_flag_square_linear_constraints_v2[j].resize(max_inner_type_size+1);
        }
                
        for (int inner_type_size = type_size; inner_type_size <= max_inner_type_size; inner_type_size++)
        {        
            int flag_size = (free_vertices-(inner_type_size-type_size))/2 + inner_type_size;
            
            //if (print_blocks)
                //cerr << "Inner type size " << inner_type_size << " flag size " << flag_size << endl;
            assert(flag_size > inner_type_size);
            
            if (flag_size == 0)
            {
                if (print_blocks)
                    cerr << "Size of big graphs is too small to use products in constraint " << j+1 << endl;
                continue;
            }
        
           // if (print_products && g_flag_square_linear_constraints_v2[j].size() == 0)
           // {
           //     continue;
           // }
        
            // generate all types of the right size
            if (g_flag_square_linear_constraints_v2_typelist[j][inner_type_size].size() == 0)
            {
                if (print_products) continue;
                generate_all_types_containing_one_subtype(inner_type_size, type, g_flag_square_linear_constraints_v2_typelist[j][inner_type_size]);
                g_flag_square_linear_constraints_v2[j][inner_type_size].resize(g_flag_square_linear_constraints_v2_typelist[j][inner_type_size].size());
                
                //if (print_blocks)
                //    cerr << "Generated " << g_flag_square_linear_constraints_v2_typelist[j][inner_type_size].size() << " inner type(s)" << endl;
                
            }
            
            for (int it = 0; it < (int)g_flag_square_linear_constraints_v2_typelist[j][inner_type_size].size(); it++)
            {
                
                flag inner_type = g_flag_square_linear_constraints_v2_typelist[j][inner_type_size][it];
                
                
                if ( g_flag_square_linear_constraints_v2[j][inner_type_size][it].size() ==  0)
                {
                    if (print_products) break;
                    get_labeled_flags_of_one_type(flag_size, inner_type, g_flag_square_linear_constraints_v2[j][inner_type_size][it]);
                    if (g_flag_square_linear_constraints_v2[j].size() == 0)
                    {
                        cerr << "Linear constraint " << j+1 << " cannot be used since generating of other flags of type " <<  type.print() << " and size " <<  flag_size  << " failed" << endl;
                        continue;
                    }
                    //if (print_blocks)
                    //    cerr << "Generated " << g_flag_square_linear_constraints_v2[j][inner_type_size][it].size() << " flags of type " << inner_type.print() << endl;
                    
                    
                }
                
                //if (verbose_output)
                //    cerr << "Constraints " << j << " is using "  << g_flag_square_linear_constraints_v2[j][inner_type_size][it].size() << " flags " << endl;
                //for (int i = 0; i < g_flag_square_linear_constraints_v2[j].size(); i++)
                //{
                //    cerr << g_flag_square_linear_constraints_v2[j][i].print() << endl;
                //}
                
                
                
                if (print_blocks)
                {
                    ostr << g_flag_square_linear_constraints_v2[j][inner_type_size][it].size() << " ";
                }
                
                if (print_products)
                {
                    //cerr <<  g_unlabeled_flags[Kn][i].print() << endl;
                    
                    for (unsigned int x = 0; x < g_flag_square_linear_constraints_v2[j][inner_type_size][it].size(); x++)
                        for (unsigned int y = x; y < g_flag_square_linear_constraints_v2[j][inner_type_size][it].size(); y++)
                        {
                            flag F1 = g_flag_square_linear_constraints_v2[j][inner_type_size][it][x];
                            flag F2 = g_flag_square_linear_constraints_v2[j][inner_type_size][it][y];
                            
                            // product
                            vector<flag_and_coefficient> F1F2 = F1_times_F2(F1, F2, type);
                            
                            // constant
                            double d = 0;
                            for (int z = 0; z < (int)F1F2.size();z++)
                            {
                                d += g_linear_constraints[j].m_constant * F1F2[z].coefficient * P_F1_IN_H(F1F2[z].g, g_unlabeled_flags[Kn][i]);
                                                       
                                for (int k = 0; k < (int)g_linear_constraints[j].m_entries.size(); k++)
                                {
                                
                                    //cerr << "Computing poduct " << F1.print() << "X" << F2.print() << " =" << F1F2[z].g.print() << " X " << g_linear_constraints[j].m_entries[k].g.print() << " in " << g_unlabeled_flags[Kn][i].print() << endl;
                                    double PF1F2g = P_F1_F2_IN_labeled_H(F1F2[z].g,g_linear_constraints[j].m_entries[k].g, g_unlabeled_flags[Kn][i]);
                                
                                    d += PF1F2g * F1F2[z].coefficient * g_linear_constraints[j].m_entries[k].coefficient; 
                                
                                }
                            }
                            if (d != 0)
                            {
                                ostr.precision(G_PRECISION);
                                 ostr << matrixID <<  " " << blockID+constraint_id << " " << x+1 << " " << y+1 << " " <<  smart_round(d) << endl;
                            }
                        }
                }                
                constraint_id++;
            }
            
        }
    }    
    
    if (print_blocks)
    {
        cerr << "Square linear constraints are working for " << constraints_applicable << "/" << g_linear_constraints.size() << " constraints" << endl;
        cerr << "Number of used blocks is " << constraint_id << endl;
        
    }
    
    //cerr << "Returning " << constraint_id << endl;
    
    return constraint_id;
}



bool has_slack(int Kn, int i)
{
    if (find_flag_in_list(g_unlabeled_flags[Kn][i],g_no_slack_flags) != -1)
    {
        cerr << "Flag " << i << " has no slack" << endl;
    }
    return (find_flag_in_list(g_unlabeled_flags[Kn][i],g_no_slack_flags) == -1);
}

// When calculating parts of the SDP in parallel, we want to have a nice ordered output
// and we do not want to wait for the things to be ordered. So here thing that orders the
// output. By increasing the id which is being written.
//
class parallel_output_linearizer
{
public:

    parallel_output_linearizer(ostream &ostr, int first_ID):
    m_ostr(ostr)
    {
        m_next_ID = first_ID;
    }

    void print_string(int id, const string &str)
    {
        if (id == m_next_ID)
        {
            m_ostr << str;
            m_next_ID++;
            try_writing_next();
        }
        else
        {
            m_to_write.insert(make_pair(id, str));
        }
    }

    void try_writing_next()
    {
        //std::unordered_map<std::string,double>::const_iterator got = mymap.find (input);
        auto next_found = m_to_write.find(m_next_ID);
        while (next_found != m_to_write.end())
        {
            m_ostr << next_found->second;
            m_next_ID++;
            m_to_write.erase(next_found);
            next_found = m_to_write.find(m_next_ID);
        }
    }

    bool is_empty()
    {
        return m_to_write.empty();
    }

    ostream &m_ostr;
    int m_next_ID;
    unordered_map<int, string> m_to_write;
};



int print_CSDP_constraints_header(int Kn, ostream &ostr = cout, int verbose_output = 1)
{   
    int blocks = 0;
    
    // next we have linear constraints blocks

    if (g_linear_constraints.size() != 0)
    {
        if (g_use_simple_linear_constraints)
        {
            // Simple linear constraints
            int simple_constraints = print_CSDP_simple_linear_constraints(ostr, Kn, 0, 0, 0, true, false, verbose_output);
            if (simple_constraints != 0)
            {
                //block_ID_simple_linear_constraints = blocks+1;
                blocks++;
            }
        }
        
        if (g_use_product_linear_constraints)
        {
            int product_constraints = print_CSDP_product_linear_constraints(ostr, Kn, 0, 0, 0, true, false, verbose_output);
            if (product_constraints != 0)
            {
                //block_ID_product_linear_constraints = blocks+1;
                blocks++;
            }
        }
        
        
        if (g_use_square_linear_constraints)
        {
            int square_constraints_v2 = print_CSDP_square_linear_constraints_v2(ostr, Kn, 0, 0, 0, true, false, verbose_output);
            if (square_constraints_v2 != 0)
            {
                //block_ID_square_linear_constraints_v2 = blocks+1;
                blocks += square_constraints_v2;
            }
        }
    }  
    
    return blocks;
}

        
int print_CSDP_constraints_blocks(ostream &ostr, int Kn, int i, int matrixID, int current_csdp_block, int verbose_output = 1)
{        
    if (g_use_simple_linear_constraints)
    {
        if (print_CSDP_simple_linear_constraints(ostr, Kn, i, matrixID, current_csdp_block, false, true, verbose_output) != 0)
        {
            current_csdp_block++;
        }
    }
    
    if (g_use_product_linear_constraints)
    {
        if (print_CSDP_product_linear_constraints(ostr, Kn, i, matrixID, current_csdp_block, false, true, verbose_output) != 0)
        {
            current_csdp_block++;
        }
    }


    if (g_use_square_linear_constraints)
    {
        current_csdp_block += print_CSDP_square_linear_constraints_v2(ostr, Kn, i, matrixID, current_csdp_block, false, true, verbose_output);
    }  
    
    return current_csdp_block;
}




int print_CSDP_additional_blocks_header(int Kn, ostream &ostr, int verbose_output = 1)
{   
    int blocks = 0;
    
    for (int i = 0; i < (int)g_additional_csdp_blocks.size(); i++)
    {
        ifstream infile;
        infile.open (g_additional_csdp_blocks[i].c_str(), ifstream::in);
        if (!infile.good())
        {
            cerr << "Failed opening file " << g_additional_csdp_blocks[i] << endl;
            assert(0);
            return 0;
        }
        int blockKn = 0;
        infile >> blockKn;

        if (blockKn != Kn)
        {
            cerr << "Loading additional CSDP blocks from " << g_additional_csdp_blocks[i] << " failed." << endl;
            cerr << "Blocks were computed n=" << blockKn << " while this program is running with n="  << Kn << endl;
            assert(0);
            return 0;
        }
        
        int blocksHere = 0;
        infile >> blocksHere;
        assert(blocksHere != 0);
        
        blocks += blocksHere;
        
        for (int j = 0; j < blocksHere; j++)
        {
            int block_size=0;
            infile >> block_size;
            assert(block_size != 0);
            ostr << " " << block_size;
        }
                
        infile.close();
    }
    
    if (verbose_output)
    {
        cerr << "Added " << blocks << " additional blocks." << endl;
    }
    return blocks;
}
        
int print_CSDP_additional_blocks(int Kn, ostream &ostr, int block_offset, int verbose_output = 1)
{   
    int blocks = 0;
    
    for (int i = 0; i < (int)g_additional_csdp_blocks.size(); i++)
    {
        ifstream infile;
        infile.open (g_additional_csdp_blocks[i].c_str(), ifstream::in);
        if (!infile.good())
        {
            cerr << "Failed opening file " << g_additional_csdp_blocks[i] << endl;
            assert(0);
            return 0;
        }
        int blockKn = 0;
        infile >> blockKn;

        if (blockKn != Kn)
        {
            cerr << "Loading additional CSDP blocks from " << g_additional_csdp_blocks[i] << " failed." << endl;
            cerr << "Blocks were computed n=" << blockKn << " while this program is running with n="  << Kn << endl;
            assert(0);
            return 0;
        }
        
        int blocksHere = 0;
        infile >> blocksHere;
        assert(blocksHere != 0);
        
        blocks += blocksHere;
        
        cerr << "Copying " << blocksHere << " dat-s blocks from " << g_additional_csdp_blocks[i] << endl;

        for (int j = 0; j < blocksHere; j++)
        {
            int block_size=0;
            infile >> block_size;
            assert(block_size != 0);
        }

        // HERE we read/write lines as long as possible...
        int matrix_id;
        int block_id;
        int x;
        int y;
        string value;

        int lines_copied = 0;
        
        while(infile.good())
        {    
            infile >> matrix_id >> block_id >> x >> y >> value;
            if (!infile.good())
            {
                break;
            }
            assert(value.size() != 0);
            
            ostr << matrix_id << " " << block_id+block_offset << " " << x << " " << " " << y << " " << value << endl;  
            lines_copied++;
        }
        
        infile.close();
        if (verbose_output)
        {
            cerr << lines_copied << " lines copied from " << g_additional_csdp_blocks[i] << endl;
        } 
        block_offset += blocksHere;
    }
    
    if (verbose_output)
    {
        cerr << "Added " << blocks << " additional blocks." << endl;
    }
    
    return block_offset;
}

        
        


// According to http://plato.asu.edu/ftp/sdpa_format.txt
//
// Solved program
//   max C.X
//   subject to A_1.X = b_1
//              A_2.X = b_2
//              A_mdim.X = b_mdim
//   where X is positive semidefinite matrix
//
// In our UPPER BOUND application
//
//      min    t
//      s.t.  D_i + [A.X]_i <= t    
//
// We add a slack variable s_i to make the <= to = and move D_i (density - constant to the other side)
// and change to maximization:
//
//      max   -t
//      s.t.  [A.X]_i + s_i - t = -D_i
//
//   C.X is just entry 0 2 1 1 -1.0 which means that in the second block,
//   the first variable is with -1. Se we are maximizing -t which corresponds to minimizing t
//
// Now LOWER BOUND application
//
//      max    t
//      s.t.  D_i - [A.X]_i  >= t
//
// We add a slack variable s_i to make the >= to = and move D_i (density - constant to the other side)
// and change to maximization:
//
//      max   t
//      s.t.  [A.X]_i + s_i + t = D_i
//
// Finally, Additional constraints applications:
// Say we want to add  a*G >= c. Then we could add
//  k(a*G(sampled in H_i) -c) >= 0
//
// For lower bound we get (max t)
//     [A.X]_i + k*(a*G_i-c) + s_i + t = D_i
//
// For upper bound we get (min t)
//     [A.X]_i + k*(a*G_i-c) + s_i -t = D_i
//
//
void print_CSDP_specific_part(int Kn, bool upper_bound = true, ostream &ostr = cout, int verbose_output = 1)
{
    int mdim = 0;
    mdim = (int)g_unlabeled_flags[Kn].size();
    
    // hack
    int extra_onstraints  = 0;
    //extra_onstraints += print_CSDP_product_linear_constraints_extra_lines_in_csdp(ostr, Kn, 0, false, verbose_output)

    ostr << mdim + extra_onstraints << endl; //<<" =mdim" << endl; // number of constraint matrices
    
    stringstream ssblocks; // String stream for all blocks in the matrix
    int blocks = 0; // number of blocks in ssblocks
 
    // First blocks are flags products:
    for (int f = 0; f < (int)g_flags.size(); f++) ssblocks <<  g_flags[f].size() << " ";
    blocks = (int)g_flags.size();
    
    // next we have one diagonal block of slack variables
    ssblocks << -1-mdim << " ";  // get_number_of_slacks() but we don't go for it here...
    
    blocks++;

    int first_constraint_block = blocks+1;    
 
    blocks += print_CSDP_constraints_header(Kn, ssblocks, verbose_output);
    
    int additional_blocks_offset = blocks+1;
    blocks += print_CSDP_additional_blocks_header(Kn, ssblocks, verbose_output);
    
    ostr << blocks << endl;
    ostr << ssblocks.str() << endl;
 
    for (unsigned int i = 0; i < g_unlabeled_flags[Kn].size(); i++)
    {
        ostr.precision(G_PRECISION);
        if (upper_bound) ostr << -1*smart_round(compute_objective_function_for_SDP(g_unlabeled_flags[Kn][i])) << " "; // outputing b_i
        else             ostr << smart_round(compute_objective_function_for_SDP(g_unlabeled_flags[Kn][i])) << " "; // outputing b_i
    }
    ostr << endl;
    
    if (upper_bound) ostr << "0 " << (int)g_flags.size()+1 << " 1 1 -1" << endl; // Objective function
    else             ostr << "0 " << (int)g_flags.size()+1 << " 1 1 1" << endl; // Objective function


    mini_timer mt;

    parallel_output_linearizer output_buffer(ostr, 0);

    // printing of subflags
#pragma omp parallel for ordered schedule(dynamic)
    for (int i = 0; i < (int)g_unlabeled_flags[Kn].size(); i++)
    {            
        //cout << "flag: "; c4free_subrgaphs[e][i].print();
        if (verbose_output) 
        {
            #pragma omp critical(cerr)
            {
                std::cerr << "Computing specific part " << i+1 << "/" << g_unlabeled_flags[Kn].size()
                          << mt.report(i, g_unlabeled_flags[Kn].size()) << endl;
            }
        }
        stringstream ss;
        
        int matrixID = i + 1;
     
        // variable for objective function
        double objective_divisor = 1;
        if (g_objective_divisor.size() > 0)
        {
            objective_divisor = smart_round( compute_linear_combination_in_g(g_unlabeled_flags[Kn][i], g_objective_divisor));
            if (objective_divisor != 0)
            {
                if (upper_bound) ss << matrixID << " " << (int)g_flags.size()+1 << " 1 1 " << -objective_divisor << endl;
                else             ss << matrixID << " " << (int)g_flags.size()+1 << " 1 1 " <<  objective_divisor << endl;
            }
        }
        else
        {
            if (upper_bound) ss << matrixID << " " << (int)g_flags.size()+1 << " 1 1 -1" << endl;
            else             ss << matrixID << " " << (int)g_flags.size()+1 << " 1 1 1" << endl;
        }
        
        // variable for slack
        if (has_slack(Kn,i))
        {
            ss << matrixID << " " << (int)g_flags.size()+1 << " " << matrixID+1 << " " <<  matrixID+1 << " 1" << endl;
        }
        
                    
        //int next_block = 
            print_CSDP_constraints_blocks(ss, Kn, i, matrixID, first_constraint_block, verbose_output);                        
// We do not back-up these so whatever order works
// For serious constraints like the cuts heavy programs, this actually may make a differnece.
//#pragma omp ordered
#pragma omp critical
{
            output_buffer.print_string(i, ss.str());
            //ostr << ss.str();
}
    }
    assert(output_buffer.is_empty());

    ostr.flush();
    
    print_CSDP_additional_blocks(Kn, ostr, additional_blocks_offset, verbose_output = 1);
}

void print_CSDP_flag_products_part(int Kn, bool upper_bound = true, ostream &ostr = cout, int verbose_output = 1, bool use_sdp_temp = false, int lowest_not_computed=0)
{
    mini_timer mt;

    parallel_output_linearizer output_buffer(ostr, lowest_not_computed);

#pragma omp parallel for ordered schedule(dynamic)
    for (int i = lowest_not_computed; i < (int)g_unlabeled_flags[Kn].size(); i++)
    {            
        if (verbose_output) 
        {
            #pragma omp critical(cerr)
            {
                std::cerr << "Computing flag products " << i+1 << "/" << g_unlabeled_flags[Kn].size()
                          << mt.report(i, g_unlabeled_flags[Kn].size()) << endl;
            }            
        }

        stringstream ss;
        
        int matrixID = i + 1;
        
#ifdef G_FLAG_PRODUCTS_SLOW_LOW_MEMORY       
         for (int f = 0; f < (int)g_flags.size(); f++)
         {
             for (unsigned int x = 0; x < g_flags[f].size(); x++)
                 for (unsigned int y = x; y < g_flags[f].size(); y++)
                 {
         
                     double PF = P_F1_F2_IN_H(g_flags[f][x], g_flags[f][y], g_unlabeled_flags[Kn][i], true);
                     //					if (x == y) PF = PF*0.5;
                     if (PF != 0)
                     {
                         ss << matrixID <<  " " << f+1 << " " << x+1 << " " << y+1 << " " <<  PF << endl;
                     }
                 }
         }
#else
        count_flag_products(ss, matrixID, g_unlabeled_flags[Kn][i]);
#endif 

    #pragma omp critical (ostr)
    {
        output_buffer.print_string(i, ss.str());
    }

    }
    assert(output_buffer.is_empty());
    ostr.flush();
}

void print_CSDP(int Kn, bool upper_bound = true, ostream &ostr = cout, int verbose_output = 1, bool use_sdp_temp=false, bool sdp_temp_up_to_date=true)
{
    cerr << "Computing specific part of the SDP..." << endl;
    print_CSDP_specific_part(Kn, upper_bound, ostr, verbose_output);
    
    if (!use_sdp_temp || sdp_temp_up_to_date == false)
    {
        cerr << "Computing flag products part of the SDP..." << endl;
        print_CSDP_flag_products_part(Kn, upper_bound, ostr, verbose_output, use_sdp_temp);
        return;
    }
    
    cerr << "Getting flag products for SDP..." << endl;    
    
    stringstream filename;
    
    filename << filename_prefix() << "__n" << Kn << "_sdp_products.txt";

    // The sdp_product_file may not contain all products if the computation
    // was accidentally interrupted. We try to go from we stopped...
    // First find were we stopped
    
    int lowest_not_computed = 0;

    ifstream sdp_product_file;
    sdp_product_file.open (filename.str().c_str(),ifstream::in);
    while (sdp_product_file.good())
    {
        int id, type, f1, f2, d;
        sdp_product_file >> id >> type >> f1 >> f2 >> d;
        if (!sdp_product_file.good()) break; // if the read failed, do not use it...
        if (id > lowest_not_computed)
        {
            lowest_not_computed = id; // note that in file the ID is +1
        }
    }
    sdp_product_file.close();
    
    if (lowest_not_computed < (int)g_unlabeled_flags[Kn].size())
    {
        if (lowest_not_computed != 0)
        {
            cerr << "Warning: file " << filename.str() << " contains only products up to " << lowest_not_computed << " out of " << (int)g_unlabeled_flags[Kn].size() << endl;
            cerr << "          We assume it happened by accident so the rows starting with " << lowest_not_computed << " will be deleted and recomputed" << endl;
            stringstream mv_command;
            mv_command << "rm -f "<< filename.str() << "~ ;  mv " << filename.str() << " " << filename.str() << "~";
            cerr << "          Executing " << mv_command.str() << endl;
            if (system(mv_command.str().c_str()) != 0)
            {
                cerr << " Execution failed" << endl;
                exit(1);
            }
            // We remove the last line since it may be incoplete and not catched by the grep
            stringstream grep_command;
            grep_command << "sed \\$d " << filename.str() << "~ |   grep -v '^" << lowest_not_computed << " ' "  << " >" << filename.str();
            cerr << "          Executing " << grep_command.str() << endl;
            if (system(grep_command.str().c_str()) != 0)
            {
                cerr << " Execution failed" << endl;
                exit(1);
            }
            lowest_not_computed--;
        }
     
        ofstream sdp_product_file_out;
        sdp_product_file_out.open (filename.str().c_str(),ofstream::out|ofstream::app);
        if (!sdp_product_file_out.good())
        {
            cerr << "Failed creating file with sdp products " << filename.str() << endl;
            exit(1);
        }
        
        cerr << "Computing flag products of the SDP to file " << filename.str() << endl;
        print_CSDP_flag_products_part(Kn, upper_bound, sdp_product_file_out, verbose_output, use_sdp_temp, lowest_not_computed);
        sdp_product_file_out.close();
    }

   
    sdp_product_file.open (filename.str().c_str(),ifstream::in);
    if (!sdp_product_file.good())
    {
        cerr << "Failed obtaining file with sdp products " << filename.str() << endl;
        exit(1);
    }
    
    cerr << "Copying flag products of the SDP from " << filename.str() << endl;
    ostr << sdp_product_file.rdbuf();
    sdp_product_file.close();
    
}



void process_csdp_solution(istream &ist, int Kn, double lower_bound, double upper_bound, bool print_density)
{
    // First line is just solution to the dual - Y
	// which mean combination of flags...
	
	double sum = 0;
    double sum_print = 0;
    double count_printed = 0;
	
	for (int i = 0; i < (int)g_unlabeled_flags[Kn].size(); i++)
	{
		double density;
		ist >> density;
		
		sum += density;
		
		if (lower_bound <= density && density < upper_bound)
		{
            count_printed++;
            sum_print += density;
            if (print_density)
            {
                std::cout.fill(' ');
                std::cout.width(G_PRECISION+1);
                std::cout.precision(G_PRECISION);
                std::cout << std::left << density << " ";
            }
            
			//cout << density << " ";
			//if (g_unlabeled_flags[Kn][i].count_crossings() == 2)
			{
                std::cout << g_unlabeled_flags[Kn][i].print() << "  # " << i+1 << endl;
			}
		}
	}
	
	cerr << "Printed: " << count_printed << "/" <<  g_unlabeled_flags[Kn].size() << "  Sum of coefficients: " << sum_print << "/" << sum  << endl;
	cerr << "Numering is 1.." << g_unlabeled_flags[Kn].size() << endl;
	return;
}








// Not paralel version - very slow and not completely correct - in particular for ordered version.
void generate_labeled_flags_old(int flag_size, int type_size, int verbose_output)
{
    int to_label_size = (int)g_unlabeled_flags[flag_size].size();
    //cerr << "Type size " << type_size << " for " << g_unlabeled_flags[flag_size][i].print() << endl;
    int mapping[flag_size];
    flag F;
    
    for (int i = 0; i < to_label_size; i++)
    {
        if (verbose_output)
            cerr << "Flag size "<<flag_size << " type size " << type_size << " labeling " << i << "/" << to_label_size << endl;

        for (int j = 0; j < flag_size; j++) mapping[j]=j;
        
        do {
            F.as_subflag(g_unlabeled_flags[flag_size][i], mapping, flag_size, type_size); // taking a subflag
            include_flag_in_list_if_new(F,g_flags);
        } while ( std::next_permutation(mapping,mapping+flag_size) );
    }
}


//vector<flag> g_types[V];
void create_types_of_size(int type_size, bool check_sensible_flags=false)
{
    if (g_types[type_size].size() != 0) return;


    // Look if there are already some flags with given types...
    for (int i = 0; i < (int)g_flags.size(); i++)
    {
        if (g_flags[i].size() == 0) continue;
        
        if (g_flags[i][0].labeled_vertices_cnt() == type_size)
        {
            flag type;
            g_flags[i][0].get_type_subflag(type);
            g_types[type_size].push_back(type);
        }
    }
    if (g_types[type_size].size() != 0) return;
    
    
    // if not, just label graphs of apropriate size
    get_unlabeled_flags_of_size(type_size);
    
    if (!check_sensible_flags)
    {
        g_types[type_size] = g_unlabeled_flags[type_size];
    }
    else
    {
        for (int i = 0; i < (int)g_unlabeled_flags[type_size].size(); i++)
        {
            g_types[type_size].push_back(g_unlabeled_flags[type_size][i]);
        }
    }
    
#pragma omp parallel for
    for (int i = 0; i < (int)g_types[type_size].size(); i++)
    {
        g_types[type_size][i].make_all_vertices_labeled();
        //
        cerr << "New type: " << g_types[type_size][i].print() << endl;
    }
    
    
    cerr << "Generated " << g_types[type_size].size()  << " types of size " << type_size << endl;
}

int get_reduced_type_id(const flag &f)
{
    int type_size = f.labeled_vertices_cnt();
    
    //cerr << "Testing size " << type_size << endl;

    
    for (int i = 0; i < (int)g_types[type_size].size(); i++)
    {
        if (f.have_same_type(g_types[type_size][i]))
        {
            return i;
        }
    }    
    
    return -1;
}


void get_types_of_size(int type_size, bool types_from_file, bool check_sensible_flags=false)
{
//    create_types_of_size(type_size);
//    return;
        
    if (!types_from_file)
    {
        create_types_of_size(type_size, check_sensible_flags);
        cerr <<  g_types[type_size].size() <<  " types of size " << type_size << " generated" << endl;
        return;
    }
    
    string filename = filename_prefix() + "__type_"+ to_string(type_size) +".txt";
    
    if (!load_flags_from_file(filename, g_types[type_size]))
    {
        create_types_of_size(type_size, check_sensible_flags);
        dump_flags_to_file(filename, g_types[type_size]);
    }
    else
    {
        cerr <<  g_types[type_size].size() <<  " types of size " << type_size << " loaded from file " << filename << endl;
    }
}


//
void generate_labeled_flags(int flag_size, int type_size, int verbose_output, bool types_from_file=false, bool check_sensible_flags=false)
{        
    // First make types
    get_types_of_size(type_size, types_from_file, check_sensible_flags);
    
    //cerr << "B:EEEE  " << g_types[type_size].size() << endl;
    //cerr <<  "flag_size=" << flag_size << " type_size=" << type_size << endl;
    //cerr << "check_sensible_flags=" << check_sensible_flags << endl;
    
    int to_label_size = (int)g_unlabeled_flags[flag_size].size();
    
    // 
    vector< vector<flag> >  *flags;
    flags = new vector< vector<flag> >[to_label_size];
    
    int types = (int)g_types[type_size].size();
    
    // Avoid copying flags
    for (int i = 0; i < to_label_size; i++)
    {
        flags[i].resize(types);
    }
    
    cerr << "Generating labeled flags for Cauchy-Swartz in parallel " << endl;
#pragma omp parallel for ordered
    for (int i = 0; i < to_label_size; i++)
    {
        if (verbose_output)
            cerr << "Flag size "<<flag_size << " type size " << type_size << " labeling " << i << "/" << to_label_size << endl;

        //cerr << "Type size " << type_size << " for " << g_unlabeled_flags[flag_size][i].print() << endl;
        // Permutations for labeling... Could be improved to distinguish only the first type_size vertices.
        int mapping[flag_size];
        flag F;
        for (int j = 0; j < flag_size; j++) mapping[j]=j;
        
        //int found_good_type = 0;
        do {
            F.as_subflag(g_unlabeled_flags[flag_size][i], mapping, flag_size, type_size); // taking a subflag
            //cerr << "Got one" << endl;
            int type_of_F = get_reduced_type_id(F);
            if (type_of_F == -1)
            {
                // We have each type only once -  3 3   1 1  2  and  3 3   1 2  1 are different types but enough to use just one
                continue;
            }
            //found_good_type++;
            if (!g_already_in_known_flags(F,flags[i][type_of_F]))
            {
                if (check_sensible_flags)
                {                    
                    bool flag_OK = true;
                    if (!flag_OK) continue;
                }

                
                    flags[i][type_of_F].push_back(F);
            }
            //include_flag_in_list_if_new(F,flags[i],false);
        } while ( std::next_permutation(mapping,mapping+flag_size) );
        //cout << "Found good type " << found_good_type << endl;
    }
    
    // Todo: this might be also maybe paralelized - there is some paralelism inside...
    for (int i = 0; i < to_label_size; i++)
    {
        for (int j = 0; j < (int)flags[i].size(); j++)
            for (int k = 0; k < (int)flags[i][j].size(); k++)
                include_flag_in_list_if_new(flags[i][j][k], g_flags);
        
        //merge_typed_flags_lists_parallel(g_flags, flags[i]);
    }
    
    for (int t = 0; t < types; t++)
    {
        int found_t = 0;
        for (int i = 0; i < to_label_size; i++)
        {
            found_t +=  (int)flags[i][t].size();
        }        
        if (found_t == 0)
        {
            cerr << "Type " << t << " has no labeled flags" << endl;
        }
    }
  

 

    delete[] flags;

    //cerr << "done " << endl;
}


//////////////////////////////////////////////////////// LOAD & DUMP
//////////////////////////////////////////////////////// LOAD & DUMP
//////////////////////////////////////////////////////// LOAD & DUMP
//////////////////////////////////////////////////////// LOAD & DUMP
//////////////////////////////////////////////////////// LOAD & DUMP
//////////////////////////////////////////////////////// LOAD & DUMP


// return true if constraint added
bool add_linear_constraint_if_new(linear_constraint &lc, bool test_with_type_permutations, vector<linear_constraint> &where_to_add = g_linear_constraints, bool use_openmp = false)
{
    if (use_openmp == false)
    {
        for (const linear_constraint &lc_known : where_to_add)
        {
            if (lc_known == lc)
            {
                //cerr << "Duplicate linear constraint..." << endl;
                return false;
            }

            if (test_with_type_permutations && lc.is_identical_after_type_permutation(lc_known))
            {
                //cerr << "Duplicate after relabeling found" << endl;
                // assert(0);
                return false;
            }
        }
    }
    else // parallel version
    {
        bool constraint_is_new = true;

        #pragma omp parallel for ordered schedule(dynamic)
        for (int i = 0; i < (int)where_to_add.size(); i++)
        {
            const linear_constraint &lc_known = where_to_add[i];
            if (constraint_is_new == false) continue;

            if (lc_known == lc)
            {
                #pragma omp critical
                constraint_is_new = false;
                continue;
            }

            if (test_with_type_permutations && lc.is_identical_after_type_permutation(lc_known))
            {
                #pragma omp critical
                constraint_is_new = false;
            }
        }
        if (constraint_is_new == false)
            return false;
    }


    where_to_add.push_back(lc);  
    return true;
}

void remove_constraints_implied_by_others(int verbose_output)
{
    int to_be_removed_counter = 0;


    std::cerr << "Testing useless constraints. To test " << g_linear_constraints.size() << endl;

    vector<linear_constraint> all_linear_constraints(std::move(g_linear_constraints));
    g_linear_constraints.clear();
    g_linear_constraints.reserve(all_linear_constraints.size());


    mini_timer mt;

//    for (const linear_constraint &lc_to_test : all_linear_constraints)
//    const linear_constraint &lc_to_test = g_linear_constraints[0];
    #pragma omp parallel for ordered schedule(dynamic)
    for (int i = 0; i < (int)all_linear_constraints.size(); i++)
    {

        if (verbose_output)
        {
            #pragma omp critical(cerr)
            {
                cerr << "Testing constraint " << i << " out of " << all_linear_constraints.size()
                     << mt.report(i, all_linear_constraints.size()) << endl;
            }
        }

        const linear_constraint &lc_to_test = all_linear_constraints[i];

        bool useful_constraint = true;

        for (const linear_constraint &lc : all_linear_constraints)
        {
            // skip icomparable constriants or is we know that
            // lc cannot be < lc_to_test
            if (lc.m_constant > lc_to_test.m_constant) continue;
            if (lc.m_labeled_vertices_in_type_cnt != lc_to_test.m_labeled_vertices_in_type_cnt) continue;
            if (lc.m_entries_max_size > lc_to_test.m_entries_max_size) continue;
            if (lc.m_entries.size() > lc_to_test.m_entries.size()) continue;
            if (!lc_to_test.m_type.is_isomorphic_to(lc.m_type)) continue;


            // For every flag in lc_to_test, its coefficient in lc is >=
            bool found_one_strict_less = false;
            bool lc_is_not_smaller = false;
            int  matches_not_found = 0;

            for (const flag_and_coefficient &entry_test : lc_to_test.m_entries)
            {
                // Test if it is the same as some other entry in lc
                bool found_match = false;
                for (const flag_and_coefficient &entry2 : lc.m_entries)
                {
                    if (entry_test.g.is_isomorphic_to(entry2.g))
                    {
                        if (entry_test.coefficient < entry2.coefficient)
                        {
                            lc_is_not_smaller = true;
                            break;
                        }

                        if (entry_test.coefficient > entry2.coefficient)
                        {
                            found_one_strict_less = true;
                        }

                        found_match = true;
                        break;
                    }
                }
                if (found_match == false) 
                {
                    found_one_strict_less = true;
                    matches_not_found++;
                }
                if (lc_is_not_smaller) 
                {
                    break;
                }
            }

            if (lc_is_not_smaller) continue;
            if (!found_one_strict_less) continue;
            if (lc.m_entries.size()+matches_not_found != lc_to_test.m_entries.size()) continue;

            if (verbose_output >= 2)
            {
                // We found and extra one
                cerr << "I found a linear constraint that is not needed" << endl;
                cerr << lc_to_test << "is weaker than \n" << lc << endl;
            }
            to_be_removed_counter++;
            useful_constraint = false;
            break;
        }

        #pragma omp ordered
        if (useful_constraint)
        {
            // This is causing a trouble since it can reorder the constraints
            //#pragma omp critical (g_linear_constraints)
            {
                g_linear_constraints.push_back(lc_to_test);
            }
        }
    }
    cerr << "Found " << to_be_removed_counter << " not needed constraints." << endl;
    cerr << "Current number of constraints is " << g_linear_constraints.size() << endl;
}

void simplify_FC(vector<flag_and_coefficient> &FC, double epsilon=0)
{
    vector<flag_and_coefficient> new_FC;
    
    for (int j = 0; j < (int)FC.size(); j++)
    {
        bool processed = false;
        for (int k = 0; k < (int)new_FC.size(); k++)
        {
            if (FC[j].g.is_isomorphic_to(new_FC[k].g))
            {
                new_FC[k].coefficient += FC[j].coefficient;
                processed = true;
                break;
            }
        }
        if (!processed)
        {
            flag_and_coefficient fc = FC[j];
            new_FC.push_back(fc);
        }
    }
    
    FC.clear();
    for (int j = 0; j < (int)new_FC.size(); j++)
    {
        if (abs(new_FC[j].coefficient) > epsilon)
        {
            new_FC[j].coefficient = smart_round(new_FC[j].coefficient, epsilon);
            FC.push_back(new_FC[j]);
        }
    }
    
//    FC = new_FC;
}

void load_linear_constraints_from_stream(istream &instream, bool test_duplicates_with_type_permutations, bool remove_implied_constraints, int verbose_output)
{
    //
    flag_and_coefficient fc;
    string lc_constat;

    int removed_duplicate_constraints = 0;
    int processed_constraints = 0;

    while (instream)
    {
        linear_constraint lc;
        instream >> lc.m_constant;
        if(!instream.good()) {
            // Not a number means comment continues...
            std::istreambuf_iterator<char> eos;
            std::string s(std::istreambuf_iterator<char>(instream), eos);
            g_program_description = s;
            break;
        }
        processed_constraints++;
        instream >> fc.coefficient;
        
        if (verbose_output >= 3)
            cerr << "Loading linear constraint 0 <= " << lc.m_constant <<  endl;
                
        assert (fc.coefficient != 0);
        
        while (instream && fc.g.load_from_stream(instream,-1,-1))
        {
            if (verbose_output >= 3)
                cerr << "..+..  " << fc.coefficient << " * " << fc.g.print() << endl;
            
            lc.m_entries.push_back(fc);
            instream >> fc.coefficient;

            if (fc.coefficient == 0)
            {
                break;
            }
        }

        // Sometimes linear constraints may contain same entry several times
        // so we just sum these all up
        simplify_FC(lc.m_entries);
        
        if (lc.m_entries.size() == 0 && lc.m_constant == 0)
        {
            if (verbose_output >= 2)
            {
                cerr << "Constraint empty after simplification" << endl;
                continue;
            }
        }
        
        lc.check_constraint();
        if (lc.m_checked == false)
        {
            cerr << "check_constraint() failed for constraint:" << endl;
            cerr << lc << endl;
        }
        assert(lc.m_checked);
        if (!add_linear_constraint_if_new(lc, test_duplicates_with_type_permutations,g_linear_constraints, true))
        {
            removed_duplicate_constraints++;
        }

        if (verbose_output >= 2)
            cerr << "Processed " << processed_constraints << "  constraints, found " << removed_duplicate_constraints << " duplicates " <<  endl;


    }
    if (removed_duplicate_constraints > 0)
    {
        cerr << "Removed " << removed_duplicate_constraints << " duplicate constriants." << endl;
    }

    if (remove_implied_constraints)
    {
        remove_constraints_implied_by_others(verbose_output);
    }
}


bool load_linear_constraints_from_file(const string &filename, bool test_duplicates_with_type_permutations, int verbose_output)
{
    ifstream infileF;
    infileF.open (filename.c_str(), ifstream::in);

    if (infileF.good())
    {
        FilteringIstream infile(infileF);
    
        cerr << "Loading additional constraints from file " << filename  <<  endl;
        int tmp;
        infile >> tmp;
        assert(tmp == 0);
        load_linear_constraints_from_stream(infile, test_duplicates_with_type_permutations, test_duplicates_with_type_permutations, verbose_output);
        infileF.close();
        return true;
    }

    cerr << "Failed opening file with objective function " << filename << endl;

    return false;
}


bool load_objective_from_file(const string &filename, bool test_duplicates_in_linear_constriants_with_type_permutations,  int verbose_output)
{
    OPEN_FILE_SMARTLY_RETURN_FALSE_ON_FAIL(istr, filename);
    /*
    ifstream infileF;
    infileF.open (filename.c_str(), ifstream::in);
    if (!infileF.good())
    {
        cerr << "Failed opening file with objective function " << filename << endl;
        return false;
    }
    
    FilteringIstream infile(infileF);
    */

    cerr << "Loading objective function from file " << filename << endl;
        
    flag_and_coefficient fc;
    
    (*istr) >> fc.coefficient;

    // THis is useful if we want to draw the file of constraints.
    if (fc.coefficient == 0)
    {
        load_linear_constraints_from_stream((*istr), test_duplicates_in_linear_constriants_with_type_permutations, test_duplicates_in_linear_constriants_with_type_permutations, verbose_output);
    } 
    else
    { 
    while ((*istr) && fc.g.load_from_stream((*istr),-1,-1))
    {
        // check if the flag is labeled
        if (fc.g.labeled_vertices_cnt() != 0)
        {
            cerr << "WARNING: You are using labeled flags in th objective function" << endl;
        }
        // check if there are no repetitions in the objective function..
        for (int i = 0; i < (int)g_objective_combination.size(); i++)
        {
            if (g_objective_combination[i].g.is_isomorphic_to(fc.g))
            {
                cerr << "WARNING: In the objecive function, flags " << g_objective_combination[i].g.print() << " and " <<
                fc.g.print() << " are the same " << endl;
                //exit(1);
            }
        }
        g_objective_combination.push_back(fc);
        (*istr) >> fc.coefficient;
        
        if (fc.coefficient == 0)
        {
            load_linear_constraints_from_stream((*istr), test_duplicates_in_linear_constriants_with_type_permutations, test_duplicates_in_linear_constriants_with_type_permutations, verbose_output);
            break;
        }
    }
    }
        
    cerr << "# of entries in objective is " << g_objective_combination.size() << endl;
    cerr << "# of linear constraints is " << g_linear_constraints.size() << endl;
    
    return true;
}





bool load_labeled_flags_from_file(int sizeKn, int verbose_output = 0)
{
    stringstream filename;
    filename << filename_prefix() << "__n" << sizeKn << "_labeled.txt";
    
    ifstream infileF;
    infileF.open (filename.str().c_str(), ifstream::in);
    if (!infileF.good())
    {
        cerr << "Failed opening file with labeled flags " << filename.str() << endl;
        return false;
    }
    
    FilteringIstream infile(infileF);

/*
    ifstream infile;
    infile.open (filename.str().c_str(), ifstream::in);
    if (!infile.good())
    {
        cerr << "Failed opening file with labeled flags " << filename.str() << endl;
        return false;
    }
*/

    if (verbose_output)
    {
        cerr << "Loading labeled flags from file " << filename.str() << endl;
    }

    
   flag f;
	while (f.load_from_stream(infile,-1,-1))
	{
		include_flag_in_list(f, g_flags);
	}
	
	infileF.close();
    
	return true;
}

bool dump_labeled_flags(int sizeKn)
{
    stringstream filename;
    filename << filename_prefix() << "__n" << sizeKn << "_labeled.txt";
    
    ofstream outfile;
    outfile.open (filename.str().c_str(), ofstream::out);
    if (!outfile.good())
    {
        cerr << "Failed opening file " << filename.str() << endl;
        return false;
    }
    
    cerr << "Writing labeled flags to file " << filename.str() << endl;
    
    for (int f = 0; f < (int)g_flags.size(); f++)
    {
        for (unsigned int x = 0; x < g_flags[f].size(); x++)
        {
            outfile << g_flags[f][x].print() << endl;
        }
        if (f < (int)g_flags.size()-1) outfile << endl;
    }
    
    outfile.close();
    
    return true;
}



bool compare_flag_sizes(const flag&f1, const flag &f2)
{
    return f1.m_vertices < f2.m_vertices;
}

void sort_flags_by_size(vector<flag> &flag_list)
{
    sort(flag_list.begin(), flag_list.end(), compare_flag_sizes);
}



void load_forbidden()
{
    stringstream filename;
    filename <<  filename_prefix() << "__forbidden.txt";
    
    ifstream infile;
    infile.open (filename.str().c_str(), ifstream::in);
    if (!infile.good())
    {
        cerr << "Failed opening file " << filename.str() << " no forbidden structures are used" << endl;
        return;
    }
    
    FilteringIstream filteredinfile(infile);


    int flags_in_file = 0;
    
    flag h;
    while (h.load_from_stream(filteredinfile,-1,0))
    {
        flags_in_file++;
        // check if h is not already there
        if (!(g_already_in_known_flags(h,g_forbidden_subflags)))
        {
            g_forbidden_subflags.push_back(h);
            g_forbidden_subflags_by_size[h.m_vertices].push_back(h);
            //cerr << "Loaded forbidden graph " << h.print() << endl;
        }
    }
    
    infile.close();
    
    // Idea is that it is faster to test smaller subflags than bigger
    // so lets first test smaller ones
    sort_flags_by_size(g_forbidden_subflags);
    
    cerr << "Loaded " << g_forbidden_subflags.size() << " forbidden graphs from " << filename.str() << endl;
    if (flags_in_file != (int)g_forbidden_subflags.size())
    {
        cerr << "File " << filename.str() << " contains " << flags_in_file-g_forbidden_subflags.size() << " duplicate  graphs." << endl;    
    }
    
    //for (int i = 0; i < g_forbidden_subflags.size(); i++)
    //{
    //    cerr << g_forbidden_subflags[i].print() << endl;
    //}
}









void add_linear_constraints_flags_equal(const flag &f1, const flag &f2)
{
    linear_constraint lc;
    lc.m_constant = 0;
    
    flag_and_coefficient gd;
    gd.coefficient = 1;
    gd.g = f1;
    lc.m_entries.push_back(gd);
    
    gd.coefficient = -1;
    gd.g = f2;
    lc.m_entries.push_back(gd);
    
    lc.check_constraint();
    assert(lc.m_checked);
    g_linear_constraints.push_back(lc);

    // swap the coefficients and also add as a linear constraint
    lc.m_entries[0].coefficient = -1;
    lc.m_entries[1].coefficient = 1;
    lc.check_constraint();
    assert(lc.m_checked);
    g_linear_constraints.push_back(lc);
}






void extensions_of_g_with_print(flag &g, bool extension_count_copies)
{
    vector<flag> flag_list;
    
    extensions_of_g(g, flag_list);
    
    for (unsigned int x = 0; x < flag_list.size(); x++)
    {
        if (extension_count_copies)
        {
            cout << flag_list[x].count_labeled_copies_of(g) << "    ";
        }
        cout << flag_list[x].print() << endl;
    }
}


void extensions_of_fc(const string &filename)
{
    cerr << "Extensions of " << filename << endl;
    vector<flag_and_coefficient> F1;
    load_flags_and_coefficients_from_file(filename, F1);

    vector<flag_and_coefficient> FE; 

    for (int i = 0; i < (int)F1.size(); i++)
    {
        vector<flag> flag_list;
        extensions_of_g(F1[i].g, flag_list);
    
        for (unsigned int x = 0; x < flag_list.size(); x++)
        {
            flag_and_coefficient fc;
            fc.g = flag_list[x];
            fc.coefficient = F1[i].coefficient;

            dump_flag_and_coefficient(fc);
            //FE.push_back(fc);
        }
    }
        
    //dump_flags_and_coefficients(FE);
}

void generate_products_of_constraints(int Kn, string objective_file_name, bool force_generating_constriants, int verbose_output)
{
    if (g_linear_constraints.size() == 0) return;
    
    
    stringstream filename;
    filename <<  filename_prefix() << "__n" << Kn << "_generated_constraints_" << objective_file_name << ".txt";
    
    if (!force_generating_constriants)
    {
        if (load_linear_constraints_from_file(filename.str(), false, verbose_output)) return;
    }
    
    
    cerr << "Generating additoinal linear constraints...." << endl;
    // 
    int original_size = (int)g_linear_constraints.size();
    for (int i = 0; i < original_size; i++)
    {
        //if (!g_linear_constraints[i].m_same_types) continue;
        for (int j = i; j < (int)g_linear_constraints.size(); j++)
        {
            
            //if (!g_linear_constraints[j].m_same_types) continue;

            linear_constraint &lci = g_linear_constraints[i];
            linear_constraint &lcj = g_linear_constraints[j];

            
            if (!lci.m_type.have_same_type(lcj.m_type)) continue;
            
            if (lci.m_entries_max_size + lcj.m_entries_max_size - lci.m_labeled_vertices_in_type_cnt > Kn) continue;
            
            cerr << "Adding product of constraints " << i+1 << " and " << j+1 << endl;
            linear_constraint lc;
            
            // Does product of (F1+c1)*(F1+c2)
            
            // c1*c2
            lc.m_constant = lci.m_constant * lcj.m_constant;
            
            flag_and_coefficient fc;
            // F1*c2
            for (int x = 0; x < (int)lci.m_entries.size(); x++)
            {
                fc = lci.m_entries[x];
                fc.coefficient *= lcj.m_constant;
                lc.add_entry(fc);
            }
            // F2*c1
            for (int y = 0; y < (int)lcj.m_entries.size(); y++)
            {
                fc = lcj.m_entries[y];
                fc.coefficient *= lci.m_constant;
                lc.add_entry(fc);
            }
            // F1*F2
            for (int x = 0; x < (int)lci.m_entries.size(); x++)
            {
                for (int y = 0; y < (int)lcj.m_entries.size(); y++)
                {
                    vector<flag_and_coefficient> F1F2 = F1_times_F2(lci.m_entries[x].g, lcj.m_entries[y].g, lci.m_type);
                    
                    double scale = lci.m_entries[x].coefficient * lcj.m_entries[y].coefficient;
                    
                    // polish the list of products and add them to lc
                    for (int z = 0; z < (int)F1F2.size(); z++)
                    {
                        F1F2[z].coefficient *= scale;
                        lc.add_entry(F1F2[z]);
                    }                    
                }              
            }

            lc.check_constraint();
            assert(lc.m_checked);
            g_linear_constraints.push_back(lc);
            if (verbose_output)
                cerr << "Adding new linear constraint " << endl << lc.print();
        }
    }
    
    cerr << "Generated " << g_linear_constraints.size()-original_size << " new linear constraints" << endl;
    cerr << "Writing new linear constraints to " << filename.str() << endl;
    
    ofstream outfile;
    outfile.open (filename.str().c_str(), ofstream::out);
    if (!outfile.good())
    {
        cerr << "Failed opening file " << filename.str() << endl;
        return;
    }
        
    for (int c = original_size; c < (int)g_linear_constraints.size(); c++)
    {
        outfile << g_linear_constraints[c].print() << endl;        
    }
    
    outfile.close();
}


void scale_linear_constraints_to_integers(int Kn)
{
    if (g_linear_constraints.size() == 0) return;
    
    cerr << "Scaling linear constraints...." << endl;
    for (int j = 0; j < (int)g_linear_constraints.size(); j++)
    {
        // Computing the right scale
        double scale=1;
        bool used[V];
        for (int i = 0; i < V; i++) used[i] = false;
        for (int k = 0; k < (int)g_linear_constraints[j].m_entries.size(); k++)
        {
            int v = g_linear_constraints[j].m_entries[k].g.m_vertices;
            if (v < Kn && used[v] == false)
            {
                used[v] = true;
                scale *= binomial(Kn,v);
            }
        }
        
        if (scale != 1)
        {
            cerr << "Scaling linear constraint " << j+1 << " by " << scale << endl;
            // Applying the right scale
            g_linear_constraints[j].m_constant *= scale;
            for (int k = 0; k < (int)g_linear_constraints[j].m_entries.size(); k++)
            {
                g_linear_constraints[j].m_entries[k].coefficient *= scale;
            }
        }
        else
        {
            cerr << "Constraint " << j+1 << " not scaled" << endl;            
        }
    }   
}


// Generates Baber's equalities
void generate_baber_equalities(int Kn, int verbose_output)
{
    cerr << "Baber's inequalities were not enabled during compile time" << endl;
}


void print_latex_header(ostream &outfile, bool color_1_nonedge = false)
{
   
    outfile << ""
    "\\documentclass{article} \n"
    "\\usepackage{tikz}\n"
    "\\usepackage{fullpage} \% Disable if equations are overflowing \n"
    "\\usepackage{multicol}\n"
    "\n \% tikz style \n\n"
    "\\newcommand{\\vc}[1]{\\ensuremath{\\vcenter{\\hbox{#1}}}}\n"
    "\\tikzset{flag_pic/.style={scale=1}}  \%  Change the scale to change all figures \n"
    "\\tikzset{unlabeled_vertex/.style={inner sep=1.7pt, outer sep=0pt, circle, fill}} \n"
    "\\tikzset{labeled_vertex/.style={inner sep=2.2pt, outer sep=0pt, rectangle, fill=yellow, draw=black}} \n"
    "\\tikzset{edge_color0/.style={color=black,line width=1.2pt,opacity=0.5}} \n"
    ;
    if (color_1_nonedge)
    {
        outfile << "\\tikzset{edge_color1/.style={color=red,  line width=1.2pt,opacity=0}} \n";
    }
    else
    {
        outfile << "\\tikzset{edge_color1/.style={color=red,  line width=1.2pt,opacity=1}} \n";
    }
    outfile <<
    "\\tikzset{edge_color2/.style={color=blue, line width=1.2pt,opacity=1}} \n"
    "\\tikzset{edge_color3/.style={color=green,line width=1.2pt}} \n"
    "\\tikzset{edge_color4/.style={color=red,  line width=1.2pt,dotted}} \n"
    "\\tikzset{edge_color5/.style={color=blue, line width=1.2pt,dotted}} \n"
    "\\tikzset{edge_color6/.style={color=green, line width=1.2pt,dotted}} \n"
    "\\tikzset{edge_color7/.style={color=orange, line width=1.2pt}} \n"
    "\\tikzset{edge_color8/.style={color=gray, line width=1.2pt}} \n"
    "\\tikzset{edge_thin/.style={color=black}} \n"    
    "\\tikzset{edge_hidden/.style={color=black,dotted,opacity=0}} \n"    
    "\\tikzset{vertex_color1/.style={inner sep=1.7pt, outer sep=0pt, draw, circle, fill=red}} \n"
    "\\tikzset{vertex_color2/.style={inner sep=1.7pt, outer sep=0pt, draw, circle, fill=blue}} \n"
    "\\tikzset{vertex_color3/.style={inner sep=1.7pt, outer sep=0pt, draw, circle, fill=green}} \n"
    "\\tikzset{vertex_color4/.style={inner sep=1.7pt, outer sep=0pt, draw, circle, fill=pink}} \n"
    "\\tikzset{labeled_vertex_color1/.style={inner sep=2.2pt, outer sep=0pt, draw, rectangle, fill=red}} \n"
    "\\tikzset{labeled_vertex_color2/.style={inner sep=2.2pt, outer sep=0pt, draw, rectangle, fill=blue}} \n"
    "\\tikzset{labeled_vertex_color3/.style={inner sep=2.2pt, outer sep=0pt, draw, rectangle, fill=green}} \n"
    "\\tikzset{labeled_vertex_color4/.style={inner sep=2.2pt, outer sep=0pt, draw, rectangle, fill=pink}} \n"
    "\n"
    "\\def\\outercycle#1#2{ \n"
    "\\pgfmathtruncatemacro{\\plusone}{#1+1} \n"
// \def\outercycle#1#2{ 
//\pgfmathtruncatemacro{\plusone}{#2+1}
//\pgfmathtruncatemacro{\shift}{270- (360/#2)/2}
//\draw \foreach \x in {0,1,...,\plusone}{(\shift+\x*360/#2:0.9) coordinate(x\x)};
    "\\pgfmathtruncatemacro{\\zeroshift}{270 - (#2-1)*360/#1/2 } \n"
//    "\\draw  \\foreach \\x in {0,1,...,#1}{(270-360/#1/2+\\x*360/#1:1) coordinate(x\\x)};} \n"
    "\\draw  \\foreach \\x in {0,1,...,#1}{(\\zeroshift+\\x*360/#1:1) coordinate(x\\x)};} \n"
//    "\\foreach \\x in {0,1,...,#1}{(270-45+\\x*360/#2:1) coordinate(x\\x)};} \n"
    "\\def\\drawhypervertex#1#2{ \\draw[edge_color2] (x#1)++(0,-0.2-0.2*#2)+(-0.2,0) -- +(0.2,0);} \n"
    "\\def\\drawhypervertexcolor#1#2#3{ \\draw[edge_color#3] (x#1)++(0,-0.2-0.2*#2)+(-0.2,0) -- +(0.2,0);} \n"
    "\\def\\drawhyperedge#1#2{ \\draw[dotted] (x0)++(0,-0.2-0.2*#1)--++(0.5*#2-0.5,0);} \n"
    "\n"
    "%\\def\\labelvertex#1{\\draw (x#1) node[below]{#1}; }\n"
    "\\def\\labelvertex#1{\\pgfmathtruncatemacro{\\vertexlabel}{#1+1 } \\draw (x#1) node{\\tiny\\vertexlabel}; }\n"
    "%\\def\\labelvertex#1{}"
    "\n\n\\begin{document}\n"
    << endl;    
}

// A human friendly printout that can be process by latex - makes it easier
// to check if the program undersood input. It is not implemented completely.
//
void print_problem_in_latex(string objective_file_name, bool color_1_nonedge = false)
{
    stringstream filename;
    if (objective_file_name != "")
        filename << objective_file_name << "__latex.tex";
    else 
        filename << filename_prefix() << "__latex.tex";
    
    ofstream outfile;
    outfile.open (filename.str().c_str(), ofstream::out);
    if (!outfile.good())
    {
        cerr << "Failed opening file " << filename.str() << endl;
        return;
    }
    
    cerr << "Writing problem formulation in latex to a file " << filename.str() << endl;
    
    print_latex_header(outfile, color_1_nonedge);


    outfile << "Problem: \\verb!" << filename_prefix() << "!\n\n";

    outfile << "Description:\n  \\begin{verbatim}" << g_program_description << "\\end{verbatim}\n";
    
    outfile << "Compile options: \\begin{itemize}\n"
#ifdef G_COLORED_EDGES
        "\\item Edges have "<< G_COLORED_EDGES << " colors\n"
        "\\item Color coding:  0 black, 1 red, 2 blue, 3 green, 4 red dotted, 5 blue dotted, 6 green dotted\n"
#endif
#ifdef G_COLORED_EDGES_BLIND 
        "\\item Edges are treated colorblind TODO write color permutations\n"
#endif
    
    "\\end{itemize}\n";
    

    if (objective_file_name == "")
        outfile << "\n\n Objective file name: default \n";
    else    
        outfile << "\n\n Objective file name: \\verb!" << objective_file_name << "!\n";

    
    outfile << "\n\n Forbidden flags: " << g_forbidden_subflags.size() << " \\begin{center}\n";
    for (int i = 0; i < (int)g_forbidden_subflags.size(); i++)
    {
        outfile << g_forbidden_subflags[i].print_latex(false,0);
    }
    outfile << "\\end{center}\n";


    outfile << "\n\n Objective is a linear combination of  " << g_objective_combination.size() << " flag(s): \\begin{center}\n";
    for (int i = 0; i < (int)g_objective_combination.size(); i++)
    {
        outfile << std::setprecision(G_PRECISION);
        if (g_objective_combination[i].coefficient >= 0)
            outfile << "+" << g_objective_combination[i].coefficient << " \\vc{" << g_objective_combination[i].g.print_latex(false,0) << "}";
        else
            outfile << g_objective_combination[i].coefficient << " \\vc{" << g_objective_combination[i].g.print_latex(false,0) << "}";            
    }
    outfile << "\\end{center}\n";

    
    outfile << "\\newpage \n\n Linear constraints count : " << g_linear_constraints.size() << "\n";
    for (int i = 0; i < (int)g_linear_constraints.size(); i++)
    {
        outfile << " \\begin{center}\n";
        outfile << " $0 \\leq $";
//        outfile << std::setprecision(G_PRECISION) << g_linear_constraints[i].m_constant;
        outfile << std::setprecision(G_PRECISION) << g_linear_constraints[i].m_constant;
        for (int j = 0; j < (int)g_linear_constraints[i].m_entries.size(); j++)
        {
            if (g_linear_constraints[i].m_entries[j].coefficient >= 0)
                outfile << " +" << g_linear_constraints[i].m_entries[j].coefficient << "~\\vc{" << g_linear_constraints[i].m_entries[j].g.print_latex(false,0) << "}";  
            else
                outfile << " " << g_linear_constraints[i].m_entries[j].coefficient << "~\\vc{" << g_linear_constraints[i].m_entries[j].g.print_latex(false,0) << "}";
        }
        outfile << "\\end{center}\n";
    }
    
//    string print_latex(bool use_label, const T &graph_label) const;

    outfile << "\\end{document}\n";
    outfile.close();
    
    execlp("pdflatex","pdflatex",filename.str().c_str(),(char *)NULL);    
}




void draw_graphs(const string& path, bool file_has_with_densities, bool color_1_nonedge = false,  bool latex_enumerated = false, bool use_cout = false)
{
    OPEN_FILE_SMARTLY(istr, path);

    
    ostream *ostr = &std::cout;

    stringstream filename;
    filename << path << "__latex.tex";
    
    ofstream outfile;

    if (use_cout == false)
    {
        outfile.open (filename.str().c_str(), ofstream::out);
        if (!outfile.good())
        {
            cerr << "Failed opening file " << filename.str() << endl;
            return;
        }
        ostr = &outfile;
    }


    
    //cerr << "Writing problem formulation in latex to a file " << filename.str() << endl;
    print_latex_header((*ostr), color_1_nonedge);    
    //outfile << "Problem: \\verb!" << filename_prefix() << "!\n\n";

    
    
    /*
    vector<flag> drawn;
    
    flag h;
    while (h.load_from_stream(infile,-1,0))
    {
        if (find_flag_in_list(h, drawn) < 0)
        {
            
            outfile << h.print_latex(false,0);
            drawn.push_back(h);
        }
    }
    */
    
    if (latex_enumerated)
    {
        (*ostr) << "%\\begin{multicols}{2}" << endl;
        (*ostr) << "\\begin{enumerate}" << endl;
    }
    
    double coefficient = 0.0;
    flag h;
    if (file_has_with_densities)
    {
        (*istr) >> coefficient;
    }
    while (h.load_from_stream((*istr),-1,-1))
    {
        //if (!is_flag_forbidden(h))
        {

            (*ostr) << "% " << h.print() << endl;

            if (latex_enumerated)
            {
                (*ostr) << "\\item" << endl;
            }

            if (file_has_with_densities)
            {
                //cerr << "Printing " << coefficient << endl;
                //if (abs(coefficient) < 0.0002)
                //    outfile << 0;
                ///else 
                    (*ostr).precision(G_PRECISION);
                    (*ostr) << coefficient << " ";
            }
            
            (*ostr) << "\\vc{"<<h.print_latex(false,0)<<"}\n" ;
            //g_forbidden_subflags.push_back(h);
        }
        if (file_has_with_densities)
        {
            (*istr) >> coefficient;
        }
    }
    

    if (latex_enumerated)
    {
        (*ostr) << "\\end{enumerate}" << endl;
        (*ostr) << "%\\end{multicols}" << endl;
    }

    (*ostr) << "\\end{document}\n";

    if (use_cout == false)
    {
        (*ostr).flush();
        execlp("pdflatex","pdflatex",filename.str().c_str(),(char *)NULL);
    }
}



void compute_densities_in(const string& path, int Kn)
{
    vector<flag_and_coefficient> big_flags;

    OPEN_FILE_SMARTLY(istr, path);


    flag_and_coefficient fc;
    
    (*istr) >> fc.coefficient;
    while (fc.g.load_from_stream((*istr),-1,0))
    {
        big_flags.push_back(fc);
        //cout << fc.coefficient << " " << fc.g.print() << endl;
        (*istr) >> fc.coefficient;
    }

    cerr << "Loaded " << big_flags.size() << " flags from " << path << endl;
    
    double all_densitysum = 0;
    
    for (int i = 0; i < (int)g_unlabeled_flags[Kn].size(); i++)
    {
        double densitysum = 0;
        for (int j = 0; j < (int)big_flags.size(); j++)
        {
            if (big_flags[j].g.m_vertices >= Kn)
            {              
                densitysum += big_flags[j].coefficient*P_F1_IN_H(g_unlabeled_flags[Kn][i], big_flags[j].g);
            }
            else
            {
                densitysum += big_flags[j].coefficient*P_F1_IN_H(big_flags[j].g, g_unlabeled_flags[Kn][i]);                
            }
        }
        std::cout.precision(G_PRECISION);
        cout << densitysum << " " <<  g_unlabeled_flags[Kn][i].print() << " # " <<  i+1 << endl;
        all_densitysum += densitysum;
    }
     std::cout.precision(G_PRECISION);
    cerr << "Sum of all densities is " << all_densitysum << endl;
}

// filter_allowed means that the pathfilter contain allowed graphs
// otherwise is contains forbidden graphs
void filter_flags(const string& pathall, const string& pathfilter, bool filter_allowed, bool coefficients_in_input)
{
    vector<flag> filter_flags;
    load_flags_from_file(pathfilter, filter_flags);
    cerr << "Loaded " << filter_flags.size() << " flags ";
    if (filter_allowed) cerr << "that are allowed" << endl;
    else cerr << "that are forbidden" << endl;
    

    OPEN_FILE_SMARTLY(istr, pathall);


    flag_and_coefficient fc;

    if (coefficients_in_input) (*istr) >> fc.coefficient;
    while (fc.g.load_from_stream((*istr),-1,-1))
    {
        if (find_flag_in_list(fc.g,filter_flags) != -1)
        {
            if (filter_allowed == true)
            {
                if (coefficients_in_input)
                {
                    std::cout.precision(G_PRECISION);
                    cout << fc.coefficient << " ";
                }
                cout << fc.g.print() << endl;
            }
        }
        else
        {
            if (filter_allowed == false)
            {
                if (coefficients_in_input)
                {
                    std::cout.precision(G_PRECISION);
                    cout << fc.coefficient << " ";
                }
                cout << fc.g.print() << endl;
            }
        }
        if (coefficients_in_input) (*istr) >> fc.coefficient;
    }    
}


void filter_flags_using_subflags_fun(const string& pathall, const string& pathsubflagfilter, bool filter_forbidden_subflag, bool coefficients_in_input)
{
    vector<flag> filter_subflags;
    load_flags_from_file(pathsubflagfilter, filter_subflags);
    cerr << "Loaded " << filter_subflags.size() << " flags that" << endl;
    if (filter_forbidden_subflag)
    {
        cerr << " are forbidden as subflags" << endl;
    }
    else
    {
        cerr << " must at least once occur as subflags" << endl;
    }
    
    OPEN_FILE_SMARTLY(istr, pathall);

    flag_and_coefficient fc;
    
    if (coefficients_in_input) (*istr) >> fc.coefficient;
    while (fc.g.load_from_stream((*istr),-1,-1))
    {
        bool contains_subflag = false;
        for (int i = 0; i < (int)filter_subflags.size(); i++)
        {
            if (fc.g.contains_as_subflag(filter_subflags[i]))
            {
                contains_subflag = true;
                break;
            }
        }
        
        if (contains_subflag == false && filter_forbidden_subflag == true)
        {
                if (coefficients_in_input)
                {
                    std::cout.precision(G_PRECISION);
                    cout << fc.coefficient << " ";
                }
                cout << fc.g.print() << endl;
        }
        
        if (coefficients_in_input) (*istr) >> fc.coefficient;
    }
    
    //infile.close();
}


inline bool EXISTS_F_IN_blowup_of_H(const flag &F1, const flag &H)
{
    
    
    // Take mappings of vertices from F1 to H .... Then count how many are isomorphic to the blow-up.
    int mapping[V];
    for (int i = 0; i < V; i++)
    {
        mapping[i] = 0;
    }
    
    while (mapping[F1.m_vertices] == 0)
    {
        // Test current mapping
        flag tmp;
        tmp.as_subflag_in_blowup(H,mapping,F1.m_vertices, F1.m_Theta);
        
        if (tmp.is_isomorphic_to(F1)) 
        {
            return true;
        }
        
        // Increment mapping
        // Keep ordering of the array such that mapping[0] >= mapping[1] >= mapping[2] >=....
        bool add_one = true;
        int  add_to = 0;  // Coordinate in the array that overflew...
        while (add_one)
        {
            if (++mapping[add_to] >= H.m_vertices)
            {
                    add_to++;                    
            }
            else
            {
                while (--add_to >= 0)
                {
                    mapping[add_to] = mapping[add_to+1];
                }
                add_one = false;
            }
        }
    }
    
    return false;
}


inline double P_F1_IN_blowup_of_H(const flag &F1, const flag &H, bool density=true)
{
    
    
    // Take mappings of vertices from F1 to H .... Then count how many are isomorphic to the blow-up.
    int mapping[V];
    for (int i = 0; i < V; i++)
    {
        mapping[i] = 0;
    }
    
    int64_t good_maps = 0;
    int64_t all_maps = 0;
    
    
    // Compute all maps as H.m_vertices ^ F1.m_vertices
    all_maps = 1;
    for (int p = 0; p < F1.m_vertices; p++) all_maps *= H.m_vertices;
    
    // We are doing this because it can run in parallel

    //cout << good_maps << " " << all_maps << endl;

    //cout << "All maps=" << all_maps << endl;
    
    all_maps = 0;
    good_maps = 0;
    while (mapping[F1.m_vertices] == 0)
    {
        //if (all_maps %100000000 == 0) cout << all_maps << endl;
        //cerr << all_maps << endl;
        
        // Test current mapping
        flag tmp;
        tmp.as_subflag_in_blowup(H,mapping,F1.m_vertices, F1.m_Theta);

        all_maps++;
        if (tmp.is_isomorphic_to(F1)) 
        {
            good_maps++;
        }
        
                
        // Increment mapping
        bool add_one = true;
        int  add_to = 0;
        while (add_one)
        {
            if (++mapping[add_to] >= H.m_vertices)
            {
                mapping[add_to] = 0;
                add_to++;
            }
            else
            {
                add_one = false;
            }
        }
    }

    //cout << good_maps << " " << all_maps << endl;

//   return good_maps;
//    cout << all_maps << endl;
    
    if (density)
        return (double)good_maps/(double)all_maps;
    else
        return good_maps;
}


void print_F1_times_F2(double coeff, const flag &F1, const flag &F2)
{
    flag type;
    F1.get_type_subflag(type);
    vector<flag_and_coefficient> F1F2 = F1_times_F2(F1, F2, type);

    for (int i = 0; i < (int)F1F2.size(); i++)
    {
        cout.precision(G_PRECISION);
        cout << coeff*F1F2[i].coefficient << " " << F1F2[i].g.print() << endl;
    }
}

    


    
    
    // Builds extremal vectors for type i and in blowup of graph B
    // Extremal vector is a vector with coordinates being (labeled) flags f of type i
    //  Now the entry in the vector is p(f,B)
    //  We try to compute all possible rootings of the type and for the rooting of the type
    //  we compute the rest of the vertices. Each rooting could give us a different flag.
    //
    // Output is a vector, that contains for each flag of a given type 'coefficient' The coefficient
    // is a sum (another vector) of products (another vector).
    // The products are indexed by the numeber of vertices in B and it is saying how many verties were mapped to each
    // vertex in b. 
    // To make a simpler interpretation of this extremal vector, see the following two functions below 
vector<vector<vector<int> > > find_extremal_vectors_in_blow_up_of_B_for_vector_of_flags_with_fixed_theta(const flag &B, vector<flag> f_for_test, int *mapping_Theta, bool print_output = true)
    {
        // We have one mapping of Theta. Now count how many ways it extends to each flag.
        
        // Entry is 
        vector<vector<vector<int> > > one_extremal_vector;

        // Go over rooting of the type part of the mapping
        if (f_for_test.size() == 0) 
        {   
//            assert(0);
            return one_extremal_vector;
        }
        
        int Theta    = f_for_test[0].m_Theta;
                           
        vector<vector<int> > one_entry;
        one_extremal_vector.resize(f_for_test.size(), one_entry);
        
        
            for (int j = 0; j < (int)f_for_test.size(); j++)
            {
                flag f;
                f = f_for_test[j];
                
                //int good_maps = 0;
                vector<vector<int> > good_maps;
                
                string good_maps_str = "0";
                
                int mapping[V+1];
                for (int i = 0; i < V+1; i++)
                {
                    mapping[i] = mapping_Theta[i];
                }
                
                while (mapping[f.m_vertices] == 0)
                {
                    // Test current mapping
                    flag tmp;
                    tmp.as_subflag_in_blowup(B,mapping,f.m_vertices, f.m_Theta);
                    if (tmp.is_isomorphic_to(f)) 
                    {
                        vector<int> mapping_list;
                        mapping_list.resize(B.m_vertices, 0);
                        

                        for (int i = f.m_Theta; i < f.m_vertices; i++)
                        {
                            mapping_list[mapping[i]]++;
                        }
                    
                        good_maps.push_back(mapping_list);
                    }
                    
                    
                    // Increment mapping
                    bool add_one = true;
                    int  add_to = Theta; // Change mapping only after theta
                    while (add_one)
                    {
                        if (++mapping[add_to] >= B.m_vertices)
                        {
                            mapping[add_to] = 0;
                            add_to++;
                        }
                        else
                        {
                            add_one = false;
                        }
                    }
                }
                one_extremal_vector[j] = good_maps;
                //cout << good_maps << " ";
            }
            //cout << endl;
            
            // Add new extremal vector if found...
        if (print_output)
        {
            cout << "[";
            cout << endl;
            for (int k = 0; k < (int)one_extremal_vector.size(); k++)
            {
                if (k != 0) cout << ",";
                for (int x = 0; x < (int)one_extremal_vector[k].size(); x++)
                {
                    cout << "+";
                    for (int y = 0; y < (int)one_extremal_vector[k][x].size(); y++)
                    {
                        cout << one_extremal_vector[k][x][y];
                    }                    
                }
                
            }
            cout << "]";
        }
        
        return one_extremal_vector;
    }



vector<int> big_extremal_vector_to_int_weighted(const vector<vector<vector<int> > > &big_vector, vector<int> &weights)
{
    vector<int> one_extremal_vector;
    
        for (int j = 0; j < (int)big_vector.size(); j++)
        {
            int weight_sum = 0;
            for (int x = 0; x < (int)big_vector[j].size(); x++)
            {
                
                int one_weight=1;
                for (int y = 0; y < (int)big_vector[j][x].size(); y++)
                {
                    for (int z = 0; z < big_vector[j][x][y]; z++)
                    {
                        one_weight *= weights[y];
                    }
                }
                weight_sum += one_weight;
            }
            
            one_extremal_vector.push_back(weight_sum);
        }
    return  one_extremal_vector;   
}


vector<string> big_extremal_vector_to_str_weighted(const vector<vector<vector<int> > > &big_vector, vector<string> &weights_str)
{
    vector<string> one_extremal_vector;
    
    for (int j = 0; j < (int)big_vector.size(); j++)
    {
        string weight_sum = "";
        for (int x = 0; x < (int)big_vector[j].size(); x++)
        {
            
            string one_weight="";
            for (int y = 0; y < (int)big_vector[j][x].size(); y++)
            {
                for (int z = 0; z < big_vector[j][x][y]; z++)
                {
                    if (one_weight.length() != 0) one_weight.append("*");
                    one_weight.append(weights_str[y]);
                }
            }
            if (one_weight.length() != 0)
            {
                if (weight_sum.length() != 0) weight_sum.append("+");
                weight_sum.append(one_weight);
            }
        }
        
        if (weight_sum.length() == 0) weight_sum = "0";
        one_extremal_vector.push_back(weight_sum);
    }
    return  one_extremal_vector;   
}


    // Builds extremal vectors for type i and in blowup of graph B
    // Extremal vector is a vector with coordinates being (labeled) flags f of type i
    //  Now the entry in the vector is p(f,B)
    //  We try to compute all possible rootings of the type and for the rooting of the type
    //  we compute the rest of the vertices. Each rooting could give us a different flag.
    //
    //
    //
    //
    void find_extremal_vectors_in_blow_up_of_B_for_vector_of_flags(const flag &B, vector<flag> f_for_test, vector<int> &weights, vector<string> &weights_str, bool python_output, int verbose_output, bool use_weights_str)
    {
        // Go over rooting of the type part of the mapping
        if (f_for_test.size() == 0) return;
        int Theta    = f_for_test[0].m_Theta;
        //    int vertices = f_for_test[0].m_vertices;
        
        // First create mapping of Theta
        int mapping_Theta[V+1];
        for (int j = 0; j < V+1; j++)
        {
            mapping_Theta[j] = 0;
        }
        
        vector<vector<vector<vector<int> > > > extremal_vectors;
        //vector< vector<int> > extremal_zero(f_for_test.size(),0);

        vector<vector<int> >  extremal_vectors_int;
        vector<vector<string> >  extremal_vectors_str;
        //vector< vector<int> > extremal_zero(f_for_test.size(),0);
        //vector< vector<int> > extremal_zero_str(f_for_test.size(),"");
        
        while (mapping_Theta[Theta] == 0)
        {
            //cout << "One rooting does: ";
            // We have one mapping of Theta. Now count how many ways it extends to each flag.
//            vector<vector<vector<int> > > one_extremal_vector;
            vector<vector<vector<int> > > one_extremal_vector;
            //one_extremal_vector.resize(f_for_test.size(),0);
            
//            one_extremal_vector = find_extremal_vectors_in_blow_up_of_B_for_vector_of_flags_with_fixed_theta(B, f_for_test, weights, weights_str, mapping_Theta, false);
            one_extremal_vector = find_extremal_vectors_in_blow_up_of_B_for_vector_of_flags_with_fixed_theta(B, f_for_test,  mapping_Theta, false);

            
            // Add new extremal vector if found...
            bool vector_empty = true;
            for (int j = 0; j < (int)one_extremal_vector.size(); j++)
            {
                if (one_extremal_vector[j].size() > 0)
                {
                    vector_empty = false;
                    break;
                }
            }
            
            if (vector_empty == false)
            {                
                // Integer weights
                if (use_weights_str == false)
                {
                    vector<int> one_extremal_vector_int =  big_extremal_vector_to_int_weighted(one_extremal_vector, weights);

                    // Test if one_extremal_vector is already among discovered vectors
                    bool match=false;
                    for (int k = 0; k < (int)extremal_vectors_int.size(); k++)
                    {
                        if (extremal_vectors_int[k] == one_extremal_vector_int)
                        {
                            match = true;
                            break;
                        }
                    }
                    if (!match) 
                    {
                        if (python_output)
                        {
                            if (!extremal_vectors_int.empty()) cout << ",";
                            cout << "[";
                            for (int j = 0; j < (int)one_extremal_vector_int.size(); j++)
                            {
                                    if (j != 0) cout << ",";                            
                                    cout << one_extremal_vector_int[j];
                            }
                            cout << "]";                        
                        }
                        else
                        {
                            assert(0);
                        }
                        extremal_vectors_int.push_back(one_extremal_vector_int);
                    }                   
                }
                else
                {
                    vector<string> one_extremal_vector_str =  big_extremal_vector_to_str_weighted(one_extremal_vector, weights_str);
                    
                    // Test if one_extremal_vector is already among discovered vectors
                    bool match=false;
                    for (int k = 0; k < (int)extremal_vectors_str.size(); k++)
                    {
                        if (extremal_vectors_str[k] == one_extremal_vector_str)
                        {
                            match = true;
                            break;
                        }
                    }
                    if (!match) 
                    {
                        if (python_output)
                        {
                            if (!extremal_vectors_str.empty()) cout << ",";
                            cout << "[";
                            for (int j = 0; j < (int)one_extremal_vector_str.size(); j++)
                            {
                                if (j != 0) cout << ",";                            
                                cout << one_extremal_vector_str[j];
                            }
                            cout << "]";                        
                        }
                        else
                        {
                            assert(0);
                        }
                        extremal_vectors_str.push_back(one_extremal_vector_str);
                    }                                       
                }
                
            }
            
            
            // Increment mapping
            bool add_one = true;
            int  add_to = 0;
            while (add_one)
            {
                if (++mapping_Theta[add_to] >= B.m_vertices)
                {
                    mapping_Theta[add_to] = 0;
                    add_to++;
                }
                else
                {
                    add_one = false;
                }
            }
        }
    }
    
    


vector<int> find_extremal_vectors_read_weights(const flag &B, const string &weights_str)
{
    vector<int> weights;

    if (weights_str == "")
    {
        for (int i = 0; i < B.m_vertices; i++)
        weights.push_back(1);
    }
    else
    {
        stringstream ss(weights_str);
        for (int i = 0; i < B.m_vertices; i++)
        {
            string tmp;
            ss >> tmp;
            int w = stol(tmp);
            if (w == 0)
            {
                cerr << "Unable to convert '" << weights_str << "' in " <<  B.m_vertices << " integer weights." << endl;
                cerr << "Failed at index " << i << endl;
                assert(0);
                return weights;
            }
            weights.push_back(w);
        }
    }
    return weights;
}


vector<string> find_extremal_vectors_read_weights_str(const flag &B, const string &weights_str)
{
    vector<string> weights;
    
    if (weights_str == "")
    {
        for (int i = 0; i < B.m_vertices; i++)
            weights.push_back("");
        //cerr << "Empty" << endl;
    }
    else
    {
        stringstream ss(weights_str);
        for (int i = 0; i < B.m_vertices; i++)
        {
            string tmp;
            ss >> tmp;
            weights.push_back(tmp);
            //cout << tmp << endl;
        }
    }
    return weights;
}


void find_extremal_vectors_process_colors(const flag &B, const string &colors_str)
{
    if (colors_str == "")
    {
        //Just do nothing...
    }
    else
    {
        stringstream ss(colors_str);
        for (int i = 0; i < B.m_vertices; i++)
        {
            string tmp;
            ss >> tmp;
            int c = stol(tmp);
            if (c == 0)
            {
                cerr << "Unable to convert '" << colors_str << "' in " <<  B.m_vertices << " integer colors." << endl;
                cerr << "Failed at index " << i << endl;
                assert(0);
                return;
            }
            g_blow_up_color_edges[i] = c; 
        }
    }
}
    
void find_extremal_vectors_in_blow_up_of(const flag &B, const string &weights_in_str, const string &weights_str_in_str, const string &colors_str, bool python_output, int verbose_output)
{
    vector<int> weights = find_extremal_vectors_read_weights(B, weights_in_str);
    vector<string> weights_str = find_extremal_vectors_read_weights_str(B, weights_str_in_str);
    find_extremal_vectors_process_colors(B, colors_str);

    if (python_output) cout << "extremalVectors=[";
    for (int i = 0; i < (int)g_flags.size(); i++)
    {
        if (verbose_output)
        {
            flag type;
            g_flags[i][0].get_type_subflag(type);
            cerr << endl << "Type " << i << " i.e. " << type.print() << endl;
        }
        if (python_output)
        {
            if (i == 0)
                cout << "[";
            else
                cout << ",[";
        }
        //find_extremal_vectors_in_blow_up_of_B_for_type(B, i, weights, verbose_output);
        find_extremal_vectors_in_blow_up_of_B_for_vector_of_flags(B, g_flags[i], weights, weights_str, python_output, verbose_output, weights_str_in_str != "");
        if (python_output) cout << "]";
    }
    if (python_output) cout << "]" << endl;
}

    

void find_projection_vectors_for_one_file(const string &find_projection_vectors_input_file)
{
    vector<flag_and_coefficient> F;
    load_flags_and_coefficients_from_file(find_projection_vectors_input_file, F);

    if (F.size() == 0)
    {
        cerr << "Input file had no flags!!" << endl;
        //return;
    }

    cerr << "Loaded linear combination(s) of " << F.size() << " flags" << endl;


    // Get all possible types
    vector<flag> F_types; 

    for (const auto &fc : F)
    {
        bool known_type = false;
        for (const auto &t : F_types)
        {
            if (fc.g.have_same_type(t))
            {
                known_type = true;
                break;
            }
        }

        if (known_type == false)
        {
            flag fc_type;
            fc.g.get_type_subflag(fc_type);
            F_types.push_back(std::move(fc_type));
        }
    }

    cerr << "Identified " << F_types.size() << " different linear combinations" << endl;

    // TODO: Projection matrices - better output.
    //cerr << g_types.size() << endl;
    cout << "projection = [[] for i in range(" << g_flags.size() << ")]" << endl;
    for (const auto &F_type : F_types)
    {
        cerr << "# Using type " << F_type.print() << endl;

        //for (const auto & fl : g_flags)
        for (int i = 0; i< (int)g_flags.size(); i++)
        {
            const vector<flag> &fl = g_flags[i];
            if (fl.size() == 0) continue;

            if (fl[0].have_same_type(F_type))
            {
                cout << "projection["<<i<<"].append([";
                for (const auto &f: fl)
                {
                    double coefficient = 0;
                    for (const auto &fc : F)
                    {
                        if (!fc.g.have_same_type(F_type))
                            continue;

                        if (fc.g.is_isomorphic_to(f))
                        {
                            coefficient += fc.coefficient;
                        }
                    }
                    cout.precision(G_PRECISION);
                    cout << coefficient << ",";
                }
                cout << "])" << endl;
                break;
            }
        }
    }
}


         
void test_tight_for_blow_up_B_func(const flag &B, const string &filename, const string &weights_str, const string &colors_str, const string &mapping_Theta_str, int verbose_output)
{
    vector<flag_and_coefficient> flag_list;
    if (!load_flags_and_coefficients_from_file(filename, flag_list))
    {
        return;
    }
    
    vector<int> weights = find_extremal_vectors_read_weights(B, weights_str);
    find_extremal_vectors_process_colors(B, colors_str);    

    vector<flag> f_for_test;
    for (int i = 0; i < (int)flag_list.size(); i++)
        f_for_test.push_back(flag_list[i].g);

    int Theta = f_for_test[0].m_Theta;

    int mapping_Theta[V+1];
    for (int j = 0; j < V+1; j++)
    {
        mapping_Theta[j] = 0;
    }
    
    // Special use, where Theta is provided
    if (mapping_Theta_str != "")
    {
        stringstream ss(mapping_Theta_str);
        for (int i = 0; i < Theta; i++)
        {
            string tmp;
            ss >> tmp;
            int T = 0;
            try {
                T = stol(tmp);
            } catch (...)
            {
                cerr << "Unable to convert '" << mapping_Theta_str << "' in " <<  Theta << " integers." << endl;
                cerr << "Failed at index " << i << endl;
                assert(0);
            }
            if (T < 0 || T >= B.m_vertices)
            {
                cerr << "At index " << i << " got matting of theta to " << T << " which is not in the blowup graph." << endl;
                assert(0);
            }
            mapping_Theta[i] = T;
        }
        
        vector<int> coefs  = big_extremal_vector_to_int_weighted(find_extremal_vectors_in_blow_up_of_B_for_vector_of_flags_with_fixed_theta(B, f_for_test, mapping_Theta, false), weights);
        
        // and now make the dot product...
        double all_sum = 0;
        bool all_zeros = true;
        for (int i = 0; i < (int)flag_list.size(); i++)
        { 
            if (coefs[i] != 0) all_zeros = false;
            //cerr << "Adding " << coefs[i] * flag_list[i].coefficient << endl;
            all_sum += coefs[i] * flag_list[i].coefficient;
        }
        
        if (all_zeros)
            cerr << "Tight since always zero coefficients" << endl;
        else
        {
        if (all_sum == 0)
            cerr << "Tight since dot product zero" << endl;
        else
            cerr << "Not tight, dot product is " << all_sum << endl;
        }
        return;
    }

    
    bool all_zeros_for_all_Thetas = true;
    int only_zeros = 0;
    int zero_sums = 0;
    int nonzero_sums = 0;
    
    // Normal use that tests all Thetas
    while (mapping_Theta[Theta] == 0)
    {
        vector<int> coefs  = big_extremal_vector_to_int_weighted(find_extremal_vectors_in_blow_up_of_B_for_vector_of_flags_with_fixed_theta(B, f_for_test, mapping_Theta, false), weights);

        // and now make the dot product...
        double all_sum = 0;
        bool all_zeros = true;
        for (int i = 0; i < (int)flag_list.size(); i++)
        { 
            if (coefs[i] != 0) all_zeros = false;
            //cerr << "Adding " << coefs[i] * flag_list[i].coefficient << endl;
            all_sum += coefs[i] * flag_list[i].coefficient;
        }
        
        if (all_zeros)
        {
            only_zeros++;
        }
        else
        {
            all_zeros_for_all_Thetas = false;
            if (all_sum == 0)
                zero_sums++;
            else
            {
                nonzero_sums++;
                if (verbose_output)
                {
                    cerr << "Theta mapping ";
                    for (int i = 0; i < Theta; i++) cerr << " " << mapping_Theta[i];
                    cerr << " gave dot product " << all_sum << endl;
                }
            }
        }
        
        
        // Increment mapping
        bool add_one = true;
        int  add_to = 0;
        while (add_one)
        {
            if (++mapping_Theta[add_to] >= B.m_vertices)
            {
                mapping_Theta[add_to] = 0;
                add_to++;
            }
            else
            {
                add_one = false;
            }
        }
    }    

    if (all_zeros_for_all_Thetas)
    {
        cerr << "All mappings had all coefficients always zero" << endl;
    }
    else
    {
        cerr << "Mappings with only zeros          : " << only_zeros << endl;
        cerr << "Mappings with     zero dot product: " << zero_sums << endl;
        cerr << "Mappings with non-zero dot product: " << nonzero_sums << endl;
    }
}
    
    
int labeled_automorphisms(const flag &G)
{
    flag Gl = G;
    Gl.make_all_vertices_labeled();
    
    int count=0;
    
    int vertex_permutation[] = {0,1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20}; // probably longer than needed
     
    flag F;
    do {
        F.as_subflag(G,vertex_permutation, G.m_vertices, G.m_vertices);
        
        
        if (F.is_isomorphic_to(Gl))
        {
            count++;
        }
         
     } while ( std::next_permutation(vertex_permutation,vertex_permutation+G.m_vertices));
    
    
    return count;
}


class fraction
{
    public:
    fraction()
    {
        m_numerator = 0;
        m_denominator = 1;
    }
    
    fraction(int i)
    {
        m_numerator = i;
        m_denominator = 1;
    }
    
    fraction(int num, int den)
    {
        m_numerator = num;
        m_denominator = den;
        simplify();
    }

    
    
    double get_double()
    {
        return (double)m_numerator/(double)m_denominator;
    }

    void simplify()
    {
        if (m_numerator == 0)
        {
            m_denominator = 1;
            return;
        }
        
        for (int t = 2; t < m_denominator; t++)
        {
            if (m_numerator%t == 0  && m_denominator%t == 0)
            {
                m_numerator /= t;
                m_denominator /= t;
                t--;
            }
        }
    }
    
    fraction operator*(const fraction& a)
    {
        fraction f;
        f.m_numerator = m_numerator * a.m_numerator;
        f.m_denominator = m_denominator * a.m_denominator;
        f.simplify();
        return f;
    }

    fraction operator+(const fraction& a)
    {
        fraction f;
        f.m_numerator = m_numerator * a.m_denominator + m_denominator*a.m_numerator;
        f.m_denominator = m_denominator * a.m_denominator;
        f.simplify();
        return f;
    }

    fraction operator-(const fraction& a)
    {
        fraction f;
        f.m_numerator = m_numerator * a.m_denominator - m_denominator*a.m_numerator;
        f.m_denominator = m_denominator * a.m_denominator;
        f.simplify();
        return f;
    }

    fraction operator/(const long l)
    {
        fraction f;
        f.m_numerator = m_numerator;
        f.m_denominator = m_denominator * l;
        f.simplify();
        return f;
    }

    fraction operator/(const fraction& a)
    {
        fraction f;
        f.m_numerator = m_numerator * a.m_denominator;
        f.m_denominator = m_denominator * a.m_numerator;
        f.simplify();
        return f;
    }
    
    long long m_numerator;
    long long m_denominator;
};

std::ostream& operator<< (std::ostream& stream, const fraction& f)
{
    stream << f.m_numerator << "/" << f.m_denominator;
    return stream;
}


void jan_to_nice_convert(const string&  filename, bool with_coefficients)
{
    cerr << "Converting flags in Jan Volec format in file " << filename << endl;

#if G_COLORED_EDGES == 2
    
    OPEN_FILE_SMARTLY(istr, filename);

    
    /*
    istream *istr = &cin;
    
    ifstream infile;
    if (filename != "cin")
    {
        infile.open (filename.c_str(), ifstream::in);
        if (!infile.good())
        {
            cerr << "Failed opening file " << filename << endl;
            return;
        }
        istr = &infile;
    }
    */
    
    // In case loabeled vertices are out of order....
    //int remapping_jan[]={2,3,0,1,4,5,6,7,8,9,10};
    int remapping_jan[]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    
    while (istr->good())
    {
        flag f;
        int vertices;
        int theta;
        int edges;
        double coefficient = 0;
        string coefficient_str;
        //TODO
        (*istr) >> vertices >> theta >> edges;
        
        if (!istr->good()) break;
        //cerr << vertices << " " << theta << " " << edges << " " << coefficient_str << " ";
  
        
        if (with_coefficients)
        {
            (*istr) >> coefficient_str;
            if (coefficient_str[0] != '[') break;
#ifdef DONT_USE_C11
            coefficient_str.erase( coefficient_str.size() - 1,1);
#else            
            coefficient_str.pop_back();
#endif
            coefficient_str.erase(0,1);
            coefficient = stod(coefficient_str);
        }
         
        f.set_vertices_and_Theta(vertices,theta);
        
        for (int u = 0; u < vertices; u++)
            for (int v = u+1; v < vertices; v++)
                f.color_edge(u,v,1);
        
        for (int e = 0; e < edges; e++)
        {
            string e_str;
            (*istr) >> e_str;
            int u = (int)e_str[0] - 'a';
            int v = (int)e_str[1] - 'a';
            
            u = remapping_jan[u];
            v = remapping_jan[v];
            
            
            f.color_edge(u,v,2);
        }
        
        if (with_coefficients) 
        {
            cout.precision(G_PRECISION);
            cout << coefficient << " ";
        }
        cout << f.print() << endl;
    }
    //jan.close();
#else
    cerr << "Not implemented." << endl;
#endif
}

void find_density_vectors_in_iterated_blow_up_of(const flag &B, int Kn)
{
    cerr << "Finding density vectors in iterated blow-up" << endl;
    

#ifdef G_COLORED_EDGES_BLIND
    cerr << "Not sure how to do iterated constructions with blind colored edges." << endl;
    return;
#endif
    

#ifdef G_COLORED_EDGES
    
    vector<double> densities[V];
    vector<double> densities_labeled[V];

    vector<fraction> f_densities[V];
    vector<fraction> f_densities_labeled[V];

    
    // We say that vertex always has density one
    // which makes life nice
    densities[1].push_back(1);
    densities_labeled[1].push_back(1);
    
    f_densities[1].push_back(1);
    f_densities_labeled[1].push_back(1);
    
    int VB = B.m_vertices;
    
    for (int VG = 2; VG <= Kn; VG++)
    {
        cerr << "****************** VG=" << VG << endl;
        double sum_density=0;
        
        for (int v = 0; v < (int)g_unlabeled_flags[VG].size(); v++)
        {
            
            flag G = g_unlabeled_flags[VG][v];
            //cerr << "Testing " << G.print() << endl;
            
            // here we generate all mappings of i
            int mapping[VG];
            for (int t = 0; t < VG; t++) mapping[t] = 0;
            
            double label_d = 0; // density counted so far - labeled
            fraction f_label_d = 0;
            
            int recursive_d_cnt = 0;
            
            // This is unrolled recursion for generating the mapping
            while (true)
            {
                //int all_maps = 0;
                
                // test mapping
                // How many mapped to each class
                vector<int> mapped_to[V];
                for (int t = 0; t < VG; t++)
                    mapped_to[mapping[t]].push_back(t); 


                // check that edges between pieces in different parts are correct
                bool invalid_map = false;
                bool recursive_map = false;
                
                for (int a = 0; a < VG-1; a++)
                    for (int b = a+1; b < VG; b++)
                    {
                        if (mapping[a] == mapping[b])
                            continue;
                        if (G.m_color_edge[a][b] != B.m_color_edge[mapping[a]][mapping[b]])
                        {
                            invalid_map = true;
                        }
                        
                    }
                
                if (!invalid_map)
                {
                    double density_map = 1;
                    fraction f_density_map = 1;
                    
                                
                    // count probability inside each part separately
                    for (int n = 0; n < VB; n++)
                    {
                        // all recursively mapped to the same blob
                        if (mapped_to[n].size() == (unsigned int)VG)
                        {
                            recursive_d_cnt++;
                            recursive_map = true;
                            break;
                        }
                    
                        // If not recursive, count what is the probability of the induced subrgaph in n
                        //TODO some scale will be needed
                        if (mapped_to[n].size() >= 2)
                        {
                            int mapping_G_in_n[V];
                            for (int j = 0; j < (int)mapped_to[n].size(); j++) mapping_G_in_n[j] = mapped_to[n][j];
                            flag G_in_n;
                            G_in_n.as_subflag(G, mapping_G_in_n, mapped_to[n].size(), 0);
                            
                            int small_ID = find_flag_in_list(G_in_n, g_unlabeled_flags[mapped_to[n].size()]);
                            assert(small_ID >= 0);
                        
                            density_map *= densities_labeled[mapped_to[n].size()][small_ID];
                            f_density_map = f_density_map * f_densities_labeled[mapped_to[n].size()][small_ID];
                        }
                    }
                    
                    if (!recursive_map)
                    {
                        label_d += density_map;
                        f_label_d = f_label_d + f_density_map;
                    }

                    //cerr << "_label_d="  << label_d << endl;
                    //cerr << "_f_label_d="  << f_label_d << endl;
                    
                    
                }
                            
                
                /*
                cerr << "Map ";
                for (int t = 0; t < i ; t++) 
                    cerr << mapping[t] << " ";
                cerr << "label_d=" << label_d  << " inv_map="  << invalid_map  <<  " rec_map=" <<  recursive_map  <<  " rec_cnt=" <<  recursive_d_cnt << endl;
                */
                
                // increment mapping
                int inc = 0;
                while (inc < VG && mapping[inc] == VB-1) mapping[inc++] = 0;
                if (inc == VG) break; // we are done
                mapping[inc]++;
            }
            
            // count number of non-isomorphic labelings on G....
            // factorial(i)/labeled_automorphisms(G)
            
            //cerr << "label_d="  << label_d << endl;
            //cerr << "f_label_d="  << f_label_d << endl;
            
            //cerr << "labeled_automorphisms(G)=" << labeled_automorphisms(G) << endl;
            double non_iso_cnt = (double)factorial(VG)/labeled_automorphisms(G);
            fraction f_non_iso_cnt = fraction(factorial(VG),labeled_automorphisms(G));
            //cerr << "non_iso_cnt=" << non_iso_cnt << endl;
            //cerr << "f_non_iso_cnt=" << f_non_iso_cnt << endl;
            
            
            // non-recursive part
            double unlabel_d = non_iso_cnt* label_d/pow((double)VB,(double)VG);
            fraction f_unlabel_d = f_non_iso_cnt*f_label_d/( (long)pow(VB,VG));
            //cerr << "non-recursive: " << unlabel_d << endl;
            //cerr << "non-recursive: " << f_unlabel_d << endl;
            
            // now recursive part
            // d \binom{n}{VG}    =  d*VB*\binom{n/VB}{VG} + C
            // d ( \frac{1}{VG!}(1- VB^(VG-1))              =  c
            //
            double recusive_scale = 1 - 1/pow((double)VB,(double)VG-1);
            fraction f_recusive_scale = fraction(1) -  fraction(1)/((long)pow((double)VB,(double)VG-1));
            //cerr << "scale: " << (recusive_scale) << endl;
            //cerr << "scale: " << (f_recusive_scale) << endl;
            double full_d = unlabel_d / recusive_scale;
            fraction f_full_d = f_unlabel_d / f_recusive_scale;
            
            //cout << std::setw(10) << full_d << "   " << f_full_d  << "     " << G.print() << endl;
            cout << std::setw(10) << full_d << "     " << G.print() << endl;
            //cout << f_full_d << " " << G.print() << endl;

            //cerr << "full_d: " << full_d << endl;
            //cerr << "label_d: " << label_d/pow((double)VB,(double)VG)/recusive_scale << endl;
            
            densities[VG].push_back(full_d);
            densities_labeled[VG].push_back(label_d/pow((double)VB,(double)VG)/recusive_scale);
            f_densities[VG].push_back(f_full_d);
            f_densities_labeled[VG].push_back(f_label_d/((long)pow((double)VB,(double)VG))/f_recusive_scale);
            sum_density += full_d;
        }
        
        cerr << "Density sum=" << sum_density << endl;
    }
    
    cerr << "WARNING: This is experimental feature and may not work 100% correctly." << endl;
#else
    cerr << "Feature not implemented for this kind of flags" << endl;
#endif
}





void extend_flag()
{
    
}

    
    
template<typename T>
void make_patterns(vector<T> &ptr,  vector<vector<T> > &patterns, int length, int maxx, T scale)
{
    if ((int)ptr.size() == length)
    {
        patterns.push_back(ptr);
        return;
    }
    for (int x = 0; x < maxx; x++)
    {
        ptr.push_back(x/scale);
        make_patterns(ptr,patterns,length,maxx,scale);
        ptr.pop_back();
    }
}


template<typename T>
void expand_pattern(vector<T> &ptr, vector<T> &template_ptr, vector<vector<T> > &patterns, int length, int maxx, T scale)
{
    if ((int)ptr.size() == length)
    {
        patterns.push_back(ptr);
        return;
    }
    
    if (template_ptr[ptr.size()] == 0 || template_ptr[ptr.size()] == 1)
    {
        ptr.push_back(template_ptr[ptr.size()]);
        expand_pattern(ptr,template_ptr,patterns,length,maxx,scale);
        ptr.pop_back();
        return;
    }
    
    for (int x = 0; x < maxx; x++)
    {
        ptr.push_back(x/scale);
        expand_pattern(ptr,template_ptr,patterns,length,maxx,scale);
        ptr.pop_back();
    }
}


void fix_part(vector<vector<int> > adjacencies, vector<double> &tmpd, int a, int b, int c, int d, int e, double part)
{
    for (int t = 0; t < (int)adjacencies.size(); t++)
    {
        if (adjacencies[t][0] == a && adjacencies[t][1] == b && adjacencies[t][2] == c && adjacencies[t][3] == d && adjacencies[t][4] == e)
        {
            tmpd[t] = part;
            return;
        }
    }
}


void fix_part(vector<vector<int> > adjacencies, vector<double> &tmpd, int a, int b, int c, int d, double part)
{
    for (int t = 0; t < (int)adjacencies.size(); t++)
    {
        if (adjacencies[t][0] == a && adjacencies[t][1] == b && adjacencies[t][2] == c && adjacencies[t][3] == d)
        {
            tmpd[t] = part;
            return;
        }
    }
}


void fix_part(vector<vector<int> > adjacencies, vector<double> &tmpd, int a, int b, int c,  double part)
{
    for (int t = 0; t < (int)adjacencies.size(); t++)
    {
        if (adjacencies[t][0] == a && adjacencies[t][1] == b && adjacencies[t][2] == c)
        {
            tmpd[t] = part;
            return;
        }
    }
}



void fix_part(vector<vector<int> > adjacencies, vector<double> &tmpd, int a, int b, double part)
{
    for (int t = 0; t < (int)adjacencies.size(); t++)
    {
        if (adjacencies[t][0] == a && adjacencies[t][1] == b)
        {
            tmpd[t] = part;
            return;
        }
    }
}


double get_part(vector<vector<int> > adjacencies, vector<double> &part, int a, int b, int c, int d)
{
    for (int t = 0; t < (int)adjacencies.size(); t++)
    {
        if (adjacencies[t][0] == a && adjacencies[t][1] == b && adjacencies[t][2] == c && adjacencies[t][3] == d)
        {
            return part[t];
        }
    }
    assert(0);
    return -1;
}

double get_part(vector<vector<int> > adjacencies, vector<double> &part, int a, int b, int c)
{
    for (int t = 0; t < (int)adjacencies.size(); t++)
    {
        if (adjacencies[t][0] == a && adjacencies[t][1] == b && adjacencies[t][2] == c)
        {
            return part[t];
        }
    }
    assert(0);
    return -1;
}

bool are_constraints_same()
{
    
    return false;
}



void generate_subcombinations(int n, int r, int shift, vector< vector<int> > &result)
{
    std::vector<bool> v(n);
    std::fill(v.begin(), v.begin() + r, true);

    do {
        vector<int>  one_result;
        for (int i = 0; i < n; ++i) {
            if (v[i]) {
                //std::cout << (i) << " ";
                one_result.push_back(i+shift);
            }
        }
        result.push_back(one_result);
        //std::cout << "\n";
    } while (std::prev_permutation(v.begin(), v.end()));
}


//void hack_cuts_generator()
void hack(string cutinput = "cutinput.txt")
{
      
} 

// Adding symmetry constraints....
void hack_cuts()
{
}

// Translates rotation systems to crossings
void hack_cross()
{
}


void extra_arguments(int extra, int processed, int argc)
{
    if (processed+extra >= argc)
    {
        cerr << "Expected more arguments..." << endl;
        exit(1);
    }
}




// THis prints extensive help
void help()
{
    cout << "Extensive help for flag.cpp program" << endl
    << "Compile flags: "
    
    
    
#ifdef G_COLORED_EDGES    
    << "G_COLORED_EDGES=" << G_COLORED_EDGES
#endif
    
    
    << endl
    << "File format:" << endl
    << "n         Number of vertices, integer >= 0" << endl
    
    << "Theta     Number of lableled vertices. First Theta vertices are\n"
    << "          considered labeled\n"

    
    
#ifdef G_COLORED_EDGES
    << "edges     Colors of edges. Listing by top right of adjacency matrix that contains color entries.\n"
    << "          For example C_5 in colors 2 and 1 would look like  2 1 1 2   1 1 2   1 2   2.\n"
    << "          Color 0 can be used instead of 'any' color. Good when asking for extensions. Use caution when\n"
    << "          actually using it in some program. It might not work right.\n"    
#endif    
    
    
 
   
    ;
    
}


void expand_linear_constraints_by_products(int Kn, int verbose_output, string dump_to_tmp_prefix="")
{
    // g_flag_product_linear_constraints[ constraint ] is a vector of flags that can be used to multiply  [constraint]
    //static vector<vector<flag> >  g_flag_product_linear_constraints;
    
    // First try to find all types
    vector<flag> types_used;
    for (int j = 0; j < (int)g_linear_constraints.size(); j++)
    {
        add_g_to_flags_list_if_new(g_linear_constraints[j].m_type, types_used);
    }
    cerr << "Using " << types_used.size() << endl;

    vector<flag> types_not_generated;
    for (const auto &f : types_used)
    {
        // try loading, if fail, add to list for generation
        if (labeled_flags_of_one_already_exist(Kn, f) == false)
        {
            types_not_generated.push_back(f);
            cerr << "Will generate  labeled flags of size " << Kn << " of type " << f.print() << endl;
        }
    }

    #pragma omp parallel for ordered schedule(dynamic)
    for (int j = 0; j < (int)types_not_generated.size(); j++)
    {
        // Generate the files here!!!!
        vector<flag> flag_list;
        generate_labeled_flags_of_one_type(Kn, types_not_generated[j], flag_list);

        string filename = get_filename_for_labeled_flags(Kn,types_not_generated[j]);
        #pragma omp critical
        {
            g_labeled_flags_of_one_type_map.insert(make_pair(filename, flag_list));
        }
    }

    vector<vector<linear_constraint> > linear_constraints_expanded;

    linear_constraints_expanded.resize(g_linear_constraints.size());

    int removed_duplicate_constraints = 0;

    mini_timer mt;

    #pragma omp parallel for ordered schedule(dynamic)
    for (int j = 0; j < (int)g_linear_constraints.size(); j++)
    {
                        
        // Size of the other flags in multiplication
        int type_size = g_linear_constraints[j].m_labeled_vertices_in_type_cnt;         
        int flag_size = (Kn - g_linear_constraints[j].m_entries_max_size) + type_size;

        if (flag_size < 2)   // multiplying by just 1 vertex does not do much (unless the vertex has a color)
        {
            cerr << "Size of big graphs is too small to use products linear constraint " << j+1 << endl;
            continue;
        }
        
                
        // getting labeled flags of the right size

        vector<flag> constraint_multiplied_by;
        flag type;
        type = g_linear_constraints[j].m_type;
        get_labeled_flags_of_one_type(flag_size, type, constraint_multiplied_by);
        
        if (constraint_multiplied_by.size() == 0)
        {
            cerr << "Linear constraint " << j+1 << " cannot be used since generating flags of type " <<  type.print() << " and size " <<  flag_size  << " gave zero flags" << endl;
            continue;
        }
        
        if (verbose_output)
        {
            #pragma omp critical(cerr)
            {
                cerr << "Constraints " << j << "/" << g_linear_constraints.size() << " can be multiplied by " <<  constraint_multiplied_by.size()  << " other flags "
                     << mt.report(j, g_linear_constraints.size()) << endl;
            }
        }

        //constraints_possible++;

//        #pragma omp parallel for ordered schedule(dynamic)
        for (int z = 0; z < (int)constraint_multiplied_by.size(); z++)
        {
            flag f = constraint_multiplied_by[z];

            //g_linear_constraints[j].m_entries;

            linear_constraint new_lc_f_i;
            //check_constraint()

            vector<flag_and_coefficient> product_labeled;

            for (int i = 0; i < (int)g_linear_constraints[j].m_entries.size(); i++)
            {
                    //cerr << i << " " << j << endl;
                    flag type = g_linear_constraints[j].m_type;
                    vector<flag_and_coefficient> F1F2 = F1_times_F2(g_linear_constraints[j].m_entries[i].g, f, type);
                                    
                    multiply_FC_by_C(F1F2, g_linear_constraints[j].m_entries[i].coefficient);
                    
                    fc_add_FC_to_first(product_labeled,F1F2);
            }

            // Unlabeling fc
            //cerr << "Unlabeling " << endl;

            for (int i = 0; i < (int)product_labeled.size(); i++)
            {
                flag_and_coefficient fc;
                
                fc.g = product_labeled[i].g;
                fc.g.m_Theta = 0;
                fc.coefficient = P_F1_IN_H(product_labeled[i].g, fc.g)*product_labeled[i].coefficient;
                
                vector<flag_and_coefficient> tmp;
                tmp.push_back(fc);
                
                fc_add_FC_to_first(new_lc_f_i.m_entries, tmp);
            }

            // Sometimes linear constraints may contain same entry several times
            // so we just sum these all up
            simplify_FC(new_lc_f_i.m_entries);

            if (!new_lc_f_i.check_constraint())
            {
                cerr << "Something went wrong!" << endl;
                exit(1);
            }



            // Each has a private one...
            //#pragma omp critical
            {
                //linear_constraints_expanded[j].push_back(new_lc_f_i);
            }
             if (!add_linear_constraint_if_new(new_lc_f_i, false,linear_constraints_expanded[j]))
             {
                #pragma omp atomic
                    removed_duplicate_constraints++;
             }
        }
    }


//
//    
//

    cerr << "Merging generated constraints" << endl;
    //vector<linear_constraint> > linear_constraints_reduced;
    int provided_constraints = g_linear_constraints.size();
    g_linear_constraints.clear();
    int generated_constraints = 0;
    for (int j = 0; j < (int)linear_constraints_expanded.size(); j++)
    {
        generated_constraints += linear_constraints_expanded[j].size();
        merge_vectors_parallel(g_linear_constraints, linear_constraints_expanded[j]);
    }

    removed_duplicate_constraints += generated_constraints - g_linear_constraints.size();

    cerr << "Got " << generated_constraints << " expanded constraints out of " << provided_constraints << " provided ones and after removing duplicates " << g_linear_constraints.size() << endl;

    remove_constraints_implied_by_others(verbose_output);


    if (removed_duplicate_constraints > 0)
    {
        cerr << "Removed " << removed_duplicate_constraints << " duplicate constriants." << endl;
    }

    cerr << "New lc " <<  g_linear_constraints.size() << endl;
}



#include <sys/time.h>
#include <sys/resource.h>
// return memory in Kbytes
long getMemoryUsage() 
{
  struct rusage usage;
  if(0 == getrusage(RUSAGE_SELF, &usage))
#ifdef __APPLE__    
    return usage.ru_maxrss/1024; // bytes on OSX, Kbytes on Linux
#else
    return usage.ru_maxrss; // bytes on OSX, Kbytes on Linux
#endif
  else
    return 0;
}


string getMemoryUsageStr()
{
    long used_K = getMemoryUsage();
    long used_M = used_K / 1024 ;
    long used_G = used_K / 1024 / 1024;

    stringstream ss;

    if (used_G > 0)
    {
        ss << used_G << "G";
    }
    else if (used_M > 0)
    {
        ss << used_M << "M";
    }
    else
    {
        ss << used_K << "K";
    }
    return ss.str();
}

string g_log_string;
bool g_log_disabled = false;
time_t g_log_start_timestamp=0;
string g_log_filename = "flag.log";
int g_log_min_time_taken = 3600; // Minimum number of seconds to run to generate a log


void write_log()
{
    if (g_log_disabled)
    {
        return;
    }

    time_t time_taken = time(NULL) - g_log_start_timestamp;
    if (g_log_start_timestamp == 0)
    {
        time_taken = 0;
    }

    char hostname[512];
    char username[512];
    gethostname(hostname, 512);
    getlogin_r(username, 512);

    // Declaring argument for time() 
    time_t tt; 
  
    // Declaring variable to store return value of 
    // localtime() 
    struct tm * ti; 
  
    // Applying time() 
    time (&tt); 
  
    // Using localtime() 
    ti = localtime(&tt); 

    string current_time_str = asctime(ti);
    current_time_str.pop_back();

    stringstream ss;
    ss << current_time_str << " " << username << "@" << hostname << " " << g_log_string << " time taken " 
          << time_to_str(time_taken) <<  " Memory=" <<  getMemoryUsageStr() <<  endl; 

    cerr << ss.str();

    // More than 1 hour of time
    if (time_taken > g_log_min_time_taken)
    {
        string logfile=getenv("HOME");
        logfile +=  "/" + g_log_filename;

        std::ofstream outfile;
        outfile.open(logfile, std::ios_base::app | std::ios::out);
        if (outfile.is_open())
        {
            outfile << ss.str(); 
            cerr << "Written to a log file " << logfile << endl;
        }
        else
        {
            cerr << "Error opening log file '" << logfile << "':" <<  strerror(errno) << endl;
        }
    }
}

void init_log(int argc, char *argv[])
{
    g_log_start_timestamp = time(NULL);

    stringstream ss;

    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) != NULL) {
        ss << cwd;
    }
    else
    {
        ss << "UNKNOWN_DIRECTORY";
    }

    for (int i  = 0; i < argc; i++)
    {
        string tmp(argv[i]);
        if (tmp.find_first_of("\t\n ") == string::npos)
            ss << " " << argv[i] << " ";
        else
            ss << " \'" << argv[i] << "\' ";
    }
    
    g_log_string = ss.str();
}

#ifdef G_EXTRNAL_HACK_CPP
#include "external_hack.cpp"
//void external_hack();
#endif

#ifndef G_DISABLE_MAIN
int main(int argc, char *argv[])
{
    init_log(argc, argv);
    atexit(write_log);

    // Parse command line
    int Kn = 0;
    string objective_file = "";
    string objective_file_divisor = "";
    string output_file = "";
    int verbose_output = 0;
    bool use_sdp_temp = false;
    bool extension = false;
    bool extension_count_copies = false;
    string extensions_str = "";
    //flag g_extensions;

    bool forbidden_test = false;
    flag g_forbidden_test;    
    bool use_sdpa = false;
    
    bool extensions_fc = false;
    string extensions_fc_file = "";

    bool dump_unlabeled = false;
    bool dump_unlabeled_while_generating = false;
    bool generate_small_unlabeled_from_large = false;
    bool quit_after_generating_unlabeled = false;
   //bool quit_after_generating_labeled = false;
    bool force_generate_flags = false;
    bool force_generate_labeled_flags = false;
    bool remove_duplicates_while_loading = false;
    bool remove_forbidden_wile_loading=false;
    bool types_from_file = false;

    bool dump_types_in_order = false;

    
    bool upper_bound = false;
    bool lower_bound = false;
    bool force_generating_constriants = false;
    
    bool don_run_csdp = false;
    bool only_compute_objective = false;
    
    bool find_extremal_vectors = false;
    flag find_extremal_vectors_construction;
    string find_extremal_vectors_weights = "";
    string find_extremal_vectors_weights_str = "";
    string find_extremal_vectors_colors = "";
    bool find_extremal_vectors_python_output = true;


    
    bool print_problem_in_latex_and_quit = false;
    bool generate_baber_equalities_and_quit = false;
    
    bool print_constraints_blocks_and_quit = false;
    int  print_constraints_blocks_and_quit_start_from=0;
    
    bool draw_graphs_in_file = false;
    string draw_graphs_in_filename = "";
    //bool draw_graphs_with_densities = false;
    bool draw_graphs_color_1_nonedge = true;
    bool draw_graphs_use_enumerate = false;
    bool draw_graphs_use_cout = false;

    bool process_solution = false;
    double process_solution_lower_bound = 0;
    double process_solution_upper_bound = 0.000001;
    bool process_solution_print_density = false;

    bool compute_denisties_in_file = false;
    string compute_denisties_in_filename = "";

    //bool no_slack_flags = false;
    //string no_slack_flags_filename = "";


    
    string csdp_binary="csdp";

    bool P_F_ISO_H_test = false;
    flag P_F_ISO_H_test_F;
    flag P_F_ISO_H_test_H;
    


    bool flag_calculator_general = false;
    bool flag_calculator_general_draw = false;    

    string flag_calculator_F1 = "";
    string flag_calculator_F2 = "";
    double flag_calculator_C = 1;
    string flag_calculator_TYPE = "";
    double flag_calculator_epsilon = 0;
    string flag_calculator_general_string = "";

    bool generate_subflags_of_size_n_switch = false;
    int generate_subflags_of_size_n_n = 0;
    string generate_subflags_of_size_n_input = "";

    vector<string> additional_constraints;
    vector<string> additional_constraints_strings;
    bool linear_constriants_remove_duplicates_with_types = false;
//    bool dump_linear_constriants_and_quite = false;
    
    string csdp_OMP = "";

    // Enables to call a python script before csdp happens
    string csdp_preprocessing = "";

    for (int i = 0; i < V; i++) g_blow_up_color_edges[i] = G_BLOW_UP_COLOR_EDGES; 

    
    bool use_hack = false;
    string use_hack_1 = "";
    
    if (argc <= 1)
    {
        cerr << endl;
        cerr << "Usage: " << argv[0] << " -n SIZE_OF_UNLABELED_FLAGS [-obj OBJECTIVE_FILENAME] [-v] [-o OUTPUT_DAT-S] [-lb | -ub]" << endl;
        cerr << " -n SIZE_OF_UNLABELED_FLAGS   Size of unlabeled flags" << endl;
        cerr << " -obj OBJECTIVE_FILENAME      File containing objective function" << endl;
        cerr << " -v                           Display progress when creating SDP and generating flags and others" << endl;
        cerr << " -ub                          Creates SDP computing UPPER BOUND on the objective" << endl;        
        cerr << " -lb                          Creates SDP computing LOWER BOUND on the objective" << endl;
        cerr << " -elcp                        Enable Linear Constraints products - multiplied with other flags." << endl;
        cerr << " -fp                          Stores flag products or tries to reuse them - usefull for repeated computations, eats more space" << endl;
        cerr << " -fc EXPRESSION               Flag Calculator: simple calculator" << endl;
        cerr << " -fcd EXPRESSION              Flag Calculator: draw simple calculator" << endl;
        cerr << endl;

   
        cerr << "Filenames:" << endl;
        cerr << "    " << filename_prefix() << "__n?_unlabeled.txt    unlabeled flags" << endl;
        cerr << "    " << filename_prefix() << "__n?_labeled.txt      labeled flags" << endl;
        cerr << "    " << filename_prefix() << "__n?_sdp_products.txt part of sdp" << endl;
        cerr << "    " << filename_prefix() << "__forbidden.txt       forbidden flags" << endl;
        cerr << "    " << filename_prefix() << "__objective.txt       default objective" << endl;
#ifdef G_COLORED_EDGES_BLIND
        cerr << "    " << filename_prefix() << "__edgeblind_permutations.txt  allowed perms" << endl;
#endif        
        cerr << endl;
        return 0;
    }  
    
    
    int i = 1;
    while (i < argc)
    {        
        if (string(argv[i]) == "-help") {
            help();
            return 0;
        } else if (string(argv[i]) == "-n") {
            extra_arguments(1,i,argc);
            Kn = atoi(argv[i+1]);
            i++;
        } else if (string(argv[i]) == "-omp" || string(argv[i]) == "-OMP") {
            extra_arguments(1,i,argc);
            int num_threads = atoi(argv[i+1]);
            if (num_threads > 0)
            {
#ifdef _USING_OMP_                
                omp_set_num_threads(num_threads);
#else
                cerr << "ERROR: Program compiled without OpenMP support. -omp is useless" << endl;
#endif                
            }
            else
            {
                cerr << "ERROR: num_threads in -omp should be an integer > 0";
                cerr << "  maybe coversion error for '" << argv[i+1] << "'" << endl;
                return 1;
            }            
            i++;
        } else if (string(argv[i]) == "-obj") {
            extra_arguments(1,i,argc);
            objective_file = argv[i + 1];
            i++;
        } else if (string(argv[i]) == "-o") {
            extra_arguments(1,i,argc);
            output_file = argv[i + 1];
            i++;
        } else if (string(argv[i]) == "-csdp") {
            extra_arguments(1,i,argc);
            csdp_binary = argv[i + 1];
            i++;
        } else if (string(argv[i]) == "-v") {
            verbose_output = 1;
        } else if (string(argv[i]) == "-vv") {
            verbose_output = 2;
        } else if (string(argv[i]) == "-vvv") {
            verbose_output = 3;
        } else if (string(argv[i]) == "-vvvv") {
            verbose_output = 4;
        } else if (string(argv[i]) == "-latex") {
            print_problem_in_latex_and_quit = true;
        } else if (string(argv[i]) == "-lb") {
            lower_bound = true;
        } else if (string(argv[i]) == "-ub") {
            upper_bound = true;
        } else if (string(argv[i]) == "-srp") {
            extra_arguments(1,i,argc);
            g_smart_round_precision = strtod(argv[i+1], NULL);
            i++;     
        } else if (string(argv[i]) == "-cdin") {
            extra_arguments(1,i,argc);
            compute_denisties_in_file = true;
            compute_denisties_in_filename = argv[i + 1];
            i++;
        } else if (string(argv[i]) == "-ps") {
            process_solution = true;
        } else if (string(argv[i]) == "-psd") {
            process_solution = true;
            process_solution_print_density = true;
        } else if (string(argv[i]) == "-gbe") {
            generate_baber_equalities_and_quit = true;
        } else if (string(argv[i]) == "-ncsdp") {
            don_run_csdp = true; 
        } else if (string(argv[i]) == "-nocsdp") {
            don_run_csdp = true; 
        } else if (string(argv[i]) == "-csdp_OMP") {
            extra_arguments(1,i,argc);
            csdp_OMP = argv[i + 1];
            i++;
        } else if (string(argv[i]) == "-csdp_pp") {
            extra_arguments(1,i,argc);
            csdp_preprocessing = argv[i + 1];
            i++;
        } else if (string(argv[i]) == "-sdpa") {
            use_sdpa = true;     
        } else if (string(argv[i]) == "-fgf") {
            force_generate_flags = true;
            force_generate_labeled_flags = true;

        } else if (string(argv[i]) == "-tff") {
            types_from_file = true;
        } else if (string(argv[i]) == "-fglf") {
            force_generate_labeled_flags = true;
        } else if (string(argv[i]) == "-fev") {
            extra_arguments(1,i,argc);
            find_extremal_vectors = true;
            find_extremal_vectors_construction.load_from_string(argv[i + 1]);
            i++;
        } else if (string(argv[i]) == "-fp") {
            use_sdp_temp = true;            
        } else if (string(argv[i]) == "-dlcs") {
            g_use_simple_linear_constraints = false;
        } else if (string(argv[i]) == "-lcsuse") {
            extra_arguments(1,i,argc);
            int uselc = (int)strtol(argv[i+1],NULL,0);
            g_use_simple_linear_constraints_list.push_back(uselc);
            i++;
        } else if (string(argv[i]) == "-elcp") {
            g_use_product_linear_constraints = true;
        } else if (string(argv[i]) == "-fc") {
            extra_arguments(1,i,argc);
            flag_calculator_general = true;
            flag_calculator_general_string = argv[i + 1];
            i += 1;  
        } else if (string(argv[i]) == "-fcd") {
            extra_arguments(1,i,argc);
            flag_calculator_general = true;
            flag_calculator_general_string = argv[i + 1];
            i += 1;  
            flag_calculator_general_draw = true;         
        } else {
            cerr << "Unsupported argument " << argv[i] << " " << endl;
            return 0;
        }
        i++;
    }
    
    if (i > argc)
    {
        cerr << "Probably some missing argument." << endl;
        return 1;
    }

    
    
    
#ifdef G_LOAD_COLORED_EDGES_BLIND_PERMUTATIONS    
    load_blind_colorededge_permutations();
#endif

        
    
    load_forbidden();


    if (flag_calculator_general)
    {
        vector<flag_and_coefficient> F;


        ostream *ostr = NULL;

        if (flag_calculator_general_draw)
        {
            ostr = &cout;
        }

        if (ostr != NULL)
        {
            print_latex_header((*ostr), draw_graphs_color_1_nonedge);
        }

        flag_calculator_simple(F, flag_calculator_general_string, ostr);
        
        if (ostr != NULL)
        {
            (*ostr) << "\\end{document}" << endl;
        }


        //simplify_FC(F);  
        dump_flags_and_coefficients(F);
        return 0;
    }

    
    //////////////////////////////// Normal use....
    
    if (objective_file == "")
    {
        stringstream filename;
        filename <<  filename_prefix() << "__objective.txt";
    
        objective_file = filename.str();
        cerr << "Objective file not provided, using default filename " << objective_file << endl;
    }
    
    if (!load_objective_from_file(objective_file, linear_constriants_remove_duplicates_with_types, verbose_output))
    {
        return 1;
    }    
//#ifdef USE_FOR_CROSSINGS
//    cerr << "Warning: using special function for crossings - not everything works." << endl;
//    return main_crossings(Kn);
//#endif

    
    assert(Kn <= V);
    
    if (print_problem_in_latex_and_quit)
    {
        print_problem_in_latex(objective_file, draw_graphs_color_1_nonedge);
        return 0;
    }
    

    get_unlabeled_flags_of_size(Kn, force_generate_flags, remove_duplicates_while_loading, remove_forbidden_wile_loading, verbose_output, dump_unlabeled_while_generating);
    if (generate_small_unlabeled_from_large) generate_all_unlabeled_subflags_from_size(Kn, verbose_output);
    if (dump_unlabeled) dump_unlabeled_flags(Kn);
    if (quit_after_generating_unlabeled)
    {
        cerr << "Quitting after generating unlabeled flags." << endl;
        return 0;
    }

    
    

    if (additional_constraints.size() > 0)
    {
        for (int i = 0; i < (int)additional_constraints.size(); i++)
        {
            if (!load_linear_constraints_from_file(additional_constraints[i], linear_constriants_remove_duplicates_with_types, verbose_output))
            {
                cerr << "Loading additional constraints from " << additional_constraints[i] << " failed" << endl;
                return 0;
            }
        }
    }
    
    if (additional_constraints_strings.size() > 0)
    {
        cerr << "Loadig " << additional_constraints_strings.size()  << " in-command line constraint(s) " << endl;
        for (int i = 0; i < (int)additional_constraints_strings.size(); i++)
        {
            //if (verbose_output)
            {
                cerr << "Loadig constraint from  '" << additional_constraints_strings[i] << "'" << endl;
            }

            
            istringstream ss(additional_constraints_strings[i]);
            
            int tmp=-1;
            ss >> tmp;
            if (tmp != 0)
            {
                cerr << "Constraint should start with 0 but it is not the case for '" << additional_constraints_strings[i] << "'" << endl;
                return 1;
            }
            
            load_linear_constraints_from_stream(ss, linear_constriants_remove_duplicates_with_types, linear_constriants_remove_duplicates_with_types, verbose_output);   
        }
    }
    
    if (g_additional_csdp_blocks.size() > 0)
    {
        for (int i = 0; i < (int)g_additional_csdp_blocks.size(); i++)
        {
            
            ifstream infile;
            infile.open (g_additional_csdp_blocks[i].c_str(), ifstream::in);
            if (!infile.good())
            {
                cerr << "Failed opening file " << g_additional_csdp_blocks[i] << endl;
                return 0;
            }
            int blockKn;
            infile >> blockKn;
            
            if (blockKn != Kn)
            {
                cerr << "Loading additional CSDP blocks from " << g_additional_csdp_blocks[i] << " failed." << endl;
                cerr << "Blocks were computed n=" << blockKn << " while this program is running with n="  << Kn << endl;
                return 0;
            }
            infile.close();
            cerr << "Using additional csdp blocks from " << g_additional_csdp_blocks[i] << endl;
        }
    }

#ifdef G_EDGE_COLOR_SYMMETRY
    add_edge_color_symmetry_constraints(Kn);
#endif

    
    // labeled flags
    g_flags.reserve(MAX_FLAG_TYPES); // Should make live much faster
    if (force_generate_labeled_flags || !load_labeled_flags_from_file(Kn,true))
    {

        // Initialize g_exact_number_of_colored_vertices - currently important just
        // when computing with crossings of bipartite graphs...
        bool check_sensible_flags = false;
        
        cerr << "Generating labeled flags..." << endl;
        // Generating all labeled flags:
        int last_types = (int)g_flags.size();
//// TODO: Make this parallel        
        for (int i = 1; i <= Kn/2; i++)
        {
#ifdef DISABLE_UNLABELED_PRODUCTS
            if (Kn-2*i == 0) continue;
#endif
            cerr << "Getting labeled flags of size " << Kn-i << ":" << Kn-2*i << endl;
            generate_labeled_flags(Kn-i,Kn-2*i, verbose_output, types_from_file, check_sensible_flags);
            int gain = (int)g_flags.size() - last_types;
            cerr << "Got "  << gain << " types ";
            for (int j = last_types; j < (int)g_flags.size(); j++)
                cerr << g_flags[j].size() << " ";
            cerr << endl;
            last_types = (int)g_flags.size();
        }
//#endif
        dump_labeled_flags(Kn);
    }
    cerr << "Labeled flag have " << (int)g_flags.size() << " types. Counts: ";
    for (int i = 0; i < (int)g_flags.size(); i++) cerr << " " << g_flags[i].size();
    cerr << endl;
    if (g_flags.size() == 0)
    {
        cerr << "WARNING: No labeled flags. Computations might not work." << endl;
    }
    
    

    if (!upper_bound && !lower_bound)
    {
        cerr << "I dont know if I should try upper bound our lower bound - you need to tell me using -lb or -ub" << endl;
        return 1;
    }
    
    if (upper_bound && lower_bound)
    {
        cerr << "I'm sorry, but I cannot do both upper and lower bound at the same time." << endl;
        return 1;
    }    

    
    // create result_file
    if (output_file == "" )
    {
        stringstream filename;
        filename << "SDP_n" << Kn << "_";
        if (lower_bound) filename << "LB_";
        else  filename << "UB_";
        filename <<  objective_file << ".dat-s"; 
        output_file = filename.str();
    }
    stringstream filename;
    filename << output_file << ".result";
    string result_file = filename.str();   
    
    
    
    if (process_solution)
    {
        if (output_file == "stdout")
        {
            process_csdp_solution(cin, Kn, process_solution_lower_bound, process_solution_upper_bound, process_solution_print_density);
        } else {
            ifstream results;
            results.open (result_file.c_str(),ifstream::in);
            if (!results.good())
            {
                cerr << "Failed opening file with CSDP result " << result_file << endl;
                return 1;
            }
            cerr << "Processing CSDP solution file " << result_file << endl;
            process_csdp_solution(results, Kn, process_solution_lower_bound, process_solution_upper_bound, process_solution_print_density);            
        }
        return 0;
    }
    
    bool sdp_temp_up_to_date = true;
    if (force_generate_flags || force_generate_labeled_flags)
    {
        sdp_temp_up_to_date = false;
        //cerr << "Todo";
        //assert(0);
    }
        
    if (output_file == "stdout")
    {
        cerr << "Generating SDP program to stdout..." << endl;
        print_CSDP(Kn,upper_bound,cout,verbose_output,use_sdp_temp, sdp_temp_up_to_date);
        
    } else {
        if (output_file == "" )
        {
            stringstream filename;
            filename << "SDP_n" << Kn << "_";
            if (lower_bound) filename << "LB_";
            else  filename << "UB_";
            filename <<  objective_file << ".dat-s"; 
            output_file = filename.str();
        }
        stringstream filename;
        filename << output_file << ".result";
        string result_file = filename.str();
               
        ofstream outfile;
        outfile.open (output_file.c_str(), ofstream::out);
        if (!outfile.good())
        {
            cerr << "Failed opening file " << output_file << endl;
            return -1;
        }
        cerr << "Generating SDP program to " << output_file << endl;
        print_CSDP(Kn,upper_bound,outfile,verbose_output,use_sdp_temp, sdp_temp_up_to_date);
        outfile.close();

        if (csdp_preprocessing != "")
        {
            stringstream ss;
            ss << "python " << csdp_preprocessing << " " << output_file;
            cerr << "Executing preprocessor " << csdp_preprocessing << endl;
            int rv = system(ss.str().c_str());
            if (rv == -1)
            {
                cerr << "Execition failed!" << endl;
                assert(0);
            }
        }

        if (!don_run_csdp)
        {
            if (!csdp_OMP.empty())
            {
                setenv("OMP_NUM_THREADS", csdp_OMP.c_str(), 1);
            }
            

            int csdp_return_value = 1;
           if (system(NULL)) 
           {
                // cout << "Command processor exists";
                assert(use_sdpa == false); 
 
                stringstream system_command;

                string log_file = output_file;
                log_file += ".log";


                system_command <<  "'" << csdp_binary.c_str() << "'  '" << output_file << "'  '" <<  result_file << "'  2>&1  | tee  '"  <<  log_file << "'";
                cerr << "Executing: " << system_command.str() << endl;
                csdp_return_value = system(system_command.str().c_str());

                 //execlp(csdp_binary.c_str(),"csdp",output_file.c_str(),result_file.c_str(),(char *)NULL);
  

           }
            if (csdp_return_value != 0)
           {
               cout << "Command processor doesn't exists or the system call failed."; 
  
                pid_t csdp_PID = 0;
                csdp_PID = fork();
                
                if (csdp_PID == 0)
                {
                    // clean memory! TODO

                    if (use_sdpa)
                    {
                        execlp("sdpa","sdpa","-ds",output_file.c_str(),"-o",result_file.c_str(),(char *)NULL);                
                    }
                    else
                    {
                        execlp(csdp_binary.c_str(),"csdp",output_file.c_str(),result_file.c_str(),(char *)NULL);
                    }

                    cerr << "Strugglig with executing: " << csdp_binary.c_str() << endl;
                    cerr << "Something went wrong: " <<  strerror(errno) << endl;
                }
                else
                {
                    int csdo_status;
                    pid_t tpid = wait(&csdo_status);
                    if (tpid != csdp_PID)
                    {
                        cerr << "Something went wrong!" << endl;
                    }

                }
           }

        }
    }
    
    cerr << "Done." << endl;
    
    
    return 0;    
}
#endif


