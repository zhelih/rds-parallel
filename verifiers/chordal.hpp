#ifndef _CHORDAL_H
#define _CHORDAL_H
#include "verifier.hpp"

class Chordal: public RegisterVerifier<Chordal> {
  private:
    static Chordal instance;
    std::vector<std::vector<unsigned int>> colors;
    std::vector<unsigned int> colors_to_update;

    unsigned int level = 0;

  public:
    bool check_pair(uint i, uint j) const {
      return true;
    }

    bool check(const std::vector<uint>& p, uint n) const {
//      printf("\nChecking if %d and (", n+1);
//      for(unsigned int i = 0; i < p.size(); ++i) {
//        printf("%d, ", p[i]+1);
//      }
//      printf("\b\b) form a chordal graph\n");
      int cnt_bottom[200];
      bool ok_left[200];
      for(int i = 0; i < 200; ++i) {
        ok_left[i] = false;
        cnt_bottom[i] = 0;
      }

      for(uint v: p) {
        for(uint u: p) {
          if (v == u) {
            continue;
          }
//          printf("Checking pair %d, %d\n", v+1, u+1);
          if (!g->is_edge(n, v) || !g->is_edge(n, u)) {
//            printf("Verdict: no edge from %d to one of them\n", n+1);
            continue;
          }
          else {
//            printf("There are edges to both of them\n");
          }

          if (g->is_edge(v, u)) {
//            printf("Verdict: there is edge (%d, %d)\n", v+1, u+1);
            cnt_bottom[v]++;
//            printf("cnt_bottom[%d] now is %d\n", v+1, cnt_bottom[v]);
            continue;
          }

          if (colors[level][v] == colors[level][u]) {
//            printf("Verdict: (%d, %d, %d) is a claw. Solve it later\n", v+1, n+1, u+1);
            continue;
          }
        }
      }

      for(uint v: p) {
        for(uint u: p) {
          if (v == u) {
            continue;
          }
          if (g->is_edge(v, u)) {
//            printf("There is an edge (%d, %d), cnt_bottom[%d] = %d\n", v+1, u+1, v+1, cnt_bottom[v]);
            if (cnt_bottom[v] > 1) {
              ok_left[u] = true;
//              printf("ok_left[%d] = %d\n", u+1, (int)ok_left[u]);
            }
          }
        }
      }

      for(uint v: p) {
        for(uint u: p) {
          if (v == u) {
            continue;
          }
          if (g->is_edge(n, v) && g->is_edge(n, u) && !g->is_edge(v, u) && colors[level][v] == colors[level][u]) {
//            printf("Claw: (%d, %d, %d). ok_left[%d] = %d\n", v+1, n+1, u+1, v+1, (int)ok_left[v]);
            if (!ok_left[v]) {
//              printf("Verdict: no\n");
              return false;
            }
          }
        }
      }

//      printf("Verdict: yes\n");
      return true;
    }

    bool check_solution(const std::vector<uint>& res) const {
      std::vector<uint> res_n = res;
/*    printf("\n Vertices in solution: {");
      for(auto v: res_n) printf("%d, ", v);
      printf("\b\b}\n");
*/
      for(uint i = 0; i < res.size()-1; ++i) {
        bool simplical_found = false;
        for(uint j = i+1; j < res.size(); ++j) {
          uint v = res_n[j];
//          printf("checking if vertex %d is simplicial \n", v);
          bool v_ok = true;

          for(uint z_1 = j+1; z_1 < res.size(); ++z_1) {
            uint z_1_v = res_n[z_1];
            if (!g->is_edge(v, z_1_v)) continue;

            for(uint z_2 = z_1 + 1; z_2 < res.size(); ++z_2) {
              uint z_2_v = res_n[z_2];
              if (!g->is_edge(v, z_2_v)) {
                continue;
              }
              if (!g->is_edge(z_1_v, z_2_v)) {
                v_ok = false;
                break;
              }
            }
            
            if (!v_ok) break;	
          }

          if (v_ok) {
//            printf("Vertex %d is simplicial\n", v);
            std::swap(res_n[i], res_n[j]);
            simplical_found = true;
            break;
          }
        }
//        printf("Vertex %d: simplicial_found is %d\n", res_n[i], simplical_found);
        if (!simplical_found) return false;
      }
      return true;
    }

    void init_aux(uint i, const std::vector<uint>& c) {
      colors.resize(g->nr_nodes); 
      for(auto& v: colors) {
        v.resize(g->nr_nodes, 0);
      }
      for(int i = 0; i < g->nr_nodes; ++i) {
        colors[0][i] = i;
      }

      colors_to_update.reserve(g->nr_nodes);
    }

    void prepare_aux(const std::vector<uint>& p, uint j, const std::vector<uint>& c)
    {
      ++level;
      colors[level] = colors[level-1];
      colors[level][j] = j;
      for(uint v: p) {
        if (g->is_edge(v, j)) {
          colors_to_update.push_back(colors[level][v]);
        }
      }
      for(uint c: colors_to_update) {
        for(uint& ccolor: colors[level]) {
          if (ccolor == c) {
            ccolor = j;
          }
        }
      }
      colors_to_update.resize(0);
    }

    void undo_aux(const std::vector<uint>& p, uint j, const std::vector<uint>& c)
    {
      --level;
    }

    void free_aux() {
      level = 0;
    }

    Chordal() {
      name = "Chordal";
      description = "Well, it's a chordal graph";
      shortcut = "-ch";
    }

    Chordal* clone() const {
        return new Chordal(*this);
    }
};

Chordal Chordal::instance = Chordal();

#endif
