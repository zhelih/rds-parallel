#include "verifier.h"

// clique
bool clique::check_pair(graph* g, uint i, uint j) const { return g->is_edge(i,j); }
bool clique::check(graph* g, const std::vector<uint>& p, uint n, void* _aux) const
{
  for(auto it = p.begin(); it != p.end(); ++it)
    if(!g->is_edge(*it, n))
      return false;
  return true;
}
bool clique::check_solution(graph* g, const std::vector<uint>& res) const
{
  for(uint i = 0; i < res.size(); ++i)
    for(uint j = i+1; j > res.size(); ++j)
      if(!g->is_edge(res[i], res[j]))
        return false;
  return true;
}

// stable
bool stable::check_pair(graph* g, uint i, uint j) const { return !g->is_edge(i,j); }
bool stable::check(graph* g, const std::vector<uint>& p, uint n, void* _aux) const
{
  for(auto it = p.begin(); it != p.end(); ++it)
    if(g->is_edge(*it, n))
      return false;
  return true;
}
bool stable::check_solution(graph* g, const std::vector<uint>& res) const
{
  for(uint i = 0; i < res.size(); ++i)
    for(uint j = i+1; j > res.size(); ++j)
      if(g->is_edge(res[i], res[j]))
        return false;
  return true;
}

// iuc
bool iuc::check_pair(graph* g, uint i, uint j) const { return true; }
bool iuc::check(graph* g, const std::vector<uint>& p, uint n, void* _aux) const
{
  for(auto it = p.begin(); it != p.end(); ++it)
    for(auto it2 = next(it); it2 != p.end(); ++it2)
    if(g->is_edge(*it, *it2) + g->is_edge(*it, n) + g->is_edge(*it2, n) == 2)
      return false;
  return true;
}
bool iuc::check_solution(graph* g, const std::vector<uint>& res) const {
  for(uint i = 0; i < res.size(); ++i)
    for(uint j = i+1; j < res.size(); ++j)
      for(uint k = j+1; k < res.size(); ++k)
      {
        uint ij = g->is_edge(res[i], res[j]);
        uint ik = g->is_edge(res[i], res[k]);
        uint jk = g->is_edge(res[j], res[k]);
        if(ij + ik + jk == 2)
          return false;
      }
  return true;
}

//s-defective clique
bool defective_clique::check_pair(graph* g, uint i, uint j) const { return (s>0)||(g->is_edge(i,j)); }
bool defective_clique::check(graph* g, const std::vector<uint>& p, uint n, void* aux) const
{
  t_aux* a = (static_cast<t_aux*>(aux));
  return a->nnv + a->nncnt[n] <= s;
}
void* defective_clique::init_aux(graph* g, uint i, const std::vector<uint>& c)
{
  t_aux *aux = new t_aux();
  aux->nnv=0;
  aux->nncnt.resize(g->nr_nodes);
  for(uint it = 0; it < c.size(); ++it)
  {
    aux->nncnt[c[it]]=!g->is_edge(c[it], i);
  }
  return aux;
}
void defective_clique::prepare_aux(graph* g, const std::vector<uint>& p, uint j, const std::vector<uint>& c, void* aux)
{
  t_aux* a = static_cast<t_aux*>(aux);

  a->nnv += a->nncnt[j];
  for(uint i = 0; i < c.size(); ++i)
  {
    if(!g->is_edge(c[i],j))
      a->nncnt[c[i]]++;
  }
}
void defective_clique::undo_aux(graph* g, const std::vector<uint>& p, uint j, const std::vector<uint>& c, void* aux)
{
  t_aux* a = static_cast<t_aux*>(aux);

  for(uint i = 0; i < c.size(); ++i)
  {
    if(!g->is_edge(c[i],j))
      a->nncnt[c[i]]--;
  }
  a->nnv -= a->nncnt[j];
}
void defective_clique::free_aux(void* aux) { t_aux* i = static_cast<t_aux*>(aux); delete i; }
bool defective_clique::check_solution(graph* g, const std::vector<uint>& res) const
{
  uint edges = 0;
  for(uint i = 0; i < res.size(); ++i)
    for(uint j = i+1; j < res.size(); ++j)
      if(g->is_edge(res[i],res[j]))
        edges++;
  return edges >= (res.size()*(res.size()-1)/2 - s);
}

//s-plex
bool plex::check_pair(graph* g, uint i, uint j) const
{
  if(s == 1)
    return g->is_edge(i,j);
  else
    return true;
}
bool plex::check_solution(graph* g, const std::vector<uint>& res) const
{
  std::vector<uint> degrees(res.size());
  for(uint i = 0; i < res.size(); ++i)
    degrees[i] = 0;
  for(uint i = 0; i < res.size(); ++i)
    for(uint j = 0; j < res.size(); ++j)
      if(i != j && g->is_edge(res[i],res[j]))
        degrees[i]++;
  uint m_degree = degrees[0];
  for(uint i = 1; i < res.size(); ++i)
    if(m_degree > degrees[i])
      m_degree = degrees[i];
  printf("min_degree = %u\n", m_degree);
  return m_degree >= (res.size() - s);
}

bool plex::check(graph* g, const std::vector<uint>& p, uint n, void* aux) const
{
  t_aux* a = static_cast<t_aux*>(aux);
  if(a->nncnt[n] >= s) // degree check
    return false;
  for(uint i = 0; i < a->nr_sat; ++i) // SAT connectivity check
    if(!g->is_edge(n, a->sat[i]))
      return false;
  return true;
}
void* plex::init_aux(graph* g, uint i, const std::vector<uint>& c)
{
  t_aux *aux = new t_aux();
  aux->nncnt.resize(g->nr_nodes);
  aux->sat.resize(g->nr_nodes);
  aux->nr_sat = 0;
  for(uint it = 0; it < g->nr_nodes; ++it)
    aux->nncnt[it] = 0;
  for(uint it = 0; it < c.size(); ++it)
    aux->nncnt[c[it]]=!g->is_edge(c[it], i);
  return aux;
}
void plex::prepare_aux(graph*g, const std::vector<uint>& p, uint j, const std::vector<uint>& c, void* aux)
{
  t_aux* a = static_cast<t_aux*>(aux);
  a->nr_sat = 0;
  for(uint i = 0; i < p.size(); ++i)
  {
    if(!g->is_edge(p[i],j))
    {
      a->nncnt[p[i]]++;
      if(a->nncnt[p[i]] == s-1)
      {
        a->sat[a->nr_sat] = p[i]; // push_back
        a->nr_sat++;
      }
    }
  }

  for(uint i = 0; i < c.size(); ++i)
  {
    if(!g->is_edge(c[i],j))
    {
      a->nncnt[c[i]]++;
    }
  }
  if(a->nncnt[j] == s-1)
  {
    a->sat[a->nr_sat] = j;
    a->nr_sat++;
  }
}
void plex::undo_aux(graph* g, const std::vector<uint>& p, uint j, const std::vector<uint>& c, void* aux)
{
  t_aux* a = static_cast<t_aux*>(aux);

  for(uint i = 0; i < c.size(); ++i)
  {
    if(!g->is_edge(c[i],j))
    {
      a->nncnt[c[i]]--;
    }
  }
  for(uint i = 0; i < p.size(); ++i)
  {
    if(!g->is_edge(p[i],j))
    {
      a->nncnt[p[i]]--;
    }
  }
}
void plex::free_aux(void* aux) { t_aux* a = static_cast<t_aux*>(aux); delete a; }
