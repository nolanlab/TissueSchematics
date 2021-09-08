#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <utility>
#include <map>
#include <vector>
#include <random>
#include <algorithm>


namespace py = pybind11;


class ColorMC{
public:
  std::map<int, int> initial_coloring;
  std::map<int, std::vector<int>> neighbors;
  std::map<int, std::vector<int>> initial_node_to_class;
  std::map<std::vector<int>, std::set<int>> initial_class_to_nodeset;

  
    // get isomorphism class function
  virtual std::vector<int> getIsom(int node, const std::map<int,int>& coloring ){
    auto neighbor_nodes = this->neighbors.at(node);

    std::vector<int> out(neighbor_nodes.size());

    for (int i = 0; i < neighbor_nodes.size(); i++){
      out[i] = coloring.at(neighbor_nodes[i]);
    }
    std::sort(out.begin(), out.end());
    return out;
  };



  
  // constructor requires graph as map of nodes to vector of nodes (neighbors), coloring
  ColorMC(std::map<int, std::vector<int>> neighbors, std::map<int, int> coloring){

    this-> neighbors = neighbors;
    this-> initial_coloring = coloring;

    // build initial_node_to_class and initial_class_to_nodeset
    
    for (auto node_and_neighbors: neighbors){
        
      auto node = node_and_neighbors.first;
      auto isom = getIsom(node, coloring);
      auto kvpair = std::pair<int, std::vector<int>>(node,isom);
      this->initial_node_to_class.insert(kvpair);
	
      if (initial_class_to_nodeset.count(isom)==0){
	auto kvpair = std::pair<std::vector<int>, std::set<int>>(isom, std::set<int>());
	this->initial_class_to_nodeset.insert(kvpair);

      };
      
      this->initial_class_to_nodeset.at(isom).insert(node);
    }
  }

 void transposeColoring(std::pair<int,int> transposition, std::map<int,int> &coloring){
   
    int node1 = transposition.first;
    int node2 = transposition.second;

    int color1 = coloring.at(node1);
    int color2 = coloring.at(node2);

    coloring.at(node1) = color2;
    coloring.at(node2) = color1;

  }
  
  void updateQuotient(std::pair<int,int> transposition,std::map<int,int> &coloring, std::map<int, std::vector<int>> &node_to_class, std::map<std::vector<int>, std::set<int>> &class_to_nodeset){

    int node1 = transposition.first;
    int node2 = transposition.second;

    std::set<int> nodes_to_update;
    for (int a: this->neighbors.at(node1)){
      nodes_to_update.insert(a);
    }

    for (int a: this->neighbors.at(node2)){
      nodes_to_update.insert(a);
    }

    transposeColoring(transposition,coloring);
    for (int node: nodes_to_update){
      auto old_isom = node_to_class.at(node);
      auto new_isom = getIsom(node, coloring);

      node_to_class.at(node) = new_isom;
      class_to_nodeset.at(old_isom).erase(node);
      
      if (class_to_nodeset.count(new_isom) ==0){
	class_to_nodeset.insert(std::pair<std::vector<int>, std::set<int>>(new_isom, std::set<int>()));
      }

      class_to_nodeset.at(new_isom).insert(node);
      
    }
      
  }

  std::map<int, std::set<int>> findCandTranspositions(std::map<int,int> &coloring, std::map<int, std::vector<int>> &node_to_class, std::map<std::vector<int>, std::set<int>> &class_to_nodeset){
    // for each node, builds a set of possible swaps (those with same isom and different color)
    std::map<int, std::set<int>> cand_swaps;
    for (auto node_class: node_to_class){
      int node = node_class.first;
      std::vector<int> isom = node_class.second;
      std::set<int> node_cand_swaps;
      for (int cand:class_to_nodeset.at(isom)){	
	if (coloring.at(cand)!=coloring.at(node)){
	  node_cand_swaps.insert(cand);
	}
      }

      cand_swaps.insert(std::pair<int,std::set<int>>(node,node_cand_swaps));
      
    }
    return cand_swaps;
  }


  std::vector<std::map<int,int>> runMC(int num_steps, int save_every, int burn_steps){

    std::vector<std::map<int,int>> out_colorings;
    std::default_random_engine gen;
    
    
    auto coloring = this->initial_coloring;
    auto node_to_class = this->initial_node_to_class;
    auto class_to_nodeset = this->initial_class_to_nodeset;



    auto cand_swaps = findCandTranspositions(coloring,node_to_class,class_to_nodeset);
    std::vector<int> weights;
    std::vector<int> nodes;
    int num_viable = 0;
    for (auto kvpair: cand_swaps){
      nodes.push_back(kvpair.first);
      int viable = kvpair.second.size()>0;
      weights.push_back(viable);
      if (viable==1){
	num_viable++;
      }
     }

    if (num_viable==0){
      py::print("no viable transpositions at step 0");
      return out_colorings;
    }

    for (int i=0; i<num_steps; i++){

      // sample a node  with viable swap
      std::discrete_distribution<> d1(weights.begin(), weights.end());
      int sampled_node1 = nodes[d1(gen)];

      //sample one in its isom class with a different color
      std::vector<int> swapset(cand_swaps.at(sampled_node1).size());
      std::vector<int> unif_probs(cand_swaps.at(sampled_node1).size());
      int j = 0;
      for (int a:cand_swaps.at(sampled_node1)){
	swapset[j] = a;
	unif_probs[j] = 1;
	j++;
      }
      std::discrete_distribution<> d2(unif_probs.begin(),unif_probs.end());      
      int sampled_node2 = swapset[d2(gen)];

      //generate the transposition
      std::pair<int,int> transposition(sampled_node1,sampled_node2);

      // update everything
      updateQuotient(transposition,coloring,node_to_class,class_to_nodeset);


      // work out going backwards for MH
      auto back_cand_swaps = findCandTranspositions(coloring,node_to_class,class_to_nodeset);
      std::vector<int> back_weights;
      std::vector<int> back_nodes;
      int back_num_viable = 0;
      for (auto kvpair: back_cand_swaps){
	back_nodes.push_back(kvpair.first);
	int viable = kvpair.second.size()>0;
	back_weights.push_back(viable);
	if (viable==1){
	back_num_viable++;
	}
      }

      // compute the MH acceptance probability
      double g_forwards = (1.0/num_viable) *((1.0/cand_swaps.at(sampled_node1).size()) + (1.0/cand_swaps.at(sampled_node2).size()));
      double g_backwards = (1.0/back_num_viable) *((1.0/back_cand_swaps.at(sampled_node1).size()) + (1.0/back_cand_swaps.at(sampled_node2).size()));
      double accept = g_backwards/g_forwards;

      
      std::uniform_real_distribution<double> unif(0.0,1.0);
      if (unif(gen) <accept){
	num_viable = back_num_viable;
	nodes = back_nodes;
	weights = back_weights;
	cand_swaps = back_cand_swaps;
      }
      else {
	// better to just re-update with transposition than copy everything
	updateQuotient(transposition,coloring,node_to_class,class_to_nodeset);
      }

      // save the coloring if appropriate
      if ((i%save_every ==0) && (i>burn_steps)){
	out_colorings.push_back(coloring);
    }
    }
      return out_colorings;
  }
    
    
  
  

  // tests
  std::vector<std::pair<std::vector<int>, std::set<int>>>  cast_eq(){
      std::vector<std::pair<std::vector<int>, std::set<int>>>  out;
    for (auto a:this->initial_class_to_nodeset){
      out.push_back(a);
    }
    return out;
  }

  void test_update(int node1, int node2){
    py::print("before");
    py::print(this->initial_coloring);
    py::print(cast_eq());
    updateQuotient(std::pair<int,int>(node1,node2), this->initial_coloring, this->initial_node_to_class, this->initial_class_to_nodeset);
    py::print("after");
    py::print(this->initial_coloring);
    py::print(cast_eq());
    
  }

  
  std::map<int, std::vector<int>> get_ntc(){return this->initial_node_to_class;}
    
};

class ColorMCSet : public ColorMC{
public:
  using ColorMC::ColorMC;
  //override isomorphism function
    std::vector<int> getIsom(int node, const std::map<int,int>& coloring ) override {
    auto neighbor_nodes = this->neighbors.at(node);

    std::vector<int> out(neighbor_nodes.size());

    for (int i = 0; i < neighbor_nodes.size(); i++){
      out[i] = coloring.at(neighbor_nodes[i]);
    }
    
    std::sort(out.begin(), out.end());
    auto last = std::unique(out.begin(),out.end());
    out.erase(last, out.end());
    
    return out;
    }

};


PYBIND11_MODULE(colormc, m) {
  py::class_<ColorMC>(m, "ColorMC").def(py::init<const std::map<int,std::vector<int>>&, const std::map<int,int>>())
    .def("runMC", &ColorMC::runMC)
    .def("getIsom", &ColorMC::getIsom)
    .def("test_update", &ColorMC::test_update)
    .def("get_ntc", &ColorMC::get_ntc);

 py::class_<ColorMCSet>(m, "ColorMCSet").def(py::init<const std::map<int,std::vector<int>>&, const std::map<int,int>>())
    .def("runMC", &ColorMCSet::runMC)
    .def("getIsom", &ColorMCSet::getIsom)
    .def("test_update", &ColorMCSet::test_update)
    .def("get_ntc", &ColorMCSet::get_ntc);

  

}
