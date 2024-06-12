#include "exODT.h"
#include "fractionMissing.h"

using namespace std;
using namespace bpp;
static scalar_type EPSILON = numeric_limits<scalar_type>::min();

void exODT_model::construct_undated(const string &Sstring,
                                    const string &fractionMissingFile) {
  daughter.clear();
  son.clear();
  name_node.clear();
  node_name.clear();
  node_ids.clear();
  id_nodes.clear();

  string_parameter["S_un"] = Sstring;
  S = TreeTemplateTools::parenthesisToTree(
      string_parameter["S_un"],
      true); //(string_parameter["BOOTSTRAP_LABELS"]=="yes")

  S_root = S->getRootNode();
  vector<Node *> nodes = TreeTemplateTools::getNodes(*S_root);

  for (vector<Node *>::iterator it = nodes.begin(); it != nodes.end(); it++)
    if ((*it)->isLeaf()) {
      name_node[(*it)->getName()] = (*it);
      node_name[(*it)] = (*it)->getName();
    } else {
      vector<string> leafnames = TreeTemplateTools::getLeavesNames(*(*it));
      sort(leafnames.begin(), leafnames.end());
      stringstream name;
      for (vector<string>::iterator st = leafnames.begin();
           st != leafnames.end(); st++)
        name << (*st) << ".";

      name_node[name.str()] = (*it);
      node_name[(*it)] = name.str();
    }
  // register species
  last_branch = 0;
  last_leaf = 0;

  set<Node *> saw;
  for (map<string, Node *>::iterator it = name_node.begin();
       it != name_node.end(); it++)
    if ((*it).second->isLeaf()) {
      Node *node = (*it).second;
      extant_species[last_branch] = node->getName();
      // stringstream name;
      // name << extant_species[last_branch] <<"("<< last_branch<<")";
      // node->setName(name);
      node_ids[node] = last_branch;
      id_nodes[last_branch] = node;
      last_branch++;
      last_leaf++;
      saw.insert(node);
      // a leaf
      daughter[last_branch] = -1;
      // a leaf
      son[last_branch] = -1;
      vector_parameter["BL_rate_multiplier"].push_back(
          node->getDistanceToFather());
      vector_parameter["rate_multiplier_tau_to"].push_back(1);
      vector_parameter["rate_multiplier_tau_from"].push_back(1);
      wT.push_back(1);
      rmD.push_back(1);
      rmT.push_back(1);
      rmL.push_back(1);
      vector_parameter["rate_multiplier_delta"].push_back(1);
      vector_parameter["rate_multiplier_lambda"].push_back(1);
      vector_parameter["rate_multiplier_O"].push_back(1);
    }

  // ad-hoc postorder
  vector<Node *> next_generation;
  for (map<string, Node *>::iterator it = name_node.begin();
       it != name_node.end(); it++)
    if ((*it).second->isLeaf()) {
      Node *node = (*it).second;
      next_generation.push_back(node);
    }
  while (next_generation.size()) {
    vector<Node *> new_generation;
    for (vector<Node *>::iterator it = next_generation.begin();
         it != next_generation.end(); it++) {
      Node *node = (*it);
      if (node->hasFather()) {
        Node *father = node->getFather();
        vector<Node *> sons = father->getSons();
        Node *sister;
        if (sons[0] == node)
          sister = sons[1];
        else
          sister = sons[0];

        if (not node_ids.count(father) and saw.count(sister)) {
          node_ids[father] = last_branch;
          id_nodes[last_branch] = father;
          stringstream name;
          name << last_branch;

          vector_parameter["BL_rate_multiplier"].push_back(
              node->getDistanceToFather());
          vector_parameter["rate_multiplier_tau_to"].push_back(1);
          vector_parameter["rate_multiplier_tau_from"].push_back(1);
          wT.push_back(1);
          rmD.push_back(1);
          rmT.push_back(1);
          rmL.push_back(1);
          vector_parameter["rate_multiplier_delta"].push_back(1);
          vector_parameter["rate_multiplier_lambda"].push_back(1);
          vector_parameter["rate_multiplier_O"].push_back(1);

          father->setBranchProperty("ID", BppString(name.str()));
          last_branch++;

          saw.insert(father);
          new_generation.push_back(father);
        }
      }
    }
    next_generation.clear();
    for (vector<Node *>::iterator it = new_generation.begin();
         it != new_generation.end(); it++)
      next_generation.push_back((*it));
  }

  // for (vector <Node * >::iterator it=nodes.begin();it!=nodes.end();it++ )
  // (*it)->setDistanceToFather(1);
  below.clear();
  for (int e = 0; e < last_branch - 1; e++) {
    for (int f = 0; f < last_branch - 1; f++) {
      // cout <<height(id_nodes[e]->getFather()) <<" e:"<<e<<" "<<" "<<
      // height(id_nodes[f])  <<" f:"<<f<<" "<<" "<< endl;
      if (height(id_nodes[e]->getFather()) < height(id_nodes[f])) {
        // cout << e << " below " << f << endl;
        below[e][f] = 1;
      } else
        below[e][f] = 0;
    }
  }

  vector_parameter["BL_rate_multiplier"][last_branch] =
      scalar_parameter["root_BL"];
  vector_parameter["rate_multiplier_tau_to"].push_back(1);
  vector_parameter["rate_multiplier_tau_from"].push_back(1);
  wT.push_back(1);
  vector_parameter["rate_multiplier_delta"].push_back(1);
  vector_parameter["rate_multiplier_lambda"].push_back(1);
  vector_parameter["rate_multiplier_O"].push_back(1);

  for (map<Node *, int>::iterator it = node_ids.begin(); it != node_ids.end();
       it++) {
    Node *node = (*it).first;
    int branch = (*it).second;
    stringstream out;
    stringstream out1;
    stringstream out2;
    out1 << t_begin[branch];
    out2 << t_end[branch];
    int rank = branch;
    out << rank;
    if (node->hasBranchProperty("bootstrap")) {
      rank2label[rank] = node->getBootstrapValue();
      // cout <<rank2label[rank]<<"->"<<rank<< endl;
    } else {
      rank2label[rank] = -1;
    }
    node->setBranchProperty("ID", BppString(out.str()));
  }

  string_parameter["S_with_ranks"] =
      TreeTemplateTools::treeToParenthesis(*S, false, "ID");

  // map <string,map<string,int> > ancestral_names;
  // map <int,map<int,int> > ancestral;
  ancestors.clear();
  for (int e = 0; e < last_branch; e++) {
    vector<int> tmp;
    ancestors.push_back(tmp);
    for (int f = 0; f < last_branch; f++)
      ancestral[e][f] = 0;
  }

  for (vector<Node *>::iterator it = nodes.begin(); it != nodes.end(); it++) {
    Node *node = (*it);
    int e = node_ids[node];
    stringstream name_from;
    if (e < last_leaf)
      name_from << node_name[node];
    else
      name_from << e;
    map<string, int> tmp;
    ancestral_names[name_from.str()] = tmp;
    while (node->hasFather()) {
      stringstream name_to;
      int f = node_ids[node];
      if (f < last_leaf)
        name_to << node_name[node];
      else
        name_to << f;
      node = node->getFather();
      ancestral_names[name_from.str()][name_to.str()] = 1;
      if (not ancestral[e][f])
        ancestors[e].push_back(f);
      ancestral[e][f] = 1;
    }
    stringstream name_to;
    int f = node_ids[node];
    name_to << f;
    ancestral_names[name_from.str()][name_to.str()] = 1;
    if (not ancestral[e][f])
      ancestors[e].push_back(f);
    ancestral[e][f] = 1;
  }
  if (scalar_parameter["reldate"] == true)
    for (int e = 0; e < last_branch; e++) {
      for (int f = 0; f < last_branch; f++)
        if (below[e][f] == 1) {
          if (not ancestral[e][f])
            ancestors[e].push_back(f);
          ancestral[e][f] = 1;
        }
    }

  for (map<string, Node *>::iterator it = name_node.begin();
       it != name_node.end(); it++)
    if (not(*it).second->isLeaf()) {
      Node *node = (*it).second;
      vector<Node *> sons = node->getSons();
      daughter[node_ids[node]] = node_ids[sons[0]];
      son[node_ids[node]] = node_ids[sons[1]];
      // cout << node_ids[node] << " => " << node_ids[sons[0]] << " & " <<
      // node_ids[sons[1]] << endl; cout << node_name[node] << " => " <<
      // node_name[sons[0]] << " & " << node_name[sons[1]] << endl;
    }
  branch_counts["Os"].clear();
  branch_counts["Ds"].clear();
  branch_counts["Ts"].clear();
  branch_counts["Tfroms"].clear();
  branch_counts["Ls"].clear();
  branch_counts["count"].clear();
  branch_counts["presence"].clear();
  branch_counts["saw"].clear();
  branch_counts["O_LL"].clear();
  branch_counts["copies"].clear();
  branch_counts["singleton"].clear();

  for (int e = 0; e < last_branch; e++) {
    branch_counts["Os"].push_back(0);
    branch_counts["Ds"].push_back(0);
    branch_counts["Ts"].push_back(0);
    branch_counts["Tfroms"].push_back(0);
    branch_counts["Ls"].push_back(0);
    branch_counts["count"].push_back(0);
    branch_counts["presence"].push_back(0);
    branch_counts["saw"].push_back(0);
    branch_counts["O_LL"].push_back(0);
    branch_counts["copies"].push_back(0);
    branch_counts["singleton"].push_back(0);
  }
  T_to_from.clear();
  for (int e = 0; e < last_branch; e++) {
    vector<scalar_type> tmp;
    T_to_from.push_back(tmp);
    for (int f = 0; f < last_branch; f++)
      T_to_from[e].push_back(0);
  }

  last_rank = last_branch;
  set_model_parameter("N", 1);

  // Put default values for the fraction of missing genes at the leaves.
  vector_parameter["fraction_missing"] = vector<scalar_type>(last_leaf, 0.0);
  // Put user-defined values, if available
  if (fractionMissingFile == "") {

  } else {
    fraction_missing = readFractionMissingFile(fractionMissingFile);
    // Now we need to fill up the vector_parameter, and we have to be careful
    // about the order.
    size_t index = 0;
    for (map<string, Node *>::iterator it = name_node.begin();
         it != name_node.end(); it++) {
      if ((*it).second->isLeaf()) {
        Node *node = (*it).second;
        string currentSpecies = node->getName();
        vector_parameter["fraction_missing"][index] =
            fraction_missing[currentSpecies];
        index++;
      }
    }
    VectorTools::print(vector_parameter["fraction_missing"]);
  }
}

void exODT_model::calculate_undatedEs() {
  uE.clear();
  fm.clear();
  mPTE_ancestral_correction.clear();
  PD.clear();
  tau_norm.clear();
  PL.clear();
  PS.clear();
  // scalar_type P_T=0;
  map<string, scalar_type> rm_norms;
  rm_norms["tau_to"] = 0;
  rm_norms["tau_from"] = 0;
  rm_norms["delta"] = 0;
  rm_norms["lambda"] = 0;
  rm_norms["O"] = 0;

  for (int e = 0; e < last_branch; e++) {
    rm_norms["tau_to"] += vector_parameter["rate_multiplier_tau_to"][e];
    rm_norms["tau_from"] += vector_parameter["rate_multiplier_tau_from"][e];
    rm_norms["delta"] += vector_parameter["rate_multiplier_delta"][e];
    rm_norms["lambda"] += vector_parameter["rate_multiplier_lambda"][e];
    rm_norms["O"] += vector_parameter["rate_multiplier_O"][e];
  }

  for (int e = 0; e < last_branch; e++) {
    if (scalar_parameter["undatedBL"] == true) {
      wT[e] = vector_parameter["rate_multiplier_tau_to"][e] *
              vector_parameter["BL_rate_multiplier"][e];
      rmT[e] = vector_parameter["tau"][e] *
               vector_parameter["rate_multiplier_tau_from"][e] *
               vector_parameter["BL_rate_multiplier"][e];
      rmD[e] = vector_parameter["delta"][e] *
               vector_parameter["rate_multiplier_delta"][e] *
               vector_parameter["BL_rate_multiplier"][e];
      rmL[e] = vector_parameter["lambda"][e] *
               vector_parameter["rate_multiplier_lambda"][e] *
               vector_parameter["BL_rate_multiplier"][e];
    } else {
      wT[e] = vector_parameter["rate_multiplier_tau_to"][e];
      rmT[e] = vector_parameter["tau"][e] *
               vector_parameter["rate_multiplier_tau_from"][e];
      rmD[e] = vector_parameter["delta"][e] *
               vector_parameter["rate_multiplier_delta"][e];
      rmL[e] = vector_parameter["lambda"][e] *
               vector_parameter["rate_multiplier_lambda"][e];
    }
  }
  scalar_type tau_sum = 0;
  for (int f = 0; f < last_branch; f++)
    tau_sum += wT[f];

  for (int e = 0; e < last_branch; e++) {
    // cout << e << " is " << node_name[id_nodes[e]] << endl;
    scalar_type P_D = rmD[e];
    scalar_type P_T = rmT[e];
    scalar_type P_L = rmL[e];
    scalar_type P_S = 1;

    scalar_type tmp = P_D + P_T + P_L + P_S;
    P_D /= tmp;
    P_T /= tmp;
    P_L /= tmp;
    P_S /= tmp;
    PD.push_back(P_D);
    tau_norm.push_back(tau_sum);
    for (vector<int>::iterator it = ancestors[e].begin();
         it != ancestors[e].end(); it++) {
      int f = (*it);
      tau_norm[e] -= wT[f];
    }
    tau_norm[e] /= P_T;

    PL.push_back(P_L);
    PS.push_back(P_S);
    uE.push_back(0);
    if (e < last_leaf) { // we are at a leaf
      fm.push_back(vector_parameter["fraction_missing"][e]);
    } else {
      fm.push_back(0);
    }
    // uE.push_back(fm[e]);

    mPTE_ancestral_correction.push_back(0);
  }

  // In the loop below with 4 iterations, we calculate the mean probability mPTE
  // for a gene to become extinct across all branches.
  mPTE = 0;
  for (int i = 0; i < 4; i++) {
    scalar_type newmPTE = 0;
    // vector<scalar_type> ancestral_correction;
    if (i > 0) // There should be no need for this loop at the first iteration,
               // because then it leaves mPTE_ancestral_correction at 0.
    {
      for (int e = 0; e < last_branch; e++) {
        mPTE_ancestral_correction[e] = 0;
        // for (map<int,int>::iterator it=ancestral[e].begin();(
        // it!=ancestral[e].end() and i>0);it++)
        for (vector<int>::iterator it = ancestors[e].begin();
             it != ancestors[e].end(); it++) {
          // int f=(*it).first;
          int f = (*it);
          // if (ancestral[e][f]==1)
          mPTE_ancestral_correction[e] +=
              (wT[f]) *
              uE[f]; // That's how we forbid transfers to ancestors of a branch
        }
      }
    }

    for (int e = 0; e < last_branch; e++) {
      if (e < last_leaf) // we are at a leaf, there cannot be a speciation
                         // event, but the gene has been lost
      {
        uE[e] = PL[e] + PD[e] * uE[e] * uE[e] +
                uE[e] * (mPTE - mPTE_ancestral_correction[e]) / tau_norm[e];
      } else // Not at a leaf: the gene was lost once on branch e, or on all
             // descendants, including after speciation
      {
        int f = daughter[e];
        int g = son[e];
        uE[e] = PL[e] + PS[e] * uE[f] * uE[g] + PD[e] * uE[e] * uE[e] +
                uE[e] * (mPTE - mPTE_ancestral_correction[e]) / tau_norm[e];
      }
      newmPTE += (wT[e]) * uE[e];
    }
    mPTE = newmPTE;
  } // End of the loop to compute mPTE

  // Now we add one more update of mPTE to take into account the fraction of
  // missing genes, which had been ignored so far.

  scalar_type newmPTE = 0;
  for (int e = 0; e < last_branch; e++) {
    if (e < last_leaf) // we are at a leaf: either the gene has been lost, or it
                       // is one of those missing genes
    {
      uE[e] = (1 - fm[e]) * uE[e] + fm[e];
    } else // Not at a leaf: the gene was lost once on branch e, or on all
           // descendants
    {
      int f = daughter[e];
      int g = son[e];
      uE[e] = PL[e] + PS[e] * uE[f] * uE[g] + PD[e] * uE[e] * uE[e] +
              uE[e] * (mPTE - mPTE_ancestral_correction[e]) / tau_norm[e];
    }
    newmPTE += (wT[e]) * uE[e];
  }
  mPTE = newmPTE;
}

scalar_type exODT_model::pun(approx_posterior *ale, bool verbose, bool no_T) {
  scalar_type survive = 0;
  scalar_type root_sum = 0;
  scalar_type O_norm = 0;
  mPTuq_ancestral_correction.clear();
  uq.clear();
  mPTuq.clear(); // XX
  ale_pointer = ale;

  for (std::map<long int,
                std::map<scalar_type, std::map<int, scalar_type>>>::iterator
           it = q.begin();
       it != q.end(); it++) {
    for (std::map<scalar_type, std::map<int, scalar_type>>::iterator jt =
             (*it).second.begin();
         jt != (*it).second.end(); jt++)
      (*jt).second.clear();
    (*it).second.clear();
  }
  q.clear();

  // directed partitions and their sizes
  // vector <long int>  g_ids;
  // vector <long int>  g_id_sizes;
  g_ids.clear();
  g_id_sizes.clear();

  for (map<int, vector<long int>>::iterator it = ale->size_ordered_bips.begin();
       it != ale->size_ordered_bips.end(); it++)
    for (vector<long int>::iterator jt = (*it).second.begin();
         jt != (*it).second.end(); jt++) {
      g_ids.push_back((*jt));
      g_id_sizes.push_back((*it).first);
    }
  // root bipartition needs to be handled separately
  g_ids.push_back(-1);
  g_id_sizes.push_back(ale->Gamma_size);

  root_i = g_ids.size() - 1;

  // gene<->species mapping
  // mpi_tree's pun() breaks if mapping is saved.. Sz.G. 12.29/
  gid_sps.clear();
  if (gid_sps.size() == 0) // If the mapping has not been done yet
  {
    // Test that the species associated to genes are really in the species tree
    std::set<string> species_set;
    for (std::map<int, string>::iterator iter = extant_species.begin();
         iter != extant_species.end(); ++iter) {
      species_set.insert(iter->second);
    }

    if (verbose)
      cout << "\nGene"
           << "\t:\t"
           << "Species" << endl;
    for (int i = 0; i < (int)g_ids.size(); i++) {
      long int g_id = g_ids[i];

      if (g_id_sizes[i] == 1) {
        int id = 0;
        for (auto i = 0; i < ale->Gamma_size + 1; ++i) {
          if (ale->id_sets[g_id][i]) {
            id = i;
            break;
          }
        }

        string gene_name = ale->id_leaves[id];
        vector<string> tokens;
        boost::split(tokens, gene_name,
                     boost::is_any_of(string_parameter["gene_name_separators"]),
                     boost::token_compress_on);
        string species_name;
        if ((int)scalar_parameter["species_field"] == -1) {
          species_name = tokens[1];
          for (int fi = 2; fi < tokens.size(); fi++)
            species_name += "_" + tokens[fi];
        }
        //	    species_name=tokens[tokens.size()-1];
        else
          species_name = tokens[(int)scalar_parameter["species_field"]];
        gid_sps[g_id] = species_name;
        if (species_set.find(species_name) == species_set.end()) {
          cout << "Error: gene name " << gene_name
               << " is associated to species name " << species_name
               << " that cannot be found in the species tree." << endl;
          exit(-1);
        }
        if (verbose)
          cout << gene_name << "\t:\t" << species_name << endl;
      }
    }
  }
  // map <long int, int> g_id2i;
  // XX ancestral_correction ..
  for (int i = 0; i < (int)g_ids.size(); i++) {
    long int g_id = g_ids[i];
    g_id2i[g_id] = i;

    if (not(i < (int)uq.size())) {
      vector<scalar_type> tmp;
      uq.push_back(tmp);
      vector<scalar_type> tmp2;
      mPTuq_ancestral_correction.push_back(tmp2);

      mPTuq.push_back(0);
    } else
      mPTuq[i] = 0;

    for (int e = 0; e < last_branch; e++)
      if (not(e < (int)uq[i].size())) {
        uq[i].push_back(0);
        mPTuq_ancestral_correction[i].push_back(0);
      } else {
        uq[i][e] = 0;
        mPTuq_ancestral_correction[i][e] = 0;
      }
  }

  for (int iter = 0; iter < 4; iter++) {
    for (int i = 0; i < (int)g_ids.size(); i++) {

      scalar_type new_mPTuq = 0;

      // directed partition (dip) gamma's id
      bool is_a_leaf = false;
      long int g_id = g_ids[i];
      if (g_id_sizes[i] == 1)
        is_a_leaf = true;
      vector<int> gp_is;
      vector<long int> gpp_is;
      vector<scalar_type> p_part;
      if (g_id != -1)
        for (unordered_map<pair<long int, long int>, scalar_type>::iterator kt =
                 ale->Dip_counts[g_id].begin();
             kt != ale->Dip_counts[g_id].end(); kt++) {
          pair<long int, long int> parts = (*kt).first;
          long int gp_id = parts.first;
          long int gpp_id = parts.second;
          int gp_i = g_id2i[parts.first];
          int gpp_i = g_id2i[parts.second];
          gp_is.push_back(gp_i);
          gpp_is.push_back(gpp_i);
          if (ale->Bip_counts[g_id] <= scalar_parameter["min_bip_count"])
            p_part.push_back(0);
          else
            p_part.push_back(
                pow((scalar_type)ale->p_dip(g_id, gp_id, gpp_id),
                    (scalar_type)scalar_parameter["seq_beta"])); // set pp
        }
      else {
        // XX
        // root bipartition needs to be handled separately
        map<set<long int>, int> bip_parts;
        for (map<long int, scalar_type>::iterator it = ale->Bip_counts.begin();
             it != ale->Bip_counts.end(); it++) {
          long int gp_id = (*it).first;
          boost::dynamic_bitset<> gamma = ale->id_sets.at(gp_id);
          boost::dynamic_bitset<> not_gamma = ~gamma;
          not_gamma[0] = 0;
          long int gpp_id = ale->set_ids.at(not_gamma);

          set<long int> parts;
          parts.insert(gp_id);
          parts.insert(gpp_id);
          bip_parts[parts] = 1;
          // gamma.clear();
          // not_gamma.clear();
        }
        for (map<set<long int>, int>::iterator kt = bip_parts.begin();
             kt != bip_parts.end(); kt++) {
          vector<long int> parts;
          for (set<long int>::iterator sit = (*kt).first.begin();
               sit != (*kt).first.end(); sit++) {
            parts.push_back((*sit));
          }
          long int gp_id = parts[0];
          // long int gpp_id=parts[1];

          int gp_i = g_id2i[parts[0]];
          int gpp_i = g_id2i[parts[1]];
          gp_is.push_back(gp_i);
          gpp_is.push_back(gpp_i);

          // Here we can create a new ale->Bip_counts[gp_id], in particular for
          // leaves. We may want to add the leaf entries for Bip_counts when
          // Bip_counts is first created.
          if (ale->Bip_counts[gp_id] <= scalar_parameter.at("min_bip_count") and
              not ale->Gamma_size < 4)
            p_part.push_back(0);
          else
            p_part.push_back(
                pow((scalar_type)ale->p_bip(gp_id),
                    (scalar_type)scalar_parameter["seq_beta"])); // set pp
        }
        bip_parts.clear();
      }
      // ######################################################################################################################
      // #########################################INNNER
      // LOOP##################################################################
      // ######################################################################################################################

      for (int e = 0; e < last_branch; e++) {
        scalar_type uq_sum = 0;
        // S leaf and G leaf
        if (e < last_leaf and is_a_leaf and
            extant_species[e] == gid_sps[g_id]) {
          // present
          uq_sum += PS[e] * 1;
        } else {
          // G internal
          if (not is_a_leaf) {
            int N_parts = gp_is.size();
            for (int i = 0; i < N_parts; i++) {
              int gp_i = gp_is[i];
              int gpp_i = gpp_is[i];
              scalar_type pp = p_part[i];
              if (not(e < last_leaf)) {
                int f = daughter[e];
                int g = son[e];
                // S event
                uq_sum +=
                    PS[e] *
                    (uq[gp_i][f] * uq[gpp_i][g] + uq[gp_i][g] * uq[gpp_i][f]) *
                    pp;
              }
              // D event
              uq_sum += PD[e] * (uq[gp_i][e] * uq[gpp_i][e]) *
                        pp; // no factor of two needed here
              // T event
              if (not no_T)
                uq_sum +=
                    (uq[gp_i][e] *
                         (mPTuq[gpp_i] - mPTuq_ancestral_correction[gpp_i][e]) /
                         tau_norm[e] +
                     uq[gpp_i][e] *
                         (mPTuq[gp_i] - mPTuq_ancestral_correction[gp_i][e]) /
                         tau_norm[e]) *
                    pp;
            }
          }
          if (not(e < last_leaf)) {
            int f = daughter[e];
            int g = son[e];
            // SL event
            uq_sum += PS[e] * (uq[i][f] * uE[g] + uq[i][g] * uE[f]);
          }
          // DL event
          uq_sum += PD[e] * (uq[i][e] * uE[e] * 2);
          // TL event
          if (not no_T)
            uq_sum += ((mPTuq[i] - mPTuq_ancestral_correction[i][e]) /
                           tau_norm[e] * uE[e] +
                       uq[i][e] * (mPTE - mPTE_ancestral_correction[e]) /
                           tau_norm[e]);
        }
        if (uq_sum < EPSILON)
          uq_sum = EPSILON;
        uq[i][e] = uq_sum;
        new_mPTuq += (wT[e]) * uq_sum;
      }
      for (int e = 0; e < last_branch; e++) {
        mPTuq_ancestral_correction[i][e] = 0;
        for (vector<int>::iterator it = ancestors[e].begin();
             it != ancestors[e].end(); it++) {
          int f = (*it);
          mPTuq_ancestral_correction[i][e] += (wT[f]) * uq[i][f];
        }
      }
      mPTuq[i] = new_mPTuq;

      // ######################################################################################################################
      // #########################################INNNER
      // LOOP##################################################################
      // ######################################################################################################################
    }
    survive = 0;
    root_sum = 0;
    bool single_O = false;
    for (int e = 0; e < last_branch; e++)
      if (vector_parameter["rate_multiplier_O"][e] < 0)
        single_O = true;
    if (single_O) {
      for (int e = 0; e < last_branch; e++) {
        if (vector_parameter["rate_multiplier_O"][e] > 0) {
          vector_parameter["rate_multiplier_O"][e] = 0;
        } else {
          vector_parameter["rate_multiplier_O"][e] = 1;
        }
      }
    }

    for (int e = 0; e < last_branch; e++) {
      scalar_type O_p = vector_parameter["rate_multiplier_O"][e];
      if (e == (last_branch - 1) and O_p == 1)
        O_p = scalar_parameter["O_R"];
      root_sum += uq[root_i][e] * O_p;
      survive += O_p * (1 - uE[e]);
    }
    for (int e = 0; e < last_branch; e++) {
      scalar_type O_p = vector_parameter["rate_multiplier_O"][e];
      if (e == (last_branch - 1) and O_p == 1)
        O_p = scalar_parameter["O_R"];
      branch_counts["O_LL"].at(e) = log(uq[root_i][e]) + log(O_p);
    }
    // cout << root_sum/survive << endl;
  }
  return root_sum * last_branch / survive;
}

string exODT_model::sample_undated(bool no_T) {

  scalar_type r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1);

  scalar_type root_sum = 0;
  scalar_type O_norm = 0;
  for (int e = 0; e < last_branch; e++) {
    branch_counts["saw"].at(e) = 0;
    scalar_type O_p = vector_parameter["rate_multiplier_O"][e];
    if (e == (last_branch - 1) and O_p == 1)
      O_p = scalar_parameter["O_R"];
    O_norm += O_p;
    root_sum += uq[g_ids.size() - 1][e] * O_p + EPSILON;
  }
  scalar_type root_resum = 0;
  for (int e = 0; e < last_branch; e++) {
    scalar_type O_p = vector_parameter["rate_multiplier_O"][e];
    if (e == (last_branch - 1) and O_p == 1)
      O_p = scalar_parameter["O_R"];
    root_resum += uq[root_i][e] * O_p + EPSILON;
    if (r * root_sum < root_resum) {
      register_O(e);
      return sample_undated(e, root_i, "O", "", no_T) + ";";
    }
  }
  return "-!=-";
}

string exODT_model::sample_undated(int e, int i, string last_event,
                                   string branch_string, bool no_T) {

  scalar_type r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1);

  bool is_a_leaf = false;
  long int g_id = g_ids[i];
  if (g_id_sizes[i] == 1)
    is_a_leaf = true;

  stringstream bl;
  if (ale_pointer->Bip_counts.count(g_id) and ale_pointer->Bip_counts[g_id] > 0)
    bl << max(ale_pointer->Bip_bls[g_id] / ale_pointer->Bip_counts[g_id],
              (scalar_type)scalar_parameter["min_branch_lenghts"]);
  else
    bl << max(ale_pointer->Bip_bls[g_id] / ale_pointer->observations,
              (scalar_type)scalar_parameter["min_branch_lenghts"]);
  string branch_length = bl.str();

  vector<int> gp_is;
  vector<long int> gpp_is;
  vector<scalar_type> p_part;
  if (g_id != -1)
    for (unordered_map<pair<long int, long int>, scalar_type>::iterator kt =
             ale_pointer->Dip_counts[g_id].begin();
         kt != ale_pointer->Dip_counts[g_id].end(); kt++) {
      pair<long int, long int> parts = (*kt).first;
      long int gp_id = parts.first;
      long int gpp_id = parts.second;
      int gp_i = g_id2i[parts.first];
      int gpp_i = g_id2i[parts.second];
      gp_is.push_back(gp_i);
      gpp_is.push_back(gpp_i);
      if (ale_pointer->Bip_counts[g_id] <= scalar_parameter["min_bip_count"])
        p_part.push_back(0);
      else
        p_part.push_back(
            pow((scalar_type)ale_pointer->p_dip(g_id, gp_id, gpp_id),
                (scalar_type)scalar_parameter["seq_beta"])); // set pp
    }
  else {
    // root bipartition needs to be handled separately
    map<set<long int>, int> bip_parts;
    for (map<long int, scalar_type>::iterator it =
             ale_pointer->Bip_counts.begin();
         it != ale_pointer->Bip_counts.end(); it++) {
      long int gp_id = (*it).first;
      boost::dynamic_bitset<> gamma = ale_pointer->id_sets.at(gp_id);
      boost::dynamic_bitset<> not_gamma = ~gamma;
      not_gamma[0] = 0;
      long int gpp_id = ale_pointer->set_ids.at(not_gamma);

      set<long int> parts;
      parts.insert(gp_id);
      parts.insert(gpp_id);
      bip_parts[parts] = 1;
      // gamma.clear();
      // not_gamma.clear();
    }
    for (map<set<long int>, int>::iterator kt = bip_parts.begin();
         kt != bip_parts.end(); kt++) {
      vector<long int> parts;
      for (set<long int>::iterator sit = (*kt).first.begin();
           sit != (*kt).first.end(); sit++) {
        parts.push_back((*sit));
      }
      long int gp_id = parts[0];
      // long int gpp_id=parts[1];

      int gp_i = g_id2i[parts[0]];
      int gpp_i = g_id2i[parts[1]];
      gp_is.push_back(gp_i);
      gpp_is.push_back(gpp_i);

      // Here we can create a new ale->Bip_counts[gp_id], in particular for
      // leaves. We may want to add the leaf entries for Bip_counts when
      // Bip_counts is first created.
      if (ale_pointer->Bip_counts[gp_id] <=
              scalar_parameter.at("min_bip_count") and
          not ale_pointer->Gamma_size < 4)
        p_part.push_back(0);
      else
        p_part.push_back(
            pow((scalar_type)ale_pointer->p_bip(gp_id),
                (scalar_type)scalar_parameter["seq_beta"])); // set pp
    }
    bip_parts.clear();
  }

  // ######################################################################################################################
  // #########################################INNNER
  // LOOP##################################################################
  // ######################################################################################################################
  scalar_type uq_sum = 0;
  // S leaf and G leaf
  if (e < last_leaf and is_a_leaf and extant_species[e] == gid_sps[g_id]) {
    // present
    uq_sum += PS[e] * 1 + EPSILON;
  } else {

    // G internal
    if (not is_a_leaf) {
      int N_parts = gp_is.size();
      for (int i = 0; i < N_parts; i++) {
        int gp_i = gp_is[i];
        int gpp_i = gpp_is[i];
        scalar_type pp = p_part[i];
        if (not(e < last_leaf)) {
          int f = daughter[e];
          int g = son[e];
          // S event
          uq_sum += PS[e] * uq[gp_i][f] * uq[gpp_i][g] * pp + EPSILON;
          uq_sum += PS[e] * uq[gp_i][g] * uq[gpp_i][f] * pp + EPSILON;
        }
        // D event
        uq_sum += PD[e] * (uq[gp_i][e] * uq[gpp_i][e] * 2) * pp + EPSILON;
        // T event
        for (int f = 0; f < last_branch; f++)
          if (not ancestral[e][f] and not no_T) {
            uq_sum += uq[gp_i][e] * (wT[f] / tau_norm[e]) * uq[gpp_i][f] * pp +
                      EPSILON;
            uq_sum += uq[gpp_i][e] * (wT[f] / tau_norm[e]) * uq[gp_i][f] * pp +
                      EPSILON;
          }
      }
    }

    if (not(e < last_leaf)) {
      int f = daughter[e];
      int g = son[e];
      // SL event
      uq_sum += PS[e] * uq[i][f] * uE[g] + EPSILON;
      uq_sum += PS[e] * uq[i][g] * uE[f] + EPSILON;
    }

    // DL event
    uq_sum += PD[e] * (uq[i][e] * uE[e] * 2) + EPSILON;
    // TL event
    for (int f = 0; f < last_branch; f++)
      if (not ancestral[e][f] and not no_T) {
        uq_sum += (wT[f] / tau_norm[e]) * uq[i][f] * uE[e] + EPSILON;
        uq_sum += (wT[f] / tau_norm[e]) * uE[f] * uq[i][e] + EPSILON;
      }
  }

  // ######################################################################################################################
  // #########################################INNNER
  // LOOP##################################################################
  // ######################################################################################################################

  // ######################################################################################################################
  // #########################################INNNER
  // LOOP##################################################################
  // ######################################################################################################################

  stringstream estring;
  if (not(e < last_leaf))
    estring << e;
  else
    estring << extant_species[e];
  string estr = estring.str();

  scalar_type uq_resum = 0;
  // S leaf and G leaf
  if (e < last_leaf and is_a_leaf and extant_species[e] == gid_sps[g_id]) {
    // present
    uq_resum += PS[e] * 1 + EPSILON;
    if (r * uq_sum < uq_resum) {
      register_leafu(e, last_event);
      return ale_pointer->set2name(ale_pointer->id_sets[g_id]) + branch_string +
             ":" + branch_length;
    }
  } else {
    // G internal
    if (not is_a_leaf) {
      int N_parts = gp_is.size();
      for (int i = 0; i < N_parts; i++) {
        int gp_i = gp_is[i];
        int gpp_i = gpp_is[i];
        scalar_type pp = p_part[i];
        if (not(e < last_leaf)) {
          int f = daughter[e];
          int g = son[e];
          // S event
          uq_resum += PS[e] * uq[gp_i][f] * uq[gpp_i][g] * pp + EPSILON;
          if (r * uq_sum < uq_resum) {
            register_Su(e, last_event);
            return "(" + sample_undated(f, gp_i, "S", "", no_T) + "," +
                   sample_undated(g, gpp_i, "S", "", no_T) + ")." + estr +
                   branch_string + ":" + branch_length;
          }
          uq_resum += PS[e] * uq[gp_i][g] * uq[gpp_i][f] * pp + EPSILON;
          if (r * uq_sum < uq_resum) {
            register_Su(e, last_event);
            return "(" + sample_undated(g, gp_i, "S", "", no_T) + "," +
                   sample_undated(f, gpp_i, "S", "", no_T) + ")." + estr +
                   branch_string + ":" + branch_length;
          }
        }
        // D event
        uq_resum += PD[e] * (uq[gp_i][e] * uq[gpp_i][e] * 2) * pp + EPSILON;
        if (r * uq_sum < uq_resum or no_T) {
          register_D(e);
          return "(" + sample_undated(e, gp_i, "D", "", no_T) + "," +
                 sample_undated(e, gpp_i, "D", "", no_T) + ").D@" + estr +
                 branch_string + ":" + branch_length;
        }

        // T event
        for (int f = 0; f < last_branch; f++)
          if (not ancestral[e][f] and not no_T) {
            stringstream fstring;
            if (not(f < last_leaf))
              fstring << f;
            else
              fstring << extant_species[f];
            string fstr = fstring.str();

            uq_resum +=
                uq[gp_i][e] * (wT[f] / tau_norm[e]) * uq[gpp_i][f] * pp +
                EPSILON;
            if (r * uq_sum < uq_resum) {
              register_Tfrom(e);
              register_Tto(f);
              register_T_to_from(e, f);
              stringstream Ttoken;
              Ttoken << estr << ">" << fstr << "|"
                     << ale_pointer->set2name(
                            ale_pointer->id_sets[g_ids[gpp_i]]);
              Ttokens.push_back(Ttoken.str());

              return "(" + sample_undated(e, gp_i, "S", "", no_T) + "," +
                     sample_undated(f, gpp_i, "T", "", no_T) + ").T@" + estr +
                     "->" + fstr + branch_string + ":" + branch_length;
            }
            uq_resum +=
                uq[gpp_i][e] * (wT[f] / tau_norm[e]) * uq[gp_i][f] * pp +
                EPSILON;
            if (r * uq_sum < uq_resum) {
              register_Tfrom(e);
              register_Tto(f);
              register_T_to_from(e, f);
              stringstream Ttoken;
              Ttoken << estr << ">" << fstr << "|"
                     << ale_pointer->set2name(
                            ale_pointer->id_sets[g_ids[gp_i]]);
              Ttokens.push_back(Ttoken.str());
              return "(" + sample_undated(e, gpp_i, "S", "", no_T) + "," +
                     sample_undated(f, gp_i, "T", "", no_T) + ").T@" + estr +
                     "->" + fstr + branch_string + ":" + branch_length;
            }
          }
      }
    }
    if (not(e < last_leaf)) {
      int f = daughter[e];
      int g = son[e];
      // SL event
      uq_resum += PS[e] * uq[i][f] * uE[g] + EPSILON;
      if (r * uq_sum < uq_resum) {
        register_Su(e, last_event);
        register_L(g);
        return sample_undated(f, i, "S", "." + estr + branch_string, no_T);
      }
      uq_resum += PS[e] * uq[i][g] * uE[f] + EPSILON;
      if (r * uq_sum < uq_resum) {
        register_Su(e, last_event);
        register_L(f);
        return sample_undated(g, i, "S", "." + estr + branch_string, no_T);
      }
    }
    // DL event
    uq_resum += PD[e] * (uq[i][e] * uE[e] * 2) + EPSILON;
    if (r * uq_sum < uq_resum) {
      return sample_undated(e, i, "S", branch_string, no_T);
    }
    // TL event
    for (int f = 0; f < last_branch; f++)
      if (not ancestral[e][f] and not no_T) {
        stringstream fstring;
        if (not(f < last_leaf))
          fstring << f;
        else
          fstring << extant_species[f];
        string fstr = fstring.str();

        uq_resum += (wT[f] / tau_norm[e]) * uq[i][f] * uE[e] + EPSILON;
        if (r * uq_sum < uq_resum) {
          register_Tfrom(e);
          register_Tto(f);
          register_T_to_from(e, f);
          /*
          stringstream Ttoken;
          Ttoken<<estr<<">"<<fstr<<"|"<<ale_pointer->set2name(ale_pointer->id_sets[g_id]);
          Ttokens.push_back(Ttoken.str());
          */
          register_L(e);
          return sample_undated(
              f, i, "T", ".T@" + estr + "->" + fstr + branch_string, no_T);
        }
        uq_resum += (wT[f] / tau_norm[e]) * uE[f] * uq[i][e] + EPSILON;
        if (r * uq_sum < uq_resum) {
          return sample_undated(e, i, "S", "", no_T);
        }
      }
  }
  // ######################################################################################################################
  // #########################################INNNER
  // LOOP##################################################################
  // ######################################################################################################################
  cout << "sum error!" << endl;
  return "-!=-";
}

string exODT_model::counts_string_undated(scalar_type samples) {

  stringstream out;
  for (int e = 0; e < last_branch; e++) {
    bool isleaf = e < last_leaf;
    stringstream named_branch;
    if (e < last_leaf)
      named_branch << extant_species[e] << "(" << e << ")";
    else
      named_branch << e;

    if (not isleaf)
      out << "S_internal_branch\t" << named_branch.str() << "\t"
          << branch_counts["Ds"][e] / samples << "\t"
          << branch_counts["Ts"][e] / samples << "\t"
          << branch_counts["Ls"][e] / samples << "\t"
          << branch_counts["Os"][e] / samples << "\t"
          << branch_counts["copies"][e] / samples << "\t"
          << branch_counts["singleton"][e] / samples << "\t" << uE[e] << "\t"
          << branch_counts["presence"][e] / samples << "\t"
          << branch_counts["O_LL"][e] << "\n";

    else
      out << "S_terminal_branch\t" << named_branch.str() << "\t"
          << branch_counts["Ds"][e] / samples << "\t"
          << branch_counts["Ts"][e] / samples << "\t"
          << branch_counts["Ls"][e] / samples << "\t"
          << branch_counts["Os"][e] / samples << "\t"
          << branch_counts["copies"][e] / samples << "\t"
          << branch_counts["singleton"][e] / samples << "\t" << uE[e] << "\t"
          << branch_counts["presence"][e] / samples << "\t"
          << branch_counts["O_LL"][e] << "\n";
  }
  return out.str();
}

void exODT_model::register_Su(int e, string last_event) {
  MLRec_events["S"] += 1;
  if (e > -1) {
    int f = daughter[e];
    int g = son[e];
    if (last_event == "S" or last_event == "O")
      branch_counts["singleton"].at(e) += 1;
    branch_counts["copies"].at(e) += 1;
    if (branch_counts["saw"].at(e) == 0)
      branch_counts["presence"].at(e) += 1;
    branch_counts["saw"].at(e) = 1;

    branch_counts["count"].at(f) += 1;
    branch_counts["count"].at(g) += 1;
  }
}

void exODT_model::register_leafu(int e, string last_event) {
  if (e > -1) {
    branch_counts["copies"].at(e) += 1;
    if (branch_counts["saw"].at(e) == 0)
      branch_counts["presence"].at(e) += 1;
    branch_counts["saw"].at(e) = 1;
    if (last_event == "S" or last_event == "O")
      branch_counts["singleton"].at(e) += 1;
  }
  // MLRec_events["genes"]+=1;
}

void exODT_model::register_T_to_from(int e, int f) { T_to_from[e][f] += 1; }

void exODT_model::reset_T_to_from() {
  for (int e = 0; e < last_branch; e++) {
    for (int f = 0; f < last_branch; f++) {
      T_to_from[e][f] = 0;
    }
  }
}

string exODT_model::feSPR(int e, int f) {
  tree_type *newS = TreeTemplateTools::parenthesisToTree(
      string_parameter["S_un"],
      (string_parameter["BOOTSTRAP_LABELS"] == "yes"));
  Node *newS_root = newS->getRootNode();
  vector<Node *> nodes = TreeTemplateTools::getNodes(*newS_root);

  string e_name = node_name[id_nodes[e]];
  string f_name = node_name[id_nodes[f]];
  ;

  Node *e_node, *f_node;

  for (vector<Node *>::iterator it = nodes.begin(); it != nodes.end(); it++) {
    string name_it;
    if ((*it)->isLeaf()) {
      name_it = (*it)->getName();
    } else {
      vector<string> leafnames = TreeTemplateTools::getLeavesNames(*(*it));
      sort(leafnames.begin(), leafnames.end());
      stringstream name;
      for (vector<string>::iterator st = leafnames.begin();
           st != leafnames.end(); st++)
        name << (*st) << ".";

      name_it = name.str();
    }
    if (name_it == e_name)
      e_node = (*it);
    if (name_it == f_name)
      f_node = (*it);
  }

  if (e == f)
    return string_parameter["S_un"];

  bool e_below_f = false;
  Node *node;
  node = e_node;
  while (node->hasFather()) {
    node = node->getFather();
    if (node == f_node)
      e_below_f = true;
  }
  if (e_below_f) {
    Node *swap_tmp = e_node;
    e_node = f_node;
    f_node = swap_tmp;
  }
  if (f_node->hasFather() and f_node->getFather() == e_node)
    return string_parameter["S_un"];

  Node *f_father = f_node->getFather();
  vector<Node *> f_sons = f_father->getSons();
  Node *f_sister;
  if (f_sons[0] == f_node)
    f_sister = f_sons[1];
  else
    f_sister = f_sons[0];
  f_father->removeSon(f_sister);
  if (f_father->hasFather()) {
    Node *f_grand_father = f_father->getFather();
    f_grand_father->removeSon(f_father);
    f_grand_father->addSon(f_sister);
  } else {
    newS->setRootNode(f_sister);
  }
  if (e_node->hasFather()) {
    Node *e_father = e_node->getFather();
    e_father->removeSon(e_node);
    e_father->addSon(f_father);
  } else
    newS->setRootNode(f_father);

  f_father->addSon(e_node);

  // for (vector <Node * >::iterator it=nodes.begin();it!=nodes.end();it++ )
  // (*it)->setDistanceToFather(1);
  //?

  return TreeTemplateTools::treeToParenthesis(*newS, false, "ID");
}

vector<string> exODT_model::NNIs(int e) {
  vector<string> NNIs;
  int left_e, right_e, f;

  Node *root = id_nodes[e];

  if (root->isLeaf())
    return NNIs;

  vector<Node *> roots_sons = root->getSons();

  right_e = node_ids[roots_sons[0]];
  left_e = node_ids[roots_sons[1]];

  if (roots_sons[0]->isLeaf())
    ;
  else {
    vector<Node *> right_sons = roots_sons[0]->getSons();
    f = node_ids[right_sons[0]];
    NNIs.push_back(feSPR(left_e, f));
    f = node_ids[right_sons[1]];
    NNIs.push_back(feSPR(left_e, f));
  }

  if (roots_sons[1]->isLeaf())
    ;
  else {
    vector<Node *> left_sons = roots_sons[1]->getSons();
    f = node_ids[left_sons[0]];
    NNIs.push_back(feSPR(right_e, f));
    f = node_ids[left_sons[1]];
    NNIs.push_back(feSPR(right_e, f));
  }
  return NNIs;
}
