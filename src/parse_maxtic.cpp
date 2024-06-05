#include "ALE.h"

using namespace std;
using namespace bpp;

string azname(Node *node) {
  vector<string> left_aznames =
      TreeTemplateTools::getLeavesNames(*(node->getSons()[0]));
  vector<string> right_aznames =
      TreeTemplateTools::getLeavesNames(*(node->getSons()[1]));

  sort(left_aznames.begin(), left_aznames.end());
  sort(right_aznames.begin(), right_aznames.end());

  string azname =
      left_aznames[0] + "-" + right_aznames[right_aznames.size() - 1];
  return azname;
}

int main(int argc, char **argv) {

  ifstream tree_stream(argv[1]);
  string fname = argv[1];

  string tree;
  getline(tree_stream, tree);
  tree_type *T = TreeTemplateTools::parenthesisToTree(tree, true);

  vector<Node *> nodes = T->getNodes();

  for (auto it = nodes.begin(); it != nodes.end(); it++)
    if ((*it)->hasFather() and not(*it)->isLeaf()) {
      scalar_type h = TreeTemplateTools::getHeight(*(*it));
      Node *node = (*it);
      string name = azname(node);
      int id = node->getBootstrapValue();
      cout << id << " " << name << endl;
    }
}
