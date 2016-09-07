#ifndef TWOLAYERGRAPH_H_
#define TWOLAYERGRAPH_H_

/*!
 \file   TwoLayerGraph.h
 \brief  This file provides an undirected graph between regions in an image
 \author Frank Moosmann (<frank.moosmann@kit.edu>),
 \date   2011

         Copyright: Institute for Measurement and Control Systems,
         Karlsruhe Institute of Technology.
         http://www.mrt.uni-karlsruhe.de
         All rights reserved
*/

#include <string>
#include <iostream>
#include <fstream>
#include <deque>
#include <vector>
#include <map>
#include <boost/tuple/tuple.hpp>
#include <boost/foreach.hpp>
#include <boost/function.hpp>

/*!
 \class TwoLayerGraph
 \brief

 */
template<typename NodeT1, typename NodeT2>
class TwoLayerGraph {
public:
  typedef unsigned int ValT; //!< internal representation of the link-count.
  typedef boost::function<void (std::ofstream &of)> KeySpecifier;
  typedef boost::function<void (const NodeT1 &node1, std::ofstream &of)> Node1Visitor;
  typedef boost::function<void (const NodeT2 &node2, std::ofstream &of)> Node2Visitor;
  typedef boost::function<void (const NodeT1 &node1, const NodeT2 &node2, ValT val, std::ofstream &of)> ConnectionVisitor;
private:
  // typedefs and data
  typedef size_t IdT; //!< internal representation of an index
  typedef std::pair<IdT, IdT> ConnectionT; //!< internal representation of a connection index
  //IdT                   _nextID;
  std::deque<NodeT1> _nodes1; //!< Storage for all 1st-level-Nodes of the graph.
  std::deque<NodeT2> _nodes2; //!< Storage for all 2nd-level-Nodes of the graph.
  std::vector< std::vector<ValT> > _connections; //!< Storage for all connections of the graph in a matrix
private:
  // functions
  IdT _getNode1(NodeT1 node1) const {
    IdT idx = 0;
    while ((idx < _nodes1.size()) && (_nodes1[idx] != node1))
      ++idx;
    return idx;
  };
  IdT _addNode1(NodeT1 node1) {
    IdT idx = _getNode1(node1);
    if (idx == _nodes1.size()) {
      _nodes1.push_back(node1);
      _connections.push_back(std::vector<ValT>( _nodes2.size() ));
    }
    return idx;
  };
  IdT _getNode2(NodeT2 node2) const {
    IdT idx = 0;
    while ((idx < _nodes2.size()) && (_nodes2[idx] != node2))
      ++idx;
    return idx;
  };
  IdT _addNode2(NodeT2 node2) {
    IdT idx = _getNode2(node2);
    if (idx == _nodes2.size()) {
      _nodes2.push_back(node2);
      BOOST_FOREACH(std::vector<ValT> &v, _connections)
        v.resize(_nodes2.size());
    }
    return idx;
  };

public:
  //! constructs an empty graph
  TwoLayerGraph() { clear(); };
  //! destructor
  virtual ~TwoLayerGraph() {};

  //! returns the number of nodes in 1st layer of the graph
  size_t size1() const { return _nodes1.size(); };
  //! returns the number of nodes in 2nd layer of the graph
  size_t size2() const { return _nodes2.size(); };
  //! reset graph to a new one
  void clear() { _nodes1.clear(); _nodes2.clear(); _connections.clear(); };
  //! tells if the graph is empty or not
  bool empty() const { return (_nodes1.empty()) && (_nodes2.empty()); };
  //! sets all 1st-layer nodes, makes adding connections more efficient?
  template <typename InputIteratorN1, typename InputIteratorN2>
  void setNodes(InputIteratorN1 begin1, InputIteratorN1 end1, InputIteratorN2 begin2, InputIteratorN2 end2) {
    clear();
    while (begin1 != end1)
      _nodes1.push_back(*begin1++);
    while (begin2 != end2)
      _nodes2.push_back(*begin2++);
    for (unsigned int i=0; i<_nodes1.size(); ++i)
      _connections.push_back(std::vector<ValT>(_nodes2.size()));
  };
  //! adds a 1st-level-node without connection
  void addNode1(NodeT1 node1) {
    _addNode1(node1);
  };
  //! adds a 2nd-level-node without connection
  void addNode2(NodeT2 node2) {
    _addNode2(node2);
  };
  //! adds a connection
  void addConnection(NodeT1 node1, NodeT2 node2) {
    IdT i1 = _addNode1(node1);
    IdT i2 = _addNode2(node2);
    _connections[i1][i2]++;
  };
  //! returns all 1st-level-nodes
  std::vector< NodeT1 > getNodes1() const {
    std::vector< NodeT1 > retVec;
    retVec.reserve(_nodes1.size());
    BOOST_FOREACH(NodeT1 &n, _nodes1)
      retVec.push_back(n);
    return retVec;
  };
  //! returns all 2nd-level-nodes
  std::vector< NodeT2 > getNodes2() const {
    std::vector< NodeT2 > retVec;
    retVec.reserve(_nodes2.size());
    BOOST_FOREACH(NodeT2 &n, _nodes2)
      retVec.push_back(n);
    return retVec;
  };
  //! returns all connections of specified node
  std::vector< std::pair<NodeT2, ValT> > getConnectionsToL2(NodeT1 node1) const {
    std::vector< std::pair<NodeT2, ValT> > retVec;
    IdT i1 = _getNode1(node1);
    if (i1 == _nodes1.size()) return retVec; // node doesn't exist
    for (IdT i2 = 0; i2 < _nodes2.size(); ++i2) {
      if (_connections[i1][i2] > 0)
        retVec.push_back(std::make_pair(_nodes2[i2], _connections[i1][i2]));
    }
    return retVec;
  };
  //! returns matrix indicating linkage: mat[l1][l2] = linkage count
  std::vector< std::vector<ValT> > getAssociationMatrix() const {
    return _connections;
  };
  //!< save graph, preferred ending: .graphml
  void saveToGraphML(std::string filename, KeySpecifier keyvisit = KeySpecifier(), Node1Visitor n1visit = Node1Visitor(), Node2Visitor n2visit = Node2Visitor(), ConnectionVisitor cvisit = ConnectionVisitor()) const {
    using namespace std;
    ofstream of(filename.c_str(), ios::out);
    of << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>" << endl;
    of << "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"" << endl;
    of << "         xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"" << endl;
    of << "         xmlns:y=\"http://www.yworks.com/xml/graphml\" xmlns:yed=\"http://www.yworks.com/xml/yed/3\"" << endl;
    of << "         xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns http://www.yworks.com/xml/schema/graphml/1.1/ygraphml.xsd\">" << endl;

    // 1) specifiy keys
    if (keyvisit) keyvisit(of);
    of << "  <graph edgedefault=\"undirected\" id=\"G\">" << endl;
    // 2) output nodes
    for (unsigned int i1 = 0; i1 < _nodes1.size(); ++i1) {
      of << "    <node id=\"n1_" << i1 << "\">" << endl;
      if (n1visit) n1visit(_nodes1[i1], of);
      of << "    </node>" << endl;
    }
    for (unsigned int i2 = 0; i2 < _nodes2.size(); ++i2) {
      of << "    <node id=\"n2_" << i2 << "\">" << endl;
      if (n2visit) n2visit(_nodes2[i2], of);
      of << "    </node>" << endl;
    }
    // 3) output connections
    for (unsigned int i1 = 0; i1 < _nodes1.size(); ++i1) {
      for (unsigned int i2 = 0; i2 < _nodes2.size(); ++i2) {
        if (_connections[i1][i2] > 0) {
          of << "    <edge source=\"n1_" << i1 << "\" target=\"n2_" << i2 << "\">" << endl;
          if (cvisit) cvisit(_nodes1[i1], _nodes2[i2], _connections[i1][i2], of);
          of << "    </edge>" << endl;
        }
      }
    }

    of << "  </graph>" << endl;
    of << "</graphml>" << endl;
  };

  //test functions
//  static void test(); //!< synthetic data to test correct graph building process


};


#endif /* TWOLAYERGRAPH_H_ */
