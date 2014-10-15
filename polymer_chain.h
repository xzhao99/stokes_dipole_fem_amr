// The libMesh Finite Element Library.
// Copyright (C) 2014  Xujun Zhao


//#ifndef POLYMER_CHAIN_H
//#define POLYMER_CHAIN_H


// C++ Includes
#include <iostream>
#include <algorithm>
#include <vector>
#include <math.h>

// LibMesh library includes
#include "libmesh/libmesh_common.h"
#include "libmesh/mesh.h"
#include "libmesh/elem.h"
//#include "libmesh/numeric_vector.h"   // difference between this and next?
//#include "libmesh/vector_value.h"
#include "libmesh/id_types.h"
#include "libmesh/point.h"

// include user defined classes
#include "GeomTools.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;



/**
 * This class defines a polymer chain class using bead-spring model
 * which is defined in the LIBMESH_DIM dimensional Real space
 */

class PolymerChain
{
public:
  // constructor
  PolymerChain();
  
  // destructor
  ~PolymerChain();
  
  // initialize the particle-mesh relation
  void init(const MeshBase& mesh);
  
  // return total number of beads in the polymer chain
  unsigned int num_beads() const;
  
  // return total number of spring in the polymer chain
  unsigned int num_springs() const;
  
  // check if this element contains particles
  bool contain_particle(const dof_id_type elem_id) const;
  
  // return the spacial coordinates of a particle
  Point particle_coordinate(dof_id_type particle_id) const;
  
  // return spring directions, pointing from bead A to B
  Point spring_direction(dof_id_type spring_id) const;
  
  // return bead forces, including magnitude and direction for each bead
  // this is a temporary function for dipoles with 2 particles!
  Point bead_force(dof_id_type spring_id) const;

  
  // if this element contains particles, return the particle id's
  // if it does NOT contain particles, the size of vector returned is 0!
  std::vector<dof_id_type> element_particle_map(const dof_id_type elem_id) const;
  
  
  // print out the information
  void print_info() const;
  
private:
  // initialize the polymer chain
  void init();
  
  // spatial dimension = LIBMISH_DIM
  //unsigned int _dim;
  
  // total number of beads
  unsigned int _num_beads;
  
  // total number of springs
  unsigned int _num_springs;
  
  // spatial coordinates of beads
  std::vector<Point> _bead_coordinates;
  
  // bead connection: id of two springs connected to a bead,
  // only one connected spring for the beads at two ends
  std::vector<std::vector<dof_id_type> > _bead_connection;
  
  // spring connection: id of two beads(A-B) to form a spring
  std::vector<std::vector<dof_id_type> > _spring_connection;
    
  // direction of spring connected by two beads. pointing from bead A to B
  std::vector<Point> _spring_direction;
  
  // bead forces, including magnitude and direction for each bead
  std::vector<Point> _bead_force;
  
  
  // map vector from finite element to particle
  // one element can contains multiple particles!
  // if an element doesn't contain particles, _map[ielm].size()==0
  std::vector<std::vector<dof_id_type> > _element_particle_map;
  
  // map vector from particle to element.
  // Assume a particle can only be in one element, even if it is on the sides or vertex
  std::vector<dof_id_type> _particle_element_map;
  
  // flag if an element contains particles
  std::vector<bool> _element_has_particle;
  
};


/// ============================================================================= ///
PolymerChain::PolymerChain()
{
  // do nothing
  init();
}

/// ============================================================================= ///
PolymerChain::~PolymerChain()
{
  // do nothing
}


/// ============================================================================= ///
void PolymerChain::init(const MeshBase& mesh)
{
  // initialize variables
  const dof_id_type num_elem = mesh.n_elem();
  _element_has_particle.resize(num_elem);
  _element_particle_map.resize(num_elem);
  _particle_element_map.resize(_num_beads);
  
  // loop over each finite element
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
  
  dof_id_type elem_id = 0;
  for ( ; el != end_el; ++el)
  {
    // Store a pointer to the element we are currently working on.
    const Elem* elem = *el;
    
    // loop over each particle (need N_particle*N_elem operations)
    for (std::size_t i=0; i<_num_beads; ++i)
    {
      const bool contain_point = elem->contains_point(_bead_coordinates[i]);
      if(contain_point)
      {
        //_element_particle_map[elem_id].push_back(i); => cause problems when bead sits on edges and nodes
        _particle_element_map[i] = elem_id;
      } // end if(contain_point)
    } // end of loop for(i)
    
    _element_has_particle[elem_id] = false; // set all false!
    
    ++elem_id;
  } // end of loop for(el)
  
  // loop over each bead to find the parent element in which they are
  for (std::size_t i=0; i<_num_beads; ++i)
  {
    elem_id = _particle_element_map[i];
    _element_particle_map[elem_id].push_back(i);
    _element_has_particle[elem_id] = true;
  }

}


/// ============================================================================= ///
unsigned int PolymerChain::num_beads() const
{
  return _num_beads;
}


/// ============================================================================= ///
unsigned int PolymerChain::num_springs() const
{
  return _num_springs;
}


/// ============================================================================= ///
bool PolymerChain::contain_particle(const dof_id_type elem_id) const
{
  return _element_has_particle[elem_id];
}


/// ============================================================================= ///
Point PolymerChain::particle_coordinate(dof_id_type p_id) const
{
  return _bead_coordinates[p_id];
}


/// ============================================================================= ///
Point PolymerChain::spring_direction(dof_id_type spring_id) const
{
  return _spring_direction[spring_id];
}


/// ============================================================================= ///
Point PolymerChain::bead_force(dof_id_type bead_id) const
{
  return _bead_force[bead_id];
}



/// ============================================================================= ///
std::vector<dof_id_type> PolymerChain::
  element_particle_map(const dof_id_type elem_id) const
{
  return _element_particle_map[elem_id];
}




/// ============================================================================= ///
void PolymerChain::print_info() const
{
  std::cout<<std::endl;
  std::cout<<"------------------------PolymerChain::print_info---------------------------"<<std::endl;
  std::cout<< "There are totally "<< _num_beads<<" beads, and "<< _num_springs
           << " springs in this chain." <<std::endl;
  
  // print particle-element map info
  for(std::size_t i=0; i<_num_beads; ++i)
  {
    std::cout<<"bead "<<i<<" is in the element "<<_particle_element_map[i]<<std::endl;
  } // end for i-loop
  
  // print element-particle map info
  for(std::size_t i=0; i<_element_particle_map.size(); ++i)
  {
    const size_t n = _element_particle_map[i].size();
    const bool has_particle = _element_has_particle[i];
    if(has_particle)  //(n>0)
    {
      //std::cout<<"element "<<i<<" contains particles:";
      std::cout<<"element "<<i<<" contains "<<n<<" particles:";
      
      for(std::size_t j=0; j<n; ++j)
        std::cout << _element_particle_map[i][j] ;
      std::cout << std::endl;
    } // end if (n>0)
    
  } // end for i-loop
  
  std::cout<<"---------------------------------------------------------------"<<std::endl;
  std::cout<<std::endl;
}


/// ============================================================================= ///
void PolymerChain::init()
{
  // simple case: 1 spring with two beads
  _num_springs = 1;
  _num_beads = 2;
  
  // initialize the particle position
  _bead_coordinates.resize(_num_beads);
  Point pt0, pt1;
  pt0(0) = -0.10; pt0(1) = -0.1;
  pt1(0) = +0.10; pt1(1) = +0.1;
  _bead_coordinates[0] = pt0;
  _bead_coordinates[1] = pt1;
  
  // initialize the spring direction
  _spring_direction.resize(_num_springs);
  Point pt_n = pt1 - pt0;
  Real pt_n_norm = GeomTools::point_norm(pt_n);
  _spring_direction[0] = pt_n/pt_n_norm;

  // bead force direction.
  _bead_force.resize(2);
  _bead_force[0] = _spring_direction[0];
  _bead_force[1] = -_spring_direction[0];
  
  // initialize the bead connections and spring connections of the chain.
  
  
}


/// ============================================================================= ///



/// ============================================================================= ///
