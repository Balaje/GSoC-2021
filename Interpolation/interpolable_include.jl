## Methods by @fverdugo
using Gridap
using Test
using Gridap.FESpaces
using Gridap.ReferenceFEs
using Gridap.Fields
using Gridap.CellData
using Gridap.Arrays
using Gridap.Geometry
using StaticArrays
using NearestNeighbors
using Gridap.Helpers

struct PushDofMap <: Map end

function Arrays.evaluate!(cache,::PushDofMap,f::Dof,m::Field)
  @abstractmethod
end

function Arrays.evaluate!(cache,::PushDofMap,f::AbstractArray{<:Dof},m::Field)
  @abstractmethod
end

# Local implementations
function Arrays.return_cache(::PushDofMap,f::LagrangianDofBasis,m::Field)
  q = f.nodes
  return_cache(m,q)
end
function replace_nodes(f::LagrangianDofBasis,x)
  LagrangianDofBasis(x, f.dof_to_node, f.dof_to_comp, f.node_and_comp_to_dof)
end
function Arrays.evaluate!(cache,::PushDofMap,f::LagrangianDofBasis,m::Field)
  q = f.nodes
  x = evaluate!(cache,m,q)
  LagrangianDofBasis(x, f.dof_to_node, f.dof_to_comp, f.node_and_comp_to_dof)
end

function Arrays.return_cache(::PushDofMap,f::MomentBasedDofBasis,m::Field)
  q = f.nodes
  return_cache(m,q)
end
function replace_nodes(f::MomentBasedDofBasis,x)
  MomentBasedDofBasis(x, f.face_moments, f.face_nodes)
end
function Arrays.evaluate!(cache,::PushDofMap,f::MomentBasedDofBasis,m::Field)
  q = f.nodes
  x = evaluate!(cache,m,q)
  MomentBasedDofBasis(x, f.face_moments, f.face_nodes)
end

# For performance reasons we also want these global implementations
# in practice, only the global implementations will be used

function Arrays.lazy_map(
  ::PushDofMap,
  cell_f::AbstractArray{<:LagrangianDofBasis},
  cell_m::AbstractArray{<:Field})

  cell_q = lazy_map(f->f.nodes,cell_f)
  cell_x = lazy_map(evaluate,cell_m,cell_q)
  lazy_map(LagrangianDofBasis,cell_f,cell_x)
end

function Arrays.lazy_map(
  ::PushDofMap,
  cell_f::AbstractArray{<:MomentBasedDofBasis},
  cell_m::AbstractArray{<:Field})

  cell_q = lazy_map(f->f.nodes,cell_f)
  cell_x = lazy_map(evaluate,cell_m,cell_q)
  lazy_map(replace_nodes,cell_f,cell_x)
end


# In CellData/CellDofs.jl

function CellData.change_domain(a::CellDof,::ReferenceDomain,::PhysicalDomain)
  trian = get_triangulation(a)
  cell_m = get_cell_map(trian)
  cell_f_ref = get_data(a)
  cell_f_phys = lazy_map(PushDofMap(),cell_f_ref,cell_m)
  CellDof(cell_f_phys,trian,DomainStyle(a))
end
