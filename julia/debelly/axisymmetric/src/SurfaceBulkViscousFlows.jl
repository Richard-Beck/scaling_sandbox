module SurfaceBulkViscousFlows

using Gridap
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.FESpaces

using GridapEmbedded
using GridapEmbedded.AlgoimUtils
using GridapEmbedded.Interfaces: IN, OUT, CUT

using SparseMatricesCSR
using GridapPETSc
using GridapPETSc: PETSC

using LinearAlgebra: diag, cond
using Symbolics

using Plots
using DelimitedFiles

using Base: @kwdef # Para permitir la sintaxis Struct(campo=valor)

import Random

include("ActivityFunctions.jl")
include("AuxiliaryFunctions.jl")
include("InitialDensityFunctions.jl")
include("TangentOperators.jl")
include("SolverFunctions.jl")
include("WeakForms.jl")

include("MechanochemicalAxisymmetricVector.jl")

export unit_density
export verification
export mechanostability
  
export run_mechanochemical_axisymmetric_vector 

export MechanicalParams, KineticParams, CouplingParams, SimControl


end # module SurfaceBulkViscousFlows
