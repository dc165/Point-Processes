using KDTrees
struct MPoxModel{T <: AbstractFloat}
    S           :: Matrix{T}
    B           :: Vector{Vector}
    KDT         :: KDTree
    mu_xy       :: Matrix{T} 
    C_xy        :: Array{T, 3}
    theta_B     :: Vector{Vector{T}}
    theta_kappa :: Matrix{T}
    d_theta     :: Vector{T}
    a_mu        :: Vector{T}
    lambda_mu   :: Vector{T}
    
end