using LinearAlgebra
using FaADE
using JLD2
#using Revise


Δt = 1e-4
t_f = 1e-2


δ = parse(Float64, ARGS[1]) # 0.01, 0.05, 0.1, 0.5
k_para = parse(Float64, ARGS[2]) # 1.0e10, 1.0e1
k_perp = 1.0

exp_k = log10(k_para/k_perp)

order = 2

# Spatial domain
nx = ny = parse(Int, ARGS[3]) # 201

Dx = [0.0,1.0]
Dy = [0.0,2π]
Dom = Grid2D(Dx,Dy,nx,ny)

# Initial condition
u₀(x, y) = x 

coord = :Cartesian

#=== BOUNDARY CONDITIONS ===# 
BoundaryLeft    = SAT_Dirichlet((y,t)->0.0, Dom.Δx, Left, order)
BoundaryRight   = SAT_Dirichlet((y,t)->1.0, Dom.Δx, Right, order)
BoundaryUp      = SAT_Periodic(Dom.Δy, Up, order)
BoundaryDown    = SAT_Periodic(Dom.Δy, Down, order)

BC = (BoundaryLeft, BoundaryRight, BoundaryUp, BoundaryDown)

println("Solving for δ = $δ")
#=== PARALLEL MAP ===#
xₛ = 0.5  # Island position

function B(X, x, p, t) # x = [x, y]; t = z; # X = B(x)
    X[1] = δ*x[1]*(1 - x[1])*sin(x[2])
    X[2] = 2x[1] - 2*xₛ + δ*(1 - x[1])*cos(x[2]) - δ*x[1]*cos(x[2])
end

#parallel grid and grid data construction
gridoptions = Dict("interpmode"=>:chs,"xbound"=>[0.0,1.0],"ybound"=>[0.0,2π])

PGrid = construct_grid(B, Dom, [-2π,2π], ymode=:period, interpmode=:chs, gridoptions=gridoptions)
PData = ParallelData(PGrid, Dom, order, κ=k_para, interpolant=:chs)

#construct problem
P = Problem2D(order, u₀, k_perp, k_perp,Dom,BC, parallel=PData)

println("Done")

println("Solving")

soln = solve(P, Dom, Δt, t_f)


# Mid-plane index
midy = floor(Int, ny / 2) + 1
midx = floor(Int, nx / 2) + 1

# Concatenate u(x) across the three blocks
x_mid = Dom.gridx[:, midy]

u_mid = vcat(soln.u[2][:, midy])


save(
    "num_results/SI_k$(exp_k)_tf$(t_f)_d$(δ)_res$(nx)_dT.jld2",
    "soln", soln,
    "x_mid", x_mid,
    "u_mid", u_mid,
    "δ", δ,
    "t_f", t_f
)