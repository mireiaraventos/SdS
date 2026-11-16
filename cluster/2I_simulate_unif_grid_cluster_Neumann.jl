using LinearAlgebra
using FaADE
using JLD2
#using Revise


Δt = 1e-4
t_f = 1e-2


δ1 = parse(Float64, ARGS[1]) # 0.01, 0.05, 0.1, 0.5
δ2 = parse(Float64, ARGS[2]) # 0.01, 0.05, 0.1, 0.5
k_para = parse(Float64, ARGS[3]) # 1.0e10, 1.0e1
k_perp = 1.0
iota = parse(Float64, ARGS[5])

exp_k = log10(k_para/k_perp)

order = 2

# Spatial domain
nx = ny = parse(Int, ARGS[4]) # 201

Dx = [0.0,2.0]
Dy = [0.0,2π]
Dom = Grid2D(Dx,Dy,nx,ny)

# Initial condition
u₀(x, y) = x 

coord = :Cartesian

#=== BOUNDARY CONDITIONS ===# 
BoundaryLeft    = SAT_Neumann((y,t)->1.0, Dom.Δx, Left, order)
BoundaryRight   = SAT_Neumann((y,t)->1.0, Dom.Δx, Right, order)
BoundaryUp      = SAT_Periodic(Dom.Δy, Up, order)
BoundaryDown    = SAT_Periodic(Dom.Δy, Down, order)

BC = (BoundaryLeft, BoundaryRight, BoundaryUp, BoundaryDown)

println("Solving for δ1 = $δ1, δ2 = $δ2")
#=== PARALLEL MAP ===#
xₛ = 1.0  # Island position

function B(X, x, p, t) # x = [x, y]; t = z; # X = B(x)
    m1 = 1
    n1 = -1
    m2 = 1
    n2 = 1

    X[1] = δ1*x[1]*(2-x[1])*0.25*m1*sin(m1*x[2]-n1*t) + δ2*x[1]*(2-x[1])*0.25*m2*sin(m2*x[2]-n2*t)
    X[2] = iota*(x[1]-xₛ) + δ1*(2-2*x[1])*0.25*cos(m1*x[2]-n1*t) + δ2*(2-2*x[1])*0.25*cos(m2*x[2]-n2*t)
end

#parallel grid and grid data construction
gridoptions = Dict("interpmode"=>:chs,"xbound"=>[0.0,2.0],"ybound"=>[0.0,2π])

PGrid = construct_grid(B, Dom, [-2π,2π], ymode=:period, interpmode=:chs, gridoptions=gridoptions)
PData = ParallelData(PGrid, Dom, order, κ=k_para, interpolant=:chs)

#construct problem
P = Problem2D(order, u₀, k_perp, k_perp, Dom,BC, parallel=PData)

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
    "num_results/2I_k$(exp_k)_tf$(t_f)_d$(δ1)_$(δ2)_res$(nx)_m1_dT_iota$(iota)_Neumann.jld2",
    "soln", soln,
    "x_mid", x_mid,
    "u_mid", u_mid,
    "δ1", δ1,
    "δ2", δ2,
    "t_f", t_f,
    "iota", iota
)