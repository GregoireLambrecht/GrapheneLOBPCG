include("lobpcg.jl ")
include("DiagBand.jl")
using FFTW
using Plots
using KrylovKit
using LinearAlgebra


O = [0,0]
K = b1
M = b1 .+ b2 ./2
nc = 10
T = [i/nc for i = 0:nc]
ϕ1 = [(1 - T[i]) .* O .+ T[i] .* K for i = 1:(nc+1)]
ϕ2 = [(1 - T[i]) .* K .+ T[i] .* M for i = 1:(nc+1)]
ϕ3 = [(1 - T[i]) .* M .+ T[i] .* O for i = 1:(nc+1)]

λ1 = solveVect.(ϕ1)
λ2 = solveVect.(ϕ2)
λ3 = solveVect.(ϕ3)

plot(norm.(ϕ1),linearize(λ1))
plot!(norm.(ϕ2),linearize(λ2))
plot!(norm.(ϕ3),linearize(λ3))
