include("lobpcg.jl ")
using FFTW
using Plots
using KrylovKit
using LinearAlgebra

mutable struct Params
	dim #La dimension : 1, 2 ou 3
	n# La précision
	k_grid
	k2_grid
	k2lin #Le laplacien linéarisé (tableau de taille n^dim)
end

function init_struct(n,dim)
	k_axis = fftfreq(n)*n
	k_grid = [[k_axis[i],k_axis[j]] for i=1:n, j=1:n]
	k2_grid = [norm(k_grid[i,j])^2 for i=1:n, j=1:n]
	k2_grid = norm.(k2_grid).^2
	k2lin = linearize(k2_grid)
	Params(dim,n,k_grid,k2_grid,k2lin)
end

linearize(X) = vcat(X...)
Reshape(Xlin,p) = reshape(Xlin,Tuple(fill(p.n,p.dim)))

convolve(X,V) = fft(ifft(X).*ifft(V))# X et V en fourier

#V doit être un tenseur d'ordre 2, de taille n
#le k du dual pour Hk
#l : nb de valeurs propres renvoyées
function solve(V,p,k,l)
	N = p.n^p.dim
	δ = [1im .* k .+ p.k_grid[i,j] for i=1:p.n,j=1:p.n]
	δ = δ'δ
	Δklin = linearize(δ)
	VLinFour = linearize(fft(V))
	H = X->-Δklin.*X + convolve(VLinFour,X) # X est en fourier et linéraire
	l = 100 # nb modes propres
	(λs,ϕs,cv) = solve_lobpcg(H,N,l,p.k2lin;tol=1e-7)
	λs
end


p = init_struct(100,2)
#le k du dual pour Hk
k = [1,1]
V = [ rand() for i=1:p.n,j = 1:p.n]
solve(V,p,k,5)
