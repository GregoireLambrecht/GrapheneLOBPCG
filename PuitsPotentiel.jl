include("lobpcg.jl ")
using FFTW
using Plots
using KrylovKit
using LinearAlgebra

mutable struct Params
	dim
	n
	k_axis
	k2_grid
	k2lin
end

function init_struct(n,dim)
	k_axis = fftfreq(n)*n
	if (dim == 1)
		k_grid = k_axis
		k2_grid = norm.(k_grid)
		k2_grid = norm.(k2_grid).^2
		k2lin = k2_grid
	elseif (dim == 2)
		k_grid = [[k_axis[i],k_axis[j]] for i=1:n, j=1:n]
		k2_grid = [norm(k_grid[i,j])^2 for i=1:n, j=1:n]
		k2_grid = norm.(k2_grid).^2
		k2lin = linearize(k2_grid)
	else
		k2_grid = [[k_axis[i],k_axis[j], k_axis[k]] for i=1:n, j=1:n, k = 1:n]
		k2_grid = norm.(k2_grid).^2
		k2lin = linearize(k2_grid)
	end
	Params(dim,n,k_axis,k2_grid,k2lin)
end

linearize(X) = vcat(X...)
Reshape(Xlin,p) = reshape(Xlin,Tuple(fill(p.n,p.dim)))

convolve(X,V) = fft(ifft(X).*ifft(V))# X et V en fourier

#V doit être un tenseur d'ordre dim, de taille n
#k : mode recherché, doit être plus petit que 5 et plus grand que 1
function solve(V,p,k)
	N = p.n^p.dim
	VLinFour = linearize(fft(V))
	H = X->-4*pi^2 .*p.k2lin.*X + convolve(VLinFour,X) # X est en fourier et linéraire
	l = 5 # nb modes propres
	(λs,ϕs,cv) = solve_lobpcg(H,N,l,p.k2lin;tol=1e-7)
	#(λs,ϕs) = eigsolve(H,p.k2lin,l)
	# k ieme etat excite
	ψ = ϕs[k]; Eψ = λs[k]
	ψ = Reshape(ψ,p)
	Mψ = ifft(ψ)
	Mψ = real.(conj.(Mψ) .* Mψ)
	Mψ ./= norm(Mψ)
	(Eψ,Mψ)
end

#Une dimension
p = init_struct(200,1)
x = [(i*1/p.n -0.5) for i=1:p.n]
V = -exp.(-(x.^2)./ 0.01)
(E,ϕ) = solve(V,p,1)
plot(x,ϕ)
#plot!(x,V)
#Deux dimensions
par = init_struct(100,2)
y = [(i*1/p.n -0.5) for i=1:p.n]
VS = [ exp(-(x[i]^2 +x[j]^2)/0.01) for i=1:par.n,j = 1:par.n]
(ES,ϕS) = solve(VS,par,1)
Plots.surface(x,y,ϕS)