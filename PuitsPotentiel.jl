include("lobpcg.jl")
using FFTW
using Plots
using LinearAlgebra

mutable struct Params
	dim
	n
	k_axis
	k2_grid
	k2lin
	solver
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
	Params(dim,n,k_axis,k2_grid,k2lin,"lobpcg")
end

linearize(X) = vcat(X...)
Reshape(Xlin,p) = reshape(Xlin,Tuple(fill(p.n,p.dim)))

convolve(X,V) = fft(ifft(X).*ifft(V))# X et V en fourier

#V doit être un tenseur d'ordre dim, de taille n
#k : mode recherché, doit être plus petit que 5 et plus grand que 1
function solve(V,p,k)
	N = p.n^p.dim
	V_four = fft(V)
	H = X->4*pi^2 .*p.k2lin.*X .+ linearize(convolve(V_four,Reshape(X,p))) # X est en fourier et linéraire
	l = k+1 # nb modes propres
	#if p.solver == "lobpcg"
	(λs,ϕs,cv) = solve_lobpcg(H,N,l,p.k2lin;tol=1e-7)
	#else
		#eis = []
		#[H(ei) for ei in ]
		#(λs,ϕs) = eigen(Array(H))
	#end
	# k ieme etat excite
	ψ = ϕs[k]; Eψ = λs[k]
	ψ = Reshape(ψ,p)
	Mψ = ifft(ψ)
	Mψ ./= norm(Mψ)
	(Eψ,Mψ)
end

#Une dimension. x pour l'abscisse, V pour le potentiel, k pour le mode
function D1(x,V,k)
	p = init_struct(100,1)
	(E,ϕ) = solve(V,p,k)
	pl = plot(x,abs2.(ϕ))
	savefig(pl,"pot1d.pdf")
end

x = [(i/p.n) for i=1:p.n]
D1(x,-200*exp.((-(x.-0.5).^2)./ 0.1),4)

#2D, y pour l'ordonnée, V pour le potentiel (tenseur d'ordre 2)
function D2(x,y,V,k)
	p = init_struct(100,2)
	(E,ϕ) = solve(V,p,k)
	pl = Plots.heatmap(x,y,abs2.(ϕ))
	surf = Plots.surface(x,y,abs2.(ϕ))
	savefig(pl,"pot2dheat.pdf")
	savefig(surf,"pot2surf.pdf")
end

function f(x,y)
	-2000*exp(-((x-0.5)^2+(y-0.5)^2)/0.1)
end

D2(x,x,[ f(x[i],y[j]) for i=1:p.n, j=1:p.n],3)
