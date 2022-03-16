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
	if p.solver == "lobpcg"
		(λs,ϕs,cv) = solve_lobpcg(H,N,l,p.k2lin;tol=1e-7)
		ψ = ϕs[k]; Eψ = λs[k]
	else
		eis = zeros(p.n,p.n,p.n,p.n)
		for i =1:p.n
			for j=1:p.n
				eis[i,j,i,j] = 1
			end
		end
		MH = [H(linearize([eis[i,j,k,l] for k=1:p.n, l=1:p.n])) for i=1:p.n, j=1:p.n]
		MH = [MH[i][j] for i=1:N,j=1:N]
		(λs,ϕs) = eigen(MH)
		ψ = [ϕs[i,k] for i=1:N]; Eψ = λs[k]
	end
	ψ = Reshape(ψ,p)
	Mψ = ifft(ψ)
	Mψ ./= norm(Mψ)
	(Eψ,Mψ)
end

#Une dimension. x pour l'abscisse, V pour le potentiel, k pour le mode, n pour la précision
function D1(x,V,n,k)
	p = init_struct(n,1)
	(E,ϕ) = solve(V,p,k)
	pl = plot(x,abs2.(ϕ),xlabel = "x", ylabel = "Densité", title = "Densité de Probabilité pour le mode $k")
	savefig(pl,"pot1d.pdf")
end


#2D, y pour l'ordonnée, V pour le potentiel (tenseur d'ordre 2)
function D2(x,y,V,n,k)
	p = init_struct(n,2)
	(E,ϕ) = solve(V,p,k)
	pl = Plots.heatmap(x,y,abs2.(ϕ), xlabel = "x",ylabel = "y", title = "Densité de Probabilité pour le mode $k")
	surf = Plots.surface(x,y,abs2.(ϕ), xlabel = "x",ylabel = "y",title = "Densité de Probabilité pour le mode $k")
	savefig(pl,"pot2dheat.pdf")
	savefig(surf,"pot2surf.pdf")
end

function f(x,y)
	-2000*exp(-((x-0.5)^2+(y-0.5)^2)/0.1)
end

function timeCompare(enu)
	N = [i for i=10:enu]
	X = [[i/N[j] for i=1:N[j]] for j=1:(enu-9)]
	V = [[ f(X[l][i],X[l][j]) for i=1:N[l], j=1:N[l]] for l=1:(enu-9)]
	P = [init_struct(N[j],2) for j=1:(enu-9)]
	TLOBPCG = [(@elapsed solve(V[i],P[i],1)) for i=1:(enu-9)]
	for i=1:(enu-9)
		P[i].solver = "eig"
	end
	TLIN = [(@elapsed solve(V[i],P[i],1)) for i=1:(enu-9)]
	pl = plot(N,TLOBPCG,xlabel = "N : taille de la discrétisation", ylabel = "Temps en seconde",title = "Comparaison des temps d'éxecution",label = "LOBPCG")
	pl = plot!(N,TLIN, label = "LinearAlgebra")
	savefig(pl,"Times.pdf")
end

n=200
x = [(i/n) for i=1:n]
D1(x,-200*exp.((-(x.-0.5).^2)./ 0.1),n,4)
D2(x,x,[ f(x[i],x[j]) for i=1:n, j=1:n],n,2)
timeCompare(25)
