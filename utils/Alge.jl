"""
	eye(n, T)
Generate a n×n matrix in type `T` (Float64 by default)
"""
function eye(n, T = Float64)
	f(i,j) = T(i == j)
	return StaticArrays.sacollect(SMatrix{n,n}, f(i,j) for i in 1:n, j in 1:n)
end

function bitget(n,i)
	D = digits(n, base = 2)
	if i > length(D)
		return 0
	else
		return D[i]
	end
end

function NearZero(near)
	norm(near)<1e-15
end

function trace(T)
	n = size(T,2)
	sum([T[i,i] for i = 1:n])
end

"""
	Sdiag(v)
Generate a diagonal SMatrix with diagonal indices comes from v.
"""
function Sdiag(v::SVector)
	n = size(v,1);
	f(i,j) = i == j ? v[i] : 0.0;
	return StaticArrays.sacollect(SMatrix{n,n}, f(i,j) for i in 1:n, j in 1:n)
end

"""
	cvtSMatrix(SE3s)
Converting 4×4 Matrix into SMatrix.
"""
function cvtSMatrix(SE3s)
	n = size(SE3s,3);
	SSE3s = Vector{SMatrix{4,4,Float64,16}}(undef,n)
	for i = 1:n
		SSE3s[i] = SMatrix{4,4,Float64,16}(SE3s[:,:,i])
	end
	return SSE3s
end