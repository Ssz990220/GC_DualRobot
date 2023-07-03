# using Grassmann
# using StaticArrays
# @basis "1,1,1,0"

function katyushaX(V,P,get_i_gard_func,get_full_grad_func, err_func,normfunc, conf)
	max_iter = conf.max_iter; m = conf.m; η₀ = conf.η; n = size(P[1],1);
	τ = conf.τ; 
	∇ᵥg(V) = get_full_grad_func(V,P,conf)
	∇ᵥgᵢ(V,i) = get_i_gard_func(V,P,i,conf);
	local yₜ = V; local yₜ₁ = V; local xₜ = V;
	local ∇ᵥgṼₜ = ∇ᵥg(V); local xₜ₁ = V;
	conf.err=="Ext" ? err = zeros(max_iter,conf.errdim) : err = zeros(max_iter,1);
	@inbounds for i = 1:max_iter
		xₜ₁ = xₜ
		xₜ = update(yₜ,yₜ₁,xₜ₁,τ)
		∇ᵥgṼₜ = ∇ᵥg(xₜ)
		yₖ = xₜ
		@inbounds for k = 1:m
			iₖ = rand(1:n)
			yₖ = yₖ .- η₀.*(∇ᵥgᵢ(yₖ,iₖ).-∇ᵥgᵢ(xₜ,iₖ).+∇ᵥgṼₜ)
		end
		yₜ₁ = yₜ
		if conf.reg
			yₜ = yₖ
		else
			yₜ = normfunc(yₖ)
		end
		if conf.err=="Ext"
			err[i,:] = err_func(yₜ,P,∇ᵥgṼₜ,0.0,conf)
		else
			err[i] = err_func(yₜ,P,∇ᵥgṼₜ,0.0,conf)
		end
	end
	return yₜ,err
end

@inline function update(yₜ,yₜ₁,xₜ₁,τ)
	(3.0/2*yₜ+1.0/2*xₜ₁-(1-τ)*yₜ₁)/(1.0+τ);
end

function GradientDescent(V,P,get_i_gard_func,get_full_grad_func, err_func,normfunc, conf)
	max_iter = conf.max_iter; m = conf.m; η₀ = conf.η; n = size(P[1],1);
	∇ᵥg(V) = get_full_grad_func(V,P,conf)
	local Ṽₜ = V; local Ṽₜ₁ = Ṽₜ; local ∇ᵥgṼₜ = ∇ᵥg(Ṽₜ); local ∇ᵥgṼₜ₁ = ∇ᵥgṼₜ;
	local ηₜ = 0;
	conf.err=="Ext" ? err = zeros(max_iter,conf.errdim) : err = zeros(max_iter,1);
	@inbounds for i = 1:max_iter
		∇ᵥgṼₜ = ∇ᵥg(Ṽₜ)
		Ṽₜ = Ṽₜ .- η₀.* ∇ᵥgṼₜ
		if conf.reg
		else
			Ṽₜ = normfunc(Ṽₜ)
		end
		if conf.err=="Ext"
			err[i,:] = err_func(Ṽₜ,P,∇ᵥgṼₜ,ηₜ,conf)
		else
			err[i] = err_func(Ṽₜ,P,∇ᵥgṼₜ,ηₜ,conf)
		end
	end
	return Ṽₜ,err
end

function SVRG(V,P,get_i_gard_func,get_full_grad_func, err_func,normfunc, conf)
	max_iter = conf.max_iter; m = conf.m; η₀ = conf.η; n = size(P[1],1);
	∇ᵥg(V) = get_full_grad_func(V,P,conf)
	∇ᵥgᵢ(V,i) = get_i_gard_func(V,P,i,conf);
	local Ṽₜ = V; local Ṽₜ₁ = Ṽₜ; local ∇ᵥgṼₜ = ∇ᵥg(Ṽₜ); local ∇ᵥgṼₜ₁ = ∇ᵥgṼₜ;
	local ηₜ = 0;
	conf.err=="Ext" ? err = zeros(max_iter,conf.errdim) : err = zeros(max_iter,1);
	@inbounds for i = 1:max_iter
		∇ᵥgṼₜ = ∇ᵥg(Ṽₜ)
		Vₖ = Ṽₜ
		@inbounds for k = 1:m
			iₖ = rand(1:n)
			Vₖ = Vₖ .- η₀.*(∇ᵥgᵢ(Vₖ,iₖ).-∇ᵥgᵢ(Ṽₜ,iₖ).+∇ᵥgṼₜ)
		end
		Ṽₜ₁ = Ṽₜ
		if conf.reg
			Ṽₜ = Vₖ
		else
			Ṽₜ = normfunc(Vₖ)
		end
		# Ṽₜ = Vₖ
		∇ᵥgṼₜ₁ = ∇ᵥgṼₜ;
		if conf.err=="Ext"
			err[i,:] = err_func(Ṽₜ,P,∇ᵥgṼₜ,ηₜ,conf)
		else
			err[i] = err_func(Ṽₜ,P,∇ᵥgṼₜ,ηₜ,conf)
		end
	end
	return Ṽₜ,err
end

function SVRG_BB(V,P,get_i_gard_func, get_full_grad_func, get_ηₜ, err_func, conf)
	max_iter = conf.max_iter; m = conf.m; η₀ = conf.η; n = size(P[1],1);
	∇ᵥg(V) = get_full_grad_func(V,P,conf)
	∇ᵥgᵢ(V,i) = get_i_gard_func(V,P,i,conf);
	local Ṽₜ = V; local Ṽₜ₁ = Ṽₜ; local ∇ᵥgṼₜ = ∇ᵥg(Ṽₜ); local ∇ᵥgṼₜ₁ = ∇ᵥgṼₜ;
	local ηₜ = 0;
	conf.err=="Ext" ? err = zeros(max_iter,conf.errdim) : err = zeros(max_iter,1);
	@inbounds for i = 1:max_iter
		∇ᵥgṼₜ = ∇ᵥg(Ṽₜ)
		i > 1 ? ηₜ =  get_ηₜ(Ṽₜ,Ṽₜ₁,∇ᵥgṼₜ,∇ᵥgṼₜ₁,conf) : ηₜ = η₀;
		if conf.stop
			if grad_norm(∇ᵥgṼₜ) < 1e-14
				return Ṽₜ,err[1:i-1]
			end
		end
		Vₖ = Ṽₜ
		@inbounds for k = 1:m
			iₖ = rand(1:n)
			Vₖ = Vₖ .- ηₜ.*(∇ᵥgᵢ(Vₖ,iₖ).-∇ᵥgᵢ(Ṽₜ,iₖ).+∇ᵥgṼₜ)
		end
		Ṽₜ₁ = Ṽₜ
		Ṽₜ = Vₖ
		∇ᵥgṼₜ₁ = ∇ᵥgṼₜ;
		if conf.err=="Ext"
			err[i,:] = err_func(Ṽₜ,P,∇ᵥgṼₜ,ηₜ,conf)
		else
			err[i] = err_func(Ṽₜ,P,∇ᵥgṼₜ,ηₜ,conf)
		end
	end
	return Ṽₜ,err
end

function grad_norm(∇ᵥgṼₜ)
	return norm(∇ᵥgṼₜ,2)
end

function grad_norm(∇ᵥgṼₜ::Union{SVector,Vector})
	return norm(∇ᵥgṼₜ)
end