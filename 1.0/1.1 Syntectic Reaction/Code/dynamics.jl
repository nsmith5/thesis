export step!

function ΔFmix(c::Float64, s::State)
	@unpack ω, σ, σ₀, ϵ₀ = s
	ϵ = -4.0 + ϵ₀*(σ - σ₀)
	return ω * (c * log(2.0 * c) + (1 - c)*log(2.0 * (1. - c)) + 0.5 * ϵ * (c - 0.5) ^ 2)
end

δΔFmixδc(c::Float64, s::State) = s.ω*(log(2.0 * c) - log(2.0 * (1 - c)))

function noise!(s::State)
	@unpack N, kbT, Mc, Mₙ, Δx, Δt = s
	c_scale = im*√(kbT*Mc/(Δx^2*Δt))
	n_scale = im*√(kbT*Mₙ/(Δx^2*Δt))
	for j in 1:N
		for i in 1:N>>1 + 1
			k = index_to_k(s, i, j)
			s.ξc[i, j] = k*c_scale*complex(randn(), randn())
			s.ξₙ[i, j] = k*n_scale*complex(randn(), randn()) 
		end
	end
	return nothing
end

λ(c::Float64, s::State) = exp(-(c-0.5)^2/(2.0*s.αc^2))

function set_nonlinear!(s::State)
	@unpack N, n, c, η, χ, αc, C₂n = s
	third = 1./3.
	for i in 1:N*N
		s.nₙₗ[i] = n[i]^2*(-0.5*η + third*χ*n[i]) + ΔFmix(c[i], s) - λ(c[i], s)*C₂n[i]
		s.cₙₗ[i] = (1+n[i])*δΔFmixδc(c[i], s) + 0.5*n[i]*((c[i]-0.5)/αc^2*λ(c[i], s))*C₂n[i]
	end
end

function calccorr!(s::State)
	@unpack N, C₂n, kC₂n, C₂, ñ, fftplan = s
	for i in 1:(s.N>>1+1)*s.N
        kC₂n[i] = C₂[i]*ñ[i]
    end
    A_ldiv_B!(C₂n, fftplan, kC₂n)
    return nothing
end

function step!(s::State)
	@unpack ζ, Δt, ∇², Mₙ, Mc, ω, c̃, ñ, n, c, Wc, ξc, ξₙ = s
	A_mul_B!(ñ, s.fftplan, n)
	A_mul_B!(c̃, s.fftplan, c)
	calccorr!(s)
	set_nonlinear!(s)
	A_mul_B!(s.ñₙₗ, s.fftplan, s.nₙₗ)
	A_mul_B!(s.c̃ₙₗ, s.fftplan, s.cₙₗ)
	noise!(s)
	ϵ = -4.0 + s.ϵ₀*(s.σ - s.σ₀)
	for i in 1:(s.N>>1 + 1)*s.N
		Λ = Mc*∇²[i]*(ω*ϵ - Wc*∇²[i])
		prefn = 1.0/(1.0 - ζ*Δt*Mₙ*∇²[i])
		prefc = 1.0/(1.0 - ζ*Δt*Λ)
		s.ñ[i] = prefn*((1.0 + (1 - ζ)*Δt*Mₙ*∇²[i])*ñ[i] + Mₙ*Δt*∇²[i]*s.ñₙₗ[i]+Δt*ξₙ[i])
		s.c̃[i] = prefc*((1.0 + (1 - ζ)*Δt*Λ)*c̃[i] + Mc*Δt*∇²[i]*s.c̃ₙₗ[i] + Δt*ξc[i])	
	end
	A_ldiv_B!(s.n, s.fftplan, s.ñ)
	A_ldiv_B!(s.c, s.fftplan, s.c̃)
	return nothing
end
