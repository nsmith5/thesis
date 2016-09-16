export step!

function ΔFmix(c::Float64, s::State)
	ϵ = -4.0 + s.ϵ₀*(s.σ - s.σ₀)
	return s.ω*(c*log(2.0c) + (1-c)*log(2.0*(1-c)) + 0.5*ϵ*(c-0.5)^2)
end

δΔFmixδc(c::Float64, s::State) = s.ω*(log(2.0c) - log(2.0*(1-c)))

function noise!(s::State)
	c_scale = sqrt(2)*s.kbT*s.Mc/(s.Δx^2*s.Δt)
	n_scale = sqrt(2)*s.kbT*s.Mₙ/(s.Δx^2*s.Δt)
	for i in 1:length(s.ξc)
		s.ξc[i] = c_scale*complex(randn(), randn())
		s.ξₙ[i] = n_scale*complex(randn(), randn()) 
	end
	return nothing
end

function set_nonlinear!(s::State)
	third = 1/3
	for i in 1:s.N*s.N
		s.nₙₗ[i] = s.n[i]^2*(-0.5*s.η + third*s.χ*s.n[i]) + ΔFmix(s.c[i], s) - exp(-(s.c[i]-0.5)^2/(2s.αc^2))*s.C₂n[i]
		s.cₙₗ[i] = (s.n[i] + 1.)*δΔFmixδc(s.c[i], s) - 0.5*s.n[i]*(-(s.c[i]-0.5)/s.αc^2*exp(-(s.c[i]-0.5)^2/(2s.αc^2))*s.C₂n[i])
	end
end

function calccorr!(s::State)
	for i in 1:(s.N>>1+1)*s.N
        s.kC₂n[i] = s.C₂[i]*s.ñ[i]
    end
    A_ldiv_B!(s.C₂n, s.fftplan, s.kC₂n)
    return nothing
end

function step!(s::State)
	A_mul_B!(s.ñ, s.fftplan, s.n)
	A_mul_B!(s.c̃, s.fftplan, s.c)
	calccorr!(s)
	set_nonlinear!(s)
	A_mul_B!(s.ñₙₗ, s.fftplan, s.nₙₗ)
	A_mul_B!(s.c̃ₙₗ, s.fftplan, s.cₙₗ)
	noise!(s)
	ϵ = -4.0 + s.ϵ₀*(s.σ - s.σ₀)
	for i in 1:(s.N>>1 + 1)*s.N
		s.ñ[i] = 1.0/(1.0-s.Δt*s.∇²[i])*(s.ñ[i] + s.Δt*s.∇²[i]*(s.ñₙₗ[i]+s.ξₙ[i]))
		s.c̃[i] = 1.0/(1.0-s.Δt*s.∇²[i]*(s.ω*ϵ-s.Wc*s.∇²[i]))*(s.c̃[i] + s.Δt*s.∇²[i]*(s.c̃ₙₗ[i] + s.ξc[i]))	
	end
	A_ldiv_B!(s.n, s.fftplan, s.ñ)
	A_ldiv_B!(s.c, s.fftplan, s.c̃)
	return nothing
end


