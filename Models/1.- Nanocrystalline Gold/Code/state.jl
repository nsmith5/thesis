import Base: show,
			 get

export State,
       set!,
       get

type State
	N_init::Bool

	# Free Energy Parameters
	η::Float64
	χ::Float64
	Wc::Float64
	ϵ₀::Float64
	σ₀::Float64
	σ::Float64
	ω::Float64

	# Dynamic Parameters
	Mₙ::Float64
	Mc::Float64

	# Correlation Function Parameters
	k′::Float64
	α::Float64
	β::Float64
	ρ::Float64
	αc::Float64

	# Numerical Parameters
	Δx::Float64
	Δt::Float64
	N::Int

	# Temporary Fields
	ñₙₗ::Array{Complex128, 2}
	c̃ₙₗ::Array{Complex128, 2}
	c̃::Array{Complex128, 2}
	ñ::Array{Complex128, 2}
	kC₂n::Array{Complex128, 2}
	ξₙ::Array{Complex128, 2}
	ξc::Array{Complex128, 2}

	C₂n::Array{Float64, 2}
	nₙₗ::Array{Float64, 2}
	cₙₗ::Array{Float64, 2}
	fftplan::Base.DFT.FFTW.FFTWPlan

	# Operators
	∇²::Array{Float64, 2}
	C₂::Array{Float64, 2}

	# Actual State
	n::Array{Float64, 2}
	c::Array{Float64, 2}

	State() = new(false)
end

function show(io::IO, s::State)
	# Pretty printing of a state
	println(io, "Simulation State")
	symbols = [:η, :χ, :Wc, :ϵ₀, :σ₀, :σ, :ω]
	println(io, "\nThermodynamic Parameters:")
	for sym in symbols
		println(io, "    ", sym, " = ", eval(:($s.$sym)))
	end
	println(io, "\nDynamic Parameters")
	println(io, "    Mₙ = ", s.Mₙ)
	println(io, "    Mc = ", s.Mc)
	println(io, "\nNumerical Parameters")
	println(io, "    Δx = ", s.Δx)
	println(io, "    Δt = ", s.Δt)
	print(io, "    N  = ", s.N)
end

function set!(s::State, sym::Symbol, val)
	if sym ∉ fieldnames(State)
		error("$sym not a field of type State")
	elseif sym ≠ :N && !s.N_init
		error("N must be initialized before all other state parameters")
	elseif sym ∈ [:σ, :k′, :α, :β, :ρ]
		@assert typeof(val) == Float64
		@eval $s.$sym = $val
		setC₂!(s)
		return nothing
	elseif sym == :Δx
		@assert typeof(val) == Float64
		s.Δx = val
		set∇²!(s)
		return nothing
	elseif sym == :N
		s.N_init = true
		s.N = val
		for arr_sym in [:ñₙₗ, :c̃ₙₗ, :c̃, :ñ, :kC₂n, :ξₙ, :ξc]
			@eval $s.$arr_sym = Array(Complex128, $val>>1 + 1, $val)
		end
		for arr_sym in [:∇², :C₂]
			@eval $s.$arr_sym = Array(Float64, $val>>1 + 1, $val)
		end
		for arr_sym in [:C₂n, :nₙₗ, :cₙₗ, :n, :c]
			@eval $s.$arr_sym = Array(Float64, $val, $val)
		end
		FFTW.set_num_threads(4)
		s.fftplan = plan_rfft(s.nₙₗ, (1, 2), flags=FFTW.MEASURE)
		return nothing
	else
		@assert @eval typeof($s.$sym) == typeof($val)
		@eval $s.$sym = $val
		return nothing
	end
end

function get(sym::Symbol, s::State)
	if sym ∈ fieldnames(State	)
		return @eval $s.$sym
	else
		error("$sym not a field of type State")
	end
end

function index_to_k(s::State, i::Int64, j::Int64)
	N = s.N
	L = s.N*s.Δx
	kx² = i < N>>1 ? (2π*(i-1)/L)^2 : (2π*(i-N-1)/L)^2
	ky² = j < N>>1 ? (2π*(j-1)/L)^2 : (2π*(j-N-1)/L)^2
	sqrt(kx² + ky²)
end

function set∇²!(s::State)
	N = s.N
	out = Array(Float64, N>>1 + 1, N)
	for j in 1:N
		for i in 1:N>>1 + 1
			out[i,j] = -index_to_k(s, i, j)^2
		end
	end
	s.∇² = out
	return
end

function setC₂!(s::State)
	symlist = [:N, :σ, :k′, :α, :β, :ρ]
	for sym in symlist
		@eval $sym = $s.$sym
	end
	C = Array(Float64, N>>1 + 1, N)
	for j in 1:N
		for i in 1:N>>1 + 1
			k = index_to_k(s, i, j)
			C[i,j] = exp(-σ^2*k′^2/(2*ρ*β))*exp(-(k-k′)^2/(2α^2))
		end
	end
	s.C₂ = C
	return
end

#function save(filename, s::State)
#	 TODO Make a function to save a state
#end

#function load!(filename)
#	# TODO make a function to load a state from file
#end
