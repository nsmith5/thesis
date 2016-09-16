using Syntectic

s = State()

# Numerical Parameters
set!(s, :N, 256)
set!(s, :Δx, 0.125)
set!(s, :Δt, 0.00125)

# Correlation Function Parameters
set!(s, :k′, 2π)
set!(s, :α , 0.8)
set!(s, :β , 6.0)
set!(s, :ρ , √3/2)
set!(s, :αc, 0.4)
set!(s, :c₀, 0.5)

# Free Energy Parameters
set!(s, :χ, 1.0)
set!(s, :η, 2.0)
set!(s, :ϵ₀, 10.0)
set!(s, :σ₀, 0.350)
set!(s, :σ, 0.30)
set!(s, :ω, 0.30)
set!(s, :Wc, 1.0)
set!(s, :kbT, 0.0001)

# Dynamic Parameters
set!(s, :Mₙ, 1.0)
set!(s, :Mc, 1.0)
