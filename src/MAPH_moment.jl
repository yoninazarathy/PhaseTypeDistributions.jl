include("MAPH.jl")

absorption_probs(maph::MAPHDist) = maph.α*maph.R
kth_moment(maph::MAPHDist, k::Real) = (-1.0)^(k+1) * factorial(k) .* maph.α * inv(maph.T)^(k+1) * maph.D ./ absorption_probs(maph)
mean(maph::MAPHDist) = kth_moment(maph, 1)
var(maph::MAPHDist) = kth_moment(maph, 2) - mean(maph).^2
scv(maph::MAPHDist) = var(maph) ./ (mean(maph).^2) 

sub_distribution(maph::MAPHDist, k::Real) = x -> maph.α * (exp(x * maph.T) - I) * inv(maph.T) * maph.D[:,k] 
sub_pdf(maph::MAPHDist, k::Real) = x-> maph.α * exp(maph.T*x) * maph.D[:, k]
mgf(maph::MAPHDist, k::Real) = z-> -maph.α * inv(z*I + maph.T) * maph.D[:, k]

