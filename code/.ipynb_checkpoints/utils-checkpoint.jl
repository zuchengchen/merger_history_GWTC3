using Cuba
include("events.jl")
using Revise
using Cuba
using HypergeometricFunctions
using SpecialFunctions
# using GSL
using BenchmarkTools
using Interpolations
using Distributions
using KernelDensity
using DelimitedFiles
using Memoize
using StatsPlots
using NestedSamplers
using StatsBase: sample
using Random
# using MCMCChains
using FileIO
using BAT
using IntervalSets
using FileIO
using LaTeXStrings
using Roots
# using GWSC
using IntervalSets
using AbstractMCMC
# AbstractMCMC.setprogress!(false)
Random.seed!(8452);

const m_min = 3
const m_max = 130
const δm = m_max - m_min

# function beta_func(x, y, a, b)
#     result = sf_beta_inc(a, b, y) - sf_beta_inc(a, b, x)
#     result * sf_beta(a, b)
# end

@memoize function beta_func2(x, y, a, b)
    δ = y - x
    int0(t) = t^(a-1) * (1-t)^(b-1)
    int(t) = δ * int0(x + δ*t)
    function integrand(x, f)
        f[1] = int(x[1])
    end
    result, err = cuhre(integrand, rtol=1e-1)
    result[1]
end

function get_post_func(event)
    post_file = "LIGO-data/masses_source_frame_posterior/GW$(event)_posterior.txt"
    post = readdlm(post_file);
    d1 = fit(MvNormal, transpose(post))
    
    μ1, μ2 = d1.μ
    σ1 = sqrt(d1.Σ[1])
    σ2 = sqrt(d1.Σ[4])
    ρ = d1.Σ[2]/σ1/σ2
    norm = 2pi * sqrt(1 - ρ^2) * σ1 * σ2
    
    post(x1, x2) = 1/norm * exp(- 1/2/(1-ρ^2) * ((x1-μ1)^2/σ1^2 - 2ρ*(x1-μ1)*(x2-μ2)/σ1/σ2 + (x2-μ2)^2/σ2^2))
end

#### 
# vt
const time_O1 = 46.1; # (* days or 48.6 days *)
const time_O2 = 117; #(* days *)
const time_O3a = 149.9 #;(* days *)

function get_VT_GWTC2()
    m1s = m2s = collect(1.0:1.0:136.0)
    vts_O1O2 = (time_O1 + time_O2)/365 * reshape(readdlm("backup/VT_1yr_m1m2_LIGO_O1.txt")[:, 3], length(m1s), length(m1s))
    vts_O3a = (time_O3a)/365 * reshape(readdlm("backup/VT_1yr_m1m2_LIGO_O3.txt")[:, 3], length(m1s), length(m1s))
    vts = vts_O1O2 + vts_O3a
    
    vt_int0 = LinearInterpolation((m1s, m2s), vts)
    
    VT(m1, m2) = vt_int0(m1, m2)
end

VT = get_VT_GWTC2()