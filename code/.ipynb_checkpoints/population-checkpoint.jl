# add Cuba HypergeometricFunctions SpecialFunctions BenchmarkTools Interpolations Distributions KernelDensity DelimitedFiles StatsPlots NestedSamplers StatsBase Random MCMCChains AbstractMCMC BAT IntervalSets FileIO JLD2

include("events.jl")
using Revise
using Cuba
using HypergeometricFunctions
using SpecialFunctions
using BenchmarkTools
using Interpolations
using Distributions
using KernelDensity
using DelimitedFiles
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

const m_min = 1.0
const m_max = 130.0
const δm = m_max - m_min
const σM = 6e-3
# const base_dir = homedir() * "/projects/working/pbh_population/code/"

############################################################
# binned mass function
function Pm(m, p1, p2, p3, p4)
    """
    binned mass function    
    """
    
    p5 = (1.0 - 2.0p1 - 7.0p2 - 30.0p3 - 40.0p4)/50.0
    if p5 < 0.0
        println("p1 = $p1, p2 = $p2, p3 = $p3, p4 = $p4")
        error("Totol Pm should <= 1")
    end
    
    Pm(m, p1, p2, p3, p4, p5)
end

function Pm(m, p1, p2, p3, p4, p5)
    """
    binned mass function    
    """
        
    if 1.0 <= m < 3.0
        p = p1
    elseif 3.0 <= m < 10.0
        p = p2
    elseif 10.0 <= m < 40.0
        p = p3
    elseif 40.0 <= m < 80.0
        p = p4
    elseif 80.0 <= m <= 130.0
        p = p5
    else
        error("PBH mass is $m, but should be in the range of [1, 130] solar masses.")        
    end
end

function get_p5(p1, p2, p3, p4) 
    p5 = (1.0 - 2.0p1 - 7.0p2 - 30.0p3 - 40.0p4)/50.0
    if p5 < 0.0
        println("p5 < 0")
        p5 = 1e-20
    end
    p5
end

Pm_over_m(m, p1, p2, p3, p4) = Pm(m, p1, p2, p3, p4)/m

hypergeometric_int(t, a, b, z) = exp(-z * t) * t^(a - 1) * (1 + t)^(b - a - 1)

function hypergeometricU(a, b, z)
    int(t, a0, b0, z0) = 1/(1 - t)^2 * hypergeometric_int(t/(1 - t), a0, b0, z0)
    
    function integrand(x, f)
        f[1] = int(x[1], a, b, z)
    end
    result, err = cuhre(integrand, rtol=1e-3)
    result[1]/gamma(a)
end

hyperU(fpbh) = hypergeometricU(21.0/74.0, 0.5, 5.0fpbh/6.0/σM^2.0)

log_fpbhs = -6:0.2:0
log_hyperUs = log10.(hyperU.(10 .^ log_fpbhs))
log_log_hyperU_inter = interpolate((log_fpbhs,), log_hyperUs, Gridded(Linear()))

hyperU_inter(x) = 10^log_log_hyperU_inter(log10(x))

function C_func(fpbh, avm, avm2)
#      gamma(29.0/37.0)/sqrt(pi) = 0.6674098580018847
#     (-74.0/21.0) = -3.5238095238095237
    tmp = (0.6674098580018847 * hyperU_inter(fpbh))^(-3.5238095238095237) - 1.0
    fpbh^2 * avm2/avm^2/σM^2 / tmp
end

function S_factor(m1, m2, p1, p2, p3, p4, fpbh)
    
    avm = 105.0 - 513.0p1 - 3637.5p2 - 3450.0p3 + 2400.0p4
    avm2 = 11233.3 - 56125.3p1 - 392875.0p2 - 428333.0p3 + 149333.0p4
    Nbar = (m1 + m2)/avm * fpbh/(fpbh + σM)
    C = C_func(fpbh, avm, avm2)
    
    #     (-21/74) = -0.28378378378378377
    S1 = 1.42 * exp(-Nbar) * (avm2/avm^2/(Nbar + C) + σM^2/fpbh^2)^(-0.28378378378378377)
    S2 = min(1.0, 9.6e-3 * fpbh^(-0.65) * exp(0.03 * log(fpbh)^2))

    S = S1 * S2
end

diff_S(m1, m2, p1, p2, p3, p4, fpbh) = S_factor(m1, m2, p1, p2, p3, p4, fpbh) * (1 + σM^2/fpbh^2)^(21/74)

function R12(m1, m2, p1, p2, p3, p4, fpbh)
    σeq = 5e-3
    f = 0.85fpbh
    
    p10 = Pm_over_m(m1, p1, p2, p3, p4)
    p20 = Pm_over_m(m2, p1, p2, p3, p4)
    
    3.9e6 * f^2 * (f^2 + σeq^2)^(-21/74) * min(p10, p20) * (p10+p20) * (m1*m2)^(3/37) * (m1+m2)^(36/37)
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

function β_func(p1, p2, p3, p4, fpbh)
    int0(m1, m2) = R12(m1, m2, p1, p2, p3, p4, fpbh) * VT(m1, m2)
    int(m1, m2) = δm^2 * int0(m_min + δm * m1, m_min + δm * m2)
    function integrand(x, f)
        f[1] = int(x[1], x[2])
    end
    result, err = cuhre(integrand, rtol=1e-1)
    0.5 * result[1]
end

function log_pR_func(p1, p2, p3, p4, fpbh, post)
    int0(m1, m2) = R12(m1, m2, p1, p2, p3, p4, fpbh) * post(m1, m2)
    int(m1, m2) = δm^2 * int0(m_min + δm * m1, m_min + δm * m2)
    function integrand(x, f)
        f[1] = int(x[1], x[2])
    end
    result, err = cuhre(integrand, rtol=1e-1)
    log(0.5 * result[1])
end


function logL(log_p1, log_p2, log_p3, log_p4, log_fpbh, posts)
    p1, p2, p3, p4, fpbh = 10.0 .^ (log_p1, log_p2, log_p3, log_p4, log_fpbh)
    if 2.0p1 + 7.0p2 + 30.0p3 + 40.0p4 > 1.0
        return -1e6
    end
    result = - β_func(p1, p2, p3, p4, fpbh) + sum(log_pR_func.(p1, p2, p3, p4, fpbh, posts))

end

function merger_rate(p1, p2, p3, p4, fpbh, m_min=m_min, m_max=m_max)
    δm = m_max - m_min
    int(m1, m2) = δm^2 * R12(m_min + δm * m1, m_min + δm * m2, p1, p2, p3, p4, fpbh)
#     println("p1=$p10, p2=$p20, p3=$p30")
    function integrand(x, f)
        f[1] = int(x[1], x[2])
    end
    result, err = cuhre(integrand, rtol=1e-1)
    0.5 * result[1]
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


# SGWB
const pc = 3.08567758149137e16 # meter; 
const Mpc = 1e6 * pc;
const km = 1e3;
const yr = 365.25 * 24 * 3600; # sec
const lC = 2.99792e8
const z0 = 1e-14;
const Gpc = 1e9 * pc;
const G = 6.67408e-11;

# TT,TE,EE+lowE+lensing+BAO of Table 2 in https://arxiv.org/pdf/1807.06209.pdf
const Om0 = 0.3111 # including neutrinos
const Or0 = 5.4020151370547025e-05
const H0 = 67.66 * km * Mpc^-1 #; (* Hubble constant *)
const OΛ0 = 1 - Om0 - Or0


function νs(i, m1, m2)
    # data from Table 1 in arXiv:0710.2335
    a = [2.9740e-1, 5.9411e-1, 8.4845e-1, 5.0801e-1];
    b = [4.4810e-2, 8.9794e-2, 1.2848e-1, 7.7515e-2];
    c = [9.5560e-2, 1.9111e-1, 2.7299e-1, 2.2369e-2];
    M = m1 + m2
    η = (m1 * m2)/M^2
#     mSun = (1.98855e30)
#     G = 6.67408e-11
#     lC = 2.99792458e8
    # lC^3/(π * G * mSun) = 64622.582946495444
    (a[i] * η^2 + b[i] * η + c[i])/M * 64622.582946495444
end

function dE_over_dν(ν, m1, m2, z)
    νz = (1 + z) * ν
    # (π*G)^(2/3)/3 * mSun^(5/3) = 3.6994212882380265e43
    result = 3.6994212882380265e43 * m1*m2 * (m1+m2)^(-1/3)
    
    ν1, ν2, ν3, ν4 = νs.(1:4, m1, m2)
    
    if 0 <= νz < ν1
        result *= νz^(-1/3)
    elseif ν1 <= νz < ν2
        result *= νz^(2/3) / ν1
    elseif ν2 <= νz <= ν3
        result *= νz^2 * ν4^4 * ν1^(-1) * ν2^(-4/3) * (ν4^2 + 4*(-ν2 + νz)^2)^(-2)
    else
        result *= 0
    end
    
    result
end

Ez(z) = sqrt(Or0*(1+z)^4 + Om0*(1+z)^3 + OΛ0)
Hz(z) = H0 * Ez(z)

function tz(z)
    int0(zz) = 1/Hz(zz)/(1+zz)
    int(zz) = int0(z + zz/(1-zz)) / (1-zz)^2
    function integrand(x, f)
        f[1] = int(x[1])
    end
    result, err = cuhre(integrand, rtol=1e-3, minevals=10^5, maxevals=10^10)
    result[1]
end
   
log_zs = cat(log10(z0):-3, -2.9:0.1:3.6, 4:8, dims=(1));
log_tzs = log10.(tz.(10 .^ log_zs));

log_log_tz_inter = interpolate((log_zs,), log_tzs, Gridded(Linear()));
tz_inter(z) = 10^log_log_tz_inter(log10(z))

t0 = tz_inter(z0)

# function Pm_log(σ, Mc, m)
#     1/m/σ/sqrt(2pi) * exp(- log(m/Mc)^2/2/σ^2)
# end

# function R12_log(m1, m2, σ, Mc, fpbh)
#     σeq = 5e-3
#     f = 0.85fpbh
    
#     p10 = Pm_log(σ, Mc, m1)/m1
#     p20 = Pm_log(σ, Mc, m2)/m2
    
#     3.9e6 * f^2 * (f^2 + σeq^2)^(-21/74) * min(p10, p20) * (p10+p20) * (m1*m2)^(3/37) * (m1+m2)^(36/37)
# end

# function Ωgw_int_log(σ, Mc, fpbh, ν, m1, m2, z)
#     if νs(3, m1, m2) >= ν * (1 + z0)
#         # (Gpc^(-3) * yr^(-1)) / lC^2 /ρc0 = 1.3955693892394134e-75
#         0.5* 1.3955693892394134e-75 * ν * R12_log(m1, m2, σ, Mc, fpbh) * (tz_inter(z)/t0)^(-34/37) / (1+z) / Hz(z) * dE_over_dν(ν, m1, m2, z)
#     else
#         0
#     end
# end

# function Ωgw2_log(σ, Mc, fpbh, ν, m1, m2)
    
#     int(z) = let δz = min((νs(3, m1, m2)/ν - 1 - z0), 1000) # min((νs(3, m1, m2)/ν - 1 - z0), 30)
#         δz * Ωgw_int_log(σ, Mc, fpbh, ν, m1, m2, z0 + δz * z)
#     end
    
# #     int(m1, m2, z) = δm^2 * (νs(3, m1, m2)/ν - 1 - z0) * Ωgw_int(p1, p2, p3, p4, fpbh, ν, m_min + δm * m1, m_min + δm * m2, z0 + (νs(3, m1, m2)/ν - 1 - z0) * z)
    
# #     println(int(0.5))
#     function integrand(x, f)
#         f[1] = int(x[1])
#     end
#     result, err = cuhre(integrand, rtol=1e-3)
#     result[1]
# end

# function Ωgw_log2(σ, Mc, fpbh, ν)
#     m_min = 5
#     δm = 100-5
#     int(m1, m2) = let δz = min(νs(3, m1, m2)/ν - 1 - z0, 1000)
#         δm^2 * Ωgw2_log(σ, Mc, fpbh, ν, m_min + δm*m1, m_min + δm*m2)
#     end
    
# #     int(m1, m2, z) = δm^2 * (νs(3, m1, m2)/ν - 1 - z0) * Ωgw_int(p1, p2, p3, p4, fpbh, ν, m_min + δm * m1, m_min + δm * m2, z0 + (νs(3, m1, m2)/ν - 1 - z0) * z)
    
# #     println(int(0.5, 0.5, 1e-5))
#     function integrand(x, f)
#         f[1] = int(x[1], x[2])
#     end
#     result, err = cuhre(integrand, rtol=1e-1)
#     result[1]
# end

# function Ωgw_log(σ, Mc, fpbh, ν)
#     m_min = 5
#     δm = 100-5
#     int(m1, m2, z) = let δz = min(νs(3, m1, m2)/ν - 1 - z0, 1000)
#         δm^2 * δz * Ωgw_int_log(σ, Mc, fpbh, ν, m_min + δm*m1, m_min + δm*m2, z0 + δz*z)
#     end
    
# #     int(m1, m2, z) = δm^2 * (νs(3, m1, m2)/ν - 1 - z0) * Ωgw_int(p1, p2, p3, p4, fpbh, ν, m_min + δm * m1, m_min + δm * m2, z0 + (νs(3, m1, m2)/ν - 1 - z0) * z)
    
# #     println(int(0.5, 0.5, 1e-5))
#     function integrand(x, f)
#         f[1] = int(x[1], x[2], x[3])
#     end
#     result, err = cuhre(integrand, 3, 1, rtol=1e-3)
#     result[1]
# end

ρc0 = (3.0*H0^2)/(8*pi*G)

function Ωgw_int(p1, p2, p3, p4, fpbh, ν, m1, m2, z)
    if νs(3, m1, m2) >= ν * (1 + z0)
        # (Gpc^(-3) * yr^(-1)) / lC^2 /ρc0 = 1.3955693892394134e-75
        1.3955693892394134e-75 * ν * R12(m1, m2, p1, p2, p3, p4, fpbh) * (tz_inter(z)/t0)^(-34/37) / (1+z) / Hz(z) * dE_over_dν(ν, m1, m2, z)
    else
        0
    end
end

function Ωgw(p1, p2, p3, p4, fpbh, ν)
#     m_min = 5
#     δm = m_max - m_min
    int(m1, m2, z) = let δz = min((νs(3, m1, m2)/ν - 1 - z0), 2000)
        δm^2 * δz * Ωgw_int(p1, p2, p3, p4, fpbh, ν, m_min + δm * m1, m_min + δm * m2, z0 + δz * z)
    end
    
    function integrand(x, f)
        f[1] = int(x[1], x[2], x[3])
    end
    result, err = cuhre(integrand, 3, 1, rtol=1e-4)
    0.5 * result[1]
end