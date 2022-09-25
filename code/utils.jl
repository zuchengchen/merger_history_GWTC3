using BAT
using BenchmarkTools
using Cuba
using DelimitedFiles
using FileIO
using Interpolations
using IntervalSets
using LaTeXStrings
using Memoize
using Revise
using StatsPlots
using Trapz

const yr = 365.25 * 24 * 3600 # seconds
const Mpc_to_km = 3.08567758128e19 # https://www.unitconverters.net/length/megaparsec-to-kilometer.htm
const lc = 2.99792458e8; # speed of light in m/s
const Gpc_in_m = 1e6 * Mpc_to_km;

abstract type Cosmology end
struct LCDM <: Cosmology
    h::Float64 # hubble const
    t0::Float64
    ΩM::Float64
    ΩΛ::Float64

    t_at_z::Function
    z_at_DL::Function # luminosity distance
    dDL_over_dz::Function
    dVC_over_dz::Function

    function LCDM(h, ΩM)
        ΩΛ = 1 - ΩM
        H0 = 100h / Mpc_to_km # 1/second

        # eq.(A2) of arxiv:1712.01168
        Ez(z) = sqrt(ΩΛ + ΩM * (1 + z)^3)

        # eq.(A5) of arxiv:1712.01168
        function ϕt0(z)
            int0(z) = 1 / Ez(z) / (1 + z)
            int(y) = int0(z + y / (1 - y)) / (1 - y)^2
            function integrand(x, f)
                f[1] = int(x[1])
            end
            result, err = cuhre(integrand, rtol=1e-3)
            result[1]
        end

        # eq.(A6) of arxiv:1712.01168
        function ϕr0(z)
            int0(z) = 1 / Ez(z)
            int(y) = z * int0(z * y)
            function integrand(x, f)
                f[1] = int(x[1])
            end
            result, err = cuhre(integrand, rtol=1e-3, minevals=10^5)
            result[1]
        end

        logzs = -4:0.002:1
        zs = 10 .^ logzs
        logϕts = log10.(ϕt0.(zs))
        logϕrs = log10.(ϕr0.(zs))

        logz_logϕt_inter = interpolate((logzs,), logϕts, Gridded(Linear()))
        ϕt(z) = 10^logz_logϕt_inter(log10(z))
        t_at_z(z) = ϕt(z) / H0 # eq.(A5) of arxiv:1712.01168, in second

        logts = log10.(t_at_z.(zs))

        logt_logz_inter = interpolate((reverse(logts),), reverse(logzs), Gridded(Linear()))
        z_at_t(t) = 10^logt_logz_inter(log10(t)) # t in second

        logz_logϕr_inter = interpolate((logzs,), logϕrs, Gridded(Linear()))
        ϕr(z) = 10^logz_logϕr_inter(log10(z))
        rz(z) = ϕr(z) / H0 # eq.(A6) of arxiv:1712.01168
        DC_at_z(z) = lc * ϕr(z) / H0 / Gpc_in_m # in Gpc
        DL_at_z(z) = (1 + z) * DC_at_z(z) # in Gpc

        dDC_over_dz(z) = lc / H0 / Ez(z) / Gpc_in_m # Gpc
        dDL_over_dz(z) = DC_at_z(z) + (1 + z) * dDC_over_dz(z) # Gpc
        dVC_over_dz(z) = 4 * π * DC_at_z(z)^2 * dDC_over_dz(z) # Gpc^3

        logDLs = log10.(DL_at_z.(zs))

        logz_logDL_inter = interpolate((logDLs,), logzs, Gridded(Linear()))
        z_at_DL(DL) = 10^logz_logDL_inter(log10(DL))

        # eq.(A7) of arxiv:1712.01168
        # ϕV(z) = 4π * ϕr(z)^2 / (1 + z)^3 / Ez(z)

        t0 = t_at_z(1e-4)

        new(h, t0, ΩM, ΩΛ, t_at_z, z_at_DL, dDL_over_dz, dVC_over_dz)
    end
end


function get_post(event)
    post_file = "/home/czc/projects/working/LIGO_posterior/masses_DL_det_posterior/GW$(event)_posterior.txt"
    post = readdlm(post_file)
end


function reduce_post(post, new_size)
    N = size(post)[1]
    if N < new_size
        return post
    end
    idx = rand(1:N, new_size)
    post[idx, :]
end

function update_cut(injections, snr_cut, ifar_cut, fraction=nothing)

    println("Selecting injections with SNR $snr_cut and IFAR $ifar_cut yr.")

    idx0 = @. (injections.snr > snr_cut) & (injections.ifar > ifar_cut)
    idx = findall(x -> x == 1, idx0)

    # Sub-sample the selected injections in order to reduce the computational load
    if fraction != nothing
        len = length(idx)
        idx = sample(1:len, float2int(len / fraction))
        ntotal = float2int(injections.ntotal / fraction)
        println("Working with a total of $(ntotal) injections")
    else
        ntotal = injections.ntotal
    end

    m1_det = injections.m1_det[idx]
    m2_det = injections.m2_det[idx]
    DL = injections.DL[idx]
    prior_det = injections.prior_det[idx]
    ifar = injections.ifar[idx]
    snr = injections.snr[idx]

    Injection(m1_det, m2_det, DL, prior_det, ifar, snr, ntotal, injections.Tobs)
end


function detector_frame_to_source_frame(cosmo, m1_det, m2_det, DL)
    z = cosmo.z_at_DL(DL)
    m1_src = m1_det / (1 + z)
    m2_src = m2_det / (1 + z)

    m1_src, m2_src, z
end

detector_to_source_jacobian(cosmo, z) = (1 + z)^2 * cosmo.dDL_over_dz(z) * 1e3 # convert GPc to Mpc

float2int(n) = Int(floor(n))

struct Injection
    m1_det
    m2_det
    DL
    prior_det
    ifar
    snr
    ntotal
    Tobs
end

abstract type MergerRateDensity end
struct Lognormal_R_Density <: MergerRateDensity
    mergerRateDensity1st::Function
    mergerRateDensity2nd::Function
    mergerRateDensity::Function

    function Lognormal_R_Density(mc, σc, f)

        function mergerRateDensity1st(cosmo, i, j, z)
            tz = cosmo.t_at_z(z)
            t0 = cosmo.t0
            (tz / t0)^(-34 / 37) * 210084.52488130186 / (i^2 * j^2 * σc^2) * exp(-(743 * σc^4 + 1369 * log(i / mc)^2 + 1369 * log(j / mc)^2) / (2738 * σc^2)) * (i * j)^(3 / 37) * (i + j)^(36 / 37) * ((exp(σc^2 / 2) * f) / mc)^(53 / 37) * mc^(53 / 37)
        end

        function mergerRateDensity2nd1(i, j)
            tmp1 = 1009.5488113544313 * f^(69 / 37) * i^(6 / 37) * (i + j)^(72 / 37) / j^(68 / 37) / σc^3
            tmp2 = exp(-(-3318 * σc^4 + 1369 * log(j / mc)^2) / (2738 * σc^2))

            int0(e) = exp(-(1369 * log(e / mc)^2 + 1369 * log((-e + i) / mc)^2) / (2738 * σc^2)) / (e^2 * (-e + i)^2)

            vx = range(1e-1, i - 1e-1, length=20)
            M = [int0(e) for e = vx]
            # @show M
            result = trapz((vx), M)

            # int(e) = i * int0(i * e)
            # function integrand(x, f)
            #     f[1] = int(x[1])
            # end
            # result, err = cuhre(integrand, rtol=1e-1)
            tmp1 * tmp2 * result[1]
        end

        function mergerRateDensity2nd(cosmo, i, j, z)
            tz = cosmo.t_at_z(z)
            t0 = cosmo.t0
            (tz / t0)^(-31 / 37) * 0.5 * (mergerRateDensity2nd1(i, j) + mergerRateDensity2nd1(j, i))
        end

        # mergerRateDensity(cosmo, i, j, z) = mergerRateDensity1st(cosmo, i, j, z) + mergerRateDensity2nd(cosmo, i, j, z)
        mergerRateDensity(cosmo, i, j, z) = mergerRateDensity1st(cosmo, i, j, z)

        new(mergerRateDensity1st, mergerRateDensity2nd, mergerRateDensity)
    end
end

function R_over_pdraw(cosmo, R_density, m1_src, m2_src, z, pdraw)
    R = R_density.mergerRateDensity(cosmo, m1_src, m2_src, z)
    jacobian_term = detector_to_source_jacobian(cosmo, z)

    source_to_detector_factor(cosmo, z) * R / pdraw / jacobian_term
end

function cal_Nexp(injections, cosmo, R_density)
    m1_det = injections.m1_det
    m2_det = injections.m2_det
    DL = injections.DL
    prior_det = injections.prior_det
    ntotal = injections.ntotal
    Tobs = injections.Tobs

    len = length(m1_det)

    R_over_pdraws = zeros(len)
    for i in 1:len
        m1_src, m2_src, z = detector_frame_to_source_frame(cosmo, m1_det[i], m2_det[i], DL[i])
        R_over_pdraws[i] = R_over_pdraw(cosmo, R_density, m1_src, m2_src, z, prior_det[i])
    end

    sum(R_over_pdraws) / ntotal * Tobs
end

function lnLike_single_event(cosmo, R_density, post)
    result = 0.0
    for p in eachrow(post)
        m1_det, m2_det, DL = p
        m1_src, m2_src, z = detector_frame_to_source_frame(cosmo, m1_det, m2_det, DL)
        R = R_density.mergerRateDensity(cosmo, m1_src, m2_src, z)
        jacobian_term = detector_to_source_jacobian(cosmo, z)
        result += source_to_detector_factor(cosmo, z) * R / DL^2 / jacobian_term * 1e-6 # DL from Gpc to Mpc
    end
    log(result / size(post)[1])
end

function lnLike_events(cosmo, R_density, posts)
    lnLs = [lnLike_single_event(cosmo, R_density, post) for post in posts]
    sum(lnLs)
end

function lnL0(injections, cosmo, R_density, posts)
    Tobs = injections.Tobs
    Nobs = length(posts)
    Nexp = cal_Nexp(injections, cosmo, R_density)

    Nobs * log(Tobs) - Nexp + lnLike_events(cosmo, R_density, posts)
end

function merger_rate(cosmo, R_density, z)
    m_min = 1
    m_max = 150
    δm = m_max - m_min
    int(m1, m2) = δm^2 * R_density.mergerRateDensity(cosmo, m_min + δm * m1, m_min + δm * m2, z)
    #     println("p1=$p10, p2=$p20, p3=$p30")
    function integrand(x, f)
        f[1] = int(x[1], x[2])
    end
    result, err = cuhre(integrand, rtol=1e-3)
    0.5 * result[1]
end

function merger_rate_1st(cosmo, R_density, z)
    m_min = 1
    m_max = 150
    δm = m_max - m_min
    int(m1, m2) = δm^2 * R_density.mergerRateDensity1st(cosmo, m_min + δm * m1, m_min + δm * m2, z)
    #     println("p1=$p10, p2=$p20, p3=$p30")
    function integrand(x, f)
        f[1] = int(x[1], x[2])
    end
    result, err = cuhre(integrand, rtol=1e-3)
    0.5 * result[1]
end

"""
    merger_rate_2nd(cosmo, R_density, z)

TBW
"""
function merger_rate_2nd(cosmo, R_density, z)
    m_min = 1
    m_max = 150
    δm = m_max - m_min
    int(m1, m2) = δm^2 * R_density.mergerRateDensity2nd(cosmo, m_min + δm * m1, m_min + δm * m2, z)
    #     println("p1=$p10, p2=$p20, p3=$p30")
    function integrand(x, f)
        f[1] = int(x[1], x[2])
    end
    result, err = cuhre(integrand, rtol=1e-3)
    0.5 * result[1]
end

function source_to_detector_factor(cosmo, z)
    cosmo.dVC_over_dz(z) / (1 + z)
end