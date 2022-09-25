const yr = 365.25 * 24 * 3600 # seconds

abstract type Cosmology end
struct LCDM <: Cosmology
    h::Float64 # hubble const
    H0::Float64
    t0::Float64

    ΩM::Float64
    ΩΛ::Float64

    HH::Function
    ϕt::Function
    ϕr::Function
    ϕV::Function
    zt::Function
    tz::Function
    rz::Function
    DLz::Function # luminosity distance


    function LCDM(h, ΩM)
        ΩΛ = 1 - ΩM
        H0 = 100h / 3.08567758128e19 # https://www.unitconverters.net/length/megaparsec-to-kilometer.htm

        # eq.(A2) of arxiv:1712.01168
        HH(z) = sqrt(ΩΛ + ΩM * (1 + z)^3)

        # eq.(A5) of arxiv:1712.01168
        function ϕt0(z)
            int0(z) = 1 / HH(z) / (1 + z)
            int(y) = int0(z + y / (1 - y)) / (1 - y)^2
            # int(y) = 1/ (1-y) / (1+(1-y)*z) / HH(y/(1-y) + z)
            function integrand(x, f)
                f[1] = int(x[1])
            end
            result, err = cuhre(integrand, rtol=1e-3)
            result[1]
        end

        # eq.(A6) of arxiv:1712.01168
        function ϕr0(z)
            int0(z) = 1 / HH(z)
            int(y) = z * int0(z * y)
            # int(y) = z / HH(z*y)
            function integrand(x, f)
                f[1] = int(x[1])
            end
            result, err = cuhre(integrand, rtol=1e-3, minevals=10^5)
            result[1]
        end

        logzs = -10:0.01:2
        logϕts = log10.(ϕt0.(10 .^ logzs))
        logts = logϕts .- log10(H0) # eq.(A5) of arxiv:1712.01168
        logϕrs = log10.(ϕr0.(10 .^ logzs))

        writedlm("logzs_logϕts.txt", zip(logzs, logϕts))
        log_log_ϕt_inter = interpolate((logzs,), logϕts, Gridded(Linear()))
        ϕt(z) = 10^log_log_ϕt_inter(log10(z))
        tz(z) = ϕt(z) / H0 # eq.(A5) of arxiv:1712.01168

        log_log_zt_inter = interpolate((reverse(logts),), reverse(logzs), Gridded(Linear()))
        zt(t) = 10^log_log_zt_inter(log10(t))

        writedlm("logzs_logϕrs.txt", zip(logzs, logϕrs))
        log_log_ϕr_inter = interpolate((logzs,), logϕrs, Gridded(Linear()))
        ϕr(z) = 10^log_log_ϕr_inter(log10(z))
        rz(z) = ϕr(z) / H0 # eq.(A6) of arxiv:1712.01168
        DLz(z) = (1 + z) * ϕr(z) / H0

        # eq.(A7) of arxiv:1712.01168
        ϕV(z) = 4π * ϕr(z)^2 / (1 + z)^3 / HH(z)

        t0 = tz(1e-10)

        new(h, H0, t0, ΩM, ΩΛ, HH, ϕt, ϕr, ϕV, zt, tz, rz, DLz)
    end
end


function get_post(event)
    post_file = "/home/czc/projects/working/LIGO_posterior/src_masses_DL_z_posterior/GW$(event)_posterior.txt"
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