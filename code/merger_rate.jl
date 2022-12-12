using Trapz, SpecialFunctions

const m_min = 1e-1
const m_max = 1e3

# log-normal mass function 
function mergerRateDensity1st0_log(mc, σc, log_fpbh, i, j)
    fpbh = 10^log_fpbh
    210084.52488130186 / (i^2 * j^2 * σc^2) * exp(-(743 * σc^4 + 1369 * log(i / mc)^2 + 1369 * log(j / mc)^2) / (2738 * σc^2)) * (i * j)^(3 / 37) * (i + j)^(36 / 37) * ((exp(σc^2 / 2) * fpbh) / mc)^(53 / 37) * mc^(53 / 37)
end

function mergerRateDensity2nd00_log(mc, σc, log_fpbh, i, j)
    fpbh = 10^log_fpbh
    tmp1 = 1009.5488113544313 * fpbh^(69 / 37) * i^(6 / 37) * (i + j)^(72 / 37) / j^(68 / 37) / σc^3
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

function mergerRateDensity2nd0_log(mc, σc, log_fpbh, i, j)
    0.5 * (mergerRateDensity2nd00_log(mc, σc, log_fpbh, i, j) + mergerRateDensity2nd00_log(mc, σc, log_fpbh, j, i))
end

mergerRateDensity1st_log(mc, σc, log_fpbh, i, j) = mergerRateDensity1st0_log.(mc, σc, log_fpbh, i, j)
mergerRateDensity2nd_log(mc, σc, log_fpbh, i, j) = mergerRateDensity2nd0_log.(mc, σc, log_fpbh, i, j)


# power-law mass function
function mergerRateDensity1st0_power(α, M, log_fpbh, i, j)
    if (i > M) & (j > M)
        fpbh = 10^log_fpbh
        numerator = 4.884e7 * (i + j)^(36 / 37) * (i / M)^(-α) * (j / M)^(-α) * (fpbh * (α - 1.0) / M / α)^(53 / 37) * α^3.0
        denominator = (i * j)^(34 / 37) * M^(21.0 / 37.0) * (21.0 + 37.0 * α)
        numerator / denominator
    else
        0.0
    end
end

function mergerRateDensity2nd00_power(α, M, log_fpbh, i, j)
    if (i > 2 * M) & (j > M)
        fpbh = 10^log_fpbh
        t1 = 588300.0 * i^(6 / 37) * (i + j)^(72 / 37) * (fpbh * (α - 1.0) / M / α)^(69 / 37) * α^4 * M^(3 * α)
        t2 = j^(31 / 37) * M^(42 / 37) * (42 + 37 * α)

        int0(e) = (e * (-e + i) * j)^(-α) / e / (-e + i)

        ve = range(M, i - M, length=20)
        result = trapz((ve), int0.(ve))

        t1 / t2 * result[1]
    else
        0.0
    end
end

function mergerRateDensity2nd0_power(α, M, log_fpbh, i, j)
    0.5 * (mergerRateDensity2nd00_power(α, M, log_fpbh, i, j) + mergerRateDensity2nd00_power(α, M, log_fpbh, j, i))
end

mergerRateDensity1st_power(α, M, log_fpbh, i, j) = mergerRateDensity1st0_power.(α, M, log_fpbh, i, j)
mergerRateDensity2nd_power(α, M, log_fpbh, i, j) = mergerRateDensity2nd0_power.(α, M, log_fpbh, i, j)


# critical collapse
function mergerRateDensity1st0_CC(α, Mf, log_fpbh, i, j)
    fpbh = 10^log_fpbh
    t1 = (i * j)^(-34 / 37) * 1.32e6 * exp(-(i / Mf)^α - (j / Mf)^α) * i^α * j^α * (i + j)^(36 / 37) * Mf^(-(21 / 37) - 2 * α) * α^2
    t2 = (fpbh * α / Mf / gamma(1 / α))^(53 / 37) * (gamma(1 - 21 / 37 / α, (10.0 * Mf)^(-α)) - gamma(1 - 21 / 37 / α, (1e-3 * Mf)^(-α)))
    t1 * t2
end

function mergerRateDensity2nd00_CC(α, Mf, log_fpbh, i, j)
    fpbh = 10^log_fpbh

    t1 = 15900.0 * i^(6 / 37) * j^(-(31 / 37) + α) * (i + j)^(72 / 37) * Mf^(-(42 / 37) - 3 * α) * α^3
    t2 = (fpbh * α / Mf / gamma(1 / α))^(69 / 37) * (gamma(1 - 42 / 37 / α, (10.0 * Mf)^-α) - gamma(1 - 42 / 37 / α, (1e-3 * Mf)^(-α)))

    int1(e) = e^(-1 + α) * exp(-(e / Mf)^α - ((-e + i) / Mf)^α - (j / Mf)^α) * (-e + i)^(-1 + α)

    ve = range(m_min, i - 1e-5, length=20)
    Ie = trapz((ve), int1.(ve))

    t1 * t2 * Ie
end

function mergerRateDensity2nd0_CC(α, Mf, log_fpbh, i, j)
    0.5 * (mergerRateDensity2nd00_CC(α, Mf, log_fpbh, i, j) + mergerRateDensity2nd00_CC(α, Mf, log_fpbh, j, i))
end

mergerRateDensity1st_CC(α, Mf, log_fpbh, i, j) = mergerRateDensity1st0_CC.(α, Mf, log_fpbh, i, j)
mergerRateDensity2nd_CC(α, Mf, log_fpbh, i, j) = mergerRateDensity2nd0_CC.(α, Mf, log_fpbh, i, j)

# broken power-law
function piecewise(i, ms, α1, α2)
    (i < ms) ? (i / ms)^α1 : (i / ms)^(-α2)
end


function mergerRateDensity1st0_bpower(ms, α1, α2, log_fpbh, i, j)
    fpbh = 10.0^log_fpbh
    t1 = 3.5832299e7 * 10^(-α1 - 3 * α2) * fpbh * (i + j)^(36 / 37) * ms^(-58 / 37 - α1) * α1^2 * (1.0 + α1) * (-1.0 + α2)
    t2 = (fpbh * (1 + α1) * (-1 + α2) / ms / α1 / α2)^(16 / 37) * α2^2
    t31 = 10^α1 * ms^(21 / 37 + α1 + α2) * (-0.567568 + α1) + 10^(α1 + 3 * α2) * ms^α1 * (-50.4316 * α1 - 50.4316 * α2)
    t32 = 1000.0^α2 * ms^(21 / 37) * (105.752 + 186.325 * α2)
    t4 = (i * j)^(34 / 37) * (-21 + 37 * α1) * (α1 + α2)^3 * (21 + 37 * α2)
    -t1 * t2 * (t31 + t32) * piecewise(i, ms, α1, α2) * piecewise(j, ms, α1, α2) / t4
end

function mergerRateDensity2nd00_bpower(ms, α1, α2, log_fpbh, i, j)

    fpbh = 10.0^log_fpbh

    t1 = 8558.45 * 10^(-α1 - 3 * α2) * fpbh * i^(6 / 37) * (i + j)^(72 / 37) * ms^(-79 / 37 - α1) * α1^3 * (1 + α1) * (-1 + α2)
    t2 = (fpbh * (1 + α1) * (-1 + α2) / ms / α1 / α2)^(32 / 37) * α2^3
    t31 = 10^α1 * ms^(42 / 37 + α1 + α2) * (-1.13514 + α1)
    t32 = 10^(α1 + 3 * α2) * ms^α1 * (-2543.35 * α1 - 2543.35 * α2)
    t33 = 1000^α2 * ms^(42 / 37) * (39408.3 + 34716.9 * α2)
    t4 = j^(31 / 37) * (-42.0 + 37.0 * α1) * (α1 + α2)^4 * (42.0 + 37.0 * α2)

    int(e) = piecewise(e, ms, α1, α2) * piecewise(i - e, ms, α1, α2) / e / (-e + i)

    ve = range(m_min, i - 1e-5, length=20)
    Ie = trapz((ve), int.(ve))

    -t1 * t2 * (t31 + t32 + t33) * piecewise(j, ms, α1, α2) / t4 * Ie
end

function mergerRateDensity2nd0_bpower(ms, α1, α2, log_fpbh, i, j)
    0.5 * (mergerRateDensity2nd00_bpower(ms, α1, α2, log_fpbh, i, j) + mergerRateDensity2nd00_bpower(ms, α1, α2, log_fpbh, j, i))
end

mergerRateDensity1st_bpower(ms, α1, α2, log_fpbh, i, j) = mergerRateDensity1st0_bpower.(ms, α1, α2, log_fpbh, i, j)
mergerRateDensity2nd_bpower(ms, α1, α2, log_fpbh, i, j) = mergerRateDensity2nd0_bpower.(ms, α1, α2, log_fpbh, i, j)