# modify model
model = "GWTC2_lognormal_heavy_mass"

include("utils.jl")

# modify post_events
posts = [get_post_func(event) for event in GWTC2_events_heavy_mass];
save("backup/posts_$(model).jld2", Dict("posts" => posts))
posts = load("backup/posts_$(model).jld2", "posts");

f = 0.0025693801260000975 # fpbh
log_fpbh = 1e-3
mc = 8.896735628637506
σc = 0.9108320383237213
i = 10
j = 30

function mergerRateDensity1st(i, j, mc, σc, f)
    210084.52488130186/(i^2 * j^2 * σc^2) * exp(-(743*σc^4 + 1369*log(i/mc)^2 + 1369*log(j/mc)^2)/(2738*σc^2)) * (i*j)^(3/37) * (i+j)^(36/37) * ((exp(σc^2/2)*f)/mc)^(53/37) * mc^(53/37)
end

function mergerRateDensity2nd1(i, j, mc, σc, f)
    tmp1 = 1009.5488113544313 * f^(69/37) * i^(6/37) * (i+j)^(72/37) / j^(68/37) / σc^3 
    tmp2 = exp(-(-3318*σc^4 + 1369*log(j/mc)^2)/(2738*σc^2))
    
    int0(e) = exp(-(1369*log(e/mc)^2 + 1369*log((-e+i)/mc)^2)/(2738*σc^2)) / (e^2 * (-e + i)^2)
    int(e) = i * int0(i*e)
    function integrand(x, f)
        f[1] = int(x[1])
    end
    result, err = cuhre(integrand, rtol=1e-1)
    tmp1 * tmp2 * result[1]
end

mergerRateDensity2nd(i, j, mc, σc, f) = 0.5 * (mergerRateDensity2nd1(i, j, mc, σc, f) + mergerRateDensity2nd1(j, i, mc, σc, f))

# function mergerRateDensity(i, j, mc, σc, f)
#     if i >= j
#         mergerRateDensity1st(i, j, mc, σc, f) + mergerRateDensity2nd(i, j, mc, σc, f)
#     else
#         mergerRateDensity(j, i, mc, σc, f)
#     end
# end

mergerRateDensity(i, j, mc, σc, f) = mergerRateDensity1st(i, j, mc, σc, f) + mergerRateDensity2nd(i, j, mc, σc, f)

function merger_rate(mc, σc, f) 
    δm = m_max - m_min
    int(m1, m2) = δm^2 * mergerRateDensity(m_min + δm*m1, m_min + δm*m2, mc, σc, f) 
#     println("p1=$p10, p2=$p20, p3=$p30")
    function integrand(x, f)
        f[1] = int(x[1], x[2])
    end
    result, err = cuhre(integrand, rtol=1e-3)
    0.5*result[1]
end

# @time mergerRateDensity(i, j, mc, σc, f)
# mergerRateDensity1st(i, j, mc, σc, f)
# @time merger_rate(mc, σc, f)
# mergerRateDensity2nd(i, j, mc, σc, f)
# @time mergerRateDensity2nd(i, j, mc, σc, f)
# @time mergerRateDensity2nd1(i, j, mc, σc, f)

# function β_func(mc, σc, f)
#     int0(m1, m2) = mergerRateDensity(m1, m2, mc, σc, f) * VT(m1, m2)
#     int(m1, m2) = δm^2 * int0(m_min + δm * m1, m_min + δm * m2)
#     function integrand(x, f)
#         f[1] = int(x[1], x[2])
#     end
#     result, err = cuhre(integrand, rtol=1e-1)
#     0.5 * result[1]
# end

function β_func_1st(mc, σc, f)
    int0(m1, m2) = mergerRateDensity1st(m1, m2, mc, σc, f) * VT(m1, m2)
    int(m1, m2) = δm^2 * int0(m_min + δm * m1, m_min + δm * m2)
    function integrand(x, f)
        f[1] = int(x[1], x[2])
    end
    result, err = cuhre(integrand, rtol=1e-1)
    0.5 * result[1]
end

function β_func_2nd(mc, σc, f)
    
    int00(i, j, e) = 1009.5488113544313 * f^(69/37) * i^(6/37) * (i+j)^(72/37) / j^(68/37) / σc^3 * exp(-(-3318*σc^4 + 1369*log(j/mc)^2)/(2738*σc^2)) * exp(-(1369*log(e/mc)^2 + 1369*log((-e+i)/mc)^2)/(2738*σc^2)) / (e^2 * (-e + i)^2) * VT(i, j)
#     int0(i, j, e) = 0.5 * (int00(i, j, e) + int00(j, i, e))
    int(i, j, e) = δm^2 * (m_min + δm * i) * int00(m_min + δm * i, m_min + δm * j, (m_min + δm * i)*e)
    function integrand(x, f)
        f[1] = int(x[1], x[2], x[3])
    end
    result, err = cuhre(integrand, 3, 1, rtol=1e-1)
    0.5 * result[1]
end

β_func(mc, σc, f) = β_func_1st(mc, σc, f) + β_func_2nd(mc, σc, f)

@btime β_func(mc, σc, f)

# function log_pR_func(mc, σc, f, post)
#     int0(m1, m2) = mergerRateDensity(m1, m2, mc, σc, f) * post(m1, m2)
#     int(m1, m2) = δm^2 * int0(m_min + δm * m1, m_min + δm * m2)
#     function integrand(x, f)
#         f[1] = int(x[1], x[2])
#     end
#     result, err = cuhre(integrand, rtol=1e-1)
#     log(0.5 * abs(result[1]))
# end

function log_pR_func_1st(mc, σc, f, post)
    int0(m1, m2) = mergerRateDensity1st(m1, m2, mc, σc, f) * post(m1, m2)
    int(m1, m2) = δm^2 * int0(m_min + δm * m1, m_min + δm * m2)
    function integrand(x, f)
        f[1] = int(x[1], x[2])
    end
    result, err = cuhre(integrand, rtol=1e-1)
    log(0.5 * abs(result[1]))
end

function log_pR_func_2nd(mc, σc, f, post)
    
    int00(i, j, e) = 1009.5488113544313 * f^(69/37) * i^(6/37) * (i+j)^(72/37) / j^(68/37) / σc^3 * exp(-(-3318*σc^4 + 1369*log(j/mc)^2)/(2738*σc^2)) * exp(-(1369*log(e/mc)^2 + 1369*log((-e+i)/mc)^2)/(2738*σc^2)) / (e^2 * (-e + i)^2) * post(i, j)
#     int0(i, j, e) = 0.5 * (int00(i, j, e) + int00(j, i, e))
    int(i, j, e) = δm^2 * (m_min + δm * i) * int00(m_min + δm * i, m_min + δm * j, (m_min + δm * i)*e)
    function integrand(x, f)
        f[1] = int(x[1], x[2], x[3])
    end
    result, err = cuhre(integrand, 3, 1, rtol=1e-1)
    result[1]
end

log_pR_func(mc, σc, f, post) = log_pR_func_1st(mc, σc, f, post) + log_pR_func_2nd(mc, σc, f, post)

@time log_pR_func_1st(mc, σc, f, posts[8])

@time log_pR_func_2nd(mc, σc, f, posts[8])

@time log_pR_func(mc, σc, f, posts[8])

function logL(mc, σc, log_fpbh, posts)
    - β_func(mc, σc, 10.0^log_fpbh) + sum(log_pR_func.(mc, σc, 10.0^log_fpbh, posts))
end

@btime logL(mc, σc, -3, posts)

likelihood = ps -> LogDVal(logL(ps.mc, ps.σc, ps.log_fpbh, posts))

prior = BAT.NamedTupleDist(
    mc = (1)..(30),
    σc = (0.1)..(2.0),
    log_fpbh = (-4)..(-0.0)
)

posterior = PosteriorDensity(likelihood, prior);

@time begin
println("Start sampling.")

burnin = MCMCMultiCycleBurnin(max_ncycles=1000)
samples, chains = bat_sample(posterior, MCMCSampling(mcalg=MetropolisHastings(), nsteps=10^5, nchains=4, burnin = burnin));
    
println(" ")
println("Finish sampling.")
end

save("backup/samples_$model.jld2", Dict("samples" => samples, "chains" => chains))