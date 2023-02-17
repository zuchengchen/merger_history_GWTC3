# modify model
model = "GWTC2_powerlaw_heavy_mass"

include("utils.jl")

# modify post_events
# posts = [get_post_func(event) for event in GWTC2_events_heavy_mass];
# save("backup/posts_$(model).jld2", Dict("posts" => posts))
posts = load("backup/posts_$(model).jld2", "posts");

fpbh = 0.0025693801260000975 # fpbh
log_fpbh = -3
M = 8.896735628637506
α = 1.9108320383237213
i = 10
j = 30

function mergerRateDensity1st(i, j, M, α, fpbh)
    4.884e7* fpbh^(53/37) * (i + j)^(36/37) * (i/M)^(-α) * (j/M)^(-α) * α^3 / (i*j)^(34/37) / M^(21/37) / (M*α/(-1 + α))^(
 53/37) / (21 + 37α)
end

function mergerRateDensity2nd0(i, j, M, α, fpbh)
    if i > 2*M
        5.883e5 * fpbh^(69/37) * i^(-2α) * j^(-α) * (i+j)^(72/37) * M^(-3+3α) * α^4 * beta_func2(M/i, 1-M/i, -α, -α) * ((α-1)/α)^(69/37) / (i*j)^(31/37) / (42 + 37α)
    else
        0
    end
end

mergerRateDensity2nd(i, j, M, α, f) = 0.5 * (mergerRateDensity2nd0(i, j, M, α, f) + mergerRateDensity2nd0(j, i, M, α, f))

# function mergerRateDensity(i, j, mc, σc, f)
#     if i >= j
#         mergerRateDensity1st(i, j, mc, σc, f) + mergerRateDensity2nd(i, j, mc, σc, f)
#     else
#         mergerRateDensity(j, i, mc, σc, f)
#     end
# end

mergerRateDensity(i, j, M, α, f) = mergerRateDensity1st(i, j, M, α, f) + mergerRateDensity2nd(i, j, M, α, f)

function merger_rate(M, α, f) 
    int(m1, m2) = δm^2 * mergerRateDensity(m_min + δm*m1, m_min + δm*m2, M, α, f) 
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

function β_func(M, α, f)
    int0(m1, m2) = mergerRateDensity(m1, m2, M, α, f) * VT(m1, m2)
    int(m1, m2) = δm^2 * int0(m_min + δm * m1, m_min + δm * m2)
    function integrand(x, f)
        f[1] = int(x[1], x[2])
    end
    result, err = cuhre(integrand, rtol=1e-1)
    0.5 * result[1]
end

@btime β_func(M, α, fpbh)

function log_pR_func(M, α, f, post)
    int0(m1, m2) = mergerRateDensity(m1, m2, M, α, f) * post(m1, m2)
    int(m1, m2) = δm^2 * int0(m_min + δm * m1, m_min + δm * m2)
    function integrand(x, f)
        f[1] = int(x[1], x[2])
    end
    result, err = cuhre(integrand, rtol=1e-1)
    log(0.5 * abs(result[1]))
end

@time log_pR_func(M, α, fpbh, posts[8])

function logL(M, α, log_fpbh, posts)
    - β_func(M, α, 10.0^log_fpbh) + sum(log_pR_func.(M, α, 10.0^log_fpbh, posts))
end

@time logL(1, 1.01, log_fpbh, posts)

logL(1, 4, 0, posts)

@btime logL(M, α, log_fpbh, posts)

likelihood = ps -> LogDVal(logL(ps.M, ps.α, ps.log_fpbh, posts))

prior = BAT.NamedTupleDist(
    M = (1)..(20),
    α = (1.01)..(4.0),
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