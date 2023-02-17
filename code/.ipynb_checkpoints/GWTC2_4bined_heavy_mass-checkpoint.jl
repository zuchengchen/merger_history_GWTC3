# modify model
model = "GWTC2_4bined_heavy_mass"

include("utils.jl")

# modify post_events
# posts = [get_post_func(event) for event in GWTC2_events_heavy_mass];
# save("backup/posts_$(model).jld2", Dict("posts" => posts))
posts = load("backup/posts_$(model).jld2", "posts");

fpbh = 1e-3 # fpbh
log_fpbh = -3
p1, p2, p3 = 0.02, 0.01, 0.005
log_p1, log_p2, log_p3 = log10.((p1, p2, p3))
i = 10
j = 30
0.02 - 0.18p1 - 0.6p2 - 0.8p3

function Pm(j, p1, p2, p3)
    if 1 <= j <= 10
        p1
    elseif 10 < j <= 40
        p2
    elseif 40 < j <= 80
        p3
    elseif 80 < j <= 130
        0.02 - 0.18p1 - 0.6p2 - 0.8p3
    else
        0
    end
end

Pm(3, 0.1, 0.2, 0.3)

function mergerRateDensity1st(i, j, p1, p2, p3, fpbh)
    
    t1 = fpbh^(53/37) * (i+j)^(36/37) / (0.485508 + 110.76p1 + 54.7495p2 + 15.237p3)^(21/37)  
    t2 = 3873.23 + 7.01778e6p1 + 1.30948e6p2 + 232650p3
    t3 = (i*j)^(34/37) * (0.00438343 + p1 + 0.494309p2 + 0.137569p3)    
    
    t1 * t2 / t3 * Pm(i, p1, p2, p3) * Pm(j, p1, p2, p3)
    
end

function mergerRateDensity2nd0(i, j, p1, p2, p3, fpbh)
#     println("$i  $j")
    denominator = i * j^(31/37)
    numerator = if 2 < i <= 11
        2p1^2 * log(i-1)
    elseif 11 < i <= 20
        2p1 * (p1 * log(10/(i-10)) + p2 * log(0.1 * (i-10) * (i-1)))
    elseif 20 < i <= 41
        p2 * (-2p2 * log(10/(i-10)) + 2p1 * log(10(i-1)/(i-10)))
    elseif 41 < i <= 50
        -2p2^2 * log(10/(i-10)) + 2p1 * p3 * log(0.025(i-40) * (i-1)) + 2p1 * p2 * log(400/(400 + (i-50)i))
    elseif 50 < i <= 80
        2(p2^2 * log(40/(i-40)) + p2 * p3 * log(0.0025(i-40) * (i-10)) + p1 * p3 * log(10(i-1)/(i-10)))
    elseif 80 < i <= 81
        p3 * (-2p3 * log(40/(i-40)) + 2p2 * log(4(i-10)/(i-40)) + 2p1 * log(10(i-1)/(i-10)))
    elseif 81 < i <= 90
        p3 * (-2p3 * log(40/(i-40)) + 2p2 * log(4(i-10)/(i-40)) + 2p1 * log(800/(800 + (i-90)i))) + p1 * (-0.04 + 0.36p1 + 1.2p2 + 1.6p3) * log(80/(80 + (i-81)i))
    elseif 90 < i <= 120
        16.1418p2 * p3 - 2p3^2 * log(40/(i-40)) - 2p2*p3 * log((i-80) * (i-40))
            + (-0.04 + 0.36p1 + 1.2p2 + 1.6p3) * (p1 * log(0.1(i-10)/(i-1)) + p2 * log(800/(800 + (i-90)i)))
    elseif 120 < i <= 130
        2p3^2 * log(80/(i-80)) 
            + (-0.04 + 0.36p1 + 1.2p2 + 1.6p3) * (
                p2 * log(0.25(i-40)/(i-10)) + p1 * log(0.1(i-10)/(i-1)) + p3 * log(3200/(3200 + (i-120)i))
                )
    else
        0
    end
    
    t1 = 4485.672695128149fpbh^(69/37) * i^(6/37) * (i+j)^(72/37)
    t2 = 0.00006325682443982787 + p1 + 0.06080232325374703p2 + 0.006400685352784427p3
    t3 = (0.4855078157817008 + 110.75968430766699p1 +54.749483582543505p2 + 15.237046396729234p3)^(5/37) 
    t4 = (0.004383434449244752 + p1 + 0.49430877240911064p2 + 0.13756852497343655p3)^2
    
    
    t1 * t2 / t3 / t4 * Pm(j, p1, p2, p3)  * numerator/denominator
    
end

mergerRateDensity2nd(i, j, p1, p2, p3, fpbh) = 0.5 * (mergerRateDensity2nd0(i, j, p1, p2, p3, fpbh) + mergerRateDensity2nd0(i, j, p1, p2, p3, fpbh)) 

mergerRateDensity(i, j, p1, p2, p3, fpbh) = mergerRateDensity1st(i, j, p1, p2, p3, fpbh) + mergerRateDensity2nd(i, j, p1, p2, p3, fpbh)

function merger_rate(p1, p2, p3, fpbh)
    int(m1, m2) = δm^2 * mergerRateDensity(m_min + δm*m1, m_min + δm*m2, p1, p2, p3, fpbh) 
#     println("p1=$p10, p2=$p20, p3=$p30")
    function integrand(x, f)
        f[1] = int(x[1], x[2])
    end
    result, err = cuhre(integrand, rtol=1e-3)
    0.5*result[1]
end

merger_rate(p1, p2, p3, fpbh)

# @time mergerRateDensity(i, j, mc, σc, f)
# mergerRateDensity1st(i, j, mc, σc, f)
# @time merger_rate(mc, σc, f)
# mergerRateDensity2nd(i, j, mc, σc, f)
# @time mergerRateDensity2nd(i, j, mc, σc, f)
# @time mergerRateDensity2nd1(i, j, mc, σc, f)

function β_func(p1, p2, p3, fpbh)
    int0(m1, m2) = mergerRateDensity(m1, m2, p1, p2, p3, fpbh) * VT(m1, m2)
    int(m1, m2) = δm^2 * int0(m_min + δm * m1, m_min + δm * m2)
    function integrand(x, f)
        f[1] = int(x[1], x[2])
    end
    result, err = cuhre(integrand, rtol=1e-1)
    0.5 * result[1]
end

@btime β_func(p1, p2, p3, fpbh)

function log_pR_func(p1, p2, p3, fpbh, post)
    int0(m1, m2) = mergerRateDensity(m1, m2, p1, p2, p3, fpbh) * post(m1, m2)
    int(m1, m2) = δm^2 * int0(m_min + δm * m1, m_min + δm * m2)
    function integrand(x, f)
        f[1] = int(x[1], x[2])
    end
    result, err = cuhre(integrand, rtol=1e-1)
    log(0.5 * abs(result[1]))
end

@btime log_pR_func(p1, p2, p3, fpbh, posts[8])

function logL(log_p1, log_p2, log_p3, log_fpbh, posts)
    p1, p2, p3, fpbh = 10.0 .^ (log_p1, log_p2, log_p3, log_fpbh)
    if 0.02 - 0.18p1 - 0.6p2 - 0.8p3 < 0.0
        return -1e6
    end
    result = - β_func(p1, p2, p3, fpbh) + sum(log_pR_func.(p1, p2, p3, fpbh, posts))
end

@time logL(log_p1, log_p2, log_p3, log_fpbh, posts)

@btime logL(log_p1, log_p2, log_p3, log_fpbh, posts)

likelihood = ps -> LogDVal(logL(ps.log_p1, ps.log_p2, ps.log_p3, ps.log_fpbh, posts))

prior = BAT.NamedTupleDist(
    log_p1 = (-2)..(1/9),
    log_p2 = (-4)..(1/30),
    log_p3 = (-6)..(1/40),
    log_fpbh = (-4)..(-0.0)
)

posterior = PosteriorDensity(likelihood, prior);

@time begin
println("Start sampling.")
flush(stdout)

burnin = MCMCMultiCycleBurnin(max_ncycles=1000)
samples, chains = bat_sample(posterior, MCMCSampling(mcalg=MetropolisHastings(), nsteps=10^5, nchains=4, burnin = burnin));
    
println(" ")
println("Finish sampling.")
end

save("backup/samples_$model.jld2", Dict("samples" => samples, "chains" => chains))