using DataFrames
using SmoothingKernels
using Gadfly

if length(ARGS)<4
    println("julia comparison.jl run_cnt bins max_time bandwidth")
    println("julia comparison.jl 1000000 2000 20 0.1")
    exit()
end
run_cnt=int(ARGS[1])
bin_cnt=int(ARGS[2])
max_time=float(ARGS[3])
bandwidth=float(ARGS[4])

function smoothed(data::Array{Float64,1}, bin_cnt::Int, max_time::Float64,
        bandwidth::Float64)
    row_cnt=size(data)[1]
    cumulative=zeros(Float64, bin_cnt)
    times=zeros(Float64, bin_cnt)
    lambda=1/bandwidth
    kernel=SmoothingKernels.kernels[:epanechnikov]
    for i in 1:bin_cnt
        dx=(i-1)*1.1*max_time/bin_cnt
        cumulative[i]=sum(lambda*kernel(lambda*(data-dx)))
        times[i]=dx
    end
    cumulative/=row_cnt
    (cumulative, times)
end


function simulated_distribution(runs)
    chart=zeros(Float64, (runs, 2))
    const two_floats=r"^([\d\.eE\-\+]+)\t([\d\.eE\-\+]+)"
    open(`./individual --runcnt $(runs)`) do f
        i=0
        for l in eachline(f)
            m=match(two_floats, l)
            if m!=nothing && length(m.captures)==2
                i+=1
                chart[i,:]=Float64[float(x) for x in m.captures]'
            end
        end
        if i!=runs
            error("missed entries. found ", i, " out of ", runs)
        end
    end
    chart
end

function check_gamma(samples::Array{Float64,1}, alpha, beta)
    println("mean ", mean(samples))
    println("expected mean ", alpha/beta)
    println("variance ", var(samples))
    println("expected variance ", alpha/(beta^2))
end

function plot_two(runs::Int, bin_cnt::Int, max_time, bandwidth::Float64)
    samples=simulated_distribution(runs)
    dist_cnt=size(samples)[2]
    binned=Array(Any,0)
    for didx in 1:dist_cnt
        push!(binned, smoothed(samples[:,didx], bin_cnt, max_time, bandwidth))
    end

    check_gamma(samples[:,2], 3.969, 1/1.107)

    df=DataFrame(Times=binned[1][2], Latent=binned[1][1],
        Recovery=binned[2][1])
    writetable("out.csv", df)
    title="Latent and Recovery Periods"
    myplot=plot(df, layer(x="Times", y="Latent", Geom.line,
            Theme(default_color=color("blue"))),
        layer(x="Times", y="Recovery", Geom.line,
            Theme(default_color=color("green"))))
    filename=join(matchall(r"[A-Za-z]", title))
    draw(PDF("$(filename).pdf", 20cm, 15cm), myplot)
end


plot_two(run_cnt, bin_cnt, max_time, bandwidth)
