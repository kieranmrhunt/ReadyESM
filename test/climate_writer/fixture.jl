function fixture(; FT=Float64, pressure=[100000.,80000.,20000.])
    sigma = FT[.1,.2,.5,.95]
    thick = FT[.15,.20,.4,.25]
    ps = FT.(pressure)
    t = [FT(250 + 3*log(Float64(ps[i])*Float64(sigma[k]))) for i in eachindex(ps), k in eachindex(sigma)]
    u = [FT(-80 + 7*log(Float64(ps[i])*Float64(sigma[k]))) for i in eachindex(ps), k in eachindex(sigma)]
    q = fill(FT(.004),length(ps),length(sigma))
    return t,q,u,ps,sigma,thick,FT(10)
end
