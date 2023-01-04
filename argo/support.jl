using SupportPoints, UnicodePlots, Loess

include("read.jl")

m = 5

sp = supportpoints(temp, m; maxit=100)

# Spaggheti plot of support points
ex = extrema(vec(sp))
plt = lineplot(pressure, sp[:, 1], xlabel="Pressure", ylabel="Temperature", ylim=ex)
for j in 2:m
    lineplot!(plt, pressure, sp[:, j])
end
println(plt)

fn = ["Latitude", "Longitude", "Day"]

for j in 1:m

    d = [norm(c-sp[:, j]) for c in eachcol(temp)]

    plt = lineplot(pressure, sp[:, j], xlabel="Pressure", ylabel="Support $(j)")
    println(plt)

    for k in 1:3
        m = loess(Y[:, k], d)
        xx = range(extrema(Y[:, k])..., length=100)
        yy = predict(m, xx)
        plt = lineplot(xx, yy, xlabel=fn[k], ylabel="Support $(j)", ylim=[0, ceil(maximum(yy))])
        println(plt)
    end

end
