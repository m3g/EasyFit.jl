
include("../src/EasyFit.jl")

x = sort(rand(10)) .- 0.5
y = sort(rand(10)) .- 0.5

plot(layout=(3,2))
plot!(framestyle=:box,grid=false)

sp=1
plot!(title="Linear fit",subplot=sp)
scatter!(x,y,label="",subplot=sp)
fit = EasyFit.fitlinear(x,y)
plot!(fit.x,fit.y,label="",linewidth=2,subplot=sp)

sp=2
y2 = y.^2
plot!(title="Quadratic fit",subplot=sp)
scatter!(x,y2,label="",subplot=sp)
fit = EasyFit.fitquad(x,y2)
plot!(fit.x,fit.y,label="",linewidth=2,subplot=sp)

sp=3
y2 = y.^3
plot!(title="Cubic fit",subplot=sp)
scatter!(x,y2,label="",subplot=sp)
fit = EasyFit.fitcubic(x,y2)
plot!(fit.x,fit.y,label="",linewidth=2,subplot=sp)

sp=4
y2 = 0.3*exp.(y)
plot!(title="Mono-exponential fit",subplot=sp)
scatter!(x,y2,label="",subplot=sp)
fit = EasyFit.fitexp(x,y2)
plot!(fit.x,fit.y,label="",linewidth=2,subplot=sp)

sp=5
y2 = 0.3*exp.(y) + 0.7*exp.(y2)
plot!(title="Bi-exponential fit",subplot=sp)
scatter!(x,y2,label="",subplot=sp)
fit = EasyFit.fitexp(x,y2,n=2)
plot!(fit.x,fit.y,label="",linewidth=2,subplot=sp)

sp=6
y2 = y.^5
plot!(title="Spline fit",subplot=sp)
scatter!(x,y2,label="",subplot=sp)
fit = EasyFit.fitexp(x,y2,n=2)
plot!(fit.x,fit.y,label="",linewidth=2,subplot=sp)

plot!(size=(600,600))



