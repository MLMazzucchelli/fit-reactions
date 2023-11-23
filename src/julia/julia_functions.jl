# Density of solid, density of fluid and weigth fraction of the non-volatile components as a function of pressure

using DelimitedFiles
using GLMakie
#using CairoMakie

P    = LinRange(0.87867,1.6429,500) #P range in scaled units

# Load the fitting parameters 
data = readdlm("fitting_parameters.txt",',',header=false)
data = data[2,:]

# Fitting functions
Rhos(P,params) = -tanh.(params[1].*(P.-params[2])).*(params[3]/2+params[4]) .+ (params[3]/2-params[4]) .+ params[5] .+ ((P.-params[6])./params[7].*params[8])
Rhof(P,params) = params[1]*log.(P.+params[2]).^params[3]
Xs(P,params)   = -tanh.(params[1]*(P.-params[2]))*params[3]/2 .+ params[3]/2 .+ params[4] 



# Calculate the density and mass fraction from the fitting equations
params = data[1:8]
rhos = Rhos(P,params)

params = data[9:11]
rhof = Rhof(P,params)

params = data[12:15]
xs     = Xs(P,params)

# Plots

fig = Figure()
lines(fig[1, 1], P,rhos)
lines(fig[2, 1], P,rhof)
lines(fig[3, 1], P, xs)
fig


#= p1 = plot(P,rhos)
#title!("Solid density [kg/m3]")
xlabel!("P")
ylabel!("Solid density [kg/m3]")
p2 = plot(P,rhof)
xlabel!("P")
ylabel!("Fluid density [kg/m3]")
p3 = plot(P,rhos)
plot(p1,p2,p3, layout=(3,1)) =#