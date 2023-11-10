# Density of solid, density of fluid and weigth fraction of the non-volatile components as a function of pressure

using DelimitedFiles
using Plots


P    = LinRange(0.87867,1.6429,500) #P range in scaled units

#load the parameters of fitting
data = readdlm("fitting_parameters.txt",',',header=false)
data = data[2,:]

#Functions
Rhos(P,params) = -tanh.(params[1].*(P.-params[2])).*(params[3]/2+params[4]) .+ (params[3]/2-params[4]) .+ params[5] .+ ((P.-params[6])./params[7].*params[8])
Rhof(P,params) = params[1]*log.(P.+params[2]).^params[3]
Xs(P,params)   = -tanh.(params[1]*(P.-params[2]))*params[3]/2 .+ params[3]/2 .+ params[4] 




params = data[1:8]
rhos = Rhos(P,params)

params = data[9:11]
rhof = Rhof(P,params)

params = data[12:15]
xs     = Xs(P,params)

plot(P,rhos)