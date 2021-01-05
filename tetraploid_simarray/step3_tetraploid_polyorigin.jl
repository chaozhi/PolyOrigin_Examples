using PolyOrigin
using Plots
cd(@__DIR__)
pwd()

dataid = "tetraploid_simarray"
genofile = string(dataid,"_geno_disturbed.csv")
pedfile = string(dataid,"_ped.csv")
outstem = string(dataid,"_output")

polygeno = readPolyGeno(genofile, pedfile)
PolyOrigin.plotdesign(polygeno)

@time polyancestry =  polyOrigin(genofile,pedfile;
    refinemap = true,
    refineorder = true,
    outstem = outstem,
)
outfiles = filter(x -> occursin(outstem, x), readdir())
println(outfiles)

# plot relative frequencies of valent configurations
polyancestry = readPolyAncestry(outstem*"_polyancestry.csv")
valentfreq = calvalentfreq(polyancestry)
plotvalentfreq(valentfreq)

# calculate estimation accuracies
truefile = string(dataid,"_true.csv")
truegeno = readTruegeno!(truefile, polyancestry)
acc = calAccuracy!(truegeno, polyancestry)
println(acc)

# plot conditional probabilities
plotCondprob(polyancestry, truegeno = truegeno, offspring = 1)
animCondprob(
    polyancestry,
    truegeno = truegeno,
    fps = 0.5,
    outfile=string(outstem,"_condprob.gif"),
)

# comparing maps
polygeno = readPolyGeno(genofile, pedfile)
fig = plotMapComp(
    truegeno.truemap,
    polygeno.markermap,
    xlabel = "True poistion (cM)",
    ylabel = "Input position (cM)",
)
fig2 = plotMapComp(
    truegeno.truemap,
    polyancestry.markermap,
    xlabel = "True poistion (cM)",
    ylabel = "Estimated position (cM)",
)
plot(fig, fig2)
