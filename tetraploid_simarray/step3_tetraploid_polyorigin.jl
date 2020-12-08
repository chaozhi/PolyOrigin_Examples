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

@time polyancestry = polyOrigin!(polygeno;
    refinemap = true,
    refineorder = true,
    outstem = outstem,
)
outfiles = filter(x -> occursin(outstem, x), readdir())
println(outfiles)

polyancestry = readPolyAncestry(string(outstem, "_genoprob.csv"), pedfile)
truefile = string(dataid,"_true.csv")
truegeno = readTruegeno!(truefile, polyancestry)
acc = calAccuracy!(truegeno, polyancestry)
println(acc)

plotCondprob(polyancestry, truegeno = truegeno, offspring = 1)
animCondprob(
    polyancestry,
    truegeno = truegeno,
    fps = 0.5,
    outfile=string(outstem,"_condprob.gif"),
)

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
