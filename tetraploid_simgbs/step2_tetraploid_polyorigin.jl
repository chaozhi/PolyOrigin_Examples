using PolyOrigin
cd(@__DIR__)
pwd()

dataid = "tetraploid_simgbs"
genofile = string(dataid,"_geno.csv")
pedfile = string(dataid,"_ped.csv")
outstem = string(dataid,"_output")

polygeno = readPolyGeno(genofile, pedfile)
PolyOrigin.plotdesign(polygeno)

@time polyancestry = polyOrigin(genofile,pedfile;outstem)
outfiles = filter(x -> occursin(outstem, x), readdir())
println(outfiles)

# plot relative frequencies of valent configurations
polyancestry = readPolyAncestry(outstem*"_polyancestry.csv")
valentfreq = calvalentfreq(polyancestry)
plotvalentfreq(valentfreq)

# calculate estimation accuracies
truegeno = readTruegeno!(string(dataid,"_true.csv"), polyancestry)
println(keys(truegeno))
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
