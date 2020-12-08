using PolyOrigin
cd(@__DIR__)
pwd()

# run polyorigin
genofile = "geno.csv"
pedfile = "ped.csv"
outstem = "example_output"
@time polyancestry = polyOrigin(genofile, pedfile;
    refinemap=true,
    refineorder=true,
    isplot=true,
    outstem,
)

# calculate the accuracy of parental phasing and ancestral inference
truefile = "true.csv"
truegeno = readTruegeno!(truefile, polyancestry)
acc = calAccuracy!(truegeno, polyancestry)
println(acc)

# plot conditional probabilities
plotCondprob(polyancestry, truegeno = truegeno, offspring = 1)
animCondprob(polyancestry;
    truegeno = truegeno,
    fps = 0.5,
    outfile = string(outstem,"_condprob.gif"),
)
