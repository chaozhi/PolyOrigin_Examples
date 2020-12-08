using Distributed
addprocs(4) # add worker processes on local machine
@info string("nworkers=", nworkers())
@everywhere  using PolyOrigin
cd(@__DIR__)
pwd()

genofile = "TableS2_dose.csv"
pedfile = "TableS3_ped.csv"
outstem = "potato_output"
chrsubset = 1:4
@time polyancestry = polyOrigin(genofile, pedfile;
    isphysmap = true,
    recomrate = 1.25,
    snpsubset = 1:5:1000,
    chrsubset = chrsubset,
    isparallel = true,
    refinemap = true,
    refineorder = false,
    outstem = outstem,
)
outfiles = filter(x -> occursin(outstem, x), readdir())
println("outfiles=", outfiles)

polygeno = readPolyGeno(genofile, pedfile;isphysmap = true)
fig = plotMapComp(
    polygeno.markermap[chrsubset],
    polyancestry.markermap,
    xlabel = "Physical poistion (Mbp)",
    ylabel = "Estimated position (cM)",
)
