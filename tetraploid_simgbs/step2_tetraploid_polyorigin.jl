using PolyOrigin
cd(@__DIR__)
pwd()

dataid = "tetraploid_simgbs"
genofile = string(dataid,"_geno.csv")
pedfile = string(dataid,"_ped.csv")
outstem = string(dataid,"_output")
@time polyancestry = polyOrigin(genofile,pedfile;outstem)
outfiles = filter(x -> occursin(outstem, x), readdir())
println(outfiles)

polyancestry = readPolyAncestry(string(outstem, "_polyancestry.csv"))
truegeno = readTruegeno!(string(dataid,"_true.csv"), polyancestry)
println(keys(truegeno))
acc = calAccuracy!(truegeno, polyancestry)
println(acc)

plotCondprob(polyancestry, truegeno = truegeno, offspring = 1)
animCondprob(
    polyancestry,
    truegeno = truegeno,
    fps = 0.5,
    outfile=string(outstem,"_condprob.gif"),
)
