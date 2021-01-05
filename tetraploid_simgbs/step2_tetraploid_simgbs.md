---
author: "Chaozhi Zheng"
title: "Haplotype reconstruction in a simulated tetrapoloid multiparental population"
---


# Setup PolyOrigin

Non-parallel computation for simulated data with one linkage group
~~~~{.julia}
using PolyOrigin
~~~~~~~~~~~~~




# Input files
Set working directory to the direcotry of this file.
~~~~{.julia}
cd(@__DIR__)
~~~~~~~~~~~~~




Set input files and outstem
~~~~{.julia}
dataid = "tetraploid_simgbs"
genofile = string(dataid,"_geno.csv")
pedfile = string(dataid,"_ped.csv")
outstem = string(dataid,"_output")
~~~~~~~~~~~~~




Check input files. In the pedfile, the parents of founders are set to 0.
~~~~{.julia}
using CSV, DataFrames
show(CSV.read(pedfile,DataFrame)[1:8,:],eltypes=false)
~~~~~~~~~~~~~

~~~~
8×5 DataFrame
 Row │ Individual    Population    MotherID    FatherID    Ploidy
─────┼────────────────────────────────────────────────────────────
   1 │ A                       0   0           0                4
   2 │ B                       0   0           0                4
   3 │ C                       0   0           0                4
   4 │ D                       0   0           0                4
   5 │ E                       0   0           0                4
   6 │ AxB0001                 1   A           B                4
   7 │ AxB0002                 1   A           B                4
   8 │ AxB0003                 1   A           B                4
~~~~




In the genofile, columns 1-3 denote the marker map, and rest
columns for the read counts of parents and offspring.
~~~~{.julia}
show(CSV.read(genofile,DataFrame)[1:5,1:6],eltypes=false)
~~~~~~~~~~~~~

~~~~
5×6 DataFrame
 Row │ marker    chromosome    position    A       B       C
─────┼───────────────────────────────────────────────────────────
   1 │ A_0001    A                  0.0    16|6    18|8    8|9
   2 │ A_0002    A                  1.94   0|22    6|14    0|21
   3 │ A_0003    A                  2.04   8|22    8|10    5|21
   4 │ A_0004    A                  2.2    19|4    12|8    1|15
   5 │ A_0005    A                  4.94   1|17    2|16    17|0
~~~~





plot crossdesign for the 5x5 diallel cross.
~~~~{.julia}
polygeno = readPolyGeno(genofile, pedfile)
PolyOrigin.plotdesign(polygeno)
~~~~~~~~~~~~~

![](figures/step2_tetraploid_simgbs_6_1.png)\ 




# Run polyOrigin

run polyOrigin
~~~~{.julia}
@time polyancestry= polyOrigin(genofile, pedfile; outstem)
~~~~~~~~~~~~~


Keyarg `outstem` specifies the stem of output files.

The returned polyancestry from polyOrigin has been saved.

~~~~{.julia}
outfiles = filter(x->occursin(outstem,x), readdir())
~~~~~~~~~~~~~

~~~~
5-element Array{String,1}:
 "tetraploid_simgbs_output.log"
 "tetraploid_simgbs_output_genoprob.csv"
 "tetraploid_simgbs_output_parentphased.csv"
 "tetraploid_simgbs_output_polyancestry.csv"
 "tetraploid_simgbs_output_postdoseprob.csv"
~~~~




[Click to view log file](tetraploid_simgbs_output.log)

# Check output

The main output file constains multiple dataframes.
~~~~{.julia}
ancestryfile = outstem*"_polyancestry.csv"
res = PolyOrigin.readdlm2dict(ancestryfile)
keys(res)
~~~~~~~~~~~~~

~~~~
Base.KeySet for a Dict{SubString{String},DataFrame} with 10 entries. Keys:
  "parentinfo"
  "offspringinfo"
  "valentprob"
  "ancestralgenotype"
  "correction"
  "designinfo"
  "genoprob"
  "valentlist"
  "delmarker"
  "parentgeno"
~~~~




The parentgeno dataframe. At each marker for each parent, the phased parental
genotype is given by e.g. 1|2|2|1, where 1 and 2 denote the two alleles.
~~~~{.julia}
show(res["parentgeno"][1:10:100,1:6],eltypes=false)
~~~~~~~~~~~~~

~~~~
10×6 DataFrame
 Row │ marker  chromosome  position  A        B        C
─────┼─────────────────────────────────────────────────────────
   1 │ A_0001  A               0.0   1|1|2|1  2|2|1|1  2|2|1|1
   2 │ A_0011  A               8.08  1|1|2|1  1|1|2|1  1|1|1|2
   3 │ A_0021  A              12.56  1|2|2|1  1|1|2|1  2|2|1|2
   4 │ A_0031  A              19.64  2|2|1|2  2|1|2|1  1|1|1|2
   5 │ A_0041  A              25.65  1|1|2|2  2|1|2|2  1|2|1|2
   6 │ A_0051  A              36.6   2|1|2|2  1|2|2|1  2|2|2|2
   7 │ A_0061  A              43.9   1|1|2|2  2|2|1|1  1|1|2|2
   8 │ A_0071  A              57.17  2|1|2|1  2|1|1|2  2|1|2|2
   9 │ A_0081  A              66.28  2|1|1|2  1|2|2|2  1|1|2|2
  10 │ A_0091  A              75.46  2|1|1|1  1|1|1|1  2|1|1|2
~~~~





The genoprob dataframe. At each marker for each offspring, the posterior probability
distribution is given by e.g. 57|87=>0.336|0.664, meaning that the posterior
probabilities of ancestral genotypes 57 and 87 are 0.336 and 0.664, respectively.
~~~~{.julia}
show(res["genoprob"][1:10:100,1:4],eltypes=false)
~~~~~~~~~~~~~

~~~~
10×4 DataFrame
 Row │ marker  chromosome  position  AxB0001
─────┼─────────────────────────────────────────────────────────────────
   1 │ A_0001  A               0.0   23|24|26|27|53|54|56|57|77|83|84…
   2 │ A_0011  A               8.08  17|27|54|57=>0.015|0.011|0.175|0…
   3 │ A_0021  A              12.56  56|57=>0.006|0.994
   4 │ A_0031  A              19.64  17|54|56|57|67|87=>0.003|0.027|0…
   5 │ A_0041  A              25.65  24|53|54|57|59|60|64|87=>0.001|0…
   6 │ A_0051  A              36.6   54|57=>0.023|0.975
   7 │ A_0061  A              43.9   27|52|54|56|57|67|77=>0.003|0.00…
   8 │ A_0071  A              57.17  16|17|27|36|37|47|52|54|55|56|57…
   9 │ A_0081  A              66.28  47|57|67|87=>0.006|0.931|0.056|0…
  10 │ A_0091  A              75.46  47|56|57=>0.001|0.001|0.997
~~~~




The ancestral genotypes are defined for each sub-population in the following.
The columns state and stateindex denote the ancestral genotypes.
~~~~{.julia}
show(res["ancestralgenotype"][sort(rand(1:100,10)),:],eltypes=false)
~~~~~~~~~~~~~

~~~~
10×5 DataFrame
 Row │ population  parentindex  parent  stateindex  state
─────┼──────────────────────────────────────────────────────
   1 │          1  1|2          A|B             10  1-1-8-8
   2 │          1  1|2          A|B             33  1-4-5-7
   3 │          1  1|2          A|B             34  1-4-5-8
   4 │          1  1|2          A|B             37  1-4-6-8
   5 │          1  1|2          A|B             38  1-4-7-7
   6 │          1  1|2          A|B             51  2-3-5-5
   7 │          1  1|2          A|B             55  2-3-6-6
   8 │          1  1|2          A|B             67  2-4-6-8
   9 │          1  1|2          A|B             89  3-4-7-8
  10 │          1  1|2          A|B             93  4-4-5-7
~~~~




where 1-4 denote the homologs in the first parent, and 5-8 for the second parent.

The valentprob dataframe. For each offpsring in each chromosome, the valent is
given eg. 1:3-2:4&5:8-6:7|1:3-2:4&5:6:7:8 where & delimits the valent configurations
between two parents, and |  delimits each possible combination of valent configurations.
The column valentprob gives the full posterior probabilties.
~~~~{.julia}
println(join(names(res["valentprob"])[[1,2,5,7]],", "))
println(join(Vector(res["valentprob"][1,[1,2,5,7]]),", "))
~~~~~~~~~~~~~

~~~~
chromosome, individual, valent, valentprob
A, AxB0001, 1:3-2:4&5:6-7:8|1:2-3:4&5:6-7:8|1:2:3:4&5:6-7:8|1:2-3:4&5:6:7:8
|1:3-2:4&5:6:7:8|1:2:3:4&5:6:7:8, 0.4043|0.3901|0.1456|0.0263|0.0243|0.0094
~~~~





# Relative frequencies of valent configurations

We can read polyancestry
~~~~{.julia}
polyancestry = readPolyAncestry(outstem*"_polyancestry.csv")
~~~~~~~~~~~~~




Relative frequencies of bi- or quadri-valent formations for each parent in
each chromosome
~~~~{.julia}
valentfreq = calvalentfreq(polyancestry)
~~~~~~~~~~~~~

~~~~
5×6 DataFrame
 Row │ chromosome  parent  1:2-3:4  1:3-2:4  1:4-2:3  1:2:3:4
     │ String      String  Float64  Float64  Float64  Float64
─────┼────────────────────────────────────────────────────────
   1 │ A           A        0.3875   0.275    0.3125   0.025
   2 │ A           B        0.3625   0.3625   0.25     0.025
   3 │ A           C        0.3      0.4      0.3      0.0
   4 │ A           D        0.2875   0.325    0.3625   0.025
   5 │ A           E        0.2875   0.325    0.375    0.0125
~~~~





plot the relative frequences averaging overage chromosomes/parents
~~~~{.julia}
plotvalentfreq(valentfreq)
~~~~~~~~~~~~~

![](figures/step2_tetraploid_simgbs_17_1.png)\ 



where 1:2:3:4 denotes quadrivalent formation and the others for bivalent formations

# Calculate estimation error probability

Read true value file used in simulating data
~~~~{.julia}
truegeno=readTruegeno!(string(dataid,"_true.csv"),polyancestry)
keys(truegeno)
~~~~~~~~~~~~~

~~~~
(:truemap, :estmap, :parentgeno, :offspringgeno)
~~~~




where the absolute phase of `polyancestry.parentgeno` is set to be consistent
with that of `truegeno.parentgeno`, and `polyancestry.genoprob` and
`polyancestry.haploprob` are re-ordered accordingly, and
* `:truemap` = true genetic map,
* `:estmap` = estimated genetic map resulting from map refinement,
* `:parentgeno` = true phased paental genotypes,
* `:offspringgeno` = true parental origins.

~~~~{.julia}
acc=calAccuracy!(truegeno,polyancestry)
~~~~~~~~~~~~~

~~~~
(ndoseerr = 1, nphaseerr = 1, nalleleerr = 1, nparentgeno = 600, assignerr 
= 0.102187, callerr = 0.0685833, delfraction = 0.0)
~~~~




where
* `ndoseerr` = number of wrongly estimated parental dosages,
* `nphaseerr` = number of wrongly estimated parental phases,
* `nalleleerr` = number of wrongly estimated alleles for  parental phased genotypes,
* `nparentgeno` = total number of parental genotypes,
* `assignerr` = one minus probability of true unphased ancestral genotype,
* `callerr` = fraction of unphased origin-genotypes being wrongly called,
* `delfraction` = fraction of markers being deleted.

# Visualize conditional probability

Visualize haplotype probabilities of single offspring
~~~~{.julia}
plotCondprob(polyancestry,truegeno=truegeno,offspring=1)
~~~~~~~~~~~~~

![](figures/step2_tetraploid_simgbs_20_1.png)\ 




Visualize haplotype probabilities of all offspring
~~~~{.julia}
animCondprob(polyancestry,truegeno=truegeno, fps=0.5,
  outfile=outstem*"_condprob.gif")
~~~~~~~~~~~~~



![](tetraploid_simgbs_output_condprob.gif)

where `fps` specifies the number of frames per seconds, for exmaple, `fps=0.5`
means one figure every two seconds.
