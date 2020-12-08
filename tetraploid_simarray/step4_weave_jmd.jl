#First add depencies for the example
# using Pkg; Pkg.add.(["Plots", "DSP"])
# using Plots
using Weave
# doctype="pandoc" (Markdown), "md2html"(HTML),  "md2pdf"(pdf)
# out_path: Path where the output is generated. Can be
# : :doc: Path of the source document,
# :pwd: Julia working directory,
# "somepath": output directory as a String e.g "/home/mpastell/weaveout"
# filename as string e.g. ~/outpath/outfile.tex.
using PolyOrigin
using Plots
cd(@__DIR__)

weave("step3_tetraploid_simarray.jmd", doctype = "md2html", out_path = :pwd)

weave("step3_tetraploid_simarray.jmd", doctype = "pandoc", out_path = :pwd)
