#First add depencies for the example
using Weave
# doctype="pandoc" (Markdown), "md2html"(HTML),  "md2pdf"(pdf)
# out_path: Path where the output is generated. Can be
# : :doc: Path of the source document,
# :pwd: Julia working directory,
# "somepath": output directory as a String e.g "/home/mpastell/weaveout"
# filename as string e.g. ~/outpath/outfile.tex.
using PolyOrigin

cd(@__DIR__)

dataid = "tetraploid_realpotato"
weave(string(dataid,".jmd"),doctype="md2html",out_path=:pwd)

weave(string(dataid,".jmd"),doctype="pandoc",out_path=:pwd)


# convert_doc(string(dataid,".jmd"),string(dataid,".ipynb"))
# notebook(string(dataid,".jmd"),out_path=:pwd)
