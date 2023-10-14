# this script is in the top directory of a package
cd(@__DIR__)

using Pkg
Pkg.activate("./")
Pkg.instantiate()

using Revise
push!(LOAD_PATH,string(@__DIR__,"/src/"))
using EasyFit
