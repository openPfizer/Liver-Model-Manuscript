# Simple setup script to activate the project.
# First time through this step can be quite long, so we separated out from main.jl

cd(dirname(@__FILE__))
using Pkg
Pkg.activate(".")
Pkg.instantiate()