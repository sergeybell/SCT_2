# one-time environment setup
using Pkg; Pkg.activate(@__DIR__); Pkg.develop(path = raw"C:\GIT_REPOS\BELL"); Pkg.add("BenchmarkTools"); Pkg.instantiate()
