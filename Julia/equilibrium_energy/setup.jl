using Pkg
Pkg.develop(path = raw"C:\GIT_REPOS\BELL\Julia")
Pkg.instantiate()
using Bell
println("OK ", pkgversion(Bell))
