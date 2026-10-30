using Pkg
Pkg.develop(path = raw"C:\GIT_REPOS\BELL")
Pkg.instantiate()
using Bell
println("OK ", pkgversion(Bell))
