import Pkg

root = @__DIR__
for project in (joinpath(root, "one_d"), joinpath(root, "axisymmetric"))
    Pkg.activate(project)
    Pkg.resolve()
    Pkg.instantiate()
end
