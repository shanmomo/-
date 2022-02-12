using Pkg

ToAddList = ["GSL",
    "HypergeometricFunctions",
    ]

for toadditem in ToAddList
    Pkg.add(toadditem)
end
