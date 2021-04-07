<div align="left">
  <h1 id="dopingfvm"><img src="/images/Logo_notex_dopingFVM.png" width="125" title="dopingFVM logo"> dopingFVM</h1>
  <p>
    <a style="text-decoration: none" href="https://TTp95.github.io/dopingFVM.jl/stable">
      <img alt="Stable" src="https://img.shields.io/badge/docs-stable-blue.svg" />
    </a>
    <a style="text-decoration: none" href="https://TTp95.github.io/dopingFVM.jl/dev">
      <img alt="Dev" src="https://img.shields.io/badge/docs-dev-blue.svg" />
    </a>
    <a style="text-decoration: none" href="https://github.com/TTp95/dopingFVM.jl/actions">
      <img alt="Build Status" src="https://github.com/TTp95/dopingFVM.jl/workflows/CI/badge.svg"/>
    </a>
    <a style="text-decoration: none" href="https://codecov.io/gh/TTp95/dopingFVM.jl">
      <img alt="Coverage" src="https://codecov.io/gh/TTp95/dopingFVM.jl/branch/master/graph/badge.svg" />
    </a>
  </p>
</div>

`dopingFVM` is a Julia package that provides Finite Volume Method (FVM) tools for numerically solving the Partial Differential Equations Systems (PDEs) that describe the transport phenomenon. The `dopingFVM` package uses structured grids with the FVM to numerically solve fluid mechanics, heat transfer, and any similar phenomenon described as a transport equation system.

---

## Table of Content
* [dopingFVM](#dopingfvm)
* [Features](#features)
* [Installation](#installation)
* [Purpose](#propurse)
* [To Do List](#to-do-list)
* [Citation](#to-do-list)

---

# Features

# Installation

Julia dopingFVM *will be* register in the Julia General Repository

```julia
#using the internal package manager
]
add dopingFVM
```

or

```julia
#calling the Pkg manager from the Pkg library
using Pkg

Pkg.add("dopingFVM")
```

To install the nightly version of dopingFVM:

```julia
#only recommended if the user required an urgent hotfix
]
add https://github.com/TTp95/dopingFVM.jl#master
```

*Note: The `master` branch will never be merged with an incomplete in-development feature. Bugs could appear in the code, but no compiling error should be merged in the `master` branch. Commits without a version update will be only for minor hotfixes.*

# Purpose



# To Do List
- [ ] Crear plots agregando bordes
- [ ] Terminar salida archivos transientes PVD
- [ ] cambiar systemconfig a controlconfig
- [ ] evaluar posibilidad de cambiar time.time1,2,3 a un vector time.dt DeltaTime
- [ ] falta sumario, función error relativo imprimir, imprimir enunciado e imprimir cambio de paso de tiempo con el error relativo escrito
- [ ] en caso que la implementación sea steady, que en pantalla diga steady y no el tiempo...


# Citation
