```@meta
CurrentModule = WaveOperators   
```

# WaveOperators 


## Overview
This package contains utility functions to setup physical design problems[^1],
using an integral equation approximation to the Helmhotlz equation.
While it was developed for computational bounds in photonics[^2], we believe
that these operators may be useful for anyone working on problems in physical
design. To get started, check out the examples!

For more details on the methodology and performance tricks, please see appendices
A & B of [Bounds on Efficiency Metrics in Photonics](https://arxiv.org/abs/2204.05243).


#### Examples:
```@contents
Pages = ["examples/create-problem-data.md", "examples/simple-example.md"]
Depth = 1
```



## References
[^1]: Angeris, G., Vučković, J., & Boyd, S. (2021). Heuristic methods and performance bounds for photonic design. Optics Express, 29(2), 2827-2854.
[^2]: Angeris, G., Diamandis, T., Vučković, J., & Boyd, S. (2022). Bounds on Efficiency Metrics in Photonics. arXiv preprint arXiv:2204.05243.
