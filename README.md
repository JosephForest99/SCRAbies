
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SCRAbies

Este es un paquete sencillo con funciones para generar una simulación de
crecimiento a nivel de rodal.

## Instalación

Puedes instalar el paquete con la siguiente instrucción:
`devtools::installl_github("JosephForest99/SCRAbies")`

``` r
# devtools::installl_github("JosephForest99/SCRAbies")
```

## Example

Este es un ejemplo para simular un sistema de crecimiento y rendimiento
con la FRA de Bailey & Ware (1983).

``` r
library(SCRAbies)
scr_fra_bailey_ware_abies(E = 1:5, Narb = 2000, IS = 28)
#> # A tibble: 5 × 8
#>   Thin  Densidad     E    IS     AD     N      AB   DCM
#>   <chr>    <dbl> <int> <dbl>  <dbl> <dbl>   <dbl> <dbl>
#> 1 SA        2000     1    28 0.0642 2000  0.00278 0.133
#> 2 SA        2000     2    28 0.189  1999. 0.0168  0.327
#> 3 SA        2000     3    28 0.355  1997. 0.0480  0.553
#> 4 SA        2000     4    28 0.554  1995. 0.101   0.802
#> 5 SA        2000     5    28 0.782  1991. 0.179   1.07
scr_isup_abies(E = 1:5, Narb = 2000, IS = 28)
#> # A tibble: 5 × 9
#>   Thin      E Densidad    IS     AD     N  ISup      AB   DCM
#>   <chr> <int>    <dbl> <dbl>  <dbl> <dbl> <dbl>   <dbl> <dbl>
#> 1 SA        1     2000    28 0.0642 2000      0 0.00235 0.122
#> 2 SA        2     2000    28 0.189  1999.     0 0.0147  0.306
#> 3 SA        3     2000    28 0.355  1997.     0 0.0431  0.524
#> 4 SA        4     2000    28 0.554  1995.     0 0.0919  0.766
#> 5 SA        5     2000    28 0.782  1991.     0 0.165   1.03
```
