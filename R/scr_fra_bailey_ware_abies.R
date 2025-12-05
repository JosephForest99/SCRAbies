#' @title Sistema de crecimiento con la FRA de Bailey y Ware (1983)
#' @description
#' Simulador de crecimiento de atributos de rodal como altura dominante, mortalidad y área basal
#' sin o con efecto de aclareo.
#'
#' @param E, Vector de valores numérico. Al menos utilizar un vector con dos valores.
#' @param Narb, Valor numérico. Indica la densidad inicial del rodal - número de árboles por ha.
#' @param IS, Valor numérico. Indica el índice de sitio del rodal.
#' @param Thinning, Valor lógico. FALSE no simula el crecimiento con efecto de aclareo. TRUE incluye el efecto de aclareo
#' @param Ets, Valor o vector de valores numérico. Indica la edad del rodal cuando se aplica el aclareo.
#' @param ABpts, Valor o vector de valores numérico. Indica la intensidad de aclareo respecto al AB, por ejemplo 0.3, equivale al 30 %.
#'
#' @return Regresa un objeto tibble.
#' @examples
#' # scr_fra_bailey_ware_abies(E = 1:5, Narb = 2000, IS = 28)
#'
#'
#' @import dplyr
#' @import tidyr

#' @export

# Función para simular un SCR con aclareo
scr_fra_bailey_ware_abies <- function(E, Narb, IS, Thinning = FALSE, Ets = NULL, ABpts = NULL){

  `%nin%` = Negate(`%in%`)

  if (isTRUE(Thinning)) {
    if (length(Ets) != length(ABpts)) {
      stop("Los vectores `Ets` y `ABpts` debe tener la misma longitud")
    } else{
      E = E |> c(Ets) |> sort()
      n = length(E)

      if (n <= 2) {
        stop("Al tratarse de un aclareo es necesario que E2>E1")
      }else{
        idthin = as.vector(rep(NA, n));
        Thin = as.vector(rep(NA, n));
        Et = as.vector(rep(NA, n));
        ABpt = as.vector(rep(NA, n))
        for (j in seq_along(Ets)) {
          for (i in length(E):1) {
            if (E[i] == Ets[j]) {
              idthin[i] = 1
              Thin[i] = paste0("A",j)
              Et[i] = Ets[j]
              ABpt[i] = ABpts[j]
              break
            }
          }
        }

        df <- data.frame(E = E, idthin = idthin, Thin = Thin, Et = Et, ABpt = ABpt) |>
          fill(c(idthin, Thin, Et), .direction = "down") |>
          mutate(idthin = ifelse(is.na(idthin), 0, idthin)) |>
          mutate(Thin = ifelse(is.na(Thin), "SA", Thin)) |>
          mutate_at(.vars = vars(c("Et")), .funs = \(x) ifelse(is.na(x), 0, x))
      }

    }
  } else{
    E = E
    n = length(E)

    df <- data.frame(E = E, idthin = 0, Thin = "SA", Et = 0, ABpt = 0)

  }

  idthin = df$idthin
  Thin = df$Thin
  Et = df$Et


  if (n <= 1) {
    stop("Se requiere que la E2 sea al menos un año mayor/menor que la E1 para realizar la proyección")
  }

  # **ALTURA DOMINANTE**
  Eb = 80
  AD2_fn <- function(AD1, E1, E2){
    B1 = 0.0191326; B2 = 1.5619704
    B0 = E1^B2/AD1 - B1*E1^B2
    AD = E2^B2 / (B0 + B1*E2^B2)
    return(AD)
  }

  AD2 = AD2_fn(AD1 = IS, E1 = Eb, E2 = E)

  # **MORTALIDAD**
  N_fn <- function(N1, E1, E2){
    B1 = 0.02025389
    N2 = (N1^-0.5 + B1 * ((E2/100)^2 - (E1/100)^2) )^(-2)
    return(N2)
  }

  Nt = vector(mode = "numeric", length = length(E));
  Nb = vector(mode = "numeric", length = length(E));
  ABpt = as.vector(df$ABpt)
  Nt_bandera = 0; Nb_bandera = 0; ABpt_bandera = 0;

  N2 = NA
  N2[1] = Narb
  for (i in 2:length(E)) {
    if (E[i] == E[i-1]) {
      Nt_bandera = N2[i-1] * ABpt[i]**(1/1.649317)
      Nb_bandera = N2[i-1]
      ABpt_bandera = ABpt[i]
      Nt[i] = N2[i-1] * ABpt[i]**(1/1.649317)
      Nb[i] = N2[i-1]
      N2[i] = N2[i-1] - N2[i-1] * ABpt[i]**(1/1.649317)
    }else{
      N2[i] = N_fn(N1 = N2[i-1], E1 = E[i-1], E2 = E[i])
      Nt[i] = Nt_bandera
      Nb[i] = Nb_bandera
      ABpt[i] = ABpt_bandera
    }
  }

  # **ÁREA BASAL**
  AB1_fn <- function(E1, N1, AD1, Xt, Et, idthin){
    B0 = -6.3708686; B1 = 1.1695811; B2 = 0.7713131; B3 = 0.4866774; B4 = -49.9570687;
    AB1 = exp(B0)*AD1**B1*E1**B2*N1**B3*
      exp(ifelse(idthin == 0 , 0, B4*Xt/(E1**2*Et)))
    return(AB1)
  }
  AB2_fn <- function(AB1, E1, E2, N1, N2, AD1, AD2, Xt, Et, idthin){
    B0 = -6.3708686; B1 = 1.1695811; B2 = 0.7713131; B3 = 0.4866774; B4 = -49.9570687;
    AB2 = AB1*(AD2/AD1)**B1*(E2/E1)**B2*(N2/N1)**B3*
      exp(ifelse(idthin == 0, 0, B4*Xt*(1/(E2**2*Et)-1/(E1**2*Et)) ))
    return(AB2)
  }

  AB2_fn90 <- function(AB1, E1, E2, AD1, AD2, Xt, Et, idthin){
    B1 = 9.066359e-01; B4 = -4.106880e+04;
    AB2 = AB1*(AD2/AD1)^B1
    exp(ifelse(idthin == 0, 0, B4*Xt*(1/(E2**2*Et)-1/(E1**2*Et)) ))
    return(AB2)
  }

  ABt = vector(mode = "numeric", length = length(E));
  ABb = vector(mode = "numeric", length = length(E));
  ABt_bandera = 0;
  ABb_bandera = 0
  DCM = vector(mode = "numeric", length = length(E));
  DCMt = vector(mode = "numeric", length = length(E));
  DCMb = vector(mode = "numeric", length = length(E));
  Xt = vector(mode = "numeric", length = length(E));
  DCMt_bandera = 0; DCMb_bandera = 0

  AB2 = NA
  AB2[1] = AB1_fn(E1 = E[1], N1 = N2[1], AD1 = AD2[1], Xt = Xt[1],
                  Et = Et[1], idthin = idthin[1])
  DCM[1] = sqrt(40000/pi * AB2[1]/N2[1])

  for (i in 2:length(E)) {
    if (Et[i] == 0) {
      if (E[i]<90) {
        AB2[i] = AB2_fn(AB1 = AB2[i-1], E1 = E[i-1], E2 = E[i],
                        N1 = N2[i-1], N2 = N2[i], AD1 = AD2[i-1], AD2 = AD2[i],
                        Xt = Xt[i], Et = Et[i], idthin = idthin[i])
      } else{
        AB2[i] = AB2_fn90(AB1 = AB2[i-1],
                          E1 = E[i-1], E2 = E[i],
                          AD1 = AD2[i-1], AD2 = AD2[i],
                          Xt = Xt[i], Et = Et[i], idthin = idthin[i])
      }
      DCM[i] = sqrt(40000/pi * AB2[i]/N2[i])
    }else if (E[i] == E[i-1]){
      ABt_bandera = AB2[i-1] * ABpt[i]
      ABt[i] = AB2[i-1] * ABpt[i]
      ABb_bandera = AB2[i-1]
      ABb[i] = AB2[i-1]
      AB2[i] = AB2[i-1] - ABt[i]
      DCM[i] = sqrt(40000/pi * AB2[i]/N2[i])
      DCMb_bandera = DCM[i-1]
      DCMb[i] = DCM[i-1]
      DCMt_bandera = sqrt(40000/pi * AB2[i]/N2[i])
      DCMt[i] = sqrt(40000/pi * AB2[i]/N2[i])
      Xt_bandera = 1-(DCMt[i]/DCMb[i])
      Xt[i] = 1-(DCMt[i]/DCMb[i])
    }else{
      ABt[i] = ABt_bandera
      ABb[i] = ABb_bandera
      DCMt[i] = DCMt_bandera
      DCMb[i] = DCMb_bandera
      Xt[i] = Xt_bandera
      if (E[i]<90) {
        AB2[i] = AB2_fn(AB1 = AB2[i-1], E1 = E[i-1], E2 = E[i],
                        N1 = N2[i-1], N2 = N2[i], AD1 = AD2[i-1], AD2 = AD2[i],
                        Xt = Xt[i], Et = Et[i], idthin = idthin[i])
      } else{
        AB2[i] = AB2_fn90(AB1 = AB2[i-1],
                          E1 = E[i-1], E2 = E[i],
                          AD1 = AD2[i-1], AD2 = AD2[i],
                          Xt = Xt[i], Et = Et[i], idthin = idthin[i])
      }
      DCM[i] = sqrt(40000/pi * AB2[i]/N2[i])
    }
  }
  df <- tibble(Thin = Thin, Densidad = Narb, E = E, IS = IS, AD = AD2,
               N = N2, AB = AB2, DCM = DCM
  )
  return(df)
}
