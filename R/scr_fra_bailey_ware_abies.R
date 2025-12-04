#' @export

# Función para simular un SCR con aclareo
scr_fra_bailey_ware_abies <- function(E, Narb, IS, Thinning = FALSE, Ets = NULL, ABpts = NULL){

  # Función negativa
  `%nin%` = Negate(`%in%`)

  if (isTRUE(Thinning)) {
    if (length(Ets) != length(ABpts)) {
      stop("Los vectores `Ets` y `ABpts` debe tener la misma longitud")
    } else{
      # Si es TRUE es necesario los valores de Et
      E = E |> c(Ets) |> sort()
      n = length(E)

      if (n <= 2) {
        stop("Al tratarse de un aclareo es necesario que E2>E1")
      }else{
        # Crear el df base
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

        # Base de datos creada
        df <- data.frame(E = E, idthin = idthin, Thin = Thin, Et = Et, ABpt = ABpt) |>
          fill(c(idthin, Thin, Et), .direction = "down") |>
          mutate(idthin = ifelse(is.na(idthin), 0, idthin)) |>
          mutate(Thin = ifelse(is.na(Thin), "SA", Thin)) |>
          mutate_at(.vars = vars(c("Et")), .funs = \(x) ifelse(is.na(x), 0, x))
      }

    }
  } else{
    # Si es FALSE no es necesario Et
    E = E
    n = length(E)

    df <- data.frame(E = E, idthin = 0, Thin = "SA", Et = 0, ABpt = 0)

  }

  # RECUPERACIÓN DE VECTORES
  idthin = df$idthin
  Thin = df$Thin
  Et = df$Et


  # CONDICIONES PARA EJECUTAR FUNCIÓN
  if (n <= 1) {
    stop("Se requiere que la E2 sea al menos un año mayor/menor que la E1 para realizar la proyección")
  }
  # if(is.null(AD1) & is.null(IS)) {
  #   stop("Se tiene que proporcionar información de AD1 o IS para realizar la proyección")
  # }

  # ESCENARIOS DE APLICACIÓN DEL MODELO DE AD E IS

  # **ALTURA DOMINANTE**

  # Edad base
  Eb = 80

  # Función de altura dominante
  AD2_fn <- function(AD1, E1, E2){
    # Parámetros estimados
    B1 = 0.0191326; B2 = 1.5619704
    # Parámetro relacionado con el sitio
    B0 = E1^B2/AD1 - B1*E1^B2
    # Modelo de hossfeld iv B0 de proyección de AD
    AD = E2^B2 / (B0 + B1*E2^B2)
    return(AD)
  }

  AD2 = AD2_fn(AD1 = IS, E1 = Eb, E2 = E)

  # **MORTALIDAD**

  # Función de mortalidad
  N_fn <- function(N1, E1, E2){
    # Parámetro
    B1 = 0.02025389
    N2 = (N1^-0.5 + B1 * ((E2/100)^2 - (E1/100)^2) )^(-2)
    return(N2)
  }

  # Algunas variables para realizar el bucle de mortalidad
  Nt = vector(mode = "numeric", length = length(E));
  Nb = vector(mode = "numeric", length = length(E));
  ABpt = as.vector(df$ABpt)
  Nt_bandera = 0; Nb_bandera = 0; ABpt_bandera = 0;

  # Bucle para estimar mortaliad
  N2 = NA
  N2[1] = Narb
  for (i in 2:length(E)) {
    # Condición para restar el número de árboles en el año de aclareo
    if (E[i] == E[i-1]) {
      # N2[i] = N2[i-1] - Nt[i]

      # Almacenar el valor del aclareo en "banderas"
      Nt_bandera = N2[i-1] * ABpt[i]**(1/1.649317)
      Nb_bandera = N2[i-1]
      ABpt_bandera = ABpt[i]

      # Valores relacionados con el aclareo
      Nt[i] = N2[i-1] * ABpt[i]**(1/1.649317)
      Nb[i] = N2[i-1]

      # Disminución de árboles de acuerdo a la proporción de AB asignado
      N2[i] = N2[i-1] - N2[i-1] * ABpt[i]**(1/1.649317)
    }else{
      N2[i] = N_fn(N1 = N2[i-1], E1 = E[i-1], E2 = E[i])

      # Banderas (para guardar el último valor)
      Nt[i] = Nt_bandera
      Nb[i] = Nb_bandera
      ABpt[i] = ABpt_bandera
    }
  }

  # **ÁREA BASAL**

  # Predicción
  AB1_fn <- function(E1, N1, AD1, Xt, Et, idthin){
    # Parámetros del modelo
    B0 = -6.3708686; B1 = 1.1695811; B2 = 0.7713131; B3 = 0.4866774; B4 = -49.9570687;
    # Modelo de predicción
    AB1 = exp(B0)*AD1**B1*E1**B2*N1**B3*
      exp(ifelse(idthin == 0 , 0, B4*Xt/(E1**2*Et)))

    return(AB1)
  }
  # Proyección
  AB2_fn <- function(AB1, E1, E2, N1, N2, AD1, AD2, Xt, Et, idthin){
    # Parámetros del modelo
    B0 = -6.3708686; B1 = 1.1695811; B2 = 0.7713131; B3 = 0.4866774; B4 = -49.9570687;
    # Modelo de proyección
    AB2 = AB1*(AD2/AD1)**B1*(E2/E1)**B2*(N2/N1)**B3*
      exp(ifelse(idthin == 0, 0, B4*Xt*(1/(E2**2*Et)-1/(E1**2*Et)) ))

    return(AB2)
  }

  # Proyección E>90 AÑOS
  AB2_fn90 <- function(AB1, E1, E2, AD1, AD2, Xt, Et, idthin){
    # Parámetros del modelo
    B1 = 9.066359e-01; B4 = -4.106880e+04;
    # Modelo de proyección
    AB2 = AB1*(AD2/AD1)^B1
    exp(ifelse(idthin == 0, 0, B4*Xt*(1/(E2**2*Et)-1/(E1**2*Et)) ))

    return(AB2)
  }

  # Algunas variables para realizar el bucle del modelo de AB
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

    # Aplicación del modelo SIN ACLAREO
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

      # AÑO (s) de aclareo (s) - Solo se considera la diferencia de ABb - ABt
    }else if (E[i] == E[i-1]){
      ABt_bandera = AB2[i-1] * ABpt[i]
      ABt[i] = AB2[i-1] * ABpt[i]

      ABb_bandera = AB2[i-1]
      ABb[i] = AB2[i-1]

      # ESTIMACIÓN DE AB COMO DIFERENCIA (ABb-ABt)
      AB2[i] = AB2[i-1] - ABt[i]

      # ESTIMACIÓN DE AB DIRECTAMENTE CON EL MODELO
      # AB2[i] = AB2_fn(AB1 = AB2[i-1], E1 = E[i-1], E2 = E[i],
      #                 N1 = N2[i-1], N2 = N2[i], AD1 = AD2[i-1], AD2 = AD2[i],
      #                 Xt = Xt[i], Et = Et[i], idthin = idthin[i])

      DCM[i] = sqrt(40000/pi * AB2[i]/N2[i])

      DCMb_bandera = DCM[i-1]
      DCMb[i] = DCM[i-1]

      DCMt_bandera = sqrt(40000/pi * AB2[i]/N2[i])
      DCMt[i] = sqrt(40000/pi * AB2[i]/N2[i])

      Xt_bandera = 1-(DCMt[i]/DCMb[i])
      Xt[i] = 1-(DCMt[i]/DCMb[i])

      # Aplicación del modelo CON ACLAREO O SIN ACLAREO
    }else{
      # Recuperar valores de ABt y ABb
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
