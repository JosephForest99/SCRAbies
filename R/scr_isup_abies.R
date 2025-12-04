#' @title Sistema de crecimiento con el Índice de Supresión de Pienaar (1979)
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
#' # scr_isup_abies(E = 1:5, Narb = 2000, IS = 28)
#'
#' @import dplyr
#' @import tidyr
#'
#' @export

scr_isup_abies <- function(E, Narb, IS, Thinning = FALSE, Ets = NULL, ABpts = NULL){

  # Función negativa
  `%nin%` = Negate(`%in%`)

  if (isTRUE(Thinning)) {
    if (length(Ets) != length(ABpts)) {
      stop("Los vectores `Ets` y `ABpts` deben tener la misma longitud")
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
  N2 = NA; Nu2 = NA
  N2[1] = Narb
  Nu2[1] = Narb
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

      # Estimación de N de una parcela sin aclareo
      Nu2[i] = N_fn(N1 = Nu2[i-1], E1 = E[i-1], E2 = E[i])
      # Disminución de árboles de acuerdo a la proporción de AB asignado
      N2[i] = N2[i-1] - N2[i-1] * ABpt[i]**(1/1.649317)
    }else{
      Nu2[i] = N_fn(N1 = Nu2[i-1], E1 = E[i-1], E2 = E[i])
      N2[i] = N_fn(N1 = N2[i-1], E1 = E[i-1], E2 = E[i])

      # Banderas (para guardar el último valor)
      Nt[i] = Nt_bandera
      Nb[i] = Nb_bandera
      ABpt[i] = ABpt_bandera
    }
  }

  # # **ÁREA BASAL**

  # Predicción
  AB1_fn <- function(E1, N1, AD1){
    # Parámetros del modelo
    B0 = -5.9580735; B1 = 1.3068509; B2 = 0.6145979; B3 = 0.4594637
    # Modelo de predicción
    AB1 = exp(B0)*AD1**B1*E1**B2*N1**B3

    return(AB1)
  }
  # Proyección
  AB2_fn <- function(AB1, E1, E2, N1, N2, AD1, AD2){
    # Parámetros del modelo
    B0 = -5.9580735; B1 = 1.3068509; B2 = 0.6145979; B3 = 0.4594637
    # Parámetros del modelo (EDAD > 90 AÑOS)
    A1 = 0.941209
    # Modelo de proyección
    # Edad < 90 años
    if (E2<90) {
      AB2 = AB1*(AD2/AD1)**B1*(E2/E1)**B2*(N2/N1)**B3
      # Edad > 90 años
    } else{
      AB2 = AB1*(AD2/AD1)^A1
    }

    return(AB2)
  }

  # Índice de supresión
  ISup_fn <- function(ISup1, E1, E2){
    # Parámetros del modelo
    B1 = -1.482653
    ISup2 = ISup1 * exp(B1/E2 * (E2 - E1))

    return(ISup2)
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

  ISup = vector(mode = "numeric", length = length(E));

  AB2 = NA; ABu2 = NA
  AB2[1] = AB1_fn(E1 = E[1], N1 = N2[1], AD1 = AD2[1])
  ABu2[1] = AB1_fn(E1 = E[1], N1 = N2[1], AD1 = AD2[1])
  DCM[1] = sqrt(40000/pi * AB2[1]/N2[1])

  for (i in 2:length(E)) {

    # Aplicación del modelo SIN ACLAREO
    if (Et[i] == 0) {
      ABu2[i] = AB2_fn(AB1 = ABu2[i-1], E1 = E[i-1], E2 = E[i],
                       N1 = Nu2[i-1], N2 = Nu2[i], AD1 = AD2[i-1], AD2 = AD2[i])
      AB2[i] = AB2_fn(AB1 = AB2[i-1], E1 = E[i-1], E2 = E[i],
                      N1 = N2[i-1], N2 = N2[i], AD1 = AD2[i-1], AD2 = AD2[i])

      DCM[i] = sqrt(40000/pi * AB2[i]/N2[i])

      # AÑO (s) de aclareo (s) - Solo se considera la diferencia de ABb - ABt
    }else if (E[i] == E[i-1]){
      ABt_bandera = AB2[i-1] * ABpt[i]
      ABt[i] = AB2[i-1] * ABpt[i]

      ABb_bandera = AB2[i-1]
      ABb[i] = AB2[i-1]

      # Estimación de ABu (sin aclareo)
      ABu2[i] = AB2_fn(AB1 = ABu2[i-1], E1 = E[i-1], E2 = E[i],
                       N1 = Nu2[i-1], N2 = Nu2[i], AD1 = AD2[i-1], AD2 = AD2[i])

      # ESTIMACIÓN DE AB COMO DIFERENCIA (ABb-ABt)
      AB2[i] = AB2[i-1] - ABt[i]

      # Estimación del Índice de supresión
      ISup[i] = (ABu2[i] - AB2[i])/ABu2[i]

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

      # Estimación de ABu (sin aclareo)
      ABu2[i] = AB2_fn(AB1 = ABu2[i-1], E1 = E[i-1], E2 = E[i],
                       N1 = Nu2[i-1], N2 = Nu2[i], AD1 = AD2[i-1], AD2 = AD2[i])

      # Estimación del índice de supresión
      ISup[i] = ISup_fn(ISup[i-1], E[i-1], E[i])

      # Estimación de AB2 (con aclareo)
      AB2[i] = ABu2[i] * (1-ISup[i])

      DCM[i] = sqrt(40000/pi * AB2[i]/N2[i])

    }
  }

  df <- tibble(Thin = Thin, E = E, Densidad = Narb, IS = IS, AD = AD2,
               N = N2, ISup = ISup, AB = AB2, DCM = DCM
  )
  return(df)

}
