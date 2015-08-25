library(RGtk2)
library(cairoDevice)
library(Rgraphviz)
           
# tipi di controlli: 1 = testo, 2 = lista, 3 = testo multiplo, 4 = checkbox
TESTO = 1
LISTA = 2
VETT = 3
CHECKBOX = 4

# etichette
lbls <- list(
  c("alpha", "vector of parameters for the Activation sigmoid function"),
  c("beta", "vector of parameters for the Activation sigmoid function"),
  c("lambda", "vector of time constants influencing both the rate of transcription and the spontaneous degradation term"),
  c("Xmin", "vector of minimum level of expression of genes"),
  c("Xmax", "vector of maximum level of expression of genes"),
  c("X0", "initial conditions, i.e. gene expression values at time 0 scaled between 0 and 1")
)

calc_dll <-
    function(nome, args)
{
    stopifnot(is.character(nome))
    stopifnot(is.vector(args))
    vett <- .Call("calcv",
                 nome,
                 args
					)
    return(vett)
}

interf_s <- new.env()

# costanti simboliche
WEIGHTS = 1
RULES = 2
FPR_S = 7
EXTF = 8
EXTIN = 9
ACTF = 10
KO = 11
ALPHA = 12
TIMES = 60
STAT_T = 64
STAT_W = 65
METHOD = 66
PARAM = 69
SIM = 67
SAVE= 68

TOT_ARG_s = 20

non.vuota_ctrl_s <- function(widget, event, userData)
{
  if (gtkEntryGetText(widget) == "") {
    message_dialog_s("the weights matrix must be specified")
    pad_w <- get("pad_w", interf_s)
    widget$setText(sprintf("weights%0*d.txt", pad_w, 1))
    widget$grabFocus()
  }
  return (FALSE)
}

pars_ctrl_s <- function(widget, event, userData)
{
  txt <- gtkEntryGetText(widget)
  ris <- try(calc_dll(txt, c(0)))
  if (inherits(ris, "try-error")) {
    message_dialog_s(geterrmessage())
    widget$setText("0.0")
    return (FALSE)
  }
  return (FALSE)
}

vett_ctrl_s <- function(widget, event, userData)
{
  v <- sprintf("c(%s)", gtkEntryGetText(widget))
  ris <- try(eval(parse(text=v)))
  if (inherits(ris, "try-error")) {
    message_dialog_s("the value must be a list of numbers separated by commas: set to '1, 2'")
    widget$setText("1, 2")
    widget$grabFocus()
  }
  return (FALSE)
}

lstval_ctrl_s <- function(widget, event, userData)
{
  txt <- gtkEntryGetText(widget)
  # "g1=0.5, g2=1"
  coppie <- try(strsplit(txt, ","))
  if (inherits(coppie, "try-error")) {
    message_dialog_s("the value must be a list like 'g1=0.5, g2=1': resetted to the example")
    widget$setText("g1=0.5, g2=1")
    widget$grabFocus()
    return (FALSE)          
  }
  for (cp in coppie[[1]]) {
    coppia <- strsplit(cp, "=")
    elem <- coppia[[1]]
    if (length(elem) != 2 || inherits(coppia, "try-error")) {
      message_dialog_s("the value must be a list like 'g1=0.5, g2=1': resetted to the example")
      widget$setText("g1=0.5, g2=1")
      widget$grabFocus()
      return (FALSE)
    }
    if (length(grep("g\\d+", elem[1], perl=TRUE)) == 0) {
      message_dialog_s("the first element of the pair must be named gi: resetted to the example")
      widget$setText("g1=0.5, g2=1")
      widget$grabFocus()
      return (FALSE)    
    }
    n <- as.double(elem[2])
    if (is.na(n)) {
      message_dialog_s("the second element of the pair must be a number: resetted to the example")
      widget$setText("g1=0.5, g2=1")
      widget$grabFocus()
      return (FALSE)    
    }
  }
  return (FALSE)
}

num_ctrl_s <- function(widget, event, userData)
{
  num <- as.double(gtkEntryGetText(widget))
  if (is.na(num)) {
    message_dialog_s("the value must be a number: set to 0")
    widget$setText(0.0)
    widget$grabFocus()
  }
  return (FALSE)
}

gz_ctrl_s <- function(widget, event, userData)
{
  num <- as.double(gtkEntryGetText(widget))
  if (is.na(num)) {
    message_dialog_s("the value must be specified: set to 5")
    widget$setText(5)
    return (FALSE)
  }
  else if (num <= 0) {
    message_dialog_s("the value must be greater than zero")
    widget$setText(abs(num))
    widget$grabFocus()
  }
  return (FALSE)
}

gt_ctrl_s <- function(widget, event, userData)
{
  rif <- as.double(gtkEntryGetText(userData))
  if (as.double(gtkEntryGetText(widget)) <= rif) {
    message_dialog_s(paste("the value must be greater than", rif))
    widget$setText(rif + 1)
    widget$grabFocus()
  }
  return (FALSE)
}

message_dialog_s <- function(testo)
{
  w <- get("main", interf_s)
  dialog <- gtkMessageDialogNew(w, c("modal", "destroy-with-parent"), "error", "ok", testo)
  dialog$run()
  dialog$destroy()
}

# valutazione di espressioni
run_s <- function(code) {
  #e <- try(parse(text=paste("x <- (0 : 100) / 100;", code)), TRUE)
#  if (inherits(e, "try-error")) {
#    message_dialog(geterrmessage())
#    return()
#  }
#  ris <- try(eval(e), TRUE)
#  if (inherits(ris, "try-error")) {
#    message_dialog(geterrmessage())
#    return()
#  }
  x <- (0 : 100) / 100
  ris <- try(calc_dll(code, x))
  if (inherits(ris, "try-error")) {
    message_dialog_s(geterrmessage())
    return()
  }
  return (ris)
}

# funzione di call-back per il pulsante "Funzione"
preview_s <- function(widget, userData)
{
  x <- (0 : 100) / 100
  y <- run_s(gtkEntryGetText(userData))
  plot(x, y, type="l")
}

act.fun_cb_s <- function(widget, event, userData)
{
  alpha <- get("alpha", interf_s)
  beta <- get("beta", interf_s)
  txt <- gtkComboBoxGetActiveText(widget)
  if (txt == "sigmoidal") {
    gtkWidgetSetSensitive(alpha, TRUE)
    gtkWidgetSetSensitive(beta, TRUE)
  }
  else {
    gtkWidgetSetSensitive(alpha, FALSE)
    gtkWidgetSetSensitive(beta, FALSE)
  }
}

toggle_cb_s <- function(widget, userData)
{
  if (gtkToggleButtonGetActive(widget)) {
    assign(paste("toggle", userData[1], sep="_"), userData[2], interf_s)
    #cat(sprintf("%s = %s\n", paste("toggle", userData[1], sep="_"), userData[2]), "\n")
  }
}

toggle1_cb_s <- function(widget, userData)
{
  if (gtkToggleButtonGetActive(widget)) {
    if (gtkToggleButtonGetActive(userData[3][[1]]))
      gtkToggleButtonSetActive(userData[3][[1]], FALSE)
    assign(paste("toggle1", userData[1], sep="_"), userData[2], interf_s)
    #cat(sprintf("%s = %s\n", paste("toggle", userData[1], sep="_"), userData[2]), "\n")
  }
}

vuota <- list("", "")
ef <- rep(list(vuota), 100)
ef[[1]] <- list("sin(x * 10)", "g1=0.5, g2=1")
ef[[2]] <- list("cos(x * 10)", "g1=1")

# variabili globali
glbls_s <- list(
  # nome, tipo, default, on_exit, descr
  list("W_file", TESTO, "", NULL, "file with the connectivity matrix defining network topology (void for all)"),
  # nome, tipo, default, valori, call_back, descr
  list("r_file", TESTO, "", NULL, "rules from a file"),
  list("", TESTO, 12, NULL, "non usato"),
  list("", TESTO, 2.2, NULL, "non usato"),
  list("", LISTA, 1, c("free", "out"), NULL, "non usato"),
  list("", TESTO, 0.4, NULL, "non usato"),
  # nome, tipo, defaults
  # 7
  list("f.pr.and", TESTO, "0.5", pars_ctrl_s, "function with domain and codomain in [0-1] that expresses the probability to obtain a cooperative rather than a synergic rule"),
  list("ext.fun", TESTO, ef[[1]][[1]], NULL, "external stimulus"),
  list("ext.in", TESTO, ef[[1]][[2]], NULL, "genes to apply the external stimulus to, as a list whose elements are in the form 'gn=weight'"),
  list("act.fun", LISTA, 2, c("linear", "sigmoidal"), act.fun_cb_s, "the activation function"),
  list("ko.experiment", TESTO, "", vett_ctrl_s, "set of genes, separated by commas"),
  # 12
  list("alpha_const", TESTO, 1, gz_ctrl_s, "constant value (ignored when the checkboxes below are selected)"), # vector of parameters of the Activation sigmoid function
  list("alpha_a", TESTO, 0, num_ctrl_s, "parameter 'a' for the uniform distribution"),
  list("alpha_b", TESTO, 1, num_ctrl_s, "parameter 'b' for the uniform distribution"),
  list("alpha_nm", TESTO, 10, gz_ctrl_s, "mean for the normal distribution"),
  list("alpha_nsd", TESTO, 0.2, gz_ctrl_s, "std. dev. for the normal distribution"),
  list("alpha_lm", TESTO, 1, gz_ctrl_s, "mean for the log-normal distribution"),
  list("alpha_lsd", TESTO, 0, gz_ctrl_s, "std. dev. for the log-normal distribution"),
  list("alpha_file", TESTO, NULL, NULL, "file from which to read the desired values"),
  # 20
  list("beta_const", TESTO, 1, gz_ctrl_s, "constant value (ignored when the checkboxes below are selected)"),
  list("beta_a", TESTO, 0, num_ctrl_s, "parameter 'a' for the uniform distribution"),
  list("beta_b", TESTO, 1, num_ctrl_s, "parameter 'b' for the uniform distribution"),
  list("beta_nm", TESTO, 0.5, gz_ctrl_s, "mean for the normal distribution"),
  list("beta_nsd", TESTO, 0.01, gz_ctrl_s, "std. dev. for the normal distribution"),
  list("beta_lm", TESTO, 1, gz_ctrl_s, "mean for the log-normal distribution"),
  list("beta_lsd", TESTO, 0, gz_ctrl_s, "std. dev. for the log-normal distribution"),
  list("beta_file", TESTO, NULL, NULL, "file from which to read the desired values"),
  # 28
  list("lambda_const", TESTO, 1, gz_ctrl_s, "constant value (ignored when the checkboxes below are selected)"),
  list("lambda_a", TESTO, 0, num_ctrl_s, "parameter 'a' for the uniform distribution"),
  list("lambda_b", TESTO, 1, num_ctrl_s, "parameter 'b' for the uniform distribution"),
  list("lambda_nm", TESTO, 1, gz_ctrl_s, "mean for the normal distribution"),
  list("lambda_nsd", TESTO, 0.1, gz_ctrl_s, "std. dev. for the normal distribution"),
  list("lambda_lm", TESTO, 1, gz_ctrl_s, "mean for the log-normal distribution"),
  list("lambda_lsd", TESTO, 0, gz_ctrl_s, "std. dev. for the log-normal distribution"),
  list("lambda_file", TESTO, NULL, NULL, "file from which to read the desired values"),
  # 36
  list("Xmin_const", TESTO, 0, gz_ctrl_s, "constant value (ignored when the checkboxes below are selected)"),
  list("Xmin_a", TESTO, 0, num_ctrl_s, "parameter 'a' for the uniform distribution"),
  list("Xmin_b", TESTO, 1, num_ctrl_s, "parameter 'b' for the uniform distribution"),
  list("Xmin_nm", TESTO, 1, gz_ctrl_s, "mean for the normal distribution"),
  list("Xmin_nsd", TESTO, 0, gz_ctrl_s, "std. dev. for the normal distribution"),
  list("Xmin_lm", TESTO, 1, gz_ctrl_s, "mean for the log-normal distribution"),
  list("Xmin_lsd", TESTO, 0, gz_ctrl_s, "std. dev. for the log-normal distribution"),
  list("Xmin_file", TESTO, NULL, NULL, "file from which to read the desired values"),
  # 44
  list("Xmax_const", TESTO, 10, gz_ctrl_s, "constant value (ignored when the checkboxes below are selected)"),
  list("Xmax_a", TESTO, 0, num_ctrl_s, "parameter 'a' for the uniform distribution"),
  list("Xmax_b", TESTO, 1, num_ctrl_s, "parameter 'b' for the uniform distribution"),
  list("Xmax_nm", TESTO, 1, gz_ctrl_s, "mean for the normal distribution"),
  list("Xmax_nsd", TESTO, 0, gz_ctrl_s, "std. dev. for the normal distribution"),
  list("Xmax_lm", TESTO, 1, gz_ctrl_s, "mean for the log-normal distribution"),
  list("Xmax_lsd", TESTO, 0, gz_ctrl_s, "std. dev. for the log-normal distribution"),
  list("Xmax_file", TESTO, NULL, NULL, "file from which to read the desired values"),
  # 52
  list("X0_const", TESTO, 1, gz_ctrl_s, "constant value (ignored when the checkboxes below are selected)"),
  list("X0_a", TESTO, 0, num_ctrl_s, "parameter 'a' for the uniform distribution"),
  list("X0_b", TESTO, 1, num_ctrl_s, "parameter 'b' for the uniform distribution"),
  list("X0_nm", TESTO, 1, gz_ctrl_s, "mean for the normal distribution"),
  list("X0_nsd", TESTO, 0, gz_ctrl_s, "std. dev. for the normal distribution"),
  list("X0_lm", TESTO, 1, gz_ctrl_s, "mean for the log-normal distribution"),
  list("X0_lsd", TESTO, 0, gz_ctrl_s, "std. dev. for the log-normal distribution"),
  list("X0_file", TESTO, NULL, NULL, "file from which to read the desired values"),
  # 60
  list("from", TESTO, "0", num_ctrl_s, "from"),
  list("to", TESTO, "5", num_ctrl_s, "to"),
  list("step", TESTO, "0.05", num_ctrl_s, "step"),
  list("file_times", TESTO, "", NULL, "file from which to read the desired values"),
  list("stat_t", TESTO, "0.001", num_ctrl_s, "absolute difference under which the signal can be considered stationary"),
  list("stat_w", TESTO, "0", num_ctrl_s, "percentage of the maximum sampling time during which the signal has to satisfy the stationarity hypothesis (zero to disable the stationarity check)"),
  # 66
  list("method", LISTA, 2, c("Euler",
  "Embedded Runge-Kutta-Fehlberg (4,5)",
  "Embedded Runge-Kutta Cash-Karp (4,5)",
  "Embedded Runge-Kutta Prince-Dormand (8,9)",
  "Implicit 2nd order Runge-Kutta at Gaussian points",
  "Implicit 4th order Runge-Kutta at Gaussian points",
  "Implicit Gear method (M = 1)",
  "Implicit Gear method (M = 2)"),
  NULL, "Method used to solve differential equations"),
  # 67
  list("simulations", TESTO, "3", gz_ctrl_s, "number of simulations"),
  list("save", CHECKBOX, TRUE, NULL, "save simulations on disk"),
  # 69
  list("param", VETT, NULL, "non visuale")
)

titoli_s <- c("Network", "Rules", "Parameters", "External stimuli", "Knock out", "Simulations", "Done")

non_vuote <- function(ls, parte)
{
  ris <- list()
  for (l in ls) {
    if (l[[parte]] != "")
      ris <- c(ris, l[[parte]])
  }
  return (ris)
}

riepilogo_s <- function(table1, num, num1, x, y)
{
	# "primo_s" mi dice se e' la prima volta che ci passo, nel qual caso devo prima creare le liste in cui memorizzera' i valori
  primo_s <- get("primo_s", interf_s)
  argom_s <- get("argom_s", interf_s)
  N <- get("N", interf_s)
  N <- 5
  if (num != PARAM)
    ctrl <- get(glbls_s[[num]][[1]], interf_s)
  txt_lbl <- glbls_s[[num]][[1]]
  if (num == WEIGHTS) {
    txt <- get("toggle_w", interf_s)
    if (txt == 1)
      txt1 <- "NULL"
    else {
      str <- gsub("\\", "/", gtkEntryGetText(get(glbls_s[[WEIGHTS]][[1]], interf_s)), fixed=TRUE)
      txt1 <- sprintf("read.table('%s')", str)
    }
    if (primo_s)
      argom_s[[num1]] <- list("espressione", txt1, NULL)
    else
      argom_s[[num1]][[2]] <- txt1
  }
  else if (num >= ALPHA && num < TIMES) {
    l <- (num - 4) / 8
    txt_lbl <- lbls[[l]][1]
    attv <- get(paste("toggle", lbls[[l]][1], sep="_"), interf_s)
    #cat(attv)
    if (attv == 1) # costante
      txt1 <- sprintf("rep(%s, %s)", gtkEntryGetText(get(glbls_s[[num]][[1]], interf_s)), N)
    else if (attv == 2) { # uniforme
      a <- gtkEntryGetText(get(glbls_s[[num + 1]][[1]], interf_s))
      b <- gtkEntryGetText(get(glbls_s[[num + 2]][[1]], interf_s))
      txt1 <- sprintf("c(%s, %s)", a, b)
    }
    else if (attv == 3) { #normale
      m <- gtkEntryGetText(get(glbls_s[[num + 3]][[1]], interf_s))
      s <- gtkEntryGetText(get(glbls_s[[num + 4]][[1]], interf_s))
      txt1 <- sprintf("c(%s, %s)", m, s)
    }
    else if (attv == 4) { # log-normale
      m <- gtkEntryGetText(get(glbls_s[[num + 5]][[1]], interf_s))
      s <- gtkEntryGetText(get(glbls_s[[num + 6]][[1]], interf_s))
      txt1 <- sprintf("c(%s, %s)", m, s)
    }
    else
      stop(sprintf("valore di 'attv' non valido (%d)!", attv))
    if (get(paste("toggle1", lbls[[l]][1], sep="_"), interf_s) == 2) {
      str <- gsub("\\", "/", gtkEntryGetText(get(glbls_s[[num + 7]][[1]], interf_s)), fixed=TRUE)
      if (attv == 1)
        txt1 <- sprintf("read.table('%s')", str)
      else
        txt1 <- sprintf("c(read.table('%s')[[1]], %s)", str, txt1)
    }
    if (primo_s)
      argom_s[[num1]] <- list("espressione", txt1, NULL)
    else
      argom_s[[num1]][[2]] <- txt1
  }
  else if (num == RULES) {
    txt_lbl <- "Rules"
    txt <- get("toggle_r", interf_s)
    if (txt == 3) {
      str <- gsub("\\", "/", gtkEntryGetText(get(glbls_s[[2]][[1]], interf_s)), fixed=TRUE)
      txt1 <- sprintf("read.table('%s')", str)
      if (primo_s)
        argom_s[[num1]] <- list("espressione", txt1, NULL)
      else
        argom_s[[num1]][[2]] <- txt1
    }
    else if (txt == 2) {
      txt1 <- "NULL"
      if (primo_s)
        argom_s[[num1]] <- list("espressione", txt1, NULL)
      else
        argom_s[[num1]][[2]] <- txt1
    }
    else {
      txt1 <- "list()"
      if (primo_s)
        argom_s[[num1]] <- list("stringa", txt1, NULL)
      else
        argom_s[[num1]][[2]] <- txt1
    }
  }
  else if (num == FPR_S) {
    txt <- get("toggle_r", interf_s)
    if (txt == 1) {
      txt1 <- gtkEntryGetText(get(glbls_s[[num]][[1]], interf_s))
      if (primo_s) {
        argom_s[[2]][[1]] <- "stringa"
        argom_s[[2]][[2]] <- "list()"
        if (length(argom_s[[2]]) < 3)
            argom_s[[2]][[3]] <- NULL
        if (is.null(argom_s[[num1]]))
          argom_s[[num1]] <- list("stringa", txt1, NULL)
        else {
          argom_s[[num1]][[1]] <- "stringa"
          argom_s[[num1]][[2]] <- txt1
        }
      }
      else {
        argom_s[[2]][[2]] <- "list()"
        argom_s[[num1]][[2]] <- txt1
      }
    }
    else if (txt == 2) {
      txt1 <- "NULL"
      if (primo_s) {
        argom_s[[2]][[1]] <- "espressione"
        argom_s[[2]][[2]] <- "NULL"
        if (length(argom_s[[2]]) < 3)
          argom_s[[2]][[3]] <- NULL
        if (is.null(argom_s[[num1]]))
          argom_s[[num1]] <- list("espressione", NULL, NULL)
        else {
          argom_s[[num1]][[1]] <- "espressione"
          argom_s[[num1]][[2]] <- NULL
        }
      }
      else {
        argom_s[[2]][[2]] <- NULL
        argom_s[[num1]][[2]] <- NULL
      }
    }
    else if (txt == 3) {
      txt1 <- "list()"
      if (primo_s) {
        argom_s[[2]][[1]] <- "stringa"
        argom_s[[2]][[2]] <- "list()"
        if (length(argom_s[[2]]) < 3)
          argom_s[[2]][[3]] <- NULL
        if (is.null(argom_s[[num1]]))
          argom_s[[num1]] <- list("espressione", NULL, NULL)
        else {
          argom_s[[num1]][[1]] <- "espressione"
          argom_s[[num1]][[2]] <- NULL
        }
      }
      else {
        argom_s[[2]][[2]] <- list()
        argom_s[[num1]][[2]] <- NULL
      }
    }
  }
  else if (num == EXTF) {
    txt <- get("toggle_ext.fun", interf_s)
    txt1 <- list()
    if (txt == 1) {
      ef <- get("ef", interf_s)
      txt1 <- non_vuote(ef, 1)
    }
    if (primo_s)
      argom_s[[num1]] <- list("stringa", txt1, NULL)
    else
      argom_s[[num1]][[2]] <- txt1
  }
  else if (num == EXTIN) {
    txt <- get("toggle_ext.fun", interf_s)
    n <- 0
    if (txt == 1) {
      ef <- get("ef", interf_s)
      l <- non_vuote(ef, 2)
      n <- length(l)
    }
    if (n == 0)
      txt1 <- NA
    else {
      txt1 <- matrix(0, nrow=N, ncol=n)
      for (i in 1:n) {
        expr1 <- sprintf("list(%s)", ef[[i]][[2]])
        if (!is.null(expr1)) # altrimenti compare il prompt, NON NULL!
          expr <- try(parse(text=expr1))
        else
          expr <- NULL
        val <- eval(expr)
        for (j in 1:N) {
          g <- paste("val$g", j, sep="")
          if (!is.null(g)) # altrimenti compare il prompt, NON NULL!
            val1 <- eval(try(parse(text=g)))
          else
            g <- NULL
          if (!is.null(val1))
            txt1[j, i] <- val1
        }
      }
    }
    if (primo_s)
      argom_s[[num1]] <- list("stringa", txt1, NULL)
    else
      argom_s[[num1]][[2]] <- txt1
  }
  else if (num == KO) {
    txt <- get("toggle_ko", interf_s)
    if (txt == 3)
      txt1 <- "c(0)"
    else if (txt == 2)
      txt1 <- sprintf("c(%s)", gtkEntryGetText(get(glbls_s[[KO]][[1]], interf_s)))
    else
      txt1 <- "NULL"
    if (primo_s)
      argom_s[[num1]] <- list("espressione", txt1, NULL)
    else
      argom_s[[num1]][[2]] <- txt1
  }
  else if (num >= TIMES && num < METHOD - 2) {
    txt_lbl <- "times"
    txt <- get("toggle_times", interf_s)
    if (txt == 1) {
      da <- gtkEntryGetText(get(glbls_s[[TIMES]][[1]], interf_s))
      a <- gtkEntryGetText(get(glbls_s[[TIMES + 1]][[1]], interf_s))
      passo <- gtkEntryGetText(get(glbls_s[[TIMES + 2]][[1]], interf_s))
      txt1 <- sprintf("seq(%s, %s, %s)", da, a, passo)
    }
    else {
      str <-gsub("\\", "/", gtkEntryGetText(get(glbls_s[[TIMES + 3]][[1]], interf_s)), fixed=TRUE)
      txt1 <- sprintf("read.table('%s')", str)
    }
    if (primo_s)
      argom_s[[num1]] <- list("espressione", txt1, NULL)
    else
      argom_s[[num1]][[2]] <- txt1
  }
  else if (num == PARAM) {
    # l'ordine vero e` lambda, alpha, beta, ...
    txt <- get(paste("toggle", lbls[[3]][1], sep="_"), interf_s)
    n <- as.integer(txt) - 1
    if (get(paste("toggle1", lbls[[3]][1], sep="_"), interf_s) == 1) {
      n <- n + 7
      if (n == 7)
        n <- 4
    }
    else if (get(paste("toggle1", lbls[[3]][1], sep="_"), interf_s) == 2) {
      n <- n + 4
      if (n == 4)
        n <- 0
    }
    txt1 <- paste("c", n, sep="(")
    for (i in c(1,2,4,5,6)) {
      txt <- get(paste("toggle", lbls[[i]][1], sep="_"), interf_s)
      n <- as.integer(txt) - 1
      if (get(paste("toggle1", lbls[[i]][1], sep="_"), interf_s) == 1) {
        n <- n + 7
         if (n == 7)
           n <- 4
      }
      else if (get(paste("toggle1", lbls[[i]][1], sep="_"), interf_s) == 2) {
        n <- n + 4
        if (n == 4)
          n <- 0
      }
      txt1 <- paste(txt1, n, sep=",")
    }
    txt1 <- paste(txt1, ")", sep="")
    if (primo_s)
      argom_s[[num1]] <- list("espressione", txt1, NULL)
    else
      argom_s[[num1]][[2]] <- txt1
  }
  else if (gtkWidgetIsSensitive(ctrl)) {
    if (glbls_s[[num]][[2]] == 1) {
      txt1 <- as.numeric(gtkEntryGetText(ctrl))
      if (primo_s)
        argom_s[[num1]] <- list("numero", sprintf('%3.3f', txt1), NULL)
      else
        argom_s[[num1]][[2]] <- sprintf('%3.3f', txt1)
    }
    else if (glbls_s[[num]][[2]] == LISTA) {
      if (num == METHOD) {
        sigle <- c("Euler", "rkf45", "rkck", "rk8pd", "rk2imp", "rk4imp", "gear1", "gear2")
        txt1 <- sigle[gtkComboBoxGetActive(ctrl) + 1]
      }
      else
        txt1 <- gtkComboBoxGetActiveText(ctrl)
      if (primo_s)
        argom_s[[num1]] <- list("stringa", as.character(txt1), NULL)
      else
        argom_s[[num1]][[2]] <- as.character(txt1)
    }
    else if (glbls_s[[num]][[2]] == CHECKBOX) {
      txt1 <- gtkToggleButtonGetActive(ctrl)
      if (primo_s)
        argom_s[[num1]] <- list("booleano", as.logical(txt1), NULL)
      else
        argom_s[[num1]][[2]] <- as.logical(txt1)
    }
    else {
      txt1 <- "-"
      if (primo_s)
        argom_s[[num1]] <- list("saltato", glbls_s[[num]][[3]], NULL)
    }
  }
  else {
    txt1 <- "-"
    if (primo_s)
      argom_s[[num1]] <- list("saltato", glbls_s[[num]][[3]], NULL)
  }
  if (primo_s) {
    label1 <- gtkLabelNew()
    label1$setAlignment(0, 0.5)
    gtkLabelSetText(label1, txt_lbl)
    table1$Attach(label1, 4 * x, 4 * x + 1, y, y + 1)
    label2 <- gtkLabelNew()
    label2$setAlignment(1, 1)
    argom_s[[num1]][[3]] <- label2
    table1$Attach(label2, 4 * x + 2, 4 * x + 3, y, y + 1)
  }
  else {
    label2 <- argom_s[[num1]][[3]]
  }
  if (num == EXTF) {
    if (get("toggle_ext.fun", interf_s) == 1)
      gtkLabelSetText(label2, "list(...)")
    else
      gtkLabelSetText(label2, "list()")
  }
  else if (num == EXTIN && get("toggle_ext.fun", interf_s) == 1)
    gtkLabelSetText(label2, "matrix")
  else
    gtkLabelSetText(label2, txt1)
  assign("argom_s", argom_s, interf_s)
  return
}

on_assistant_apply_s <- function(widget, data)
{
  a <- list()
  argom_s <- get("argom_s", interf_s)
  #argom[[16]][[2]] <- sprintf("%s / %s", argom_s[[16]][[2]], argom_s[[15]][[2]])
  for (l in 1:TOT_ARG_s) {
    if (argom_s[[l]][[1]] == "espressione" || argom_s[[l]][[1]] == "numero") {
       if (!is.null(argom_s[[l]][[2]])) # altrimenti compare il prompt, NON NULL!
        expr <- try(parse(text=argom_s[[l]][[2]]))
       else       
        expr <- NULL
       a[[l]] <- eval(expr)
    }
    #else if (argom_s[[l]][[1]] == "funzione") {
    #   a[[l]] <- eval(try(parse(text=paste("function(x) { return(", argom_s[[l]][[2]], ")}"))))
    #}
    else
      a[[l]] <- argom_s[[l]][[2]]
  }
  salvati <- as.logical(a[[TOT_ARG_s - 1]])
  if (!is.null(a[[1]]))
    a[[1]] <- as.matrix(a[[1]]) # weights
  if (!is.null(a[[2]])) {
    if (a[[2]] == "list()")
      a[[2]] <- list() # rnd
    else
      a[[2]] <- as.list(a[[2]]) # rules
  }
  if (is.list(a[[4]]))
    a[[4]] <- a[[4]][[1]]
  if (is.list(a[[5]]))
    a[[5]] <- a[[5]][[1]]
  if (is.list(a[[6]]))
    a[[6]] <- a[[6]][[1]]
  if (is.list(a[[7]]))
    a[[7]] <- a[[7]][[1]]
  if (is.list(a[[9]]))
    a[[9]] <- a[[9]][[1]]
  if (is.list(a[[10]]))
    a[[10]] <- a[[10]][[1]]
  if (is.list(a[[11]]))
    a[[11]] <- a[[11]][[1]]
  cat("\nArgs: \n", paste(a, "\n"), "\n")
  d <- dir()
  file.remove(d[grep("SIMdata_simulateprofiles.*\\.txt", d, perl=TRUE)])
  assign("res_netsim_simulated", do.call("simulateprofiles", a), globalenv())
  if (salvati)
    cat("\nResults have been saved in variable 'res_netsim_simulated' and in folder '", getwd(), "'\n")
  else
    cat("\nResults have been saved in variable 'res_netsim_simulated'\n")
}

on_assistant_prepare_s <- function(widget, page, data)
{
  current_page <- widget$getCurrentPage()
  n_pages <- widget$getNPages()
  if (current_page == 1) {
    pesi <- gtkEntryGetText(get(glbls_s[[1]][[1]], interf_s))
    if (pesi == "") {
    	pad_w <- get("pad_w", interf_s)
      str <- sprintf("weights%0*d.txt", pad_w, 1)
      if (!file.exists(str)) {
         message_dialog_s(paste("there are no files named 'weights*.txt' in the working directory!", get("main", interf_s)))
         get("main", interf_s)$destroy()
      }
    }
    else
      str <- gsub("\\", "/", pesi, fixed=TRUE)
    W <- read.table(str, sep="\t", na="NA", dec=".")
    assign("N", dim(W)[1] - 1, interf_s)
    gtkWidgetSetSensitive(get("rules_radio", interf_s), (get("toggle_w", interf_s) == 1))
  }
  if (current_page == 2) {
    pesi <- gtkEntryGetText(get(glbls_s[[1]][[1]], interf_s))
    for (sr in get("sn_radio", interf_s))
      gtkWidgetSetSensitive(sr, (get("toggle_w", interf_s) == 1))
  }
  if (current_page == 6) {
    k <- 1
    arg <- c(WEIGHTS, RULES, FPR_S, ACTF, ALPHA, ALPHA + 8, ALPHA + 2 * 8, ALPHA + 3 * 8, ALPHA + 4 * 8, ALPHA + 5 * 8, PARAM, KO, TIMES, STAT_T, STAT_W, METHOD, EXTIN, EXTF, SIM, SAVE)
    table1 <- get("riep", interf_s)
    for (j in 0:1) {
      for (i in 0:9) { # 9 righe
        if (k <= length(arg)) {
          riepilogo_s(table1, arg[k], k, j, i)
          k <- k + 1
        }
      }
    }
    assign("primo_s", FALSE, interf_s)
  }
  title <- sprintf("Netsim assistant - %s (%d of %d)", titoli_s[current_page + 1], current_page + 1, n_pages)
  widget$setTitle(title)
}

crea_testo_s <- function(l, tips_s)
{
  ctrl <- gtkEntryNew()
  gtkEntrySetAlignment(ctrl, 1)
  if (!is.null(l[[3]]))
  gtkEntrySetText(ctrl, l[[3]])
  if (!is.null(l[[4]]))
    gSignalConnect(ctrl, "focus-out-event", l[[4]])
  gtkTooltipsSetTip(tips_s, ctrl, l[[5]])
  assign(l[[1]], ctrl, interf_s)
}

crea_chk_s <- function(l, tips_s)
{
  ctrl <- gtkCheckButtonNew()
  assign(l[[1]], ctrl, interf_s)
  gtkToggleButtonSetActive(ctrl, l[[3]])
  gtkTooltipsSetTip(tips_s, ctrl, l[[4]])
  return (ctrl)
}

crea_vett_s <- function(l, tips_s)
{
  frm <- gtkVBoxNew(FALSE, 0)
  frm1 <- gtkHBoxNew(FALSE, 0)
  # etichetta
  lbl <- gtkFrameNew(l[[1]])
  frm$packStart(lbl, TRUE, TRUE, 0)
  # testi
  cz <- vector()
  i <- 1
  for (v in l[[3]]) {
    ctrl <- gtkEntryNew()
    gtkEntrySetAlignment(ctrl, 1)
    gtkEntrySetText(ctrl, v)
    gSignalConnect(ctrl, "focus-out-event", num_ctrl_s)
    gtkTooltipsSetTip(tips_s, ctrl, l[[4]][i])
    algn <- gtkAlignmentNew(0, 0, 0, 0)
    gtkAlignmentSetPadding(algn, 5, 5, 5, 5)
    algn$add(ctrl)
    frm1$packStart(algn, FALSE, TRUE, 0)
    cz <- c(cz, ctrl)
    i <- i + 1
  }
  lbl$add(frm1)
  assign(l[[1]], cz, interf_s)
  return (frm)
}

crea_lista_s <- function(l, tips_s)
{
  # lista combinata
  ctrl <- gtkComboBoxNewText()
  for (t in l[[4]])
    gtkComboBoxAppendText(ctrl, t)
  gtkComboBoxSetActive(ctrl, l[[3]] - 1)
  if (!is.null(l[[5]])) {
     gSignalConnect(ctrl, "changed", l[[5]])
  }
  gtkTooltipsSetTip(tips_s, ctrl, l[[6]])
  assign(l[[1]], ctrl, interf_s)
  return (ctrl)
}

crea_controllo_s <- function(num, etichetta, attivo)
{
  tips_s <- get("tips_s", interf_s)
  l <- glbls_s[[num]]
  if (etichetta) {
    frm <- gtkHBoxNew(FALSE, 0)
    # etichetta
    lbl <- gtkLabelNew(l[[1]])
    algn1 <- gtkAlignmentNew(0, 0, 0, 0)
    gtkAlignmentSetPadding(algn1, 5, 5, 5, 5)
    algn1$add(lbl)
  }
  if (l[[2]] == TESTO) {
    ctrl <- crea_testo_s(l, tips_s)
  }
  else if (l[[2]] == LISTA) {
    ctrl <- crea_lista_s(l, tips_s)
  }
  else if (l[[2]] == VETT) {
    ctrl <- crea_vett_s(l, tips_s)
  }
  else if (l[[2]] == CHECKBOX) {
    ctrl <- crea_chk_s(l, tips_s)
  }
  gtkWidgetSetSensitive(ctrl, attivo)
  if (etichetta) {
    algn2 <- gtkAlignmentNew(1, 0, 0, 0)
    gtkAlignmentSetPadding(algn2, 5, 5, 5, 5)
    algn2$add(ctrl)
    frm$packStart(algn1, TRUE, TRUE, 0)
    frm$packStart(algn2, TRUE, TRUE, 0)
    return (frm)
  }
  else {
    return (ctrl)
  }
}

preview1_s <- function(widget, userData)
{
  win <- get("superf_s", interf_s)
  w <- win[["allocation"]]$width
  h <- win[["allocation"]]$height
  nome <- get("W_file", interf_s)$getText()
  if (nome == "")
    return
  W <- read.table(nome, header=TRUE)
  n <- nrow(W)
  nAttrs <- list()
  eAttrs <- list()
  nomi <- as.character(lapply(seq(1, n), function(x) {sprintf("g%d", x)}))
  colori <- rep("white", n)
  stili <- vector()
  nAttrs$fillcolor <- colori
  names(nAttrs$fillcolor) <- nomi
  nomi1 <- vector()
  rEG  <-  new("graphNEL", nodes=nomi, edgemode="directed")
  for (i in 1:n) {
    for (j in 1:n) {
      if (W[i, j] != 0) {
        rEG <- addEdge(nomi[j], nomi[i], rEG, 1)
        if (W[i, j] > 0) {
          stili <- c(stili, "black")
      }
        else {
          stili <- c(stili, "red")
        }
        nomi1 <- c(nomi1, sprintf("%s~%s", nomi[j], nomi[i]))
    }
  }
  }
  eAttrs$color <- stili
  names(eAttrs$color) <- nomi1
  require(graph)  
  da2 <- win[["window"]]
  gc <- gdkGCNew(da2)
  pixmap <- gdkPixmapNew(NULL, w, h, depth=24)
  asCairoDevice(pixmap)
  plot(rEG, recipEdges="distinct", nodeAttrs=nAttrs, edgeAttrs=eAttrs)
  img <- gdkPixbufGetFromDrawable(NULL, pixmap, pixmap$getColormap(), 0, 0, 0, 0, w, h)
  assign("img_s", img, interf_s)
  gdkDrawPixbuf(da2, gc, pixbuf=img, 0, 0, 0, 0, w, h)
}

preview2_s <- function(widget, userData)
{
  nome <- get("W_file", interf_s)$getText()
  if (nome == "")
    return
  W <- read.table(nome, header=TRUE)
  n <- nrow(W)
  nomi <- as.character(lapply(seq(1, n), function(x) {sprintf("g%d", x)}))
  f <- file("tmp.dot", "w")
  win <- get("superf_s", interf_s)
  w <- win[["allocation"]]$width
  h <- win[["allocation"]]$height
  cat(file=f, sprintf("digraph G {\n\tratio=fill;\n\tsize=\"%3.3f,%3.3f\";\n", w / 75, h / 75))
  for (i in 1:n) {
    for (j in 1:n) {
      if (W[i, j] > 0) {
        cat(file=f, sprintf("\t%s -> %s;\n", nomi[j], nomi[i]))
      }
      else if (W[i, j] < 0) {
        cat(file=f, sprintf("\t%s -> %s [style=dashed];\n", nomi[j], nomi[i]))
      }
    }
  }
  cat(file=f, "}\n")
  close(f)
  system(sprintf("dot -Tpng -otmp.png tmp.dot"), ignore.stderr=FALSE, show.output.on.console=FALSE)
  img <- gdkPixbufNewFromFile("tmp.png")$retval
  da2 <- win[["window"]]
  gc <- gdkGCNew(da2)
  img1 <- gdkPixbufScaleSimple(img, w, h, "bilinear")
  assign("img_s", img1, interf_s)
  gdkDrawPixbuf(da2, gc, pixbuf=img1, 0, 0, 0, 0, w, h)
}

expose_event_cb <- function(widget, event, userData)
{
  #img <- gdkPixbufNewFromFile("tmp.png")$retval
  if (exists("img_s", interf_s)) {
    img <- get("img_s", interf_s)
    da2 <- widget[["window"]]
    w <- widget[["allocation"]]$width
    h <- widget[["allocation"]]$height
    gc <- gdkGCNew(da2)
    img1 <- gdkPixbufScaleSimple(img, w, h, "bilinear")
    gdkDrawPixbuf(da2, gc, pixbuf=img1, 0, 0, 0, 0, w, h) 
  }
}

create_page1_s <- function(assistant)
{
  box1 <- gtkVBox(FALSE, 12)
  box <- gtkHBox(FALSE, 2)
  box1$setBorderWidth(12)

  radiobutton1 <- gtkRadioButtonNewWithLabel(NULL, "Follow simulatenet networks")
  gSignalConnect(radiobutton1, "toggled", toggle_cb_s, c("w", 1))
  assign(paste("toggle", "w", sep="_"), 2, interf_s)
  box1$packStart(radiobutton1, FALSE, TRUE, 10)
  radiobutton1_group <- gtkRadioButtonGetGroup(radiobutton1)

  radiobutton2 <- gtkRadioButtonNewWithLabel(radiobutton1_group, "weights")
  gSignalConnect(radiobutton2, "toggled", toggle_cb_s, c("w", 2))
  gtkToggleButtonSetActive(radiobutton2, TRUE)
  box$packStart(radiobutton2, TRUE, TRUE, 0)

  entry1 <- crea_controllo_s(WEIGHTS, FALSE, TRUE)
  entry1$setEditable(FALSE)
  box$packStart(entry1, TRUE, TRUE, 0)
  assign("file_w", entry1, interf_s)

  button1 <- gtkButtonNewWithLabel("...")
  gSignalConnect(button1, "clicked", select_file_s, "W")
  box$packStart(button1, FALSE, TRUE, 0)

  button2 <- gtkButton("Draw")
  box$packStart(button2, FALSE, TRUE, 0)

  box1$packStart(box, FALSE, FALSE, 0)

  entry1 <- gtkDrawingArea()
  asCairoDevice(entry1)
  box1$packStart(entry1, TRUE, TRUE, 0)
  entry1$show()
  gtkWidgetSetSizeRequest(entry1, 300, 200)
  gSignalConnect(entry1, "expose-event", expose_event_cb)
  assign("superf_s", entry1, interf_s)
  if (get("vista_s", interf_s)) {
    gSignalConnect(button2, "clicked", preview2_s)
  }
  else {
    gSignalConnect(button2, "clicked", preview1_s)
  }

  assistant$appendPage(box1)
  assistant$setPageTitle(box1, "Topology: common parameters")
  #assistant$setPageType(box, "intro")
  assistant$setPageComplete(box1, TRUE)

  pixbuf <- assistant$renderIcon(GTK_STOCK_DIALOG_INFO, "dialog")
  assistant$setPageHeaderImage(box1, pixbuf)
}

create_page2_s <- function(assistant)
{
  box <- gtkVBox(FALSE, 12)
  box$setBorderWidth(12)

  entry1 <- gtkDrawingArea()
  asCairoDevice(entry1)
  entry1$show()
  gtkWidgetSetSizeRequest(entry1, 200, 300)
  #scale_cb <- function(range) { plot(1:20,(1:20)^range$getValue()) }
  box$packStart(entry1, FALSE, FALSE, 0)

  box2 <- gtkVBox(FALSE, 12)
  box1 <- gtkHBox(FALSE, 2)
  radiobutton1 <- gtkRadioButtonNewWithLabel(NULL, "random")
  gSignalConnect(radiobutton1, "toggled", toggle_cb_s, c("r", 1))
  box1$packStart(radiobutton1, TRUE, TRUE, 0)
  radiobutton1_group <- gtkRadioButtonGetGroup(radiobutton1)
  #radiobutton1$setGroup(radiobutton1_group)
  gtkToggleButtonSetActive(radiobutton1, TRUE)
  assign(paste("toggle", "r", sep="_"), 1, interf_s)
  #cat(sprintf("%s = %s\n", paste("toggle", "r", sep="_"), 1), "\n")

  entry2 <- crea_controllo_s(FPR_S, TRUE, TRUE)
  box1$packStart(entry2, FALSE, FALSE, 0)
  entry3 <- gtkButton("Draw")
  gSignalConnect(entry3, "clicked", preview_s, get("f.pr.and", interf_s))
  box1$packStart(entry3, FALSE, FALSE, 0)
  box$packStart(box1, FALSE, FALSE, 0)

  radiobutton2 <- gtkRadioButtonNewWithLabel(radiobutton1_group, "follow simulatenet rules")
  assign("rules_radio", radiobutton2, interf_s)
  gSignalConnect(radiobutton2, "toggled", toggle_cb_s, c("r", 2))
  box$packStart(radiobutton2, TRUE, TRUE, 0)
#  gtkRadioButtonSetGroup(radiobutton2, radiobutton1_group)

  box1 <- gtkHBox(FALSE, 12)
  radiobutton3 <- gtkRadioButtonNewWithLabel(radiobutton1_group, "file")
  gSignalConnect(radiobutton3, "toggled", toggle_cb_s, c("r", 3))
  box1$packStart(radiobutton3, TRUE, TRUE, 0)
#  gtkRadioButtonSetGroup(radiobutton3, radiobutton1_group)

  entry1 <- crea_controllo_s(RULES, FALSE, TRUE)
  gtkEntrySetWidthChars(entry1, 25)
  entry1$setEditable(FALSE)
  box1$packStart(entry1, TRUE, TRUE, 0)

  button1 <- gtkButtonNewWithLabel("...")
  gSignalConnect(button1, "clicked", select_file_s, "r")
  box1$packStart(button1, FALSE, FALSE, 0)

  box2$packStart(box1, TRUE, TRUE, 0)

  box$packStart(box2, TRUE, TRUE, 0)

  assistant$appendPage(box)
  #assistant$setPageType(box, "confirm")
  assistant$setPageComplete(box, TRUE)
  assistant$setPageTitle(box, "Rules")

  pixbuf <- assistant$renderIcon(GTK_STOCK_DIALOG_INFO, "dialog")
  assistant$setPageHeaderImage(box, pixbuf)
}

select_file_s <- function(widget, userData)
{
  dialog <- gtkFileChooserDialog("Select file", NULL, "open",
                                    "gtk-cancel", GtkResponseType["cancel"],
                                    "gtk-open", GtkResponseType["accept"])
  fn <- paste(userData, "file", sep="_")
  #cat(fn)
  entry <- get(fn, interf_s)
  response <- dialog$run()
  if (response == GtkResponseType["accept"]) {
    nomefile <- dialog$getFilename()
    gtkEntrySetText(entry, nomefile)
  }
  dialog$destroy()
}

create_page3_s <- function(assistant)
{
  box <- gtkVBoxNew(FALSE, 0)
  hbox1 <- crea_controllo_s(ACTF, TRUE, TRUE)
  box$packStart(hbox1, FALSE, FALSE, 10)

  notebook1 <- gtkNotebook()
  notebook1$show()
  box$packStart(notebook1, FALSE, FALSE, 0)

  tips_s <- get("tips_s", interf_s)
  def <- c(3, 3, 3, 1, 1, 2)
  sn_radio <- vector()
  for (i in 0:5) {

    table1 <- gtkTableNew(5, 5, FALSE)

    tab1 <- gtkLabel(lbls[[i + 1]][1])
    assign(lbls[[i + 1]][1], table1, interf_s)
    gtkTooltipsSetTip(tips_s, tab1, lbls[[i + 1]][2])
    notebook1$AppendPage(table1, tab1)
    table1$setRowSpacings(5)
    table1$setColSpacings(5)
    #box$packStart(table1, FALSE, FALSE, 0)

    radiobutton1 <- gtkRadioButtonNewWithLabel(NULL, "constant")
    gSignalConnect(radiobutton1, "toggled", toggle_cb_s, c(lbls[[i + 1]][1], 1))
    assign(paste("toggle", lbls[[i + 1]][1], sep="_"), 1, interf_s)
    table1$Attach(radiobutton1, 0, 1, 0, 1, xpadding=5, ypadding=5)
    radiobutton1_group <- gtkRadioButtonGetGroup(radiobutton1)
    #radiobutton1$setGroup(radiobutton1_group)
    if (def[i + 1] == 1) {
      gtkToggleButtonSetActive(radiobutton1, TRUE)
      assign(paste("toggle", lbls[[i + 1]][1], sep="_"), 1, interf_s)
      #cat(sprintf("%s = %s\n", paste("toggle", lbls[[i + 1]][1], sep="_"), 1), "\n")
    }

    radiobutton2 <- gtkRadioButtonNewWithLabel(radiobutton1_group, "uniform")
    gSignalConnect(radiobutton2, "toggled", toggle_cb_s, c(lbls[[i + 1]][1], 2))
    table1$Attach(radiobutton2, 0, 1, 1, 2, xpadding=5, ypadding=5)
#    gtkRadioButtonSetGroup(radiobutton2, radiobutton1_group)
    if (def[i + 1] == 2)
      gtkToggleButtonSetActive(radiobutton2, TRUE)

    checkbutton1 <- gtkCheckButtonNewWithLabel("follow simulatenet outputs")
    sn_radio <- c(sn_radio, checkbutton1)
    assign(paste("toggle1", lbls[[i + 1]][1], sep="_"), 0, interf_s)
    table1$Attach(checkbutton1, 0, 1, 4, 5, xpadding=5, ypadding=5)

    checkbutton2 <- gtkCheckButtonNewWithLabel("file")
    gSignalConnect(checkbutton1, "toggled", toggle1_cb_s, c(lbls[[i + 1]][1], 1, checkbutton2))
    gSignalConnect(checkbutton2, "toggled", toggle1_cb_s, c(lbls[[i + 1]][1], 2, checkbutton1))
    table1$Attach(checkbutton2, 0, 1, 5, 6, xpadding=5, ypadding=5)

    entry1 <- crea_controllo_s(ALPHA + 8 * i, FALSE, TRUE)
    gtkEntrySetWidthChars(entry1, 5)
    table1$Attach(entry1, 1, 2, 0, 1)

    entry2 <- crea_controllo_s(ALPHA + 1 + 8 * i, FALSE, TRUE)
    gtkEntrySetWidthChars(entry2, 5)
    table1$Attach(entry2, 1, 2, 1, 2)

    entry3 <- crea_controllo_s(ALPHA + 2 + 8 * i, FALSE, TRUE)
    gtkEntrySetWidthChars(entry3, 5)
    table1$Attach(entry3, 3, 4, 1, 2, xpadding=5, ypadding=5)

    radiobutton3 <- gtkRadioButtonNewWithLabel(radiobutton1_group, "normal")
    gSignalConnect(radiobutton3, "toggled", toggle_cb_s, c(lbls[[i + 1]][1], 3))
    table1$Attach(radiobutton3, 0, 1, 2, 3, xpadding=5, ypadding=5)
#    radiobutton3$setGroup(radiobutton1_group)
    #radiobutton1_group <- gtk_radio_button_get_group(radiobutton3)
    if (def[i + 1] == 3)
      gtkToggleButtonSetActive(radiobutton3, TRUE)

    radiobutton4 <- gtkRadioButtonNewWithLabel(radiobutton1_group, "log-normal")
    gSignalConnect(radiobutton4, "toggled", toggle_cb_s, c(lbls[[i + 1]][1], 4))
    table1$Attach(radiobutton4, 0, 1, 3, 4, xpadding=5, ypadding=5)
#    radiobutton4$setGroup(radiobutton1_group)
    #radiobutton1_group <- gtk_radio_button_get_group(radiobutton4)

    entry4 <- crea_controllo_s(ALPHA + 3 + 8 * i, FALSE, TRUE)
    gtkEntrySetWidthChars(entry4, 5)
    table1$Attach(entry4, 1, 2, 2, 3)
    entry6 <- crea_controllo_s(ALPHA + 4 + 8 * i, FALSE, TRUE)
    gtkEntrySetWidthChars(entry6, 5)
    table1$Attach(entry6, 3, 4, 2, 3, xpadding=5, ypadding=5)

    entry7 <- crea_controllo_s(ALPHA + 5 + 8 * i, FALSE, TRUE)
    gtkEntrySetWidthChars(entry7, 5)
    table1$Attach(entry7, 1, 2, 3, 4)

    entry8 <- crea_controllo_s(ALPHA + 6 + 8 * i, FALSE, TRUE)
    gtkEntrySetWidthChars(entry8, 5)
    table1$Attach(entry8, 3, 4, 3, 4, xpadding=5, ypadding=5)

    entry9 <- crea_controllo_s(ALPHA + 7 + 8 * i, FALSE, TRUE)
    gtkEntrySetWidthChars(entry9, 5)
    entry9$setEditable(FALSE)
    table1$Attach(entry9, 1, 2, 5, 6)
    assign(paste("file", lbls[[i + 1]][1], sep="_"), entry9, interf_s)

    button1 <- gtkButtonNewWithLabel("...")
    gSignalConnect(button1, "clicked", select_file_s, lbls[[i + 1]][1])
    table1$Attach(button1, 2, 3, 5, 6, 0, 0, 0, 0)
  }
  assign("sn_radio", sn_radio, interf_s)

  assistant$appendPage(box)
  assistant$setPageTitle(box, "Dynamics")
  assistant$setPageComplete(box, TRUE)

  pixbuf <- assistant$renderIcon(GTK_STOCK_DIALOG_INFO, "dialog")
  assistant$setPageHeaderImage(box, pixbuf)
}

spin_changed_s <- function(widget, event, userData)
{
  val <- widget$getValueAsInt()
  ef1 <- get("ef", interf_s)
  if (val <= length(ef1)) {
    gtkEntrySetText(get("ext.fun", interf_s), ef1[[val]][[1]])
    gtkEntrySetText(get("ext.in", interf_s), ef1[[val]][[2]])
  }
  else {
    gtkEntrySetText(get("ext.fun", interf_s), "")
    gtkEntrySetText(get("ext.in", interf_s), "")
  }
}

ef_aggiorna_s <- function(widget, event, userData)
{
  ef1 <- get("ef", interf_s)
  val <- userData$getValueAsInt()
  if (val <= length(ef1)) {
    ef1[[val]][[1]] <- gtkEntryGetText(get("ext.fun", interf_s))
    pars_ctrl_s(widget, event, userData)
  }
  assign("ef", ef1, interf_s)
  return (FALSE)
}

ei_aggiorna_s <- function(widget, event, userData)
{
  ef1 <- get("ef", interf_s)
  val <- userData$getValueAsInt()
  if (val <= length(ef1)) {
    lstval_ctrl_s(widget, event, userData)
    ef1[[val]][[2]] <- gtkEntryGetText(get("ext.in", interf_s))
  }
  assign("ef", ef1, interf_s)
  return (FALSE)
}

create_page4_s <- function(assistant)
{
  box <- gtkVBox(FALSE, 12)
  box$setBorderWidth(12)
  entry1 <- gtkDrawingArea()
  asCairoDevice(entry1)
  entry1$show()
  gtkWidgetSetSizeRequest(entry1, 200, 300)
  #scale_cb <- function(range) { plot(1:20,(1:20)^range$getValue()) }
  box$packStart(entry1, FALSE, FALSE, 0)
  box1 <- gtkHBox(FALSE, 2)
  radiobutton1 <- gtkRadioButtonNew(NULL)
  gSignalConnect(radiobutton1, "toggled", toggle_cb_s, c("ext.fun", 1))
  box1$packStart(radiobutton1, FALSE, FALSE, 0)
  radiobutton1_group <- gtkRadioButtonGetGroup(radiobutton1)
  spn <- gtkSpinButtonNewWithRange(1, 100, 1)
  spn$SetWidthChars(5)
  gSignalConnect(spn, "value-changed", spin_changed_s)
  box1$packStart(spn, FALSE, FALSE, 0)
  entry2 <- crea_controllo_s(EXTF, TRUE, TRUE)
  assign("ef", ef, interf_s)
  gSignalConnect(get("ext.fun", interf_s), "focus-out-event", ef_aggiorna_s, spn)
  box1$packStart(entry2, FALSE, FALSE, 0)
  entry21 <- crea_controllo_s(EXTIN, TRUE, TRUE)
  gSignalConnect(get("ext.in", interf_s), "focus-out-event", ei_aggiorna_s, spn)
  box1$packStart(entry21, FALSE, FALSE, 0)
  entry3 <- gtkButton("Draw")
  gSignalConnect(entry3, "clicked", preview_s, get("ext.fun", interf_s))
  box1$packStart(entry3, FALSE, FALSE, 0)
  box$packStart(box1, FALSE, FALSE, 0)
  radiobutton2 <- gtkRadioButtonNewWithLabel(radiobutton1_group, "none")
  gSignalConnect(radiobutton2, "toggled", toggle_cb_s, c("ext.fun", 2))
  box$packStart(radiobutton2, TRUE, TRUE, 0)
  assign("toggle_ext.fun", 2, interf_s)
  gtkToggleButtonSetActive(radiobutton2, TRUE)
#  gtkRadioButtonSetGroup(radiobutton2, radiobutton1_group)

  assistant$appendPage(box)
  #assistant$setPageType(box, "confirm")
  assistant$setPageComplete(box, TRUE)
  assistant$setPageTitle(box, "External stimuli")

  pixbuf <- assistant$renderIcon(GTK_STOCK_DIALOG_INFO, "dialog")
  assistant$setPageHeaderImage(box, pixbuf)
}

select_file1_s <- function(widget, userData)
{
  dialog <- gtkFileChooserDialog("Select file", NULL, "open",
                                    "gtk-cancel", GtkResponseType["cancel"],
                                    "gtk-open", GtkResponseType["accept"])
  entry <- get("file_times", interf_s)
  response <- dialog$run()
  if (response == GtkResponseType["accept"]) {
    nomefile <- dialog$getFilename()
    gtkEntrySetText(entry, nomefile)
  }
  dialog$destroy()
}

create_page5_s <- function(assistant)
{
  box <- gtkVBox(FALSE, 12)
  box$setBorderWidth(12)

  radiobutton1 <- gtkRadioButtonNewWithLabel(NULL, "none")
  gSignalConnect(radiobutton1, "toggled", toggle_cb_s, c("ko", 1))
  assign(paste("toggle", "ko", sep="_"), 1, interf_s)
  box$packStart(radiobutton1, FALSE, FALSE, 0)
  radiobutton1_group <- gtkRadioButtonGetGroup(radiobutton1)

  radiobutton2 <- gtkRadioButtonNewWithLabel(radiobutton1_group, "set")
  gSignalConnect(radiobutton2, "toggled", toggle_cb_s, c("ko", 2))
  box1 <- gtkHBox(FALSE, 0)
  box1$packStart(radiobutton2, TRUE, TRUE, 0)
  entry1 <- crea_controllo_s(KO, FALSE, TRUE)
  gtkEntrySetWidthChars(entry1, 5)
  box1$packStart(entry1, TRUE, TRUE, 0)
  box$packStart(box1, FALSE, FALSE, 0)

  radiobutton3 <- gtkRadioButtonNewWithLabel(radiobutton1_group, "each in turn")
  gSignalConnect(radiobutton3, "toggled", toggle_cb_s, c("ko", 3))
  box$packStart(radiobutton3, FALSE, FALSE, 0)

  assistant$appendPage(box)
  assistant$setPageTitle(box, "Knock-out")
  assistant$setPageComplete(box, TRUE)

  pixbuf <- assistant$renderIcon(GTK_STOCK_DIALOG_INFO, "dialog")
  assistant$setPageHeaderImage(box, pixbuf)
}

create_page6_s <- function(assistant)
{
  box <- gtkVBoxNew(FALSE, 0)


  table1 <- gtkTableNew(3, 8, FALSE)
  box$packStart(table1, FALSE, FALSE, 10)
  table1$setBorderWidth(5)
  table1$setRowSpacings(5)
  table1$setColSpacings(5)

  box1 <- gtkHBox(FALSE, 12)
  label3 <- gtkLabelNew ("from")
  label3$setAlignment(0, 0.5)
  box1$packStart(label3, FALSE, FALSE, 10)
  entry1 <- crea_controllo_s(TIMES, FALSE, TRUE)
  gtkEntrySetWidthChars(entry1, 5)
  box1$packStart(entry1, FALSE, FALSE, 10)

  label4 <- gtkLabelNew("step")
  label4$setAlignment(0, 0.5)
  box1$packStart(label4, FALSE, FALSE, 10)
  entry2 <- crea_controllo_s(TIMES + 2, FALSE, TRUE)
  gtkEntrySetWidthChars(entry2, 5)
  box1$packStart(entry2, FALSE, FALSE, 10)

  label5 <- gtkLabelNew ("to")
  label5$setAlignment(0, 0.5)
  box1$packStart(label5, FALSE, FALSE, 10)
  entry3 <- crea_controllo_s(TIMES + 1, FALSE, TRUE)
  gtkEntrySetWidthChars(entry3, 5)
  gSignalConnect(entry3, "focus-out-event", gt_ctrl_s, entry1)
  box1$packStart(entry3, FALSE, FALSE, 10)
  table1$Attach(box1, 2, 3, 0, 1)

  label6 <- gtkLabelNew(glbls_s[[METHOD]][[1]])
  label6$setAlignment(0, 0.5)
  table1$Attach(label6, 0, 1, 3, 4)
  comboboxentry1 <- crea_controllo_s(METHOD, FALSE, TRUE)
  table1$Attach(comboboxentry1, 2, 3, 3, 4)

  label6 <- gtkLabelNew(glbls_s[[SIM]][[1]])
  label6$setAlignment(0, 0.5)
  table1$Attach(label6, 0, 1, 4, 5)
  entry <- crea_controllo_s(SIM, FALSE, TRUE)
  table1$Attach(entry, 2, 3, 4, 5)

  label6 <- gtkLabelNew(glbls_s[[SAVE]][[1]])
  label6$setAlignment(0, 0.5)
  table1$Attach(label6, 0, 1, 5, 6)
  chk <- crea_controllo_s(SAVE, FALSE, TRUE)
  table1$Attach(chk, 2, 3, 5, 6)

  button1 <- gtkButtonNewWithLabel("...")
  gSignalConnect(button1, "clicked", select_file1_s)
  table1$Attach(button1, 3, 4, 1, 2)

  entry4 <- crea_controllo_s(TIMES + 3, FALSE, TRUE)
  entry4$setEditable(FALSE)
  gtkEntrySetWidthChars(entry4, 5)
  table1$Attach(entry4, 2, 3, 1, 2)

  radiobutton1 <- gtkRadioButtonNewWithLabel(NULL, "sequence")
  gSignalConnect(radiobutton1, "toggled", toggle_cb_s, c("times", 1))
  assign("toggle_times", 1, interf_s)
  table1$Attach(radiobutton1, 1, 2, 0, 1)
  radiobutton1_group <- gtkRadioButtonGetGroup(radiobutton1)

  radiobutton2 <- gtkRadioButtonNewWithLabel(radiobutton1_group, "file")
  gSignalConnect(radiobutton2, "toggled", toggle_cb_s, c("times", 2))
  table1$Attach(radiobutton2, 1, 2, 1, 2)

  label1 <- gtkLabelNew("times")
  label1$setAlignment(0, 0.5)
  tips_s <- get("tips_s", interf_s)
  gtkTooltipsSetTip(tips_s, label1, "time samples at which explicit estimates of gene expression are desired")
  table1$Attach(label1, 0, 1, 0, 1)

  box2 <- gtkHBox(FALSE, 12)
  label61 <- gtkLabelNew("stationarity")
  label61$setAlignment(0, 0.5)
  table1$Attach(label61, 0, 1, 2, 3)
  label51 <- gtkLabelNew ("stat. threshold")
  label51$setAlignment(0, 0.5)
  box2$packStart(label51, FALSE, FALSE, 10)
  entry31 <- crea_controllo_s(TIMES + 4, FALSE, TRUE)
  gtkEntrySetWidthChars(entry31, 5)
  box2$packStart(entry31, FALSE, FALSE, 10)
  table1$Attach(box2, 2, 3, 2, 3)
  
  label52 <- gtkLabelNew ("stat. width")
  label52$setAlignment(0, 0.5)
  box2$packStart(label52, FALSE, FALSE, 10)
  entry32 <- crea_controllo_s(TIMES + 5, FALSE, TRUE)
  gtkEntrySetWidthChars(entry32, 5)
  box2$packStart(entry32, FALSE, FALSE, 10)
  table1$Attach(box2, 2, 3, 2, 3)

  assistant$appendPage(box)
  assistant$setPageTitle(box, "Solutions and Simulations")
  assistant$setPageComplete(box, TRUE)

  pixbuf <- assistant$renderIcon(GTK_STOCK_DIALOG_INFO, "dialog")
  assistant$setPageHeaderImage(box, pixbuf)
}

create_page7_s <- function(assistant)
{
  box <- gtkVBox(FALSE, 12)
  box$setBorderWidth(12)

  label <- gtkLabel("All parameters have been set: press 'Apply' to start the simulation")
  box$packStart(label, FALSE, FALSE, 10)

  table1 <- gtkTableNew(10, 4, FALSE)
  box$packStart(table1, FALSE, FALSE, 10)
  table1$setBorderWidth(5)
  table1$setRowSpacings(5)
  table1$setColSpacings(5)

  assign("riep", table1, interf_s)

  assistant$appendPage(box)
  assistant$setPageType(box, "confirm")
  assistant$setPageComplete(box, TRUE)
  assistant$setPageTitle(box, "Confirmation")

  pixbuf <- assistant$renderIcon(GTK_STOCK_DIALOG_INFO, "dialog")
  assistant$setPageHeaderImage(box, pixbuf)
}

onExit <- function(widget, event, userData)
{
  gtkWidgetDestroy(widget)
  setwd("..")
  return (FALSE)
}

netsim_simulate <- function(vista=FALSE)
{

  if (!file.exists("netsim"))
    stop("To simulate the profiles you must first run 'netsim_generate()' and save at least a network!")
  setwd("netsim")

  mx <- length(i <- grep("weights", dir()))
  pad_w <- ceiling(log10(mx + 1))
  assign("pad_w", pad_w, interf_s)
    
  
  assign("vista_s", vista, interf_s)
  assistant <- gtkAssistant(show=F)
  gSignalConnect(assistant, "destroy-event", onExit)
  assign("main", assistant, interf_s)

  assistant$setDefaultSize(-1, 300)

  assign("primo_s", TRUE, interf_s)

  # creo un oggetto per i Tooltips
  tips_s <- gtkTooltipsNew()
  assign("tips_s", tips_s, interf_s)
  assign("argom_s", vector("list", TOT_ARG_s), interf_s)

  create_page1_s(assistant)
  create_page2_s(assistant)
  create_page3_s(assistant)
  create_page4_s(assistant)
  create_page5_s(assistant)
  create_page6_s(assistant)
  create_page7_s(assistant)

  gSignalConnect(assistant, "cancel", onExit)
  gSignalConnect(assistant, "close", onExit)
  gSignalConnect(assistant, "apply", on_assistant_apply_s)
  gSignalConnect(assistant, "prepare", on_assistant_prepare_s)

  assistant$showAll()
}

