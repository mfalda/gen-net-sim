library(RGtk2)
library(cairoDevice)

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
                 as.character(nome),
                 as.vector(args)
					)
    return(vett)
}

interf_g <- new.env()

# costanti simboliche
NGENI = 1
CONN = 2
MAXR = 3
GAMMA = 4
IND = 5
CF = 6
NUMSUB = 7
KAPPA = 8
FPR_G = 9
ACTF = 10
WP = 11
PARAM = 67
ALPHA = 12
TIMES = 60
METHOD = 64
SIM = 65
SAVE = 66

TOT_ARG_g = 22


num_ctrl_g <- function(widget, event, userData)
{
  num <- as.double(gtkEntryGetText(widget))
  if (is.na(num)) {
    message_dialog_g("the value must be a number: set to 0")
    widget$setText(0.0)
    widget$grabFocus()
  }
  return (FALSE)
}

pars_ctrl_g <- function(widget, event, userData)
{
  txt <- gtkEntryGetText(widget)
  ris <- try(calc_dll(txt, c(0)))
  if (inherits(ris, "try-error")) {
    message_dialog_g(geterrmessage())
    widget$setText("0.0")
  }
  return (FALSE)
}

gz_ctrl_g <- function(widget, event, userData)
{
  num <- as.double(gtkEntryGetText(widget))
  if (is.na(num)) {
    message_dialog_g("the value must be specified: set to 5")
    widget$setText(5)
  }
  else if (num < 0) {
    message_dialog_g("the value must be greater or equal than zero")
    widget$setText(abs(num))
    widget$grabFocus()
  }
  return (FALSE)
}

gt_ctrl_g <- function(widget, event, userData)
{
  rif <- as.double(gtkEntryGetText(userData))
  if (as.double(gtkEntryGetText(widget)) <= rif) {
    message_dialog_g(paste("the value must be greater than", rif))
    widget$setText(rif + 1)
    widget$grabFocus()
  }
  return (FALSE)
}

message_dialog_g <- function(testo)
{
  w <- get("main", interf_g)
  dialog <- gtkMessageDialogNew(w, c("modal", "destroy-with-parent"), "error", "ok", testo)
  dialog$run()
  dialog$destroy()
}

# valutazione di espressioni
run_g <- function() {
  code <- gtkEntryGetText(get("f.pr.and", interf_g))
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
    message_dialog_g(geterrmessage())
    return()
  }
  return (ris)
}

# funzione di call-back per il pulsante "Funzione"
preview_g <- function(widget, userData)
{
  x <- (0 : 100) / 100
  y <- run_g()
  plot(x, y, type="l")
}

act.fun_cb_g <- function(widget, event, userData)
{
  alpha <- get("alpha", interf_g)
  beta <- get("beta", interf_g)
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

toggle_cb_g <- function(widget, userData)
{
  if (gtkToggleButtonGetActive(widget)) {
    assign(paste("toggle", userData[1], sep="_"), userData[2], interf_g)
    #cat(sprintf("%s = %s\n", paste("toggle", userData[1], sep="_"), userData[2]), "\n")
  }
}

# variabili globali
glbls_g <- list(
  # nome, tipo, default, on_exit, descr
  list("N", TESTO, 5, gz_ctrl_g, "number of genes in the network"),
  # nome, tipo, default, valori, call_back, descr
  list("connectivity", LISTA, 3, c("random", "scale free", "MTM", "geometric"), NULL, "type of connectivity in the network"),
  list("max.reg", TESTO, 12, gz_ctrl_g, "maximum number of regulators that each node (gene) in the network can have"),
  list("gamma", TESTO, 2.2, num_ctrl_g, "the parameter of the power law distribution of the degree of connectivity of the nodes in the graph"),
  list("INdegree", LISTA, 1, c("free", "out"), NULL, "'free': the in-degree distribution is not constrained to follow any distribution
  'out': it is constrained to follow the same distribution of the out-degree"),
  list("Cf.cl", TESTO, 0.4, num_ctrl_g, "average clustering coefficient of each sub-network in the graph"),
  # nome, tipo, defaults
  list("num.subnet", VETT, list(5, 5, 10), c("the maximum number of nodes in motif of type 1",
  "the maximum number of nodes in motif of type 2", "the maximum number of nodes in motif of type 3")),
  list("kappa", TESTO, 3, num_ctrl_g, "average number of regulators that each node (gene) in the network has with random topology"),
  list("f.pr.and", TESTO, "0.5", pars_ctrl_g, "function with domain and codomain in [0-1] that expresses the probability to obtain a cooperative rather than a synergic rule"),
  list("act.fun", LISTA, 2, c("linear", "sigmoidal"), act.fun_cb_g, "the activation function"),
  list("weight.par", VETT, list(1, 0), c("mean of the Gaussian distribution used to sample regulatory efficiency",
  "sd of the Gaussian distribution used to sample regulatory efficiency")),
  # 12
  list("alpha_const", TESTO, 1, gz_ctrl_g, "constant value"), # vector of parameters of the Activation sigmoid function
  list("alpha_a", TESTO, 0, num_ctrl_g, "parameter 'a' for the uniform distribution"),
  list("alpha_b", TESTO, 1, num_ctrl_g, "parameter 'b' for the uniform distribution"),
  list("alpha_nm", TESTO, 10, gz_ctrl_g, "mean for the normal distribution"),
  list("alpha_nsd", TESTO, 0.2, gz_ctrl_g, "std. dev. for the normal distribution"),
  list("alpha_lm", TESTO, 1, gz_ctrl_g, "mean for the log-normal distribution"),
  list("alpha_lsd", TESTO, 0, gz_ctrl_g, "std. dev. for the log-normal distribution"),
  list("alpha_file", TESTO, NULL, NULL, "file from which to read the desired values"),
  # 20
  list("beta_const", TESTO, 1, gz_ctrl_g, "constant value"),
  list("beta_a", TESTO, 0, num_ctrl_g, "parameter 'a' for the uniform distribution"),
  list("beta_b", TESTO, 1, num_ctrl_g, "parameter 'b' for the uniform distribution"),
  list("beta_nm", TESTO, 0.5, gz_ctrl_g, "mean for the normal distribution"),
  list("beta_nsd", TESTO, 0.01, gz_ctrl_g, "std. dev. for the normal distribution"),
  list("beta_lm", TESTO, 1, gz_ctrl_g, "mean for the log-normal distribution"),
  list("beta_lsd", TESTO, 0, gz_ctrl_g, "std. dev. for the log-normal distribution"),
  list("beta_file", TESTO, NULL, NULL, "file from which to read the desired values"),
  # 28
  list("lambda_const", TESTO, 1, gz_ctrl_g, "constant value"),
  list("lambda_a", TESTO, 0, num_ctrl_g, "parameter 'a' for the uniform distribution"),
  list("lambda_b", TESTO, 1, num_ctrl_g, "parameter 'b' for the uniform distribution"),
  list("lambda_nm", TESTO, 1, gz_ctrl_g, "mean for the normal distribution"),
  list("lambda_nsd", TESTO, 0.1, gz_ctrl_g, "std. dev. for the normal distribution"),
  list("lambda_lm", TESTO, 1, gz_ctrl_g, "mean for the log-normal distribution"),
  list("lambda_lsd", TESTO, 0, gz_ctrl_g, "std. dev. for the log-normal distribution"),
  list("lambda_file", TESTO, NULL, NULL, "file from which to read the desired values"),
  # 36
  list("Xmin_const", TESTO, 0, gz_ctrl_g, "constant value"),
  list("Xmin_a", TESTO, 0, num_ctrl_g, "parameter 'a' for the uniform distribution"),
  list("Xmin_b", TESTO, 1, num_ctrl_g, "parameter 'b' for the uniform distribution"),
  list("Xmin_nm", TESTO, 1, gz_ctrl_g, "mean for the normal distribution"),
  list("Xmin_nsd", TESTO, 0, gz_ctrl_g, "std. dev. for the normal distribution"),
  list("Xmin_lm", TESTO, 1, gz_ctrl_g, "mean for the log-normal distribution"),
  list("Xmin_lsd", TESTO, 0, gz_ctrl_g, "std. dev. for the log-normal distribution"),
  list("Xmin_file", TESTO, NULL, NULL, "file from which to read the desired values"),
  # 44
  list("Xmax_const", TESTO, 10, gz_ctrl_g, "constant value"),
  list("Xmax_a", TESTO, 0, num_ctrl_g, "parameter 'a' for the uniform distribution"),
  list("Xmax_b", TESTO, 1, num_ctrl_g, "parameter 'b' for the uniform distribution"),
  list("Xmax_nm", TESTO, 1, gz_ctrl_g, "mean for the normal distribution"),
  list("Xmax_nsd", TESTO, 0, gz_ctrl_g, "std. dev. for the normal distribution"),
  list("Xmax_lm", TESTO, 1, gz_ctrl_g, "mean for the log-normal distribution"),
  list("Xmax_lsd", TESTO, 0, gz_ctrl_g, "std. dev. for the log-normal distribution"),
  list("Xmax_file", TESTO, NULL, NULL, "file from which to read the desired values"),
  # 52
  list("X0_const", TESTO, 1, gz_ctrl_g, "constant value"),
  list("X0_a", TESTO, 0, num_ctrl_g, "parameter 'a' for the uniform distribution"),
  list("X0_b", TESTO, 1, num_ctrl_g, "parameter 'b' for the uniform distribution"),
  list("X0_nm", TESTO, 1, gz_ctrl_g, "mean for the normal distribution"),
  list("X0_nsd", TESTO, 0, gz_ctrl_g, "std. dev. for the normal distribution"),
  list("X0_lm", TESTO, 1, gz_ctrl_g, "mean for the log-normal distribution"),
  list("X0_lsd", TESTO, 0, gz_ctrl_g, "std. dev. for the log-normal distribution"),
  list("X0_file", TESTO, NULL, NULL, "file from which to read the desired values"),
  # 60
  list("from", TESTO, "0", num_ctrl_g, "from"),
  list("to", TESTO, "5", num_ctrl_g, "to"),
  list("step", TESTO, "0.05", num_ctrl_g, "step"),
  list("file_times", TESTO, "", NULL, "file from which to read the desired values"),
  list("method", LISTA, 2, c("Euler",
  "Embedded Runge-Kutta-Fehlberg (4,5)",
  "Embedded Runge-Kutta Cash-Karp (4,5)",
  "Embedded Runge-Kutta Prince-Dormand (8,9)",
  "Implicit 2nd order Runge-Kutta at Gaussian points",
  "Implicit 4th order Runge-Kutta at Gaussian points",
  "Implicit Gear method (M = 1)",
  "Implicit Gear method (M = 2)"),
  NULL, "Method used to solve differential equations"),
  # 65
  list("simulations", TESTO, "3", gz_ctrl_g, "number of networks to be created"),
  list("save", CHECKBOX, TRUE, NULL, "save simulations on disk"),
  # 67
  list("param", VETT, NULL, "non visuale")
)

titoli_g <- c("Topology", "Topology", "Rules", "Dynamics", "Solutions and Simulations", "Done")

riepilogo_g <- function(table1, num, num1, x, y)
{
  primo_g <- get("primo_g", interf_g)
  argom_g <- get("argom_g", interf_g)
  if (num != PARAM)
    ctrl <- get(glbls_g[[num]][[1]], interf_g)
  txt_lbl <- glbls_g[[num]][[1]]
  if (num >= ALPHA && num < TIMES) {
    l <- (num - 4) / 8
    txt_lbl <- lbls[[l]][1]
    attv <- get(paste("toggle", lbls[[l]][1], sep="_"), interf_g)
    N <- gtkEntryGetText(get(glbls_g[[NGENI]][[1]], interf_g))
    volte <- as.integer(gtkEntryGetText(get(glbls_g[[SIM]][[1]], interf_g)))
    if (attv == 1)
      txt1 <- sprintf("rep(%s, %s)", gtkEntryGetText(get(glbls_g[[num]][[1]], interf_g)), N)
    else if (attv == 2) {
      a <- gtkEntryGetText(get(glbls_g[[num + 1]][[1]], interf_g))
      b <- gtkEntryGetText(get(glbls_g[[num + 2]][[1]], interf_g))
      txt1 <- sprintf("c(%s, %s)", a, b)
    }
    else if (attv == 3) {
      m <- gtkEntryGetText(get(glbls_g[[num + 3]][[1]], interf_g))
      s <- gtkEntryGetText(get(glbls_g[[num + 4]][[1]], interf_g))
      txt1 <- sprintf("c(%s, %s)", m, s)
    }
    else if (attv == 4) {
      m <- gtkEntryGetText(get(glbls_g[[num + 5]][[1]], interf_g))
      s <- gtkEntryGetText(get(glbls_g[[num + 6]][[1]], interf_g))
      txt1 <- sprintf("c(%s, %s)", m, s)
    }
    else if (attv == 5) {
      str <- gsub("\\", "/", gtkEntryGetText(get(glbls_g[[num + 7]][[1]], interf_g)), fixed=TRUE)
      txt1 <- sprintf("read.table('%s')", str)
    }
    else
      stop(sprintf("valore di 'attv' non valido (%d)!", attv))
    if (primo_g)
      argom_g[[num1]] <- list("espressione", txt1, NULL)
    else
      argom_g[[num1]][[2]] <- txt1
  }
  else if (num == NUMSUB) {
    if (gtkWidgetIsSensitive(ctrl[[1]]))
      txt1 <- sprintf("c(%3.3f, %3.3f, %3.3f)", as.double(gtkEntryGetText(ctrl[[1]])),
        as.double(gtkEntryGetText(ctrl[[2]])), as.double(gtkEntryGetText(ctrl[[3]])))
    else
      txt1 <- "c(0, 0, 0)"
    if (primo_g)
      argom_g[[num1]] <- list("espressione", txt1, NULL)
    else
      argom_g[[num1]][[2]] <- txt1
  }
  else if (num == WP && gtkWidgetIsSensitive(ctrl[[1]])) {
    txt1 <- sprintf("c(%3.3f, %3.3f)", as.double(gtkEntryGetText(ctrl[[1]])),
      as.double(gtkEntryGetText(ctrl[[2]])))
    if (primo_g)
      argom_g[[num1]] <- list("espressione", txt1, NULL)
    else
      argom_g[[num1]][[2]] <- txt1
  }
  else if (num == FPR_G) {
    txt1 <- gtkEntryGetText(ctrl)
    if (primo_g)
      argom_g[[num1]] <- list("stringa", txt1, NULL)
    else
      argom_g[[num1]][[2]] <- txt1
  }
  else if (num >= TIMES && num < METHOD) {
    txt_lbl <- "times"
    txt <- get("toggle_times", interf_g)
    if (txt == 1) {
      da <- gtkEntryGetText(get(glbls_g[[TIMES]][[1]], interf_g))
      a <- gtkEntryGetText(get(glbls_g[[TIMES + 1]][[1]], interf_g))
      passo <- gtkEntryGetText(get(glbls_g[[TIMES + 2]][[1]], interf_g))
      txt1 <- sprintf("seq(%s, %s, %s)", da, a, passo)
    }
    else {
      str <-gsub("\\", "/", gtkEntryGetText(get(glbls_g[[TIMES + 3]][[1]], interf_g)), fixed=TRUE)
      txt1 <- sprintf("read.table('%s')", str)
    }
    if (primo_g)
      argom_g[[num1]] <- list("espressione", txt1, NULL)
    else
      argom_g[[num1]][[2]] <- txt1
  }
  else if (num == PARAM) {
    txt <- get(paste("toggle", lbls[[3]][1], sep="_"), interf_g)
    n <- as.integer(txt) - 1
    if (n == 4)
      n <- 0
    txt1 <- paste("c", n, sep="(")
    for (i in c(1,2,4,5,6)) {
      txt <- get(paste("toggle", lbls[[i]][1], sep="_"), interf_g)
      n <- as.integer(txt) - 1
      if (n == 4)
        n <- 0
      txt1 <- paste(txt1, n, sep=",")
    }
    txt1 <- paste(txt1, ")", sep="")
    if (primo_g)
      argom_g[[num1]] <- list("espressione", txt1, NULL)
    else
      argom_g[[num1]][[2]] <- txt1
  }
  else if (gtkWidgetIsSensitive(ctrl)) {
    if (glbls_g[[num]][[2]] == 1) {
      txt1 <- gtkEntryGetText(ctrl)
      if (primo_g)
        argom_g[[num1]] <- list("numero", sprintf('%3.3f', as.double(txt1)), NULL)
      else
        argom_g[[num1]][[2]] <- sprintf('%3.3f', as.double(txt1))
    }
    else if (glbls_g[[num]][[2]] == LISTA) {
      if (num == METHOD) {
        sigle <- c("Euler", "rkf45", "rkck", "rk8pd", "rk2imp", "rk4imp", "gear1", "gear2")
        txt1 <- sigle[gtkComboBoxGetActive(ctrl) + 1]
      }
      else
        txt1 <- gtkComboBoxGetActiveText(ctrl)
      if (primo_g)
        argom_g[[num1]] <- list("stringa", as.character(txt1), NULL)
      else
        argom_g[[num1]][[2]] <- as.character(txt1)
    }
    else if (glbls_g[[num]][[2]] == CHECKBOX) {
      txt1 <- gtkToggleButtonGetActive(ctrl)
      if (primo_g)
        argom_g[[num1]] <- list("booleano", as.logical(txt1), NULL)
      else
        argom_g[[num1]][[2]] <- as.logical(txt1)
    }
    else {
      txt1 <- "-"
      if (primo_g)
        argom_g[[num1]] <- list("saltato", glbls_g[[num]][[3]], NULL)
    }
  }
  else {
    txt1 <- "-"
    if (primo_g)
      argom_g[[num1]] <- list("saltato", glbls_g[[num]][[3]], NULL)
  }
  if (primo_g) {
    label1 <- gtkLabelNew()
    label1$setAlignment(0, 0.5)
    gtkLabelSetText(label1, txt_lbl)
    table1$Attach(label1, 4 * x, 4 * x + 1, y, y + 1)
    label2 <- gtkLabelNew()
    label2$setAlignment(1, 1)
    argom_g[[num1]][[3]] <- label2
    table1$Attach(label2, 4 * x + 2, 4 * x + 3, y, y + 1)
  }
  else {
    label2 <- argom_g[[num1]][[3]]
  }
  gtkLabelSetText(label2, txt1)
  assign("argom_g", argom_g, interf_g)
  return
}

on_assistant_apply_g <- function(widget, data)
{
  a <- list()
  argom_g <- get("argom_g", interf_g)
  #argom[[16]][[2]] <- sprintf("%s / %s", argom_g[[16]][[2]], argom_g[[15]][[2]])
  for (l in 1:TOT_ARG_g) {
    if (argom_g[[l]][[1]] == "espressione" || argom_g[[l]][[1]] == "numero") {
      if (!is.null(argom_g[[l]][[2]])) # altrimenti compare il prompt, NON NULL!
       expr <- try(parse(text=argom_g[[l]][[2]]))
      else
        expr <- NULL
       a[[l]] <- eval(expr)
    }
    #else if (argom_g[[l]][[1]] == "funzione") {
    #   a[[l]] <- eval(try(parse(text=paste("function(x) { return(", argom_g[[l]][[2]], ")}"))))
    #}
    else
      a[[l]] <- argom_g[[l]][[2]]
  }
  if (is.list(a[[10]]))
    a[[10]] <- a[[10]][[1]]
  if (is.list(a[[11]]))
    a[[11]] <- a[[11]][[1]]
  if (is.list(a[[12]]))
    a[[12]] <- a[[12]][[1]]
  if (is.list(a[[13]]))
    a[[13]] <- a[[13]][[1]]
  if (is.list(a[[16]]))
    a[[16]] <- a[[16]][[1]]
  if (is.list(a[[17]]))
    a[[17]] <- a[[17]][[1]]
  if (is.list(a[[18]]))
    a[[18]] <- a[[18]][[1]]
  salvati <- as.logical(a[[TOT_ARG_g - 1]])
  cat("\nArgs: \n", paste(a, "\n"), "\n")
  d <- dir()
  file.remove(d[grep(".*\\.txt", d, perl=TRUE)])
  assign("res_netsim_generated", do.call("simulatenet", a), globalenv())
  if (salvati)
    cat("\nResults have been saved in variable 'res_netsim_generated' and in folder '", getwd(), "'\n")
  else
    cat("\nResults have been saved in variable 'res_netsim_generated'\n")
}

on_assistant_prepare_g <- function(widget, page, data)
{
  current_page <- widget$getCurrentPage()
  n_pages <- widget$getNPages()
  if (current_page == 1) {
    indeg <- get("INdegree", interf_g)
    gamma <- get("gamma", interf_g)
    k <- get("kappa", interf_g)
    txt <- gtkComboBoxGetActiveText(get("connectivity", interf_g))
    cl <- get("Cf.cl", interf_g)
    for (c in get("num.subnet", interf_g))
      gtkWidgetSetSensitive(c, FALSE)
    gtkWidgetSetSensitive(indeg, FALSE)
    gtkWidgetSetSensitive(gamma, FALSE)
    gtkWidgetSetSensitive(k, FALSE)
    gtkWidgetSetSensitive(cl, FALSE)
    if (txt == "MTM") {
      gtkWidgetSetSensitive(indeg, TRUE)
      gtkWidgetSetSensitive(gamma, TRUE)
      gtkWidgetSetSensitive(cl, TRUE)
       for (c in get("num.subnet", interf_g))
        gtkWidgetSetSensitive(c, TRUE)
    }
    else if (txt == "scale free")
      gtkWidgetSetSensitive(gamma, TRUE)
    else if (txt == "random")
      gtkWidgetSetSensitive(k, TRUE)
    else if (txt == "geometric")
      gtkWidgetSetSensitive(k, TRUE)
  }
  else if (current_page == 5) {
    k <- 1
    arg <- c(NGENI, CONN, MAXR, GAMMA, IND, CF, NUMSUB, KAPPA, FPR_G, ALPHA + 3 * 8, ALPHA + 4 * 8, ALPHA + 2 * 8, ALPHA + 5 * 8, WP, ACTF, ALPHA, ALPHA + 8, TIMES, METHOD, PARAM, SAVE, SIM)
    table1 <- get("riep", interf_g)
    for (j in 0:1) {
      for (i in 0:10) {
        if (k <= length(arg)) {
          riepilogo_g(table1, arg[k], k, j, i)
          k <- k + 1
        }
      }
    }
    assign("primo_g", FALSE, interf_g)
  }
  title <- sprintf("Netsim assistant - %s (%d of %d)", titoli_g[current_page + 1], current_page + 1, n_pages)
  widget$setTitle(title)
}

crea_testo_g <- function(l, tips_g)
{
  ctrl <- gtkEntryNew()
  gtkEntrySetAlignment(ctrl, 1)
  if (!is.null(l[[3]]))
	 gtkEntrySetText(ctrl, l[[3]])
  if (!is.null(l[[4]]))
    gSignalConnect(ctrl, "focus-out-event", l[[4]])
  gtkTooltipsSetTip(tips_g, ctrl, l[[5]])
  assign(l[[1]], ctrl, interf_g)
}

crea_chk_g <- function(l, tips_g)
{
  ctrl <- gtkCheckButtonNew()
  assign(l[[1]], ctrl, interf_g)
  gtkToggleButtonSetActive(ctrl, l[[3]])
  gtkTooltipsSetTip(tips_g, ctrl, l[[4]])
  return (ctrl)
}

crea_vett_g <- function(l, tips_g)
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
    gSignalConnect(ctrl, "focus-out-event", num_ctrl_g)
    gtkTooltipsSetTip(tips_g, ctrl, l[[4]][i])
    algn <- gtkAlignmentNew(0, 0, 0, 0)
    gtkAlignmentSetPadding(algn, 5, 5, 5, 5)
    algn$add(ctrl)
    frm1$packStart(algn, FALSE, TRUE, 0)
    cz <- c(cz, ctrl)
    i <- i + 1
  }
  lbl$add(frm1)
  assign(l[[1]], cz, interf_g)
  return (frm)
}

crea_lista_g <- function(l, tips_g)
{
  # lista combinata
  ctrl <- gtkComboBoxNewText()
  for (t in l[[4]])
    gtkComboBoxAppendText(ctrl, t)
  gtkComboBoxSetActive(ctrl, l[[3]] - 1)
  if (!is.null(l[[5]])) {
     gSignalConnect(ctrl, "changed", l[[5]])
  }
  gtkTooltipsSetTip(tips_g, ctrl, l[[6]])
  assign(l[[1]], ctrl, interf_g)
  return (ctrl)
}

crea_controllo_g <- function(num, etichetta, attivo)
{
  tips_g <- get("tips_g", interf_g)
  l <- glbls_g[[num]]
  if (etichetta) {
    frm <- gtkHBoxNew(FALSE, 0)
    # etichetta
    lbl <- gtkLabelNew(l[[1]])
    algn1 <- gtkAlignmentNew(0, 0, 0, 0)
    gtkAlignmentSetPadding(algn1, 5, 5, 5, 5)
    algn1$add(lbl)
  }
  if (l[[2]] == TESTO) {
    ctrl <- crea_testo_g(l, tips_g)
  }
  else if (l[[2]] == LISTA) {
    ctrl <- crea_lista_g(l, tips_g)
  }
  else if (l[[2]] == VETT) {
    ctrl <- crea_vett_g(l, tips_g)
  }
  else if (l[[2]] == CHECKBOX) {
    ctrl <- crea_chk_g(l, tips_g)
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

create_page1_g <- function(assistant)
{
  box <- gtkVBox(FALSE, 12)
  box$setBorderWidth(12)

  entry1 <- crea_controllo_g(NGENI, TRUE, TRUE)
  box$packStart(entry1, FALSE, FALSE, 0)
  entry2 <- crea_controllo_g(CONN, TRUE, TRUE)
  box$packStart(entry2, FALSE, FALSE, 0)
  entry3 <- crea_controllo_g(MAXR, TRUE, TRUE)
  box$packStart(entry3, FALSE, FALSE, 0)

  assistant$appendPage(box)
  assistant$setPageTitle(box, "Topology: common parameters")
  #assistant$setPageType(box, "intro")
  assistant$setPageComplete(box, TRUE)

  pixbuf <- assistant$renderIcon(GTK_STOCK_DIALOG_INFO, "dialog")
  assistant$setPageHeaderImage(box, pixbuf)
}

create_page2_g <- function(assistant)
{
  box <- gtkVBox(FALSE, 12)
  box$setBorderWidth(12)

  entry1 <- crea_controllo_g(GAMMA, TRUE, TRUE)
  box$packStart(entry1, FALSE, FALSE, 0)
  entry2 <- crea_controllo_g(IND, TRUE, TRUE)
  box$packStart(entry2, FALSE, FALSE, 0)
  entry3 <- crea_controllo_g(CF, TRUE, TRUE)
  box$packStart(entry3, FALSE, FALSE, 0)
  entry3 <- crea_controllo_g(NUMSUB, FALSE, TRUE)
  box$packStart(entry3, FALSE, FALSE, 0)
  entry3 <- crea_controllo_g(KAPPA, TRUE, FALSE)
  box$packStart(entry3, FALSE, FALSE, 0)

  assistant$appendPage(box)
  #assistant$setPageType(box, "content")
  assistant$setPageComplete(box, TRUE)
  assistant$setPageTitle(box, "Topology: specific parameters")

  pixbuf <- assistant$renderIcon(GTK_STOCK_DIALOG_INFO, "dialog")
  assistant$setPageHeaderImage(box, pixbuf)
}

create_page3_g <- function(assistant)
{
  box <- gtkVBox(FALSE, 12)
  box$setBorderWidth(12)

  entry1 <- gtkDrawingArea()
  asCairoDevice(entry1)
  gtkWidgetSetSizeRequest(entry1, 200, 300)
  #scale_cb <- function(range) { plot(1:20,(1:20)^range$getValue()) }
  box$packStart(entry1, FALSE, FALSE, 0)
  box1 <- gtkHBox(FALSE, 2)
  entry2 <- crea_controllo_g(FPR_G, FALSE, TRUE)
  box1$packStart(entry2, FALSE, FALSE, 0)
  entry3 <- gtkButton("Draw")
  gSignalConnect(entry3, "clicked", preview_g)
  box1$packStart(entry3, FALSE, FALSE, 0)
  box$packStart(box1, FALSE, FALSE, 0)

  assistant$appendPage(box)
  #assistant$setPageType(box, "confirm")
  assistant$setPageComplete(box, TRUE)
  assistant$setPageTitle(box, "Rules")

  pixbuf <- assistant$renderIcon(GTK_STOCK_DIALOG_INFO, "dialog")
  assistant$setPageHeaderImage(box, pixbuf)
}

select_file_g <- function(widget, userData)
{
  dialog <- gtkFileChooserDialog("Select file", NULL, "open",
                                    "gtk-cancel", GtkResponseType["cancel"],
                                    "gtk-open", GtkResponseType["accept"])
  fn <- paste(userData, "file", sep="_")
  #cat(fn)
  entry <- get(fn, interf_g)
  response <- dialog$run()
  if (response == GtkResponseType["accept"]) {
    nomefile <- dialog$getFilename()
    gtkEntrySetText(entry, nomefile)
  }
  dialog$destroy()
}

create_page4_g <- function(assistant)
{
  box <- gtkVBoxNew(FALSE, 0)
  hbox1 <- crea_controllo_g(ACTF, TRUE, TRUE)
  box$packStart(hbox1, FALSE, FALSE, 10)

  notebook1 <- gtkNotebook()
  notebook1$show()
  box$packStart(notebook1, FALSE, FALSE, 0)

  tips_g <- get("tips_g", interf_g)
  def <- c(3, 3, 3, 1, 1, 2)
  for (i in 0:5) {

    table1 <- gtkTableNew(5, 5, FALSE)

    tab1 <- gtkLabel(lbls[[i + 1]][1])
    assign(lbls[[i + 1]][1], table1, interf_g)
    gtkTooltipsSetTip(tips_g, tab1, lbls[[i + 1]][2])
    notebook1$AppendPage(table1, tab1)
    table1$setRowSpacings(5)
    table1$setColSpacings(5)

    radiobutton1 <- gtkRadioButtonNewWithLabel(NULL, "constant")
    gSignalConnect(radiobutton1, "toggled", toggle_cb_g, c(lbls[[i + 1]][1], 1))
    assign(paste("toggle", lbls[[i + 1]][1], sep="_"), 1, interf_g)
    table1$Attach(radiobutton1, 0, 1, 0, 1, xpadding=5, ypadding=5)
    radiobutton1_group <- gtkRadioButtonGetGroup(radiobutton1)
    #radiobutton1$setGroup(radiobutton1_group)
    if (def[i + 1] == 1) {
      gtkToggleButtonSetActive(radiobutton1, TRUE)
      assign(paste("toggle", lbls[[i + 1]][1], sep="_"), 1, interf_g)
      #cat(sprintf("%s = %s\n", paste("toggle", lbls[[i + 1]][1], sep="_"), 1), "\n")
    }

    radiobutton2 <- gtkRadioButtonNewWithLabel(radiobutton1_group, "uniform")
    gSignalConnect(radiobutton2, "toggled", toggle_cb_g, c(lbls[[i + 1]][1], 2))
    table1$Attach(radiobutton2, 0, 1, 1, 2, xpadding=5, ypadding=5)
    #gtkRadioButtonSetGroup(radiobutton2, radiobutton1_group)
    if (def[i + 1] == 2)
      gtkToggleButtonSetActive(radiobutton2, TRUE)

    radiobutton5 <- gtkRadioButtonNewWithLabel(radiobutton1_group, "file")
    gSignalConnect(radiobutton5, "toggled", toggle_cb_g, c(lbls[[i + 1]][1], 5))
    table1$Attach(radiobutton5, 0, 1, 4, 5, xpadding=5, ypadding=5)
    #radiobutton5$setGroup(radiobutton1_group)
    #radiobutton1_group <- gtk_radio_button_get_group(radiobutton5)

    entry1 <- crea_controllo_g(ALPHA + 8 * i, FALSE, TRUE)
    gtkEntrySetWidthChars(entry1, 5)
    table1$Attach(entry1, 1, 2, 0, 1)

    entry2 <- crea_controllo_g(ALPHA + 1 + 8 * i, FALSE, TRUE)
    gtkEntrySetWidthChars(entry2, 5)
    table1$Attach(entry2, 1, 2, 1, 2)

    entry3 <- crea_controllo_g(ALPHA + 2 + 8 * i, FALSE, TRUE)
    gtkEntrySetWidthChars(entry3, 5)
    table1$Attach(entry3, 3, 4, 1, 2, xpadding=5, ypadding=5)

    if (i != 5) {

      radiobutton3 <- gtkRadioButtonNewWithLabel(radiobutton1_group, "normal")
      gSignalConnect(radiobutton3, "toggled", toggle_cb_g, c(lbls[[i + 1]][1], 3))
      table1$Attach(radiobutton3, 0, 1, 2, 3, xpadding=5, ypadding=5)
      #radiobutton3$setGroup(radiobutton1_group)
      #radiobutton1_group <- gtk_radio_button_get_group(radiobutton3)
      if (def[i + 1] == 3)
        gtkToggleButtonSetActive(radiobutton3, TRUE)

      radiobutton4 <- gtkRadioButtonNewWithLabel(radiobutton1_group, "log-normal")
      gSignalConnect(radiobutton4, "toggled", toggle_cb_g, c(lbls[[i + 1]][1], 4))
      table1$Attach(radiobutton4, 0, 1, 3, 4, xpadding=5, ypadding=5)
      #radiobutton4$setGroup(radiobutton1_group)
      #radiobutton1_group <- gtk_radio_button_get_group(radiobutton4)

      entry4 <- crea_controllo_g(ALPHA + 3 + 8 * i, FALSE, TRUE)
      gtkEntrySetWidthChars(entry4, 5)
      table1$Attach(entry4, 1, 2, 2, 3)
      entry6 <- crea_controllo_g(ALPHA + 4 + 8 * i, FALSE, TRUE)
      gtkEntrySetWidthChars(entry6, 5)
      table1$Attach(entry6, 3, 4, 2, 3, xpadding=5, ypadding=5)

      entry7 <- crea_controllo_g(ALPHA + 5 + 8 * i, FALSE, TRUE)
      gtkEntrySetWidthChars(entry7, 5)
      table1$Attach(entry7, 1, 2, 3, 4)

      entry8 <- crea_controllo_g(ALPHA + 6 + 8 * i, FALSE, TRUE)
      gtkEntrySetWidthChars(entry8, 5)
      table1$Attach(entry8, 3, 4, 3, 4, xpadding=5, ypadding=5)

     }

    entry9 <- crea_controllo_g(ALPHA + 7 + 8 * i, FALSE, TRUE)
    gtkEntrySetWidthChars(entry9, 5)
    entry9$setEditable(FALSE)
    table1$Attach(entry9, 1, 2, 4, 5)
    assign(paste("file", lbls[[i + 1]][1], sep="_"), entry9, interf_g)

    button1 <- gtkButtonNewWithLabel("...")
    gSignalConnect(button1, "clicked", select_file_g, lbls[[i + 1]][1])
    table1$Attach(button1, 2, 3, 4, 5, 0, 0, 0, 0)
  }

  hbox1 <- crea_controllo_g(WP, FALSE, TRUE)
  box$packStart(hbox1, FALSE, FALSE, 10)

  assistant$appendPage(box)
  assistant$setPageTitle(box, "Dynamics")
  assistant$setPageComplete(box, TRUE)

  pixbuf <- assistant$renderIcon(GTK_STOCK_DIALOG_INFO, "dialog")
  assistant$setPageHeaderImage(box, pixbuf)
}


select_file1_g <- function(widget, userData)
{
  dialog <- gtkFileChooserDialog("Select file", NULL, "open",
                                    "gtk-cancel", GtkResponseType["cancel"],
                                    "gtk-open", GtkResponseType["accept"])
  entry <- get("file_times", interf_g)
  response <- dialog$run()
  if (response == GtkResponseType["accept"]) {
    nomefile <- dialog$getFilename()
    gtkEntrySetText(entry, nomefile)
  }
  dialog$destroy()
}

create_page5_g <- function(assistant)
{
  box <- gtkVBox(FALSE, 12)
  box$setBorderWidth(12)

  table1 <- gtkTableNew(3, 8, FALSE)
  box$packStart(table1, FALSE, FALSE, 10)
  table1$setBorderWidth(5)
  table1$setRowSpacings(5)
  table1$setColSpacings(5)

  box1 <- gtkHBox(FALSE, 12)
  label3 <- gtkLabelNew ("from")
  label3$setAlignment(0, 0.5)
  box1$packStart(label3, FALSE, FALSE, 10)
  entry1 <- crea_controllo_g(TIMES, FALSE, TRUE)
  gtkEntrySetWidthChars(entry1, 5)
  box1$packStart(entry1, FALSE, FALSE, 10)

  label4 <- gtkLabelNew("step")
  label4$setAlignment(0, 0.5)
  box1$packStart(label4, FALSE, FALSE, 10)
  entry2 <- crea_controllo_g(TIMES + 2, FALSE, TRUE)
  gtkEntrySetWidthChars(entry2, 5)
  box1$packStart(entry2, FALSE, FALSE, 10)

  label5 <- gtkLabelNew ("to")
  label5$setAlignment(0, 0.5)
  box1$packStart(label5, FALSE, FALSE, 10)
  entry3 <- crea_controllo_g(TIMES + 1, FALSE, TRUE)
  gtkEntrySetWidthChars(entry3, 5)
  gSignalConnect(entry3, "focus-out-event", gt_ctrl_g, entry1)
  box1$packStart(entry3, FALSE, FALSE, 10)
  table1$Attach(box1, 2, 3, 0, 1)

  label6 <- gtkLabelNew(glbls_g[[METHOD]][[1]])
  label6$setAlignment(0, 0.5)
  table1$Attach(label6, 0, 1, 2, 3)
  comboboxentry1 <- crea_controllo_g(METHOD, FALSE, TRUE)
  table1$Attach(comboboxentry1, 2, 3, 2, 3)

  button1 <- gtkButtonNewWithLabel("...")
  gSignalConnect(button1, "clicked", select_file1_g)
  table1$Attach(button1, 3, 4, 1, 2)

  entry4 <- crea_controllo_g(TIMES + 3, FALSE, TRUE)
  entry4$setEditable(FALSE)
  gtkEntrySetWidthChars(entry4, 5)
  table1$Attach(entry4, 2, 3, 1, 2)

  radiobutton1 <- gtkRadioButtonNewWithLabel(NULL, "sequence")
  gSignalConnect(radiobutton1, "toggled", toggle_cb_g, c("times", 1))
  assign("toggle_times", 1, interf_g)
  table1$Attach(radiobutton1, 1, 2, 0, 1)
  radiobutton1_group <- gtkRadioButtonGetGroup(radiobutton1)

  radiobutton2 <- gtkRadioButtonNewWithLabel(radiobutton1_group, "file")
  gSignalConnect(radiobutton2, "toggled", toggle_cb_g, c("times", 2))
  table1$Attach(radiobutton2, 1, 2, 1, 2)

  label1 <- gtkLabelNew("times")
  label1$setAlignment(0, 0.5)
  tips_g <- get("tips_g", interf_g)
  gtkTooltipsSetTip(tips_g, label1, "time samples at which explicit estimates of gene expression are desired")
  table1$Attach(label1, 0, 1, 0, 1)

  label6 <- gtkLabelNew(glbls_g[[SIM]][[1]])
  label6$setAlignment(0, 0.5)
  table1$Attach(label6, 0, 1, 3, 4)
  entry <- crea_controllo_g(SIM, FALSE, TRUE)
  table1$Attach(entry, 2, 3, 3, 4)

  label6 <- gtkLabelNew(glbls_g[[SAVE]][[1]])
  label6$setAlignment(0, 0.5)
  table1$Attach(label6, 0, 1, 4, 5)
  chk <- crea_controllo_g(SAVE, FALSE, TRUE)
  table1$Attach(chk, 2, 3, 4, 5)

  assistant$appendPage(box)
  assistant$setPageTitle(box, "Solutions and Simulations")
  assistant$setPageComplete(box, TRUE)

  pixbuf <- assistant$renderIcon(GTK_STOCK_DIALOG_INFO, "dialog")
  assistant$setPageHeaderImage(box, pixbuf)
}

create_page6_g <- function(assistant)
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

  assign("riep", table1, interf_g)

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

netsim_generate <- function()
{
  assistant <- gtkAssistant(show=F)
  assign("main", assistant, interf_g)
  if (!file.exists(paste(getwd(), "netsim", sep="/")))
    dir.create("netsim")
  setwd(paste(getwd(), "netsim", sep="/"))
  unlink("weights*.txt")
  unlink("Rules*.txt")
  unlink("SIMdata*.txt")
  unlink("parameters*.txt")
  unlink("tmp.dot")
  unlink("tmp.png")

  assistant$setDefaultSize(-1, 300)
  gSignalConnect(assistant, "delete-event", onExit)

  # creo un oggetto per i Tooltips
  tips_g <- gtkTooltipsNew()
  assign("tips_g", tips_g, interf_g)
  assign("argom_g", vector("list", TOT_ARG_g), interf_g)

  assign("primo_g", TRUE, interf_g)
  create_page1_g(assistant)
  create_page2_g(assistant)
  create_page3_g(assistant)
  create_page4_g(assistant)
  create_page5_g(assistant)
  create_page6_g(assistant)

  gSignalConnect(assistant, "cancel", onExit)
  gSignalConnect(assistant, "close", onExit)
  gSignalConnect(assistant, "apply", on_assistant_apply_g)
  gSignalConnect(assistant, "prepare", on_assistant_prepare_g)

  assistant$showAll()
}
