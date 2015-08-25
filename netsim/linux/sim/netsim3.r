library(RGtk2)
library(netsim1all)
require(cairoDevice)

ambiente <- new.env()
assign("param_writetable", 0) # ok
assign("param_checkconn", 0) # chi la usa?
assign("param_createNEG", 0) # ok
assign("param_createRules", 1) # ok
assign("param_target", 0) # ok
assign("param_connectivitygeometric", 0) # ok
assign("param_connectivityscalefree", 0) #ok
  assign("param_Score", 0) # ok
  assign("param_Scoresf", 0) # ok
  assign("param_connectivitymodular", 0) # ok
    assign("param_module1", 0) # ok
    assign("param_module2", 0) # ok
    assign("param_module3", 0) # ok
    assign("param_assignnodes", 0) # ok
    assign("param_clustercoeff", 0) # ok
assign("param_connectivityrandom", 0) # ok

interf <- new.env()

dyn.load("calc.dll")

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

gz_ctrl <- function(widget, event, userData)
{
  if (as.double(gtkEntryGetText(widget)) <= 0) {
    message_dialog("the value must be greater than zero")
    widget$grabFocus()
  }
  return (FALSE)
}

gt_ctrl <- function(widget, event, userData)
{
  rif <- as.double(gtkEntryGetText(userData))
  if (as.double(gtkEntryGetText(widget)) <= rif) {
    message_dialog(paste("the value must be greater than", rif))
    widget$grabFocus()
  }
  return (FALSE)
}

act.fun_cb <- function(widget, event, userData)
{
  alpha = get("alpha", envir = interf)
  beta = get("beta", envir = interf)
  txt = gtkComboBoxGetActiveText(widget)
  if (txt == "sigmoidal") {
    gtkWidgetSetSensitive(alpha, TRUE)
    gtkWidgetSetSensitive(beta, TRUE)
  }
  else {
    gtkWidgetSetSensitive(alpha, FALSE)
    gtkWidgetSetSensitive(beta, FALSE)
  }
}

toggle_cb <- function(widget, userData)
{
  if (gtkToggleButtonGetActive(widget)) {
    assign(paste("toggle", userData, sep="_"), gtkButtonGetLabel(widget), envir = interf)
    #cat(sprintf("%s = %s\n", paste("toggle", userData, sep="_"), gtkButtonGetLabel(widget)))
  }
}

# etichette
lbls <- list(
  c("alpha", "vector of parameters for the Activation sigmoid function"),
  c("beta", "vector of parameters for the Activation sigmoid function"),
  c("lambda", "vector of time constants influencing both the rate of transcription and the spontaneous degradation term"),
  c("Xmax", "vector of maximum level of expression of genes"),
  c("X0", "initial conditions, i.e. gene expression values at time 0 scaled between 0 and 1")
)

# variabili globali (1 = testo, 2 = lista, 3 = testo multiplo, 4 = checkbox)
glbls <- list(
  # nome, tipo, default, on_exit, descr
  list("N", 1, 5, gz_ctrl, "number of genes in the network"),
  # nome, tipo, default, valori, call_back, descr
  list("connectivity", 2, 3, c("random", "scale free", "MTM", "geometric"), NULL, "type of connectivity in the network"),
  list("max.reg", 1, 12, NULL, "maximum number of regulators that each node (gene) in the network can have"),
  list("gamma", 1, 2.2, NULL, "the parameter of the power law distribution of the degree of connectivity of the nodes in the graph"),
  list("INdegree", 2, 1, c("free", "out"), NULL, "'free': the in-degree distribution is not constrained to follow any distribution
  'out': it is constrained to follow the same distribution of the out-degree"),
  list("Cf.cl", 1, 0.4, NULL, "average clustering coefficient of each sub-network in the graph"),
  # nome, tipo, defaults
  list("num.subnet", 3, list(5, 5, 10), c("the maximum number of nodes in motif of type 1",
  "the maximum number of nodes in motif of type 2", "the maximum number of nodes in motif of type 3")),
  list("kappa", 1, 3, NULL, "average number of regulators that each node (gene) in the network has with random topology"),
  list("f.pr.and", 1, "0.5", NULL, "function with domain and codomain in [0-1] that expresses the probability to obtain a cooperative rather than a synergic rule"),
  list("act.fun", 2, 2, c("linear", "sigmoidal"), act.fun_cb, "the activation function"),
  list("weight.par", 3, list(1, 0.2), c("mean of the Gaussian distribution used to sample regulatory efficiency",
  "sd of the Gaussian distribution used to sample regulatory efficiency")),
  # 12
  list("alpha_const", 1, 1, gz_ctrl, "constant value"), # vector of parameters of the Activation sigmoid function
  list("alpha_a", 1, 0, NULL, "parameter 'a' for the uniform distribution"),
  list("alpha_b", 1, 1, NULL, "parameter 'b' for the uniform distribution"),
  list("alpha_nm", 1, 10, gz_ctrl, "mean for the normal distribution"),
  list("alpha_nsd", 1, 0.2, gz_ctrl, "std. dev. for the normal distribution"),
  list("alpha_lm", 1, 1, gz_ctrl, "mean for the log-normal distribution"),
  list("alpha_lsd", 1, 0, gz_ctrl, "std. dev. for the log-normal distribution"),
  list("alpha_file", 1, NULL, NULL, "filename"),
  # 20
  list("beta_const", 1, 1, gz_ctrl, "constant value"),
  list("beta_a", 1, 0, NULL, "parameter 'a' for the uniform distribution"),
  list("beta_b", 1, 1, NULL, "parameter 'b' for the uniform distribution"),
  list("beta_nm", 1, 0.5, gz_ctrl, "mean for the normal distribution"),
  list("beta_nsd", 1, 0.01, gz_ctrl, "std. dev. for the normal distribution"),
  list("beta_lm", 1, 1, gz_ctrl, "mean for the log-normal distribution"),
  list("beta_lsd", 1, 0, gz_ctrl, "std. dev. for the log-normal distribution"),
  list("beta_file", 1, NULL, NULL, "filename"),
  # 28
  list("lambda_const", 1, 1, gz_ctrl, "constant value"),
  list("lambda_a", 1, 0, NULL, "parameter 'a' for the uniform distribution"),
  list("lambda_b", 1, 1, NULL, "parameter 'b' for the uniform distribution"),
  list("lambda_nm", 1, 1, gz_ctrl, "mean for the normal distribution"),
  list("lambda_nsd", 1, 0.1, gz_ctrl, "std. dev. for the normal distribution"),
  list("lambda_lm", 1, 1, gz_ctrl, "mean for the log-normal distribution"),
  list("lambda_lsd", 1, 0, gz_ctrl, "std. dev. for the log-normal distribution"),
  list("lambda_file", 1, NULL, NULL, "filename"),
  # 36
  list("Xmax_const", 1, 10, gz_ctrl, "constant value"),
  list("Xmax_a", 1, 0, NULL, "parameter 'a' for the uniform distribution"),
  list("Xmax_b", 1, 1, NULL, "parameter 'b' for the uniform distribution"),
  list("Xmax_nm", 1, 1, gz_ctrl, "mean for the normal distribution"),
  list("Xmax_nsd", 1, 0, gz_ctrl, "std. dev. for the normal distribution"),
  list("Xmax_lm", 1, 1, gz_ctrl, "mean for the log-normal distribution"),
  list("Xmax_lsd", 1, 0, gz_ctrl, "std. dev. for the log-normal distribution"),
  list("Xmax_file", 1, NULL, NULL, "filename"),
  # 44
  list("X0_const", 1, 1, gz_ctrl, "constant value"),
  list("X0_a", 1, 0, NULL, "parameter 'a' for the uniform distribution"),
  list("X0_b", 1, 1, NULL, "parameter 'b' for the uniform distribution"),
  list("X0_nm", 1, 1, gz_ctrl, "mean for the normal distribution"),
  list("X0_nsd", 1, 0, gz_ctrl, "std. dev. for the normal distribution"),
  list("X0_lm", 1, 1, gz_ctrl, "mean for the log-normal distribution"),
  list("X0_lsd", 1, 0, gz_ctrl, "std. dev. for the log-normal distribution"),
  list("X0_file", 1, NULL, NULL, "filename"),
  # 52
  list("from", 1, "0", NULL, "from"),
  list("to", 1, "5", NULL, "to"),
  list("step", 1, "0.05", NULL, "step"),
  list("file_times", 1, "0", NULL, "filename"),
  list("method", 2, 1, c("lsoda", "euler"), NULL, "Method used to solve differential equations"),
  # 57
  list("simulations", 1, "3", NULL, "number of simulations")
)

titoli <- c("Topology", "Topology", "Rules", "Dynamics", "Solutions", "Simulations", "Done")

riepilogo <- function(table1, num, num1, x, y)
{
  argom <- get("argom", interf)
  ctrl <- get(glbls[[num]][[1]], envir = interf)
  txt_lbl <- glbls[[num]][[1]]
  if (num >= 12 && num <= 51) {
    l <- (num - 4) / 8
    txt_lbl <- lbls[[l]]
    radiob <- get(paste("toggle", lbls[[l]], sep="_"), envir = interf)
    txt <- radiob
    N <- gtkEntryGetText(get(glbls[[1]][[1]], envir = interf))
    if (txt == "constant")
      txt1 <- sprintf("rep(%s, %s)", gtkEntryGetText(get(glbls[[num]][[1]], envir = interf)), N)
    else if (txt == "uniform") {
      a <- gtkEntryGetText(get(glbls[[num + 2]][[1]], envir = interf))
      b <- gtkEntryGetText(get(glbls[[num + 3]][[1]], envir = interf))
      txt1 <- sprintf("runif(%s, %s, %s)", N, a, b)
    }
    else if (txt == "normal") {
      m <- gtkEntryGetText(get(glbls[[num + 4]][[1]], envir = interf))
      s <- gtkEntryGetText(get(glbls[[num + 5]][[1]], envir = interf))
      txt1 <- sprintf("abs(rnorm(%s, %s, %s))", N, m, s)
    }
    else if (txt == "log-normal") {
      m <- gtkEntryGetText(get(glbls[[num + 6]][[1]], envir = interf))
      s <- gtkEntryGetText(get(glbls[[num + 7]][[1]], envir = interf))
      txt1 <- sprintf("abs(rlnorm(%s, %s, %s))", N, m, s)
    }
    else if (txt == "file") {
      txt1 <- sprintf("read.table('%s')", gtkEntryGetText(get(glbls[[num + 7]][[1]], envir = interf)))
    }
    if (get("primo", interf))
      argom[[num1]] <- list("espressione", txt1, NULL)
    else
      argom[[num1]][[2]] <- txt1
  }
  else if (num == 7 && gtkWidgetIsSensitive(ctrl[[1]])) {
    txt1 <- sprintf("c(%3.3f, %3.3f, %3.3f)", as.double(gtkEntryGetText(ctrl[[1]])),
      as.double(gtkEntryGetText(ctrl[[2]])), as.double(gtkEntryGetText(ctrl[[3]])))
    if (get("primo", interf))
      argom[[num1]] <- list("espressione", txt1, NULL)
    else
      argom[[num1]][[2]] <- txt1
  }
  else if (num == 11 && gtkWidgetIsSensitive(ctrl[[1]])) {
    txt1 <- sprintf("c(%3.3f, %3.3f)", as.double(gtkEntryGetText(ctrl[[1]])),
      as.double(gtkEntryGetText(ctrl[[2]])))
    if (get("primo", interf))
      argom[[num1]] <- list("espressione", txt1, NULL)
    else
      argom[[num1]][[2]] <- txt1
  }
  else if (num == 9) {
    txt1 <- gtkEntryGetText(ctrl)
    if (get("primo", interf))
      argom[[num1]] <- list("funzione", txt1, NULL)
    else
      argom[[num1]][[2]] <- txt1
  }
  else if (num >= 52 && num <= 55) {
    txt_lbl <- "times"
    txt <- get("toggle_times", envir = interf)
    if (txt == "sequence") {
      da <- gtkEntryGetText(get(glbls[[52]][[1]], envir = interf))
      a <- gtkEntryGetText(get(glbls[[53]][[1]], envir = interf))
      passo <- gtkEntryGetText(get(glbls[[54]][[1]], envir = interf))
      txt1 <- sprintf("seq(%s, %s, %s)", da, a, passo)
    }
    else
      txt1 <- sprintf("read.table('%s')", gtkEntryGetText(get(glbls[[55]][[1]], envir = interf)))
    if (get("primo", interf))
      argom[[num1]] <- list("espressione", txt1, NULL)
    else
      argom[[num1]][[2]] = txt1
  }
  else if (gtkWidgetIsSensitive(ctrl)) {
    if (glbls[[num]][[2]] == 1) {
      txt1 <- gtkEntryGetText(ctrl)
      if (get("primo", interf))
        argom[[num1]] <- list("numero", sprintf('%3.3f', txt1), NULL)
      else
        argom[[num1]][[2]] = sprintf('%3.3f', txt1)
    }
    else if (glbls[[num]][[2]] == 2) {
      txt1 <- gtkComboBoxGetActiveText(ctrl)
      if (get("primo", interf))
        argom[[num1]] <- list("stringa", as.character(txt1), NULL)
      else
        argom[[num1]][[2]] = as.character(txt1)
    }
    else {
      txt1 <- "-"
      if (get("primo", interf))
        argom[[num1]] <- list("saltato", glbls[[num]][[3]], NULL)
    }
  }
  else {
    txt1 <- "-"
    if (get("primo", interf))
      argom[[num1]] <- list("saltato", glbls[[num]][[3]], NULL)
  }
  if (get("primo", interf)) {
    label1 <- gtkLabelNew()
    label1$setAlignment(0, 0.5)
    gtkLabelSetText(label1, txt_lbl)
    table1$Attach(label1, 4 * x, 4 * x + 1, y, y + 1)
    label2 <- gtkLabelNew()
    label2$setAlignment(1, 1)
    argom[[num1]][[3]] <- label2
    table1$Attach(label2, 4 * x + 2, 4 * x + 3, y, y + 1)
  }
  else {
    label2 <- argom[[num1]][[3]]
  }
  assign("argom", argom, interf)
  gtkLabelSetText(label2, txt1)
  return
}

on_assistant_apply <- function(widget, data)
{
  a <- list()
  argom <- get("argom", interf)
  argom[[16]][[2]] <- sprintf("%s / %s", argom[[16]][[2]], argom[[15]][[2]])
  for (l in 1:18) {
    if (argom[[l]][[1]] == "espressione" || argom[[l]][[1]] == "numero") {
       expr <- try(parse(text=argom[[l]][[2]]))
       a[[l]] <- eval(expr)
    }
    #else if (argom[[l]][[1]] == "funzione") {
    #   a[[l]] <- eval(try(parse(text=paste("function(x) { return(", argom[[l]][[2]], ")}"))))
    #}
    else
      a[[l]] <- argom[[l]][[2]]
  }
  cat("\nArgomenti: \n", paste(a, "\n"), "\n")
  volte <- as.integer(gtkEntryGetText(get(glbls[[57]][[1]], interf)))
  for (i in seq(volte)) {
    cat(sprintf("\tsimulatenet(..., TRUE, %d)\n", i))
    a1 <- c(a, TRUE, i)
    do.call("simulatenet", a1)
  }
  gtkWidgetDestroy(get("main", interf))
  cat("\nBene, i risultati sono stati memorizzati nella cartella '", getwd(), "'")
}

on_assistant_prepare <- function(widget, page, data)
{
  current_page <- widget$getCurrentPage()
  n_pages <- widget$getNPages()
  if (current_page == 1) {
    indeg <- get("INdegree", envir = interf)
    gamma <- get("gamma", envir = interf)
    k <- get("kappa", envir = interf)
    txt = gtkComboBoxGetActiveText(get("connectivity", envir = interf))
    cl <- get("Cf.cl", envir = interf)
    for (c in get("num.subnet", envir = interf))
      gtkWidgetSetSensitive(c, FALSE)
    gtkWidgetSetSensitive(indeg, FALSE)
    gtkWidgetSetSensitive(gamma, FALSE)
    gtkWidgetSetSensitive(k, FALSE)
    gtkWidgetSetSensitive(cl, FALSE)
    if (txt == "MTM") {
      gtkWidgetSetSensitive(indeg, TRUE)
      gtkWidgetSetSensitive(gamma, TRUE)
      gtkWidgetSetSensitive(cl, TRUE)
       for (c in get("num.subnet", envir = interf))
        gtkWidgetSetSensitive(c, TRUE)
    }
    else if (txt == "scale free")
      gtkWidgetSetSensitive(gamma, TRUE)
    else if (txt == "random")
      gtkWidgetSetSensitive(k, TRUE)
    else if (txt == "geometric")
      gtkWidgetSetSensitive(k, TRUE)
  }
  else if (current_page == 6) {
    k <- 1
    arg <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 36, 28, 44, 11, 10, 12, 20, 52, 56)
    table1 <- get("riep", envir = interf)
    for (j in 0:1) {
      for (i in 0:8) {
        if (k <= length(arg)) {
          riepilogo(table1, arg[k], k, j, i)
          k <- k + 1
        }
      }
    }
    assign("primo", FALSE, interf)
  }
  title <- sprintf("Netsim assistant - %s (%d of %d)", titoli[current_page + 1], current_page + 1, n_pages)
  widget$setTitle(title)
}

crea_testo <- function(l, tips)
{
  ctrl <- gtkEntryNew()
  gtkEntrySetAlignment(ctrl, 1)
  gtkEntrySetText(ctrl, l[[3]])
  if (!is.null(l[[4]]))
    gSignalConnect(ctrl, "focus-out-event", l[[4]])
  gtkTooltipsSetTip(tips, ctrl, l[[5]])
  assign(l[[1]], ctrl, envir = interf)
}

crea_chk <- function(num, tips)
{
  ctrl <- gtkCheckButtonNew()
  assign(l[[1]], ctrl, envir = interf)
  gtkToggleButtonSetActive(ctrl, l[[3]])
  gtkTooltipsSetTip(tips, ctrl, l[[4]])
  return (ctrl)
}

crea_vett <- function(l, tips)
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
    gtkTooltipsSetTip(tips, ctrl, l[[4]][i])
    algn <- gtkAlignmentNew(0, 0, 0, 0)
    gtkAlignmentSetPadding(algn, 5, 5, 5, 5)
    algn$add(ctrl)
    frm1$packStart(algn, FALSE, TRUE, 0)
    cz <- c(cz, ctrl)
    i <- i + 1
  }
  lbl$add(frm1)
  assign(l[[1]], cz, envir = interf)
  return (frm)
}

crea_lista <- function(l, tips)
{
  # lista combinata
  ctrl <- gtkComboBoxNewText()
  for (t in l[[4]])
    gtkComboBoxAppendText(ctrl, t)
  gtkComboBoxSetActive(ctrl, l[[3]] - 1)
  if (!is.null(l[[5]])) {
     gSignalConnect(ctrl, "changed", l[[5]])
  }
  gtkTooltipsSetTip(tips, ctrl, l[[6]])
  assign(l[[1]], ctrl, envir = interf)
  return (ctrl)
}

crea_controllo <- function(num, etichetta, attivo)
{
  tips <- get("tips", envir = interf)
  l <- glbls[[num]]
  if (etichetta) {
    frm <- gtkHBoxNew(FALSE, 0)
    # etichetta
    lbl <- gtkLabelNew(l[[1]])
    algn1 <- gtkAlignmentNew(0, 0, 0, 0)
    gtkAlignmentSetPadding(algn1, 5, 5, 5, 5)
    algn1$add(lbl)
  }
  if (l[[2]] == 1) {
    ctrl <- crea_testo(l, tips)
  }
  else if (l[[2]] == 2) {
    ctrl <- crea_lista(l, tips)
  }
  else if (l[[2]] == 3) {
    ctrl <- crea_vett(l, tips)
  }
  else if (l[[2]] == 4) {
    ctrl <- crea_chk(l, tips)
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

create_page1 <- function(assistant)
{
  box <- gtkVBox(FALSE, 12)
  box$setBorderWidth(12)

  entry1 <- crea_controllo(1, TRUE, TRUE)
  box$packStart(entry1, FALSE, FALSE, 0)
  entry2 <- crea_controllo(2, TRUE, TRUE)
  box$packStart(entry2, FALSE, FALSE, 0)
  entry3 <- crea_controllo(3, TRUE, TRUE)
  box$packStart(entry3, FALSE, FALSE, 0)

  assistant$appendPage(box)
  assistant$setPageTitle(box, "Topology: common parameters")
  assistant$setPageType(box, "intro")
  assistant$setPageComplete(box, TRUE)

  pixbuf <- assistant$renderIcon(GTK_STOCK_DIALOG_INFO, "dialog")
  assistant$setPageHeaderImage(box, pixbuf)
}

create_page2 <- function(assistant)
{
  box <- gtkVBox(FALSE, 12)
  box$setBorderWidth(12)

  entry1 <- crea_controllo(4, TRUE, TRUE)
  box$packStart(entry1, FALSE, FALSE, 0)
  entry2 <- crea_controllo(5, TRUE, TRUE)
  box$packStart(entry2, FALSE, FALSE, 0)
  entry3 <- crea_controllo(6, TRUE, TRUE)
  box$packStart(entry3, FALSE, FALSE, 0)
  entry3 <- crea_controllo(7, FALSE, TRUE)
  box$packStart(entry3, FALSE, FALSE, 0)
  entry3 <- crea_controllo(8, TRUE, FALSE)
  box$packStart(entry3, FALSE, FALSE, 0)

  assistant$appendPage(box)
  #assistant$setPageType(box, "content")
  assistant$setPageComplete(box, TRUE)
  assistant$setPageTitle(box, "Topology: specific parameters")

  pixbuf <- assistant$renderIcon(GTK_STOCK_DIALOG_INFO, "dialog")
  assistant$setPageHeaderImage(box, pixbuf)
}

message_dialog <- function(testo)
{
  w <- get("main", envir = interf)
  dialog <- gtkMessageDialogNew(w, c("modal", "destroy-with-parent"), "error", "", testo)
  gSignalConnectSwapped(dialog, "response", gtkWidgetDestroy)
  dialog$run()
  dialog$destroy()
}

# valutazione di espressioni
run <- function() {
  code <- gtkEntryGetText(get("f.pr.and", envir = interf))
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
    message_dialog(geterrmessage())
    return()
  }
  return (ris)
}

# funzione di call-back per il pulsante "Funzione"
preview <- function(widget, userData)
{
  x <- (0 : 100) / 100
  y <- run()
  plot(x, y, type="l")
}

create_page3 <- function(assistant)
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
  entry2 <- crea_controllo(9, TRUE, TRUE)
  box1$packStart(entry2, FALSE, FALSE, 0)
  entry3 <- gtkButton("Disegna")
  gSignalConnect(entry3, "clicked", preview)
  box1$packStart(entry3, FALSE, FALSE, 0)
  box$packStart(box1, FALSE, FALSE, 0)

  assistant$appendPage(box)
  #assistant$setPageType(box, "confirm")
  assistant$setPageComplete(box, TRUE)
  assistant$setPageTitle(box, "Rules")

  pixbuf <- assistant$renderIcon(GTK_STOCK_DIALOG_INFO, "dialog")
  assistant$setPageHeaderImage(box, pixbuf)
}

select_file <- function(widget, userData)
{
  dialog <- gtkFileChooserDialog("Select file", NULL, "open",
                                    "gtk-cancel", GtkResponseType["cancel"],
                                    "gtk-open", GtkResponseType["accept"])
  fn <- paste(userData, "file", sep="_")
  cat(fn)
  entry <- get(fn, envir = interf)
  response <- dialog$run()
  if (response == GtkResponseType["accept"]) {
    nomefile <- dialog$getFilename()
    gtkEntrySetText(entry, nomefile)
  }
  dialog$destroy()
}

create_page4 <- function(assistant)
{
  box <- gtkVBoxNew(FALSE, 0)
  hbox1 <- crea_controllo(10, TRUE, TRUE)
  box$packStart(hbox1, FALSE, FALSE, 10)

  notebook1 <- gtkNotebook()
  notebook1$show()
  box$packStart(notebook1, FALSE, FALSE, 0)

  tips <- get("tips", envir = interf)
  def <- c(2, 2, 2, 0, 1)
  for (i in 0:4) {

    table1 <- gtkTableNew(5, 5, FALSE)

    tab1 <- gtkLabel(lbls[[i + 1]][1])
    assign(lbls[[i + 1]][1], table1, envir = interf)
    gtkTooltipsSetTip(tips, tab1, lbls[[i + 1]][2])
    notebook1$AppendPage(table1, tab1)
    table1$setRowSpacings(5)
    table1$setColSpacings(5)
    box$packStart(table1, FALSE, FALSE, 0)

    radiobutton1 <- gtkRadioButtonNewWithLabel(NULL, "constant")
    gSignalConnect(radiobutton1, "toggled", toggle_cb, lbls[[i + 1]][1])
    assign(paste("toggle", lbls[[i + 1]][1], sep="_"), "constant", envir = interf)
    table1$Attach(radiobutton1, 0, 1, 0, 1, xpadding = 5, ypadding = 5)
    radiobutton1_group <- gtkRadioButtonGetGroup(radiobutton1)
    #radiobutton1$setGroup(radiobutton1_group)
    if (def[i + 1] == 0)
      gtkToggleButtonSetActive(radiobutton1, TRUE)

    radiobutton2 <- gtkRadioButtonNewWithLabel(radiobutton1_group, "uniform")
    gSignalConnect(radiobutton2, "toggled", toggle_cb, lbls[[i + 1]][1])
    table1$Attach(radiobutton2, 0, 1, 1, 2, xpadding = 5, ypadding = 5)
    gtkRadioButtonSetGroup(radiobutton2, radiobutton1_group)
    if (def[i + 1] == 1)
      gtkToggleButtonSetActive(radiobutton2, TRUE)

    radiobutton3 <- gtkRadioButtonNewWithLabel(radiobutton1_group, "normal")
    gSignalConnect(radiobutton3, "toggled", toggle_cb, lbls[[i + 1]][1])
    table1$Attach(radiobutton3, 0, 1, 2, 3, xpadding = 5, ypadding = 5)
    radiobutton3$setGroup(radiobutton1_group)
    #radiobutton1_group <- gtk_radio_button_get_group(radiobutton3)
    if (def[i + 1] == 2)
      gtkToggleButtonSetActive(radiobutton3, TRUE)

    radiobutton4 <- gtkRadioButtonNewWithLabel(radiobutton1_group, "log-normal")
    gSignalConnect(radiobutton4, "toggled", toggle_cb, lbls[[i + 1]][1])
    table1$Attach(radiobutton4, 0, 1, 3, 4, xpadding = 5, ypadding = 5)
    radiobutton4$setGroup(radiobutton1_group)
    #radiobutton1_group <- gtk_radio_button_get_group(radiobutton4)

    radiobutton5 <- gtkRadioButtonNewWithLabel(radiobutton1_group, "file")
    gSignalConnect(radiobutton5, "toggled", toggle_cb, lbls[[i + 1]][1])
    table1$Attach(radiobutton5, 0, 1, 4, 5, xpadding = 5, ypadding = 5)
    radiobutton5$setGroup(radiobutton1_group)
    #radiobutton1_group <- gtk_radio_button_get_group(radiobutton5)

    entry1 <- crea_controllo(12 + 8 * i, FALSE, TRUE)
    gtkEntrySetWidthChars(entry1, 5)
    table1$Attach(entry1, 1, 2, 0, 1)

    entry2 <- crea_controllo(13 + 8 * i, FALSE, TRUE)
    gtkEntrySetWidthChars(entry2, 5)
    table1$Attach(entry2, 1, 2, 1, 2)

    entry3 <- crea_controllo(14 + 8 * i, FALSE, TRUE)
    gtkEntrySetWidthChars(entry3, 5)
    table1$Attach(entry3, 3, 4, 1, 2, xpadding = 5, ypadding = 5)

    entry4 <- crea_controllo(15 + 8 * i, FALSE, TRUE)
    gtkEntrySetWidthChars(entry4, 5)
    table1$Attach(entry4, 1, 2, 2, 3)
    entry6 <- crea_controllo(16 + 8 * i, FALSE, TRUE)
    gtkEntrySetWidthChars(entry6, 5)
    table1$Attach(entry6, 3, 4, 2, 3, xpadding = 5, ypadding = 5)

    entry7 <- crea_controllo(17 + 8 * i, FALSE, TRUE)
    gtkEntrySetWidthChars(entry7, 5)
    table1$Attach(entry7, 1, 2, 3, 4)

    entry8 <- crea_controllo(18 + 8 * i, FALSE, TRUE)
    gtkEntrySetWidthChars(entry8, 5)
    table1$Attach(entry8, 3, 4, 3, 4, xpadding = 5, ypadding = 5)

    entry9 <- crea_controllo(19 + 8 * i, FALSE, TRUE)
    gtkEntrySetWidthChars(entry9, 5)
    table1$Attach(entry9, 1, 2, 4, 5)
    assign(paste("file", lbls[[i + 1]][1], sep="_"), entry9, envir = interf)

    button1 <- gtkButtonNewWithLabel("...")
    gSignalConnect(button1, "clicked", select_file, lbls[[i + 1]][1])
    table1$Attach(button1, 2, 3, 4, 5, 0, 0, 0, 0)
  }

  hbox1 <- crea_controllo(11, FALSE, TRUE)
  box$packStart(hbox1, FALSE, FALSE, 10)

  assistant$appendPage(box)
  assistant$setPageTitle(box, "Dynamics")
  assistant$setPageComplete(box, TRUE)

  pixbuf <- assistant$renderIcon(GTK_STOCK_DIALOG_INFO, "dialog")
  assistant$setPageHeaderImage(box, pixbuf)
}


select_file1 <- function(widget, userData)
{
  dialog <- gtkFileChooserDialog("Select file", NULL, "open",
                                    "gtk-cancel", GtkResponseType["cancel"],
                                    "gtk-open", GtkResponseType["accept"])
  entry <- get("file_times", envir = interf)
  response <- dialog$run()
  if (response == GtkResponseType["accept"]) {
    nomefile <- dialog$getFilename()
    gtkEntrySetText(entry, nomefile)
  }
  dialog$destroy()
}

create_page5 <- function(assistant)
{
  box <- gtkVBox(FALSE, 12)
  box$setBorderWidth(12)

  table1 <- gtkTableNew(3, 8, FALSE)
  box$packStart(table1, FALSE, FALSE, 10)
  table1$setBorderWidth(5)
  table1$setRowSpacings(5)
  table1$setColSpacings(5)

  label4 <- gtkLabelNew("to")
  label4$setAlignment(0, 0.5)
  table1$Attach(label4, 4, 5, 0, 1)

  label5 <- gtkLabelNew ("step")
  label5$setAlignment(0, 0.5)
  table1$Attach(label5, 6, 7, 0, 1)

  label3 <- gtkLabelNew ("from")
  label3$setAlignment(0, 0.5)
  table1$Attach(label3, 2, 3, 0, 1)

  label6 <- gtkLabelNew(glbls[[56]][[1]])
  label6$setAlignment(0, 0.5)
  table1$Attach(label6, 0, 1, 2, 3)
  comboboxentry1 <- crea_controllo(56, FALSE, TRUE)
  table1$Attach(comboboxentry1, 1, 2, 2, 3)

  entry1 <- crea_controllo(52, FALSE, TRUE)
  gtkEntrySetWidthChars(entry1, 5)
  table1$Attach(entry1, 3, 4, 0, 1)

  entry3 <- crea_controllo(53, FALSE, TRUE)
  gtkEntrySetWidthChars(entry3, 5)
  gSignalConnect(entry3, "focus-out-event", gt_ctrl, entry1)
  table1$Attach(entry3, 5, 6, 0, 1)

  entry2 <- crea_controllo(54, FALSE, TRUE)
  gtkEntrySetWidthChars(entry2, 5)
  table1$Attach(entry2, 7, 8, 0, 1)

  button1 <- gtkButtonNewWithLabel("...")
  gSignalConnect(button1, "clicked", select_file1)
  table1$Attach(button1, 3, 4, 1, 2)

  entry4 <- crea_controllo(55, FALSE, TRUE)
  gtkEntrySetWidthChars(entry4, 5)
  table1$Attach(entry4, 2, 3, 1, 2)

  radiobutton1 <- gtkRadioButtonNewWithLabel(NULL, "sequence")
  gSignalConnect(radiobutton1, "toggled", toggle_cb, "time")
  assign("toggle_times", "sequence", envir = interf)
  table1$Attach(radiobutton1, 1, 2, 0, 1)
  radiobutton1_group <- gtkRadioButtonGetGroup(radiobutton1)

  radiobutton2 <- gtkRadioButtonNewWithLabel(radiobutton1_group, "file")
  gSignalConnect(radiobutton2, "toggled", toggle_cb, "time")
  table1$Attach(radiobutton2, 1, 2, 1, 2)

  label1 <- gtkLabelNew("times")
  label1$setAlignment(0, 0.5)
  tips <- get("tips", envir = interf)
  gtkTooltipsSetTip(tips, label1, "time samples at which explicit estimates of gene expression are desired")
  table1$Attach(label1, 0, 1, 0, 1)


  assistant$appendPage(box)
  assistant$setPageTitle(box, "Solutions")
  assistant$setPageComplete(box, TRUE)

  pixbuf <- assistant$renderIcon(GTK_STOCK_DIALOG_INFO, "dialog")
  assistant$setPageHeaderImage(box, pixbuf)
}

create_page6 <- function(assistant)
{
  table1 <- gtkTableNew(1, 1, FALSE)
  table1$setBorderWidth(5)
  table1$setRowSpacings(5)
  table1$setColSpacings(5)
  entry1 <- gtkLabel("How many simulations do you want?")
  table1$Attach(entry1, 0, 1, 0, 1)
  entry2 <- crea_controllo(57, FALSE, TRUE)
  table1$Attach(entry2, 1, 2, 0, 1)

  assistant$appendPage(table1)
  assistant$setPageComplete(table1, TRUE)
  assistant$setPageTitle(table1, "Simulations")

  pixbuf <- assistant$renderIcon(GTK_STOCK_DIALOG_INFO, "dialog")
  assistant$setPageHeaderImage(table1, pixbuf)
}

create_page7 <- function(assistant)
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

  assign("riep", table1, envir = interf)

  assistant$appendPage(box)
  assistant$setPageType(box, "confirm")
  assistant$setPageComplete(box, TRUE)
  assistant$setPageTitle(box, "Confirmation")

  pixbuf <- assistant$renderIcon(GTK_STOCK_DIALOG_INFO, "dialog")
  assistant$setPageHeaderImage(box, pixbuf)
}

start <- function()
{
  assistant <- gtkAssistant(show = F)
  assign("main", assistant, envir = interf)

  assistant$setDefaultSize(-1, 300)

  # creo un oggetto per i Tooltips
  tips <- gtkTooltipsNew()
  assign("tips", tips, envir = interf)
  assign("primo", TRUE, interf)
  assign("argom", vector("list", 19), interf)

  create_page1(assistant)
  create_page2(assistant)
  create_page3(assistant)
  create_page4(assistant)
  create_page5(assistant)
  create_page6(assistant)
  create_page7(assistant)

  gSignalConnect(assistant, "cancel", gtkWidgetDestroy)
  gSignalConnect(assistant, "close", gtkWidgetDestroy)
  gSignalConnect(assistant, "apply", on_assistant_apply)
  gSignalConnect(assistant, "prepare", on_assistant_prepare)

  assistant$showAll()
}

#debug(on_assistant_apply)
start()