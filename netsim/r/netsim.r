library(RGtk2)
library(netsim)

# variabili globali (1 = testo, 2 = lista, 3 = testo multiplo, 4 = checkbox)
glbls <- list(
  # nome, tipo, default, descr
  list("N", 1, 5, "number of genes in the network"),
  # nome, tipo, default, valori, call_back, descr
  list("connectivity", 2, 3, c("random", "scale free", "MTM", "geometric"), conn_cb, "type of connectivity in the network"),
  list("max.reg", 1, 12, "maximum number of regulators that each node (gene) in the network can have"),
  list("gamma", 1, 2.2, "the parameter of the power law distribution of the degree of connectivity of the nodes in the graph"),
  list("INdegree", 2, 1, c("free", "out"), NULL, "'free': the in-degree distribution is not constrained to follow any distribution
'out': it is constrained to follow the same distribution of the out-degree"),
  list("Cf.cl  ", 1, 0.4, "average clustering coefficient of each sub-network in the graph"),
  # nome, tipo, defaults
  list("num.subnet", 3, list(5, 5, 10), c("the maximum number of nodes in motif of type 1",
"the maximum number of nodes in motif of type 2", "the maximum number of nodes in motif of type 3")),
  list("kappa", 1, 3, "average number of regulators that each node (gene) in the network has with random topology"),
  list("f.pr.and", 1, NULL, "function with domain and codomain in [0-1] that expresses the probability to obtain a cooperative rather than a synergic rule"),
  list("Xmax", 1, NULL, "vector of maximum level of expression of genes"),
  list("lambda", 1, NULL, "vector of time constants influencing both the rate of transcription and the spontaneous degradation term"),
  list("x0", 1, NULL, "initial conditions, i.e. gene expression values at time 0 scaled between 0 and 1"),
  list("weight.par", 3, list(1, 0.2), c("mean of the Gaussian distribution used to sample regulatory efficiency",
 "sd of the Gaussian distribution used to sample regulatory efficiency")),
  list("act.fun", 2, 2, c("linear", "sigmoidal"), act.fun_cb, "the activation function"),
  list("alpha", 1, NULL, "vector of parameters of the Activation sigmoid function"),
  list("beta", 1, NULL, "vector of parameters of the Activation sigmoid function"),
  list("times", 1, "seq(0, 0.5, 0.05)", " time samples at which explicit estimates of gene expression are desired"),
  list("method", 2, 1, c("lsoda", "euler"), NULL, "Method used to solve differential equations"),
  list("save", 4, TRUE, "if TRUE the number specified in 'ind.itera' is appended to the output file names"),
  list("ind.itera", 1, 1, "number useful to differentiate saved output files in case of multiple runs of the function 'simulatenet'")
)

# funzione di call-back
conn_cb <- function(widget, event)
{
  indeg <- get("INdegree", envir = .GlobalEnv)
  gamma <- get("gamma", envir = .GlobalEnv)
  k <- get("kappa", envir = .GlobalEnv)
  txt = gtkComboBoxGetActiveText(widget)
  gtkWidgetSetSensitive(indeg, FALSE)
  gtkWidgetSetSensitive(gamma, FALSE)
  gtkWidgetSetSensitive(k, FALSE)
  if (txt == "MTM") {
	gtkWidgetSetSensitive(indeg, TRUE)
	gtkWidgetSetSensitive(gamma, TRUE)
  }
  else if (txt == "scale free") {
	gtkWidgetSetSensitive(gamma, TRUE)
  }
  else if (txt == "random") {
	gtkWidgetSetSensitive(k, TRUE)
  }
}

# funzione di call-back
act.fun_cb <- function(widget, event)
{
  alpha <- get("alpha", envir = .GlobalEnv)
  beta <- get("beta", envir = .GlobalEnv)
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

# funzione di call-back per il pulsante "Annulla"
Annulla <- function(widget, event)
{
  gtkWidgetDestroy(get("main", envir = .GlobalEnv))
}

# funzione di call-back per il pulsante "Bene"
on_click <- function(widget, event)
{
  args <- list()
  i <- 1
  for (l in glbls) {
    ctrl <- get(l[[1]], envir = .GlobalEnv)
    if (l[[2]] == 1) {
      if (gtkEntryGetText(ctrl) == "")
         args[[i]] <- NULL
      else if (l[[1]] == "f.pr.and")
         args[[i]] <- gtkEntryGetText(expr)
      else if (l[[1]] == "times") {
         expr <- try(parse(text=gtkEntryGetText(ctrl)))
         args[[i]] <- eval(expr)
      }
	else
         args[[i]] <- as.integer(gtkEntryGetText(ctrl))
    }
    else if (l[[2]] == 2)
      args[[i]] <- as.character(gtkComboBoxGetActiveText(ctrl))
    else if (l[[2]] == 3) {
      testo <- vector()
      for (c in ctrl)
        testo <- c(testo, as.integer(gtkEntryGetText(c)))
      args[[i]] <- testo
    }
    else if (l[[2]] == 4)
      args[[i]] <- gtkToggleButtonGetActive(ctrl)
    i <- i + 1
  }
  cat("\nArgomenti: ", paste(args, ", "))
  do.call("simulatenet", args)
  gtkWidgetDestroy(get("main", envir = .GlobalEnv))
  cat("\nBene, i risultati sono stati memorizzati nella cartella '", getwd(), "'")
}

variabile <- function(num)
{
  l <- glbls[[num]]
  tips <- get("tips", envir = .GlobalEnv)
  if (l[[2]] == 1) {
    frm <- gtkHBoxNew(FALSE, 10)
    # etichetta
    lbl <- gtkLabelNew(l[[1]])
    # testo
    ctrl <- gtkEntryNew()
    assign(l[[1]], ctrl, envir = .GlobalEnv)
    gtkEntrySetAlignment(ctrl, 1)
    gtkEntrySetText(ctrl, l[[3]])
    gtkTooltipsSetTip(tips, ctrl, l[[4]])
    if (l[[1]] == "kappa")
        gtkWidgetSetSensitive(ctrl, FALSE)
    algn1 <- gtkAlignmentNew(0, 0, 0, 0)
    gtkAlignmentSetPadding(algn1, 5, 5, 5, 5)
    algn1$add(lbl)
    algn2 <- gtkAlignmentNew(0, 0, 0, 0)
    gtkAlignmentSetPadding(algn2, 5, 5, 5, 5)
    algn2$add(ctrl)
    frm$packStart(algn1, TRUE, TRUE, 0)
    frm$packStart(algn2, FALSE, FALSE, 0)

  }
  else if (l[[2]] == 4) {
    frm <- gtkHBoxNew(FALSE, 0)
    # etichetta
    lbl <- gtkLabelNew(l[[1]])
    # testo
    ctrl <- gtkCheckButtonNew()
    assign(l[[1]], ctrl, envir = .GlobalEnv)
    gtkToggleButtonSetActive(ctrl, l[[3]])
    gtkTooltipsSetTip(tips, ctrl, l[[4]])
    algn1 <- gtkAlignmentNew(0, 0, 0, 0)
    gtkAlignmentSetPadding(algn1, 5, 5, 5, 5)
    algn1$add(lbl)
    algn2 <- gtkAlignmentNew(0, 0, 0, 0)
    gtkAlignmentSetPadding(algn2, 5, 5, 5, 5)
    algn2$add(ctrl)
    frm$packStart(algn1, TRUE, TRUE, 0)
    frm$packStart(algn2, TRUE, TRUE, 0)
  }
  else if (l[[2]] == 3) {
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
	frm1$packStart(algn, TRUE, TRUE, 0)
      cz <- c(cz, ctrl)
      i <- i + 1
    }
    lbl$add(frm1)
    assign(l[[1]], cz, envir = .GlobalEnv)
  }
  else if (l[[2]] == 2) {
    frm <- gtkHBoxNew(FALSE, 0)
    # etichetta
    lbl <- gtkLabelNew(l[[1]])
    # lista combinata
    ctrl <- gtkComboBoxNewText()
    assign(l[[1]], ctrl, envir = .GlobalEnv)
    for (t in l[[4]])
      gtkComboBoxAppendText(ctrl, t)
    gtkComboBoxSetActive(ctrl, l[[3]] - 1)
    gtkTooltipsSetTip(tips, ctrl, l[[6]])
    if (!is.null(l[[5]])) {
       gSignalConnect(ctrl, "changed", l[[5]])
    }
    algn1 <- gtkAlignmentNew(0, 0, 0, 0)
    gtkAlignmentSetPadding(algn1, 5, 5, 5, 5)
    algn1$add(lbl)
    algn2 <- gtkAlignmentNew(0, 0, 0, 0)
    gtkAlignmentSetPadding(algn2, 5, 5, 5, 5)
    algn2$add(ctrl)
    frm$packStart(algn1, TRUE, TRUE, 0)
    frm$packStart(algn2, TRUE, TRUE, 0)
  }
  return (frm)
}

main <- function()
{
  # creo la finestra
  win <- gtkWindowNew("toplevel", show = F)
  assign("main", win, envir = .GlobalEnv)
  # titolo della finestra
  gtkWindowSetTitle(win, "NetSim")

  # creo un oggetto per i Tooltips
  tips <- gtkTooltipsNew()
  assign("tips", tips, envir = .GlobalEnv)

  frm <- gtkHBoxNew(FALSE, 10)
  win$add(frm)
  frm1 <- gtkVBoxNew(FALSE, 10)
  for (i in seq(1, 10, 1)) {
     box <- variabile(i)
     frm1$packStart(box, TRUE, TRUE, 0)
  }
  frm$packStart(frm1, FALSE, FALSE, 0)
  frm2 <- gtkVBoxNew(FALSE, 10)
  for (i in seq(10, 20, 1)) {
     box <- variabile(i)
     frm2$packStart(box, TRUE, TRUE, 0)
  }
  frm$packStart(frm2, FALSE, FALSE, 0)

  button.widget <- gtkButton("Bene")
  gSignalConnect(button.widget, "clicked", on_click)
  frm$packStart(button.widget, FALSE, FALSE, 0)

  cancel.widget <- gtkButton("Annulla")
  gSignalConnect(cancel.widget, "clicked", Annulla)
  frm$packStart(cancel.widget, FALSE, FALSE, 0)

  win$showAll()
}

main()
