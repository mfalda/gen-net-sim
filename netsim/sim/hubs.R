library(RGtk2)
library(cairoDevice)
library(Rgraphviz)
library(graph)


interf_h <- new.env()

col_toggled1 <- function(widget, userData)
{
  spin_changed1_h(userData, widget)
}

col_toggled2 <- function(widget, userData)
{
  spin_changed2_h(userData, widget)
}

spin_changed1_h <- function(widget, userData)
{
  N <- widget$getValueAsInt()
  W <- read.table(sprintf("weights%d.txt", N), row.names=1, header=TRUE)
  n <- nrow(W)
  W1 <- abs(W)
  outd <- apply(W1, 2, sum)
  nAttrs <- list()
  nomi <- as.character(lapply(seq(1, n), function(x) {sprintf("g%d", x)}))
  colori <- vector()
  mx <- max(outd)
  if (userData$getActive()) {
    for (i in 1:n) {
      colore <- sprintf("#%02X%02X%02X", as.integer(255 - 127 / mx * outd[i]), as.integer(255 - 255 / mx * outd[i]), 255)
      colori <- c(colori, colore)
    }
  }
  else
    colori <- rep("white", n)
  nAttrs$fillcolor <- colori
  names(nAttrs$fillcolor) <- nomi
  rEG  <-  new("graphNEL", nodes=nomi, edgemode="directed")
  for (i in 1:n) {
    for (j in 1:n) {
      if (W[i, j] != 0) {
        rEG <- addEdge(nomi[j], nomi[i], rEG, 1)
      }
    }
  }
  require(graph)
  plot(rEG, recipEdges="distinct", nodeAttrs=nAttrs)
}

spin_changed2_h <- function(widget, userData)
{
  pad_w <- get("pad_w", interf_h)
  N <- widget$getValueAsInt()
  W <- read.table(sprintf("weights%0*d.txt", pad_w, N), row.names=1, header=TRUE)
  n <- nrow(W)
  W1 <- abs(W)
  outd <- apply(W1, 2, sum)
  nomi <- as.character(lapply(seq(1, n), function(x) {sprintf("g%d", x)}))
  colori <- vector()
  mx <- max(outd)
  if (userData$getActive()) {
    for (i in 1:n) {
      colore <- sprintf("#%02X%02X%02X", as.integer(255 - 127 / mx * outd[i]), as.integer(255 - 255 / mx * outd[i]), 255)
      colori <- c(colori, colore)
    }
  }
  else
    colori <- rep("white", n)
  f <- file("tmp.dot", "w")
  win <- get("win_h", interf_h)
  w <- win[["allocation"]]$width
  h <- win[["allocation"]]$height
  cat(file=f, sprintf("digraph G {\n\tratio=fill;\n\tsize=\"%3.3f,%3.3f\";\n", w / 75, h / 75))
  for (i in 1:n) {
    cat(file=f, sprintf("\t%s [style=filled,fillcolor=\"%s\"];\n", nomi[i], colori[i]))
  }
  for (i in 1:n) {
    for (j in 1:n) {
      if (W[i, j] != 0) {
        cat(file=f, sprintf("\t%s -> %s;\n", nomi[j], nomi[i], colori[i]))
      }
    }
  }
  cat(file=f, "}\n")
  close(f)
  system(sprintf("dot -Tpng -otmp.png tmp.dot"), ignore.stderr=FALSE, show.output.on.console=FALSE)
  img <- gdkPixbufNewFromFile("tmp.png")$retval
  assign("img_h", img, interf_h)
  da2 <- win[["window"]]
  gc <- gdkGCNew(da2)
  img1 <- gdkPixbufScaleSimple(img, w, h, "bilinear")
  gdkDrawPixbuf(da2, gc, pixbuf=img1, 0, 0, 0, 0, w, h)
}

expose_event_cb <- function(widget, event, userData)
{
  #img <- gdkPixbufNewFromFile("tmp.png")$retval
  if (exists("img_h", interf_h)) {
    img <- get("img_h", interf_h)
    da2 <- widget[["window"]]
    w <- widget[["allocation"]]$width
    h <- widget[["allocation"]]$height
    gc <- gdkGCNew(da2)
    img1 <- gdkPixbufScaleSimple(img, w, h, "bilinear")
    gdkDrawPixbuf(da2, gc, pixbuf=img1, 0, 0, 0, 0, w, h) 
  }
}

onExit <- function(widget, event, userData)
{
  gtkWidgetDestroy(widget)
  setwd("..")
  return (FALSE)
}

netsim_graph <- function(vista=FALSE)
{
  if (!file.exists("netsim")
    stop("To plot the profiles you must first generate them using 'netsim_simulate()' and save at least a network\n(also 'netsim_generate()', but in that case use 'netsim_profiles(sn=TRUE)')!")
  else {
  	 setwd("netsim")
    assign("vista_h", vista, interf_h)
    win <- gtkWindow(show=F)
    win$setTitle("Graph and hubs")
    gSignalConnect(win, "delete-event", onExit)
    da <- gtkDrawingArea()
    label1 <- gtkLabelNew("Network: ")
    mx <- length(i <- grep("weights", dir()))
    pad_w <- ceiling(log10(mx + 1))
    assign("pad_w", pad_w, interf_h)
    spn <- gtkSpinButtonNewWithRange(1, mx, 1)
    spn$SetWidthChars(5)
    hbox <- gtkHBox()
    hbox$packStart(label1, FALSE, TRUE, 5)
    hbox$packStart(spn, FALSE, TRUE, 5)
    colori <- gtkCheckButtonNewWithLabel("hubs")
    if (vista) {
      gSignalConnect(colori, "toggled", col_toggled2, spn)
      gSignalConnect(spn, "value-changed", spin_changed2_h, colori)
      gSignalConnect(da, "expose-event", expose_event_cb)
    }
    else {
      gSignalConnect(colori, "toggled", col_toggled1, spn)
      gSignalConnect(spn, "value-changed", spin_changed1_h, colori)
    }
    hbox$packStart(colori, FALSE, FALSE, 5)
    vbox <- gtkVBox()
    vbox$packStart(da)
    vbox$packStart(hbox, FALSE, FALSE, 5)
    win$add(vbox)
    win$setDefaultSize(650, 500)
    win$showAll()
    asCairoDevice(da)
    da$show()
    assign("win_h", da, interf_h)
    gtkWidgetSetSizeRequest(da, 500, 400)
    par(pty = "s")
    if (vista)
      spin_changed2_h(spn, colori)
    else
      spin_changed1_h(spn, colori)
  }
}
