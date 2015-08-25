library("RGtk2")
library("cairoDevice")

# costanti simboliche
SEMPLICE = 0 # SIMdata?.txt
MULTIPLO = 1 # SIMdata?_?.txt
KO_SEMPLICE = 2 # SIMdata?_ko?.txt
KO_MULTIPLO = 3 # SIMdata?_?_ko?.txt

interf_p <- new.env()

ncl_changed_p <- function(widget, userData)
{
  if (get("clust", interf_p)$getActive()) {
    assign("ricalcola_c", TRUE, interf_p)
    spin_changed_p(get("spn1", interf_p))
  }
  return (FALSE)
}

rete_changed_p <- function(widget, userData)
{
  tipo <- get("tipo", interf_p)
  pads <- get("pads", interf_p)
  prof <- list()
  R <- gtkSpinButtonGetValueAsInt(widget)
  if (tipo == SEMPLICE || tipo == MULTIPLO) {
    if (get("sn", interf_p))
      prof[[1]] <- read.table(sprintf("SIMdata%0*d.txt", pads[1], R), row.names=1, header=TRUE)
    else
      prof[[1]] <- read.table(sprintf("SIMdata_simulateprofiles%0*d_%0*d.txt", pad[1], R, pad[2], 1), row.names=1, header=TRUE)
    media <- as.matrix(prof[[1]])
    j <- 2
    while (j <= userData[[1]]) {
      prof[[j]] <- read.table(sprintf("SIMdata_simulateprofile%0*d_%0*d.txt", pad[1], R, pad[2], j), row.names=1, header=TRUE)
      media <- media + as.matrix(prof[[j]])
      j <- j + 1
    }
    media <- media / as.integer(userData[[1]])
  }
  else {
    media <- list()
    for (i in 1:length(userData[[3]])) {
      prof[[i]] <- list()
      prof[[i]][[1]] <- read.table(sprintf("SIMdata_simulateprofiles%0*d_%0*d_ko%0*d.txt", pad[1], R, pad[2], 1, pad[3], as.integer(userData[[3]][i])), row.names=1, header=TRUE)
      media[[i]] <- as.matrix(prof[[i]][[1]])
      j <- 2
      while (j <= userData[[2]]) {
        prof[[i]][[j]] <- read.table(sprintf("SIMdata_simulateprofiles%0*d_%0*d_ko%0*d.txt", pad[1], R, pad[2], j, pad[3], as.integer(userData[[3]][i])), row.names=1, header=TRUE)
        media[[i]] <- media[[i]] + as.matrix(prof[[i]][[j]])
        j <- j + 1
      }
      media[[i]] <- media[[i]] / as.integer(userData[[2]])
    }
  }
  assign("prof", prof, interf_p)
  assign("medie", media, interf_p)
  assign("weights", read.table(sprintf("weights%0*d.txt", pad[1], R), row.names=1, header=TRUE), interf_p)
  spin_changed_p(get("spn1", interf_p))
}

chk_toggled <- function(widget, userData)
{
  if (widget$name == "clust") {
    if (widget$getActive()) {
      assign("ricalcola_c", TRUE, interf_p)
    }
    else {
      gtkSpinButtonSetRange(get("spn5", interf_p), 0, 0)
    }
  }
  spin_changed_p(userData, widget)
  return (FALSE)
}

message_dialog_p <- function(testo)
{
  w <- get("main", interf_p)
  dialog <- gtkMessageDialogNew(w, c("modal", "destroy-with-parent"), "error", "ok", testo)
  #gSignalConnectSwapped(dialog, "response", gtkWidgetDestroy)
  dialog$run()
  dialog$destroy()
}

spin_changed_p <- function(widget, userData)
{
  N <- gtkSpinButtonGetValueAsInt(get("spn1", interf_p))
  R <- gtkSpinButtonGetValueAsInt(get("spn2", interf_p))
  tipo <- get("tipo", interf_p)
  if (widget$name == "ko") {
    K <- gtkSpinButtonGetValueAsInt(get("spn4", interf_p))
    if (K == 0)
      gtkLabelSetText(get("gn", interf_p), "(none)")
    else {
      ko <- get("ko", interf_p)
      gtkLabelSetText(get("gn", interf_p), sprintf("(g%s)", ko[K + 1]))
    }
  }
  differ <- FALSE
  if (exists("media", interf_p) && get("media", interf_p)$getActive()) {
    if (tipo == 1) {
      P <- get("medie", interf_p)
    }
    else if (tipo == KO_MULTIPLO) {
      K <- gtkSpinButtonGetValueAsInt(get("spn4", interf_p))
      differ <- get("diff", interf_p)$getActive()
      P0 <- get("medie", interf_p)[[1]]
      P <- get("medie", interf_p)[[K + 1]]
    }
  }
  else {
    if (tipo == SEMPLICE)
      P <- get("prof", interf_p)[[1]]
    else if (tipo == MULTIPLO) {
      E <- gtkSpinButtonGetValueAsInt(get("spn3", interf_p))
      P <- get("prof", interf_p)[[E]]
    }
    else if (tipo == KO_SEMPLICE) {
      K <- gtkSpinButtonGetValueAsInt(get("spn4", interf_p))
      differ <- get("diff", interf_p)$getActive()
      P0 <- get("prof", interf_p)[[1]][[1]]
      P <- get("prof", interf_p)[[K + 1]][[1]]
    }
    else if (tipo == KO_MULTIPLO) {
      E <- gtkSpinButtonGetValueAsInt(get("spn3", interf_p))
      K <- gtkSpinButtonGetValueAsInt(get("spn4", interf_p))
      P0 <- get("prof", interf_p)[[1]][[E]]
      P <- get("prof", interf_p)[[K + 1]][[E]]
      differ <- get("diff", interf_p)$getActive()
    }
  }
  W <- get("weights", interf_p)
  regol <- W[N,]
  n <- ncol(P)
  t <- seq(1, n)
  if (differ)
    P1 <- P - P0
  else
    P1 <- P
  if (get("regol", interf_p)$getActive()){
    y <- prof <- as.numeric(P1[N,])
    regol <- which(regol != 0)
    for (r in regol) {
      cat(sprintf("Gene %d is regulated by gene %d\n", N, r))
      y <- cbind(y, as.numeric(P1[r,]))
    }
    if (length(regol) == 0)
      cat(sprintf("Gene %d has no regulators\n", N))
  }
  else {
    y <- as.numeric(P[N,])
    regol <- vector()
    y <- cbind(y, as.numeric(P1[N,]))
  }
  linee <- rep("l", n)
  par(bg = "black")
  par(col.main = "white")
  par(col.axis = "white")
  par(col.lab = "white")
  plot(t, type="n")
  if (!get("clust", interf_p)$getActive()) {
  y <- as.numeric(y)
  y <- matrix(y, nrow=n)
    line_weights <- c(2, rep(1, n - 1))
    colori <- c("yellow", rainbow(length(regol)))
  matplot(t, y, type=linee, lwd=line_weights, col=colori, xlab="time samples", ylab="expression level", xlim=c(0, n * 6 / 5))
  legend("bottomright", c("Regul.", regol), lwd=line_weights, col=colori, text.col="white")
  }
  else {
    y <- matrix(nrow=n)
    sel_c <- gtkSpinButtonGetValue(get("spn5", interf_p))
    if (get("ricalcola_c", interf_p)) {
      ncl <- gtkSpinButtonGetValue(get("ncl", interf_p))
      ris <- try(cl <- kmeans(P1, ncl))
      if (inherits(ris, "try-error")) {
        message_dialog_p(geterrmessage())
        return()
      }
      assign("ricalcola_c", FALSE, interf_p)
      p <- cl$cluster
      indxNA <- which(is.na(p))
      p[indxNA] <- 1
      un <- length(p[!duplicated(p)])
      gtkSpinButtonSetRange(get("spn5", interf_p), 0, un)
      assign("p", p, interf_p)
      colori <- rainbow(un)
      colori <- sample(colori)
      assign("colori", colori, interf_p)
    }
    else {
      p <- get("p", interf_p)
    }
    p1 <- c()
    for (c in 1:nrow(P1)) {
      if (sel_c > 0) {
        if (p[c] == sel_c) {
          cat(sprintf("Gene %d belongs to cluster %d\n", c, p[c]))
          y <- cbind(y, as.numeric(P1[c,]))
          p1 <- c(p1, p[c])
        }
      }
      else {
        cat(sprintf("Gene %d belongs to cluster %d\n", c, p[c]))
        y <- cbind(y, as.numeric(P1[c,]))
        p1 <- c(p1, p[c])
      }
    }
    cat("---------------------------\n")
    line_weights <- rep(1, length(p1))
    colori <- get("colori", interf_p)
    unici <- p1[!duplicated(p1)]
    colori1 <- colori[p1]
    colori2 <- colori[unici]
    if (length(p1) > 0) {
      matplot(t, y, type=linee, lty=rep(1, length(p1)), lwd=line_weights, col=colori1, xlab="time samples", ylab="expression level", xlim=c(0, n * 6 / 5))
      legend("bottomright", as.character(unici), lwd=line_weights, col=colori2, text.col="white")
    }
  }
  axis(1, col="white")
  axis(2, col="white")
}

onExit <- function(widget, event, userData)
{
  gtkWidgetDestroy(widget)
  setwd("..")
  return (FALSE)
}

netsim_profiles <- function(sn=FALSE)
{
  if (!file.exists(paste(getwd(), "netsim", sep="/")))
    stop("To plot the profiles you must first generate them using 'netsim_simulate()' and save at least a network\n(also 'netsim_generate()', but in that case use 'netsim_profiles(sn=TRUE)')!")
  if (!file.exists("netsim")) {
    stop("There are no suitable profiles to analyze (sub-folder 'netsim' does not exists!)\n")
  }
  else {
  	 setwd("netsim")
    win <- gtkWindow(show=F)
    assign("main", win, interf_p)
    da <- gtkDrawingArea()
    gSignalConnect(win, "delete-event", onExit)
    assign("sn", sn, interf_p)
    mx <- length(i <- grep("weights", dir()))
    pad_w <- ceiling(log10(mx + 1))
    assign("pad_w", pad_w, interf_h)
    W <- read.table(sprintf("weights%0*d.txt", pad_w, 1), row.names=1, header=TRUE)
    assign("weights", W, interf_p)
    # differenziare tra: ko (SIMdata?_ko?.txt), ko ripetuto (SIMdata?_?_ko?.txt) ripetuto (SIMdata?_?.txt) oppure semplice (SIMdata?.txt)
    # se ko c'è anche il wild type, se ripetuto anche l'esperimento e la media, regolatori opzionali
    dirs <- dir()
    base <- "SIMdata_simulateprofiles"
    if(length(grep(paste(base, "\\d+\\.txt", dirs, perl=TRUE)) > 0)
    	tipo = SEMPLICE
    else if (length(grep(paste(base, "\\d+_\\d+\\.txt", dirs, perl=TRUE)) > 0)
    	tipo = MULTIPLO
    else if (length(grep(paste(base, "\\d+_ko\\d+\\.txt", dirs, perl=TRUE)) > 0)
    	tipo = KO_SEMPLICE
    else if (length(grep(paste(base, "\\d+_\\d+_ko\\d+\\.txt", dirs, perl=TRUE)) > 0)
    	tipo = KO_MULTIPLO
    pads <- c(0, 0, 0) # reti, exp, ko
    if (sn) {
      tipo <- 0
      n_exp <- 1
      n_exp1 <- 1
      ko <- ko1 <- vector()
      n_reti <- length(grep("SIMdata\\d+\\.txt", dirs, perl=TRUE))
      pads[1] <- ceiling(log10(n_reti + 1))
      n_prof <- nrow(read.table(sprintf("SIMdata%0*d.txt", pads[1], 1), row.names=1, header=TRUE))
      win$setTitle("Random initial conditions profiles")
    }
    else if (tipo == KO_MULTIPLO) {
      n_reti <- length(grep("SIMdata_simulateprofiles\\d+_[0]*1_ko[0]+\\.txt", dirs, perl=TRUE))
      pads[1] <- ceiling(log10(n_reti + 1))
		n_exp <- length(grep("SIMdata_simulateprofiles[0]+1_\\d+_ko[0]+\\.txt", dirs, perl=TRUE))      
		pads[2] <- ceiling(log10(n_exp + 1))
		n_ko <- length(grep("SIMdata_simulateprofiles[0]+1_[0]*1_ko\\d+\\.txt", dirs, perl=TRUE))      
		pads[3] <- ceiling(log10(n_ko + 1))
      n_prof <- nrow(read.table(sprintf("SIMdata_simulateprofiles%0*d_%0*d_ko*0*d.txt", pads[1], 1, pads[2], 1, pads[3], 0), row.names=1, header=TRUE))
      pads[4] <- ceiling(log10(n_prof + 1))
      ko <- gsub(sprintf("SIMdata_simulateprofiles[0]*1_[0]*1_ko(\\d+)\\.txt", pads[1], 1), "\\1", ko1, perl=TRUE)
      win$setTitle("Multiple knock-out experiments")
    }
    else if (tipo == KO_SEMPLICE) {
      n_reti <- length(grep("SIMdata_simulateprofiles\\d+_ko[0]+\\.txt", dirs, perl=TRUE))
      pads[1] <- ceiling(log10(n_reti + 1))
		pads[2] <- 0
		n_ko <- length(grep("SIMdata_simulateprofiles[0]+1_ko\\d+\\.txt", dirs, perl=TRUE))      
		pads[3] <- ceiling(log10(n_ko + 1))
      n_prof <- nrow(read.table(sprintf("SIMdata_simulateprofiles%0*d_ko*0*d.txt", pads[1], 1, pads[3], 0), row.names=1, header=TRUE))
      pads[4] <- ceiling(log10(n_prof + 1))
      ko <- gsub(sprintf("SIMdata_simulateprofiles[0]*1_ko(\\d+)\\.txt", pads[1], 1), "\\1", ko1, perl=TRUE)
      win$setTitle("Single knock-out experiments")
    }
    else if (tipo == MULTIPLO) {
      n_reti <- length(grep("SIMdata_simulateprofiles\\d+_[0]*1\\.txt", dirs, perl=TRUE))
      pads[1] <- ceiling(log10(n_reti + 1))
		n_exp <- length(grep("SIMdata_simulateprofiles[0]+1_\\d+_\\.txt", dirs, perl=TRUE))      
		pads[2] <- ceiling(log10(n_exp + 1))  
		pads[3] <- 0
      n_prof <- nrow(read.table(sprintf("SIMdata_simulateprofiles%0*d_%0*d.txt", pads[1], 1, pads[2], 1), row.names=1, header=TRUE))
      pads[4] <- ceiling(log10(n_prof + 1))
      ko <- ko1 <- vector()
      win$setTitle("Multiple experiments")
    }
    else if (tipo == SEMPLICE) {
      n_reti <- length(grep("SIMdata_simulateprofiles\\d+\\.txt", dirs, perl=TRUE))
      pads[1] <- ceiling(log10(n_reti + 1))
		pads[2] <- 0   
		pads[3] <- 0
      n_prof <- nrow(read.table(sprintf("SIMdata_simulateprofiles%0*d.txt", pads[1], 1), row.names=1, header=TRUE))
      pads[4] <- ceiling(log10(n_prof + 1))
      ko <- ko1 <- vector()
      win$setTitle("Single experiments")
    }
    else {
      setwd("..")
      stop("There are no suitable profiles to analyze!\n")
    }
    assign("tipo", tipo, interf_p)
    assign("ricalcola_c", TRUE, interf_p)
    assign("pads", pads, interf_p)
    label1 <- gtkLabelNew("Profile: ")
    label2 <- gtkLabelNew("Network: ")
    hbox <- gtkHBox()
    spn1 <- gtkSpinButtonNewWithRange(1, n_prof, 1)
    assign("spn1", spn1, interf_p)
    spn1$SetWidthChars(5)
    spn2 <- gtkSpinButtonNewWithRange(1, n_reti, 1)
    assign("spn2", spn2, interf_p)
    spn2$SetWidthChars(5)
    vbox <- gtkVBox(FALSE, 5)
    vbox$packStart(da, FALSE, FALSE, 5)
    hbox <- gtkHBox()
    hbox$packStart(label1, FALSE, FALSE, 5)
    hbox$packStart(spn1, FALSE, FALSE, 5)
    hbox$packStart(label2, FALSE, FALSE, 5)
    hbox$packStart(spn2, FALSE, FALSE, 5)
    clust <- gtkCheckButtonNewWithLabel("Clusters:")
    assign("clust", clust, interf_p)
    ncl <- gtkSpinButtonNewWithRange(2, n_prof - 1, 1)
    ncl$setValue(3)
    gSignalConnect(ncl, "value-changed", ncl_changed_p)
    regol <- gtkCheckButtonNewWithLabel("regulators")
    regol$name <- "regol"
    assign("regol", regol, interf_p)
    if (tipo == MULTIPLO || tipo == KO_MULTIPLO) {
      lbl <- gtkLabelNew("experiment")
      hbox$packStart(lbl, FALSE, FALSE, 5)
      n <- n_exp
      spn3 <- gtkSpinButtonNewWithRange(1, n, 1)
      hbox$packStart(spn3, FALSE, FALSE, 5)
      gSignalConnect(spn3, "value-changed", spin_changed_p)
      assign("spn3", spn3, interf_p)
    }
    if (tipo == KO_SEMPLICE || tipo == KO_MULTIPLO) {
      lbl <- gtkLabelNew("knocked gene")
      n <- length(ko) - 1
      spn4 <- gtkSpinButtonNewWithRange(0, n, 1)
      spn4$name <- "ko"
      gn <- gtkLabelNew("(none)")
      hbox$packStart(lbl, FALSE, FALSE, 5)
      hbox$packStart(spn4, FALSE, FALSE, 5)
      gSignalConnect(spn4, "value-changed", spin_changed_p, ko)
      assign("spn4", spn4, interf_p)
      hbox$packStart(gn, FALSE, FALSE, 5)
      assign("ko", ko, interf_p)
      assign("gn", gn, interf_p)
    }
    vbox$packStart(hbox, FALSE, FALSE, 5)
    hbox1 <- gtkHBox()
    gSignalConnect(clust, "toggled", chk_toggled, spn1)
    clust$name <- "clust"
    hbox1$packStart(clust, FALSE, FALSE, 5)
    hbox1$packStart(ncl, FALSE, FALSE, 5)
    assign("ncl", ncl, interf_p)
    lcl <- gtkLabelNew("cluster:")
    hbox1$packStart(lcl, FALSE, FALSE, 5)
    spn5 <- gtkSpinButtonNewWithRange(0, 0, 1)
    hbox1$packStart(spn5, FALSE, FALSE, 5)
    gSignalConnect(spn5, "value-changed", spin_changed_p)
    assign("spn5", spn5, interf_p)
    gSignalConnect(regol, "toggled", chk_toggled, spn1)
    hbox1$packStart(regol, FALSE, FALSE, 5)
    if (tipo == MULTIPLO || tipo == KO_MULTIPL) {
      media <- gtkCheckButtonNewWithLabel("mean")
      gSignalConnect(media, "toggled", chk_toggled, spn1)
      hbox1$packStart(media, FALSE, FALSE, 5)
      assign("media", media, interf_p)
    }
    if (tipo == KO_SEMPLICE || tipo == KO_MULTIPLO) {
      differ <- gtkCheckButtonNewWithLabel("difference w.r.t. WT")
      gSignalConnect(differ, "toggled", chk_toggled, spn1)
      hbox1$packStart(differ, FALSE, FALSE, 5)
      assign("diff", differ, interf_p)
    }
    vbox$packStart(hbox1, FALSE, FALSE, 5)
    win$add(vbox)
    gSignalConnect(spn1, "value-changed", spin_changed_p)
    gSignalConnect(spn2, "value-changed", rete_changed_p, list(n_exp, n_exp1, ko))
    win$setDefaultSize(600,400)
    gtkWidgetSetSizeRequest(da, 400, 400)
    win$showAll()
    require(cairoDevice)
    asCairoDevice(da)
    rete_changed_p(spn1, list(n_exp, n_exp1, ko))
  }
}

