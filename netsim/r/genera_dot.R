expose_event_callback <- function(widget, event, data) {
  cr <- gdkCairoCreate(widget[["window"]])
image = cairoImageSurfaceCreateFromPng("eg.png")
w = widget[["allocation"]]$width
h = widget[["allocation"]]$height
w1 <- cairoImageSurfaceGetWidth(image)
h1 <- cairoImageSurfaceGetHeight(image)
#cairoScale(cr, w, h)

cairoSetSourceSurface(cr, image, 0, 0)
cairoPaint(cr)
  return(TRUE)
}

f <- file("eg.dot", "w")
cat(file=f, "digraph G {\n\tviewport=\"200,200,1,100,100\"\n\t1 -> 2;\n\t2 -> 3;\n\t1 -> 3;\n}")
close(f)
system(sprintf("dot -Tpng -oeg.png eg.dot"))
require(RGtk2)
win <- gtkWindow()
da <- gtkDrawingAreaNew()
da$setSizeRequest(200, 200)
gSignalConnect(da, "expose_event", expose_event_callback)
win$add(da)