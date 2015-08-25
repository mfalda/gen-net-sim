escape <- function(pos, julia, n)
{
  z <- pos + julia  
  c <- z
	for (i in 0:n) {
	  z <- z ^ 2 + c
		if(Mod(z) > 2.0)  
			return (i)        
	}
	return (0)
}

min_x <- 0.32
max_x <- 0.38
min_y <- -0.38
max_y <- -0.32
risol <- 50
julia <- 0 + 0i
n <- 50
j <- 1

passo_x <- (max_x - min_x) / risol
passo_y <- (max_y - min_y) / risol
m <- matrix(ncol=risol + 1, nrow=risol + 1) 
for (y in seq(min_y, max_y, passo_y)) {
  i <- 1
  for (x in seq(min_x, max_x, passo_x)) {
    m[i, j] <- x + y * 1i
    i <- i + 1
  }
  j <- j + 1
}             
cat("Calcolo la matrice ... ")
t <- system.time(m1 <- apply(m, 1:2, escape, julia, n))
cat(sprintf("fatto: (%3.3f secondi)\n", t[1][[1]]))
image(seq(min_x, max_x, passo_x), seq(min_y, max_y, passo_y), m1, xlab="Real", ylab="Imag", col=c("black", rainbow(n)))                                                                   