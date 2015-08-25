init <- function()
{
	id <<- 0
	dir.create("distrib")
	salva <- FALSE
	setwd("distrib")
	fp <<- file("moduli.txt", "w")
}

scrivi_prob <- function()
{
	id <<- 0
	salva <<- FALSE
	close(fp)
	unlink("sample*.txt")
	unlink("normali*.txt")
	unlink("uniformi*.txt")
	unlink("log-normali*.txt")
	fp <<- file("moduli.txt", "w")
}

leggi_prob <- function()
{
	id <<- 0
	salva <<- TRUE
	close(fp)
	fp <<- file("moduli.txt", "w")
}

chiudi <- function()
{
	close(fp)
}

print_list <-function(l)
{
	for (i in l)
		cat(file=fp, "\t", i, "\n")
	cat(file=fp, "\n")
}