init <- function()
{
	dir.create("distrib")
	setwd("distrib")
	scrivi_prob()
}

scrivi_prob <- function()
{
	id <<- 0
	salva <<- T
	unlink("sample*.txt")
	unlink("normali*.txt")
	unlink("uniformi*.txt")
	unlink("log-normali*.txt")
	fp <<- file("moduli_R.txt", "w")
}

leggi_prob <- function()
{
	id <<- 0
	salva <<- F
	fp <<- file("moduli_R.txt", "w")
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
