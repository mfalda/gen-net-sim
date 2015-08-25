rnorm_s <- function(n, mean=0, sd=1, chi)
{
	id <<- id + 1
	nome <- sprintf("normali_%d.txt", id)
	if (exists("salva") && get("salva") == T) {
		if (!file.exists(nome)) {
			stop("Il file '" + nome + "' non esiste!\n")
		}
		v <- scan(nome)
		cat(file=fp, sprintf("+++Letti %d valore/i per '%s' da '%s'\n", length(v), chi, nome))
	}
	else {
		v <- rnorm(n, mean, sd)
		write(file=nome, v, ncolumns=length(v))
	}
	return (v)
}

rlnorm_s <- function(n, mean=0, sd=1, chi)
{
	id <<- id + 1
	nome <- sprintf("log-normali_%d.txt", id)
	if (exists("salva") && get("salva") == T) {
		if (!file.exists(nome)) {
			stop("Il file '" + nome + "' non esiste!\n")
		}
		v <- scan(nome)
		cat(file=fp, sprintf("+++Letti %d valore/i per '%s' da '%s'\n", length(v), chi, nome))
	}
	else {
		v <- rlnorm(n, mean, sd)
		write(file=nome, v, ncolumns=length(v))
	}
	return (v)
}

runif_s <- function(n, min=0, max=1, chi)
{
	id <<- id + 1
	nome <- sprintf("unif_%d.txt", id)
	if (exists("salva") && get("salva") == T) {
		if (!file.exists(nome)) {
			stop("Il file '" + nome + "' non esiste!\n")
		}
		v <- scan(nome)
		cat(file=fp, sprintf("+++Letti %d valore/i per '%s' da '%s'\n", length(v), chi, nome))
	}
	else {
		v <- runif(n, min, max)
		write(file=nome, v, ncolumns=length(v))
	}
	return (v)
}

sample_s <- function(x, size, replace=FALSE, prob=NULL, chi="sampleB")
{
	id <<- id + 1
	nome <- sprintf("sample_%d", id)
	if (exists("salva") && get("salva") == T) {
		if (!file.exists(nome, ".txt")) {
			if (!file.exists(nome, "-vuoto.txt")) {
				stop("Il file n. " + id + "per sample non esiste!\n")
			}
			cat(file=fp, sprintf("+++Vettore nullo per '%s' da '%s'\n", chi, nome))
			ris <- integer(0)
		}
		ris <- scan(nome)
		cat(file=fp, sprintf("+++Letti %d valore/i per '%s' da '%s'\n", length(ris), chi, nome))
	}
	else {
		ris <- sample(x, size, replace, prob)
		if (length(ris) == 0) {
			nome <- sprintf("%s-vuoto.txt", nome)
			write(file=nome, ris)
		}
		else {
			nome <- sprintf("%s.txt", nome)
			write(file=nome, ris, ncolumns=length(ris))
		}
	}
	return (ris)
}
