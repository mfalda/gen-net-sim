setwd("/home/marco/R2Clib")

dyn.load("write_table.so")

wtable_int <-
    function(nome, m)
{
    stopifnot(is.character(nome))
    stopifnot(is.matrix(m))
    vett <- .Call("write_m_int",
                 as.character(nome),
                 as.matrix(m)
					)
}

wtable_double <-
    function(nome, m)
{
    stopifnot(is.character(nome))
    stopifnot(is.matrix(m))
    vett <- .Call("write_m_double",
                 as.character(nome),
                 as.matrix(m)
					)
}

wtable_param_double <-
    function(nome, m, nomi)
{
    stopifnot(is.character(nome))
    stopifnot(is.matrix(m))
    stopifnot(is.list(nomi))
    vett <- .Call("write_param_double",
                 as.character(nome),
                 as.matrix(m),
                 as.list(nomi)
					)
}

m <- matrix(seq(0, 1 - 1 / 6, 1/6), 2, 3)
wtable_double("matrice_d.txt", m)

m <- matrix(seq(1,6, 1), 2, 3)
wtable_int("matrice_i.txt", m)

m <- matrix(seq(0, 1 - 1 / 6, 1/6), 2, 3)
wtable_param_double("matrice_p.txt", m, list("a", "b", "c"))