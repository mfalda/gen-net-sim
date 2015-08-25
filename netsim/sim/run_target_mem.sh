R -d "valgrind -v --xml=yes --log-file=errori.xml --tool=memcheck --leak-check=full" --vanilla < $1

valkyrie --view-log=errori.xml