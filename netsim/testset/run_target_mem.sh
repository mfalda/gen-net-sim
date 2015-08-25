R -d "valgrind -v  --xml=yes --log-file-exactly=errori.xml --tool=memcheck --leak-check=full" --vanilla < target.R

valkyrie --view-log=errori.xml