sudo opcontrol --reset
sudo opcontrol --no-vmlinux
sudo opcontrol --start
R -f prova_HMM.R
sudo opcontrol --shutdown
opreport --demangle=smart --long-filenames --symbols /usr/lib/R/bin/exec/R  --debug-info > out.txt
