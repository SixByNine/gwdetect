set log
set term png size 800,600 large
set out "J0437-4715.png" ; set title "J0437-4715"
plot 'J0437-4715-J0613-0200.ss' u 1:3 w lp lt 1, 'J0437-4715-J0613-0200.ss' u 1:5 w l lt 2, 'J0437-4715-J0613-0200.ss' u 1:7 w l lt 3, 'J0437-4715-J0613-0200.ss' u 1:($5+$7) w l lt 4

set out "J0613-0200.png" ; set title "J0613-0200"
plot 'J0437-4715-J0613-0200.ss' u 1:4 w lp lt 1, 'J0437-4715-J0613-0200.ss' u 1:6 w l lt 2, 'J0437-4715-J0613-0200.ss' u 1:7 w l lt 3, 'J0437-4715-J0613-0200.ss' u 1:($6+$7) w l lt 4
set out "J0711-6830.png" ; set title "J0711-6830"
plot 'J0437-4715-J0711-6830.ss' u 1:4 w lp lt 1, 'J0437-4715-J0711-6830.ss' u 1:6 w l lt 2, 'J0437-4715-J0711-6830.ss' u 1:7 w l lt 3, 'J0437-4715-J0711-6830.ss' u 1:($6+$7) w l lt 4
set out "J1022+1001.png" ; set title "J1022+1001"
plot 'J0437-4715-J1022+1001.ss' u 1:4 w lp lt 1, 'J0437-4715-J1022+1001.ss' u 1:6 w l lt 2, 'J0437-4715-J1022+1001.ss' u 1:7 w l lt 3, 'J0437-4715-J1022+1001.ss' u 1:($6+$7) w l lt 4
set out "J1024-0719.png" ; set title "J1024-0719"
plot 'J0437-4715-J1024-0719.ss' u 1:4 w lp lt 1, 'J0437-4715-J1024-0719.ss' u 1:6 w l lt 2, 'J0437-4715-J1024-0719.ss' u 1:7 w l lt 3, 'J0437-4715-J1024-0719.ss' u 1:($6+$7) w l lt 4
set out "J1045-4509.png" ; set title "J1045-4509"
plot 'J0437-4715-J1045-4509.ss' u 1:4 w lp lt 1, 'J0437-4715-J1045-4509.ss' u 1:6 w l lt 2, 'J0437-4715-J1045-4509.ss' u 1:7 w l lt 3, 'J0437-4715-J1045-4509.ss' u 1:($6+$7) w l lt 4
set out "J1600-3053.png" ; set title "J1600-3053"
plot 'J0437-4715-J1600-3053.ss' u 1:4 w lp lt 1, 'J0437-4715-J1600-3053.ss' u 1:6 w l lt 2, 'J0437-4715-J1600-3053.ss' u 1:7 w l lt 3, 'J0437-4715-J1600-3053.ss' u 1:($6+$7) w l lt 4
set out "J1603-7202.png" ; set title "J1603-7202"
plot 'J0437-4715-J1603-7202.ss' u 1:4 w lp lt 1, 'J0437-4715-J1603-7202.ss' u 1:6 w l lt 2, 'J0437-4715-J1603-7202.ss' u 1:7 w l lt 3, 'J0437-4715-J1603-7202.ss' u 1:($6+$7) w l lt 4
set out "J1643-1224.png" ; set title "J1643-1224"
plot 'J0437-4715-J1643-1224.ss' u 1:4 w lp lt 1, 'J0437-4715-J1643-1224.ss' u 1:6 w l lt 2, 'J0437-4715-J1643-1224.ss' u 1:7 w l lt 3, 'J0437-4715-J1643-1224.ss' u 1:($6+$7) w l lt 4
set out "J1713+0747.png" ; set title "J1713+0747"
plot 'J0437-4715-J1713+0747.ss' u 1:4 w lp lt 1, 'J0437-4715-J1713+0747.ss' u 1:6 w l lt 2, 'J0437-4715-J1713+0747.ss' u 1:7 w l lt 3, 'J0437-4715-J1713+0747.ss' u 1:($6+$7) w l lt 4
set out "J1730-2304.png" ; set title "J1730-2304"
plot 'J0437-4715-J1730-2304.ss' u 1:4 w lp lt 1, 'J0437-4715-J1730-2304.ss' u 1:6 w l lt 2, 'J0437-4715-J1730-2304.ss' u 1:7 w l lt 3, 'J0437-4715-J1730-2304.ss' u 1:($6+$7) w l lt 4
set out "J1732-5049.png" ; set title "J1732-5049"
plot 'J0437-4715-J1732-5049.ss' u 1:4 w lp lt 1, 'J0437-4715-J1732-5049.ss' u 1:6 w l lt 2, 'J0437-4715-J1732-5049.ss' u 1:7 w l lt 3, 'J0437-4715-J1732-5049.ss' u 1:($6+$7) w l lt 4
set out "J1744-1134.png" ; set title "J1744-1134"
plot 'J0437-4715-J1744-1134.ss' u 1:4 w lp lt 1, 'J0437-4715-J1744-1134.ss' u 1:6 w l lt 2, 'J0437-4715-J1744-1134.ss' u 1:7 w l lt 3, 'J0437-4715-J1744-1134.ss' u 1:($6+$7) w l lt 4
set out "J1857+0943.png" ; set title "J1857+0943"
plot 'J0437-4715-J1857+0943.ss' u 1:4 w lp lt 1, 'J0437-4715-J1857+0943.ss' u 1:6 w l lt 2, 'J0437-4715-J1857+0943.ss' u 1:7 w l lt 3, 'J0437-4715-J1857+0943.ss' u 1:($6+$7) w l lt 4
set out "J1909-3744.png" ; set title "J1909-3744"
plot 'J0437-4715-J1909-3744.ss' u 1:4 w lp lt 1, 'J0437-4715-J1909-3744.ss' u 1:6 w l lt 2, 'J0437-4715-J1909-3744.ss' u 1:7 w l lt 3, 'J0437-4715-J1909-3744.ss' u 1:($6+$7) w l lt 4
set out "J2124-3358.png" ; set title "J2124-3358"
plot 'J0437-4715-J2124-3358.ss' u 1:4 w lp lt 1, 'J0437-4715-J2124-3358.ss' u 1:6 w l lt 2, 'J0437-4715-J2124-3358.ss' u 1:7 w l lt 3, 'J0437-4715-J2124-3358.ss' u 1:($6+$7) w l lt 4
set out "J2129-5721.png" ; set title "J2129-5721"
plot 'J0437-4715-J2129-5721.ss' u 1:4 w lp lt 1, 'J0437-4715-J2129-5721.ss' u 1:6 w l lt 2, 'J0437-4715-J2129-5721.ss' u 1:7 w l lt 3, 'J0437-4715-J2129-5721.ss' u 1:($6+$7) w l lt 4
set out "J2145-0750.png" ; set title "J2145-0750"
plot 'J0437-4715-J2145-0750.ss' u 1:4 w lp lt 1, 'J0437-4715-J2145-0750.ss' u 1:6 w l lt 2, 'J0437-4715-J2145-0750.ss' u 1:7 w l lt 3, 'J0437-4715-J2145-0750.ss' u 1:($6+$7) w l lt 4

set out
