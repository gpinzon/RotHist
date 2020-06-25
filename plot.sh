echo " "

echo "   PMS Rotational Model     "
echo " "
date
echo " "
rm fome
rm last*
rm output*


g++ history_v15_7Mo.c
./a.out

xterm +ah -e sm inp OuT7.0
ps2pdf model7Mo.ps

mv output output7.0



g++ history_v15_6Mo.c
./a.out

xterm +ah -e sm inp OuT6.0
ps2pdf model6Mo.ps

mv output output6.0



g++ history_v15_5Mo.c
./a.out

xterm +ah -e sm inp OuT5.0
ps2pdf model5Mo.ps

mv output output5.0





g++ history_v15_4Mo.c
./a.out

xterm +ah -e sm inp OuT4.0
ps2pdf model4Mo.ps

mv output output4.0



g++ history_v15_35Mo.c
./a.out

xterm +ah -e sm inp OuT3.5
ps2pdf model3.5Mo.ps

mv output output3.5


g++ history_v15_3Mo.c
./a.out

xterm +ah -e sm inp OuT3.0
ps2pdf model3Mo.ps

mv output output3.0


g++ history_v15_27Mo.c
./a.out

xterm +ah -e sm inp OuT2.7
ps2pdf model2.7Mo.ps

mv output output2.7


g++ history_v15_25Mo.c
./a.out

xterm +ah -e sm inp OuT2.5
ps2pdf model2.5Mo.ps

mv output output2.5


g++ history_v15_22Mo.c
./a.out

xterm +ah -e sm inp OuT2.2
ps2pdf model2.2Mo.ps

mv output output2.2


g++ history_v15_2Mo.c
./a.out

xterm +ah -e sm inp OuT2.0
ps2pdf model2Mo.ps

mv output output2.0


g++ history_v15_19Mo.c
./a.out

xterm +ah -e sm inp OuT1.9
ps2pdf model1.9Mo.ps

mv output output1.9


g++ history_v15_18Mo.c
./a.out

xterm +ah -e sm inp OuT1.8
ps2pdf model1.8Mo.ps

mv output output1.8

g++ history_v15_17Mo.c
./a.out

xterm +ah -e sm inp OuT1.7
ps2pdf model1.7Mo.ps

mv output output1.7


g++ history_v15_16Mo.c
./a.out

xterm +ah -e sm inp OuT1.6
ps2pdf model1.6Mo.ps

mv output output1.6


g++ history_v15_15Mo.c
./a.out

xterm +ah -e sm inp OuT1.5
ps2pdf model1.5Mo.ps

mv output output1.5


g++ history_v15_14Mo.c
./a.out

 xterm +ah -e sm inp OuT1.4
 ps2pdf model1.4Mo.ps

 mv output output1.4


 g++ history_v15_13Mo.c
 ./a.out

 xterm +ah -e sm inp OuT1.3
 ps2pdf model1.3Mo.ps

 mv output output1.3


 g++ history_v15_12Mo.c
 ./a.out

 xterm +ah -e sm inp OuT1.2
 ps2pdf model1.2Mo.ps

 mv output output1.2

 g++ history_v15_11Mo.c
 ./a.out

 xterm +ah -e sm inp OuT1.1
 ps2pdf model1.1Mo.ps

 mv output output1.1


 g++ history_v15_1Mo.c
 ./a.out

 xterm +ah -e sm inp OuT1
 ps2pdf model1Mo.ps

 mv output output1.0

 g++ history_v15_08Mo.c
 ./a.out

 xterm +ah -e sm inp OuT0.8
 ps2pdf model0.8Mo.ps


 mv output output0.8


 g++ history_v15_06Mo.c
 ./a.out

 xterm +ah -e sm inp OuT0.6
 ps2pdf model0.6Mo.ps


 mv output output0.6



 g++ history_v15_05Mo.c
 ./a.out



 mv output output0.5





 g++ history_v15_04Mo.c
 ./a.out

 xterm +ah -e sm inp OuT0.4
 ps2pdf model0.4Mo.ps


 mv output output0.4


 g++ history_v15_03Mo.c
 ./a.out



 mv output output0.3




 g++ history_v15_02Mo.c
 ./a.out

 xterm +ah -e sm inp OuT0.2
 ps2pdf model0.2Mo.ps


 mv output output0.2


 g++ history_v15_013Mo.c
 ./a.out

 xterm +ah -e sm inp OuT0.13
 ps2pdf model0.13Mo.ps

 mv output output0.13

tail -1 output7.0 > last7.0
tail -1 output6.0 > last6.0
tail -1 output5.0 > last5.0
tail -1 output4.0 > last4.0
tail -1 output3.5 > last3.5
tail -1 output3.0 > last3.0
tail -1 output2.7 > last2.7
tail -1 output2.5 > last2.5
tail -1 output2.2 > last2.2
tail -1 output2.0 > last2.0
tail -1 output1.9 > last1.9
tail -1 output1.8 > last1.8
tail -1 output1.7 > last1.7
tail -1 output1.6 > last1.6
tail -1 output1.5 > last1.5
tail -1 output1.4 > last1.4
tail -1 output1.3 > last1.3
tail -1 output1.2 > last1.2
tail -1 output1.1 > last1.1
tail -1 output1.0 > last1.0
tail -1 output0.8 > last0.8
tail -1 output0.6 > last0.6
tail -1 output0.5 > last0.5
tail -1 output0.4 > last0.4
tail -1 output0.3 > last0.3
tail -1 output0.2 > last0.2
tail -1 output0.13 > last0.13

cat last* > fome
echo " "
echo "        "
date
echo " "










