set xlabel "Frequency (Hz)"
set ylabel "Function value"

plot "function_data.fout" u 1:2 title "Re{f}" w p ,\
     "function_data.fout" u 1:3 title "Im{f}" w p ,\
     "function_data.fout" u 1:4 title "|f|" w p ,\
     "function_data.fout" u 1:5 title "Re{f_fit}" w l ,\
     "function_data.fout" u 1:6 title "Im{f_fit}" w l,\
     "function_data.fout" u 1:7 title "|f_fit|" w l
pause -1

set autoscale x
set autoscale y
set logscale x
plot "function_data.fout" u 1:2 title "Re{f}" w p ,\
     "function_data.fout" u 1:3 title "Im{f}" w p ,\
     "function_data.fout" u 1:4 title "|f|" w p ,\
     "function_data.fout" u 1:5 title "Re{f_fit}" w l ,\
     "function_data.fout" u 1:6 title "Im{f_fit}" w l,\
     "function_data.fout" u 1:7 title "|f_fit|" w l
pause -1

set autoscale x
set autoscale y
set logscale x
set logscale y
plot  "function_data.fout" u 1:4 title "|f|" w p ,\
     "function_data.fout" u 1:7 title "|f_fit|" w l
pause -1

set autoscale x
set autoscale y
set nologscale x
set logscale y
plot  "function_data.fout" u 1:4 title "|f|" w p ,\
     "function_data.fout" u 1:7 title "|f_fit|" w l
pause -1

     
