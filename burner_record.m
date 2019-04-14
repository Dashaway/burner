%190414
%Dashaway
%数据计算和保存


clear;
close all;
%写入excel

%算例1
%装药参数
D = 0.300;         %药柱外直径(m)
r = 0.010;        %微圆弧半径(m)
l = 0.130;        %弧心距(m)
R = 0.01;       %中心孔(m)
m_s = 8;        %弧孔药比
n_s = 3;        %弧数量
Lp = 0.500;      %药柱长度(m)
burner_DA_calculate;
close all;

FileName = 'data_burner.xlsx';
SheetName = 'sheet1';

Str = '算例1';
xlswrite(FileName,{Str},SheetName,'C3');
Str = 't_max';
xlswrite(FileName,{Str},SheetName,'C4');
xlswrite(FileName,t_max,SheetName,'C5');
Str = 'F_max';
xlswrite(FileName,{Str},SheetName,'C6');
xlswrite(FileName,F_max,SheetName,'C7');
Str = 'F_min';
xlswrite(FileName,{Str},SheetName,'C8');
xlswrite(FileName,F_min,SheetName,'C9');
Str = 'i';
xlswrite(FileName,{Str},SheetName,'C10');
xlswrite(FileName,i,SheetName,'C11');
Str = 'F';
xlswrite(FileName,{Str},SheetName,'D3');
xlswrite(FileName,F',SheetName,'D4');

clear;

%算例2
%装药参数
D = 0.300;         %药柱外直径(m)
r = 0.010;        %微圆弧半径(m)
l = 0.120;        %弧心距(m)
R = 0.01;       %中心孔(m)
m_s = 8;        %弧孔药比
n_s = 3;        %弧数量
Lp = 0.500;      %药柱长度(m)
burner_DA_calculate;
close all;

FileName = 'data_burner.xlsx';
SheetName = 'sheet2';

Str = '算例2';
xlswrite(FileName,{Str},SheetName,'C3');
Str = 't_max';
xlswrite(FileName,{Str},SheetName,'C4');
xlswrite(FileName,t_max,SheetName,'C5');
Str = 'F_max';
xlswrite(FileName,{Str},SheetName,'C6');
xlswrite(FileName,F_max,SheetName,'C7');
Str = 'F_min';
xlswrite(FileName,{Str},SheetName,'C8');
xlswrite(FileName,F_min,SheetName,'C9');
Str = 'i';
xlswrite(FileName,{Str},SheetName,'C10');
xlswrite(FileName,i,SheetName,'C11');
Str = 'F';
xlswrite(FileName,{Str},SheetName,'D3');
xlswrite(FileName,F',SheetName,'D4');

clear;

%算例3
%装药参数
D = 0.300;         %药柱外直径(m)
r = 0.010;        %微圆弧半径(m)
l = 0.130;        %弧心距(m)
R = 0.05;       %中心孔(m)
m_s = 8;        %弧孔药比
n_s = 3;        %弧数量
Lp = 0.500;      %药柱长度(m)
burner_DA_calculate;
close all;

FileName = 'data_burner.xlsx';
SheetName = 'sheet3';

Str = '算例3';
xlswrite(FileName,{Str},SheetName,'C3');
Str = 't_max';
xlswrite(FileName,{Str},SheetName,'C4');
xlswrite(FileName,t_max,SheetName,'C5');
Str = 'F_max';
xlswrite(FileName,{Str},SheetName,'C6');
xlswrite(FileName,F_max,SheetName,'C7');
Str = 'F_min';
xlswrite(FileName,{Str},SheetName,'C8');
xlswrite(FileName,F_min,SheetName,'C9');
Str = 'i';
xlswrite(FileName,{Str},SheetName,'C10');
xlswrite(FileName,i,SheetName,'C11');
Str = 'F';
xlswrite(FileName,{Str},SheetName,'D3');
xlswrite(FileName,F',SheetName,'D4');

clear;

%算例4
%装药参数
D = 0.300;         %药柱外直径(m)
r = 0.005;        %微圆弧半径(m)
l = 0.130;        %弧心距(m)
R = 0.01;       %中心孔(m)
m_s = 8;        %弧孔药比
n_s = 3;        %弧数量
Lp = 0.500;      %药柱长度(m)
burner_DA_calculate;
close all;

FileName = 'data_burner.xlsx';
SheetName = 'sheet4';

Str = '算例4';
xlswrite(FileName,{Str},SheetName,'C3');
Str = 't_max';
xlswrite(FileName,{Str},SheetName,'C4');
xlswrite(FileName,t_max,SheetName,'C5');
Str = 'F_max';
xlswrite(FileName,{Str},SheetName,'C6');
xlswrite(FileName,F_max,SheetName,'C7');
Str = 'F_min';
xlswrite(FileName,{Str},SheetName,'C8');
xlswrite(FileName,F_min,SheetName,'C9');
Str = 'i';
xlswrite(FileName,{Str},SheetName,'C10');
xlswrite(FileName,i,SheetName,'C11');
Str = 'F';
xlswrite(FileName,{Str},SheetName,'D3');
xlswrite(FileName,F',SheetName,'D4');

clear;

%算例5
%装药参数
D = 0.300;         %药柱外直径(m)
r = 0.010;        %微圆弧半径(m)
l = 0.130;        %弧心距(m)
R = 0.01;       %中心孔(m)
m_s = 6;        %弧孔药比
n_s = 3;        %弧数量
Lp = 0.500;      %药柱长度(m)
burner_DA_calculate;
close all;

FileName = 'data_burner.xlsx';
SheetName = 'sheet5';

Str = '算例5';
xlswrite(FileName,{Str},SheetName,'C3');
Str = 't_max';
xlswrite(FileName,{Str},SheetName,'C4');
xlswrite(FileName,t_max,SheetName,'C5');
Str = 'F_max';
xlswrite(FileName,{Str},SheetName,'C6');
xlswrite(FileName,F_max,SheetName,'C7');
Str = 'F_min';
xlswrite(FileName,{Str},SheetName,'C8');
xlswrite(FileName,F_min,SheetName,'C9');
Str = 'i';
xlswrite(FileName,{Str},SheetName,'C10');
xlswrite(FileName,i,SheetName,'C11');
Str = 'F';
xlswrite(FileName,{Str},SheetName,'D3');
xlswrite(FileName,F',SheetName,'D4');

clear;

%算例6
%装药参数
D = 0.300;         %药柱外直径(m)
r = 0.010;        %微圆弧半径(m)
l = 0.130;        %弧心距(m)
R = 0.01;       %中心孔(m)
m_s = 6;        %弧孔药比
n_s = 5;        %弧数量
Lp = 0.500;      %药柱长度(m)
burner_DA_calculate;
close all;

FileName = 'data_burner.xlsx';
SheetName = 'sheet6';

Str = '算例6';
xlswrite(FileName,{Str},SheetName,'C3');
Str = 't_max';
xlswrite(FileName,{Str},SheetName,'C4');
xlswrite(FileName,t_max,SheetName,'C5');
Str = 'F_max';
xlswrite(FileName,{Str},SheetName,'C6');
xlswrite(FileName,F_max,SheetName,'C7');
Str = 'F_min';
xlswrite(FileName,{Str},SheetName,'C8');
xlswrite(FileName,F_min,SheetName,'C9');
Str = 'i';
xlswrite(FileName,{Str},SheetName,'C10');
xlswrite(FileName,i,SheetName,'C11');
Str = 'F';
xlswrite(FileName,{Str},SheetName,'D3');
xlswrite(FileName,F',SheetName,'D4');

clear;

disp('done');
