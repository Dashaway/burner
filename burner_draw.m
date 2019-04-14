%190414
%Dashaway
%数据读取保存及绘图



close all;
%写入excel

%文件1
burner_DA_wheel;
close all;

FileName = 'data_burner.xlsx';
SheetName = 'sheet1';

Str = 'burner_DA_wheel';
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

%文件2
burner_DA_mutihole;
close all;

FileName = 'data_burner.xlsx';
SheetName = 'sheet1';

Str = 'burner_DA_mutihole';
xlswrite(FileName,{Str},SheetName,'F3');
Str = 't_max';
xlswrite(FileName,{Str},SheetName,'F4');
xlswrite(FileName,t_max,SheetName,'F5');
Str = 'F_max';
xlswrite(FileName,{Str},SheetName,'F6');
xlswrite(FileName,F_max,SheetName,'F7');
Str = 'F_min';
xlswrite(FileName,{Str},SheetName,'F8');
xlswrite(FileName,F_min,SheetName,'F9');
Str = 'i';
xlswrite(FileName,{Str},SheetName,'F10');
xlswrite(FileName,i,SheetName,'F11');
Str = 'F';
xlswrite(FileName,{Str},SheetName,'G3');
xlswrite(FileName,F',SheetName,'G4');

%文件3
burner_DA_mutihole;
close all;

FileName = 'data_burner.xlsx';
SheetName = 'sheet1';

Str = 'burner_DA_mutihole';
xlswrite(FileName,{Str},SheetName,'F3');
Str = 't_max';
xlswrite(FileName,{Str},SheetName,'F4');
xlswrite(FileName,t_max,SheetName,'F5');
Str = 'F_max';
xlswrite(FileName,{Str},SheetName,'F6');
xlswrite(FileName,F_max,SheetName,'F7');
Str = 'F_min';
xlswrite(FileName,{Str},SheetName,'F8');
xlswrite(FileName,F_min,SheetName,'F9');
Str = 'i';
xlswrite(FileName,{Str},SheetName,'F10');
xlswrite(FileName,i,SheetName,'F11');
Str = 'F';
xlswrite(FileName,{Str},SheetName,'G3');
xlswrite(FileName,F',SheetName,'G4');

disp('done');
