%190414
%Dashaway
%���ݼ���ͱ���


clear;
close all;
%д��excel

%����1
%װҩ����
D = 0.300;         %ҩ����ֱ��(m)
r = 0.010;        %΢Բ���뾶(m)
l = 0.130;        %���ľ�(m)
R = 0.01;       %���Ŀ�(m)
m_s = 8;        %����ҩ��
n_s = 3;        %������
Lp = 0.500;      %ҩ������(m)
burner_DA_calculate;
close all;

FileName = 'data_burner.xlsx';
SheetName = 'sheet1';

Str = '����1';
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

%����2
%װҩ����
D = 0.300;         %ҩ����ֱ��(m)
r = 0.010;        %΢Բ���뾶(m)
l = 0.120;        %���ľ�(m)
R = 0.01;       %���Ŀ�(m)
m_s = 8;        %����ҩ��
n_s = 3;        %������
Lp = 0.500;      %ҩ������(m)
burner_DA_calculate;
close all;

FileName = 'data_burner.xlsx';
SheetName = 'sheet2';

Str = '����2';
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

%����3
%װҩ����
D = 0.300;         %ҩ����ֱ��(m)
r = 0.010;        %΢Բ���뾶(m)
l = 0.130;        %���ľ�(m)
R = 0.05;       %���Ŀ�(m)
m_s = 8;        %����ҩ��
n_s = 3;        %������
Lp = 0.500;      %ҩ������(m)
burner_DA_calculate;
close all;

FileName = 'data_burner.xlsx';
SheetName = 'sheet3';

Str = '����3';
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

%����4
%װҩ����
D = 0.300;         %ҩ����ֱ��(m)
r = 0.005;        %΢Բ���뾶(m)
l = 0.130;        %���ľ�(m)
R = 0.01;       %���Ŀ�(m)
m_s = 8;        %����ҩ��
n_s = 3;        %������
Lp = 0.500;      %ҩ������(m)
burner_DA_calculate;
close all;

FileName = 'data_burner.xlsx';
SheetName = 'sheet4';

Str = '����4';
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

%����5
%װҩ����
D = 0.300;         %ҩ����ֱ��(m)
r = 0.010;        %΢Բ���뾶(m)
l = 0.130;        %���ľ�(m)
R = 0.01;       %���Ŀ�(m)
m_s = 6;        %����ҩ��
n_s = 3;        %������
Lp = 0.500;      %ҩ������(m)
burner_DA_calculate;
close all;

FileName = 'data_burner.xlsx';
SheetName = 'sheet5';

Str = '����5';
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

%����6
%װҩ����
D = 0.300;         %ҩ����ֱ��(m)
r = 0.010;        %΢Բ���뾶(m)
l = 0.130;        %���ľ�(m)
R = 0.01;       %���Ŀ�(m)
m_s = 6;        %����ҩ��
n_s = 5;        %������
Lp = 0.500;      %ҩ������(m)
burner_DA_calculate;
close all;

FileName = 'data_burner.xlsx';
SheetName = 'sheet6';

Str = '����6';
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
