%190414
%Dashaway
%���ݶ�ȡ����ͼ


clear;
close all;
%��ȡexcel
FileName = 'data_burner.xlsx';

t_max = 0;
F_max = 0;

%����1
SheetName = 'sheet1';
i = xlsread(FileName,SheetName,'C11') - 1;
t = xlsread(FileName,SheetName,'C5');
if(t_max < t)
    t_max = t;
end
f = xlsread(FileName,SheetName,'C7');
if(F_max < f)
    F_max = f;
end
t_1 = (5e-5)*(0:1:i - 1);
la = strcat('D3:D',num2str(i + 3));
F_1 = xlsread(FileName,SheetName,la);
F_1 = F_1';

%����2
SheetName = 'sheet2';
i = xlsread(FileName,SheetName,'C11') - 1;
t = xlsread(FileName,SheetName,'C5');
if(t_max < t)
    t_max = t;
end
f = xlsread(FileName,SheetName,'C7');
if(F_max < f)
    F_max = f;
end
t_2 = (5e-5)*(0:1:i - 1);
la = strcat('D3:D',num2str(i + 3));
F_2 = xlsread(FileName,SheetName,la);
F_2 = F_2';

%����3
SheetName = 'sheet3';
i = xlsread(FileName,SheetName,'C11') - 1;
t = xlsread(FileName,SheetName,'C5');
if(t_max < t)
    t_max = t;
end
f = xlsread(FileName,SheetName,'C7');
if(F_max < f)
    F_max = f;
end
t_3 = (5e-5)*(0:1:i - 1);
la = strcat('D3:D',num2str(i + 3));
F_3 = xlsread(FileName,SheetName,la);
F_3 = F_3';

%����4
SheetName = 'sheet4';
i = xlsread(FileName,SheetName,'C11') - 1;
t = xlsread(FileName,SheetName,'C5');
if(t_max < t)
    t_max = t;
end
f = xlsread(FileName,SheetName,'C7');
if(F_max < f)
    F_max = f;
end
t_4 = (5e-5)*(0:1:i - 1);
la = strcat('D3:D',num2str(i + 3));
F_4 = xlsread(FileName,SheetName,la);
F_4 = F_4';

%����5
SheetName = 'sheet5';
i = xlsread(FileName,SheetName,'C11') - 1;
t = xlsread(FileName,SheetName,'C5');
if(t_max < t)
    t_max = t;
end
f = xlsread(FileName,SheetName,'C7');
if(F_max < f)
    F_max = f;
end
t_5 = (5e-5)*(0:1:i - 1);
la = strcat('D3:D',num2str(i + 3));
F_5 = xlsread(FileName,SheetName,la);
F_5 = F_5';

%����6
SheetName = 'sheet6';
i = xlsread(FileName,SheetName,'C11') - 1;
t = xlsread(FileName,SheetName,'C5');
if(t_max < t)
    t_max = t;
end
f = xlsread(FileName,SheetName,'C7');
if(F_max < f)
    F_max = f;
end
t_6 = (5e-5)*(0:1:i - 1);
la = strcat('D3:D',num2str(i + 3));
F_6 = xlsread(FileName,SheetName,la);
F_6 = F_6';




%����ͼ��
close all;
figure;
hold on;
box on;

plot(t_1,F_1,'color','k','LineStyle','-','LineWidth',0.8);
plot(t_2,F_2,'color','k','LineStyle',':','LineWidth',1.5);

plot(t_3,F_3,'color','k','LineStyle','-','LineWidth',1.5);
plot(t_4,F_4,'color','k','LineStyle',':','LineWidth',1);

plot(t_5,F_5,'color','k','LineStyle','-.','LineWidth',1);
plot(t_6,F_6,'color','k','LineStyle','--','LineWidth',1);

legend('����1','����2','����3','����4','����5','����6');

axis ([-0.05*t_max,1.05*t_max,-0.05*F_max,1.05*F_max ]);

title('����');
xlabel('ʱ��(s)');
ylabel('��(N)');


disp('done');
