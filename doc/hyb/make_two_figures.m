close all
bdwidth = 5;
topbdwidth = 30;
set(0,'Units','pixels') 
scnsize = get(0,'ScreenSize');
pos1  = [bdwidth,... 
	1/8*scnsize(4) + bdwidth,...
	scnsize(3)/2 - 2*bdwidth,...
	scnsize(4)*3/4 - (topbdwidth + bdwidth)];
pos2 = [pos1(1) + scnsize(3)/2,...
	pos1(2),...
	pos1(3),...
	pos1(4)];
figure('Position',pos1);
figure('Position',pos2);
