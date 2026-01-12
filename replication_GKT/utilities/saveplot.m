function saveplot(NAME,DIR)
ratio=16/9;
trim_hor=0.6;
trim_top=0.05;

name=NAME;
width_plot=6;
height_plot=width_plot/ratio;

fac=1+trim_hor/width_plot;

set(gcf, 'PaperPosition', [0-trim_hor*fac 0+trim_top  width_plot+trim_hor*2 height_plot+trim_top]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [ width_plot height_plot]); %Set the paper to have width 5 and height 5.
%saveas(gcf, name, 'pdf',) %Save figure
if strcmp(DIR,'')
namefull=[name];
else
namefull=[DIR '\' name];
end

print(gcf, namefull, '-dpdf'  , '-r600', '-cmyk');
print(gcf, namefull, '-depsc' , '-r600', '-cmyk');
end