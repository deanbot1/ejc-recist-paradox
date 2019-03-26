function h = scaledarrows(x,y,varargin)


axpos = get(gca,'Position');ax = axpos(1); ay = axpos(2); dax = axpos(3); day = axpos(4);
dxl = diff(get(gca,'Xlim')); mxl = min(get(gca,'Xlim'));
dyl = diff(get(gca,'Ylim')); myl = min(get(gca,'Ylim'));

for j = 1:length(x)-1
    h(j)=annotation('arrow',ax + dax*(x([j j+1])-mxl)/dxl, ay + day*(y([j j+1])-myl)/dyl,varargin{:}); hold on;
end
