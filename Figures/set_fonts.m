function [ output_args ] = set_fonts( fig, gca, fontsize )
%SET_FONTS Set font in plot to given size
% Sebastian Wallot, Dan Mønster, 2016
    set(gca,'FontSize',fontsize,'fontWeight','normal')
    set(findall(fig,'type','text'),'FontSize',fontsize,'fontWeight','normal')
end

