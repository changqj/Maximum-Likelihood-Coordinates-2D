function showbase(n,coors,f,xx,yy,cagepoly)
load('colorbar_mlc.mat')
for i = 1:n

    figure
    hold on,axis equal off
%     title('gbc  ');

    h = trisurf(f,xx,yy, coors(i,:));
    set(h,'EdgeAlpha',0)
    shading interp
    caxis([-0.1 1.1]);
    zz = zeros(size(cagepoly(1,:)));
    zz(i) = 1;
    plot3([cagepoly(1,:) cagepoly(1,1)],[cagepoly(2,:) cagepoly(2,1)],...
        [zz zz(1)],'b-o', 'MarkerSize',8,'MarkerFaceColor','b')
    %     light('Position',[1 1 1])
    caxis([-0.1 1.1]);

    colormap(colorbar_mlc)
end
end
