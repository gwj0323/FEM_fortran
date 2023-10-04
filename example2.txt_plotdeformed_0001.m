 % Plotting Un-Deformed and Deformed Structure
 close all
 clear all
 example2.txt_plotdeformed_data0001;
 % Make plot
 figure
 set(gcf,'Name','Deformed')
 subplot(2,1,2)
 hold on
 for e = 1:size(IX,1)
    edof = [2*IX(e,1)-1 2*IX(e,1) 2*IX(e,2)-1 2*IX(e,2)];
    xx = X(IX(e,1:2),1) + D(edof(1:2:4));
    yy = X(IX(e,1:2),2) + D(edof(2:2:4));
    plot(xx,yy,'b','LineWidth',1.5)
 end
 title('Deformed')
 axis equal
 xaxes = get(gca,'xlim');
 yaxes = get(gca,'ylim');
 axis off;
 subplot(2,1,1)
 hold on
 for e = 1:size(IX,1)
    xx = X(IX(e,1:2),1);
    yy = X(IX(e,1:2),2);
    plot(xx,yy,'b','LineWidth',1.5)
 end
 title('Undeformed')
 axis([min(xaxes(1),min(X(:,1))) max(xaxes(2),max(X(:,1)))...
  min(yaxes(1),min(X(:,2))) max(yaxes(2),max(X(:,2))) ]);
 axis off;
 hold off
 set(gcf,'color',[ 1  1 1]);
