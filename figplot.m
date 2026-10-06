function f=figplot(X,Y,Z,F,vecx,vecy,vecz,axval,iso,no,LR,type,cm,al)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%      Routine to Plot Dyson Orbitals at specified iso value      %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%              Last Modified:    Jan 29 2022                      %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=gcf;
clf
f.Color='white';
colormap(cm);
if type==1
    hold on
subplot(1,no,LR);
isoval=iso;
 axis([-axval axval -axval axval -axval axval], 'off');
 h=gca;
 h.XAxis.TickLength=[0 0];
 h.YAxis.TickLength=[0 0];
 h.ZAxis.TickLength=[0 0];
 h.XAxis.TickLabel=[];
 h.YAxis.TickLabel=[];
 h.ZAxis.TickLabel=[];
 p=patch(isosurface(X,Y,Z,F,isoval));
  p.FaceColor = 'blue'; 
  p.EdgeColor = 'none';
 p.FaceAlpha = al;
 p.EdgeAlpha = 0.0;
 p.SpecularStrength=1.0;
%     purple = [150 37 255]/255;
%     peach  = [255 150 120]/255;
%     set(p, 'FaceColor', peach, 'EdgeColor', 'none','FaceAlpha', '0.5');
   
 hold on
 q=patch(isosurface(X,Y,Z,F,-isoval));
 q.FaceColor = 'red';
 q.EdgeColor = 'none';
 q.FaceAlpha = al;
 q.EdgeAlpha = 0.0;
 q.SpecularStrength=1.0;
%   set(q, 'FaceColor', purple, 'EdgeColor', 'none','FaceAlpha', '0.5');
axis tight;
 grid on
 daspect([1 1 1])
 view([3 -90])
%  x1=[0,0];
%  y1=[0,0];
%  z1=[-axval,axval];
% plot3(x1,y1,z1,'.-black','LineWidth',3);
%  x1=[-axval,axval];
%  y1=[0,0];
%  z1=[0,0];
%  plot3(x1,y1,z1,'.-black','LineWidth',3);
%  x1=[0,0];
%  y1=[-axval,axval];
%  z1=[0,0];
%  plot3(x1,y1,z1,'.-black','LineWidth',3);
rotate3d;
hold off
end
 hold on

% ax=gca;
% [spx,spy,spz]=sphere;
% sps=surf(spx*axval,spy*axval,spz*axval);
% sps.FaceAlpha=0.0;
% sps.FaceColor=[0.7,0.7,0.7];
% sps.EdgeColor=[0.5,0.5,0.5];
% %sps.FaceColor=[1,1,1];
% %sps.EdgeColor=[1,1,1];
% sps.EdgeAlpha=0.0;
% %p=plot3(vec(1,1),vec(2,1),vec(3,1),'o');
% %p1=arrow([-vecz(1,1),vecz(1,1)],[-vecz(2,1),vecz(2,1)],[-vecz(3,1),vecz(3,1)],'LineWidth',3,'Color','black');
% % p1=quiver3(0,0,0,vecz(1,1),vecz(2,1),vecz(3,1));
% % %arrow3([-vecz(1,1),vecz(1,1)],[-vecz(2,1),vecz(2,1)],[-vecz(3,1),vecz(3,1)],'LineWidth',3,'Color','black');
% % p1.Color='red';
% % p1.LineWidth=2;
% % p1.ShowArrowHead='on';
% % %annotation('textbox',[vecz(1,1) vecz(2,1) vecz(3,1)],'String','z');
% % 
% % 
% % p2=quiver3(0,0,0,vecx(1,1),vecx(2,1),vecx(3,1));
% % p3=quiver3(0,0,0,vecy(1,1),vecy(2,1),vecy(3,1));
% % p2.Color='green';
% % p2.LineWidth=2;
% % p2.ShowArrowHead='on';
% % p3.Color='blue';
% % p3.LineWidth=2;
% % p3.ShowArrowHead='on';
% %colormap(cm);
%  %view(45,45);
% view([90,90,90]);
% %campos=[0,0,10];
% % ax.XLim=[-axval/7 axval/7];
% % ax.YLim=[-axval/7 axval/7];
% % ax.ZLim=[-axval/7 axval/7];
  drawnow

end