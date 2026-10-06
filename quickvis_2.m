function quickvis(filename1,X,Y,Z,WF,npts,cm,iso,A,alpha,beta,gamma,gridmax,me_et,LR,hal,orb,MO)
col=[0,1,0;0.85,0.325,0.098;0.4940,0.1840,0.5560];


%Quick visualization of Cansu's plots of 6/26/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%I-CH3I quickvis(X,Y,Z,WF,npts,cm,0.001,A,pi,pi/2,pi,20,1,2,3,1)
%%%I-CH3Br quickvis(X,Y,Z,WF,npts,cm,0.001,A,0,pi/2,pi/2,20,1,2,2,0)
%%%I-CH3Cl quickvis(X,Y,Z,WF,npts,cm,0.001,A,0,pi/2,pi/2,20,1,2,1,0)
%%%I-C2H5I quickvis(X,Y,Z,WF,npts,cm,0.001,A,0,0,0,20,2,2,3,0)
%%%I-C2H5Br quickvis(X,Y,Z,WF,npts,cm,0.001,A,0,0,0,20,2,2,2,0)
%%%I-C2H5Cl  quickvis(X,Y,Z,WF,npts,cm,0.001,A,0,pi,0,20,2,2,1,0)
%me_et=2;
% alpha = 0;
% beta=pi/2;
% gamma =pi/2;
 Rmat=[(cos(alpha)*cos(beta)*cos(gamma)-sin(alpha)*sin(gamma)),...
         (-cos(gamma)*sin(alpha)-cos(alpha)*cos(beta)*sin(gamma)),cos(alpha)*sin(beta);...
         (cos(alpha)*sin(gamma)+cos(beta)*cos(gamma)*sin(alpha)),...
         (cos(alpha)*cos(gamma)-cos(beta)*sin(alpha)*sin(gamma)),sin(alpha)*sin(beta);...
         -cos(gamma)*sin(beta),sin(beta)*sin(gamma),cos(beta)];

    

xLF=X*Rmat(1,1)+Y*Rmat(1,2)+Z*Rmat(1,3);
yLF=X*Rmat(2,1)+Y*Rmat(2,2)+Z*Rmat(2,3);
zLF=X*Rmat(3,1)+Y*Rmat(3,2)+Z*Rmat(3,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%% Setting up each individual parent and excited orientations
%%%%%%%%%% (this is the time consuming bit. By defining the original
%%%%%%%%%% wavefunction within a spherical volume of the grid, linear
%%%%%%%%%% or cubic interpolation is possible, quicker than splines)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars WFr
WFr=zeros(npts,npts,npts,2);
%gridmax=30;
if orb==1 || orb ==2
  WFr(:,:,:,1)=interp3(X,Y,Z,WF(:,:,:,1),xLF,yLF,zLF,'linear',0);
  WFr(:,:,:,2)=interp3(X,Y,Z,WF(:,:,:,2),xLF,yLF,zLF,'linear',0);
end
vecz(:,1)=transpose(Rmat)*[0;0;gridmax+1];
 vecy(:,1)=transpose(Rmat)*[0;gridmax+1;0];
 vecx(:,1)=transpose(Rmat)*[gridmax+1;0;0];
fig=figplot(X,Y,Z,WFr(:,:,:,LR),vecx(:,1),vecy(:,1),vecz(:,1),gridmax,iso,1,1,1,cm,0.5);
hold on
Af(:,:)=A(:,:)*Rmat(:,:);
ax=gca;
light(ax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%% Plotting the Methyls
%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if me_et==1
%%%H atoms
plot3(Af(1:3,1),Af(1:3,2),Af(1:3,3),"o",'Color','black','MarkerSize',6,...
    'MarkerFaceColor','white');

%%% C atoms
plot3(Af(4,1),Af(4,2),Af(4,3),"o",'Color','black','MarkerSize',6,...
    'MarkerFaceColor','black');

%%% atom X
plot3(Af(5,1),Af(5,2),Af(5,3),"o",'Color','black','MarkerSize',10,...
    'MarkerFaceColor',col(hal,:));

%%% iodine
plot3(Af(6,1),Af(6,2),Af(6,3),"o",'Color','black','MarkerSize',10,...
    'MarkerFaceColor',col(3,:));

%%% covalent bonds
lw=2;
plot3([Af(1,1) Af(4,1)],[Af(1,2) Af(4,2)],[Af(1,3) Af(4,3)],'-','Color','black','LineWidth',lw);
plot3([Af(2,1) Af(4,1)],[Af(2,2) Af(4,2)],[Af(2,3) Af(4,3)],'-','Color','black','LineWidth',lw);
plot3([Af(3,1) Af(4,1)],[Af(3,2) Af(4,2)],[Af(3,3) Af(4,3)],'-','Color','black','LineWidth',lw);
plot3([Af(4,1) Af(5,1)],[Af(4,2) Af(5,2)],[Af(4,3) Af(5,3)],'-','Color','black','LineWidth',lw);

%%% draw axes
if orb == 0 || orb == 2
p1=quiver3(0,0,0,(vecz(1,1)/2),(vecz(2,1)/2),(vecz(3,1))/2);
p1m=quiver3(0,0,0,-(vecz(1,1)/2),-(vecz(2,1)/2),-(vecz(3,1)/2));
p2m=quiver3(0,0,0,-(vecx(1,1)/2),-(vecx(2,1)/2),-(vecx(3,1)/2));
if hal==3
 p2m=quiver3(0,0,0,vecx(1,1)/2,vecx(2,1)/2,vecx(3,1)/2);
end
 p2m.Color=[0 0 0];
 p2m.LineWidth=2;
 p2m.LineStyle="--"';
 p2m.ShowArrowHead='off';
p3m=quiver3(0,0,0,-vecy(1,1)/2,-vecy(2,1)/2,-vecy(3,1)/2);
p1.Color=[0 0 0];
p1.LineStyle="-"';
p1.LineWidth=2;

p1.ShowArrowHead='on';
p1.MaxHeadSize=0.5;
p2=quiver3(0,0,0,(vecx(1,1)/2),(vecx(2,1)/2),(vecx(3,1)/2));
if hal==3
p2=quiver3(0,0,0,-(vecx(1,1)/2),-(vecx(2,1)/2),-(vecx(3,1))/2);
end
p2.ShowArrowHead='on';
p2.MaxHeadSize=0.5;
p2.LineStyle="-"';
p2.LineWidth=2;

p3=quiver3(0,0,0,vecy(1,1)/2,vecy(2,1)/2,vecy(3,1)/2);
p3.ShowArrowHead='on';
p3.LineStyle="-"';
p3.LineWidth=2;

p3.MaxHeadSize=0.5;
p2.Color=[0 0 0];
p2.LineWidth=2;
p2.ShowArrowHead='on';
% p2m.Color=[0 0 0];
% p2m.LineWidth=2;
% p2m.LineStyle="--"';
% p2m.ShowArrowHead='off';
 p3m.Color=[0 0 0];
 p3m.LineWidth=2;
 p3m.LineStyle="--"';
 p3m.ShowArrowHead='off';
p1m.Color=[0 0 0];
p1m.LineWidth=2;
p1m.LineStyle="--"';
p1m.ShowArrowHead='off';
p3.Color=[0 0 0];
p3.LineWidth=2;
p3.ShowArrowHead='on';
% if ne(hal,3) 
 text(vecz(1,1)/2,vecz(2,1)/2,vecz(3,1)/2+1,'y','FontSize',16);
 text(vecx(1,1)/2+1,vecx(2,1)/2,vecx(3,1)/2+1,'x','FontSize',16);
 text(vecy(1,1)/2,vecy(2,1)/2+1,vecy(3,1)/2,'z','FontSize',16);
% end
% if hal == 3
% text(vecz(1,1)/2+2,vecz(2,1)/2,vecz(3,1)/2,'x','FontSize',16);
% text((vecx(1,1))/2,vecx(2,1)/2,-(vecx(3,1)/2-2),'z','FontSize',16);
% text(vecy(1,1)/2,vecy(2,1)/2+1,vecy(3,1)/2,'y','FontSize',16);
% end
end
end
lw=2;
if me_et==2

plot3(Af(1:5,1),Af(1:5,2),Af(1:5,3),"o",'Color','black','MarkerSize',6,...
     'MarkerFaceColor','white');
plot3(Af(6:7,1),Af(6:7,2),Af(6:7,3),"o",'Color','black','MarkerSize',8,...
    'MarkerFaceColor','black');
plot3(Af(8,1),Af(8,2),Af(8,3),"o",'Color','black','MarkerSize',10,...
    'MarkerFaceColor',col(hal,:));
plot3(Af(9,1),Af(9,2),Af(9,3),"o",'Color','black','MarkerSize',10,...
    'MarkerFaceColor',col(3,:));

plot3([Af(1,1) Af(6,1)],[Af(1,2) Af(6,2)],[Af(1,3) Af(6,3)],'-','Color','black','LineWidth',lw);
plot3([Af(2,1) Af(6,1)],[Af(2,2) Af(6,2)],[Af(2,3) Af(6,3)],'-','Color','black','LineWidth',lw);
plot3([Af(3,1) Af(7,1)],[Af(3,2) Af(7,2)],[Af(3,3) Af(7,3)],'-','Color','black','LineWidth',lw);
plot3([Af(4,1) Af(7,1)],[Af(4,2) Af(7,2)],[Af(4,3) Af(7,3)],'-','Color','black','LineWidth',lw);
plot3([Af(5,1) Af(7,1)],[Af(5,2) Af(7,2)],[Af(5,3) Af(7,3)],'-','Color','black','LineWidth',lw);
plot3([Af(8,1) Af(6,1)],[Af(8,2) Af(6,2)],[Af(8,3) Af(6,3)],'-','Color','black','LineWidth',lw);
plot3([Af(6,1) Af(7,1)],[Af(6,2) Af(7,2)],[Af(6,3) Af(7,3)],'-','Color','black','LineWidth',lw);
%end
 
%  if hal==1
% plot3(Af(3,1),Af(3,2),Af(3,3),"o",'Color','black','MarkerSize',8,...
%     'MarkerFaceColor','black');
% plot3(Af(6,1),Af(6,2),Af(6,3),"o",'Color','black','MarkerSize',8,...
%     'MarkerFaceColor','black');
% plot3(Af(9,1),Af(9,2),Af(9,3),"o",'Color','black','MarkerSize',10,...
%     'MarkerFaceColor',col(hal,:));
% plot3(Af(1,1),Af(1,2),Af(1,3),"o",'Color','black','MarkerSize',10,...
%     'MarkerFaceColor',col(3,:));
% plot3(Af(2,1),Af(2,2),Af(2,3),"o",'Color','black','MarkerSize',6,...
%      'MarkerFaceColor','white');
% plot3(Af(4:5,1),Af(4:5,2),Af(4:5,3),"o",'Color','black','MarkerSize',6,...
%      'MarkerFaceColor','white');
% plot3(Af(7:8,1),Af(7:8,2),Af(7:8,3),"o",'Color','black','MarkerSize',6,...
%      'MarkerFaceColor','white'); 
% plot3([Af(3,1) Af(2,1)],[Af(3,2) Af(2,2)],[Af(3,3) Af(2,3)],'-','Color','black','LineWidth',lw);
% plot3([Af(3,1) Af(4,1)],[Af(3,2) Af(4,2)],[Af(3,3) Af(4,3)],'-','Color','black','LineWidth',lw);
% plot3([Af(3,1) Af(5,1)],[Af(3,2) Af(5,2)],[Af(3,3) Af(5,3)],'-','Color','black','LineWidth',lw);
% plot3([Af(3,1) Af(6,1)],[Af(3,2) Af(6,2)],[Af(3,3) Af(6,3)],'-','Color','black','LineWidth',lw);
% plot3([Af(6,1) Af(7,1)],[Af(6,2) Af(7,2)],[Af(6,3) Af(7,3)],'-','Color','black','LineWidth',lw);
% plot3([Af(6,1) Af(8,1)],[Af(6,2) Af(8,2)],[Af(6,3) Af(8,3)],'-','Color','black','LineWidth',lw);
% plot3([Af(6,1) Af(9,1)],[Af(6,2) Af(9,2)],[Af(6,3) Af(9,3)],'-','Color','black','LineWidth',lw);
% p1=quiver3(0,0,0,-(vecz(1,1)/2),-(vecz(2,1)/2),-(vecz(3,1))/2);
% p2=quiver3(0,0,0,-(vecx(1,1)/2),-(vecx(2,1)/2),-(vecx(3,1)/2));
% p3=quiver3(0,0,0,(vecy(1,1)/2),(vecy(2,1)/2),(vecy(3,1)/2));
% p1m=quiver3(0,0,0,(vecz(1,1)/2),(vecz(2,1)/2),(vecz(3,1)/2));
% p2m=quiver3(0,0,0,(vecx(1,1)/2),(vecx(2,1)/2),(vecx(3,1)/2));
% p3m=quiver3(0,0,0,-(vecy(1,1)/2),-(vecy(2,1)/2),-(vecy(3,1)/2));
% p1.Color=[0.5 0.5 0.5];
% p1.LineWidth=2;
% p1.ShowArrowHead='on';
% p2.Color=[0.5 0.5 0.5];
% p2.LineWidth=2;
% p2.ShowArrowHead='on';
% p2m.Color=[0.5 0.5 0.5];
% p2m.LineWidth=2;
% p2m.LineStyle="--"';
% p2m.ShowArrowHead='off';
% p3m.Color=[0.5 0.5 0.5];
% p3m.LineWidth=2;
% p3m.LineStyle="--"';
% p3m.ShowArrowHead='off';
% p1m.Color=[0.5 0.5 0.5];
% p1m.LineWidth=2;
% p1m.LineStyle="--"';
% p1m.ShowArrowHead='off';
% p3.Color=[0.5 0.5 0.5];
% p3.LineWidth=2;
% p3.ShowArrowHead='on';
% text(-vecz(1,1),-vecz(2,1),-(vecz(3,1)-1),'z','FontSize',16);
% text(-(vecx(1,1)+1),-vecx(2,1),-(vecx(3,1)+1),'x','FontSize',16);
% text(vecy(1,1),(vecy(2,1)+1),vecy(3,1),'y','FontSize',16);
% end

% if ne(hal,1)
%     %Ethyls
% plot3(Af(1,1),Af(1,2),Af(1,3),"o",'Color','black','MarkerSize',8,...
%     'MarkerFaceColor','black');
% plot3(Af(2,1),Af(2,2),Af(2,3),"o",'Color','black','MarkerSize',8,...
%     'MarkerFaceColor','black');
% plot3(Af(3,1),Af(3,2),Af(3,3),"o",'Color','black','MarkerSize',10,...
%     'MarkerFaceColor',col(hal,:));
% plot3(Af(4,1),Af(4,2),Af(4,3),"o",'Color','black','MarkerSize',10,...
%     'MarkerFaceColor',col(3,:));
% plot3(Af(5:9,1),Af(5:9,2),Af(5:9,3),"o",'Color','black','MarkerSize',6,...
%      'MarkerFaceColor','white');
% plot3([Af(1,1) Af(2,1)],[Af(1,2) Af(2,2)],[Af(1,3) Af(2,3)],'-','Color','black','LineWidth',lw);
% plot3([Af(1,1) Af(5,1)],[Af(1,2) Af(5,2)],[Af(1,3) Af(5,3)],'-','Color','black','LineWidth',lw);
% plot3([Af(1,1) Af(6,1)],[Af(1,2) Af(6,2)],[Af(1,3) Af(6,3)],'-','Color','black','LineWidth',lw);
% plot3([Af(1,1) Af(3,1)],[Af(1,2) Af(3,2)],[Af(1,3) Af(3,3)],'-','Color','black','LineWidth',lw);
% plot3([Af(2,1) Af(7,1)],[Af(2,2) Af(7,2)],[Af(2,3) Af(7,3)],'-','Color','black','LineWidth',lw);
% plot3([Af(2,1) Af(8,1)],[Af(2,2) Af(8,2)],[Af(2,3) Af(8,3)],'-','Color','black','LineWidth',lw);
% plot3([Af(2,1) Af(9,1)],[Af(2,2) Af(9,2)],[Af(2,3) Af(9,3)],'-','Color','black','LineWidth',lw);
if orb == 0 || orb == 2
% p1=quiver3(0,0,0,(vecz(1,1)/2),(vecz(2,1)/2),(vecz(3,1)/2));
% p1m=quiver3(0,0,0,-(vecz(1,1)/2),-(vecz(2,1)/2),-(vecz(3,1)/2));
% p2m=quiver3(0,0,0,-(vecx(1,1)/2),-(vecx(2,1)/2),-(vecx(3,1)/2));
% p3m=quiver3(0,0,0,-(vecy(1,1)/2),-(vecy(2,1)/2),-(vecy(3,1)/2));
% p1.Color=[0.5 0.5 0.5];
% p1.LineWidth=2;
% p1.ShowArrowHead='on';
% p1.MaxHeadSize=0.5;
% p2=quiver3(0,0,0,(vecx(1,1)/2),(vecx(2,1)/2),(vecx(3,1)/2));
% p3=quiver3(0,0,0,(vecy(1,1)/2),(vecy(2,1)/2),(vecy(3,1)/2));
% p2.Color=[0.5 0.5 0.5];
% p2.LineWidth=2;
% p2.ShowArrowHead='on';
% p2.MaxHeadSize=0.5;
% p2m.Color=[0.5 0.5 0.5];
% p2m.LineWidth=2;
% p2m.LineStyle="--"';
% p2m.ShowArrowHead='off';
% p3m.Color=[0.5 0.5 0.5];
% p3m.LineWidth=2;
% p3m.LineStyle="--"';
% p3m.ShowArrowHead='off';
% p1m.Color=[0.5 0.5 0.5];
% p1m.LineWidth=2;
% p1m.LineStyle="--"';
% p1m.ShowArrowHead='off';
% p3.Color=[0.5 0.5 0.5];
% p3.LineWidth=2;
% p3.ShowArrowHead='on';
% p3.MaxHeadSize=0.5;
% text((vecz(1,1)/2),(vecz(2,1)/2),(vecz(3,1)/2)+2,'z','FontSize',16);
% text((vecx(1,1)/2)+2,(vecx(2,1)/2),(vecx(3,1)/2),'x','FontSize',16);
% text((vecy(1,1)/2),(vecy(2,1)/2)+1,(vecy(3,1)/2),'y','FontSize',16);
end
%end
% if hal == 1
% text(vecz(1,1),vecz(2,1),vecz(3,1)+1,'x','FontSize',12);
% text((vecx(1,1)),vecx(2,1),-(vecx(3,1)+1),'z','FontSize',12);
% text(vecy(1,1),vecy(2,1)+1,vecy(3,1),'y','FontSize',12);
% end
end


if orb==0
% saveas(fig,strcat(filename1(1,1:(size(filename1,2)-4)),num2str(MO),'_axes.tif'));
end
if orb==1 
saveas(fig,strcat(filename1(1,1:(size(filename1,2)-4)),num2str(MO),'_HF.tif'));
end
if orb==2 
saveas(fig,strcat(filename1(1,1:(size(filename1,2)-4)),num2str(MO),'_HF_axes.tif'));
end
end
