% function quickvis(filename1,X,Y,Z,WF,npts,cm,iso,A,alpha,beta,gamma,gridmax,me_et,LR,hal,orb,MO)
function quickvis(filename1,X,Y,Z,WF,npts,cm,iso,A,alpha,beta,gamma,gridmax,me_et,LR,hal,orb,MO)
col=[0,0.69,0.31;0.75,0.000,0.000;0.4340,0.1880,0.6270];

alpha = 0;
beta=pi;
gamma =pi/2;
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
fig=figplot(X,Y,Z,WFr(:,:,:,LR),vecx(:,1),vecy(:,1),vecz(:,1),gridmax,iso,1,1,1,cm,0.3);
hold on
Af(:,:)=A(:,:)*Rmat(:,:);
ax=gca;
camlight('headlight');
lighting gouraud;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%% Plotting the Methyls
%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if me_et==1
%%%H atoms
plot3(Af(1:3,1),Af(1:3,2),Af(1:3,3),"o",'Color','black','MarkerSize',10,...
    'MarkerFaceColor','white');
% text(Af(1,1),Af(1,2),Af(1,3), 'Ha',  'FontSize', 16, 'FontWeight', 'bold', 'Color', 'black');
%% C atoms
plot3(Af(4,1),Af(4,2),Af(4,3),"o",'Color','black','MarkerSize',13,...
    'MarkerFaceColor','black');
%% atom X
plot3(Af(5,1),Af(5,2),Af(5,3),"o",'Color','black','MarkerSize',17,...
    'MarkerFaceColor',col(hal,:));
%% iodine
plot3(Af(6,1),Af(6,2),Af(6,3),"o",'Color','black','MarkerSize',17,...
    'MarkerFaceColor',col(3,:));
% % plot3(Af(1:3,1),Af(1:3,2),Af(1:3,3),"o",'Color','black','MarkerSize',6,...
% %     'MarkerFaceColor','white');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%CH3X plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot3(Af(1:3,1),Af(1:3,2),Af(1:3,3),"o",'Color','black','MarkerSize',13,...
%     'MarkerFaceColor','white');
% % % text(Af(1,1),Af(1,2),Af(1,3), 'Ha',  'FontSize', 16, 'FontWeight', 'bold', 'Color', 'black');
% % %% C atoms
% plot3(Af(4,1),Af(4,2),Af(4,3),"o",'Color','black','MarkerSize',17,...
%     'MarkerFaceColor','black');
% % %% atom X
% plot3(Af(5,1),Af(5,2),Af(5,3),"o",'Color','black','MarkerSize',17,...
%     'MarkerFaceColor',col(hal,:));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% C atoms
% plot3(Af(4,1),Af(4,2),Af(4,3),"o",'Color','black','MarkerSize',6,...
%     'MarkerFaceColor','black');
% %%% atom X
% plot3(Af(5,1),Af(5,2),Af(5,3),"o",'Color','black','MarkerSize',10,...
%     'MarkerFaceColor',col(hal,:));
% 
% %%% iodine
% plot3(Af(6,1),Af(6,2),Af(6,3),"o",'Color','black','MarkerSize',10,...
%     'MarkerFaceColor',col(3,:));
%%% covalent bonds
% lw=2;
% plot3([Af(1,1) Af(4,1)],[Af(1,2) Af(4,2)],[Af(1,3) Af(4,3)],'-','Color','black','LineWidth',lw);
% plot3([Af(2,1) Af(4,1)],[Af(2,2) Af(4,2)],[Af(2,3) Af(4,3)],'-','Color','black','LineWidth',lw);
% plot3([Af(3,1) Af(4,1)],[Af(3,2) Af(4,2)],[Af(3,3) Af(4,3)],'-','Color','black','LineWidth',lw);
% plot3([Af(4,1) Af(5,1)],[Af(4,2) Af(5,2)],[Af(4,3) Af(5,3)],'-','Color','black','LineWidth',lw);
% text(Af(2,1)+2, Af(2,2)-0.5, Af(2,3), 'H_a','FontSize', 16,'Color', 'black');
% text(Af(3,1)+3, Af(3,2)+0.9, Af(3,3), 'H_b','FontSize', 16, 'Color', 'black');
% text(Af(5,1), Af(5,2)+0.5, Af(5,3), 'X','FontSize', 16, 'Color', 'black');
%%% draw axes
if orb == 0 || orb == 2
% p1=quiver3(0,0,0,(vecz(1,1)/2),(vecz(2,1)/2),(vecz(3,1))/2);
% p1m=quiver3(0,0,0,-(vecz(1,1)/2),-(vecz(2,1)/2),-(vecz(3,1)/2));
% p2m=quiver3(0,0,0,-(vecx(1,1)/2),-(vecx(2,1)/2),-(vecx(3,1)/2));
% % if hal==3
% %  p2m=quiver3(0,0,0,vecx(1,1)/2,vecx(2,1)/2,vecx(3,1)/2);
% % end
% p3m=quiver3(0,0,0,-vecy(1,1)/2,-vecy(2,1)/2,-vecy(3,1)/2);
% p2=quiver3(0,0,0,(vecx(1,1)/2),(vecx(2,1)/2),(vecx(3,1))/2);
% % if hal==3
% % p2=quiver3(0,0,0,-(vecx(1,1)/2),-(vecx(2,1)/2),-(vecx(3,1))/2);
% % end
% p3=quiver3(0,0,0,vecy(1,1)/2,vecy(2,1)/2,vecy(3,1)/2);
% set(p1m, 'Color', [0 0 0], 'LineWidth', 2, 'LineStyle', '--', 'ShowArrowHead', 'off');
% set(p2m, 'Color', [0 0 0], 'LineWidth', 2, 'LineStyle', '--', 'ShowArrowHead', 'off');
% set(p3m, 'Color', [0 0 0], 'LineWidth', 2, 'LineStyle', '--', 'ShowArrowHead', 'off');
% set(p1, 'Color', [0 0 0], 'LineWidth', 2, 'LineStyle', '-', 'ShowArrowHead', 'off');
% set(p2, 'Color', [0 0 0], 'LineWidth', 2, 'LineStyle', '-', 'ShowArrowHead', 'off');
% set(p3, 'Color', [0 0 0], 'LineWidth', 2, 'LineStyle', '-', 'ShowArrowHead', 'off','MaxHeadSize',0.5);
% axis equal;
% if ne(hal,3) 
%  text(vecz(1,1)/2,vecz(2,1)/2,vecz(3,1)/2,'y','FontSize',16);
%  text(vecx(1,1)/2+1,vecx(2,1)/2+1,vecx(3,1)/2+1,'z','FontSize',16);
%  text(vecy(1,1)/2-1,vecy(2,1)/2-1,vecy(3,1)/2-1,'x','FontSize',16);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% ethyls %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lw=2;
if me_et==2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%C2H5X plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% H atoms (1 to 5)
%     plot3(Af(1:5,1),Af(1:5,2),Af(1:5,3),"o",'Color','black','MarkerSize',10,...
%          'MarkerFaceColor','white');
%     %%% C atoms (6 and 7)
%     plot3(Af(6:7,1),Af(6:7,2),Af(6:7,3),"o",'Color','black','MarkerSize',13,...
%         'MarkerFaceColor','black');
%     %%% Halogen atom X (8)
%     plot3(Af(8,1),Af(8,2),Af(8,3),"o",'Color','black','MarkerSize',17,...
%         'MarkerFaceColor',col(hal,:));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% plot3(Af(1:5,1),Af(1:5,2),Af(1:5,3),"o",'Color','black','MarkerSize',6,...
%      'MarkerFaceColor','white');
% plot3(Af(6:7,1),Af(6:7,2),Af(6:7,3),"o",'Color','black','MarkerSize',8,...
%     'MarkerFaceColor','black');
% plot3(Af(8,1),Af(8,2),Af(8,3),"o",'Color','black','MarkerSize',10,...
%     'MarkerFaceColor',col(hal,:));
% plot3(Af(9,1),Af(9,2),Af(9,3),"o",'Color','black','MarkerSize',10,...
%     'MarkerFaceColor',col(3,:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot3(Af(1:5,1),Af(1:5,2),Af(1:5,3),"o",'Color','black','MarkerSize',10,...
     'MarkerFaceColor','white');
plot3(Af(6:7,1),Af(6:7,2),Af(6:7,3),"o",'Color','black','MarkerSize',13,...
    'MarkerFaceColor','black');
plot3(Af(8,1),Af(8,2),Af(8,3),"o",'Color','black','MarkerSize',17,...
    'MarkerFaceColor',col(hal,:));
plot3(Af(9,1),Af(9,2),Af(9,3),"o",'Color','black','MarkerSize',17,...
    'MarkerFaceColor',col(3,:));
% text(Af(6,1) + 0.5, Af(6,2), Af(6,3)+ 1.0, 'C_1','FontSize', 15,'Color', 'black');
% text(Af(7,1)-2.7, Af(7,2)- 0.8, Af(7,3), 'C_2','FontSize', 15, 'Color', 'black');
% text(Af(8,1) + 1, Af(8,2)+ 0.5, Af(8,3)+ 0.8, 'X','FontSize', 15, 'Color', 'black');
% text(Af(2,1) , Af(2,2)+0.5, Af(2,3), 'H_a','FontSize', 15,'Color', 'black');
% text(Af(3,1) + 1, Af(3,2)+0.2, Af(3,3)+ 1.0, 'H_b','FontSize', 15,'Color', 'black');
% text(Af(5,1) + 1, Af(5,2), Af(5,3)+ 1.0, 'H_c','FontSize', 15,'Color', 'black');

% plot3([Af(1,1) Af(6,1)],[Af(1,2) Af(6,2)],[Af(1,3) Af(6,3)],'-','Color','black','LineWidth',lw);
% plot3([Af(2,1) Af(6,1)],[Af(2,2) Af(6,2)],[Af(2,3) Af(6,3)],'-','Color','black','LineWidth',lw);
% plot3([Af(3,1) Af(7,1)],[Af(3,2) Af(7,2)],[Af(3,3) Af(7,3)],'-','Color','black','LineWidth',lw);
% plot3([Af(4,1) Af(7,1)],[Af(4,2) Af(7,2)],[Af(4,3) Af(7,3)],'-','Color','black','LineWidth',lw);
% plot3([Af(5,1) Af(7,1)],[Af(5,2) Af(7,2)],[Af(5,3) Af(7,3)],'-','Color','black','LineWidth',lw);
% plot3([Af(8,1) Af(6,1)],[Af(8,2) Af(6,2)],[Af(8,3) Af(6,3)],'-','Color','black','LineWidth',lw);
% plot3([Af(6,1) Af(7,1)],[Af(6,2) Af(7,2)],[Af(6,3) Af(7,3)],'-','Color','black','LineWidth',lw);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot IC2H5X NTOs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot3([Af(1,1) Af(6,1)],[Af(1,2) Af(6,2)],[Af(1,3) Af(6,3)],'-','Color','black','LineWidth',lw);
% plot3([Af(2,1) Af(6,1)],[Af(2,2) Af(6,2)],[Af(2,3) Af(6,3)],'-','Color','black','LineWidth',lw);
% plot3([Af(3,1) Af(6,1)],[Af(3,2) Af(6,2)],[Af(3,3) Af(6,3)],'-','Color','black','LineWidth',lw);
% plot3([Af(4,1) Af(7,1)],[Af(4,2) Af(7,2)],[Af(4,3) Af(7,3)],'-','Color','black','LineWidth',lw);
% plot3([Af(5,1) Af(7,1)],[Af(5,2) Af(7,2)],[Af(5,3) Af(7,3)],'-','Color','black','LineWidth',lw);
% plot3([Af(8,1) Af(7,1)],[Af(8,2) Af(7,2)],[Af(8,3) Af(7,3)],'-','Color','black','LineWidth',lw);
% plot3([Af(6,1) Af(7,1)],[Af(6,2) Af(7,2)],[Af(6,3) Af(7,3)],'-','Color','black','LineWidth',lw);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if orb == 0 || orb == 2
% p1=quiver3(0,0,0,(vecz(1,1)/2),(vecz(2,1)/2),(vecz(3,1)/2));
% p1m=quiver3(0,0,0,-(vecz(1,1)/2),-(vecz(2,1)/2),-(vecz(3,1)/2));
% p2m=quiver3(0,0,0,-(vecx(1,1)/2),-(vecx(2,1)/2),-(vecx(3,1)/2));
% p3m=quiver3(0,0,0,-(vecy(1,1)/2),-(vecy(2,1)/2),-(vecy(3,1)/2));
% p2=quiver3(0,0,0,(vecx(1,1)/2),(vecx(2,1)/2),(vecx(3,1)/2));
% p3=quiver3(0,0,0,(vecy(1,1)/2),(vecy(2,1)/2),(vecy(3,1)/2));
% set(p1m, 'Color', [0 0 0], 'LineWidth', 2, 'LineStyle', '--', 'ShowArrowHead', 'off');
% set(p2m, 'Color', [0 0 0], 'LineWidth', 2, 'LineStyle', '--', 'ShowArrowHead', 'off');
% set(p3m, 'Color', [0 0 0], 'LineWidth', 2, 'LineStyle', '--', 'ShowArrowHead', 'off');
% set(p1, 'Color', [0 0 0], 'LineWidth', 2, 'LineStyle', '-', 'ShowArrowHead', 'off');
% set(p2, 'Color', [0 0 0], 'LineWidth', 2, 'LineStyle', '-', 'ShowArrowHead', 'off');
% set(p3, 'Color', [0 0 0], 'LineWidth', 2, 'LineStyle', '-', 'ShowArrowHead', 'off');
% axis equal;
% text((-vecz(1,1)/2)+3,(vecz(2,1)/2),(vecz(3,1)/2),'z','FontSize',14);
% text((vecx(1,1)/3),(vecx(2,1)/3),(vecx(3,1)/3),'x','FontSize',14);
% text((vecy(1,1)/3),(vecy(2,1)/3),(vecy(3,1)/3),'y','FontSize',14);
end
end
% if hal == 1
% text(vecz(1,1),vecz(2,1),vecz(3,1)+1,'x','FontSize',12);
% text((vecx(1,1)),vecx(2,1),-(vecx(3,1)+1),'z','FontSize',12);
% text(vecy(1,1),vecy(2,1)+1,vecy(3,1),'y','FontSize',12);
% end
% end
% figure('Units', 'pixels', 'Position', [100, 100, 1200, 800]);
figure('Units', 'pixels', 'Position', [100, 100, 5250, 3938]);
%%%%%%%% set to get good axes figure %%%%%%
% set(gcf, 'GraphicsSmoothing', 'on', 'Renderer', 'painters');
% set(gcf, 'Units', 'inches', 'Position', [1, 1, 2, 2.5]);
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperPosition', [1, 1, 3.5, 3]);
%%%%%% save the files %%%%%%%%
baseName = filename1(1, 1:end-4); 
% if orb==0
%  fullName = [baseName, num2str(MO), '_label.tif']; % Add the .tif extension
% print(fig, fullName, '-dtiff', '-r600');    
% saveas(fig,strcat(filename1(1,1:(size(filename1,2)-4)),num2str(MO),'geo.tif'));
% end
if orb==1 
fullName = [baseName, num2str(MO), '_HF.tif']; % Add the .tif extension
% print(fig, fullName, '-dtiff', '-r900');    
% saveas(fig,strcat(filename1(1,1:(size(filename1,2)-4)),num2str(MO),'_HF.tif'));
end
% if orb==2 
% saveas(fig,strcat(filename1(1,1:(size(filename1,2)-4)),num2str(MO),'_HF.tif'));
% end
end
