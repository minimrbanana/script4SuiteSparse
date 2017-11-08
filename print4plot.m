function print4plot(m,matClass,N)
%% read the results saved and plot the figure
%  m=2 or 3 means b2 or b3

% cvx and ncvx
dataPath=['/home/yu/datasets/SuiteSparse/' matClass '/ncvx/'];
%dataPath=['/home/yu/datasets/SuiteSparse/' matClass '/cvx/'];
matlist=dir(dataPath);
matlist(1:2)=[];

% cvx and ncvx
saveDir=['/home/yu/bcd/SuiteSparse/result/' matClass '/ncvx/' matlist(N).name(1:end-4) '/b' num2str(m) '/'];
%saveDir=['/home/yu/bcd/SuiteSparse/result/' matClass '/cvx/' matlist(N).name(1:end-4) '/b' num2str(m) '/'];

if ~exist(dataPath,'dir')
    error('matClass input error');
end

if m~=4
% before reorder
    load ([saveDir 'figure_A5.mat']);
    figure(3),
    clf;
    semilogy(0:size(c5y1,1)-1,c5y1,'r','LineWidth',2.5);hold on;
    semilogy(0:size(c5y2,1)-1,c5y2,'g','LineWidth',2.5);hold on;
    semilogy(0:size(c5y3,1)-1,c5y3,'b','LineWidth',2.5);hold on;
    semilogy(0:size(r5y1,1)-1,r5y1,'r--','LineWidth',2.5);hold on;
    semilogy(0:size(r5y2,1)-1,r5y2,'g--','LineWidth',2.5);hold on;
    semilogy(0:size(r5y3,1)-1,r5y3,'b--','LineWidth',2.5);hold on;
    l1=sprintf('CBCD1, %.4f s, #%d',tC01,size(c5y1,1)-1);
    l2=sprintf('CBCD2, %.4f s, #%d',tC02,size(c5y2,1)-1);
    l3=sprintf('CBCD3, %.4f s, #%d',tC03,size(c5y3,1)-1);
    l4=sprintf('RBCD1, %.4f s, #%d',tR01,size(r5y1,1)-1);
    l5=sprintf('RBCD2, %.4f s, #%d',tR02,size(r5y2,1)-1);
    l6=sprintf('RBCD3, %.4f s, #%d',tR03,size(r5y3,1)-1);
    legend(l1,l2,l3,l4,l5,l6);
    grid on;
    axis([0 max([size(r5y1,1),size(r5y2,1),size(r5y3,1)]) 1E-10 1E5]);
    % setting
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    
    xlabel('#epoch');ylabel('KKT Condition');
    title('Before Reordering');
    set(gca,'fontsize',20,'fontweight', 'bold');
    saveas(gca,[saveDir 'figure_A5.png']);
    saveas(gca,[saveDir 'figure_A5.pdf']);

    disp('A5 fertig');
% amd
    load ([saveDir 'figure_A5amd.mat']);
    figure(6),
    clf;
    semilogy(0:size(c5ay1,1)-1,c5ay1,'r','LineWidth',2.5);hold on;
    semilogy(0:size(c5ay2,1)-1,c5ay2,'g','LineWidth',2.5);hold on;
    semilogy(0:size(c5ay3,1)-1,c5ay3,'b','LineWidth',2.5);hold on;
    semilogy(0:size(r5ay1,1)-1,r5ay1,'r--','LineWidth',2.5);hold on;
    semilogy(0:size(r5ay2,1)-1,r5ay2,'g--','LineWidth',2.5);hold on;
    semilogy(0:size(r5ay3,1)-1,r5ay3,'b--','LineWidth',2.5);hold on;
    l1=sprintf('CBCD1, %.4f s, #%d',tC11,size(c5ay1,1)-1);
    l2=sprintf('CBCD2, %.4f s, #%d',tC12,size(c5ay2,1)-1);
    l3=sprintf('CBCD3, %.4f s, #%d',tC13,size(c5ay3,1)-1);
    l4=sprintf('RBCD1, %.4f s, #%d',tR11,size(r5ay1,1)-1);
    l5=sprintf('RBCD2, %.4f s, #%d',tR12,size(r5ay2,1)-1);
    l6=sprintf('RBCD3, %.4f s, #%d',tR13,size(r5ay3,1)-1);
    legend(l1,l2,l3,l4,l5,l6);
    grid on;
    axis([0 max([size(r5ay1,1),size(r5ay2,1),size(r5ay3,1)]) 1E-10 1E5]);
    % setting
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    
    xlabel('#epoch');ylabel('KKT Condition');
    title(['After AMD, AMD time=' num2str(tAMD)]);
    set(gca,'fontsize',20,'fontweight', 'bold');
    saveas(gca,[saveDir 'figure_A5amd.png']);
    saveas(gca,[saveDir 'figure_A5amd.pdf']);

    disp('A5amd fertig');
% rcm
    load ([saveDir 'figure_A5rcm.mat']);
    figure(9),
    clf;
    semilogy(0:size(c5ry1,1)-1,c5ry1,'r','LineWidth',2.5);hold on;
    semilogy(0:size(c5ry2,1)-1,c5ry2,'g','LineWidth',2.5);hold on;
    semilogy(0:size(c5ry3,1)-1,c5ry3,'b','LineWidth',2.5);hold on;
    semilogy(0:size(r5ry1,1)-1,r5ry1,'r--','LineWidth',2.5);hold on;
    semilogy(0:size(r5ry2,1)-1,r5ry2,'g--','LineWidth',2.5);hold on;
    semilogy(0:size(r5ry3,1)-1,r5ry3,'b--','LineWidth',2.5);hold on;
    l1=sprintf('CBCD1, %.4f s, #%d',tC21,size(c5ry1,1)-1);
    l2=sprintf('CBCD2, %.4f s, #%d',tC22,size(c5ry2,1)-1);
    l3=sprintf('CBCD3, %.4f s, #%d',tC23,size(c5ry3,1)-1);
    l4=sprintf('RBCD1, %.4f s, #%d',tR21,size(r5ry1,1)-1);
    l5=sprintf('RBCD2, %.4f s, #%d',tR22,size(r5ry2,1)-1);
    l6=sprintf('RBCD3, %.4f s, #%d',tR23,size(r5ry3,1)-1);
    legend(l1,l2,l3,l4,l5,l6);
    grid on;
    axis([0 max([size(r5ry1,1),size(r5ry2,1),size(r5ry3,1)]) 1E-10 1E5]);
    % setting
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    
    xlabel('#epoch');ylabel('KKT Condition');
    title(['After RCM, RCM time=' num2str(tRCM)]);
    set(gca,'fontsize',20,'fontweight', 'bold');
    saveas(gca,[saveDir 'figure_A5rcm.png']); 
    saveas(gca,[saveDir 'figure_A5rcm.pdf']);

    disp('A5rcm fertig');


end

%% plot structure
if m==4
    load ([dataPath matlist(N).name]);
    A=Problem.A;
    d = size(A,1);
    D=spdiags(ones(d,1)*1E-5,0,d,d);
    A=A+D;
    % first amd reordering
    p1=symamd(A);
    A_amd=A(p1,p1);
    % second rcm reordering
    p2=symrcm(A);
    A_rcm=A(p2,p2);
    
    %% compute the cover of blocks
    % plate size 2
    e = ones(d,1);
    e0 = [1;0];
    e1 = e0(:,ones(ceil(d/2),1));
    e1 = reshape(e1 ,numel(e1),1);
    e1 = e1(1:d);
    P2 = spdiags([e1,e,[0;e1(1:end-1)]],[-1,0,1],d,d);
    % plate size 3
    e0 = [1;1;0];
    e1 = e0(:,ones(ceil(d/3),1));
    e1 = reshape(e1 ,numel(e1),1);
    e1 = e1(1:d);
    e0 = [1;0;0];
    e2 = e0(:,ones(ceil(d/3),1));
    e2 = reshape(e2 ,numel(e2),1);
    e2 = e2(1:d);
    e = e2*0+1;
    P3 = spdiags([e2,e1,e,[0;e1(1:end-1)],[0;0;e2(1:end-2)]],...
        [-2,-1,0,1,2],d,d);
    COVER=zeros(3,2);
    COVER(1,1)=nnz(A.*P2)/nnz(A);%ratioA2
    COVER(1,2)=nnz(A.*P3)/nnz(A);%ratioA3
    COVER(2,1)=nnz(A_amd.*P2)/nnz(A);%ratioAMD2
    COVER(2,2)=nnz(A_amd.*P3)/nnz(A);%ratioAMD3
    COVER(3,1)=nnz(A_rcm.*P2)/nnz(A);%ratioRCM2
    COVER(3,2)=nnz(A_rcm.*P3)/nnz(A);%ratioRCM3
    title1=sprintf('NNZ of A, nnz in box:%.3f & %.3f',COVER(1,1),COVER(1,2));
    figure(100),
    spy(A);
    title(title1,'Position',[d/2, -d/18]);
    set(gcf,'Position',[232  246  520  540]);
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    set(gca,'FontSize',15,'fontweight', 'bold');
    saveas(gca,[saveDir(1:end-3) 'figure_A.pdf']);
    title2=sprintf('After AMD, nnz in box:%.3f & %.3f',...
        COVER(2,1),COVER(2,2));
    figure(101),
    spy(A_amd);
    title(title2,'Position',[d/2, -d/18]);
    set(gcf,'Position',[232  246  520  540]);
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    set(gca,'FontSize',15,'fontweight', 'bold');
    saveas(gca,[saveDir(1:end-3) 'figure_AMD.pdf']);
    title3=sprintf('After RCM, nnz in box:%.3f & %.3f',...
        COVER(3,1),COVER(3,2));
    figure(102),
    spy(A_rcm);
    title(title3,'Position',[d/2, -d/18]);
    set(gcf,'Position',[232  246  520  540]);
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    set(gca,'FontSize',15,'fontweight', 'bold');
    saveas(gca,[saveDir(1:end-3) 'figure_RCM.pdf']);
end

end