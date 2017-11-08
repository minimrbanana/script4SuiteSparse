function test_suite_cvx1(matClass)
%% test_cvxbqp1
% with amd and rcm, rand and randn of b
%% set path
addpath ../../BCD
dataPath=['/home/yu/datasets/SuiteSparse/' matClass '/cvx/'];
savePath=['/home/yu/bcd/SuiteSparse/result/' matClass '/cvx/'];
if ~exist(dataPath,'dir')
    error('matClass input error');
end
matList = dir(dataPath);
matList(1:2)=[];
% refine matList
%% set parameters
l=0;
u=1;
precision=1E-4;
iters = 20000;
%%
rng(1);
for i=1:size(matList,1)
    %% load data and reordering
    load ([dataPath matList(i).name]);
    saveDir = [savePath matList(i).name(1:end-4) '/'];
    if ~exist(saveDir,'dir')
        mkdir(saveDir);
    end
    if ~exist([saveDir 'b1/'],'dir')
        mkdir([saveDir 'b1/']);
    end
    A=Problem.A;
    d = size(A,1);
    % set b
    x=rand(d,1);
    x1=0.8*x+0.1;
    b=A*x1;
    % first do amd reordering
    tstart=tic;
    p1=symamd(A);
    A_amd=A(p1,p1);
    b_amd=b(p1);
    tAMD=toc(tstart);
    % second rcm reordering
    tstart=tic;
    p2=symrcm(A);
    A_rcm=A(p2,p2);
    b_rcm=b(p2);
    tRCM=toc(tstart);
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
    P3 = spdiags([e2,e1,e,[0;e1(1:end-1)],[0;0;e2(1:end-2)]],...
        [-2,-1,0,1,2],d,d);
    COVER=zeros(3,2);
    COVER(1,1)=nnz(A.*P2)/nnz(P2);%ratioA2
    COVER(1,2)=nnz(A.*P3)/nnz(P3);%ratioA3
    COVER(2,1)=nnz(A_amd.*P2)/nnz(P2);%ratioAMD2
    COVER(2,2)=nnz(A_amd.*P3)/nnz(P3);%ratioAMD3
    COVER(3,1)=nnz(A_rcm.*P2)/nnz(P2);%ratioRCM2
    COVER(3,2)=nnz(A_rcm.*P3)/nnz(P3);%ratioRCM3
    %% compute 0, 0.5 and 1 KKT condition
    % init=0
    x0=zeros(d,1);grad=x0;
    index_l = find(x0<=l+2*eps);
    index_u = find(x0>=u-2*eps);
    index = find(x0>l+2*eps & x0<u-2*eps);
    KKT0 = norm([grad(index)-b(index);min(0,grad(index_l)-b(index_l));...
        max(0,grad(index_u)-b(index_u))],2);
    % init=1
    x1=ones(d,1);grad=A*x1;
    index_l = find(x1<=l+2*eps);
    index_u = find(x1>=u-2*eps);
    index = find(x1>l+2*eps & x1<u-2*eps);
    KKT1 = norm([grad(index)-b(index);min(0,grad(index_l)-b(index_l));...
        max(0,grad(index_u)-b(index_u))],2);
    % init=0.5
    x5=ones(d,1)/2;grad=A*x5;
    index_l = find(x5<=l+2*eps);
    index_u = find(x5>=u-2*eps);
    index = find(x5>l+2*eps & x5<u-2*eps);
    KKT5 = norm([grad(index)-b(index);min(0,grad(index_l)-b(index_l));...
        max(0,grad(index_u)-b(index_u))],2);
    %% first run with b-> [0.1,0.9] compare block size123
    %% run original matrix with init 0
    init=l;
    t0=tic;
    [c0x1, c0y1] = CBCD_size1_gc(A, b, d, iters,precision,l,u,init);
    tC01=toc(t0);t0=tic;
    [~, c0y2] = CBCD_size2_gc(A, b, d, iters,precision,l,u,init);
    tC02=toc(t0);t0=tic;
    [~, c0y3] = CBCD_size3_gc(A, b, d, iters,precision,l,u,init);
    tC03=toc(t0);t0=tic;
    [~, r0y3] = RBCD3(A, b, d, iters,precision,l,u,init,1.0);
    tR03=toc(t0);t0=tic;
    [~, r0y2] = RBCD2(A, b, d, iters,precision,l,u,init,1.0);
    tR02=toc(t0);t0=tic;
    [~, r0y1] = RBCD_size1_gc_u(A, b, d, iters,precision,l,u,init,1.0);
    tR01=toc(t0);
    % plot original matrix with init 0
    figure(1),
    clf;
    semilogy(0:size(c0y1,1)-1,c0y1,'r','LineWidth',2.5);hold on;
    semilogy(0:size(c0y2,1)-1,c0y2,'g','LineWidth',2.5);hold on;
    semilogy(0:size(c0y3,1)-1,c0y3,'b','LineWidth',2.5);hold on;
    semilogy(0:size(r0y1,1)-1,r0y1,'r--','LineWidth',2.5);hold on;
    semilogy(0:size(r0y2,1)-1,r0y2,'g--','LineWidth',2.5);hold on;
    semilogy(0:size(r0y3,1)-1,r0y3,'b--','LineWidth',2.5);hold on;
    l1=sprintf('CBCD1, %.4f s, #%d',tC01,size(c0y1,1)-1);
    l2=sprintf('CBCD2, %.4f s, #%d',tC02,size(c0y2,1)-1);
    l3=sprintf('CBCD3, %.4f s, #%d',tC03,size(c0y3,1)-1);
    l4=sprintf('RBCD1, %.4f s, #%d',tR01,size(r0y1,1)-1);
    l5=sprintf('RBCD2, %.4f s, #%d',tR02,size(r0y2,1)-1);
    l6=sprintf('RBCD3, %.4f s, #%d',tR03,size(r0y3,1)-1);
    legend(l1,l2,l3,l4,l5,l6);
    grid on;
    set(gca,'fontsize',14);
    xlabel('#epoch');ylabel('KKT Condition');
    title(['Convergence before Reordering, KKT0=' num2str(KKT0)]);
    saveas(gca,[saveDir 'b1/figure_A0.png']);
    %% run original matrix with init 1
    init=u;
    t0=tic;
    [c1x1, c1y1] = CBCD_size1_gc(A, b, d, iters,precision,l,u,init);
    tC01=toc(t0);t0=tic;
    [~, c1y2] = CBCD_size2_gc(A, b, d, iters,precision,l,u,init);
    tC02=toc(t0);t0=tic;
    [~, c1y3] = CBCD_size3_gc(A, b, d, iters,precision,l,u,init);
    tC03=toc(t0);t0=tic;
    [~, r1y3] = RBCD3(A, b, d, iters,precision,l,u,init,1.0);
    tR03=toc(t0);t0=tic;
    [~, r1y2] = RBCD2(A, b, d, iters,precision,l,u,init,1.0);
    tR02=toc(t0);t0=tic;
    [~, r1y1] = RBCD_size1_gc_u(A, b, d, iters,precision,l,u,init,1.0);
    tR01=toc(t0);
    % plot original matrix with init 1
    figure(2),
    clf;
    semilogy(0:size(c1y1,1)-1,c1y1,'r','LineWidth',2.5);hold on;
    semilogy(0:size(c1y2,1)-1,c1y2,'g','LineWidth',2.5);hold on;
    semilogy(0:size(c1y3,1)-1,c1y3,'b','LineWidth',2.5);hold on;
    semilogy(0:size(r1y1,1)-1,r1y1,'r--','LineWidth',2.5);hold on;
    semilogy(0:size(r1y2,1)-1,r1y2,'g--','LineWidth',2.5);hold on;
    semilogy(0:size(r1y3,1)-1,r1y3,'b--','LineWidth',2.5);hold on;
    l1=sprintf('CBCD1, %.4f s, #%d',tC01,size(c1y1,1)-1);
    l2=sprintf('CBCD2, %.4f s, #%d',tC02,size(c1y2,1)-1);
    l3=sprintf('CBCD3, %.4f s, #%d',tC03,size(c1y3,1)-1);
    l4=sprintf('RBCD1, %.4f s, #%d',tR01,size(r1y1,1)-1);
    l5=sprintf('RBCD2, %.4f s, #%d',tR02,size(r1y2,1)-1);
    l6=sprintf('RBCD3, %.4f s, #%d',tR03,size(r1y3,1)-1);
    legend(l1,l2,l3,l4,l5,l6);
    grid on;
    set(gca,'fontsize',14);
    xlabel('#epoch');ylabel('KKT Condition');
    title(['Convergence before Reordering, KKT1=' num2str(KKT1)]);
    saveas(gca,[saveDir 'b1/figure_A1.png']);
    %% run original matrix with init 0.5
    init=(l+u)/2;
    t0=tic;
    [c5x1, c5y1] = CBCD_size1_gc(A, b, d, iters,precision,l,u,init);
    tC01=toc(t0);t0=tic;
    [~, c5y2] = CBCD_size2_gc(A, b, d, iters,precision,l,u,init);
    tC02=toc(t0);t0=tic;
    [~, c5y3] = CBCD_size3_gc(A, b, d, iters,precision,l,u,init);
    tC03=toc(t0);t0=tic;
    [~, r5y3] = RBCD3(A, b, d, iters,precision,l,u,init,1.0);
    tR03=toc(t0);t0=tic;
    [~, r5y2] = RBCD2(A, b, d, iters,precision,l,u,init,1.0);
    tR02=toc(t0);t0=tic;
    [~, r5y1] = RBCD_size1_gc_u(A, b, d, iters,precision,l,u,init,1.0);
    tR01=toc(t0);
    % plot original matrix with init 0.5
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
    set(gca,'fontsize',14);
    xlabel('#epoch');ylabel('KKT Condition');
    title(['Convergence before Reordering, KKT0.5=' num2str(KKT5)]);
    saveas(gca,[saveDir 'b1/figure_A5.png']);
    %% run matrix after amd with init 0
    init=l;
    t0=tic;
    [c0ax1, c0ay1] = CBCD_size1_gc(A_amd, b_amd, d, iters,precision,l,u,init);
    tC11=toc(t0);t0=tic;
    [~, c0ay2] = CBCD_size2_gc(A_amd, b_amd, d, iters,precision,l,u,init);
    tC12=toc(t0);t0=tic;
    [~, c0ay3] = CBCD_size3_gc(A_amd, b_amd, d, iters,precision,l,u,init);
    tC13=toc(t0);t0=tic;
    [~, r0ay3] = RBCD3(A_amd, b_amd, d, iters,precision,l,u,init,1.0);
    tR13=toc(t0);t0=tic;
    [~, r0ay2] = RBCD2(A_amd, b_amd, d, iters,precision,l,u,init,1.0);
    tR12=toc(t0);t0=tic;
    [~, r0ay1] = RBCD_size1_gc_u(A_amd, b_amd, d, iters,precision,l,u,init,1.0);
    tR11=toc(t0);
    % plot original matrix with init 0
    figure(4),
    clf;
    semilogy(0:size(c0ay1,1)-1,c0ay1,'r','LineWidth',2.5);hold on;
    semilogy(0:size(c0ay2,1)-1,c0ay2,'g','LineWidth',2.5);hold on;
    semilogy(0:size(c0ay3,1)-1,c0ay3,'b','LineWidth',2.5);hold on;
    semilogy(0:size(r0ay1,1)-1,r0ay1,'r--','LineWidth',2.5);hold on;
    semilogy(0:size(r0ay2,1)-1,r0ay2,'g--','LineWidth',2.5);hold on;
    semilogy(0:size(r0ay3,1)-1,r0ay3,'b--','LineWidth',2.5);hold on;
    l1=sprintf('CBCD1, %.4f s, #%d',tC11,size(c0ay1,1)-1);
    l2=sprintf('CBCD2, %.4f s, #%d',tC12,size(c0ay2,1)-1);
    l3=sprintf('CBCD3, %.4f s, #%d',tC13,size(c0ay3,1)-1);
    l4=sprintf('RBCD1, %.4f s, #%d',tR11,size(r0ay1,1)-1);
    l5=sprintf('RBCD2, %.4f s, #%d',tR12,size(r0ay2,1)-1);
    l6=sprintf('RBCD3, %.4f s, #%d',tR13,size(r0ay3,1)-1);
    legend(l1,l2,l3,l4,l5,l6);
    grid on;
    set(gca,'fontsize',14);
    xlabel('#epoch');ylabel('KKT Condition');
    title(['Convergence after AMD, AMD time=' num2str(tAMD)]);
    saveas(gca,[saveDir 'b1/figure_A0amd.png']);
    %% run matrix after amd with init 1
    init=u;
    t0=tic;
    [c1ax1, c1ay1] = CBCD_size1_gc(A_amd, b_amd, d, iters,precision,l,u,init);
    tC11=toc(t0);t0=tic;
    [~, c1ay2] = CBCD_size2_gc(A_amd, b_amd, d, iters,precision,l,u,init);
    tC12=toc(t0);t0=tic;
    [~, c1ay3] = CBCD_size3_gc(A_amd, b_amd, d, iters,precision,l,u,init);
    tC13=toc(t0);t0=tic;
    [~, r1ay3] = RBCD3(A_amd, b_amd, d, iters,precision,l,u,init,1.0);
    tR13=toc(t0);t0=tic;
    [~, r1ay2] = RBCD2(A_amd, b_amd, d, iters,precision,l,u,init,1.0);
    tR12=toc(t0);t0=tic;
    [~, r1ay1] = RBCD_size1_gc_u(A_amd, b_amd, d, iters,precision,l,u,init,1.0);
    tR11=toc(t0);
    % plot original matrix with init 1
    figure(5),
    clf;
    semilogy(0:size(c1ay1,1)-1,c1ay1,'r','LineWidth',2.5);hold on;
    semilogy(0:size(c1ay2,1)-1,c1ay2,'g','LineWidth',2.5);hold on;
    semilogy(0:size(c1ay3,1)-1,c1ay3,'b','LineWidth',2.5);hold on;
    semilogy(0:size(r1ay1,1)-1,r1ay1,'r--','LineWidth',2.5);hold on;
    semilogy(0:size(r1ay2,1)-1,r1ay2,'g--','LineWidth',2.5);hold on;
    semilogy(0:size(r1ay3,1)-1,r1ay3,'b--','LineWidth',2.5);hold on;
    l1=sprintf('CBCD1, %.4f s, #%d',tC11,size(c1ay1,1)-1);
    l2=sprintf('CBCD2, %.4f s, #%d',tC12,size(c1ay2,1)-1);
    l3=sprintf('CBCD3, %.4f s, #%d',tC13,size(c1ay3,1)-1);
    l4=sprintf('RBCD1, %.4f s, #%d',tR11,size(r1ay1,1)-1);
    l5=sprintf('RBCD2, %.4f s, #%d',tR12,size(r1ay2,1)-1);
    l6=sprintf('RBCD3, %.4f s, #%d',tR13,size(r1ay3,1)-1);
    legend(l1,l2,l3,l4,l5,l6);
    grid on;
    set(gca,'fontsize',14);
    xlabel('#epoch');ylabel('KKT Condition');
    title(['Convergence after AMD, AMD time=' num2str(tAMD)]);
    saveas(gca,[saveDir 'b1/figure_A1amd.png']);
    %% run matrix after amd with init 0.5
    init=(u+l)/2;
    t0=tic;
    [c5ax1, c5ay1] = CBCD_size1_gc(A_amd, b_amd, d, iters,precision,l,u,init);
    tC11=toc(t0);t0=tic;
    [~, c5ay2] = CBCD_size2_gc(A_amd, b_amd, d, iters,precision,l,u,init);
    tC12=toc(t0);t0=tic;
    [~, c5ay3] = CBCD_size3_gc(A_amd, b_amd, d, iters,precision,l,u,init);
    tC13=toc(t0);t0=tic;
    [~, r5ay3] = RBCD3(A_amd, b_amd, d, iters,precision,l,u,init,1.0);
    tR13=toc(t0);t0=tic;
    [~, r5ay2] = RBCD2(A_amd, b_amd, d, iters,precision,l,u,init,1.0);
    tR12=toc(t0);t0=tic;
    [~, r5ay1] = RBCD_size1_gc_u(A_amd, b_amd, d, iters,precision,l,u,init,1.0);
    tR11=toc(t0);
    % plot original matrix with init 0.5
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
    set(gca,'fontsize',14);
    xlabel('#epoch');ylabel('KKT Condition');
    title(['Convergence after AMD, AMD time=' num2str(tAMD)]);
    saveas(gca,[saveDir 'b1/figure_A5amd.png']);
    %% run matrix after rcm with init 0
    init=l;
    t0=tic;
    [c0rx1, c0ry1] = CBCD_size1_gc(A_rcm, b_rcm, d, iters,precision,l,u,init);
    tC21=toc(t0);t0=tic;
    [~, c0ry2] = CBCD_size2_gc(A_rcm, b_rcm, d, iters,precision,l,u,init);
    tC22=toc(t0);t0=tic;
    [~, c0ry3] = CBCD_size3_gc(A_rcm, b_rcm, d, iters,precision,l,u,init);
    tC23=toc(t0);t0=tic;
    [~, r0ry3] = RBCD3(A_rcm, b_rcm, d, iters,precision,l,u,init,1.0);
    tR23=toc(t0);t0=tic;
    [~, r0ry2] = RBCD2(A_rcm, b_rcm, d, iters,precision,l,u,init,1.0);
    tR22=toc(t0);t0=tic;
    [~, r0ry1] = RBCD_size1_gc_u(A_rcm, b_rcm, d, iters,precision,l,u,init,1.0);
    tR21=toc(t0);
    % plot original matrix with init 0
    figure(7),
    clf;
    semilogy(0:size(c0ry1,1)-1,c0ry1,'r','LineWidth',2.5);hold on;
    semilogy(0:size(c0ry2,1)-1,c0ry2,'g','LineWidth',2.5);hold on;
    semilogy(0:size(c0ry3,1)-1,c0ry3,'b','LineWidth',2.5);hold on;
    semilogy(0:size(r0ry1,1)-1,r0ry1,'r--','LineWidth',2.5);hold on;
    semilogy(0:size(r0ry2,1)-1,r0ry2,'g--','LineWidth',2.5);hold on;
    semilogy(0:size(r0ry3,1)-1,r0ry3,'b--','LineWidth',2.5);hold on;
    l1=sprintf('CBCD1, %.4f s, #%d',tC21,size(c0ry1,1)-1);
    l2=sprintf('CBCD2, %.4f s, #%d',tC22,size(c0ry2,1)-1);
    l3=sprintf('CBCD3, %.4f s, #%d',tC23,size(c0ry3,1)-1);
    l4=sprintf('RBCD1, %.4f s, #%d',tR21,size(r0ry1,1)-1);
    l5=sprintf('RBCD2, %.4f s, #%d',tR22,size(r0ry2,1)-1);
    l6=sprintf('RBCD3, %.4f s, #%d',tR23,size(r0ry3,1)-1);
    legend(l1,l2,l3,l4,l5,l6);
    grid on;
    set(gca,'fontsize',14);
    xlabel('#epoch');ylabel('KKT Condition');
    title(['Convergence after RCM; RCM time=' num2str(tRCM)]);
    saveas(gca,[saveDir 'b1/figure_A0rcm.png']);
    %% run matrix after rcm with init 1
    init=u;
    t0=tic;
    [c1rx1, c1ry1] = CBCD_size1_gc(A_rcm, b_rcm, d, iters,precision,l,u,init);
    tC21=toc(t0);t0=tic;
    [~, c1ry2] = CBCD_size2_gc(A_rcm, b_rcm, d, iters,precision,l,u,init);
    tC22=toc(t0);t0=tic;
    [~, c1ry3] = CBCD_size3_gc(A_rcm, b_rcm, d, iters,precision,l,u,init);
    tC23=toc(t0);t0=tic;
    [~, r1ry3] = RBCD3(A_rcm, b_rcm, d, iters,precision,l,u,init,1.0);
    tR23=toc(t0);t0=tic;
    [~, r1ry2] = RBCD2(A_rcm, b_rcm, d, iters,precision,l,u,init,1.0);
    tR22=toc(t0);t0=tic;
    [~, r1ry1] = RBCD_size1_gc_u(A_rcm, b_rcm, d, iters,precision,l,u,init,1.0);
    tR21=toc(t0);
    % plot original matrix with init 1
    figure(8),
    clf;
    semilogy(0:size(c1ry1,1)-1,c1ry1,'r','LineWidth',2.5);hold on;
    semilogy(0:size(c1ry2,1)-1,c1ry2,'g','LineWidth',2.5);hold on;
    semilogy(0:size(c1ry3,1)-1,c1ry3,'b','LineWidth',2.5);hold on;
    semilogy(0:size(r1ry1,1)-1,r1ry1,'r--','LineWidth',2.5);hold on;
    semilogy(0:size(r1ry2,1)-1,r1ry2,'g--','LineWidth',2.5);hold on;
    semilogy(0:size(r1ry3,1)-1,r1ry3,'b--','LineWidth',2.5);hold on;
    l1=sprintf('CBCD1, %.4f s, #%d',tC21,size(c1ry1,1)-1);
    l2=sprintf('CBCD2, %.4f s, #%d',tC22,size(c1ry2,1)-1);
    l3=sprintf('CBCD3, %.4f s, #%d',tC23,size(c1ry3,1)-1);
    l4=sprintf('RBCD1, %.4f s, #%d',tR21,size(r1ry1,1)-1);
    l5=sprintf('RBCD2, %.4f s, #%d',tR22,size(r1ry2,1)-1);
    l6=sprintf('RBCD3, %.4f s, #%d',tR23,size(r1ry3,1)-1);
    legend(l1,l2,l3,l4,l5,l6);
    grid on;
    set(gca,'fontsize',14);
    xlabel('#epoch');ylabel('KKT Condition');
    title(['Convergence after RCM; RCM time=' num2str(tRCM)]);
    saveas(gca,[saveDir 'b1/figure_A1rcm.png']);    
    %% run matrix after rcm with init 0.5
    init=(u+l)/2;
    t0=tic;
    [c5rx1, c5ry1] = CBCD_size1_gc(A_rcm, b_rcm, d, iters,precision,l,u,init);
    tC21=toc(t0);t0=tic;
    [~, c5ry2] = CBCD_size2_gc(A_rcm, b_rcm, d, iters,precision,l,u,init);
    tC22=toc(t0);t0=tic;
    [~, c5ry3] = CBCD_size3_gc(A_rcm, b_rcm, d, iters,precision,l,u,init);
    tC23=toc(t0);t0=tic;
    [~, r5ry3] = RBCD3(A_rcm, b_rcm, d, iters,precision,l,u,init,1.0);
    tR23=toc(t0);t0=tic;
    [~, r5ry2] = RBCD2(A_rcm, b_rcm, d, iters,precision,l,u,init,1.0);
    tR22=toc(t0);t0=tic;
    [~, r5ry1] = RBCD_size1_gc_u(A_rcm, b_rcm, d, iters,precision,l,u,init,1.0);
    tR21=toc(t0);
    % plot original matrix with init 1
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
    set(gca,'fontsize',14);
    xlabel('#epoch');ylabel('KKT Condition');
    title(['Convergence after RCM; RCM time=' num2str(tRCM)]);
    saveas(gca,[saveDir 'b1/figure_A5rcm.png']); 



    %% plot the structure of matrix and save the results
    title1=sprintf('NNZ of A, nnz in box:%.3f & %.3f',COVER(1,1),COVER(1,2));
    figure(100),spy(A);
    title(title1);
    saveas(gca,[saveDir 'figure_A.png']);
    title2=sprintf('NNZ of A after AMD, nnz in box:%.3f & %.3f',...
        COVER(2,1),COVER(2,2));
    figure(101),spy(A_amd);
    title(title2);
    saveas(gca,[saveDir 'figure_AMD.png']);
    title3=sprintf('NNZ of A after RCM, nnz in box:%.3f & %.3f',...
        COVER(3,1),COVER(3,2));
    figure(102),spy(A_rcm);
    title(title3);
    saveas(gca,[saveDir 'figure_RCM.png']);
    
    x_opt=zeros(d,9);
    x_opt(:,1)=c0x1;x_opt(:,2)=c1x1;x_opt(:,3)=c5x1;
    x_opt(:,4)=c0ax1;x_opt(:,5)=c1ax1;x_opt(:,6)=c5ax1;
    x_opt(:,7)=c0rx1;x_opt(:,8)=c1rx1;x_opt(:,9)=c5rx1;
    save([saveDir 'b1/exp.mat'],'KKT0','KKT1','KKT5','COVER','x_opt');

end




end
