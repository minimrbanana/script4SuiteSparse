function test_suite_ncvx(matClass,lambda,gamma)
%% test_cvxbqp1
% with amd and rcm, rand and randn of b
%% set path
addpath ../../BCD
dataPath=['/home/yu/datasets/SuiteSparse/' matClass '_problem/ncvx/'];
savePath=['/home/yu/bcd/SuiteSparse/result/' matClass '_problem/ncvx/'];
if ~exist(dataPath,'dir')
    error('matClass input error');
end
matList = dir(dataPath);
matList(1:2)=[];
% refine matList
%% set parameters
l=0;
u=1;
precision=1E-10;

%%
rng(1);
for i=1:size(matList,1)
    %% load data and reordering
    load ([dataPath matList(i).name]);
    saveDir = [savePath matList(i).name(1:end-4) '_b=' num2str(gamma) '_lambda=' num2str(lambda) '/'];
    if ~exist(saveDir,'dir')
        mkdir(saveDir);
    end
    A=Problem.A;
    d = size(A,1);
    A((A~=0))=-1;
    A=spdiags(zeros(d,1),0,A);
    diagonal = -sum(A)+lambda;
    diagonal(diagonal==0)=1;% if sum of row/colomn is 0, set diagonal as 1
    A = spdiags(diagonal',0,A);
    iters = 2000000;
    b = randn(d,1)*gamma;
    % first do amd reordering
    tstart=tic;
    p1=symamd(A);
    A1=A(p1,p1);
    b1=b(p1);
    tAMD=toc(tstart);
    % second rcm reordering
    tstart=tic;
    p2=symrcm(A);
    A2=A(p2,p2);
    b2=b(p2);
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
    COVER(2,1)=nnz(A1.*P2)/nnz(P2);%ratioAMD2
    COVER(2,2)=nnz(A1.*P3)/nnz(P3);%ratioAMD3
    COVER(3,1)=nnz(A2.*P2)/nnz(P2);%ratioRCM2
    COVER(3,2)=nnz(A2.*P3)/nnz(P3);%ratioRCM3
    %% run original matrix with init 0
    t0=tic;
    [c0x1, c0y1] = CBCD_size1_gc(A, b, d, iters,precision,l,u,l);
    tC01=toc(t0);t0=tic;
    [c0x2, c0y2] = CBCD_size2_gc(A, b, d, iters,precision,l,u,l);
    tC02=toc(t0);t0=tic;
    [c0x3, c0y3] = CBCD_size3_gc(A, b, d, iters,precision,l,u,l);
    tC03=toc(t0);t0=tic;
    [r0x3, r0y3] = RBCD3(A, b, d, iters,precision,l,u,l,1.0);
    tR03=toc(t0);t0=tic;
    [r0x2, r0y2] = RBCD2(A, b, d, iters,precision,l,u,l,1.0);
    tR02=toc(t0);t0=tic;
    [r0x1, r0y1] = RBCD_size1_gc_u(A, b, d, iters,precision,l,u,l,1.0);
    tR01=toc(t0);
    % plot original matrix with init 0
    figure(1),
    clf;
    semilogy(0:size(c0y1,1)-1,c0y1,'r','LineWidth',2.5);
    hold on;
    semilogy(0:size(c0y2,1)-1,c0y2,'g','LineWidth',2.5);
    hold on;
    semilogy(0:size(c0y3,1)-1,c0y3,'b','LineWidth',2.5);
    hold on;
    semilogy(0:size(r0y1,1)-1,r0y1,'r--','LineWidth',2.5);
    hold on;
    semilogy(0:size(r0y2,1)-1,r0y2,'g--','LineWidth',2.5);
    hold on;
    semilogy(0:size(r0y3,1)-1,r0y3,'b--','LineWidth',2.5);
    hold on;
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
    title('Convergence before Reordering');
    saveas(gca,[saveDir 'figure_A0.png']);
    %% run original matrix with init 1
    t0=tic;
    [c0x1, c0y1] = CBCD_size1_gc(A, b, d, iters,precision,l,u,u);
    tC01=toc(t0);t0=tic;
    [c0x2, c0y2] = CBCD_size2_gc(A, b, d, iters,precision,l,u,u);
    tC02=toc(t0);t0=tic;
    [c0x3, c0y3] = CBCD_size3_gc(A, b, d, iters,precision,l,u,u);
    tC03=toc(t0);t0=tic;
    [r0x3, r0y3] = RBCD3(A, b, d, iters,precision,l,u,u,1.0);
    tR03=toc(t0);t0=tic;
    [r0x2, r0y2] = RBCD2(A, b, d, iters,precision,l,u,u,1.0);
    tR02=toc(t0);t0=tic;
    [r0x1, r0y1] = RBCD_size1_gc_u(A, b, d, iters,precision,l,u,u,1.0);
    tR01=toc(t0);
    % plot original matrix with init 1
    figure(2),
    clf;
    semilogy(0:size(c0y1,1)-1,c0y1,'r','LineWidth',2.5);
    hold on;
    semilogy(0:size(c0y2,1)-1,c0y2,'g','LineWidth',2.5);
    hold on;
    semilogy(0:size(c0y3,1)-1,c0y3,'b','LineWidth',2.5);
    hold on;
    semilogy(0:size(r0y1,1)-1,r0y1,'r--','LineWidth',2.5);
    hold on;
    semilogy(0:size(r0y2,1)-1,r0y2,'g--','LineWidth',2.5);
    hold on;
    semilogy(0:size(r0y3,1)-1,r0y3,'b--','LineWidth',2.5);
    hold on;
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
    title('Convergence before Reordering');
    saveas(gca,[saveDir 'figure_A1.png']);
    %% run matrix after amd with init 0
    t0=tic;
    [c1x1, c1y1] = CBCD_size1_gc(A1, b1, d, iters,precision,l,u,l);
    tC11=toc(t0);t0=tic;
    [c1x2, c1y2] = CBCD_size2_gc(A1, b1, d, iters,precision,l,u,l);
    tC12=toc(t0);t0=tic;
    [c1x3, c1y3] = CBCD_size3_gc(A1, b1, d, iters,precision,l,u,l);
    tC13=toc(t0);t0=tic;
    [r1x3, r1y3] = RBCD3(A1, b1, d, iters,precision,l,u,l,1.0);
    tR13=toc(t0);t0=tic;
    [r1x2, r1y2] = RBCD2(A1, b1, d, iters,precision,l,u,l,1.0);
    tR12=toc(t0);t0=tic;
    [r1x1, r1y1] = RBCD_size1_gc_u(A1, b1, d, iters,precision,l,u,l,1.0);
    tR11=toc(t0);
    % plot original matrix with init 0
    figure(3),
    clf;
    semilogy(0:size(c1y1,1)-1,c1y1,'r','LineWidth',2.5);
    hold on;
    semilogy(0:size(c1y2,1)-1,c1y2,'g','LineWidth',2.5);
    hold on;
    semilogy(0:size(c1y3,1)-1,c1y3,'b','LineWidth',2.5);
    hold on;
    semilogy(0:size(r1y1,1)-1,r1y1,'r--','LineWidth',2.5);
    hold on;
    semilogy(0:size(r1y2,1)-1,r1y2,'g--','LineWidth',2.5);
    hold on;
    semilogy(0:size(r1y3,1)-1,r1y3,'b--','LineWidth',2.5);
    hold on;
    l1=sprintf('CBCD1, %.4f s, #%d',tC11,size(c1y1,1)-1);
    l2=sprintf('CBCD2, %.4f s, #%d',tC12,size(c1y2,1)-1);
    l3=sprintf('CBCD3, %.4f s, #%d',tC13,size(c1y3,1)-1);
    l4=sprintf('RBCD1, %.4f s, #%d',tR11,size(r1y1,1)-1);
    l5=sprintf('RBCD2, %.4f s, #%d',tR12,size(r1y2,1)-1);
    l6=sprintf('RBCD3, %.4f s, #%d',tR13,size(r1y3,1)-1);
    legend(l1,l2,l3,l4,l5,l6);
    grid on;
    set(gca,'fontsize',14);
    xlabel('#epoch');ylabel('KKT Condition');
    title(['Convergence after AMD, AMD time=' num2str(tAMD)]);
    saveas(gca,[saveDir 'figure_A0amd.png']);
    %% run matrix after amd with init 1
    t0=tic;
    [c1x1, c1y1] = CBCD_size1_gc(A1, b1, d, iters,precision,l,u,u);
    tC11=toc(t0);t0=tic;
    [c1x2, c1y2] = CBCD_size2_gc(A1, b1, d, iters,precision,l,u,u);
    tC12=toc(t0);t0=tic;
    [c1x3, c1y3] = CBCD_size3_gc(A1, b1, d, iters,precision,l,u,u);
    tC13=toc(t0);t0=tic;
    [r1x3, r1y3] = RBCD3(A1, b1, d, iters,precision,l,u,u,1.0);
    tR13=toc(t0);t0=tic;
    [r1x2, r1y2] = RBCD2(A1, b1, d, iters,precision,l,u,u,1.0);
    tR12=toc(t0);t0=tic;
    [r1x1, r1y1] = RBCD_size1_gc_u(A1, b1, d, iters,precision,l,u,u,1.0);
    tR11=toc(t0);
    % plot original matrix with init 1
    figure(4),
    clf;
    semilogy(0:size(c1y1,1)-1,c1y1,'r','LineWidth',2.5);
    hold on;
    semilogy(0:size(c1y2,1)-1,c1y2,'g','LineWidth',2.5);
    hold on;
    semilogy(0:size(c1y3,1)-1,c1y3,'b','LineWidth',2.5);
    hold on;
    semilogy(0:size(r1y1,1)-1,r1y1,'r--','LineWidth',2.5);
    hold on;
    semilogy(0:size(r1y2,1)-1,r1y2,'g--','LineWidth',2.5);
    hold on;
    semilogy(0:size(r1y3,1)-1,r1y3,'b--','LineWidth',2.5);
    hold on;
    l1=sprintf('CBCD1, %.4f s, #%d',tC11,size(c1y1,1)-1);
    l2=sprintf('CBCD2, %.4f s, #%d',tC12,size(c1y2,1)-1);
    l3=sprintf('CBCD3, %.4f s, #%d',tC13,size(c1y3,1)-1);
    l4=sprintf('RBCD1, %.4f s, #%d',tR11,size(r1y1,1)-1);
    l5=sprintf('RBCD2, %.4f s, #%d',tR12,size(r1y2,1)-1);
    l6=sprintf('RBCD3, %.4f s, #%d',tR13,size(r1y3,1)-1);
    legend(l1,l2,l3,l4,l5,l6);
    grid on;
    set(gca,'fontsize',14);
    xlabel('#epoch');ylabel('KKT Condition');
    title(['Convergence after AMD, AMD time=' num2str(tAMD)]);
    saveas(gca,[saveDir 'figure_A1amd.png']);
    %% run matrix after rcm with init 0
    t0=tic;
    [c2x1, c2y1] = CBCD_size1_gc(A2, b2, d, iters,precision,l,u,l);
    tC21=toc(t0);t0=tic;
    [c2x2, c2y2] = CBCD_size2_gc(A2, b2, d, iters,precision,l,u,l);
    tC22=toc(t0);t0=tic;
    [c2x3, c2y3] = CBCD_size3_gc(A2, b2, d, iters,precision,l,u,l);
    tC23=toc(t0);t0=tic;
    [r2x3, r2y3] = RBCD3(A2, b2, d, iters,precision,l,u,l,1.0);
    tR23=toc(t0);t0=tic;
    [r2x2, r2y2] = RBCD2(A2, b2, d, iters,precision,l,u,l,1.0);
    tR22=toc(t0);t0=tic;
    [r2x1, r2y1] = RBCD_size1_gc_u(A2, b2, d, iters,precision,l,u,l,1.0);
    tR21=toc(t0);
    % plot original matrix with init 0
    figure(5),
    clf;
    semilogy(0:size(c2y1,1)-1,c2y1,'r','LineWidth',2.5);
    hold on;
    semilogy(0:size(c2y2,1)-1,c2y2,'g','LineWidth',2.5);
    hold on;
    semilogy(0:size(c2y3,1)-1,c2y3,'b','LineWidth',2.5);
    hold on;
    semilogy(0:size(r2y1,1)-1,r2y1,'r--','LineWidth',2.5);
    hold on;
    semilogy(0:size(r2y2,1)-1,r2y2,'g--','LineWidth',2.5);
    hold on;
    semilogy(0:size(r2y3,1)-1,r2y3,'b--','LineWidth',2.5);
    hold on;
    l1=sprintf('CBCD1, %.4f s, #%d',tC21,size(c2y1,1)-1);
    l2=sprintf('CBCD2, %.4f s, #%d',tC22,size(c2y2,1)-1);
    l3=sprintf('CBCD3, %.4f s, #%d',tC23,size(c2y3,1)-1);
    l4=sprintf('RBCD1, %.4f s, #%d',tR21,size(r2y1,1)-1);
    l5=sprintf('RBCD2, %.4f s, #%d',tR22,size(r2y2,1)-1);
    l6=sprintf('RBCD3, %.4f s, #%d',tR23,size(r2y3,1)-1);
    legend(l1,l2,l3,l4,l5,l6);
    grid on;
    set(gca,'fontsize',14);
    xlabel('#epoch');ylabel('KKT Condition');
    title(['Convergence after RCM; RCM time=' num2str(tRCM)]);
    saveas(gca,[saveDir 'figure_A0rcm.png']);
    %% run matrix after rcm with init 1
    t0=tic;
    [c2x1, c2y1] = CBCD_size1_gc(A2, b2, d, iters,precision,l,u,u);
    tC21=toc(t0);t0=tic;
    [c2x2, c2y2] = CBCD_size2_gc(A2, b2, d, iters,precision,l,u,u);
    tC22=toc(t0);t0=tic;
    [c2x3, c2y3] = CBCD_size3_gc(A2, b2, d, iters,precision,l,u,u);
    tC23=toc(t0);t0=tic;
    [r2x3, r2y3] = RBCD3(A2, b2, d, iters,precision,l,u,u,1.0);
    tR23=toc(t0);t0=tic;
    [r2x2, r2y2] = RBCD2(A2, b2, d, iters,precision,l,u,u,1.0);
    tR22=toc(t0);t0=tic;
    [r2x1, r2y1] = RBCD_size1_gc_u(A2, b2, d, iters,precision,l,u,u,1.0);
    tR21=toc(t0);
    % plot original matrix with init 1
    figure(6),
    clf;
    semilogy(0:size(c2y1,1)-1,c2y1,'r','LineWidth',2.5);
    hold on;
    semilogy(0:size(c2y2,1)-1,c2y2,'g','LineWidth',2.5);
    hold on;
    semilogy(0:size(c2y3,1)-1,c2y3,'b','LineWidth',2.5);
    hold on;
    semilogy(0:size(r2y1,1)-1,r2y1,'r--','LineWidth',2.5);
    hold on;
    semilogy(0:size(r2y2,1)-1,r2y2,'g--','LineWidth',2.5);
    hold on;
    semilogy(0:size(r2y3,1)-1,r2y3,'b--','LineWidth',2.5);
    hold on;
    l1=sprintf('CBCD1, %.4f s, #%d',tC21,size(c2y1,1)-1);
    l2=sprintf('CBCD2, %.4f s, #%d',tC22,size(c2y2,1)-1);
    l3=sprintf('CBCD3, %.4f s, #%d',tC23,size(c2y3,1)-1);
    l4=sprintf('RBCD1, %.4f s, #%d',tR21,size(r2y1,1)-1);
    l5=sprintf('RBCD2, %.4f s, #%d',tR22,size(r2y2,1)-1);
    l6=sprintf('RBCD3, %.4f s, #%d',tR23,size(r2y3,1)-1);
    legend(l1,l2,l3,l4,l5,l6);
    grid on;
    set(gca,'fontsize',14);
    xlabel('#epoch');ylabel('KKT Condition');
    title(['Convergence after RCM; RCM time=' num2str(tRCM)]);
    saveas(gca,[saveDir 'figure_A1rcm.png']);    




    % plot the structure of matrix
    figure(7),spy(A);
    title('Non-zero Elements of A');
    saveas(gca,[saveDir 'figure_A.png']);
    figure(8),spy(A1);
    title('Non-zero Elements of A after AMD');
    saveas(gca,[saveDir 'figure_AMD.png']);
    figure(9),spy(A2);
    title('Non-zero Elements of A after RCM');
    saveas(gca,[saveDir 'figure_RCM.png']);
    
    % compute 0 and 1 KKT condition
    % init=0
    x0=zeros(d,1);grad=x0;
    index_l = find(x0<=l+2*eps);
    index_u = find(x0>=u-2*eps);
    index = find(x0>l+2*eps & x0<u-2*eps);
    KKT0 = norm([grad(index)-b(index);min(0,grad(index_l)-b(index_l));max(0,grad(index_u)-b(index_u))],2);
    % init=1
    x1=ones(d,1);grad=A*x1;
    index_l = find(x1<=l+2*eps);
    index_u = find(x1>=u-2*eps);
    index = find(x1>l+2*eps & x1<u-2*eps);
    KKT1 = norm([grad(index)-b(index);min(0,grad(index_l)-b(index_l));max(0,grad(index_u)-b(index_u))],2);
    %cond_A=condest(A);
    save([saveDir 'exp.mat'],'KKT0','KKT1','COVER');

end




end
