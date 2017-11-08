%% get link
% to download the dataset
fid = fopen('structural_ncvx.txt');
C=textscan(fid,'%s');
fclose(fid);
text = C{1,1};
link = cell(216,1);
j=1;
for i=1:size(text,1)
    str = text{i};
    if length(str)>20 && strcmp(str(end-5:end),'.mat">')
        link{j}=str(7:end-2);
        system(['wget -P /home/yu/datasets/SuiteSparse/structural_problem/ncvx/ ' link{j}]);
        j=j+1;
    end
    
end
