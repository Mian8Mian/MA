%Coding

clear;

%set file path
quantityoffiles=5;
filelist=dir('D:\Masterarbeit\Quelle\ACE_AUTO');
dest='D:\Masterarbeit\matlab\debug';  

%define block size, set block number
rownr=[1 1];        %[startrownr endrownr]
%columnsize=200;    %comment if data is the whole matrix
blocknumber=1;      %set 1 if data is the whole matrix
orderofconditionalentropy=9;
format long;

%select random files
qty=quantityoffiles;
% randf=randperm(numel(filelist),qty);
% 
% for ii=1:qty
%     while (randf(ii)==1)||(randf(ii)==2)
%         randf(ii)=randperm(numel(filelist),1);
%     end
% end
% 
% %copy the selected files to destination folder
% for kk=1:qty
%     source=fullfile(filelist(randf(kk)).folder,filelist(randf(kk)).name);
%     copyfile(source,dest);
% end

%get filelist of the current destination folder
cfiles=dir(dest);

%initialization for clone struct
clone=struct();

%initialization for blockp struct
blockp=struct();

%cfiles(1) and cfiles(2) are '.' and '..', file name is stored from cfiles(3)
jj=3;
while jj<=qty+2

    %load data file
    load(fullfile(cfiles(jj).folder,cfiles(jj).name));

    %initialization, convert NaN to 0
    clone(jj-2).file=cfiles(jj).name;
    clone(jj-2).mat=mat;

    %next data file
    jj=jj+1;
end

%cfiles(1) and cfiles(2) are '.' and '..', file name is stored from cfiles(3)
jj=3;
while jj<=qty+2
    
    %initialization, convert NaN to 0
    clone(jj-2).mat(isnan(clone(jj-2).mat))=0;

    %define block size and split the matrix {clone} into smaller blocks 
    %nob=number of blocks
    row=rownr;
     
    %column=columnsize;  %use this if data is only a part of the matrix
    buffer=size(mat);   %use this if data is the whole matrix
    column=buffer(2);   %use this if data is the whole matrix
    nob=fix(length(clone(jj-2).mat)/column);

    if nob~=length(clone(jj-2).mat)/column
        nob=nob+1;
    end
    
    %preallocate array(part) size
    part=zeros(nob,2);    
    
    %determine if there is more than one block
    if nob<=1
        part(1,1)=1;
        part(1,2)=length(clone(jj-2).mat);
    end

    if nob>1
        for i=1:nob
            if i~=nob 
                part(i,1)=(i-1)*column+1;
                part(i,2)=column*i;
            else
                part(i,1)=(i-1)*column+1;
                part(i,2)=length(clone(jj-2).mat);
            end
        end
    end

    %set block number
    b=blocknumber;
    
    %determine if the block number too large or negative
    if (b>nob)
        disp("Maximum Block Number"),disp(nob); 
    return
    end
    if (b<0)
        disp("Block Number Can Not Be Negative");
    return
    end

    %get the required block
    blockog=clone(jj-2).mat(row(1,1):row(1,2),part(b,1):part(b,2));
    rowsize=row(1,2)-row(1,1)+1;    

    %delete unneeded zeros
    blockp(jj-2).file=clone(jj-2).file;
    blockp(jj-2).mat=blockog;
    remainzero=zeros(1,11);

    c=1;
    del=0;
    while c<=numel(blockp(jj-2).mat)
        if ((c+numel(remainzero)-1)<numel(blockp(jj-2).mat))
            if (blockp(jj-2).mat(c:c+numel(remainzero)-1)==remainzero)
                c=c+numel(remainzero);
                del=1;
            end
        end
    while c<=numel(blockp(jj-2).mat)
        if (del==1)&&(blockp(jj-2).mat(c)==0) 
            blockp(jj-2).mat(c)=[];    
        else
            break
        end
    end
    del=0;
    c=c+1;
    end
    
    %next data file
    jj=jj+1;
end

block=[];
%combine the different block
for i=1:qty
    block=horzcat(block,blockp(i).mat);
end

%block=[0 0 8 1 2 1 7 2 7 0 1 10 3 3 26 81 73458 239 12 23 349 12902 2349 230 120 2931 29 9348 93 943 124 34 0450 9435 3 21 2 4 5 6 1 10 10 21 2 73 3 7 2 8 2 38 2 9 3 4 9 2 100 0 8 1 2 1 7 2 7 0 1 10 3 3 26 81 73 21 2 4 5 6 1 10 10 21 2 73 3 7 2 8 2 ];

%readjust column size to the new combined block 
column=numel(block);

%find unique symbols of the matrix
symbol=unique(block);


%=========================================================================%
ordofe=orderofconditionalentropy;
H_ALL=zeros(ordofe+1,1);

for n=1:ordofe
  
    %find every symbol's position
    position=zeros(1,numel(block));
    for i=1:numel(block)
        position(i)=find(symbol==block(i));
    end    

    %build symbol combination of n.order
    combination_n=struct;

    k=1;
    while k<=numel(block)-n
        combination_n(k).combination=block(1,k:k+n);
        combination_n(k).position=0;
        for i=1:n
            combination_n(k).position=combination_n(k).position+((position(k+i-1)-1)*numel(symbol)^(n-i+1));
        end
        combination_n(k).position=combination_n(k).position+position(k+n);
        k=k+1;
    end

    buffermat_n=[combination_n.position];
    [symbol_n,~,pos_n]=unique(buffermat_n);
    occurrences_n=[symbol_n' accumarray(pos_n,1)];
    occurrences_n_spec=occurrences_n;

    %build symbol combination of (n-1).order
    combination_nm=struct;

    k=1;
    while k<=numel(block)-n+1
        combination_nm(k).combination=block(1,k:k+n-1);
        combination_nm(k).position=0;
        for i=1:n-1
            combination_nm(k).position=combination_nm(k).position+((position(k+i-1)-1)*numel(symbol)^(n-i));
        end
        combination_nm(k).position=combination_nm(k).position+position(k+n-1);
        k=k+1;
    end

    buffermat_nm=[combination_nm.position];
    [symbol_nm,~,pos_nm]=unique(buffermat_nm);
    occurrences_nm=[symbol_nm' accumarray(pos_nm,1)];

    %calculate Hy for (n-1).order
    occnmsize=size(occurrences_nm);
    Hy=zeros(occnmsize(1,1),1);
    if mod(occurrences_n(1,1),numel(symbol))~=0
        k=idivide(uint64(occurrences_n(1,1)),uint64(numel(symbol)))+1;
    else
        k=idivide(uint64(occurrences_n(1,1)),uint64(numel(symbol)));
    end
    j=1;
    while ~isempty(occurrences_n)
        occnr=find(occurrences_n(:,1)<=numel(symbol)*k,1,'last');
        if ~isempty(occnr)     
            occsum=sum(occurrences_n(1:occnr,2));
            p_n=zeros(occnr,1);
            for i=1:occnr
                p_n(i)=occurrences_n(i,2)/occsum;
            end
            Hy(j)=-dot(p_n,log2(p_n));
            j=j+1;
            restsize=size(occurrences_n);
            occurrences_n=occurrences_n(occnr+1:restsize(1,1),:);
            if ~isempty(occurrences_n)
                if mod(occurrences_n(1,1),numel(symbol))~=0
                    k=idivide(uint64(occurrences_n(1,1)),uint64(numel(symbol)))+1;
                else
                    k=idivide(uint64(occurrences_n(1,1)),uint64(numel(symbol)));
                end
            end
        end
    end

    %calculate probabilities of (n-1).order
    p_nm=zeros(occnmsize(1,1),2);
    for k=1:occnmsize(1,1)
        p_nm(k,1)=occurrences_nm(k,1);
        p_nm(k,2)=occurrences_nm(k,2)/numel(buffermat_nm);
    end
    
    %final calculation of conditional entropy
    if n==1
        H_ALL(1)=-dot(p_nm(:,2),log2(p_nm(:,2)));
    end
    H_ALL(n+1)=sum(Hy.*p_nm(:,2));
end

% %store result
% entropy=fullfile(dest,'conditional_entropy');
% save(entropy,'H_ALL');
% 
% allvar=fullfile(dest,'workspace');
% save(allvar);

plot(0:1:ordofe,H_ALL,'--*');
xticks(0:1:ordofe);
title('conditional entropy');
xlabel('n.order');
ylabel('value');
