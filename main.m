%Calculate conditional entropy

clear;

%set file path
quantityoffiles=1;
filelist=dir('/home/muyan/Documents/masterarbeit/ACE_AUTO');
dest='/home/muyan/Documents/%%%%%/';  

%define block size, set block number
rownr=[1 1];        %[startrownr endrownr]
%columnsize=100;    %comment if data is the whole matrix
blocknumber=1;      %set 1 if data is the whole matrix
orderofconditionalentropy=0;
format long;

%select random files
qty=quantityoffiles;
randf=randperm(numel(filelist),qty);

for ii=1:qty
    while (randf(ii)==1)||(randf(ii)==2)
        randf(ii)=randperm(numel(filelist),1);
    end
end

%copy the selected files to destination folder
for kk=1:qty
    source=fullfile(filelist(randf(kk)).folder,filelist(randf(kk)).name);
    copyfile(source,dest);
end

%initialization for result struct
result=struct();

%get filelist of the current destination folder
cfiles=dir(dest);

%cfiles(1) and cfiles(2) are '.' and '..', file name is stored from cfiles(3)
jj=3;
while jj<=qty+2

    %load data file
    load(fullfile(cfiles(jj).folder,cfiles(jj).name));
  
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
    
    %get order of conditional entrendopy
    nn=orderofconditionalentropy;
    
    %reset all variables(which are maybe necessary)
    H=zeros(1,nn+1);
    Hcy=[];
    comb=[];
    comby=[];

    %initialization, convert NaN to 0
    clone=m0filesat;
    clone(isnan(clone))=0;

    %define block size and split the matrix {clone} into smaller blocks 
    %nob=number of blocks
    row=rownr;
    
    %column=columnsize;  %use this if data is only a part of the matrix
    buffer=size(mat);   %use this if data is the whole matrix
    column=buffer(2);   %use this if data is the whole matrix
    nob=fix(length(clone)/column);

    if nob~=length(clone)/column
        nob=nob+1;
    end
    
    %preallocate array(part) size
    part=zeros(nob,2);
    
    %determine if there is more than one block
    if nob<=1
        part(1,1)=1;
        part(1,2)=length(clone);
    end

    if nob>1
        for i=1:nob
            if i~=nob 
                part(i,1)=(i-1)*column+1;
                part(i,2)=column*i;
            else
                part(i,1)=(i-1)*column+1;
                part(i,2)=length(clone);
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
    blockog=clone(row(1,1):row(1,2),part(b,1):part(b,2));
    rowsize=row(1,2)-row(1,1)+1;

    %delete unneeded zeros
    block=blockog;
    remainzero=zeros(1,11);

    c=1;
    del=0;
    while c<=numel(block)
        if ((c+numel(remainzero)-1)<numel(block))
            if (block(c:c+numel(remainzero)-1)==remainzero)
                c=c+numel(remainzero);
                del=1;
            end
        end
    while c<=numel(block)
        if (del==1)&&(block(c)==0) 
            block(c)=[];
        else
            break
        end
    end
    del=0;
    c=c+1;
    end
    
    %readjust column size to the new block which dosen't contain unneeded
    %zero
    column=numel(block);

    %find unique symbols of the matrix
    symbol=unique(block);

    %-----------------------------------------%
    %              Main(entropy)              %
    %-----------------------------------------%
    if nn<0
        disp("order of conditional entropy must greater than or equal to 0"); 
    end
    
    %Calculate 0.order entropy
    if nn>=0 
        
        %Calculate probability of every different symbol, store in pa(k)
        counter=zeros(1,numel(symbol));
        pa=zeros(numel(symbol));
        for k=1:numel(symbol)
            for i=1:rowsize
                for j=1:column
                    if block(i,j)==symbol(k)
                        counter(k)=counter(k)+1;
                    end
                end
            end
            blocksize=size(block);
            pa(k)=counter(k)/(blocksize(1)*blocksize(2));
        end
  
        %Calculate 0.order entropy of the block
        H(1)=-dot(pa,log2(pa));
    end

    %Calculate n.order conditional entropy of the block
    if nn>0
        n=1;
        while n<=nn 
            comby=combination(n,symbol);
            comb=combination(n+1,symbol); 
    
            %Calculate the probabilities of n-1 and n symbols combinations
            %Occurrences are stored in comby(i,n) and comb(i,n+1)
            %Probabilities are stored in comby(i,n+1) and comb(i,n+2)   
            for l=1:numel(symbol)^n
              comby(l,n+1)=0;
              for m=1:column-(n-1)
                pblock=block(:,m:m+n-1);
                for i=1:rowsize 
                  if comby(l,1:n)==pblock(i,1:n)
                    comby(l,n+1)=comby(l,n+1)+1;
                  end
                end
              end
            end

            for i=1:numel(symbol)^n
              comby(i,n+2)=comby(i,n+1)/(rowsize*(column-(n-1)));
            end
    
    
    
    
            for l=1:numel(symbol)^(n+1)
              comb(l,n+2)=0;
              for m=1:column-n
                pblock=block(:,m:m+n);
                for i=1:rowsize 
                  if comb(l,1:n+1)==pblock(i,1:n+1)
                    comb(l,n+2)=comb(l,n+2)+1;
                  end
                end
              end
            end

            sn=1;
            while sn<=numel(symbol)^(n+1)
              su=sum(comb(sn:sn+numel(symbol)-1,n+2)); 
              for i=sn:sn+numel(symbol)-1
                if su~=0
                  comb(i,n+3)=comb(i,n+2)/su;
                end
              end
              su=0;
              sn=sn+numel(symbol);
            end
            comb(isnan(comb))=0; 

    
            %Calculate the n.order conditional entropy
            sn=1;
            f=1;
            Hcy=zeros(numel(symbol)^n);
            while sn<=numel(symbol)^(n+1)
                pc=zeros(numel(symbol));
                for i=sn:sn+numel(symbol)-1
                    pc(i)=comb(i,n+3);
                end
                pc=nonzeros(pc);
                Hcy(f)=-dot(pc,log2(pc));
                pc=0;
                f=f+1;
                sn=sn+numel(symbol);
            end
    
            H(n+1)=0;
            for i=1:numel(symbol)^n
                H(n+1)=H(n+1)+comby(i,n+2)*Hcy(i);
            end
            
            %next higher order
            n=n+1;
        end
    end

    %store result
    result(jj-2).name=cfiles(jj).name;
    result(jj-2).entropy=H;
    
    %next data file
    jj=jj+1;
end

%save result in destination folder
matfile=fullfile(dest,'result');
save(matfile,'result');
  

%Get all combinations of unique symbols (n and n-1 symbols)
function comb=combination(n,symbol)
    m=n;
    e=0;
    i=1;
    j=1;
    c=1;
    comb=zeros(numel(symbol)^n,n);
    while m~=0
        comb(i,m)=symbol(j);
        i=i+1;
        c=c+1;   
        if c>numel(symbol)^e 
            j=j+1;
            c=1;
        end
    
        if j>numel(symbol)
            j=1;
        end
    
        if i>numel(symbol)^n
            i=1;
            e=e+1;
            m=m-1;
        end   
    end
end
