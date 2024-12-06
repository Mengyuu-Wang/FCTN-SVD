function[ticktock,node,R] = FCTN_svd(X,Tol)
tstart=tic;
dims = size(X);
C = X;
n = numel(dims);
ep=Tol/sqrt(n);
node=cell(1,n);
R = ones(n);

for i=1:n-1
    if i==1
        C = reshape(C,dims(i),numel(C)/dims(i));
        [U,S,V] = svd(C,'econ');
        S = diag(S);
        rc = my_chop2(S,sqrt(n-1)*ep*norm(S));
        R(1,2:n) = somerule(rc,1,n);
        U = U(:,1:prod(R(1,2:n)));
        dd = horzcat(dims(1), R(1,2:n));
        node{1} = reshape(U,dd);
        S = S(1:prod(R(1,2:n)));
        V = V(:,1:prod(R(1,2:n)));
        V = V*diag(S);
        posi = horzcat(dims(2:end) , R(1,2:end));
        C = reshape(V,posi);
        dd = horzcat(n:2*n-2,1:n-1);
        for j = 1:2*n-2
            if mod(j,2) == 1
                ddd(j) = dd((j+1)/2);
            else
                ddd(j) = dd(j/2-1+n);
            end
        end
        C = permute(C,ddd);
    else
        m = prod(R(1:i-1,i))*dims(i);
        C = reshape(C,m,numel(C)/m);
        [U,S,V] = svd(C,'econ');
        S = diag(S);
        rc = my_chop2(S,sqrt(n-i)*ep*norm(S));
        R(i,i+1:n) = somerule(rc,i,n); 
        U = U(:,1:prod(R(i,i+1:n)));
        posi = horzcat(R(1:i-1,i)',dims(i),R(i,i+1:n));
        node{i} = reshape(U,posi);
        S = S(1:prod(R(i,i+1:n)));
        V = V(:,1:prod(R(i,i+1:n)));
        V = V*diag(S);
        for k = 1:n-i
            dddd(k) = prod(R(1:i-1,i+k));
        end
        dddd = dddd(1:n-i);
        dd = horzcat(dddd,dims(i+1:end));
        for l = 1:2*n-2*i
            if mod(l,2) == 1
                ddd(l) = dd((l+1)/2);
            else
                ddd(l) = dd(l/2+n-i);
            end
        end
        ddd = ddd(1:2*n-2*i);
        ddd = horzcat(ddd,R(i,i+1:n));
        C = reshape(V,ddd);
        dd = 1:3*n-3*i;
        for h = 1:3*n-3*i
            if mod(h,3) == 1
                ddd(h) = dd(2*h/3+1/3);
            elseif mod(h,3) == 2
                ddd(h) = dd((h+1)/3+2*n-2*i);
            else 
                ddd(h) = dd(2*h/3);
            end
        end
        ddd = ddd(1:3*n-3*i);
        C = permute(C,ddd);
    end
end
D = horzcat(R(1:n-1,n)',dims(n));
node{n} = reshape(C,D);
ticktock = toc(tstart);

end

