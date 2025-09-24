function fval=NLDhopf2a_noglobal(x,P)



% NSUB=32;
% TR=2;
% 
% omega = repmat(2*pi*f_diff',1,2); omega(:,1) = -omega(:,1);
dt=0.1*P.TR/2;
Tmax=P.NSUB*P.TSmax;
sig=0.02;
dsig = sqrt(dt)*P.sig; % to avoid sqrt(dt) at each time step
% % 
% we=0.6030;
% C=squeeze(Coptim(find(abs(WE-we)<0.0001),:,:));
% wC = we*C;
% sumC = repmat(sum(wC,2),1,2); 



a=repmat(x(:),1,2);
FC_simul_iter=zeros(20,P.N,P.N);
for iter=1:20
%     display(iter)

        xs=zeros(P.Tmax,P.N);
        z = 0.1*ones(P.N,2); % --> x = z(:,1), y = z(:,2)
        nn=0;
        % discard first 3000 time steps
        for t=0:dt:3000
            suma = P.wC*z - P.sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
            zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
            z = z + dt*(a.*z + zz.*P.omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(P.N,2);
        end
        % actual modeling (x=BOLD signal (Interpretation), y some other oscillation)
        for t=0:dt:((Tmax-1)*P.TR)
            suma = P.wC*z - P.sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
            zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
            z = z + dt*(a.*z + zz.*P.omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(P.N,2);
            if abs(mod(t,P.TR))<0.01
                nn=nn+1;
                xs(nn,:)=z(:,1)';
            end
        end

        FC_simul_iter(iter,:,:)=corrcoef(xs(1:nn,:));   

end

FC_simul=squeeze(mean(FC_simul_iter,1));
xco=1-corrcoef(atanh(P.FCemp(P.Isubdiag)),atanh(FC_simul(P.Isubdiag)));
display(xco(2))
fval=xco(2);
