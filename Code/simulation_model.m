function [Bank,Emiss,MF, indx1, indx2] = simulation_model(Prod, DE, RF, BankSize1, N, Nyears,MF1,LTsamps,ppt_to_tonnes)

t = 1;

indx1 = randi(size(DE,2),N,1);
DEsamps = DE(t,indx1)';
RFsamps = RF(t,indx1)';

indx2 = randi(size(Prod,1),N,1);
Prodsamps =Prod(indx2,t);

[bank_samps,emiss_samps] = Bank_Emiss(Prodsamps,RFsamps,DEsamps,BankSize1*ones(size(Prodsamps)));
Bank(:,t) = bank_samps;
Emiss(:,t) = emiss_samps;
MF(:,t+1) = MF1*exp(-LTsamps(:,t))+(1/ppt_to_tonnes).*Emiss(:,t);
for t = 2:Nyears
    DEsamps = DE(t,indx1)';
    RFsamps = RF(t,indx1)';
    Prodsamps = Prod(indx2,t);
    [bank_samps,emiss_samps] = Bank_Emiss(Prodsamps,RFsamps,DEsamps,bank_samps);
    Bank(:,t) = bank_samps;
    Emiss(:,t) = emiss_samps; 
    MF(:,t+1) = MF(:,t).*exp(-LTsamps(:,t))+(1/ppt_to_tonnes).*Emiss(:,t);
end