function [Xs, Ps] = RTS(Phi, Xf, Pf, Xf1, Pf1)
    Xs = Xf;  Ps = Pf;
    for k=length(Phi)-1:-1:1
        Ks = Pf{k}*Phi{k}'*Pf1{k}^-1;
        Xs{k} = Xf{k} + Ks*(Xs{k+1} - Xf1{k+1});
        Ps{k} = Pf{k} + Ks*(Ps{k+1} - Pf1{k+1})*Ks';
    end
