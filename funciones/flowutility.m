function [U,U2]=flowutility(S)
U =log(S.C ) - (S.C+S.C2+S.Gvec)/S.rho;
U2=log(S.C2) - (S.C+S.C2+S.Gvec)/S.rho;
end