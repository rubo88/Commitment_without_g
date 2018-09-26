function SS=linearrules_eff(S)
    S.Gvec= (S.g0 + 4*S.g2*(S.Pvec-0.5).^2)*S.Geff; 
    S.Tvec = S.Geff.*( 0.5 + S.that*(S.Pvec-0.5) )+(S.Gvec-S.Geff)/2;
    S.T2vec= S.Gvec - S.Tvec;
    S.fracWP0=-S.Tvec(1)/S.Cwp;
    SS=S;
end