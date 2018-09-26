function SS=frontierspendingrules(S)
    S.Gvec= S.g0*S.Geff; 
    S.Gvec(1)= S.g2*S.Geff;S.Gvec(end)= S.g2*S.Geff;
    S.Tvec = S.Gvec.*( 0.5 + S.that*(S.Pvec-0.5) );
    S.T2vec= S.Gvec - S.Tvec;
    S.fracWP0=-S.Tvec(1)/S.Cwp;
    SS=S;
end