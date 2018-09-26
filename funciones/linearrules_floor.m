function SS=linearrules_floor(S)
    S.Gvec= (S.g0 + 4*S.g2*(S.Pvec-0.5).^2)*S.Geff; 
    S.Tvec = S.Gvec.*( 0.5 + S.that*(S.Pvec-0.5) );
    S.T2vec= S.Gvec - S.Tvec;
    S.bailouted =0;
    aux     =S.Tvec  >-S.bail*S.Cwp;
    aux2    =S.T2vec >-S.bail*S.Cwp;
    S.Tvec(aux)     =-S.bail*S.Cwp;
    S.T2vec(aux2)   =-S.bail*S.Cwp;
    S.T2vec(aux)    =S.bail*S.Cwp+S.Gvec(aux);
    S.Tvec(aux2)    =S.bail*S.Cwp+S.Gvec(aux2);
    S.bailouted     =aux+aux2;
    S.fracWP0=-S.Tvec(1)/S.Cwp;
    
    S.aux=aux;S.aux2=aux2;
    SS=S;
end