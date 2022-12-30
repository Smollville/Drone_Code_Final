function [Fbi,Gbii,Gbij,Qbi,Rbi] = GetBatchYdistMatrices(Ai,Bii,Bij,Ci,N,Pi,Qi,Ri,alphai)
% Returns matrices for batch computation of state sequence and cost 
% functional

    F = [];
    Gii = [];
    Gij = [];
    H = [];
    Qbi = [];
    Rbi = [];
    for n = 0:N
        % state matrices
        F = [F ; Ai^n];
        Gaii = [];
        Gaij = [];
        for m = 0:N-1
            ni = n-m-1;
            Gauxii = Ai^ni*Bii;
            Gauxij = Ai^ni*Bij;
            if ni < 0
                Gauxii = Gauxii*0;
                Gauxij = Gauxij*0;
            end
            Gaii = [Gaii , Gauxii];
            Gaij = [Gaij , Gauxij];
        end
        Gii = [Gii ; Gaii];
        Gij = [Gij ; Gaij];
        H = blkdiag(H,Ci);

        % cost matrices
        if exist('Pi','var') && exist('Qi','var') && exist('Ri','var')
            if n < N
                Qbi = alphai * blkdiag(Qbi,Qi);
                Rbi = alphai * blkdiag(Rbi,Ri);
            else
                Qbi = blkdiag(Qbi,Pi);
            end
        end
    end

    % compute output matrices
    Fbi = H*F;
    Gbii = H*Gii;
    Gbij = H*Gij;

end

