function S = DinosaurBoneSVD_v4(A);
% S = DinosaurBoneSVD_v4(A)
%
% Input parameter:
%     A          MxMxL1  polynomial matrix
%
% Output parameter:
%     S          L2xM matrix containing one time-domain singular value per column

%--------------------------------------------------------------
%  parameters --- potentially to be defined externally later
%--------------------------------------------------------------
Nfft = 1024; N = size(A,1);
MinRunLength = 16;
Support = 2;

%-----------------------------------------------------------------------
% (1) bin-wise singular values of A  
%-----------------------------------------------------------------------
Af = fft(A, Nfft, 3);
Sf = zeros(N, Nfft);
for k = 1:Nfft
    [~, dummy, ~] = svd(Af(:, :, k));
    Sf(:, k) = diag(dummy);
end  

%-----------------------------------------------------------------------
% (2) divide singular values into well-separated segments of minimum length
%-----------------------------------------------------------------------
D = min(-diff([Sf; -flipud(Sf)], 1), [], 1); 
Thresh = max(D)/5;

save DinosaurBoneInterimFig6.mat D Thresh

DD = (D > Thresh);
SegEnds = find(diff([DD DD(1)]) < -eps);
SegStarts = find(diff([DD(end) DD]) > eps);
if SegEnds(1) < SegStarts(1)
    SegmentIndices = [SegStarts' [SegEnds(2:end) SegEnds(1)+Nfft]'];
else
    SegmentIndices = [SegStarts' SegEnds'];
end

% purge short segments
SegmentLengths = diff(SegmentIndices, [], 2);
SegmentIndices = SegmentIndices(SegmentLengths > MinRunLength, :);
Q = size(SegmentIndices, 1);
disp(sprintf('retained %d segments of sufficient length', Q));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%-----------------------------------------------------------------------
% (3) individual reconstruction of segments and matching via error norms
%-----------------------------------------------------------------------
M = size(Sf, 1); 
N = size(Sf, 2);
t = (-Support:Support)';
Omega = SegmentIndices;
St = zeros(M, length(t), Q);
for q = 1:Q
    % extract qth segment 
    Sf_q = Sf(:, mod((Omega(q, 1):Omega(q, 2)) - 1, Nfft) + 1);
    % frequency bins of qth segment
    omega_q = ((Omega(q, 1):Omega(q, 2))' - 1) / Nfft * 2 * pi;
    % time domain reconstruction
    St(:, :, q) = inverseFourierTransformVector(Sf_q.', omega_q, t).';
    % alternative time domain reconstruction with a half sample delay
    Sf_q = Sf_q.* ( ones(2,1)*exp(-1i*omega_q'/2) );
    Sf_HalfShiftCheck = inverseFourierTransformVector(Sf_q.', omega_q, t).';
    % check if there are any reconstructions in St2 which are symmetric to (Support+.5);
    % these belong to segments of singular values with an odd number of zero crossings
    for m = 1:M,
        % check symmetry
        %
        % this needs to be further looked at --- I only consider even symmetry, but 
        % could it also be odd?
        dummy = Sf_HalfShiftCheck(m,1:end);
        if norm( dummy(2:end) - conj(dummy(end:-1:2)) )/norm(dummy) < 0.1,
           % approximate symmetry established, replace in St(:,:,q)
           St(m,:,q) = dummy;
        end;
    end;       
end

save DinosaurBoneInterimResults.mat Sf St Omega

% Align segments using error norm 
t2 = (-2:2) + Support + 1;
for q = 2:Q
   A1 = St(:, t2, q - 1);
   A2 = St(:, t2, q);
   Cost = zeros(M);
   for i = 1:M,
      for j = 1:M,
         % the following min() takes the sign ambiguity into account 
         Cost(i, j) = min( norm(A1(i,:) - A2(j,:))^2 , norm(A1(i,:) + A2(j,:))^2 );
      end;
   end;

   % [AssignVec, ~] = munkres(Cost);     
   % replaced by matlab built-in functions thanks to Sebastian Schlecht, 28/7/25
   matched = matchpairs(Cost, max(Cost(:)));
   AssignVec = accumarray(matched, 1, size(Cost));
   
   PermIndex = zeros(1, M);

   %  assignment matrix to index vector
   AssignMtx = zeros(1, size(AssignVec, 2));
   for j = 1:size(AssignMtx, 2),
      row = find(AssignVec(:, j));
      if ~isempty(row),
         AssignMtx(j) = row; % Column j assigned to row
      else
         AssignMtx(j) = 0;   % No assignment found
      end;
   end;

    for i = 1:M,
        PermIndex(AssignMtx(i)) = i;
    end;

    Ps = sign(diag(real(A1 * A2(PermIndex,:)')));  
    % need to check,I assumed after permutation the diagonal contain
    %  the sign of the aligned
    St(:, :, q) = diag(Ps) * St(PermIndex, :, q);    
end;

%-----------------------------------------------------------------------
% (4) reconstruction --- segment-size weighted average
%-----------------------------------------------------------------------
W = (Omega(:, 2) - Omega(:, 1));
W = W / sum(W);
S = zeros(M, 2 * Support + 1);
for q = 1:Q
    S = S + St(:, :, q) * W(q);
end

