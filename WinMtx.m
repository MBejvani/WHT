function [ W ] = WinMtx( N,sd )
%WINMTX creates a N*N matrix encompassing nomalized adaptive Gaussian windows in each columns  
% Input:
% N          - length of the signal
% sd         - vector of length or standard deviation of the window in each culumn 
% Output:
% W          - the Matrix of adaptive the windows

lsd = length(sd);
if lsd==1
    sd = ones(1,N)*sd;
elseif lsd~=1 && lsd~=N
    error('window-length vector error');
end
if nargin<2
    sd = ones(1,N)*15;
end
W = zeros(N);
for i = 1:N
    W(:,i) = exp(-1/2*(((1:N)-i+1)/sd(i)).^2);
end
SW = 1./sum(W);
W  = W .* SW; % Normalization
end

% function [ W ] = winmtx0( N,sd )
% %WINMTX0 produces redundant coeficients at the corners
% W = zeros(N);
% if length(sd) == 1
%     W = toeplitz(fftshift(1/(sqrt(2*pi)*sd) * ...
%         exp(-1/2*((-fix(N/2):fix(N/2)-1)/sd).^2)));
% else 
%     for i = 1:N
%         W(:,i) = circshift(1/(sqrt(2*pi)*sd(i)) * ...
%             exp(-1/2*((-fix(N/2):fix(N/2)-1)/sd(i)).^2)',fix(N/2)+i);
%     end
% end
