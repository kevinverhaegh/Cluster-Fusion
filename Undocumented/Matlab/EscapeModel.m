function output=EscapeModel(input)
output=zeros(numel(input(1,:)),1);
for i=1:numel(input(1,:))
    [S,~,T]=spectrogram(input(:,i)/max(input(:,i)), hann(16),0,32,[1e5]);
    ind=ReleaseTimer(S(end,:),100,1e-2);
    output(i) = int64(T(ind)*1e5);
end