clc; clear; close all;

%% Parameters
n = 100;       
p = 0.3;       
rng(1);        

%% Generate ER graph adjacency matrix
A = rand(n) < p;     
A = triu(A,1);        
A = A + A';           

X_ER = A(:).';   

counts = [sum(X_ER==0), sum(X_ER==1)];
prob = counts / sum(counts);

fprintf('\nER graph symbol probabilities:\n');
fprintf('P(0) = %.4f\n', prob(1));
fprintf('P(1) = %.4f\n', prob(2));

% Entropy
H_ER = -sum(prob(prob>0).*log2(prob(prob>0)));
fprintf('Entropy per bit: %.4f bits\n', H_ER);


% Even sequence length
if mod(length(X_ER),2) ~= 0
    X_ER = X_ER(1:end-1);
end

% Form pairs
pairs = reshape(X_ER,2,[]).';
symbols = pairs(:,1)*2 + pairs(:,2);   % 00=0, 01=1, 10=2, 11=3

% Probabilities
counts = histcounts(symbols,0:4);
prob = counts / sum(counts);

fprintf('\nER graph pair probabilities:\n');
fprintf('P(00)=%.4f, P(01)=%.4f, P(10)=%.4f, P(11)=%.4f\n', prob);

% Entropy of pairs
H_pair = -sum(prob(prob>0).*log2(prob(prob>0)));
fprintf('Entropy per pair: %.4f bits, per bit: %.4f bits\n', H_pair, H_pair/2);


% Sort probabilities
[prob_sorted, idx] = sort(prob,'descend');

% Defined Huffman codes
codes = cell(1,4);
codes{idx(1)} = '0';
codes{idx(2)} = '10';
codes{idx(3)} = '110';
codes{idx(4)} = '111';

% Encode
encoded_bits = '';
for i = 1:length(symbols)
    encoded_bits = [encoded_bits codes{symbols(i)+1}];
end

compressed_bits = length(encoded_bits);
compression_ratio = length(X_ER)/compressed_bits;

fprintf('Average code length = %.4f bits/pair\n', sum(prob.*cellfun(@length,codes)));
fprintf('Compression ratio = %.4f\n', compression_ratio);


% Decode map
decodeMap = containers.Map(codes, 0:3);

decoded_symbols = [];
buffer = '';
for i = 1:length(encoded_bits)
    buffer = [buffer encoded_bits(i)];
    if isKey(decodeMap,buffer)
        decoded_symbols(end+1) = decodeMap(buffer); 
        buffer = '';
    end
end

% Recovering adjacency sequence
decoded_pairs = zeros(length(decoded_symbols),2);
for i = 1:length(decoded_symbols)
    val = decoded_symbols(i);
    decoded_pairs(i,:) = [floor(val/2), mod(val,2)];
end
X_decoded = reshape(decoded_pairs.',1,[]);

if isequal(X_ER, X_decoded)
    disp('ER graph lossless decoding verified.');
else
    disp('Decoding error!');
end

G = graph(A);        
figure;
plot(G);
title('Erdős–Rényi Random Graph');

symbols = {'00','01','10','11'};
probabilities = prob;

figure;
bar(probabilities);
set(gca, 'XTickLabel', symbols);
xlabel('Symbol pair');
ylabel('Probability');
title('ER Graph Symbol Pair Probabilities');
grid on;


p_values = 0.1:0.1:0.9;
compression_ratios = zeros(size(p_values));

for k = 1:length(p_values)
    p = p_values(k);
    A = rand(n) < p;
    A = triu(A,1);
    A = A + A';
    
   
    X_ER = A(:).';
    % ...Huffman code...
    % assign compression_ratios(k) = computed ratio
end

figure;
plot(p_values, compression_ratios,'-o');
xlabel('ER Graph Edge Probability p');
ylabel('Compression Ratio');
title('Compression Ratio vs ER Graph Density');
grid on;

