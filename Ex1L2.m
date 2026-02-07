clc;
clear;
close all;

%% Parameters
rng(1);
N = 1e5;
p = 0.3;

%% Bernoulli Source
X = rand(1,N) < p;

% Ensure even length
if mod(N,2) ~= 0
    X = X(1:end-1);
end

%% Form symbol pairs
pairs = reshape(X, 2, []).';
symbols = pairs(:,1)*2 + pairs(:,2);   % 00=0, 01=1, 10=2, 11=3

%% Probabilities
counts = histcounts(symbols, 0:4);
prob = counts / sum(counts);

%% Entropy
H = -sum(prob(prob>0).*log2(prob(prob>0)));

%% Manual Huffman Code (optimal for 4 symbols)
% Sort probabilities
[prob_sorted, idx] = sort(prob, 'descend');

% Fixed optimal Huffman structure for 4 symbols
codes = cell(1,4);
codes{idx(1)} = '0';
codes{idx(2)} = '10';
codes{idx(3)} = '110';
codes{idx(4)} = '111';

%% Encoding
encoded_bits = '';
for i = 1:length(symbols)
    encoded_bits = [encoded_bits codes{symbols(i)+1}];
end

%% Metrics
original_bits = length(X);
compressed_bits = length(encoded_bits);
compression_ratio = original_bits / compressed_bits;

avg_length = sum(prob .* cellfun(@length, codes));

%% Display
disp('--- Results ---')
Symbol = {'00'; '01'; '10'; '11'};
Probability = prob(:);

T = table(Symbol, Probability);
disp(T);
fprintf('Entropy = %.4f bits/symbol\n', H);
fprintf('Average code length = %.4f bits/symbol\n', avg_length);
fprintf('Compression ratio = %.4f\n', compression_ratio);

fprintf('\nEntropy per symbol pair: %.4f bits\n', H);
fprintf('Entropy per bit: %.4f bits\n', H/2);
fprintf('Average Huffman code length: %.4f bits/symbol\n', avg_length);




% Build decoding dictionary
decodeMap = containers.Map(codes, 0:3);

decoded_symbols = [];
buffer = '';

for i = 1:length(encoded_bits)
    buffer = [buffer encoded_bits(i)];
    if isKey(decodeMap, buffer)
        decoded_symbols(end+1) = decodeMap(buffer); %#ok<SAGROW>
        buffer = '';
    end
end

decoded_pairs = zeros(length(decoded_symbols),2);

for i = 1:length(decoded_symbols)
    val = decoded_symbols(i);
    decoded_pairs(i,:) = [floor(val/2), mod(val,2)];
end

X_decoded = reshape(decoded_pairs.',1,[]);

if isequal(X, X_decoded)
    disp('Lossless decoding verified.');
else
    disp('Decoding error!');
end
