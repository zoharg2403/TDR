clear all;
close all;

load('ORF.mat');
load('RC_profiles.mat');

% get NFC values for each codon on each gene sequence
NFCval = NFC(ORF,RC_profiles);

% validation of NFC with hist 
NFCmap = codon_NFC_map(ORF, NFCval);

% removing empty field values from map
map_fields = fieldnames(NFCmap);
for i = 1:length(map_fields)
    codonNT = char(map_fields(i));
    if isempty(NFCmap.(codonNT))
       NFCmap = rmfield(NFCmap,(char(map_fields(i))));
    end 
end

% plotting histograms
map_fields = fieldnames(NFCmap);
for i = 1:length(map_fields)
    codon_NFC = NFCmap.(char(map_fields(i)));
    subplot (8,8,i);
    histogram (codon_NFC,'Normalization','probability')
    title (strcat(char(map_fields(i)),' NFC histogram'));
    xlabel('NFC');
end

% get MLE for each codon , mu and sigma accordingly
[MLE_map, mu_map, sigma_map] = NFC_MLE(NFCmap);

% get TDR value for each codon from mu of mle dist
TDR_map = codon_TDR (MLE_map);

% corrolation plot between tAI and TDR
TDR = cell2mat(struct2cell(TDR_map));
A = xlsread('tAI_tCN.xlsx','B3:B66');
A_fixed = A(~isnan(A));
vec = ones(61,1);
tAI = vec./A_fixed;
r = corr(tAI,TDR);
figure;
scatter(tAI, TDR);
xlabel('tAI');
ylabel('TDR');




    
        
        
  