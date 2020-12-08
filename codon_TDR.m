function [TDR] = codon_TDR(MLE_map)
codon_names = fieldnames(MLE_map);

% initializing TDR map with field names
TDR = ([]);

for i = 1:length(MLE_map)
    codonNT = char(codon_names(i));
    TDR.(codonNT) = [];
end

% calc TDR for each codon
for i = 1:length(codon_names)
    codon = char(codon_names(i));
    codonMU = MLE_map.(codon)(1);
    codonSIGMA = MLE_map.(codon)(2);
    codonX = 0:0.001:1;
    log_normal_value = log_normal_pdf(codonX, codonMU, codonSIGMA);
    TDR.(codon) = skewness(log_normal_value);  
end

% normalize TDR values
max_TDR_val = max(cell2mat(struct2cell(TDR)));
for i = 1:length(fieldnames(TDR))
    codon = char(codon_names(i));
    new_TDR = (TDR.(codon))/max_TDR_val;
    TDR.(codon) = new_TDR;
end

end

